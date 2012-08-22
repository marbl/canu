/*
Copyright (C) 2011, Battelle National Biodefense Institute (BNBI);
all rights reserved. Authored by: Sergey Koren

This Software was prepared for the Department of Homeland Security
(DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
part of contract HSHQDC-07-C-00020 to manage and operate the National
Biodefense Analysis and Countermeasures Center (NBACC), a Federally
Funded Research and Development Center.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

 * Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

 * Neither the name of the Battelle National Biodefense Institute nor
  the names of its contributors may be used to endorse or promote
  products derived from this software without specific prior written
  permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

using namespace std;

#include "AS_PBR_filter.hh"
#include "AS_PBR_store.hh"

#include "AS_OVS_overlapStore.h"

#include <set>
#include <map>

static const char *rcsid_AS_PBR_FILTER_C = "$Id: AS_PBR_filter.cc,v 1.3 2012-08-22 14:41:00 skoren Exp $";

/**
 * This function is responsible for computing the repeat threshold based on the number of PacBio sequences each short-read maps to
 */
void computeRepeatThreshold(PBRThreadGlobals &thread_globals, AS_IID firstFrag) {
    double mean = 0;
    double sd = 0;
    if (thread_globals.globalRepeats == TRUE) {
        // compute global coverage of the reads we're correcting
        // the number of mappings each high-identity read has should be equal to the coverage of our data we're correcting, except for repeats
        // identify the cut off the peak
        // we do this by finding the inflection point at the end of the curve (that is where we are past the peak and no longer decreasing)
        double N = 0;
        double variance = 0;
        uint32* covHist = new uint32[MAX_COV_HIST + 1];
        memset(covHist, 0, (MAX_COV_HIST + 1) * sizeof(uint32));
        double prevScore = 0;

        for (map<AS_IID, uint8>::const_iterator iter = thread_globals.readsToPrint.begin(); iter != thread_globals.readsToPrint.end(); iter++) {
            if (iter->second != 0) {
                uint32 val = iter->second;
                if (val == MAX_COV) {
                    val = thread_globals.longReadsToPrint[iter->first];
                }
                N++;
                double delta = val - mean;
                mean += delta / N;
                variance += delta * (val - mean);

                covHist[MIN(MAX_COV_HIST, val)]++;
            }
        }
        variance /= N;
        sd = sqrt(variance);

        double prevRatio = 0;
        double runningTotal = covHist[0] + covHist[1];
        for (uint16 iter = 2; iter < MAX_COV_HIST; iter++) {
            if (covHist[iter] <= covHist[iter-1] && (double)covHist[iter]/covHist[iter-1] > prevRatio && (runningTotal / N > CUMULATIVE_SUM + ONE_SD_PERCENT)) {
                thread_globals.covCutoff = MAX(iter - 1, thread_globals.covCutoff);
                break;
            }
            prevRatio = (double)covHist[iter]/covHist[iter-1];
            runningTotal += covHist[iter];
        }
        if (thread_globals.covCutoff == 0) thread_globals.covCutoff = MAX_COV;

        // output histogram we used to compute this value
        char outputName[FILENAME_MAX] = {0};
        sprintf(outputName, "%s.layout.hist", thread_globals.prefix);
        errno = 0;
        FILE *histF =  fopen(outputName, "w");
        if (errno) {
            fprintf(stderr, "Couldn't open '%s' for write: %s\n", outputName, strerror(errno));
        } else {
            for (uint16 iter = 1; iter < MAX_COV_HIST; iter++) {
                fprintf(histF, "%d\t%d\n", iter, covHist[iter]);
            }
            fclose(histF);
        }
        delete[] covHist;
        if (thread_globals.verboseLevel >= VERBOSE_OFF) fprintf(stderr, "Picking cutoff as %d mean would be %f +- %f (%d)\n", thread_globals.covCutoff, mean, sd, (int)(ceil(mean + sd)));
    }
}

void filterRepeatReads(PBRThreadGlobals &thread_globals, AS_IID firstFrag) {
    fprintf(stderr, "Storing "F_SIZE_T" fragments out of "F_SIZE_T"\n", thread_globals.gappedReadSet->count(), thread_globals.gappedReadSet->size());

    // create files to partition the scores
    FILE **partitionedScores = new FILE*[thread_globals.partitions];
    ShortMapStore **partitionedStores = new ShortMapStore*[thread_globals.partitions];
    char outputName[FILENAME_MAX] = {0};
    for (int i = 0; i < thread_globals.partitions; i++) {
        sprintf(outputName, "%s.%d.rank", thread_globals.prefix, i+1);
        errno = 0;
        if (thread_globals.verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "Trying to open file %d named %s holding range %d-%d\n", i+1, outputName, thread_globals.partitionStarts[i].first, thread_globals.partitionStarts[i].second);
        partitionedScores[i] = fopen(outputName, "w");
        if (errno) {
            fprintf(stderr, "Couldn't open '%s' for write: %s\n", outputName, strerror(errno)); exit(1);
        }

        // create store as well
        sprintf(outputName, "%s.%d", thread_globals.prefix, i+1);

        if (thread_globals.maxUncorrectedGap > 0) {
            partitionedStores[i] = new ShortMapStore(outputName, true, true);
        }
    }

    if (thread_globals.globalRepeats == TRUE) {
        // now that we have a cutoff, stream the store and record which pacbio sequences the high-identity sequences should correct
        OverlapStore *ovs = AS_OVS_openOverlapStore(thread_globals.ovlStoreUniqPath);
        uint64 olapCount = 0;
        uint64 ovlPosition = 0;
        OVSoverlap *olaps = NULL;

        ShortMapRecord record(firstFrag, thread_globals.gkp->gkStore_getNumFragments(), thread_globals.covCutoff);

        for (map<uint32, uint8>::const_iterator iter = thread_globals.readsToPrint.begin(); iter != thread_globals.readsToPrint.end(); iter++) {
            // this read was not used, ignore it
            if (iter->second == 0) {
                continue;
            }

            // skip to the read we are interested in, the overlaps are sequential but we may have skipped some when we did the mapping
            while (olaps != NULL && ovlPosition < olapCount && olaps[ovlPosition].a_iid < iter->first) {
                ovlPosition++;
            }

            // again, read the store in batches
            if (ovlPosition >= olapCount) {
                delete[] olaps;

                AS_OVS_setRangeOverlapStore(ovs, iter->first,thread_globals.readsToPrint.rbegin()->first);
                olapCount = MIN(thread_globals.numThreads*MAX_TO_READ, AS_OVS_numOverlapsInRange(ovs));
                olaps = new OVSoverlap[olapCount];
                uint64 read = 0;
                uint64 last = olapCount;
                while (read < olapCount && last > 0) {
                    if (AS_OVS_readOverlapsFromStore(ovs, NULL, 0, AS_OVS_TYPE_ANY) <= olapCount - read) {
                        last = AS_OVS_readOverlapsFromStore(ovs, olaps+read, olapCount-read, AS_OVS_TYPE_ANY);
                        read+= last;
                    } else {
                        break;
                    }
                }
                olapCount = read;
                ovlPosition = 0;

                if (thread_globals.verboseLevel >= VERBOSE_OFF) fprintf(stderr, "Loaded "F_U64" overlaps\n", olapCount);
            }

            // store information on the short-read sequence
            SeqInterval bclr;
            record.readIID = iter->first;

            // build a sorted by score map of all mapping a short-read sequence has
            multimap<uint64, OverlapPos> scoreToReads;
            map<AS_IID, uint64> readsToScore;
            map<AS_IID, uint32> readsToBgn;

            uint32 alen = thread_globals.frgToLen[iter->first];
            for (uint64 rank = ovlPosition; rank < olapCount; rank++, ovlPosition++) {
                assert(ovlPosition == rank);
                if (olaps[ovlPosition].a_iid > iter->first) {
                    break;
                }
                uint32 blen = thread_globals.frgToLen[olaps[rank].b_iid];
                if (isOlapBad(olaps[ovlPosition], alen, blen, thread_globals.erate, thread_globals.elimit, thread_globals.maxErate)) {
                    continue;
                }
                if (olaps[ovlPosition].dat.ovl.type == AS_OVS_TYPE_OVL && (olaps[ovlPosition].dat.ovl.a_hang > 0 || olaps[ovlPosition].dat.ovl.b_hang < 0)) {
                    continue;
                }
                uint64 currScore = scoreOverlap(olaps[rank], alen, blen, thread_globals.erate, thread_globals.elimit, thread_globals.maxErate);
                OverlapPos pos;
                pos.ident = olaps[rank].b_iid;
                convertOverlapToPosition(olaps[rank], pos.position, bclr, alen, blen, true);

                // check for duplicate overlaps
                // TODO: this currently assumes ties are non-randomly broken with greater positions winning which is not ideal
                map<AS_IID, uint64>::iterator ovlIter = readsToScore.find(pos.ident);
                if (ovlIter != readsToScore.end()) {
                    if (currScore == ovlIter->second && readsToBgn[pos.ident] >= pos.position.bgn) {
                        continue;
                    }
                }
                readsToScore[pos.ident] = currScore;
                readsToBgn[pos.ident] = pos.position.bgn;
                scoreToReads.insert(pair<uint64, OverlapPos>(currScore, pos));
            }

            // now keep only the best cutoff of those and store them for easy access
            uint64 lastScore = 0;
            uint16 position = 0;
            bool storeShortMap = thread_globals.gappedReadSet->test(getBitSetID(record.readIID, &thread_globals));;
            set<AS_IID> readRanking;
            for (multimap<uint64, OverlapPos>::reverse_iterator rank = scoreToReads.rbegin(); rank != scoreToReads.rend(); rank++) {
                if (readRanking.find(rank->second.ident) != readRanking.end()) {
                    continue;  // only record each pac bio read at most once
                }
                position = (position == MAX_COV_HIST ? MAX_COV_HIST : (lastScore == rank->first ? position+1 : position+1));
                if (position > thread_globals.covCutoff) {
                    break;
                }
                readRanking.insert(rank->second.ident);
                lastScore = rank->first;

                // record this read for the short read
                record.addMapping(rank->second);
            }

            // output the score to the appropriate partition
            // at the same time, we store the info on the short read
            for (set<AS_IID>::iterator j = readRanking.begin(); j != readRanking.end(); j++) {
                // guess where we think the PacBio sequence ended up based on the size of the partitions
                // use two values due to rounding errors
                uint32 mpLow = MIN(thread_globals.partitions-1, MAX(0, (uint32) ceil((double)((*j)-firstFrag+1) / thread_globals.perFile)-1));
                uint32 mpHigh= MIN(thread_globals.partitions-1, MAX(0, (uint32) floor((double)((*j)-firstFrag+1) / thread_globals.perFile)-1));

                // make sure the smaller one is always called mpLow. This can be reversed when #fragments = #partitions
                if (mpLow == mpHigh && mpLow > 0) { mpLow--; }
                uint32 min = MIN(thread_globals.partitionStarts[mpLow].first, thread_globals.partitionStarts[mpHigh].first);
                uint32 max = MAX(thread_globals.partitionStarts[mpLow].second, thread_globals.partitionStarts[mpHigh].second);
                uint32 tmp = MIN(mpLow, mpHigh);
                mpHigh = MAX(mpLow, mpHigh);
                mpLow = tmp;

                // if we should go in the low partition, write it here
                if (thread_globals.partitionStarts[mpLow].first <= (*j) && thread_globals.partitionStarts[mpLow].second > (*j)) {
                    if (storeShortMap && thread_globals.maxUncorrectedGap > 0) {
                        partitionedStores[mpLow]->appendRecord(record);
                    }
                    fprintf(partitionedScores[mpLow], F_IID"\t"F_IID"\n", iter->first, (*j));
                    // if our second guess was right, write it here
                } else if (thread_globals.partitionStarts[mpHigh].first <= (*j) && thread_globals.partitionStarts[mpHigh].second > (*j)) {
                    if (storeShortMap && thread_globals.maxUncorrectedGap > 0) {
                        partitionedStores[mpHigh]->appendRecord(record);
                    }
                    fprintf(partitionedScores[mpHigh], F_IID"\t"F_IID"\n", iter->first, (*j));
                    // last ditch effort, our guesses were both wrong, search for the partition containing the pacbio sequence we need
                } else {
                    bool error = true;
                    if ((*j) < min) {
                        while (mpLow >= 0 && thread_globals.partitionStarts[mpLow].first > (*j)) {
                            mpLow--;
                        }
                        if (thread_globals.partitionStarts[mpLow].first < (*j) && thread_globals.partitionStarts[mpLow].second > (*j)) {
                            if (storeShortMap && thread_globals.maxUncorrectedGap > 0) {
                                partitionedStores[mpLow]->appendRecord(record);
                            }
                            fprintf(partitionedScores[mpLow], F_IID"\t"F_IID"\n", iter->first, (*j));

                            error = false;
                        }
                    } else {
                        while (mpHigh < thread_globals.partitions && thread_globals.partitionStarts[mpHigh].second <= (*j)) {
                            mpHigh++;
                        }
                        if (thread_globals.partitionStarts[mpHigh].first <= (*j) && thread_globals.partitionStarts[mpHigh].second > (*j)) {
                            if (storeShortMap && thread_globals.maxUncorrectedGap > 0) {
                                partitionedStores[mpHigh]->appendRecord(record);
                            }
                            fprintf(partitionedScores[mpHigh], F_IID"\t"F_IID"\n", iter->first, (*j));
                            error = false;
                        }
                    }
                    if (error) {
                        fprintf(stderr, "ERROR: Could not find appropriate partition for %d though it was either %d or %d (%d-%d) and (%d-%d) but it was neither\n", (*j), mpLow, mpHigh, thread_globals.partitionStarts[mpLow].first, thread_globals.partitionStarts[mpLow].second, thread_globals.partitionStarts[mpHigh].first, thread_globals.partitionStarts[mpHigh].second); exit(1);
                    }
                }
            }
            record.clearMappings();
        }
        delete[] olaps;
        AS_OVS_closeOverlapStore(ovs);
    }
    for (int i = 0; i < thread_globals.partitions; i++) {
        fclose(partitionedScores[i]);

        if (thread_globals.maxUncorrectedGap > 0) {
            delete partitionedStores[i];
        }
    }
    delete[] partitionedScores;
    delete[] partitionedStores;
}
