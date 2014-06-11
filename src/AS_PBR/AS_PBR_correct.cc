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

#include "AS_PBR_correct.hh"
#include "AS_PBR_util.hh"
#include "AS_PBR_store.hh"

#include <vector>
#include <algorithm>

static const char *rcsid_AS_PBR_CORRECT_C = "$Id$";

map<AS_IID, uint64> *globalFrgToScore;

/*
 - we track for each PacBio fragment if each bid maps more than once in a different orientation.
   - afterwards, if a large percentage of fragments does have multiple alignments to the same PacBio read in expected orientation flip, mark this PacBio sequence as containing a chimera
   - try to pick a point to split the sequence
      The fragments look like:
         17795  A:  230  338 (1077)  B:  100    0 ( 100)  15.00%                              <---------
	   	   	  17795  A:  528  629 (1077)  B:    0  100 ( 100)  21.00%                                                          -------->
	   	   	  18258  A:  306  406 (1077)  B:    0  100 ( 100)  16.00%                                     -------->
	   	   	  18258  A:  452  560 (1077)  B:  100    0 ( 100)  10.00%                                                  <---------

	   	   	  We want to pick the max position of the lower sequence and min position of the greater to pick our breakpoints
 */
// WARNING: for now this only handles one chimera per fragment. If there are more than one this code will only break at the first, it should be recursive
static SeqInterval findChimera(AS_IID aid, uint32 alen, const map<AS_IID, SeqInterval> &fwdMatches, const map<AS_IID, SeqInterval> &revMatches, uint8 verbosity) {
    uint32 totalInBothOri = 0;
    SeqInterval chimeraJunction;
    chimeraJunction.bgn = 0;
    chimeraJunction.end = alen;
    for (map<AS_IID, SeqInterval>::const_iterator iter = fwdMatches.begin(); iter != fwdMatches.end(); iter++) {
        map<AS_IID, SeqInterval>::const_iterator revMatch = revMatches.find(iter->first);
        if (revMatch != revMatches.end()) {
            totalInBothOri++;

            // rev read is first, take both begin sequences
            if (MIN(revMatch->second.bgn, revMatch->second.end) < MIN(iter->second.bgn, iter->second.end)) {
                chimeraJunction.bgn = MAX(chimeraJunction.bgn, revMatch->second.bgn);
                chimeraJunction.end = MIN(chimeraJunction.end, iter->second.bgn);
            } else{
                // fwd read is first, take both end sequences
                chimeraJunction.bgn = MAX(chimeraJunction.bgn, iter->second.end);
                chimeraJunction.end = MIN(chimeraJunction.end, revMatch->second.end);
            }
        }
    }
    if (chimeraJunction.bgn > chimeraJunction.end) {
        uint32 tmp = chimeraJunction.bgn;
        chimeraJunction.bgn = chimeraJunction.end;
        chimeraJunction.end = tmp;
    }
    if ((totalInBothOri > 1 && (chimeraJunction.end - chimeraJunction.bgn) < CHIMERA_MAX_SIZE) || ((double)totalInBothOri / MIN(fwdMatches.size(), revMatches.size()) > 0.5)) {
        if (verbosity >= VERBOSE_DEBUG) fprintf(stderr, "The pac read %d has a chimera with %d reads mapping both fwd/rev (out of %d and %d) at positions %d-%d\n", aid, totalInBothOri, fwdMatches.size(), revMatches.size(), chimeraJunction.bgn, chimeraJunction.end);
    } else {
        chimeraJunction.bgn = chimeraJunction.end = 0;
    }

    return chimeraJunction;
}

bool compare_by_identity (const OverlapPos &a, const OverlapPos &b) {
    if (globalFrgToScore == NULL || (*globalFrgToScore)[a.ident] == 0 || (*globalFrgToScore)[b.ident] == 0) {
        fprintf(stderr, "Error: fragment %d or %d does not have a defined score\n", a.ident, b.ident);
        assert(0);
    }

    return (*globalFrgToScore)[a.ident] < (*globalFrgToScore)[b.ident];
}


/**
 * This function is responsible for streaming through a subset of PacBio fragments and generating a tiling for them
 */
void *  correctFragments(void *ptr) {
    // get the global and thread-specific settings. The thread-specific settings tells us what data we are responsible for
    PBRThreadWorkArea *wa = (PBRThreadWorkArea *) (ptr);
    PBRThreadGlobals* waGlobal = wa->globals;
    uint32 counter = 0;

    fprintf(stderr, "Starting in thread %d from %d to %d\n", wa->id, wa->start, wa->end);

    OverlapStore *ovs = AS_OVS_openOverlapStore(waGlobal->ovlStoreUniqPath);

    // local copies of global variables
    map<AS_IID, uint8> readsToPrint;
    map<AS_IID, uint32> longReadsToPrint;

    // we keep a count of the mean number of mappings per sequence so we can identify tentative repetitive regions for patching
    double readsMean = 0;
    int readsCount = 0;

    boost::dynamic_bitset<> *gappedReadSet = initGappedReadSet(waGlobal);
    boost::dynamic_bitset<> *workingReadSet = initGappedReadSet(waGlobal);

    uint16  *readCoverage = new uint16[AS_READ_MAX_NORMAL_LEN];
    map<AS_IID, OVSoverlap> frgToBest;
    map<AS_IID, uint64> frgToScore;
    OVSoverlap olap;
    uint64 olapCount = 0;
    uint64 ovlPosition = 0;
    OVSoverlap *olaps = NULL;

    // get the number of files we should be generating
    char outputName[FILENAME_MAX] = {0};
    uint32 fileStart = wa->fileStart;
    uint32 lastFile = wa->fileEnd;
    uint32 numFiles = lastFile - fileStart + 1;
    uint32 perFile = waGlobal->perFile;

    // create the requested output partition
    pair<AS_IID, AS_IID> partitionStartEnd;
    uint32 currOpenID = fileStart;
    sprintf(outputName, "%s.%d.olaps", waGlobal->prefix, fileStart);
    errno = 0;
    LayRecordStore *outFile = createLayFile(outputName);
    fprintf(stderr, "In thread %d going to output files %d-%d with %d\n", wa->id, fileStart, fileStart + numFiles - 1, perFile);
    partitionStartEnd.first = wa->start;

    // finally stream through the fragments we're correcting and map the best overlaps to them
    for (AS_IID i = wa->start; i <= wa->end; i++) {
        if (waGlobal->libToInclude[waGlobal->frgToLib[i]] != TRUE) {
            continue;
        }

        if (counter % 10000 == 0) {
            if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "Thread %d done with %d fragments\n", wa->id, counter);
        }

        uint32 alen = waGlobal->frgToLen[i];
        memset(readCoverage, 0, alen * sizeof(uint16));

        frgToBest.clear();
        frgToScore.clear();

        map<uint32, OverlapPos> tile;
        LayRecord layout;
        layout.iid = i;

        // always output self match
        if (waGlobal->allowLong == TRUE) {
            SeqInterval aClr;
            aClr.bgn = 0;
            aClr.end = alen;

            SeqInterval pos;
            pos.bgn = aClr.bgn;
            pos.end = aClr.end;

            pos.bgn += TRIM_BACK_AMOUNT;
            pos.end -= TRIM_BACK_AMOUNT;
            aClr.bgn += TRIM_BACK_AMOUNT;
            aClr.end -= TRIM_BACK_AMOUNT;

            OVSoverlap olap;
            olap.a_iid = olap.b_iid = i;
            olap.dat.ovl.type = AS_OVS_TYPE_OVL;
            olap.dat.ovl.flipped = 0;
            olap.dat.ovl.corr_erate = olap.dat.ovl.orig_erate = AS_OVS_encodeQuality(0);
            olap.dat.ovl.a_hang = olap.dat.ovl.b_hang = 0;

            layout.bClrs.insert(pair<AS_IID, SeqInterval>(i, aClr));
            layout.bOvls.insert(pair<AS_IID, OVSoverlap>(i, olap));
            OverlapPos tileStr;
            tileStr.position = pos;
            tileStr.ident = i;
            tile[i] = tileStr;
        }

        map<AS_IID, uint8> toSkip;

        // read next batch of records from the overlap store, since this requires locking, we try to do it in batches
        if (ovlPosition >= olapCount) {
            delete[] olaps;
            if (waGlobal->numThreads > 1) {
                pthread_mutex_lock (&waGlobal->overlapMutex);
            }
            AS_OVS_setRangeOverlapStore(ovs, i, wa->end);
            olapCount = MIN(MAX_TO_READ, AS_OVS_numOverlapsInRange(ovs));
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
            if (waGlobal->numThreads > 1) {
                pthread_mutex_unlock(&waGlobal->overlapMutex);
            }
            olapCount = read;
            ovlPosition = 0;
            if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "Thread %d loaded "F_U64" overlaps\n", wa->id, olapCount);
        }

        // track fwd/reverse matches for a PacBio sequence
        map<AS_IID, SeqInterval> fwdMatches;
        map<AS_IID, SeqInterval> revMatches;

        // go through our block of loaded overlaps until we hit the end of this PacBio sequence (the a sequence)
        for (; ovlPosition < olapCount; ovlPosition++) {
            if (olaps[ovlPosition].a_iid != i) {
                break;
            }
            olap = olaps[ovlPosition];

            AS_IID aid = olap.a_iid;
            AS_IID bid = olap.b_iid;
            uint32 blen = waGlobal->frgToLen[bid];
            uint64 score = scoreOverlap(olap, alen, blen, waGlobal->erate, waGlobal->elimit, waGlobal->maxErate);

            // record the reads whether they are fwd or reverse to detect chimeric fragments
            SeqInterval pos;
            pos.bgn = pos.end = 0;
            SeqInterval bClr;
            bClr.bgn = bClr.end = 0;

            if (waGlobal->allowLong != FALSE) {
               allowCloseDovetail(olap);
            }
            convertOverlapToPosition(olap, pos, bClr, alen, blen);

            // when we have an overlap that is coming from the long reads, trim it back
            if (waGlobal->libToInclude[waGlobal->frgToLib[bid]] == TRUE) {
                // dont use library itself to correct fragments unless requested
                if (waGlobal->allowLong == FALSE) {
                    continue;
                }
                // skip self match, we always output it
                if (bid == aid) { continue; }
                int32 min = MAX(0,MIN(pos.bgn, pos.end));
                int32 max = MIN(alen, MAX(pos.bgn, pos.end));
                // bad overlap if it is too small to trim
                if (min + TRIM_BACK_AMOUNT > alen || max < TRIM_BACK_AMOUNT || bClr.bgn + TRIM_BACK_AMOUNT > blen || bClr.end < TRIM_BACK_AMOUNT) {
                 if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Fragment %d (%d) (%d-%d) versus %d (%d) (%d-%d) with error rate %f is short 1 , skipping it\n", aid, alen, min, max, bid, blen, bClr.bgn, bClr.end, AS_OVS_decodeQuality(olap.dat.ovl.corr_erate));
                    continue;
                }
                if (min + TRIM_BACK_AMOUNT < 0 || max - TRIM_BACK_AMOUNT > alen) {
                 if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Fragment %d (%d) (%d-%d) versus %d (%d) (%d-%d) with error rate %f is short 2, skipping it\n", aid, alen, min, max, bid, blen, bClr.bgn, bClr.end, AS_OVS_decodeQuality(olap.dat.ovl.corr_erate));
                    continue;
                }
                if (pos.bgn < pos.end) {
                    pos.bgn += TRIM_BACK_AMOUNT;
                    pos.end -= TRIM_BACK_AMOUNT;
                } else {
                    pos.bgn -= TRIM_BACK_AMOUNT;
                    pos.end += TRIM_BACK_AMOUNT;
                }
                bClr.bgn += TRIM_BACK_AMOUNT;
                bClr.end -= TRIM_BACK_AMOUNT;
            }

            // should this only use contained overlaps or not?
            if (isOvlForward(olap)) {
                fwdMatches[bid] = pos;
            } else {
                revMatches[bid] = pos;
            }

            // skip overlaps below our quality criteria
            if (isOlapBad(olap, alen, blen, waGlobal->erate, waGlobal->elimit, waGlobal->maxErate)) {
                 if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Fragment %d (%d) versus %d (%d) with error rate %f is bad, skipping it\n", aid, alen, bid, blen, AS_OVS_decodeQuality(olap.dat.ovl.corr_erate));
                continue;
            }

            // non contained overlap, dont use these fragments for correction, if close then ok
            if (olap.dat.ovl.type == AS_OVS_TYPE_OVL && (olap.dat.ovl.a_hang < 0 || olap.dat.ovl.b_hang > 0)) {
                // if there is a better non-contained overlap, then do not trust the contained overlap for this same fragment
                if (frgToScore[bid] != 0) {
if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Fragment %d being kicked out of %d because it has a better non-contained overlap\n", bid, aid);
                    OVSoverlap best = frgToBest[bid];
                    uint32 min = MIN(tile[best.b_iid].position.bgn, tile[best.b_iid].position.end);
                    uint32 max = MAX(tile[best.b_iid].position.bgn, tile[best.b_iid].position.end);
                    for (uint32 iter = min; iter <= max && waGlobal->globalRepeats == FALSE; iter++) {
                        readCoverage[iter]--;
                    }
                    if (longReadsToPrint[best.b_iid] != 0) {
                        readsMean -= longReadsToPrint[best.b_iid];
                        longReadsToPrint[best.b_iid]--;
                        if (longReadsToPrint[best.b_iid] < MAX_COV) {
                            readsToPrint[best.b_iid] = longReadsToPrint[best.b_iid];
                            longReadsToPrint[best.b_iid] = 0;

                            readsMean += readsToPrint[best.b_iid];
                        } else {
                            readsMean += longReadsToPrint[best.b_iid];
                        }
                    } else {
                        readsMean -= readsToPrint[best.b_iid];
                        readsToPrint[best.b_iid]--;
                        readsMean += readsToPrint[best.b_iid];
                    }
                    tile[best.b_iid].position.bgn = tile[best.b_iid].position.end = 0;
                }
                toSkip[olap.b_iid] = 1;
                continue;
            }

            if (toSkip[bid] != 0) {
                continue;
            }

            // if this overlap is worse than one we already found, skip it
            if (score < frgToScore[bid]) {
                continue;
            }

            // remove a previous match if this one is better
            if (frgToScore[bid] <= score && frgToScore[bid] != 0) {
                OVSoverlap best = frgToBest[bid];
                uint32 min = MIN(tile[best.b_iid].position.bgn, tile[best.b_iid].position.end);
                uint32 max = MAX(tile[best.b_iid].position.bgn, tile[best.b_iid].position.end);
                for (uint32 iter = min; iter <= max && waGlobal->globalRepeats == FALSE; iter++) {
                    readCoverage[iter]--;
                }
                if (longReadsToPrint[best.b_iid] != 0) {
                    readsMean -= longReadsToPrint[best.b_iid];
                    longReadsToPrint[best.b_iid]--;
                    if (longReadsToPrint[best.b_iid] < MAX_COV) {
                        readsToPrint[best.b_iid] = longReadsToPrint[best.b_iid];
                        longReadsToPrint[best.b_iid] = MAX_COV;
                        readsMean += readsToPrint[best.b_iid];
                    } else {
                        readsMean += longReadsToPrint[best.b_iid];
                    }
                } else {
                    readsMean -= readsToPrint[best.b_iid];
                    readsToPrint[best.b_iid]--;
                    readsMean += readsToPrint[best.b_iid];
                }
                tile[best.b_iid].position.bgn = tile[best.b_iid].position.end = 0;
            }
            frgToBest[bid] = olap;
            frgToScore[bid] = score;

            layout.bClrs.insert(pair<AS_IID, SeqInterval>(bid, bClr));
            map<AS_IID, OVSoverlap>::iterator ovlIter = layout.bOvls.find(bid);
            if (ovlIter == layout.bOvls.end()) {
               layout.bOvls.insert(pair<AS_IID, OVSoverlap>(bid, olap));
            } else {
               ovlIter->second = olap;
            }

            OverlapPos tileStr;
            tileStr.position = pos;
            tileStr.ident = bid;
            tile[bid] = tileStr;

            if (readsToPrint.find(bid) == readsToPrint.end()) {
                readsToPrint[bid] = 0;
                readsCount++;
            }
            if (readsToPrint[bid] == MAX_COV) {
                if (longReadsToPrint.find(bid) == longReadsToPrint.end()) { longReadsToPrint[bid] = MAX_COV; }
                readsMean -= longReadsToPrint[bid];
                longReadsToPrint[bid]++;
                readsMean += longReadsToPrint[bid];
            } else {
                readsMean -= readsToPrint[bid];
                readsToPrint[bid]++;
                if (readsToPrint[bid] == MAX_COV) { longReadsToPrint[bid] = MAX_COV; }
                readsMean += readsToPrint[bid];
            }

            uint32 min = MIN(pos.bgn, pos.end);
            uint32 max = MAX(pos.bgn, pos.end);
            for (uint32 iter = min; iter <= max && waGlobal->globalRepeats == FALSE; iter++) {
                readCoverage[iter] = (readCoverage[iter] == MAX_COV ? MAX_COV : readCoverage[iter]+1);
            }
        }

        // check for a chimera and update this sequence if it is chimeric
        SeqInterval chimeraJunction = findChimera(i, alen, fwdMatches, revMatches, waGlobal->verboseLevel);
        fwdMatches.clear();
        revMatches.clear();

        // update our tiling based on the detected chimera
        for (map<uint32, OverlapPos>::const_iterator iter = tile.begin(); iter != tile.end(); iter++) {
            if (iter->second.position.bgn > 0 || iter->second.position.end > 0) {
                bool toAdd = true;

                // pick only the largest chunk of non-chimeric fragment to print
                if (chimeraJunction.bgn != 0 && chimeraJunction.end != 0) {
                    if (rangesOverlap(chimeraJunction, iter->second.position)) {
                        toAdd = false;
                        //fprintf(stderr, "Skipping read %d in %d because it intersects chimera %d ( %d %d )\n", iter->first, i, chimeraJunction.bgn, iter->second.position.bgn, iter->second.position.end);
                    }
                    else if (alen - chimeraJunction.end > chimeraJunction.bgn) {
                        if (MAX(iter->second.position.bgn, iter->second.position.end) < chimeraJunction.bgn) {
                            toAdd = false;
                            //fprintf(stderr, "Skipping read %d in %d because it is on the wrong (start) size of chimera %d ( %d %d )\n", iter->first, i, chimeraJunction.bgn, iter->second.position.bgn, iter->second.position.end);
                        }
                    } else {
                        if (MIN(iter->second.position.bgn, iter->second.position.end) > chimeraJunction.end) {
                            toAdd = false;
                            //fprintf(stderr, "Skipping read %d in %d because it is on the wrong (long) size of chimera %d ( %d %d )\n", iter->first, i, chimeraJunction.end, iter->second.position.bgn, iter->second.position.end);
                        }
                    }
                }

                if (toAdd == true) {
                    layout.mp.push_back(iter->second);
                }
            }
        }

        if (layout.mp.size() > 0) {
            // check for repeats within this sequence, if we were asked
            if (waGlobal->globalRepeats == FALSE) {
                double mean = 0;
                double N = 0;

                for (uint32 iter = 0; iter < alen; iter++) {
                    if (readCoverage[iter] > 0) {
                        N++;
                        double delta = readCoverage[iter] - mean;
                        mean += delta / N;
                    }
                }
                mean = MAX(2, mean);

                memset(readCoverage, 0, alen * sizeof(uint16));

                if (waGlobal->numThreads > 1) {
                    pthread_mutex_lock( &waGlobal->globalDataMutex);
                }
                globalFrgToScore = &frgToScore;
                stable_sort(layout.mp.begin(), layout.mp.end(), compare_by_identity);
                globalFrgToScore = NULL;
                if (waGlobal->numThreads > 1) {
                    pthread_mutex_unlock( &waGlobal->globalDataMutex);
                }

                for (vector<OverlapPos>::iterator iter = layout.mp.begin(); iter != layout.mp.end(); ) {
                    uint32 min = MIN(iter->position.bgn, iter->position.end);
                    uint32 max = MAX(iter->position.bgn, iter->position.end);
                    uint16  max_cov = 0;

                    for (uint32 coverage = min; coverage <= max; coverage++) {
                        if (max_cov < readCoverage[coverage]) {
                            max_cov = readCoverage[coverage];
                        }
                        readCoverage[coverage] = (readCoverage[coverage] == MAX_COV ? MAX_COV : readCoverage[coverage]+1);
                    }
                    if (max_cov > (mean * waGlobal->repeatMultiplier)) {
                        iter = layout.mp.erase(iter);
                    } else {
                        iter++;
                    }
                }
            }

            // sort the tiling by position
            stable_sort(layout.mp.begin(), layout.mp.end(), compare_tile);
            // now save the tiling, if we have filled up a file, move onto the next file before writing
            uint32 fileId = fileStart + MIN(numFiles - 1, floor((double)(i - wa->start) / perFile));
            if (fileId != currOpenID) {
                partitionStartEnd.second = i;
                if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Thread %d set partition %d to be %d-%d\n", wa->id, currOpenID-1, partitionStartEnd.first, partitionStartEnd.second);
                waGlobal->partitionStarts[currOpenID-1] = partitionStartEnd;
                partitionStartEnd.first = partitionStartEnd.second = i;
                closeLayFile(outFile);
                sprintf(outputName, "%s.%d.olaps", waGlobal->prefix, fileId);
                errno = 0;
                outFile = createLayFile(outputName);
                if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "For thread %d I'm currently on file %d and need to move to %d\n", wa->id, currOpenID, fileId);
                currOpenID = fileId;
            }
            writeLayRecord(outFile, layout, workingReadSet, waGlobal->percentShortReadsToStore, readsToPrint, longReadsToPrint, readsMean / readsCount);
            (*gappedReadSet) |= (*workingReadSet);
            workingReadSet->reset();
        }

        counter++;
    }
    partitionStartEnd.second = wa->end + 1;
    waGlobal->partitionStarts[currOpenID-1] = partitionStartEnd;
    if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Thread %d set partition %d to be %d-%d\n", wa->id, currOpenID-1, partitionStartEnd.first, partitionStartEnd.second);
    closeLayFile(outFile);

    // finally update the global data on the mappings by each short-read
    if (waGlobal->numThreads > 1) {
        pthread_mutex_lock( &waGlobal->globalDataMutex);
    }
    for (map<AS_IID, uint8>::iterator iter = readsToPrint.begin(); iter != readsToPrint.end(); iter++) {
        uint32 cov = 0;
        if (iter->second == MAX_COV) {
            assert(longReadsToPrint.find(iter->first) != longReadsToPrint.end());
            cov = longReadsToPrint.find(iter->first)->second;
        } else {
            cov = iter->second;
        }
        if (waGlobal->readsToPrint.find(iter->first) == waGlobal->readsToPrint.end()) { waGlobal->readsToPrint[iter->first] = 0; }
        if (((uint32)waGlobal->readsToPrint[iter->first] + cov) > MAX_COV) {
            waGlobal->readsToPrint[iter->first] = MAX_COV;
            if (waGlobal->longReadsToPrint.find(iter->first) == waGlobal->longReadsToPrint.end()) { waGlobal->longReadsToPrint[iter->first] = 0; }
            waGlobal->longReadsToPrint[iter->first] += cov;
        } else {
            waGlobal->readsToPrint[iter->first] += cov;
        }
    }
    if (waGlobal->hasMates == false) {
        (*waGlobal->gappedReadSet) |= (*gappedReadSet);
    }
    if (waGlobal->numThreads > 1) {
        pthread_mutex_unlock( &waGlobal->globalDataMutex);
    }

    if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "Thread shutting down after finishing %d fragments\n", counter);
    delete gappedReadSet;
    delete workingReadSet;
    delete[] readCoverage;
    delete[] olaps;
    readsToPrint.clear();
    longReadsToPrint.clear();
    AS_OVS_closeOverlapStore(ovs);

    return 0;
}
