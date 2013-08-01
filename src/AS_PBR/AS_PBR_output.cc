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

#include "AS_PBR_output.hh"
#include "AS_PBR_store.hh"

#include "AS_OVS_overlapStore.H"
#include "AS_UTL_reverseComplement.H"
#include "AS_PER_encodeSequenceQuality.H"

#include <sstream>
#include <map>
#include <vector>
#include <set>

static const char *rcsid_AS_PBR_OUTPUT_C = "$Id: AS_PBR_output.cc,v 1.9 2012-12-27 15:17:49 skoren Exp $";

static const uint32 FUDGE_BP = 5;

// search other pacbio sequences for shared short-reads and recruit their sequences to help our gaps
static void getCandidateOverlaps(PBRThreadGlobals *waGlobal, boost::dynamic_bitset<> &bits, map<AS_IID, SeqInterval> &matchingSequencePositions, map<AS_IID, set<AS_IID> > &readRanking, vector<OverlapPos>::const_iterator &iter, LayRecord &layRecord, ShortMapStore *inStore) {
    // get other candidate pacbio reads that overlap this one
    map<AS_IID, bool> matchingSequenceOrientation;
    map<AS_IID, SeqInterval> matchingSequenceLastFwd;
    map<AS_IID, SeqInterval> matchingSequenceLastRev;

    vector<OverlapPos>::const_iterator fwd = iter;
    vector<OverlapPos>::const_iterator rev = iter;

    bool fwdDistanceSatisfied = false;
    bool revDistanceSatisfied = false;
    uint32 revGapStart = 0;
    uint32 lastFwd = 0;
    uint32 lastRev = 0;

    while (!fwdDistanceSatisfied || !revDistanceSatisfied) {
        ShortMapRecord *record = NULL;
        // stop when we've eliminated all supporters
        if (bits.size() != 0 && bits.count() == 1) {
            if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Eliminated all candidates in sequence %d while pulling from %d\n", layRecord.iid, iter->ident);
            break;
        }

        // first we check the forward direction (in front of the gap) for those who can help us
        if (!fwdDistanceSatisfied && fwd != layRecord.mp.end()) {
            if (waGlobal->globalRepeats == TRUE && (readRanking[fwd->ident].find(layRecord.iid) == readRanking[fwd->ident].end())) {
                // ignore
            } else {
                if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Pulling sequence %d forward of %d in %d \n", fwd->ident, iter->ident, layRecord.iid);
                record = inStore->getRecord(fwd->ident);
                // skip sequences that dont advance us
                if (lastFwd != 0 && lastFwd < MIN(fwd->position.bgn, fwd->position.end)) {
                    fwdDistanceSatisfied = true;
                }
                else if ((fwd+1) != layRecord.mp.end() && MAX(fwd->position.bgn, fwd->position.end) <= lastFwd+FUDGE_BP) {
                    // skip
                    lastFwd = MAX(fwd->position.bgn, fwd->position.end);
                }
                else {
                    if (record != NULL) {
                        lastFwd = MAX(fwd->position.bgn, fwd->position.end);

                        if (bits.size() == 0) {
                            bits = (*record->getMappedReads());
                        } else {
                            bits &= (*record->getMappedReads());
                        }

                        // check mapped candidates to see if they are good
                        bool iterFwd = fwd->position.bgn < fwd->position.end;
                        for(int i = bits.find_first(); i != boost::dynamic_bitset<>::npos; i = bits.find_next(i)) {
                            AS_IID id = inStore->getMappedIID(i);
                            SeqInterval *search = record->getMapping(i);
                            if (search != NULL) {
                                bool isFwd = search->bgn < search->end;
                                uint32 min = MIN(search->bgn, search->end);
                                bool isContained = false;
                                if (matchingSequenceLastFwd.find(id) != matchingSequenceLastFwd.end()) {
                                    min = MIN(matchingSequenceLastFwd[id].bgn, matchingSequenceLastFwd[id].end);
                                    uint32 max = MAX(matchingSequenceLastFwd[id].bgn, matchingSequenceLastFwd[id].end);
                                    uint32 myMin = MIN(search->bgn, search->end);
                                    uint32 myMax = MAX(search->bgn, search->end);
                                    isContained =(myMin >= min && myMax <= max);
                                }

                                if (matchingSequenceOrientation.find(id) == matchingSequenceOrientation.end()) {
                                    matchingSequenceOrientation[id] = iterFwd == isFwd;
                                }
                                if (matchingSequenceOrientation[id] && isFwd == iterFwd && (MIN(search->bgn, search->end) >= MIN(0, min - FUDGE_BP) || (isContained))) {
                                    if (matchingSequencePositions.find(id) == matchingSequencePositions.end()) {
                                        matchingSequencePositions[id].bgn = 0;
                                        matchingSequencePositions[id].end = INT32_MAX;
                                    }
                                    if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "For candidate sequence %d set it to be forward is %d and current sequence is %d based on short-read %d and will check end for %d in %d\n", id, isFwd, iterFwd, fwd->ident, MIN(search->bgn, search->end), layRecord.iid);
                                    matchingSequencePositions[id].end = MIN(matchingSequencePositions[id].end, MIN(search->bgn, search->end));
                                } else if (!matchingSequenceOrientation[id] && isFwd != iterFwd && (MIN(search->bgn, search->end) <= (min + FUDGE_BP) || (isContained))) {
                                    if (matchingSequencePositions.find(id) == matchingSequencePositions.end()) {
                                        matchingSequencePositions[id].bgn = INT32_MAX;
                                        matchingSequencePositions[id].end = 0;
                                    }
                                    if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "For candidate sequence %d set it to be forward is %d and current sequence is %d based on short-read %d and will check end for %d in %d\n", id, isFwd, iterFwd, fwd->ident, MAX(search->bgn, search->end), layRecord.iid);
                                    matchingSequencePositions[id].end = MAX(matchingSequencePositions[id].end, MAX(search->bgn, search->end));
                                } else {
                                    bits.set(i, false);
                                    if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Eliminating %s candidate read %d with orientation %d vs %d (expected %d) min %d vs %d for gap in %d\n", (matchingSequenceOrientation[id] == true ? "fwd" : "rev"), id, isFwd, iterFwd, matchingSequenceOrientation[id], MIN(search->bgn, search->end), min, layRecord.iid);
                                }
                                matchingSequenceLastFwd[id] = *search;
                                matchingSequenceOrientation[id] = iterFwd == isFwd;
                            } else {
                                bits.set(i, false);
                                if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Eliminating %s candidate read %d for gap in %d it had no mapping for %d\n", (matchingSequenceOrientation[id] == true ? "fwd" : "rev"), id, layRecord.iid, record->readIID);
                            }
                        }
                    }
                }
                if (MAX(fwd->position.bgn, fwd->position.end) - MIN(iter->position.bgn, iter->position.end) > MIN_DIST_TO_RECRUIT) {
                    fwdDistanceSatisfied = true;
                }
                delete record;
            }
            // move forward to next good sequence
            if (fwd != layRecord.mp.end()) {
                fwd++;
                while (waGlobal->globalRepeats == TRUE && (fwd != layRecord.mp.end() && readRanking[fwd->ident].find(layRecord.iid) == readRanking[fwd->ident].end())) {
                    // if we can keep going, skip this sequence
                    fwd++;
                }
                if (fwd == layRecord.mp.end()) { fwdDistanceSatisfied = true; }
            }
        }

        // next, check the sequences before the gap to find helpers
        // TODO: this is almost exactly duplicating the code above, clean up
        if (!revDistanceSatisfied) {
            // move back to the next good sequence
            if (rev != layRecord.mp.begin()) {
                rev--;
                while (waGlobal->globalRepeats == TRUE && (readRanking[rev->ident].find(layRecord.iid) == readRanking[rev->ident].end())) {
                    if (rev == layRecord.mp.begin()) { revDistanceSatisfied = true; break;}
                    rev--;
                }
            }
            // if we stopped because we ran off the list
            if (revDistanceSatisfied == true) {
                continue;
            }

            if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Pulling sequence %d reverse of %d in %d \n", rev->ident, iter->ident, layRecord.iid);
            if (lastRev != 0 && lastRev > MAX(rev->position.bgn, rev->position.end)) {
                if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Hit a gap before sequence in %d\n", layRecord.iid);
                revDistanceSatisfied = true;
            }
            else if (rev != layRecord.mp.begin() && MIN(rev->position.bgn, rev->position.end) <= lastRev+FUDGE_BP) {
                if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Skipping because it didnt advance\n");
                lastRev = MIN(rev->position.bgn, rev->position.end);
            }
            else {
                record = inStore->getRecord(rev->ident);

                if (record != NULL) {
                    lastRev = MIN(rev->position.bgn, rev->position.end);

                    if (bits.size() == 0) {
                        bits = (*record->getMappedReads());
                    } else {
                        bits &= (*record->getMappedReads());
                    }

                    // check mapped candidates to see if they are good
                    bool iterFwd = rev->position.bgn < rev->position.end;
                    for(int i = bits.find_first(); i != boost::dynamic_bitset<>::npos; i = bits.find_next(i)) {
                        AS_IID id = inStore->getMappedIID(i);
                        SeqInterval *search = record->getMapping(i);
                        if (search != NULL) {
                            uint32 min = MIN(search->bgn, search->end);
                            bool isContained = false;
                            if (matchingSequenceLastRev.find(id) != matchingSequenceLastRev.end()) {
                                min = MIN(matchingSequenceLastRev[id].bgn, matchingSequenceLastRev[id].end);
                                uint32 max = MAX(matchingSequenceLastRev[id].bgn, matchingSequenceLastRev[id].end);
                                uint32 myMin = MIN(search->bgn, search->end);
                                uint32 myMax = MAX(search->bgn, search->end);
                                isContained = (myMin >= min && myMax <= max);
                            }
                            bool isFwd = search->bgn < search->end;
                            if (matchingSequenceOrientation.find(id) == matchingSequenceOrientation.end()) {
                                matchingSequenceOrientation[id] = iterFwd == isFwd;
                            }
                            if (matchingSequenceOrientation[id] && isFwd == iterFwd && (MIN(search->bgn, search->end) <= (min + FUDGE_BP) || (isContained))) {
                                if (matchingSequencePositions.find(id) == matchingSequencePositions.end()) {
                                    matchingSequencePositions[id].bgn = 0;
                                    matchingSequencePositions[id].end = INT32_MAX;
                                }
                                if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "For candidate sequence %d set it to be forward is %d and current sequence is %d based on short-read %d and will check bgn for %d in %d \n", id, isFwd, iterFwd, rev->ident, MAX(search->bgn, search->end), layRecord.iid);
                                matchingSequencePositions[id].bgn = MAX(matchingSequencePositions[id].bgn, MAX(search->bgn, search->end)-1);
                            } else if (!matchingSequenceOrientation[id] && isFwd != iterFwd && (MIN(search->bgn, search->end) >= MIN(0, min - FUDGE_BP) || (isContained))) {
                                if (matchingSequencePositions.find(id) == matchingSequencePositions.end()) {
                                    matchingSequencePositions[id].bgn = INT32_MAX;
                                    matchingSequencePositions[id].end = 0;
                                }
                                if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "For candidate sequence %d set it to be forward is %d and current sequence is %d based on short-read %d and will check bgn for %d in %d\n", id, isFwd, iterFwd, rev->ident, MIN(search->bgn, search->end), layRecord.iid);
                                matchingSequencePositions[id].bgn = MIN(matchingSequencePositions[id].bgn, MIN(search->bgn, search->end)-1);
                            } else {
                                if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Eliminating %s candidate read %d with orientation %d vs %d (expected %d) min %d vs %d for gap in %d\n", (matchingSequenceOrientation[id] == true ? "fwd" : "rev"), id, isFwd, iterFwd, matchingSequenceOrientation[id], MIN(search->bgn, search->end), min, layRecord.iid);
                                bits.set(i, false);
                            }
                            matchingSequenceLastRev[id] = *search;
                            matchingSequenceOrientation[id] = iterFwd == isFwd;
                        } else {
                            bits.set(i, false);
                            if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Eliminating %s candidate read %d for gap in %d it had no mapping for %d\n", (matchingSequenceOrientation[id] == true ? "fwd" : "rev"), id, layRecord.iid, record->readIID);
                        }
                    }
                }
                delete record;
            }

            if (revGapStart == 0) {revGapStart = MAX(rev->position.bgn, rev->position.end); }
            else {
                if (revGapStart - MIN(rev->position.bgn, rev->position.end) > MIN_DIST_TO_RECRUIT) {
                    if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "in read %d rev distance was satisfied by read %d at %d %d\n", layRecord.iid, rev->ident, rev->position.bgn, rev->position.end);
                    revDistanceSatisfied = true;
                }
            }
            if (rev == layRecord.mp.begin()) { revDistanceSatisfied = true; }
        }
    }

    // finally, go through our list of helper sequences and make sure they are properly initialized and have valid positions
    for(int i = bits.find_first(); i != boost::dynamic_bitset<>::npos; i = bits.find_next(i)) {
        AS_IID iid = inStore->getMappedIID(i);
        map<AS_IID, SeqInterval>::iterator iter = matchingSequencePositions.find(iid);

        if (matchingSequenceLastFwd.find(iid) == matchingSequenceLastFwd.end() ||
                matchingSequenceLastRev.find(iid) == matchingSequenceLastRev.end()) {
            if (waGlobal->verboseLevel >= VERBOSE_DEVELOPER) fprintf(stderr, "Uninitialized fragment %d supporter of %d on either fwd or rev end with positions %d %d\n", iter->first, layRecord.iid, iter->second.bgn, iter->second.end);
            bits.set(inStore->getStoreIID(iid), false);
        }

        int32 tmp;
        if (matchingSequenceOrientation[iid]) {
            tmp = MIN(iter->second.bgn, iter->second.end);
            iter->second.end = MAX(iter->second.bgn, iter->second.end);
            iter->second.bgn = tmp;
        } else {
            tmp = MIN(iter->second.bgn, iter->second.end);
            iter->second.bgn = MAX(iter->second.bgn, iter->second.end);
            iter->second.end = tmp;
        }
    }
}

/**
 * Output a single AMOS layout record
 */
static void closeRecord(PBRThreadGlobals *waGlobal, FILE *outFile, FILE *logFile, stringstream &layout, LayRecord &layRecord, uint32 lastEnd, int32 &offset, uint32 &readIID, uint32 &readSubID) {
    // close the layout and start a new one because we have a coverage gap
    if (lastEnd - offset >= waGlobal->minLength) {
        fprintf(outFile, "%s}\n", layout.str().c_str());
        fprintf(logFile, "%d\t%s_%d_%d\t%d\t%d\t%d\n", layRecord.iid, waGlobal->libName, layRecord.iid, readSubID, readSubID, offset, lastEnd);
    }
    readIID++;
    readSubID++;
    offset = -1;
    layout.str("");
    layout << "{LAY\neid:" << waGlobal->libName << "_" << layRecord.iid << "_" << readSubID << "\niid:" << readIID << "\n";
}

/**
 *output AMOS-style layout messaged based on our internal structures
 */
void outputResults(PBRThreadGlobals *waGlobal, uint32 part) {
        assert(part > 0);
        map<AS_IID, uint8> readsToPrint;
        map<AS_IID, uint8> readsWithGaps;
        map<AS_IID, vector<pair<AS_IID, pair<uint32, uint32> > > > gaps;
        map<AS_IID, set<AS_IID> > readRanking;
        map<AS_IID, char*> frgToEnc;

        char reportName[FILENAME_MAX] = {0};
        sprintf(reportName, "%s.%d.log", waGlobal->prefix, part);
        errno = 0;
        FILE *reportFile = fopen(reportName, "w");
        if (errno) {
            fprintf(stderr, "Couldn't open '%s' for write: %s\n", reportName, strerror(errno)); exit(1);
        }

        char outputName[FILENAME_MAX] = {0};
        sprintf(outputName, "%s.%d.lay", waGlobal->prefix, part);
        errno = 0;
        FILE *outFile = fopen(outputName, "w");
        if (errno) {
            fprintf(stderr, "Couldn't open '%s' for write: %s\n", outputName, strerror(errno)); exit(1);
        }

        char inName[FILENAME_MAX] = {0};
        if (waGlobal->hasMates) {
            sprintf(inName, "%s.%d.paired.olaps", waGlobal->prefix, part);
        } else {
            sprintf(inName, "%s.%d.olaps", waGlobal->prefix, part);
        }
        errno = 0;
        LayRecordStore *inFile = openLayFile(inName);
        if (errno) {
            fprintf(stderr, "Couldn't open '%s' for read %s\n", inName, strerror(errno));
            return;
        }

        char inRankName[FILENAME_MAX] = {0};

        ShortMapStore *inStore = NULL;
        if (waGlobal->maxUncorrectedGap > 0) {
            sprintf(inRankName, "%s.%d", waGlobal->prefix, part);
            inStore = new ShortMapStore(inRankName, false, false, true);
        }
        sprintf(inRankName, "%s.%d.rank", waGlobal->prefix, part);
        errno = 0;
        FILE *inRankFile = fopen(inRankName, "r");
        if (errno) {
            fprintf(stderr, "Couldn't open %s for read %s\n", inRankName, strerror(errno)); exit(1);
        }
        while (!feof(inRankFile)) {
            AS_IID illumina;
            AS_IID corrected;

            fscanf(inRankFile, F_IID"\t"F_IID"\n", &illumina, &corrected);
            readRanking[illumina].insert(corrected);
        }
        fclose(inRankFile);
        AS_UTL_unlink(inRankName);

        uint32 readIID = 0;
        char seq[AS_READ_MAX_NORMAL_LEN];
        char qlt[AS_READ_MAX_NORMAL_LEN];
        AS_IID gapIID = waGlobal->gkp->gkStore_getNumFragments()+1;

        LayRecord layRecord;
        while (readLayRecord(inFile, layRecord)) {
            uint32 readSubID = 1;

            stringstream layout (stringstream::in | stringstream::out);
            layout << "{LAY\neid:" << waGlobal->libName << "_" << layRecord.iid << "_" << readSubID << "\niid:" << readIID << "\n";
            uint32 lastEnd = 0;
            int32 offset = -1;

            // process record
            for (vector<OverlapPos>::const_iterator iter = layRecord.mp.begin(); iter != layRecord.mp.end(); iter++) {
                // skip reads over coverage
                if (waGlobal->globalRepeats == TRUE && (readRanking[iter->ident].find(layRecord.iid) == readRanking[iter->ident].end())) {
                    //fprintf(stderr, "Skipping read %d to correct %d it was at cutoff %d true %d\n", iter->ident, i, waGlobal->readRanking[iter->ident][i].first, waGlobal->readRanking[iter->ident][i].second);
                    continue;
                }
                // if the last fragment ended before the current one starts, we have a gap
                if (lastEnd != 0 && lastEnd <= MIN(iter->position.bgn, iter->position.end)) {
                    // if we were asked, instead of skipping, we will output the uncorrected sequence
                    if (MIN(iter->position.bgn, iter->position.end) - lastEnd < waGlobal->maxUncorrectedGap) {
                        boost::dynamic_bitset<> bits;
                        map<AS_IID, SeqInterval> matchingSequencePositions;
                        getCandidateOverlaps(waGlobal, bits, matchingSequencePositions, readRanking, iter, layRecord, inStore);
                        uint32 gapStart = lastEnd;
                        uint32 gapEnd = MIN(iter->position.bgn, iter->position.end);
                        if (bits.size() != 0 && bits.test(inStore->getStoreIID(layRecord.iid)) == true) {
                            gapStart = matchingSequencePositions[layRecord.iid].bgn;
                            gapEnd = matchingSequencePositions[layRecord.iid].end;
                        }
                        assert(gapEnd >= gapStart);

                        if (bits.count() > 1) { // we found some candidates, use them
                            uint32 count = 0;

                            for(int i = bits.find_first(); i != boost::dynamic_bitset<>::npos; i = bits.find_next(i)) {
                                AS_IID iid = inStore->getMappedIID(i);

                                if (iid != layRecord.iid) {
                                    uint32 min = MIN(matchingSequencePositions[iid].bgn,matchingSequencePositions[iid].end);
                                    uint32 max = MAX(matchingSequencePositions[iid].bgn,matchingSequencePositions[iid].end);
                                    uint32 gapSize = (int32)(MIN(iter->position.bgn, iter->position.end) - gapStart + 1);
                                    uint32 diff = abs((int32)(max-min+1) - (int32)gapSize);
                                    if ((double)diff / gapSize <= (ERATE_ADJUST*waGlobal->erate)) {	// if the gap size difference is within our error rate, it is OK
                                        uint32 tmpoff = (offset < 0 ? 0 : offset);
                                        uint32 gapPatch = 0;
                                        if (matchingSequencePositions[iid].bgn < matchingSequencePositions[iid].end) {
                                            if (matchingSequencePositions[iid].bgn < MIN_DIST_TO_RECRUIT) {
                                                gapPatch = MIN_DIST_TO_RECRUIT - matchingSequencePositions[iid].bgn;
                                            }
                                            min = (matchingSequencePositions[iid].bgn >= (MIN_DIST_TO_RECRUIT) ? matchingSequencePositions[iid].bgn - (MIN_DIST_TO_RECRUIT): 0);
                                            max = (matchingSequencePositions[iid].end+(MIN_DIST_TO_RECRUIT) < waGlobal->frgToLen[iid] ? matchingSequencePositions[iid].end+(MIN_DIST_TO_RECRUIT) : waGlobal->frgToLen[iid]);
                                        } else {
                                            if (matchingSequencePositions[iid].bgn+MIN_DIST_TO_RECRUIT >= waGlobal->frgToLen[iid]) {
                                                gapPatch = MIN_DIST_TO_RECRUIT - (waGlobal->frgToLen[iid] - matchingSequencePositions[iid].bgn);
                                            }
                                            min = (matchingSequencePositions[iid].bgn+(MIN_DIST_TO_RECRUIT) < waGlobal->frgToLen[iid] ? matchingSequencePositions[iid].bgn+(MIN_DIST_TO_RECRUIT) : waGlobal->frgToLen[iid]);
                                            max = (matchingSequencePositions[iid].end >= (MIN_DIST_TO_RECRUIT) ? matchingSequencePositions[iid].end - (MIN_DIST_TO_RECRUIT) : 0);
                                        }

                                        if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "Found read %d (will be called %d) sequence %d - %d (min: %d, max: %d) that could help gap %d to %d in %d (gap diff is %d)\n", iid, gapIID, matchingSequencePositions[iid].bgn, matchingSequencePositions[iid].end, min, max, gapStart, MIN(iter->position.bgn, iter->position.end), layRecord.iid, diff);

                                        uint32 MIN_FRM_START = 50;
                                        // ran off the end of the patchee
                                        if (gapStart + gapPatch - tmpoff < MIN_DIST_TO_RECRUIT) {
                                            if (min < max) {
                                                min += (MIN_DIST_TO_RECRUIT - (gapStart + gapPatch - tmpoff)) + MIN_FRM_START;
                                            } else {
                                                min -= (MIN_DIST_TO_RECRUIT - (gapStart +gapPatch - tmpoff)) + MIN_FRM_START;
                                            }
                                        }
                                        layout << "{TLE\nclr:"
                                                << 0 /*min*/
                                                << ","
                                                << (MAX(max, min) - MIN(max, min) + 1)
                                                << "\noff:" << (gapStart+gapPatch-tmpoff >= (MIN_DIST_TO_RECRUIT) ? gapStart + gapPatch - (MIN_DIST_TO_RECRUIT) - tmpoff: MIN_FRM_START)
                                                << "\nsrc:" << gapIID
                                                << "\n}\n";
                                        readsWithGaps[iid] = 1;
                                        pair<AS_IID, pair<uint32, uint32> > gapInfo(gapIID++, pair<uint32, uint32>(min, max));
                                        gaps[iid].push_back(gapInfo);
                                        count++;
                                    } else {
                                        if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "Found read %d sequence %d - %d that could help gap %d to %d in %d (but diff %d is too big (vs %d aka %f)\n", iid, matchingSequencePositions[iid].bgn, matchingSequencePositions[iid].end, gapStart, MIN(iter->position.bgn, iter->position.end), layRecord.iid, diff, gapSize, (double)diff/gapSize);
                                    }
                                }
                            }
                            if (count > 0) {
                                if (offset < 0) {
                                    offset = 0;
                                }
                                uint32 overlappingStart = (lastEnd-offset >= (MIN_DIST_TO_RECRUIT / 3) ? lastEnd - (MIN_DIST_TO_RECRUIT / 3) - offset: 0);
                                uint32 overlappingEnd = MIN(iter->position.bgn, iter->position.end) + (MIN_DIST_TO_RECRUIT / 3) - offset;
                                // record this gap
                                if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "For fragment %d with %d supporters had a gap from %d to %d inserting range %d from %d %d with offset %d so in original read positions are %d %d\n", layRecord.iid, count, lastEnd, MIN(iter->position.bgn, iter->position.end),gapIID, overlappingStart, overlappingEnd, offset, overlappingStart+offset, overlappingEnd+offset);
                                layout << "{TLE\nclr:"
                                        << 0
                                        << ","
                                        << overlappingEnd - overlappingStart
                                        << "\noff:" << overlappingStart
                                        << "\nsrc:" << gapIID
                                        << "\n}\n";
                                readsWithGaps[layRecord.iid] = 1;
                                pair<AS_IID, pair<uint32, uint32> > gapInfo(gapIID++, pair<uint32, uint32>(overlappingStart+offset, MIN(waGlobal->frgToLen[layRecord.iid], overlappingEnd+offset)));
                                gaps[layRecord.iid].push_back(gapInfo);
                            } else {
                                if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "For fragment %d had a gap from %d to %d but no one had a matching gap size it so breaking it\n", layRecord.iid, lastEnd, MIN(iter->position.bgn, iter->position.end));
                                closeRecord(waGlobal, outFile, reportFile, layout, layRecord, lastEnd, offset, readIID, readSubID);
                            }
                        } else { // no one could agree with this read, break it
                            /*
                            if (MIN(iter->position.bgn, iter->position.end) - lastEnd < 2*FUDGE_BP) {
                                if (offset < 0) {
                                    offset = 0;
                                }
                                uint32 overlappingStart = (lastEnd-offset >= (MIN_DIST_TO_RECRUIT / 3) ? lastEnd - (MIN_DIST_TO_RECRUIT / 3) - offset: 0);
                                uint32 overlappingEnd = MIN(iter->position.bgn, iter->position.end) + (MIN_DIST_TO_RECRUIT / 3) - offset;
                                // record this gap
                                if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "For fragment %d with 0 supporters had a gap from %d to %d inserting range %d from %d %d with offset %d so in original read positions are %d %d\n", layRecord.iid, lastEnd, MIN(iter->position.bgn, iter->position.end),gapIID, overlappingStart, overlappingEnd, offset, overlappingStart+offset, overlappingEnd+offset);
                                layout << "{TLE\nclr:"
                                        << 0
                                        << ","
                                        << overlappingEnd - overlappingStart
                                        << "\noff:" << overlappingStart
                                        << "\nsrc:" << gapIID
                                        << "\n}\n";
                                readsWithGaps[layRecord.iid] = 1;
                                pair<AS_IID, pair<uint32, uint32> > gapInfo(gapIID++, pair<uint32, uint32>(overlappingStart+offset, MIN(waGlobal->frgToLen[layRecord.iid], overlappingEnd+offset)));
                                gaps[layRecord.iid].push_back(gapInfo);
                                if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "For fragment %d had a gap from %d to %d, no one believed it but it was very small so we allowed.\n", layRecord.iid, lastEnd, MIN(iter->position.bgn, iter->position.end));
                            } else {
                            */
                            if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "For fragment %d had a gap from %d to %d but no one believed it so breaking it\n", layRecord.iid, lastEnd, MIN(iter->position.bgn, iter->position.end));
                            closeRecord(waGlobal, outFile, reportFile, layout, layRecord, lastEnd, offset, readIID, readSubID);
                            //}
                        }
                    } else {
                        // close the layout and start a new one because we have a coverage gap
                        if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "Gap in sequence %d between positions %d - %d is greater than allowed %d\n", layRecord.iid, lastEnd, MIN(iter->position.bgn, iter->position.end), waGlobal->maxUncorrectedGap);
                        closeRecord(waGlobal, outFile, reportFile, layout, layRecord, lastEnd, offset, readIID, readSubID);
                    }
                }
                if (offset < 0) {
                    offset = MIN(iter->position.bgn, iter->position.end);
                    //fprintf(stderr, "Updating offset in layout "F_IID" to be "F_U32"\n", i, offset);
                }
                SeqInterval bClr = layRecord.bClrs[iter->ident];
                uint32 min = MIN(iter->position.bgn, iter->position.end);
                uint32 max = MAX(iter->position.bgn, iter->position.end);
                uint32 length = max - min;

                //fprintf(stderr, "Writing layout for frg "F_IID" "F_IID" "F_U32" "F_U32" "F_U32"\n", i, iter->ident, iter->position.bgn, iter->position.end, offset);
                layout << "{TLE\nclr:"
                        << (iter->position.bgn < iter->position.end ? bClr.bgn : bClr.end)
                        << ","
                        << (iter->position.bgn < iter->position.end ? bClr.end : bClr.bgn)
                        << "\noff:" <<  MIN(iter->position.bgn, iter->position.end)-offset
                        << "\nsrc:" << iter->ident
                        << "\n}\n";
                readsToPrint[iter->ident]=1;
                if (lastEnd < (min + length - 1)) {
                    lastEnd = min + length - 1;
                }
            }
            if (lastEnd - offset >= waGlobal->minLength) {
                fprintf(outFile, "%s}\n", layout.str().c_str());
                fprintf(reportFile, "%d\t%s_%d_%d\t%d\t%d\t%d\n", layRecord.iid, waGlobal->libName, layRecord.iid, readSubID, readSubID, offset, lastEnd);
            }
            readIID++;
            if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "Finished processing read %d subsegments %d\n", layRecord.iid, readSubID);
        }

        if (waGlobal->verboseLevel >= VERBOSE_OFF) fprintf(stderr, "Beginning output of "F_SIZE_T" reads\n", readsToPrint.size());
        loadSequence(waGlobal->gkp, readsToPrint, frgToEnc);
        map<AS_IID, char*> *theFrgs = &frgToEnc;
        for (map<AS_IID, uint8>::const_iterator iter = readsToPrint.begin(); iter != readsToPrint.end(); iter++) {
            if (iter->second == 0) {
                continue;
            }
            if ((*theFrgs)[iter->first] == 0) {
                fprintf(stderr, "Error no ID for read %d\n",iter->first);
            }
            decodeSequenceQuality((*theFrgs)[iter->first], (char*) &seq, (char *) &qlt);
            assert(iter->first < gapIID);
            fprintf(outFile, "{RED\nclr:%d,%d\neid:%d\niid:%d\nqlt:\n%s\n.\nseq:\n%s\n.\n}\n", 0, waGlobal->frgToLen[iter->first], iter->first, iter->first, qlt, seq);
        }
        // output uncorrected sequences
        loadSequence(waGlobal->gkp, readsWithGaps, frgToEnc);
        for (map<AS_IID, uint8>::const_iterator iter = readsWithGaps.begin(); iter != readsWithGaps.end(); iter++) {
            if (iter->second == 0) {
                continue;
            }
            if (frgToEnc[iter->first] == 0) {
                fprintf(stderr, "Error no ID for read %d\n",iter->first);
            }
            if (gaps.find(iter->first) == gaps.end()) {
                fprintf(stderr, "No gap list for read %d\n", iter->first);
            }
            decodeSequenceQuality(frgToEnc[iter->first], (char*) &seq, (char *) &qlt);
            for (vector<pair<AS_IID, pair<uint32, uint32> > >::const_iterator j = gaps[iter->first].begin(); j != gaps[iter->first].end(); j++) {
                fprintf(outFile, "{RED\nclr:%d,%d\neid:%d\niid:%d\nqlt:\n%s\n.\nseq:\n%s\n.\n}\n", j->second.first, j->second.second, j->first, j->first, qlt, seq);
            }
        }

        for (map<AS_IID, char*>::iterator iter = frgToEnc.begin(); iter != frgToEnc.end(); iter++) {
            delete[] iter->second;
        }
        frgToEnc.clear();
        fclose(outFile);
        fclose(reportFile);
        closeLayFile(inFile);

        if (waGlobal->maxUncorrectedGap > 0) {
            inStore->unlink();
            delete inStore;
        }
        AS_UTL_unlink(inName);
    fprintf(stderr, "Done output\n");
}

