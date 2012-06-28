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

#include "AS_OVS_overlapStore.h"
#include "AS_UTL_reverseComplement.h"
#include "AS_PER_encodeSequenceQuality.h"

#include <sstream>
#include <map>
#include <vector>
#include <set>

static const char *rcsid_AS_PBR_OUTPUT_C = "$Id: AS_PBR_output.cc,v 1.1 2012-06-28 20:02:45 skoren Exp $";

void *outputResults(void *ptr) {
  PBRThreadWorkArea *wa = (PBRThreadWorkArea *) ptr;
  PBRThreadGlobals *waGlobal = wa->globals;

  //drand48_data rstate;
  //srand48_r(1, &rstate);
  pair<AS_IID, AS_IID> bounds(0,0);
  int part = 0;

  while (true) {
     if (waGlobal->numThreads > 1) {
        pthread_mutex_lock(&waGlobal->countMutex);
     }
     if (waGlobal->toOutput.size() == 0) {
        if (waGlobal->numThreads > 1) {
           pthread_mutex_unlock(&waGlobal->countMutex);
        }
        break;
     } else {
        part = waGlobal->toOutput.size();
fprintf(stderr, "THe thread %d has to output size of "F_SIZE_T" and partitions %d\n", wa->id, waGlobal->toOutput.size(), waGlobal->partitions);
        bounds = waGlobal->toOutput.top();
        waGlobal->toOutput.pop();
        if (waGlobal->numThreads > 1) {
           pthread_mutex_unlock(&waGlobal->countMutex);
        }
     }
     part = bounds.first;
     map<AS_IID, uint8> readsToPrint;
     map<AS_IID, uint8> readsWithGaps;
     map<AS_IID, vector<pair<AS_IID, pair<uint32, uint32> > > > gaps;
     map<AS_IID, set<AS_IID> > readRanking;
     map<AS_IID, char*> frgToEnc;

     char outputName[FILENAME_MAX] = {0};
     sprintf(outputName, "%s.%d.lay", waGlobal->prefix, part);
     errno = 0;
     FILE *outFile = fopen(outputName, "w");
     if (errno) {
        fprintf(stderr, "Couldn't open '%s' for write: %s\n", outputName, strerror(errno)); exit(1);
     }

     char inName[FILENAME_MAX] = {0};
     sprintf(inName, "%s.%d.olaps", waGlobal->prefix, part);
     errno = 0;
     FILE *inFile = fopen(inName, "r");
     if (errno) {
        fprintf(stderr, "Couldn't open '%s' for write: %s from %d-%d\n", inName, strerror(errno), waGlobal->partitionStarts[part].first, waGlobal->partitionStarts[part].second);
        assert(waGlobal->partitionStarts[part-1].first == waGlobal->partitionStarts[part-1].second  && waGlobal->partitionStarts[part-1].first == 0);
        fclose(outFile);
        continue;
     }

     char inRankName[FILENAME_MAX] = {0};
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

     fprintf(stderr, "Thread %d is running and output to file %s range %d-%d\n", wa->id, outputName, bounds.first, bounds.second);
     uint32 readIID = 0;
     char seq[AS_READ_MAX_NORMAL_LEN];
     char qlt[AS_READ_MAX_NORMAL_LEN];
     AS_IID gapIID = MAX(1, waGlobal->partitionStarts[bounds.second].second + 1);
     gapIID = MAX(gapIID, waGlobal->partitionStarts[bounds.first].second + 1);
     fprintf(stderr, "Thread %d is output starting with gap sequence %d\n", wa->id, gapIID);

     while (!feof(inFile)) {
        uint32 readSubID = 1;

        AS_IID i;
        vector<OverlapPos> mp;
        map<AS_IID, SeqInterval> bclrs;

        // read in a record
        uint32 count = 0;
        fscanf(inFile, "LAY\t"F_IID"\t"F_U32"\n", &i, &count);
        for (uint32 iter = 0; iter < count; iter++) {
           OverlapPos o;
           SeqInterval bclr;
           fscanf(inFile, "TLE\t"F_IID"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\n", &o.ident, &o.position.bgn, &o.position.end, &bclr.bgn, &bclr.end);
//fprintf(stderr, "Read in record "F_IID" with start "F_U32" "F_U32"\n", o.ident, o.position.bgn, o.position.end);
           mp.push_back(o);
           bclrs[o.ident] = bclr;
        }

        stringstream layout (stringstream::in | stringstream::out);
        layout << "{LAY\neid:" << i << "_" << readSubID << "\niid:" << readIID << "\n";
        uint32 lastEnd = 0;
        int32 offset = -1;

        // process record
        for (vector<OverlapPos>::const_iterator iter = mp.begin(); iter != mp.end(); iter++) {
           // skip reads over coverage
           if (waGlobal->globalRepeats == TRUE && (readRanking[iter->ident].find(i) == readRanking[iter->ident].end())) {
              //fprintf(stderr, "Skipping read %d to correct %d it was at cutoff %d true %d\n", iter->ident, i, waGlobal->readRanking[iter->ident][i].first, waGlobal->readRanking[iter->ident][i].second);
              continue;
           }
           if (lastEnd != 0 && lastEnd < MIN(iter->position.bgn, iter->position.end)) {
              // instead of  skipping, we will output the uncorrected sequence
              if (MIN(iter->position.bgn, iter->position.end) - lastEnd < waGlobal->maxUncorrectedGap) {
                 if (offset < 0) {
                    offset = 0;
                 }
                 uint32 overlappingStart = (lastEnd-offset >= (waGlobal->maxUncorrectedGap / 3) ? lastEnd - (waGlobal->maxUncorrectedGap / 3) - offset: 0);
                uint32 overlappingEnd = MIN(iter->position.bgn, iter->position.end) + (waGlobal->maxUncorrectedGap / 3) - offset;
                // record this gap
                fprintf(stderr, "For fragment %d had a gap from %d to %d inserting range %d from %d %d\n", i, lastEnd, MIN(iter->position.bgn, iter->position.end),gapIID, overlappingStart, overlappingEnd);
                layout << "{TLE\nclr:"
                       << 0
                       << ","
                       << overlappingEnd - overlappingStart
                       << "\noff:" << overlappingStart
                       << "\nsrc:" << gapIID
                       << "\n}\n";
               readsWithGaps[i] = 1;
               pair<AS_IID, pair<uint32, uint32> > gapInfo(gapIID++, pair<uint32, uint32>(overlappingStart, overlappingEnd));
               gaps[i].push_back(gapInfo);
              } else {
                 // close the layout and start a new one because we have a coverage gap
                 if (lastEnd - offset >= waGlobal->minLength) {
                    fprintf(outFile, "%s}\n", layout.str().c_str());
                 }
                 readIID++;
                 readSubID++;
                 offset = -1;
                 layout.str("");
                 layout << "{LAY\neid:" << i << "_" << readSubID << "\niid:" << readIID << "\n";
              }
           }
           if (offset < 0) {
              offset = MIN(iter->position.bgn, iter->position.end);
//fprintf(stderr, "Updating offset in layout "F_IID" to be "F_U32"\n", i, offset);
           }
           SeqInterval bClr;
           bClr.bgn = 0;
           bClr.end = waGlobal->frgToLen[iter->ident];
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
        }
        readIID++;
        fprintf(stderr, "Finished processing read %d subsegments %d\n", i, readSubID);
     }

     fprintf(stderr, "Thread %d beginning output of "F_SIZE_T" reads\n", wa->id, readsToPrint.size());
     if (waGlobal->fixedMemory == TRUE) {
        if (waGlobal->numThreads > 1) {
           pthread_mutex_lock(&waGlobal->globalDataMutex);
        }
        loadSequence(waGlobal->gkp, readsToPrint, frgToEnc);
        if (waGlobal->numThreads > 1) {
           pthread_mutex_unlock(&waGlobal->globalDataMutex);
        }
     }
     map<AS_IID, char*> *theFrgs = (waGlobal->fixedMemory == TRUE ? &frgToEnc : &waGlobal->frgToEnc);
     for (map<AS_IID, uint8>::const_iterator iter = readsToPrint.begin(); iter != readsToPrint.end(); iter++) {
        if (iter->second == 0) {
           continue;
        }
        if ((*theFrgs)[iter->first] == 0) {
           fprintf(stderr, "Error no ID for read %d\n",iter->first);
        }
        decodeSequenceQuality((*theFrgs)[iter->first], (char*) &seq, (char *) &qlt);
        fprintf(outFile, "{RED\nclr:%d,%d\neid:%d\niid:%d\nqlt:\n%s\n.\nseq:\n%s\n.\n}\n", 0, waGlobal->frgToLen[iter->first], iter->first, iter->first, qlt, seq);
     }
     // output uncorrected sequences
     if (waGlobal->numThreads > 1) {
        pthread_mutex_lock(&waGlobal->globalDataMutex);
     }
     loadSequence(waGlobal->gkp, readsWithGaps, frgToEnc);
     if (waGlobal->numThreads > 1) {
        pthread_mutex_unlock(&waGlobal->globalDataMutex);
     }
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
     fclose(inFile);
  }
  fprintf(stderr, "Done output of thread %d\n", wa->id);
  return(NULL);
}

