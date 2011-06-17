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
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id: CorrectPacBio.cc,v 1.1 2011-06-17 12:58:42 skoren Exp $";

#include <map>
#include <set>
#include <list>
#include <vector>
#include <stack>
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <sstream>
#include <pthread.h>

using namespace std;

#include "AS_global.h"
#include "AS_UTL_reverseComplement.h"
#include "AS_MSG_pmesg.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_PBR_util.hh"

#define  THREAD_STACKSIZE        (16 * 512 * 512)
const uint8 MAX_COV     = 255;
const double CUMULATIVE_SUM = 0.95;
const uint64 MAX_TO_READ = 100000;

map<AS_IID, uint64> *globalFrgToScore;

// global variables shared between all threads
struct PBRThreadGlobals {
   // mutex to controll access to overlap store
   pthread_mutex_t  overlapMutex;

   // mutex to controll access to gkp store/gkpStore itself
   pthread_mutex_t gkpMutex;
   gkStore *gkp;

   // writable global data (access controlled by globalDataMutex)
   pthread_mutex_t globalDataMutex;

   map<AS_IID, uint32> readsToPrint;
   map<AS_IID, vector<IntMultiPos> *> correctedTilings;
   map<AS_IID, map<AS_IID, SeqInterval> *> readBClrs;
   map<AS_IID, multimap<uint64, AS_IID> > readsToBestMatches;

   // track number of active threads for output of layouts
   stack<pair<AS_IID, AS_IID> > toOutput;
   pthread_mutex_t countMutex;

   // global parameters
   char *    ovlStoreUniqPath;
   int       numThreads;
   int       partitions;
   double    maxErate;
   double    erate;
   double    elimit;
   int       globalRepeats;
   double    repeatMultiplier;
   int       minLength;
   char      prefix[FILENAME_MAX];

   // read-only variables for thread
   uint8 covCutoff;
   map<AS_IID, uint32> frgToLen;
   map<AS_IID, uint32> frgToLib;
   map<AS_IID, map<AS_IID, pair<uint8, uint8> > > readRanking;
   map<AS_IID, string> frgToSeq;
   map<AS_IID, string> frgToQlt;
};

// this holds thread-specfic variables (plus a pointer to the global data)
struct PBRThreadWorkArea {
   uint32 id;
   AS_IID start;
   AS_IID end;
   PBRThreadGlobals *globals;
};

bool compare_by_identity (const IntMultiPos &a, const IntMultiPos &b) {
   if (globalFrgToScore == NULL || (*globalFrgToScore)[a.ident] == 0 || (*globalFrgToScore)[b.ident] == 0) {
      fprintf(stderr, "Error: fragment %d or %d does not have a defined score\n", a.ident, b.ident);
      assert(0);
   }

   return (*globalFrgToScore)[a.ident] < (*globalFrgToScore)[b.ident];
}

bool compare_tile ( const IntMultiPos& a, const IntMultiPos& b )
{
  return MIN(a.position.bgn, a.position.end) < MIN(b.position.bgn, b.position.end);
}

static uint32 loadFragments(gkStream *fs, uint32* includeLib, map<AS_IID, uint32>& frgToLen, map<AS_IID, uint32>& frgToLib, map<AS_IID, string> &frgToSeq, map<AS_IID, string> &frgToQlt) {
  gkFragment  fr;
  uint32 counter = 0;
  
  fprintf(stderr, "Streaming fragments\n");
  // figure out which libraries we want to use and store fragment info
  while (fs->next(&fr)) {
     if (includeLib[fr.gkFragment_getLibraryIID()] == TRUE) {
        frgToLib[fr.gkFragment_getReadIID()] = TRUE;
        counter++;
     } else {
        uint32 clrL, clrR;
        fr.gkFragment_getClearRegion(clrL, clrR);
        frgToSeq[fr.gkFragment_getReadIID()] = string(fr.gkFragment_getSequence()).substr(clrL, clrR);
        frgToQlt[fr.gkFragment_getReadIID()] = string(fr.gkFragment_getQuality()).substr(clrL, clrR);
     }
     frgToLen[fr.gkFragment_getReadIID()] = fr.gkFragment_getClearRegionLength();
  }

  return counter;
}

// partition the work
static AS_IID partitionWork( uint32 counter, map<AS_IID, uint32>& frgToLib, int numThreads, PBRThreadWorkArea *wa) { 
  uint32 lastEnd = 0;
  uint32 currThread = 0;
  AS_IID lastFrag = 0;
  PBRThreadGlobals* waGlobal = wa[0].globals;

  uint32 perThread = (uint32) floor((double)counter / numThreads);
  fprintf(stderr, "Each thread responsible for %d (%d threads) of %d fragments\n", perThread, numThreads, counter);

  counter = 0;

  for(map<AS_IID, uint32>::const_iterator iter = frgToLib.begin(); iter != frgToLib.end(); iter++) {
     if (iter->second == TRUE) { 
        if (counter == 0) { lastEnd = iter->first; }
        counter++;
     }
     if (currThread < waGlobal->numThreads -1 ) {
        if ((counter > 0 && counter % perThread == 0)) {
           wa[currThread].start = lastEnd;
           wa[currThread].end = iter->first-1;
           wa[currThread].id = currThread;
           lastEnd = iter->first;
           currThread++;
         }
      }
     lastFrag = (iter->first > lastFrag ? iter->first : lastFrag);
  }
  // give the rest to the last thread
  wa[currThread].start = lastEnd;
  wa[currThread].end = lastFrag;
  wa[currThread].id = currThread;

  return lastFrag;
}

static void *  correctFragments(void *ptr) {
  PBRThreadWorkArea *wa = (PBRThreadWorkArea *) (ptr);
  PBRThreadGlobals* waGlobal = wa->globals;
  uint32 counter = 0;

  fprintf(stderr, "Starting in thread %d from %d to %d\n", wa->id, wa->start, wa->end);

  OverlapStore *ovs = AS_OVS_openOverlapStore(waGlobal->ovlStoreUniqPath);

  // local copies of global variables
  map<AS_IID, uint32> readsToPrint;
  map<AS_IID, vector<IntMultiPos> *> correctedTilings;
  map<AS_IID, map<AS_IID, SeqInterval> *> readBClrs;
  map<AS_IID, multimap<uint64, AS_IID> > readsToBestMatches;

  uint8  *readCoverage = new uint8[AS_READ_MAX_NORMAL_LEN];
  map<AS_IID, OVSoverlap> frgToBest;
  map<AS_IID, uint64> frgToScore;
  OVSoverlap olap;
  uint64 olapCount = 0;
  uint64 ovlPosition = 0;
  OVSoverlap *olaps = NULL;

  // finally stream through the fragments we're correcting and map the best overlaps to them
  // figure out our block of responsibility
  for (AS_IID i = wa->start; i <= wa->end; i++) {
     if (waGlobal->frgToLib[i] != TRUE) {
       continue;
    }

    if (counter % 10000 == 0) {
       fprintf(stderr, "Thread %d done with %d fragments\n", wa->id, counter);
     }

     uint32 alen = waGlobal->frgToLen[i];
     memset(readCoverage, 0, alen * sizeof(uint8));

     frgToBest.clear();
     frgToScore.clear();

     map<uint32, IntMultiPos> tile;
     map<AS_IID, SeqInterval> *bClrs = new map<AS_IID, SeqInterval>();

     // read next batch
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
        fprintf(stderr, "Thread %d loaded %llu overlaps\n", wa->id, olapCount);
     }

     for (; ovlPosition < olapCount; ovlPosition++) {
        if (olaps[ovlPosition].a_iid != i) {
           break;
        }
        olap = olaps[ovlPosition];

        AS_IID aid = olap.a_iid;
        AS_IID bid = olap.b_iid;
        uint32 blen = waGlobal->frgToLen[bid];
        uint64 score = scoreOverlap(olap, alen, blen, waGlobal->erate, waGlobal->elimit, waGlobal->maxErate);

        if (isOlapBad(olap, alen, blen, waGlobal->erate, waGlobal->elimit, waGlobal->maxErate)) {
           continue;
        }

        // figure out what bases the bfrag covers
        if (olap.dat.ovl.type == AS_OVS_TYPE_OVL && (olap.dat.ovl.a_hang < 0 || olap.dat.ovl.b_hang > 0)) {
           continue;
        }

        // dont use library itself to correct fragments
        if (waGlobal->frgToLib[bid] == TRUE) {
           continue;
        }

        // remove mapping if we find a better version of it to this read
        if (score  < frgToScore[bid]) {
           continue;
        }

        if (frgToScore[bid] < score && frgToScore[bid] != 0) {
           OVSoverlap best = frgToBest[bid];
           uint32 min = MIN(tile[best.b_iid].position.bgn, tile[best.b_iid].position.end);
           uint32 max = MAX(tile[best.b_iid].position.bgn, tile[best.b_iid].position.end);
           for (uint32 iter = min; iter <= max && waGlobal->globalRepeats == FALSE; iter++) {
              readCoverage[iter]--;
           }
           readsToPrint[best.b_iid]--;
           tile[best.b_iid].position.bgn = tile[best.b_iid].position.end = 0;
        }
        frgToBest[bid] = olap;
        frgToScore[bid] = score; 
        
        SeqInterval pos;
        SeqInterval bClr;
        if (olap.dat.ovl.type == AS_OVS_TYPE_OVL) {
           pos.bgn = olap.dat.ovl.a_hang;
           pos.end = alen + olap.dat.ovl.b_hang;
           bClr.bgn = 0;
           bClr.end = blen;

           if (olap.dat.ovl.flipped) {
              uint32 x = pos.end;
              pos.end = pos.bgn;
              pos.bgn = x;
           }
        } else if (olap.dat.ovl.type == AS_OVS_TYPE_OBT) {
           pos.bgn = olap.dat.obt.a_beg;
           pos.end = olap.dat.obt.a_end;

           if (!olap.dat.obt.fwd) {
              pos.bgn = olap.dat.obt.a_end;
              pos.end = olap.dat.obt.a_beg;
           }
           uint32 bend = olap.dat.obt.b_end_hi << 9 | olap.dat.obt.b_end_lo;
           bClr.bgn = MIN(olap.dat.obt.b_beg, bend);
           bClr.end = MAX(olap.dat.obt.b_beg, bend);
        }

        // update end points if necessary
        uint32 len = MAX(pos.bgn, pos.end) - MIN(pos.bgn, pos.end);
        if (len > waGlobal->frgToLen[bid]) {
           if (pos.bgn > pos.end) {
              pos.bgn = pos.end + waGlobal->frgToLen[bid];
           } else {
              pos.end = pos.bgn + waGlobal->frgToLen[bid];
          }
        }
        IntMultiPos tileStr;
        tileStr.position = pos;
        tileStr.ident = bid;
        tileStr.parent = aid;
        tile[bid] = tileStr;

        bClrs->insert(pair<AS_IID, SeqInterval>(bid, bClr));
        readsToPrint[bid]++;

        uint32 min = MIN(pos.bgn, pos.end);
        uint32 max = MAX(pos.bgn, pos.end);
        for (uint32 iter = min; iter <= max && waGlobal->globalRepeats == FALSE; iter++) {
           readCoverage[iter] = (readCoverage[iter] == MAX_COV ? MAX_COV : readCoverage[iter]+1);
        }
     }

     vector<IntMultiPos>* mp = new vector<IntMultiPos>();
     for (map<uint32, IntMultiPos>::const_iterator iter = tile.begin(); iter != tile.end(); iter++) {
        if (iter->second.position.bgn > 0 || iter->second.position.end > 0) {
           mp->push_back(iter->second);

           // save the best matches for our ranking
           uint64 score = frgToScore[iter->second.ident];
           readsToBestMatches[iter->second.ident].insert(pair<uint64, AS_IID>(score, i));
        }
     }

     if (mp->size() > 0) {
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

           memset(readCoverage, 0, alen * sizeof(uint8));

           if (waGlobal->numThreads > 1) {
              pthread_mutex_lock( &waGlobal->globalDataMutex);
           }
           globalFrgToScore = &frgToScore;
           stable_sort(mp->begin(), mp->end(), compare_by_identity);
           globalFrgToScore = NULL;
           if (waGlobal->numThreads > 1) {
              pthread_mutex_unlock( &waGlobal->globalDataMutex);
           }

           for (vector<IntMultiPos>::iterator iter = mp->begin(); iter != mp->end(); ) {
              uint32 min = MIN(iter->position.bgn, iter->position.end);
              uint32 max = MAX(iter->position.bgn, iter->position.end);
              uint8  max_cov = 0;

              for (uint32 coverage = min; coverage <= max; coverage++) {
                 if (max_cov < readCoverage[coverage]) {
                    max_cov = readCoverage[coverage];
                 }
                 readCoverage[coverage] = (readCoverage[coverage] == MAX_COV ? MAX_COV : readCoverage[coverage]+1);
              }
              if (max_cov > (mean * waGlobal->repeatMultiplier)) {
                 iter = mp->erase(iter);
                 readsToPrint[iter->ident]--;
              } else {
                 iter++;
              }
           }
        }

        // now save the tiling
        stable_sort(mp->begin(), mp->end(), compare_tile);
        correctedTilings[i] = mp;
        readBClrs[i] = bClrs;
     }

     counter++;
   }

   // finally update the global data
   if (waGlobal->numThreads > 1) {
      pthread_mutex_lock( &waGlobal->globalDataMutex);
   }
   for (map<AS_IID, uint32>::const_iterator iter = readsToPrint.begin(); iter != readsToPrint.end(); iter++) {
      waGlobal->readsToPrint[iter->first] += iter->second;
   }
   if (waGlobal->numThreads > 1) {
      pthread_mutex_unlock( &waGlobal->globalDataMutex);
   }

   if (waGlobal->numThreads > 1) {
      pthread_mutex_lock( &waGlobal->globalDataMutex);
   }
   for (AS_IID i = wa->start; i <= wa->end; i++) {
      waGlobal->correctedTilings[i] = correctedTilings[i];
      waGlobal->readBClrs[i] = readBClrs[i];
   }
   if (waGlobal->numThreads > 1) {
      pthread_mutex_unlock( &waGlobal->globalDataMutex);
   }

   if (waGlobal->numThreads > 1) {
      pthread_mutex_lock( &waGlobal->globalDataMutex);
   }
   for (map<AS_IID, multimap<uint64, AS_IID> >::iterator iter = readsToBestMatches.begin(); iter != readsToBestMatches.end(); iter++) {
      multimap<uint64, AS_IID> *multMap = &(iter->second);
      multimap<uint64, AS_IID> *globalMap = &waGlobal->readsToBestMatches[iter->first];
      for (multimap<uint64, AS_IID>::const_iterator j = multMap->begin(); j != multMap->end(); j++) {
         globalMap->insert(pair<uint64, AS_IID>(j->first, j->second));
      }
   }
   if (waGlobal->numThreads > 1) {
      pthread_mutex_unlock( &waGlobal->globalDataMutex);
   }

   fprintf(stderr, "Thread shutting down after finishing %d fragments\n", counter);
   delete[] readCoverage;
   readsToPrint.clear();
   correctedTilings.clear();
   readBClrs.clear();
   readsToBestMatches.clear();
   AS_OVS_closeOverlapStore(ovs);

   return 0;
}

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
fprintf(stderr, "THe thread %d has tootput size of %d and partitiosn %d\n", wa->id, waGlobal->toOutput.size(), waGlobal->partitions);
        bounds = waGlobal->toOutput.top();
        waGlobal->toOutput.pop();
        if (waGlobal->numThreads > 1) {
           pthread_mutex_unlock(&waGlobal->countMutex);
        }
     }
     char outputName[FILENAME_MAX] = {0};
     sprintf(outputName, "%s.%d.lay", waGlobal->prefix, part);
     errno = 0;
     FILE *outFile = fopen(outputName, "w");
     if (errno) {
        fprintf(stderr, "Couldn't open '%s' for write: %s\n", outputName, strerror(errno)); exit(1);
     }

     fprintf(stderr, "Thread %d is running and output to file %s range %d-%d\n", wa->id, outputName, bounds.first, bounds.second);
     map<AS_IID, uint32> readsToPrint;
     uint32 readIID = 0;

     for (AS_IID i = bounds.first; i <= bounds.second; i++) {
        uint32 readSubID = 1;
        stringstream layout (stringstream::in | stringstream::out);
        layout << "{LAY\neid:" << i << "_" << readSubID << "\niid:" << readIID << "\n";
        uint32 lastEnd = 0;        int32 offset = -1;

        vector<IntMultiPos>* mp = waGlobal->correctedTilings[i];
        if (mp == NULL) continue;
        for (vector<IntMultiPos>::const_iterator iter = mp->begin(); iter != mp->end(); iter++) {
           // skip reads over coverage
           if (waGlobal->globalRepeats == TRUE && waGlobal->readRanking[iter->ident][i].first > waGlobal->covCutoff) {
              //fprintf(stderr, "Skipping read %d to correct %d it was at cutoff %d true %d\n", iter->ident, i, waGlobal->readRanking[iter->ident][i].first, waGlobal->readRanking[iter->ident][i].second);
              continue;
           } else if (waGlobal->globalRepeats == TRUE && waGlobal->readRanking[iter->ident][i].second > waGlobal->covCutoff) {
              //fprintf(stderr, "Read %d to correct %d is at equal position %d but at true position %d", iter->ident, i, waGlobal->readRanking[iter->ident][i].first, waGlobal->readRanking[iter->ident][i].second);
              int numInstances = (uint32) floor((double)waGlobal->readsToPrint[iter->ident] / waGlobal->covCutoff);
              double randDouble = 0;
              //drand48_r(&rstate, &randDouble);
              int randVal = 1 + (int) (numInstances * (randDouble / (RAND_MAX + 1.0)));
              //if (randVal != 1) { continue; }
              //fprintf(stderr, "\n");
           }
           if (lastEnd != 0 && lastEnd < MIN(iter->position.bgn, iter->position.end)) {
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

           if (offset < 0) {
              offset = MIN(iter->position.bgn, iter->position.end);
           }
           SeqInterval bClr = waGlobal->readBClrs[i]->find(iter->ident)->second; //waGlobal->readBClrs[i][iter->ident];
           uint32 min = MIN(iter->position.bgn, iter->position.end);
           uint32 max = MAX(iter->position.bgn, iter->position.end);
           uint32 length = max - min;

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

     fprintf(stderr, "Thread %d beginning output of %d reads\n", wa->id, readsToPrint.size());
     for (map<uint32, uint32>::const_iterator iter = readsToPrint.begin(); iter != readsToPrint.end(); iter++) {
        if (iter->second == 0) {
           continue;
        }
        string seq = waGlobal->frgToSeq[iter->first]; //seqCache[cacheRead].substr(clrL, clrR);
        string qlt = waGlobal->frgToQlt[iter->first]; // qltCache[cacheRead].substr(clrL, clrR);

        fprintf(outFile, "{RED\nclr:%d,%d\neid:%d\niid:%d\nqlt:\n%s\n.\nseq:\n%s\n.\n}\n", 0, waGlobal->frgToLen[iter->first], iter->first, iter->first, qlt.c_str(), seq.c_str());
     }

     fclose(outFile);
  }
  fprintf(stderr, "Done output of thread %d\n", wa->id);
}

int
main (int argc, char * argv []) {
  PBRThreadGlobals  thread_globals;

  char      *gkpStorePath           = NULL;
  thread_globals.ovlStoreUniqPath   = NULL;

  thread_globals.numThreads        = 2;
  thread_globals.maxErate          = 0.25;
  thread_globals.erate             = 0.15;
  thread_globals.elimit            = 4.5;
  thread_globals.globalRepeats     = TRUE;
  thread_globals.repeatMultiplier  = 2.0;
  thread_globals.partitions        = 100;
  thread_globals.minLength         = 500;
  strcpy(thread_globals.prefix, "asm");

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if (strcmp(argv[arg], "-G") == 0) {
      gkpStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      if      (thread_globals.ovlStoreUniqPath == NULL)
        thread_globals.ovlStoreUniqPath = argv[++arg];
      else
        err++;

    } else if (strcmp(argv[arg], "-p") == 0) {
       thread_globals.partitions = atoi(argv[++arg]);
       if (thread_globals.partitions <= 0) { thread_globals.partitions = 1; }

    } else if (strcmp(argv[arg], "-o") == 0) {
       strncpy(thread_globals.prefix, argv[++arg], FILENAME_MAX);

    } else if (strcmp(argv[arg], "-l") == 0) {
       thread_globals.minLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      thread_globals.erate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      thread_globals.elimit = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-c") == 0) {
      thread_globals.maxErate = atof(argv[++arg]);
    
    } else if (strcmp(argv[arg], "-R") == 0) {
       thread_globals.globalRepeats = FALSE;
       thread_globals.repeatMultiplier = atof(argv[++arg]);
    
    } else if (strcmp(argv[arg], "-t") == 0) {
       thread_globals.numThreads = atoi(argv[++arg]);
       if (thread_globals.numThreads <= 0) { thread_globals.numThreads = 1; }
 
    } else {
      err++;
    }

    arg++;
  }

  if ((thread_globals.maxErate < 0.0) || (AS_MAX_ERROR_RATE < thread_globals.maxErate))
    err++;
  if ((thread_globals.erate < 0.0) || (AS_MAX_ERROR_RATE < thread_globals.erate))
    err++;
  if (thread_globals.elimit < 0.0)
    err++;
  if (gkpStorePath == NULL)
    err++;
  if (thread_globals.ovlStoreUniqPath == NULL)
    err++;
  if (thread_globals.minLength <= 0) 
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -O ovlStore -G gkpStore [options]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -O         Mandatory path to an ovlStore.\n");
    fprintf(stderr, "  -G         Mandatory path to a gkpStore.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e 0.15   no more than 0.015 fraction (1.5%%) error\n");
    fprintf(stderr, "  -E 0      no more than 0 errors\n");
    fprintf(stderr, "  -c 0.25   ignore overlaps over this rate before correction\n");
    fprintf(stderr, "  -l %d     ignore corrected fragments less than %d bp\n", thread_globals.minLength, thread_globals.minLength);
    fprintf(stderr, "  -t %d     use %d threads to process correction in parallel\n", thread_globals.numThreads, thread_globals.numThreads);
    fprintf(stderr, "  -p %d     output %d results files, corresponds to #of parallel consensus jobs desired\n", thread_globals.partitions, thread_globals.partitions);
    fprintf(stderr, "  -p %s     output prefix of %s\n", thread_globals.prefix, thread_globals.prefix);
    fprintf(stderr, "  -R %f     consider a pileup of %f times the mean for a single corrected read to be a repeat and distribute reads to their best locations (this is only useful for metagenomic or non-even coverage datasets. Otherwise the global repeat estimate is used instead)\n", thread_globals.repeatMultiplier, thread_globals.repeatMultiplier);
    fprintf(stderr, "\n");

    if ((thread_globals.maxErate < 0.0) || (AS_MAX_ERROR_RATE < thread_globals.maxErate))
      fprintf(stderr, "Invalid overlap error threshold (-c option); must be between 0.00 and %.2f.\n", AS_MAX_ERROR_RATE);

    if ((thread_globals.erate < 0.0) || (AS_MAX_ERROR_RATE < thread_globals.erate))
      fprintf(stderr, "Invalid overlap error threshold (-e option); must be between 0.00 and %.2f.\n",
              AS_MAX_ERROR_RATE);

    if (thread_globals.elimit < 0.0)
      fprintf(stderr, "Invalid overlap error limit (-E option); must be above 0.\n");

    if (gkpStorePath == NULL)
      fprintf(stderr, "No gatekeeper store (-G option) supplied.\n");

    if (thread_globals.ovlStoreUniqPath == NULL)
      fprintf(stderr, "No overlap store (-O option) supplied.\n");

    if (thread_globals.minLength <= 0) 
       fprintf(stderr, "Invalid min length, must be positive: %d\n", thread_globals.minLength);

    exit(1);
  }

  // loop through all the pacBio reads
  thread_globals.gkp = new gkStore(gkpStorePath, FALSE, FALSE, TRUE);
  uint32 i;
  uint32 *includeLib = new uint32[thread_globals.gkp->gkStore_getNumLibraries() + 1];
  memset(includeLib, 0, (thread_globals.gkp->gkStore_getNumLibraries() + 1) * sizeof(uint32));

  for (i=1; i<=thread_globals.gkp->gkStore_getNumLibraries(); i++) {
     gkLibrary      *gkpl = thread_globals.gkp->gkStore_getLibrary(i);
     includeLib[i] |= gkpl->doConsensusCorrection;
  }

  fprintf(stderr, "Correcting fragments\n");
  pthread_attr_t    attr;
  pthread_t  *      thread_id;
  PBRThreadWorkArea* thread_wa;
  thread_id = new pthread_t[thread_globals.numThreads];
  thread_wa = new PBRThreadWorkArea[thread_globals.numThreads];

  // partition data
  thread_wa[0].globals = &thread_globals;
  gkStream   *fs = new gkStream(thread_globals.gkp, 1, thread_globals.gkp->gkStore_getNumFragments(), GKFRAGMENT_QLT);//INF);
  uint32 numFrags = loadFragments(fs, includeLib, thread_globals.frgToLen, thread_globals.frgToLib, thread_globals.frgToSeq, thread_globals.frgToQlt);
  AS_IID lastFrag = partitionWork(numFrags, thread_globals.frgToLib, thread_globals.numThreads, thread_wa);
  delete fs;

  // create parallel threads 
  if  (thread_globals.numThreads > 1) {
    pthread_attr_init (& attr);
    pthread_attr_setstacksize (& attr, THREAD_STACKSIZE);
    pthread_mutex_init (& thread_globals.overlapMutex, NULL);
    pthread_mutex_init (& thread_globals.globalDataMutex, NULL);
    pthread_mutex_init (& thread_globals.gkpMutex, NULL);
    pthread_mutex_init (& thread_globals.countMutex, NULL);
  }

  for (int i = 1; i < thread_globals.numThreads; i++) {
     thread_wa[i].globals = &thread_globals;
     int status = pthread_create(&thread_id[i], &attr, correctFragments, &(thread_wa[i]));
     if (status != 0) {
        fprintf (stderr, "pthread_create error at line %d:  %s\n",  __LINE__, strerror (status));
         exit (-3);
       }
   }
   // create the 0 thread (this thread is always created since we have at least one thread always)
   correctFragments(&(thread_wa[0]));
   for (int i = 1; i < thread_globals.numThreads; i++) {
      void *ptr;
      int status = pthread_join(thread_id[i], &ptr);
      if (status != 0) {
         fprintf (stderr, "pthread_create error at line %d:  %s\n",  __LINE__, strerror (status));
         exit (-3);
      }
   }

   // filter repeat reads out 
   fprintf(stderr, "Filtering repeats\n");

   if (thread_globals.globalRepeats == TRUE) {
      // compute global coverage of the reads we're correcting
      // the number of mappings each high-identity read has should be equal to the coverage of our data we're correcting, except for repeats
      // identify the cut off the peak
      // we do this by finding the inflection point at the end of the curve (that is where we are past the peak and no longer decreasing)
      double mean = 0;
      double N = 0;
      uint32* covHist = new uint32[MAX_COV + 1];
      memset(covHist, 0, (MAX_COV + 1) * sizeof(uint32));
      double prevScore = 0;

      for (map<uint32, uint32>::const_iterator iter = thread_globals.readsToPrint.begin(); iter != thread_globals.readsToPrint.end(); iter++) {
         if (iter->second > 0) {
             N++;
             double delta = iter->second - mean;
             mean += delta / N;

             covHist[MIN(MAX_COV, iter->second)]++;

             uint8 position = 0;
             uint8 positionAbsolute = 0;
             // also build the map for this read of its ranking of corrected reads
             for (multimap<uint64, AS_IID>::const_iterator rank = thread_globals.readsToBestMatches[iter->first].begin(); rank != thread_globals.readsToBestMatches[iter->first].end(); rank++) {
                position = (position == MAX_COV ? MAX_COV : (prevScore == rank->first ? position : position+1));
                positionAbsolute = (positionAbsolute == MAX_COV ? MAX_COV : positionAbsolute+1);
                thread_globals.readRanking[iter->first][rank->second] = pair<uint8, uint8>(position, positionAbsolute);
                prevScore = rank->first;
             }   
         }
      }
      double prevRatio = 0;
      double runningTotal = covHist[0] + covHist[1];
      for (uint8 iter = 2; iter < MAX_COV; iter++) {
         if (covHist[iter] <= covHist[iter-1] && (double)covHist[iter]/covHist[iter-1] > prevRatio && (runningTotal / N > CUMULATIVE_SUM)) {
            thread_globals.covCutoff = iter - 1;
            break;
         }
         prevRatio = (double)covHist[iter]/covHist[iter-1];
         runningTotal += covHist[iter];
       }
       if (thread_globals.covCutoff == 0) thread_globals.covCutoff = MAX_COV;
       thread_globals.readsToBestMatches.clear();
       delete[] covHist;
       fprintf(stderr, "Picking cutoff as %d mean would be %f\n", thread_globals.covCutoff, mean * 2.0);
    }

    // output our tiling
    delete[] thread_wa;
    thread_wa = new PBRThreadWorkArea[thread_globals.numThreads];
    delete[] thread_id;
    thread_id = new pthread_t[thread_globals.numThreads];
    thread_wa[0].globals = &thread_globals;
    thread_wa[0].id = 0;

    fprintf(stderr, "Processing a total of %d tilings\n", thread_globals.correctedTilings.size());
    AS_IID lastEnd = thread_globals.correctedTilings.begin()->first;
    uint32 perThread = (uint32) floor((double)(lastFrag - lastEnd + 1) / thread_globals.partitions);
    for  (uint32 i = 0;  i < thread_globals.partitions; i++ ) {
       pair<AS_IID, AS_IID> bounds;
       bounds.first = lastEnd;
       bounds.second = lastEnd + perThread - 1;
       if (i == thread_globals.partitions - 1) bounds.second = lastFrag;
       thread_globals.toOutput.push(bounds);
       lastEnd+=perThread;
    }
    for  (i = 1;  i < thread_globals.numThreads ;  i ++) {
       thread_wa[i].globals = &thread_globals;
       thread_wa[i].id = i;
       int status = pthread_create(&thread_id[i], &attr, outputResults, &thread_wa[i]);
        if  (status != 0) {
           fprintf (stderr, "pthread_create error at line %d:  %s\n",  __LINE__, strerror (status));
           exit (-3);
        }
    }
    // create the 0 thread
    outputResults(&(thread_wa[0]));
    for  (i = 1;  i < thread_globals.numThreads;  i ++) {
       void  * ptr;
        int status = pthread_join  (thread_id [i], & ptr);
    } 

   // clean up
   fprintf(stderr, "Cleaning up\n");
   if (thread_globals.numThreads > 1) {
      pthread_mutex_destroy(& thread_globals.overlapMutex);
      pthread_mutex_destroy (& thread_globals.globalDataMutex);
      pthread_mutex_destroy (& thread_globals.gkpMutex);
      pthread_mutex_destroy (& thread_globals.countMutex);
   }

   for (map<AS_IID, vector<IntMultiPos>* >::iterator i = thread_globals.correctedTilings.begin(); i != thread_globals.correctedTilings.end(); i++) {
     delete thread_globals.correctedTilings[i->first];
     delete thread_globals.readBClrs[i->first];
   }

   delete[] includeLib;
   delete[] thread_id;
   delete[] thread_wa;
   delete thread_globals.gkp;
}
