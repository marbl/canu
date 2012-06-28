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

#include <vector>
#include <algorithm>

static const char *rcsid_AS_PBR_CORRECT_C = "$Id: AS_PBR_correct.cc,v 1.1 2012-06-28 20:02:45 skoren Exp $";

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
static SeqInterval findChimera(AS_IID aid, uint32 alen, const map<AS_IID, SeqInterval> &fwdMatches, const map<AS_IID, SeqInterval> &revMatches) {
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
	if (totalInBothOri > 0 && (chimeraJunction.end - chimeraJunction.end) < CHIMERA_MAX_SIZE) {
		fprintf(stderr, "The pac read %d has a chimera with %d reads mapping both fwd/rev (out of %d and %d) at positions %d-%d\n", aid, totalInBothOri, fwdMatches.size(), revMatches.size(), chimeraJunction.bgn, chimeraJunction.end);
	} else {
		fprintf(stderr, "The pac read %d is not a chimera with %d reads mapping both, ( %d %d ) \n", aid, totalInBothOri, chimeraJunction.bgn, chimeraJunction.end);
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

void *  correctFragments(void *ptr) {
  PBRThreadWorkArea *wa = (PBRThreadWorkArea *) (ptr);
  PBRThreadGlobals* waGlobal = wa->globals;
  uint32 counter = 0;

  fprintf(stderr, "Starting in thread %d from %d to %d\n", wa->id, wa->start, wa->end);

  OverlapStore *ovs = AS_OVS_openOverlapStore(waGlobal->ovlStoreUniqPath);

  // local copies of global variables
  map<AS_IID, uint8> readsToPrint;
  map<AS_IID, uint32> longReadsToPrint;

  uint16  *readCoverage = new uint16[AS_READ_MAX_NORMAL_LEN];
  map<AS_IID, OVSoverlap> frgToBest;
  map<AS_IID, uint64> frgToScore;
  OVSoverlap olap;
  uint64 olapCount = 0;
  uint64 ovlPosition = 0;
  OVSoverlap *olaps = NULL;

  // compute the number of files we should be generating
  char outputName[FILENAME_MAX] = {0};
  uint32 fileStart = wa->fileStart;
  uint32 lastFile = wa->fileEnd;
  uint32 numFiles = lastFile - fileStart + 1;
  uint32 perFile = waGlobal->perFile;

  pair<AS_IID, AS_IID> partitionStartEnd;
  uint32 currOpenID = fileStart;
  sprintf(outputName, "%s.%d.olaps", waGlobal->prefix, fileStart);
  errno = 0;
  FILE *outFile = fopen(outputName, "w");
  if (errno) {
     fprintf(stderr, "Couldn't open '%s' for write: %s\n", outputName, strerror(errno)); exit(1);
  }
  fprintf(stderr, "In thread %d going to output files %d-%d with %d\n", wa->id, fileStart, fileStart + numFiles - 1, perFile);
  partitionStartEnd.first = wa->start;

  // finally stream through the fragments we're correcting and map the best overlaps to them
  // figure out our block of responsibility
  for (AS_IID i = wa->start; i <= wa->end; i++) {
     if (waGlobal->libToInclude[waGlobal->frgToLib[i]] != TRUE) {
       continue;
    }

    if (counter % 10000 == 0) {
       fprintf(stderr, "Thread %d done with %d fragments\n", wa->id, counter);
     }

     uint32 alen = waGlobal->frgToLen[i];
     memset(readCoverage, 0, alen * sizeof(uint16));

     frgToBest.clear();
     frgToScore.clear();

     map<uint32, OverlapPos> tile;
     map<AS_IID, SeqInterval> bClrs;

     map<AS_IID, uint8> toSkip;

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
        fprintf(stderr, "Thread %d loaded "F_U64" overlaps\n", wa->id, olapCount);
     }

     map<AS_IID, SeqInterval> fwdMatches;
     map<AS_IID, SeqInterval> revMatches;

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
        convertOverlapToPosition(olap, pos, bClr, alen, blen);

        if (isOvlForward(olap)) {
        	fwdMatches[bid] = pos;
        } else {
        	revMatches[bid] = pos;
        }

        if (isOlapBad(olap, alen, blen, waGlobal->erate, waGlobal->elimit, waGlobal->maxErate)) {
           continue;
        }

        // figure out what bases the bfrag covers
        if (olap.dat.ovl.type == AS_OVS_TYPE_OVL && (olap.dat.ovl.a_hang < 0 || olap.dat.ovl.b_hang > 0)) {
           // non contained overlap, dont use these fragments for correction
           if (frgToScore[bid] != 0) {
              OVSoverlap best = frgToBest[bid];
              uint32 min = MIN(tile[best.b_iid].position.bgn, tile[best.b_iid].position.end);
              uint32 max = MAX(tile[best.b_iid].position.bgn, tile[best.b_iid].position.end);
              for (uint32 iter = min; iter <= max && waGlobal->globalRepeats == FALSE; iter++) {
                 readCoverage[iter]--;
              }
              if (longReadsToPrint[best.b_iid] != 0) {
                 longReadsToPrint[best.b_iid]--;
                 if (longReadsToPrint[best.b_iid] < MAX_COV) {
                    readsToPrint[best.b_iid] = longReadsToPrint[best.b_iid];
                    longReadsToPrint[best.b_iid] = 0;
                 }
              } else {
                 readsToPrint[best.b_iid]--;
              }
              tile[best.b_iid].position.bgn = tile[best.b_iid].position.end = 0;
           }
           toSkip[olap.b_iid] = 1;
           continue;
        }

        if (toSkip[bid] != 0) {
           continue;
        }

        // dont use library itself to correct fragments
        if (waGlobal->libToInclude[waGlobal->frgToLib[bid]] == TRUE) {
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
           if (longReadsToPrint[best.b_iid] != 0) {
              longReadsToPrint[best.b_iid]--;
              if (longReadsToPrint[best.b_iid] < MAX_COV) {
                 readsToPrint[best.b_iid] = longReadsToPrint[best.b_iid];
                 longReadsToPrint[best.b_iid] = 0;
              }
           } else {
              readsToPrint[best.b_iid]--;
           }
           tile[best.b_iid].position.bgn = tile[best.b_iid].position.end = 0;
        }
        frgToBest[bid] = olap;
        frgToScore[bid] = score;

        bClrs.insert(pair<AS_IID, SeqInterval>(bid, bClr));

        OverlapPos tileStr;
        tileStr.position = pos;
        tileStr.ident = bid;
        tile[bid] = tileStr;

        if (readsToPrint[bid] == MAX_COV) {
           if (longReadsToPrint.find(bid) == longReadsToPrint.end()) { longReadsToPrint[bid] = MAX_COV; }
           longReadsToPrint[bid]++;
        } else {
           readsToPrint[bid]++;
        }

        uint32 min = MIN(pos.bgn, pos.end);
        uint32 max = MAX(pos.bgn, pos.end);
        for (uint32 iter = min; iter <= max && waGlobal->globalRepeats == FALSE; iter++) {
           readCoverage[iter] = (readCoverage[iter] == MAX_COV ? MAX_COV : readCoverage[iter]+1);
        }
     }

     SeqInterval chimeraJunction = findChimera(i, alen, fwdMatches, revMatches);
     fwdMatches.clear();
     revMatches.clear();

     vector<OverlapPos> mp; // = new vector<OverlapPos>();
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
        		mp.push_back(iter->second);
        	}
        }
     }

     if (mp.size() > 0) {
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
           stable_sort(mp.begin(), mp.end(), compare_by_identity);
           globalFrgToScore = NULL;
           if (waGlobal->numThreads > 1) {
              pthread_mutex_unlock( &waGlobal->globalDataMutex);
           }

           for (vector<OverlapPos>::iterator iter = mp.begin(); iter != mp.end(); ) {
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
                 iter = mp.erase(iter);
              } else {
                 iter++;
              }
           }
        }

        // now save the tiling
        stable_sort(mp.begin(), mp.end(), compare_tile);
        // write to my file
        uint32 fileId = fileStart + MIN(numFiles - 1, floor((double)(i - wa->start) / perFile));
        if (fileId != currOpenID) {
           partitionStartEnd.second = i;
fprintf(stderr, "Thread %d set partition %d to be %d-%d\n", wa->id, currOpenID-1, partitionStartEnd.first, partitionStartEnd.second);
           waGlobal->partitionStarts[currOpenID-1] = partitionStartEnd;
           partitionStartEnd.first = partitionStartEnd.second = i;
           fclose(outFile);
           sprintf(outputName, "%s.%d.olaps", waGlobal->prefix, fileId);
           errno = 0;
           outFile = fopen(outputName, "w");
           if (errno) {
              fprintf(stderr, "Couldn't open '%s' for write: %s\n", outputName, strerror(errno)); exit(1);
           }
fprintf(stderr, "For thread %d I'm currently on file %d and need to move to %d\n", wa->id, currOpenID, fileId);
           currOpenID = fileId;
        }
        fprintf(outFile, "LAY\t"F_IID"\t"F_SIZE_T"\n", i, mp.size());
        for (vector<OverlapPos>::const_iterator iter = mp.begin(); iter != mp.end(); iter++) {
        	fprintf(outFile, "TLE\t"F_IID"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\n", iter->ident, iter->position.bgn, iter->position.end, bClrs[iter->ident].bgn, bClrs[iter->ident].end);
        }
     }

     counter++;
   }
   partitionStartEnd.second = wa->end + 1;
   waGlobal->partitionStarts[currOpenID-1] = partitionStartEnd;
fprintf(stderr, "Thread %d set partition %d to be %d-%d\n", wa->id, currOpenID-1, partitionStartEnd.first, partitionStartEnd.second);
   fclose(outFile);

   // finally update the global data
   if (waGlobal->numThreads > 1) {
      pthread_mutex_lock( &waGlobal->globalDataMutex);
   }
   for (map<AS_IID, uint8>::iterator iter = readsToPrint.begin(); iter != readsToPrint.end(); iter++) {
      if (((uint32)waGlobal->readsToPrint[iter->first] + iter->second) > MAX_COV) {
         waGlobal->readsToPrint[iter->first] = MAX_COV;
         if (waGlobal->longReadsToPrint.find(iter->first) == waGlobal->longReadsToPrint.end()) { waGlobal->longReadsToPrint[iter->first] = waGlobal->readsToPrint[iter->first]; }
         waGlobal->longReadsToPrint[iter->first] += iter->second;
      } else {
         waGlobal->readsToPrint[iter->first] += iter->second;
      }
   }
   if (waGlobal->numThreads > 1) {
      pthread_mutex_unlock( &waGlobal->globalDataMutex);
   }

   fprintf(stderr, "Thread shutting down after finishing %d fragments\n", counter);
   delete[] readCoverage;
   delete[] olaps;
   readsToPrint.clear();
   AS_OVS_closeOverlapStore(ovs);

   return 0;
}
