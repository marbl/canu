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

#include "AS_OVS_overlapStore.h"

#include <set>
#include <map>

static const char *rcsid_AS_PBR_FILTER_C = "$Id: AS_PBR_filter.cc,v 1.1 2012-06-28 20:02:45 skoren Exp $";

void computeRepeatThreshold(PBRThreadGlobals &thread_globals, AS_IID firstFrag) {
   // create files to partition the scores
   FILE **partitionedScores = new FILE*[thread_globals.partitions];
   char outputName[FILENAME_MAX] = {0};
   for (int i = 0; i < thread_globals.partitions; i++) {
	  sprintf(outputName, "%s.%d.rank", thread_globals.prefix, i+1);
	  errno = 0;
	  fprintf(stderr, "Trying to open file %d named %s holding range %d-%d\n", i+1, outputName, thread_globals.partitionStarts[i].first, thread_globals.partitionStarts[i].second);
	  partitionedScores[i] = fopen(outputName, "w");
	  if (errno) {
		 fprintf(stderr, "Couldn't open '%s' for write: %s\n", outputName, strerror(errno)); exit(1);
	  }
   }

	double mean = 0;
	if (thread_globals.globalRepeats == TRUE) {
	   if (thread_globals.covCutoff == 0) {
		  // compute global coverage of the reads we're correcting
		  // the number of mappings each high-identity read has should be equal to the coverage of our data we're correcting, except for repeats
		  // identify the cut off the peak
		  // we do this by finding the inflection point at the end of the curve (that is where we are past the peak and no longer decreasing)
		  double N = 0;
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

				covHist[MIN(MAX_COV_HIST, val)]++;
			 }
		  }

		  double prevRatio = 0;
		  double runningTotal = covHist[0] + covHist[1];
		  for (uint16 iter = 2; iter < MAX_COV_HIST; iter++) {
			 if (covHist[iter] <= covHist[iter-1] && (double)covHist[iter]/covHist[iter-1] > prevRatio && (runningTotal / N > CUMULATIVE_SUM)) {
				thread_globals.covCutoff = iter - 1;
				break;
			 }
			 prevRatio = (double)covHist[iter]/covHist[iter-1];
			 runningTotal += covHist[iter];
		   }
		   if (thread_globals.covCutoff == 0) thread_globals.covCutoff = MAX_COV;
		   delete[] covHist;
		}
		fprintf(stderr, "Picking cutoff as %d mean would be %f\n", thread_globals.covCutoff, mean);

		// now that we have a cutoff, stream the store and record which pacbio sequences the high-identity sequences should correct
		OverlapStore *ovs = AS_OVS_openOverlapStore(thread_globals.ovlStoreUniqPath);
		uint64 olapCount = 0;
		uint64 ovlPosition = 0;
		OVSoverlap *olaps = NULL;

		for (map<uint32, uint8>::const_iterator iter = thread_globals.readsToPrint.begin(); iter != thread_globals.readsToPrint.end(); iter++) {
		  if (iter->second == 0) {
			 continue;
		  }

		  while (olaps != NULL && ovlPosition < olapCount && olaps[ovlPosition].a_iid < iter->first) {
			 ovlPosition++;
		  }

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

			 fprintf(stderr, "Loaded "F_U64" overlaps\n", olapCount);
		  }

		  // build a sorted by score map of all mapping an illumina sequence has
		  multimap<uint64, AS_IID> scoreToReads;
		  uint32 alen = thread_globals.frgToLen[iter->first];
		  for (uint64 rank = ovlPosition; rank < olapCount; rank++, ovlPosition++) {
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
			 scoreToReads.insert(pair<uint64, AS_IID>(currScore, olaps[rank].b_iid));
		  }

		  // now keep only the best cuttoff of those and store them for easy access
		  uint64 lastScore = 0;
		  uint16 position = 0;
		  set<AS_IID> readRanking;
		  for (multimap<uint64, AS_IID>::reverse_iterator rank = scoreToReads.rbegin(); rank != scoreToReads.rend(); rank++) {
			 if (readRanking.find(rank->second) != readRanking.end()) {
				continue;  // only record each pac bio read at most once
			 }
			 position = (position == MAX_COV ? MAX_COV : (lastScore == rank->first ? position+1 : position+1));
			 if (position > thread_globals.covCutoff) {
				break;
			 }
			 readRanking.insert(rank->second);
			 lastScore = rank->first;
		   }

		   // output the score to the appropriate partition
		   for (set<AS_IID>::iterator j = readRanking.begin(); j != readRanking.end(); j++) {
			  uint32 mpLow = MIN(thread_globals.partitions-1, MAX(0, (uint32) ceil((double)((*j)-firstFrag+1) / thread_globals.perFile)-1));
			  uint32 mpHigh= MIN(thread_globals.partitions-1, MAX(0, (uint32) floor((double)((*j)-firstFrag+1) / thread_globals.perFile)-1));
			  if (mpLow == mpHigh && mpLow > 0) { mpLow--; }
			  uint32 min = MIN(thread_globals.partitionStarts[mpLow].first, thread_globals.partitionStarts[mpHigh].first);
			  uint32 max = MAX(thread_globals.partitionStarts[mpLow].second, thread_globals.partitionStarts[mpHigh].second);
			  uint32 tmp = MIN(mpLow, mpHigh);
			  mpHigh = MAX(mpLow, mpHigh);
			  mpLow = tmp;

			  if (thread_globals.partitionStarts[mpLow].first <= (*j) && thread_globals.partitionStarts[mpLow].second > (*j)) {
				 fprintf(partitionedScores[mpLow], F_IID"\t"F_IID"\n", iter->first, (*j));
			  } else if (thread_globals.partitionStarts[mpHigh].first <= (*j) && thread_globals.partitionStarts[mpHigh].second > (*j)) {
				 fprintf(partitionedScores[mpHigh], F_IID"\t"F_IID"\n", iter->first, (*j));
			  } else {
				 bool error = true;
				 if ((*j) < min) {
					while (mpLow >= 0 && thread_globals.partitionStarts[mpLow].first > (*j)) {
					   mpLow--;
					}
					if (thread_globals.partitionStarts[mpLow].first < (*j) && thread_globals.partitionStarts[mpLow].second > (*j)) {
					   fprintf(partitionedScores[mpLow], F_IID"\t"F_IID"\n", iter->first, (*j));
					   error = false;
					}
				 } else {
					while (mpHigh < thread_globals.partitions && thread_globals.partitionStarts[mpHigh].second <= (*j)) {
					   mpHigh++;
					}
					if (thread_globals.partitionStarts[mpHigh].first <= (*j) && thread_globals.partitionStarts[mpHigh].second > (*j)) {
					   fprintf(partitionedScores[mpHigh], F_IID"\t"F_IID"\n", iter->first, (*j));
					   error = false;
					}
				 }
				 if (error) {
					fprintf(stderr, "ERROR: Could not find appropriate partition for %d though it was either %d or %d (%d-%d) and (%d-%d) but it was neither\n", (*j), mpLow, mpHigh, thread_globals.partitionStarts[mpLow].first, thread_globals.partitionStarts[mpLow].second, thread_globals.partitionStarts[mpHigh].first, thread_globals.partitionStarts[mpHigh].second); exit(1);
				 }
			  }
		   }
		}
		delete[] olaps;
		AS_OVS_closeOverlapStore(ovs);
		if (thread_globals.fixedMemory == FALSE) {
		   loadSequence(thread_globals.gkp, thread_globals.readsToPrint, thread_globals.frgToEnc);
		}
		thread_globals.readsToPrint.clear();
	 }
	 for (int i = 0; i < thread_globals.partitions; i++) {
	   fclose(partitionedScores[i]);
	 }
	 delete[] partitionedScores;
}
