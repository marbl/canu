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

#include "AS_PBR_mates.hh"
#include "MultiAlignMatePairAnalysis.H"
#include "MultiAlign.h"

#include <map>
#include <vector>

const double 	Z_MAX	= 6;

/*  The following JavaScript functions for calculating normal and
    chi-square probabilities and critical values were adapted by
    John Walker from C implementations
    written by Gary Perlman of Wang Institute, Tyngsboro, MA
    01879.  Both the original C code and this JavaScript edition
    are in the public domain.  */

double poz(double z) {
	double y, x, w;

    if (z == 0.0) {
        x = 0.0;
    } else {
        y = 0.5 * fabs(z);
        if (y > (Z_MAX * 0.5)) {
            x = 1.0;
            z = 0.0;
        } else if (y < 1.0) {
            w = y * y;
            x = ((((((((0.000124818987 * w
                     - 0.001075204047) * w + 0.005198775019) * w
                     - 0.019198292004) * w + 0.059054035642) * w
                     - 0.151968751364) * w + 0.319152932694) * w
                     - 0.531923007300) * w + 0.797884560593) * y * 2.0;
        } else {
            y -= 2.0;
            x = (((((((((((((-0.000045255659 * y
                           + 0.000152529290) * y - 0.000019538132) * y
                           - 0.000676904986) * y + 0.001390604284) * y
                           - 0.000794620820) * y - 0.002034254874) * y
                           + 0.006549791214) * y - 0.010557625006) * y
                           + 0.011630447319) * y - 0.009279453341) * y
                           + 0.005353579108) * y - 0.002141268741) * y
                           + 0.000535310849) * y + 0.999936657524;
        }
    }
    return z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5);
}

void *estimateInsertSizes(void *ptr) {
	PBRThreadWorkArea *wa = (PBRThreadWorkArea *) ptr;
	PBRThreadGlobals *waGlobal = wa->globals;

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
			bounds = waGlobal->toOutput.top();
			waGlobal->toOutput.pop();
			if (waGlobal->numThreads > 1) {
				pthread_mutex_unlock(&waGlobal->countMutex);
			}
		}
		part = bounds.first;

		char inName[FILENAME_MAX] = {0};
		sprintf(inName, "%s.%d.olaps", waGlobal->prefix, part);
		errno = 0;
		FILE *inFile = fopen(inName, "r");
		if (errno) {
			fprintf(stderr, "Couldn't open '%s' for write: %s from %d-%d\n", inName, strerror(errno), waGlobal->partitionStarts[part].first, waGlobal->partitionStarts[part].second);
			assert(waGlobal->partitionStarts[part-1].first == waGlobal->partitionStarts[part-1].second  && waGlobal->partitionStarts[part-1].first == 0);
			continue;
		}

		MultiAlignT *ma = CreateEmptyMultiAlignT();
		while (!feof(inFile)) {
			uint32 readSubID = 1;

			AS_IID i;

			// read in a record
			uint32 count = 0;
			fscanf(inFile, "LAY\t"F_IID"\t"F_U32"\n", &i, &count);
			ma->maID = i;
			ma->data.num_frags = count;
			ResetToRange_VA(ma->f_list, ma->data.num_frags);

			for (uint32 iter = 0; iter < count; iter++) {
				OverlapPos o;
				SeqInterval bclr;
				fscanf(inFile, "TLE\t"F_IID"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\n", &o.ident, &o.position.bgn, &o.position.end, &bclr.bgn, &bclr.end);

				IntMultiPos  *imp = GetIntMultiPos(ma->f_list, iter);
				imp->ident        = o.ident;
				imp->contained    = false;
				imp->parent       = 0;
				imp->ahang        = 0;
				imp->bhang        = 0;
				imp->position.bgn = o.position.bgn;
				imp->position.end = o.position.end;
			}

			if (waGlobal->numThreads > 1) {
				waGlobal->mpa->evaluateTig(ma, &waGlobal->globalDataMutex);
			} else {
				waGlobal->mpa->evaluateTig(ma);
			}
			ClearMultiAlignT(ma);
		}
		fclose(inFile);
	}

	return (NULL);
}

void *screenBadMates(void *ptr) {
  PBRThreadWorkArea *wa = (PBRThreadWorkArea *) ptr;
  PBRThreadGlobals *waGlobal = wa->globals;

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
        bounds = waGlobal->toOutput.top();
        waGlobal->toOutput.pop();
        if (waGlobal->numThreads > 1) {
           pthread_mutex_unlock(&waGlobal->countMutex);
        }
     }
     part = bounds.first;

     char outputName[FILENAME_MAX] = {0};
     sprintf(outputName, "%s.%d.paired.olaps", waGlobal->prefix, part);
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

     fprintf(stderr, "Thread %d is running and output to file %s range %d-%d\n", wa->id, outputName, bounds.first, bounds.second);

     while (!feof(inFile)) {
        AS_IID i;
        vector<OverlapPos> mp;
        map<AS_IID, SeqInterval> bclrs;
        map<AS_IID, int32> readStarts;
        map<AS_IID, bool> readOri;
        uint32 currLength = 0;

        // read in a record
        uint32 count = 0;
        fscanf(inFile, "LAY\t"F_IID"\t"F_U32"\n", &i, &count);
        for (uint32 iter = 0; iter < count; iter++) {
           OverlapPos o;
           SeqInterval bclr;
           fscanf(inFile, "TLE\t"F_IID"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\n", &o.ident, &o.position.bgn, &o.position.end, &bclr.bgn, &bclr.end);
           mp.push_back(o);
           bclrs[o.ident] = bclr;
           if (currLength < MAX(o.position.bgn, o.position.end)) {
        	   currLength = MAX(o.position.bgn, o.position.end);
           }

           AS_IID libID = waGlobal->frgToLib[o.ident];
           if (waGlobal->libToOrientation[libID] != AS_READ_ORIENT_UNKNOWN) {
			   // save the beginning position of each sequence
			   // we can save the beginning position because if the read is fwd it looks like:
			   //	PacRead: -------------------------->
			   //	fwd read:		--->
			   //	rev read:				<---
			   // so the begin is always the outer-most positions of the sequences
			   readStarts[o.ident] = o.position.bgn;
			   readOri[o.ident] = (o.position.bgn < o.position.end);
           }
        }

        map<AS_IID, bool> processed;
        map<AS_IID, bool> good;
        for (vector<OverlapPos>::const_iterator iter = mp.begin(); iter != mp.end(); iter++) {
        	if (processed.find(iter->ident) != processed.end()) {
        		continue;
        	}
        	processed[iter->ident] = true;
        	int32 mypos = readStarts[iter->ident];
        	AS_IID otherRead = waGlobal->frgToMate[iter->ident];
        	AS_IID libID = waGlobal->frgToLib[iter->ident];
        	pair<double, double> libSize = waGlobal->libToSize[libID];
        	double mean = (libSize.second + libSize.first) / 2;
        	double stdev = (libSize.second - libSize.first) / ( 2 * CGW_CUTOFF);

    		// unmated keep
        	if (otherRead == 0) {
        		good[iter->ident] = true;
        		continue;
        	}
        	if (readStarts.find(otherRead) == readStarts.end()) {
        		// mate not in same read, keep with some probability
        		// first, we need to compute the distance from the end of the read to this sequence
        		// next, we compute the probability of observing a paired-end that has distance > this value given our library distribution
        		// we drop the sequence with P(1-probability of observing pair)

        		// compute distance from end of containing read
        		uint32 dist = 0;
        		if (readOri[iter->ident] == 0) {
        			dist = iter->position.bgn;
        		} else {
        			dist = currLength - iter->position.bgn + 1;
        		}

        		// calculate the probability of having a sequence longer than this distance, given the mean/stdev
        		double zScore = (dist - mean) / stdev;
        		double twoTailedProb = poz(zScore);
        		assert(twoTailedProb >= 0 && twoTailedProb <= 1);
        		if (drand48() < twoTailedProb) {
        			good[iter->ident] = true;
        		}
        		continue;
        	}

        	// mark the mate as processed too
        	processed[otherRead] = true;

        	// make sure this fragment comes first, it should because we sorted our fragments by start position before
                if (readStarts[otherRead] <= readStarts[iter->ident]) {
                   if (MIN(iter->position.bgn, iter->position.end) <= readStarts[otherRead]) {
                      // one read is contained in the other or they've gone past each other, ignore it
                      continue;
                    }
                    assert(0);
                }

        	uint32 dist = readStarts[otherRead] - readStarts[iter->ident];
        	// check the orientation of the reads
        	switch(waGlobal->libToOrientation[libID]) {
        		case AS_READ_ORIENT_INNIE:
        			if (readOri[iter->ident] == true && readOri[otherRead] == false && libSize.first < dist && dist < libSize.second) {
        				// good pair, keep
        				good[iter->ident] = true;
        				good[otherRead] = true;
        			} else {
        				// bad mate;
        			}
        			break;
        		case AS_READ_ORIENT_OUTTIE:
        			if (readOri[iter->ident] == false && readOri[otherRead] == true && libSize.first < dist && dist < libSize.second) {
        				// good pair, keep
        				good[iter->ident] = true;
        				good[otherRead] = true;
        			} else {
        				// bad mate
        			}
        			break;
        		case AS_READ_ORIENT_NORMAL:
        		case AS_READ_ORIENT_ANTINORMAL:
        			fprintf(stderr, "Not supported yet!");
        			assert(0);
        			break;
        	}
        }

        if (good.size() > 0) {
			// finally, output the sequences we decided to keep
			fprintf(outFile, "LAY\t"F_IID"\t"F_SIZE_T"\n", i, good.size());
	        for (vector<OverlapPos>::const_iterator iter = mp.begin(); iter != mp.end(); iter++) {
			   AS_IID ident = iter->ident;

			   if (good.find(ident) != good.end()) {
				   SeqInterval bclr = bclrs[ident];
				   fprintf(outFile, "TLE\t"F_IID"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\n", ident, iter->position.bgn, iter->position.end, bclr.bgn, bclr.end);
			   }
			}
        }
     }
     fclose(inFile);
     fclose(outFile);
  }

  return (NULL);
}
