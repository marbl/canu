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
#include "AS_PBR_store.hh"
#include "MultiAlignMatePairAnalysis.H"
#include "MultiAlign.h"

#include <map>
#include <vector>

static const char *rcsid_AS_PBR_MATES_C = "$Id: AS_PBR_mates.cc,v 1.2 2012-08-20 13:10:37 skoren Exp $";

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

/**
 * Parallel function to estimate insert-sizes for each library in our short-read sequences
 *
 */
void *estimateInsertSizes(void *ptr) {
	PBRThreadWorkArea *wa = (PBRThreadWorkArea *) ptr;
	PBRThreadGlobals *waGlobal = wa->globals;

	pair<AS_IID, AS_IID> bounds(0,0);
	int part = 0;

	// we have a queue of work, grab a work unit (partition) from the queue and work on it, keep going until the queue is empty
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
		LayRecordStore *inFile = openLayFile(inName);
		if (inFile == NULL) {
			fprintf(stderr, "Couldn't open '%s' for read: %s from %d-%d\n", inName, strerror(errno), waGlobal->partitionStarts[part].first, waGlobal->partitionStarts[part].second);
			assert(waGlobal->partitionStarts[part-1].first == waGlobal->partitionStarts[part-1].second  && waGlobal->partitionStarts[part-1].first == 0);
			continue;
		}

		// build a multialignment for each tiling (pacbio read) in our partition and compute its mate values
		MultiAlignT *ma = CreateEmptyMultiAlignT();
		LayRecord layout;
		while (readLayRecord(inFile, layout, ma)) {
			uint32 readSubID = 1;

			if (waGlobal->numThreads > 1) {
				waGlobal->mpa->evaluateTig(ma, &waGlobal->globalDataMutex);
			} else {
				waGlobal->mpa->evaluateTig(ma);
			}
			ClearMultiAlignT(ma);
		}
		closeLayFile(inFile);
	}

	return (NULL);
}

/**
 * Function to screen out short-reads who are not satifisfied by their mappings
 */
void *screenBadMates(void *ptr) {
  PBRThreadWorkArea *wa = (PBRThreadWorkArea *) ptr;
  PBRThreadGlobals *waGlobal = wa->globals;

  pair<AS_IID, AS_IID> bounds(0,0);
  boost::dynamic_bitset<> *gappedReadSet = initGappedReadSet(waGlobal);
  boost::dynamic_bitset<> *workingReadSet = initGappedReadSet(waGlobal);

  int part = 0;

  // we have a queue of work, grab a work unit (partition) and work on it. Keep going until the queue is empty
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

     // open our input partition and a coressponding output file to record updated tiling with only satisfied mates
     char outputName[FILENAME_MAX] = {0};
     sprintf(outputName, "%s.%d.paired.olaps", waGlobal->prefix, part);
     errno = 0;
     LayRecordStore *outFile = createLayFile(outputName);
     if (errno) {
        fprintf(stderr, "Couldn't open '%s' for write: %s\n", outputName, strerror(errno)); exit(1);
     }

     char inName[FILENAME_MAX] = {0};
     sprintf(inName, "%s.%d.olaps", waGlobal->prefix, part);
     errno = 0;
     LayRecordStore *inFile = openLayFile(inName);
     if (inFile == NULL) {
        fprintf(stderr, "Couldn't open '%s' for write: %s from %d-%d\n", inName, strerror(errno), waGlobal->partitionStarts[part].first, waGlobal->partitionStarts[part].second);
        assert(waGlobal->partitionStarts[part-1].first == waGlobal->partitionStarts[part-1].second  && waGlobal->partitionStarts[part-1].first == 0);
        continue;
     }
     char outputOverlaps[FILENAME_MAX] = {0};
     sprintf(outputOverlaps, "%s.%d.ovb", waGlobal->prefix, part);
     BinaryOverlapFile *bof = AS_OVS_createBinaryOverlapFile(outputOverlaps, FALSE);

     if (waGlobal->verboseLevel >= VERBOSE_DEBUG) fprintf(stderr, "Thread %d is running and output to file %s range %d-%d\n", wa->id, outputName, bounds.first, bounds.second);

     LayRecord layout;
     // work on one PacBio sequence at a time, keeping a subset of its mapped short-read sequences
     while (readLayRecord(inFile, layout, waGlobal)) {
        map<AS_IID, bool> processed;
        map<AS_IID, bool> good;
        for (vector<OverlapPos>::const_iterator iter = layout.mp.begin(); iter != layout.mp.end(); iter++) {
        	if (processed.find(iter->ident) != processed.end()) {
        		continue;
        	}
        	processed[iter->ident] = true;
        	int32 mypos = layout.readStarts[iter->ident];
        	AS_IID otherRead = waGlobal->frgToMate[iter->ident];
        	AS_IID libID = waGlobal->frgToLib[iter->ident];
        	pair<double, double> libSize = waGlobal->libToSize[libID];
        	double mean = (libSize.second + libSize.first) / 2;
        	double stdev = (libSize.second - libSize.first) / ( 2 * CGW_CUTOFF);

    		// unmated sequence, keep
        	if (otherRead == 0) {
        		good[iter->ident] = true;
        		continue;
        	}
        	if (layout.readStarts.find(otherRead) == layout.readStarts.end()) {
        		// mate not in same read, keep with some probability
        		// first, we need to compute the distance from the end of the read to this sequence
        		// next, we compute the probability of observing a paired-end that has distance > this value given our library distribution
        		// we drop the sequence with P(1-probability of observing pair)

        		// compute distance from end of containing read
        		uint32 dist = 0;
        		if (layout.readOri[iter->ident] == 0) {
        			dist = iter->position.bgn;
        		} else {
        			dist = layout.length - iter->position.bgn + 1;
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
            if (layout.readStarts[otherRead] <= layout.readStarts[iter->ident]) {
               if (MIN(iter->position.bgn, iter->position.end) <= layout.readStarts[otherRead]) {
                  // one read is contained in the other or they've gone past each other, ignore it
                  continue;
               }
               fprintf(stderr, "Comparing mated reads %d and %d in %d and first ones position is %d and second is %d which is bad\n", layout.iid, iter->ident, otherRead, layout.readStarts[iter->ident], layout.readStarts[otherRead]);
               assert(0);
            }

        	uint32 dist = layout.readStarts[otherRead] - layout.readStarts[iter->ident];
        	// check the orientation of the reads
        	switch(waGlobal->libToOrientation[libID]) {
        		case AS_READ_ORIENT_INNIE:
        			if (layout.readOri[iter->ident] == true && layout.readOri[otherRead] == false && libSize.first < dist && dist < libSize.second) {
        				// good pair, keep
        				good[iter->ident] = true;
        				good[otherRead] = true;
        			} else {
        				// bad mate;
        			}
        			break;
        		case AS_READ_ORIENT_OUTTIE:
        			if (layout.readOri[iter->ident] == false && layout.readOri[otherRead] == true && libSize.first < dist && dist < libSize.second) {
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

        // finally, record the good subset of our data
        if (good.size() > 0) {
        	writeLayRecord(outFile, layout, workingReadSet, waGlobal->percentShortReadsToStore, &good, bof);
      	   (*gappedReadSet) |= (*workingReadSet);
      	   workingReadSet->reset();
        }
     }
     closeLayFile(inFile);
     closeLayFile(outFile);
     AS_OVS_closeBinaryOverlapFile(bof);
     AS_UTL_unlink(inName);
  }

  if (waGlobal->numThreads > 1) {
     pthread_mutex_lock( &waGlobal->globalDataMutex);
  }
  (*waGlobal->gappedReadSet) |= (*gappedReadSet);
  if (waGlobal->numThreads > 1) {
     pthread_mutex_unlock( &waGlobal->globalDataMutex);
  }
  delete gappedReadSet;
  delete workingReadSet;
  return (NULL);
}
