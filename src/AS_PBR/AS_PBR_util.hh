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

#ifndef AS_PBR_UTIL_H
#define AS_PBR_UTIL_H

static const char *rcsid_AS_PBR_UTIL_H = "$Id$";

#include "AS_global.H"
#include "AS_OVS_overlapStore.H"
#include "AS_PER_gkpStore.H"
#include "MultiAlignMatePairAnalysis.H"

#include <map>
#include <stack>
#include "boost/dynamic_bitset.hpp"

#include <pthread.h>

#define	THREAD_STACKSIZE	(16 * 512 * 512)
#define CGW_CUTOFF 			5
#define	NUM_BITS			8

const uint16 	MAX_COV_HIST					= 	65535;
const uint8 	MAX_COV							= 	255;
const double    FIRST_QUARTILE			=   0.25;
const double    THIRD_QUARTILE                  =   0.75;
const double    ONE_SD_PERCENT                  =   0.341;
const double    TWO_SD_PERCENT                  =   0.136;
const double    CUMULATIVE_SUM                  =   (2*ONE_SD_PERCENT) + (2*TWO_SD_PERCENT);
const uint64 	MAX_TO_READ						= 	100000;
const double	DEFAULT_SAMPLE_SIZE 			= 	0.05;
const double	DEFAULT_SHORT_READ_STORE_SIZE 	= 	0.01;
const uint32	CHIMERA_MAX_SIZE 				= 	150;
const uint32	MIN_DIST_TO_RECRUIT				=	500;
const double    ERATE_ADJUST                    =   1.5;
const uint8     TRIM_BACK_AMOUNT                =   50;
const uint8	MIN_NUM_SAMPLES			=   100;

const uint8		VERBOSE_OFF 					= 0;
const uint8		VERBOSE_DEBUG					= 1;
const uint8		VERBOSE_DEVELOPER				= 2;

struct OverlapPos {
   SeqInterval position;
   IntFragment_ID ident;
};

// global variables shared between all threads
struct PBRThreadGlobals {
   // mutex to controll access to overlap store
   pthread_mutex_t  overlapMutex;

   // mutex to control access to gkp store/gkpStore itself
   pthread_mutex_t gkpMutex;
   gkStore *gkp;

   // writable global data (access controlled by globalDataMutex)
   pthread_mutex_t globalDataMutex;

   map<AS_IID, uint8> readsToPrint;
   map<AS_IID, uint32> longReadsToPrint;
   double percentShortReadsToStore;
   boost::dynamic_bitset<> *gappedReadSet;
   uint32 bitMax;

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
   char		 exePath[FILENAME_MAX];
   double	 percentToEstimateInserts;
   uint8	 verboseLevel;
   int       allowLong;

   // read-only variables for thread
   uint16 covCutoff;
   uint16 maxUncorrectedGap;

   // read-only fragment information
   map<AS_IID, uint32> frgToLen;
   map<AS_IID, uint8> frgToLib;
   map<AS_IID, AS_IID> frgToMate;

   // read-only lib information
   uint32 *libToInclude;
   map<AS_IID, pair<double, double> > libToSize;
   map<AS_IID, uint8 > libToOrientation;
   bool hasMates;
   char libName[LIBRARY_NAME_SIZE];

   matePairAnalysis *mpa;

   int fixedMemory;
   map<AS_IID, char*> frgToEnc;
   pair<AS_IID, AS_IID> *partitionStarts;
   uint32 perFile;
};

// this holds thread-specfic variables (plus a pointer to the global data)
struct PBRThreadWorkArea {
   uint32 id;
   AS_IID start;
   AS_IID end;
   uint32 fileStart;
   uint32 fileEnd;
   PBRThreadGlobals *globals;
};

static boost::dynamic_bitset<>* initGappedReadSet(PBRThreadGlobals *globals, AS_IID max = 0) {
	if (globals->bitMax == 0) {
		assert(max != 0);
		globals->bitMax = max;
	}
	return new boost::dynamic_bitset<>(globals->bitMax + 1);
}

static AS_IID getBitSetID(AS_IID &id, PBRThreadGlobals *globals) {
	assert(globals != NULL);
	assert(id < globals->bitMax);
	return id;
}

static uint32 olapLengthOVL(OVSoverlap ovl, uint32 alen, uint32 blen) {
  int32   ah = ovl.dat.ovl.a_hang;
  int32   bh = ovl.dat.ovl.b_hang;
  uint32  le = 0;

  if (ah < 0) {
    if (bh < 0)
      le = alen + bh;
    else
      le = blen + ah - bh;
  } else {
    if (bh < 0)
      le = alen + bh - ah;
    else
      le = alen - ah;
  }

  return(le);
}

static uint32 olapLengthOBT(OVSoverlap ovl, uint32 alen, uint32 blen) {
   uint32 aovl = ovl.dat.obt.a_end - ovl.dat.obt.a_beg;
   uint32 bend = ovl.dat.obt.b_end_hi >> 9 | ovl.dat.obt.b_end_lo;
   uint32 bbgn = MIN(ovl.dat.obt.b_beg, bend);
   bend = MAX(ovl.dat.obt.b_beg, bend);

   return MIN(aovl, (bend-bbgn));
}

static uint32
olapLength(OVSoverlap ovl, uint32 alen, uint32 blen) {
   if (ovl.dat.ovl.type == AS_OVS_TYPE_OVL) {
      return olapLengthOVL(ovl, alen, blen);
   } else if (ovl.dat.ovl.type == AS_OVS_TYPE_OBT) {
      return olapLengthOBT(ovl, alen, blen);
   }
   return 0;
}

static bool isOlapBadOVL(OVSoverlap olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
    //  The overlap is ALWAYS bad if the original error rate is above what we initially required
    //  overlaps to be at.  We shouldn't have even seen this overlap.  This is a bug in the
    //  overlapper.
    //
    if (olap.dat.ovl.orig_erate > AS_OVS_encodeQuality(AS_CNS_ERROR_RATE))
       return true;

    //  The overlap is GOOD (false == not bad) if the corrected error rate is below the requested
    //  erate.
    //
    if (olap.dat.ovl.corr_erate <= AS_OVS_encodeQuality(AS_UTG_ERROR_RATE)) {
      return false;
    }

    //  If we didn't allow fixed-number-of-errors, the overlap is now bad.  Just a slight
    //  optimization.
    //
    if (AS_UTG_ERROR_LIMIT <= 0)
      return true;

    double olen = olapLength(olap, alen, blen);
    double nerr = olen * AS_OVS_decodeQuality(olap.dat.ovl.corr_erate);

    assert(nerr >= 0);

    if (nerr <= AS_UTG_ERROR_LIMIT) {
       return false;
    }

    return true;
}

static bool isOlapBadOBT(OVSoverlap olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
   if (olap.dat.obt.erate > AS_OVS_encodeQuality(AS_CNS_ERROR_RATE))
      return true;

   if (olap.dat.obt.erate <= AS_OVS_encodeQuality(AS_UTG_ERROR_RATE))
      return false;

   if (AS_UTG_ERROR_LIMIT <= 0) 
      return true;

   double olen = olapLength(olap, alen, blen);
   double nerr = olen * AS_OVS_decodeQuality(olap.dat.obt.erate);
   if (nerr <= AS_UTG_ERROR_LIMIT)
      return false;

   return true;
}

static bool isOlapBad(OVSoverlap olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
    if (olap.dat.ovl.type == AS_OVS_TYPE_OVL) {
       return isOlapBadOVL(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE);
    } else if (olap.dat.ovl.type == AS_OVS_TYPE_OBT) {
       return isOlapBadOBT(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE);
    }
    return true;
}

static uint64  scoreOverlapOVL(const OVSoverlap& olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {

    //  BPW's newer new score.  For the most part, we use the length of the overlap, but we also
    //  want to break ties with the higher quality overlap.
    //
    //  The high 20 bits are the length of the overlap.
    //  The next 12 are the corrected error rate.
    //  The last 12 are the original error rate.
    //
    //  (Well, 12 == AS_OVS_ERRBITS)
    if (isOlapBad(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE)) {
       return 0;
    }

    uint64  leng = 0;
    uint64  corr = (AS_OVS_MAX_ERATE - olap.dat.ovl.corr_erate);
    uint64  orig = (AS_OVS_MAX_ERATE - olap.dat.ovl.orig_erate);

    //  Shift AFTER assigning to a 64-bit value to avoid overflows.
    corr <<= AS_OVS_ERRBITS;

    //  Containments - the length of the overlaps are all the same.  We return the quality.
    //
    if (((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) ||
        ((olap.dat.ovl.a_hang <= 0) && (olap.dat.ovl.b_hang >= 0)))
      return(corr | orig);

    //  Dovetails - the length of the overlap is the score, but we bias towards lower error.
    //  (again, shift AFTER assigning to avoid overflows)
    //
    leng   = olapLength(olap, alen, blen);
    leng <<= (2 * AS_OVS_ERRBITS);

    return(leng | corr | orig);
}

static uint64  scoreOverlapOBT(const OVSoverlap& olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
    if (isOlapBad(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE)) {
       return 0;
    }

    uint64  leng = 0;
    uint64  orig = (AS_OVS_MAX_ERATE - olap.dat.obt.erate);

    leng   = olapLength(olap, alen, blen);
    leng <<= AS_OVS_ERRBITS;

    return(leng | orig);
}

static
  uint64  scoreOverlap(const OVSoverlap& olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
   if (olap.dat.ovl.type == AS_OVS_TYPE_OVL) {
      return scoreOverlapOVL(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE);
   } else if (olap.dat.ovl.type == AS_OVS_TYPE_OBT) {
      return scoreOverlapOBT(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE);
   }
   return 0;
}

static bool compare_tile ( const OverlapPos& a, const OverlapPos& b ) {
  return MIN(a.position.bgn, a.position.end) < MIN(b.position.bgn, b.position.end);
}

extern uint32 loadSequence(gkStore *fs, map<AS_IID, uint8> &readsToPrint, map<AS_IID, char*> &frgToEnc);
extern int32 loadOneSequence(gkStore *fs, AS_IID readIID, char* seq);

static bool isOvlForward(const OVSoverlap& olap) {
	if (olap.dat.ovl.type == AS_OVS_TYPE_OVL) {
	   return !olap.dat.ovl.flipped;
	} else if (olap.dat.ovl.type == AS_OVS_TYPE_OBT) {
		return olap.dat.obt.fwd;
	}

	return false;
}

static bool rangesOverlap(const SeqInterval &first, const SeqInterval &second) {
   uint32 minA = MIN(first.bgn, first.end);
   uint32 minB = MIN(second.bgn, second.end);
   uint32 maxA = MAX(first.bgn, first.end);
   uint32 maxB = MAX(second.bgn, second.end);

   int start = MAX(minA, minB);
   int end = MIN(maxA, maxB);

   return (end-start+1) > 0;
}

extern void convertOverlapToPosition(const OVSoverlap& olap, SeqInterval &pos, SeqInterval &bClr, uint32 alen, uint32 blen, bool forB = false);
#endif
