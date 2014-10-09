
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
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

const char *mainid = "$Id$";

#include "AS_global.H"

#include "AS_PER_gkpStore.H"
#include "intervalList.H"
#include "AS_UTL_decodeRange.H"
#include "AS_OBT_overlaps.H"

#include <vector>
#include <map>

using namespace std;


//  When searching for chimer, each overlap end is trimed back by this amount.
#define OVLTRIM  ((AS_OVERLAP_MIN_LEN / 2) - 1)


//  But then use only overlaps larger than this for some of the more questionable overlaps.
#define MIN_INTERVAL_OVERLAP  60

//  Same orient overlaps closer than this are evidence of PacBio subreads.
#define SUBREAD_LOOP_MAX_SIZE  500
#define SUBREAD_LOOP_EXT_SIZE  2000

//  Original
//#define SUBREAD_LOOP_MAX_SIZE  100
//#define SUBREAD_LOOP_EXT_SIZE  0


//  WITH_REPORT_FULL will ALL overlap evidence.
//  REPORT_OVERLAPS  will print the incoming overlaps in the log.
//
#undef WITH_REPORT_FULL
#undef REPORT_OVERLAPS

#undef DEBUG_ISLINKER
#undef DEBUG_INTERVAL
#undef DEBUG_FILTER
#undef DEBUG_MATES


//  Yuck.  Parameters shouldn't be global.
//
//  minInniePair - a gap with at least this many innie pointing fragment pairs indicating a chimera,
//  is trimmed to the longest confirmed sequence.
//
//  minOverhang - a gap with at least this many overhanging fragments indicating a chimera, is
//  trimmed to the longest confirmed sequence.  Low coverage sequence next to a repeat will have
//  this pattern.
//
//  minInniePair=1  minOverhang=infinity     The historical defaults.
//  minInniePair=1  minOverhang=4            Trim when any four fragments indicate a chimera.
//  minInniePair=0  minOverhang=0            Trim on ANY overlap coverage gap.
//
uint32   minInniePair          = 1;
uint32   minOverhang           = 1000000000;
bool     confirmWithMatedReads = false;

uint32  *libMean    = NULL;
uint32  *libStdDev  = NULL;
uint32  *libOrient  = NULL;

//  Stats for the summary file.

uint32   readsProcessed      = 0;

//  Read fates

uint32   noCoverage          = 0;
uint32   fullCoverage        = 0;
uint32   noSignalNoGap       = 0;
uint32   noSignalButGap      = 0;
uint32   gapConfirmedMate    = 0;

uint32   bothFixed           = 0;
uint32   chimeraFixed        = 0;
uint32   spurFixed           = 0;

uint32   badOvlDeleted       = 0;

uint32   bothDeletedSmall    = 0;
uint32   chimeraDeletedSmall = 0;
uint32   spurDeletedSmall    = 0;

//  Types of fixes

uint32   spurDetectedNormal = 0;
uint32   spurDetectedLinker = 0;

uint32   chimeraDetectedInnie     = 0;   //  Detected chimera
uint32   chimeraDetectedOverhang  = 0;
uint32   chimeraDetectedGap       = 0;
uint32   chimeraDetectedLinker    = 0;
uint32   chimeraDetectedGapNoMate = 0;


#define F_U32W(X)  "%" #X F_U32P
#define F_U64W(X)  "%" #X F_U64P


class chimeraClear {
public:
  uint64   deleted        : 1;
  uint64   doFixChimera   : 1;
  uint64   doCheckSubRead : 1;

  uint64   length      : AS_READ_MAX_NORMAL_LEN_BITS;

  uint64   initL       : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   initR       : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   mergL       : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   mergR       : AS_READ_MAX_NORMAL_LEN_BITS;

  uint64   tntBeg      : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   tntEnd      : AS_READ_MAX_NORMAL_LEN_BITS;

  AS_IID   libraryIID;
  AS_IID   mateIID;
};


class chimeraBad {
public:
  chimeraBad() {
    iid = 0;
    bgn = 0;
    end = 0;
  };
  chimeraBad(uint32 iid_, uint32 bgn_, uint32 end_) {
    iid = iid_;
    bgn = bgn_;
    end = end_;
  };

public:
  uint32    iid;
  uint32    bgn;
  uint32    end;
};


class chimeraOvl {
public:
  chimeraOvl() {
    flipped   = false;
    ignore    = false;
    isCrap    = false;
    style     = 0;
    Alhang    = 0;
    Abeg      = 0;
    Aend      = 0;
    Arhang    = 0;
    Blhang    = 0;
    Bbeg      = 0;
    Bend      = 0;
    Brhang    = 0;
    Aiid      = 0;
    Biid      = 0;
  };

  chimeraOvl(uint64 Aiid_, uint32 Alhang_, uint32 Abeg_, uint32 Aend_, uint32 Arhang_,
             uint64 Biid_, uint32 Blhang_, uint32 Bbeg_, uint32 Bend_, uint32 Brhang_, char ori_, bool isCrap_) {

    flipped   = (ori_ == 'r');
    ignore    = false;
    isCrap    = isCrap_;
    style     = 0;
    Alhang    = Alhang_;
    Abeg      = Abeg_;
    Aend      = Aend_;
    Arhang    = Arhang_;
    Blhang    = Blhang_;
    Bbeg      = Bbeg_;
    Bend      = Bend_;
    Brhang    = Brhang_;
    Aiid      = Aiid_;
    Biid      = Biid_;

    assert(Alhang    == Alhang_);
    assert(Abeg      == Abeg_);
    assert(Aend      == Aend_);
    assert(Arhang    == Arhang_);
    assert(Blhang    == Blhang_);
    assert(Bbeg      == Bbeg_);
    assert(Bend      == Bend_);
    assert(Brhang    == Brhang_);
    assert(Aiid      == Aiid_);
    assert(Biid      == Biid_);
    assert(isCrap    == isCrap_);

    //  We used to filter short overlaps in the adjust() function, but we want to use some of them
    //  for subread detection.  If we don't filter them out, we get invalid intervals in the chimera
    //  processing.  We could check for length at the correct time, or just set the ignore flag
    //  here.
    if ((Aend - Abeg < AS_OVERLAP_MIN_LEN) ||
        (Bend - Bbeg < AS_OVERLAP_MIN_LEN))
      ignore = true;

    if (Alhang > 0)  style |= 0x08;
    if (Arhang > 0)  style |= 0x04;
    if (Blhang > 0)  style |= 0x02;
    if (Brhang > 0)  style |= 0x01;

    switch (style) {
      case 5:
      case 7:
        //  A is chimeric anchored on the left.
      case 13:
        break;

      case 10:
      case 11:
        //  A is chimeric anchored on the right.
      case 14:
        break;

      case 6:
      case 9:
        //  Dovetail overlap
      case 1:
      case 2:
      case 3:
      case 4:
      case 8:
      case 12:
        //  Containment overlap
        if ((flipped == false) && (((Abeg >= Bbeg) && (Abeg - Bbeg < 30)) ||
                                   ((Abeg <  Bbeg) && (Bbeg - Abeg < 30)))) {
          style = 0;
        }
        break;

      case 0:
        //  Duplicate read
        break;

      case 15:
        //  Repeat overlap
        break;

      default:
        fprintf(stderr, "UNCLASSIFIED OVERLAP TYPE "F_U64"\n", style);
        exit(1);
        break;
    }
  };


  void    print(FILE *out) const {
    bool  mark = false;

    switch (style) {
      case 5:
      case 7:
      case 10:
      case 11:
      case 13:
      case 14:
        mark = true;
        break;
    }

    fprintf(out, F_U32W(6)" "F_U32W(6)" "F_U64W(2)" "F_U64W(4)" "F_U64W(4)"-"F_U64W(4)" "F_U64W(4)"  "F_U64W(4)" "F_U64W(4)"-"F_U64W(4)" "F_U64W(4)"%s\n",
            Aiid, Biid, style,
            Alhang, Abeg, Aend, Arhang,
            Blhang, Bbeg, Bend, Brhang,
            (mark) ? " ****" : "");
  };


public:
  uint64   flipped   : 1;
  uint64   ignore    : 1;
  uint64   isCrap    : 1;
  uint64   style     : 4;
  uint64   Alhang    : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Abeg      : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Aend      : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Arhang    : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Blhang    : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Bbeg      : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Bend      : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Brhang    : AS_READ_MAX_NORMAL_LEN_BITS;
  uint32   Aiid;
  uint32   Biid;
};



class chimeraRes {
public:
  chimeraRes(uint32 iid_=0, uint32 origBeg_=0, uint32 origEnd_=0) {
    iid         = iid_;

    isSpur      = false;
    isChimera   = false;
    isBadTrim   = false;
    isGood      = true;
    isBad       = false;
    deleteMe    = false;

    origBeg     = origBeg_;
    origEnd     = origEnd_;

    intervalBeg = origBeg_;
    intervalEnd = origEnd_;

    badBeg      = 0;
    badEnd      = 0;
  };

  chimeraRes(uint32 iid_,
             bool   isSpur_,
             bool   isChimera_,
             uint32 origBeg_,
             uint32 origEnd_,
             uint32 intervalBeg_,
             uint32 intervalEnd_) {
    iid         = iid_;

    isSpur      = isSpur_;
    isChimera   = isChimera_;
    isBadTrim   = false;
    isGood      = false;
    isBad       = false;
    deleteMe    = false;

    origBeg     = origBeg_;
    origEnd     = origEnd_;

    intervalBeg = intervalBeg_;
    intervalEnd = intervalEnd_;

    badBeg      = 0;
    badEnd      = 0;
  };

public:
  uint32  iid;

  bool    isSpur;
  bool    isChimera;
  bool    isBadTrim;
  bool    isGood;
  bool    isBad;
  bool    deleteMe;

  uint32  origBeg;
  uint32  origEnd;
  uint32  intervalBeg;
  uint32  intervalEnd;

  uint32  badBeg;
  uint32  badEnd;

  char    type[32];
};







chimeraClear *
readClearRanges(gkStore *gkp) {
  gkStream       *fs    = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);
  gkFragment      fr;
  chimeraClear   *clear = new chimeraClear [gkp->gkStore_getNumFragments() + 1];

  //  Enable OBTCHIMERA, this will copy the current clear range when we access it.
  //
  //  Do NOT enable TNT or the others; we expect those to be invalid unless they're set.  By
  //  enabling them, we would instead get the current clear range.
  //
  gkp->gkStore_enableClearRange(AS_READ_CLEAR_OBTCHIMERA);

  while (fs->next(&fr)) {
    AS_IID       iid  = fr.gkFragment_getReadIID();
    gkLibrary   *lr   = gkp->gkStore_getLibrary(fr.gkFragment_getLibraryIID());

    clear[iid].deleted        = fr.gkFragment_getIsDeleted() ? 1 : 0;
    clear[iid].doFixChimera   = lr->doRemoveChimericReads;
    clear[iid].doCheckSubRead = lr->doCheckForSubReads;

    clear[iid].length        = fr.gkFragment_getSequenceLength();
    clear[iid].initL         = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTINITIAL);
    clear[iid].initR         = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_OBTINITIAL);
    clear[iid].mergL         = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTMERGE);
    clear[iid].mergR         = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_OBTMERGE);

    if ((clear[iid].initL == 1) && (clear[iid].initR == 0)) {
      clear[iid].initL = 0;
      clear[iid].initR = clear[iid].length;
    }

    if ((clear[iid].mergL == 1) && (clear[iid].mergR == 0)) {
      clear[iid].mergL = 0;
      clear[iid].mergR = clear[iid].length;
    }

    clear[iid].tntBeg        = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT);
    clear[iid].tntEnd        = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT);

    clear[iid].libraryIID    = fr.gkFragment_getLibraryIID();
    clear[iid].mateIID       = fr.gkFragment_getMateIID();
  }

  delete fs;

  return(clear);
}





void
adjust(OVSoverlap         *ovl,
       uint32              ovlLen,
       const chimeraClear *clear,
       vector<chimeraOvl> &olist,
       double              errorRate,
       double              errorLimit,
       FILE               *reportFile) {

  int32  minOverlap = MAX(40, AS_OVERLAP_MIN_LEN / 2);

  for (uint32 o=0; o<ovlLen; o++) {
    int32 idA    = ovl[o].a_iid;
    int32 idB    = ovl[o].b_iid;

    if ((clear[idA].deleted) || (clear[idA].deleted))
      continue;

    //  Grab the clear ranges for each read.

    int32  ori    = (ovl[o].dat.obt.fwd) ? 'f' : 'r';

    int32 initLa = clear[idA].initL;
    int32 initRa = clear[idA].initR;
    int32 mergLa = clear[idA].mergL;
    int32 mergRa = clear[idA].mergR;

    int32 initLb = clear[idB].initL;
    int32 initRb = clear[idB].initR;
    int32 mergLb = clear[idB].mergL;
    int32 mergRb = clear[idB].mergR;

    //  OBT overlaps are relative to the OBTINITIAL clear range.  We convert them to be relative to
    //  the start of the fragment.

    int32 leftA = initLa + ovl[o].dat.obt.a_beg;
    int32 righA = initLa + ovl[o].dat.obt.a_end;
    int32 lenA  = clear[idA].length;

    int32 leftB = initLb + ovl[o].dat.obt.b_beg;
    int32 righB = initLb + ((ovl[o].dat.obt.b_end_hi << 9) | (ovl[o].dat.obt.b_end_lo));
    int32 lenB  = clear[idB].length;

    double error  = AS_OVS_decodeQuality(ovl[o].dat.obt.erate);

    if (ori == 'r') {
      int32 t = leftB;
      leftB   = righB;
      righB   = t;
    }

    //if (leftA >= righA)
    //  fprintf(stderr, "ERROR: leftA=%d !< rightA=%d\n", leftA, righA);
    //if (leftB >= righB)
    //  fprintf(stderr, "ERROR: leftB=%d !< rightB=%d\n", leftB, righB);

    if ((leftA >= righA) || (leftB >= righB))
      fprintf(stderr, "A="F_S32"\tB="F_S32"\tori=%c\tA="F_S32"-"F_S32" len "F_S32" B="F_S32"-"F_S32" len "F_S32"\t%5.3f%%\n",
              idA, idB, ori, leftA, righA, lenA, leftB, righB, lenB, error);

    assert(leftA < righA);
    assert(leftB < righB);

#ifdef REPORT_OVERLAPS
    if (reportFile)
      fprintf(reportFile, F_S32"\t"F_S32"\t%c\t"F_S32"\t"F_S32"\t"F_S32"\t"F_S32"\t"F_S32"\t"F_S32"\t%5.3f -- ",
              idA, idB, ori, leftA, righA, lenA, leftB, righB, lenB, error);
#endif

    //  Decide how much we already trimmed off.  This could be described as taking the intersection
    //  of the two OBTMERGE clear ranges for each read, with the overlap.

    int32  trimLa = (leftA  < mergLa) ? (mergLa - leftA)  : 0;
    int32  trimRa = (mergRa < righA)  ? (righA  - mergRa) : 0;
    int32  trimLb = (leftB  < mergLb) ? (mergLb - leftB)  : 0;
    int32  trimRb = (mergRb < righB)  ? (righB  - mergRb) : 0;

    assert(trimLa >= 0);
    assert(trimRa >= 0);
    assert(trimLb >= 0);
    assert(trimRb >= 0);

    //  The overlap is claiming bases a-b in one read align to bases c-d in the other.  If we have
    //  trimmed off the start of the first read (a increases) then c should increase by the same
    //  amount.  Whichever read had the most trimmed out is what we adjust by.

    int32  trimL = MAX(trimLa, trimLb);
    int32  trimR = MAX(trimRa, trimRb);

#ifdef REPORT_OVERLAPS
    if (reportFile)
      fprintf(reportFile, "CLR A %d,%d %d,%d  B %d,%d %d,%d  TRIM %d %d %d %d -- ",
              initLa, initRa, mergLa, mergRa,
              initLb, initRb, mergLb, mergRb,
              trimLa, trimRa,
              trimLb, trimRb);
#endif

    leftA += trimL;
    righA -= trimR;

    leftB += trimL;
    righB -= trimR;

#ifdef REPORT_OVERLAPS
    if (reportFile)
      fprintf(reportFile, F_S32"\t"F_S32"\t%c\t"F_S32"\t"F_S32"\t"F_S32"\t"F_S32"\t"F_S32"\t"F_S32"\t%5.3f\n",
              idA, idB, ori, leftA, righA, lenA, leftB, righB, lenB, error);
#endif

    if ((leftA < 0) ||
        (leftB < 0) ||
        (righA < 0) ||
        (righB < 0))
      //  Trimmed out.
      continue;

    if ((leftA + minOverlap > righA) ||
        (leftB + minOverlap > righB))
      //  Too small.
      continue;

    int32 lhangA = leftA  - mergLa;
    int32 rhangA = mergRa - righA;
    int32 lhangB = leftB  - mergLb;
    int32 rhangB = mergRb - righB;

    //  Third attempt.
    //
    //  The first version used hard and fast cutoffs of (>=35bp and <= 0.02 error) or (>= 70bp).
    //  These were not fair to short reads.
    //
    //  The second version used compilacted rules (see cvs), but based the decision off of the
    //  length of the A-read.  Fair, but only for assemblies of similarily sized reads.  Totally
    //  unfair for Sanger-vs-Illumin overlaps.
    //
    //  The third version sets a minimum length and identity based on the shorter fragment. These
    //  mirror closely what the second version was doing.  It was extended to allow even shorter
    //  lengths if either read is aligned to an end.

    int32  minLen = MIN(lenA, lenB);

    int32  minLenHQ = 0.10 * minLen;
    int32  minLenMQ = 0.33 * minLen;
    int32  minLenLQ = 0.50 * minLen;

    if ((lhangA == 0) ||
        (rhangA == 0) ||
        (lhangB == 0) ||
        (rhangB == 0)) {
      minLenHQ = 0.05 * minLen;
      minLenMQ = 0.15 * minLen;
      minLenLQ = 0.25 * minLen;
    }

    minLenHQ = MAX(minLenHQ, AS_OVERLAP_MIN_LEN);
    minLenMQ = MAX(minLenMQ, 100);
    minLenLQ = MAX(minLenLQ, 200);

    double maxErrHQ = 0.333 * errorRate;
    double maxErrMQ = 0.750 * errorRate;
    double maxErrLQ = 1.000 * errorRate;

    bool isNotCrap = (((righA - leftA >= minLenHQ) && (error <= maxErrHQ)) ||
                      ((righA - leftA >= minLenMQ) && (error <= maxErrMQ)) ||
                      ((righA - leftA >= minLenLQ) && (error <= maxErrLQ)));

#ifdef DEBUG_FILTER
    if (reportFile)
      fprintf(reportFile, "OVERLAP "F_IID" %3d %3d-%3d %3d vs "F_IID" %3d %3d-%3d %3d ori %c len/err %d/%f minLen/maxErr HQ:%d/%f MQ:%d/%f LQ:%d/%f%s\n",
              idA, lhangA, leftA, righA, rhangA,
              idB, lhangB, leftB, righB, rhangB, ori,
              righA - leftA, error,
              minLenHQ, maxErrHQ,
              minLenMQ, maxErrMQ,
              minLenLQ, maxErrLQ,
              (isNotCrap) ? "" : " CRAP");
#endif

    if (lhangA < 15)  lhangA = 0;
    if (rhangA < 15)  rhangA = 0;
    if (lhangB < 15)  lhangB = 0;
    if (rhangB < 15)  rhangB = 0;

    if (ori == 'r') {
      int32 t = lhangB;
      lhangB  = rhangB;
      rhangB  = t;
    }

    olist.push_back(chimeraOvl(idA, lhangA, leftA, righA, rhangA,
                               idB, lhangB, leftB, righB, rhangB, ori, !isNotCrap));
  }
}




//  Returns true if mate pairs link all regions in IL together (not chimeric).
//
//  For all mate pairs internal to the fragment, a second intervalList is built of the span
//  if the pair.  The two intervalLists are merged.
//
bool
checkSpanningMates(const AS_IID           iid,
                   const chimeraClear    *clear,
                   gkStore               *gkp,
                   vector<chimeraOvl>    &olist,
                   intervalList<int32>   &IL) {

  //  Build a list of mated reads internal to this fragment.  We have choice of algorithm:
  //    N log N -- build a set<> of the IDs for lookup
  //    N^2 / 2 -- compare all-vs-all
  //    log N   -- sort, sort, compare
  //
  //  The log N solution is to build a list of <mateID, readID, postiion>, sort that by mate id.  A
  //  second list of <readID, mateID, position>, also sorted.  The two lists are scanned, and
  //  whenever we find mateID (list 1) == readID (list 2) we have a mate internal to this fragment.

  intervalList<int32>   CC;

  //fprintf(reportFile, "checkSpanningMates()-- for read "F_IID"\n", iid);

  for (uint32 aa=0; aa<olist.size(); aa++) {
    chimeraOvl  *ovlaa   = &olist[aa];
    AS_IID       iidaa   =  ovlaa->Biid;

    for (uint32 bb=aa+1; bb<olist.size(); bb++) {
      chimeraOvl  *ovlbb   = &olist[bb];
      AS_IID       iidbb   =  ovlbb->Biid;

#if 0
      fprintf(reportFile, "A "F_IID" ("F_IID") -- B "F_IID" ("F_IID")\n",
              iidaa, clear[iidaa].mateIID,
              iidbb, clear[iidbb].mateIID);
#endif

      if (clear[iidaa].mateIID != iidbb) {
        assert(clear[iidbb].mateIID != iidaa);
        continue;
      }

      assert(clear[iidaa].mateIID == iidbb);
      assert(clear[iidbb].mateIID == iidaa);

      uint32   mean   = libMean[clear[iidaa].libraryIID];
      uint32   stddev = libStdDev[clear[iidaa].libraryIID];
      uint32   orient = libOrient[clear[iidaa].libraryIID];

      uint32   ovlaalo = ovlaa->Abeg;
      uint32   ovlaahi = ovlaa->Aend;

      uint32   ovlbblo = ovlbb->Abeg;
      uint32   ovlbbhi = ovlbb->Aend;

      assert(ovlaalo < ovlaahi);
      assert(ovlbblo < ovlbbhi);

      if ((orient == AS_READ_ORIENT_INNIE) &&
          (((ovlaalo < ovlbblo) && (ovlaa->flipped == true)) ||
           ((ovlbblo < ovlaalo) && (ovlbb->flipped == true)) ||
           (ovlaa->flipped == ovlbb->flipped)))
        continue;

      if ((orient == AS_READ_ORIENT_OUTTIE) &&
          (((ovlaalo < ovlbblo) && (ovlaa->flipped == false)) ||
           ((ovlbblo < ovlaalo) && (ovlbb->flipped == false)) ||
           (ovlaa->flipped == ovlbb->flipped)))
        continue;

      if (orient == AS_READ_ORIENT_NORMAL)
        continue;

      if (orient == AS_READ_ORIENT_ANTINORMAL)
        continue;

      uint32  bgn = MIN(ovlaalo, ovlbblo);
      uint32  end = MAX(ovlaahi, ovlbbhi);

      if ((end - bgn < mean - 3 * stddev) ||
          (mean + 3 * stddev < end - bgn))
        continue;
          
      bool contained = false;

      for (uint32 i=0; i<IL.numberOfIntervals(); i++)
        if ((IL.lo(i) < bgn) && (end < IL.hi(i))) {
#ifdef DEBUG_MATES
          fprintf(reportFile, "MATE: "F_U32","F_U32" -- "F_U32","F_U32" contained in "F_S64","F_S64"\n",
                  ovlaalo, ovlaahi,
                  ovlbblo, ovlbbhi,
                  IL.lo(i), IL.hi(i));
#endif
          contained = true;
        }

      if (contained)
        continue;

#ifdef DEBUG_MATES
      fprintf(reportFile, "MATE: "F_U32","F_U32" -- "F_U32","F_U32"\n",
              ovlaalo, ovlaahi,
              ovlbblo, ovlbbhi);
#endif

      CC.add(bgn, end - bgn);
    }
  }

  if (CC.numberOfIntervals() == 0)
    return(false);

  CC.merge();

  for (uint32 i=0; i<CC.numberOfIntervals(); i++) {
    if (CC.ct(i) < 2)
      continue;

#ifdef DEBUG_MATES
    fprintf(reportFile, "ADD: "F_S64","F_S64"\n", CC.lo(i), CC.hi(i));
#endif

    IL.add(CC.lo(i), CC.hi(i) - CC.lo(i));
  }

  IL.merge();

  if (IL.numberOfIntervals() == 1)
    return(true);
  else
    return(false);
}




chimeraRes
processChimera(const AS_IID           iid,
               const chimeraClear    *clear,
               gkStore               *gkp,
               vector<chimeraOvl>    &olist,
               FILE                  *reportFile) {
  int32   ola = clear[iid].mergL;
  int32   ora = clear[iid].mergR;

  if (olist.size() == 0)
    return(chimeraRes(iid, ola, ora));

  if (clear[iid].doFixChimera == false)
    return(chimeraRes(iid, ola, ora));

  readsProcessed++;

  uint32  loLinker = clear[iid].tntBeg;
  uint32  hiLinker = clear[iid].tntEnd;
  bool    isLinker = false;

  //  If this read has a region marked as a potential chimeric join (Illumina reads from merTrim),
  //  or as potential linker (454 mated reads from sffToCA), decide if that region is supported by
  //  overlaps.  If it isn't we need to remove overlaps from here so that it is properly detected as
  //  chimeric later on.
  //
  if (loLinker <= hiLinker) {
    uint32  isectbefore = 0;
    uint32  isect       = 0;
    uint32  isectafter  = 0;

    //  Count the number of overlaps intersecting this region, compare to the number of overlaps in
    //  the surrounding areas.
    //
    //               ---this---
    //               ----------        not isect; doesn't span more than the region
    //     --------------------------- isect
    //        ------                   isectbefore
    //   ------                        not before, too far away
    //                          -----  isectafter
    //        -------------            **
    //
    //  ** - This is trouble.  If this region is genomic dna and not linker, we expect to have more
    //  than enough true isect to notice it.  Ideally, we'd be using the overlap types (as above) to
    //  notice that this overlap has more sequence spurring off because it diagrees with the linker.
    //
    for (uint32 i=0; i<olist.size(); i++) {
      chimeraOvl  *ovl   = &olist[i];
      uint32       ovllo = ovl->Abeg;
      uint32       ovlhi = ovl->Aend;

      if ((ovl->ignore == true) ||
          (ovl->isCrap == true) ||
          (ovl->style  == 0))
        continue;

      //  Overlap ends within 5bp of the beginning of the region
      if ((ovlhi <= loLinker + OVLTRIM) && (loLinker <= ovlhi + OVLTRIM))
        isectbefore++;

      //  Overlap spans 5bp more than both ends of the region
      if ((ovllo + OVLTRIM <= loLinker) && (hiLinker + OVLTRIM <= ovlhi))
        isect++;

      //  Overlap begins within 5bp of the end of the region
      if ((ovllo <= hiLinker + OVLTRIM) && (hiLinker <= ovllo + OVLTRIM))
        isectafter++;
    }


    if (isect == 0)
      isLinker = true;

    if ((isect == 1) && ((isectbefore >= 4) || (isectafter >= 4)))
      isLinker = true;

    if ((isect == 1) && (isectbefore >= 2) && (isectafter == 0))
      isLinker = true;

    if ((isect == 1) && (isectbefore == 0) && (isectafter >= 2))
      isLinker = true;


    //  Truncate the overlaps to linker bits.
    //
    if (isLinker == true) {
#ifdef DEBUG_ISLINKER
      fprintf(reportFile, "frag "F_IID" region "F_U32"-"F_U32" isectbefore "F_U32" isect "F_U32" isectafter "F_U32"\n",
              iid,
              loLinker, hiLinker, isectbefore, isect, isectafter);
#endif

      for (uint32 i=0; i<olist.size(); i++) {
        chimeraOvl  *ovl   = &olist[i];
        uint32       ovllo =  ovl->Abeg;
        uint32       ovlhi =  ovl->Aend;

        if ((ovl->ignore == true) ||
            (ovl->isCrap == true) ||
            (ovl->style  == 0))
          continue;

        //  Overlap spans the region
        if ((ovllo <= loLinker) && (hiLinker <= ovlhi)) {
          ovl->style = 0x0f;  //  Invalid style, will be ignored
#ifdef DEBUG_ISLINKER
          fprintf(reportFile, "  overlap "F_U32"-"F_U32" --> "F_U32"-"F_U32" delete\n", ovllo, ovlhi, 0, 0);
#endif
          continue;
        }

        //  Overlap ends inside the region - trim it
        if ((ovllo < loLinker) && (loLinker < ovlhi)) {
          uint32 trim = ovlhi - loLinker;
          ovl->Aend   -= trim;
          ovl->Arhang += trim;
          ovl->Bend   -= trim;
          ovl->Brhang += trim;
          ovl->style  |= 0x04 | 0x01;
          if ((loLinker > ovlhi) || (ovl->Abeg > ovl->Aend) || (ovl->Aend - ovl->Abeg < 40)) {
            ovl->style = 0x0f;
#ifdef DEBUG_ISLINKER
            fprintf(reportFile, "  overlap "F_U32"-"F_U32" --> "F_U64"-"F_U64" begin delete\n", ovllo, ovlhi, ovl->Abeg, ovl->Aend);
          } else {
            fprintf(reportFile, "  overlap "F_U32"-"F_U32" --> "F_U64"-"F_U64" begin\n", ovllo, ovlhi, ovl->Abeg, ovl->Aend);
#endif
          }
          continue;
        }

        //  Overlap starts inside the region - trim it
        if ((ovllo < hiLinker) && (hiLinker < ovlhi)) {
          uint32 trim = hiLinker - ovllo;
          ovl->Abeg   += trim;
          ovl->Alhang += trim;
          ovl->Bbeg   += trim;
          ovl->Blhang += trim;
          ovl->style  |= 0x08 | 0x02;
          if ((ovllo > hiLinker) || (ovl->Abeg > ovl->Aend) || (ovl->Aend - ovl->Abeg < 40)) {
            ovl->style = 0x0f;
#ifdef DEBUG_ISLINKER
            fprintf(reportFile, "  overlap "F_U32"-"F_U32" --> "F_U64"-"F_U64" end delete\n", ovllo, ovlhi, ovl->Abeg, ovl->Aend);
          } else {
            fprintf(reportFile, "  overlap "F_U32"-"F_U32" --> "F_U64"-"F_U64" end\n", ovllo, ovlhi, ovl->Abeg, ovl->Aend);
#endif
          }
          continue;
        }
      }
    }
  }



  intervalList<int32>   IL;
  uint32                hasPotentialChimera = 0;
  uint32                hasInniePair        = 0;
  uint32                hasOverhang         = 0;

  //  These types are described in Bri's VI Notebook #2 page 33.

  for (uint32 i=0; i<olist.size(); i++) {
    chimeraOvl  *ovl   = &olist[i];
    bool         valid = false;
    int32        bgn   = 0;
    int32        end   = 0;
    bool         ipc   = false;

    if ((ovl->ignore == true) ||
        (ovl->isCrap == true))
      continue;

    switch (ovl->style) {
      case 5:
      case 7:
        valid = true;
        bgn   = ovl->Abeg;
        end   = ovl->Aend - OVLTRIM;
        ipc   = true;
        break;

      case 13:
        if ((ovl->Aend - ovl->Abeg) >= MIN_INTERVAL_OVERLAP) {
          valid = true;
          bgn   = ovl->Abeg + OVLTRIM;
          end   = ovl->Aend - OVLTRIM;
          ipc   = true;
        }
        break;

      case 10:
      case 11:
        valid = true;
        bgn   = ovl->Abeg + OVLTRIM;
        end   = ovl->Aend;
        ipc   = true;
        break;

      case 14:
        if ((ovl->Aend - ovl->Abeg) >= MIN_INTERVAL_OVERLAP) {
          valid = true;
          bgn   = ovl->Abeg + OVLTRIM;
          end   = ovl->Aend - OVLTRIM;
          ipc   = true;
        }
        break;

      case 6:
        valid = true;
        bgn   = ovl->Abeg;
        end   = ovl->Aend - OVLTRIM;
        break;
      case 9:
        valid = true;
        bgn   = ovl->Abeg + OVLTRIM;
        end   = ovl->Aend;
        break;

      case 1:
      case 2:
      case 3:
        valid = true;
        bgn   = ovl->Abeg;
        end   = ovl->Aend;
        break;

      case 4:
        valid = true;
        bgn   = ovl->Abeg;
        end   = ovl->Aend - OVLTRIM;
        break;
      case 8:
        valid = true;
        bgn   = ovl->Abeg + OVLTRIM;
        end   = ovl->Aend;
        break;

      case 12:
        valid = true;
        bgn   = ovl->Abeg + OVLTRIM;
        end   = ovl->Aend - OVLTRIM;
        break;

      case 0:
        break;

      case 15:
        if ((ovl->Aend - ovl->Abeg) >= MIN_INTERVAL_OVERLAP) {
          valid = true;
          bgn   = ovl->Abeg + OVLTRIM;
          end   = ovl->Aend - OVLTRIM;
        }
        break;
    }

    if (ipc)
      hasPotentialChimera++;

    if (valid) {
      IL.add(bgn, end - bgn);

#ifdef REPORT_OVERLAPS
      fprintf(reportFile, "%6d "F_U32W(6)" "F_U64W(2)" "F_U64W(4)" "F_U64W(4)"-"F_U64W(4)" "F_U64W(4)"  "F_U64W(4)" "F_U64W(4)"-"F_U64W(4)" "F_U64W(4)" interval "F_U32W(4)"-"F_U32W(4)"%s\n",
              iid,
              ovl->Biid,
              ovl->style,
              ovl->Alhang, ovl->Abeg, ovl->Aend, ovl->Arhang,
              ovl->Blhang, ovl->Bbeg, ovl->Bend, ovl->Brhang,
              bgn, end,
              (ipc) ? " ****" : "");
#endif
    }
  }


  IL.merge();

#ifdef DEBUG_INTERVAL
  for (uint32 interval=0; interval<IL.numberOfIntervals(); interval++)
    fprintf(reportFile, "interval[%d] = "F_U64"-"F_U64"\n", interval, IL.lo(interval), IL.hi(interval));
#endif


  //  Having a single style 0 overlap (no hangs on either side on
  //  either fragment) can do this.
  //
  if (IL.numberOfIntervals() == 0) {
    noCoverage++;
    return(chimeraRes(iid, ola, ora));
  }

#if 0
  if (IL.numberOfIntervals() == 1)
    fprintf(reportFile, "FULL "F_S64","F_S64" -- iid "F_IID" length "F_U64"\n",
            IL.lo(0), IL.hi(0), iid, clear[iid].length);
#endif

  if ((IL.numberOfIntervals() == 1) &&
      (IL.lo(0) == 0) &&
      (IL.hi(0) == clear[iid].length)) {
    fullCoverage++;
    return(chimeraRes(iid, ola, ora));
  }

  //  Run through the overlaps again, counting the number of innie pairs across each gap in the
  //  intervals.  Also mark hangs at the ends of intervals.
  // 
  //  The two intervalHang arrays show, for each interval, if there are bases from the other
  //  reads hanging into the gap on that side of the interval.
  //
  bool           leftIntervalHang[1025];
  bool           rightIntervalHang[1025];

  for (uint32 interval=0; interval<=IL.numberOfIntervals(); interval++) {
    int32  begGap = (interval == 0)                      ? ola : IL.hi(interval-1);
    int32  endGap = (interval == IL.numberOfIntervals()) ? ora : IL.lo(interval);

    int32  l = 0;
    int32  r = 0;

    //  initialize interval hang marks
    //
    assert(interval < 1025);

    leftIntervalHang[interval] = false;
    rightIntervalHang[interval] = false;

#ifdef DEBUG_INTERVAL
    fprintf(reportFile, "intervalHang[%d] begGap=%d endGap=%d\n", interval, begGap, endGap);
#endif

    if (begGap != endGap) {

      //  Count the number of overlaps with hangs that are
      //  on the correct side to be chimeric.
      //
      for (uint32 i=0; i<olist.size(); i++) {
        chimeraOvl  *ovl = &olist[i];

        if ((ovl->ignore == true) ||
            (ovl->isCrap == true))
          continue;

        switch (ovl->style) {
          case 5:
          case 7:
            //  These should be to the left of the endGap to count.
            if (((ovl->Aend - OVLTRIM) < endGap) && (ovl->Aend >= begGap)) {
              l++;
              assert(interval > 0);
              rightIntervalHang[interval-1] = true;
            }
            break;

          case 13:
            if ((ovl->Aend - ovl->Abeg) >= MIN_INTERVAL_OVERLAP) {
              //  These should be to the left of the endGap to count.
              if (((ovl->Aend - OVLTRIM) < endGap) && (ovl->Aend >= begGap)) {
                l++;
                assert(interval > 0);
                rightIntervalHang[interval-1] = true;
              }
            }
            break;

          case 10:
          case 11:
            //  These should be to the right of the begGap to count.
            if (((ovl->Abeg + OVLTRIM) > begGap) && (ovl->Abeg <= endGap)) {
              r++;
              assert(interval < IL.numberOfIntervals());
              leftIntervalHang[interval] = true;
            }
            break;

          case 14:
            if ((ovl->Aend - ovl->Abeg) >= MIN_INTERVAL_OVERLAP) {
              //  These should be to the right of the begGap to count.
              if (((ovl->Abeg + OVLTRIM) > begGap) && (ovl->Abeg <= endGap)) {
                r++;
                assert(interval < IL.numberOfIntervals());
                leftIntervalHang[interval] = true;
              }
            }
            break;

          case 15:
            //  Repeats.
            if ((ovl->Aend - ovl->Abeg) >= MIN_INTERVAL_OVERLAP) {
              //  These should be to the left of the endGap to count.
              if (((ovl->Aend - OVLTRIM) < endGap) && (ovl->Aend >= begGap)) {
                l++;
                assert(interval  > 0);
                rightIntervalHang[interval-1] = true;
              }
              //  These should be to the right of the begGap to count.
              if (((ovl->Abeg + OVLTRIM) > begGap) && (ovl->Abeg <= endGap)) {
                r++;
                assert(interval < IL.numberOfIntervals());
                leftIntervalHang[interval] = true;
              }
            }
            break;
        }
      }

      hasInniePair += MIN(l, r);
      hasOverhang  += l + r;
    }
  }


  bool   isChimera      = false;
  bool   isGapConfirmed = false;

  if (IL.numberOfIntervals() > 1) {

    //  Chimera induced by having linker in the middle.
    if (isLinker == true) {
      chimeraDetectedLinker++;
      isChimera = true;
    }

    //  The classic chimera pattern.
    else if ((hasPotentialChimera > 0) &&
             (hasInniePair >= minInniePair)) {
      chimeraDetectedInnie++;
      isChimera = true;
    }

    //  The aggressive chimera pattern.
    else if ((hasPotentialChimera > 0) &&
             (hasOverhang >= minOverhang)) {
      chimeraDetectedOverhang++;
      isChimera = true;
    }

    //  The super aggressive 'any gap is chimeric' pattern.
    else if ((minInniePair == 0) &&
             (minOverhang  == 0)) {
      chimeraDetectedGap++;
      isChimera = true;
    }

    else if (confirmWithMatedReads == false) {
      isChimera = false;

    } else if (checkSpanningMates(iid, clear, gkp, olist, IL) == false) {
      chimeraDetectedGapNoMate++;
      isChimera = true;

    } else {
      isChimera      = false;
      isGapConfirmed = true;
    }
  }



  //  Check if the overlaps on the left or on the right are spurs.  If so, use these min/max
  //  positions as the final clear range.
  //
  //  The general idea is to pick the overlap that reaches the end.  If that overlap indicates a
  //  spur (based on the type of overlap it is) we mark the end as containing a spur.  This is
  //  complicated by wanting to allow dovetail overlaps to reset the spur-mark, but only if they are
  //  strong dovetail overlaps compared to the spur overlaps.
  //
  //  -------------------------------------------------------------------------------------------->
  //  +212 <-----------------------------------------------
  //  +209 <-------------------------------------------------------
  //   +70 -------------------------------------------------------------->
  //  +180 <-----------------------------------------------------------------
  //  +145 <---------------------------------------------------------------------------
  //  +124 ---------------------------------------------------------------------------> +26
  //  +117 <--------------------------------------------------------------------------- +39
  //       <--------------------------------------------------------------------------- +142
  //   +64 ---------------------------------------------------------------------------> +100
  //   +100 -----------------------------------------------------------> +89
  //      +24 <------------------------------------------------------------------------ +70
  //               <------------------------------------------------------------------ +51
  //                    +1 <----------------------------------------------------------- +174
  //                               <--------------------------------------------------- +161
  //                                       <------------------------------------------- +201
  //
  //  The overlap with no left overhang is such a dovetail overlap.  We could allow this overlap to
  //  reset the spur-mark, but it's quite clear this is a spur.  We thus require dovetails to stick
  //  out past any of the spur overlaps.
  //  

#define SPUR_PADDING 5

  int32  minOvl = AS_READ_MAX_NORMAL_LEN + 1;
  int32  maxOvl = 0;
  bool    isLeftSpur  = false;
  bool    isRightSpur = false;
  bool    isSpur      = false;

  for (uint32 i=0; i<olist.size(); i++) {
    chimeraOvl  *ovl = &olist[i];

    if ((ovl->ignore == true) ||
        (ovl->isCrap == true))
      continue;

    switch (ovl->style) {
      case 5:
      case 7:
      case 13:
        if (ovl->Aend > maxOvl) {
          maxOvl = ovl->Aend;
          isRightSpur = true;
        }
        if (((isLeftSpur == true)  && (ovl->Abeg + SPUR_PADDING <= minOvl)) ||
            ((isLeftSpur == false) && (ovl->Abeg                <= minOvl))) {
          minOvl = ovl->Abeg;
          isLeftSpur = false;
        }
        break;

      case 10:
      case 11:
      case 14:
        if (ovl->Abeg + SPUR_PADDING < minOvl) {
          minOvl = ovl->Abeg;
          isLeftSpur = true;
        }
        if (((isRightSpur == true)  && (ovl->Aend - SPUR_PADDING >= maxOvl)) ||
            ((isRightSpur == false) && (ovl->Aend                >= maxOvl))) {
          maxOvl = ovl->Aend;
          isRightSpur = false;
        }
        break;

      case 6:
      case 9:
      case 1:
      case 2:
      case 3:
      case 4:
      case 8:
      case 12:
        if (((isLeftSpur == true)  && (ovl->Abeg + SPUR_PADDING <= minOvl)) ||
            ((isLeftSpur == false) && (ovl->Abeg                <= minOvl))) {
          minOvl = ovl->Abeg;
          isLeftSpur = false;
        }
        if (((isRightSpur == true)  && (ovl->Aend - SPUR_PADDING >= maxOvl)) ||
            ((isRightSpur == false) && (ovl->Aend                >= maxOvl))) {
          maxOvl = ovl->Aend;
          isRightSpur = false;
        }
        break;

      case 0:
        break;

      case 15:
        if (ovl->Abeg < minOvl) {
          minOvl = ovl->Abeg;
          isLeftSpur = true;
        }
        if (ovl->Aend > maxOvl) {
          maxOvl = ovl->Aend;
          isRightSpur = true;
        }
        break;
    }
  }

  isSpur = (isLeftSpur || isRightSpur);

  if (isSpur) {
    if (isLinker)
      spurDetectedLinker++;
    else
      spurDetectedNormal++;
  }

  if ((isChimera == false) &&
      (isSpur    == false)) {
    if      (isGapConfirmed)
      gapConfirmedMate++;

    else if (IL.numberOfIntervals() == 1)
      noSignalNoGap++;

    else
      noSignalButGap++;

    return(chimeraRes(iid, ola, ora));
  }

  //if (isSpur && gapConfirmedMate)
  //  count it how?  here?  below?  too much detail?


  //  And another pass, this time to find the largest block.  If the fragment is chimeric or a spur,
  //  we reset the clear to this interval, rather than killing the fragment completely.

  uint32  currentBeg  = ola;
  uint32  intervalBeg = 0;
  uint32  intervalEnd = 0;
  uint32  intervalMax = 0;

  bool    ignoreLast = false;

  //  A quirk is that if we are told to split on any gap, we treat every gap as having an interval
  //  hang, and so we end up picking just the largest region covered by overlaps.
  //
  if ((minInniePair == 0) &&
      (minOverhang  == 0)) {
    for (uint32 i=0; i<IL.numberOfIntervals(); i++) {
      leftIntervalHang[i] = true;
      rightIntervalHang[i] = true;
    }
  }

  //  Bad 454 reads -- usually they're too long to begin with -- tend to have only a few overlaps to
  //  the second half of the chimer, or just spurious overlaps on the end, and since the end is too
  //  long to begin with, we can't simply pick the longest interval.  Instead, we compute the
  //  average number of overlaps per interval (excluding the last interval), and ignore the last
  //  interval if it is way below average.
  //
  //  In practice this had basically no impact.  It did fix a handful of bad choices.
  //
  if (IL.numberOfIntervals() > 1) {
    uint32  aveOlaps = 0;
    for (uint32 interval=0; interval<IL.numberOfIntervals()-1; interval++)
      aveOlaps += IL.ct(interval);
    aveOlaps /= (IL.numberOfIntervals() - 1);

    if (IL.ct(IL.numberOfIntervals()-1) < 0.33 * aveOlaps)
      ignoreLast = true;
  }


  //  If there is an intervalHang covering a region with no overlaps, we do not want to include it
  //  in the result.  If the region with no overlaps is not covered by an interval hang, we will
  //  include, especially if it is at the end of the fragment.
  //
  //     [   --------^    ^-------   ]  (chimer)
  //      111111111111222223333333333
  //
  //     [   ---   --^    ^-------   ]  (chimer, with a chunk of unique sequence in the middle)
  //      111111111111222223333333333
  //
  //     [  ^-------   ]  (left spur)
  //      1112222222222
  //
  //     [   -------^  ]  (right spur)
  //      1111111111      (see below!)
  //
  for (uint32 interval=0; interval<IL.numberOfIntervals(); interval++) {
#ifdef DEBUG_INTERVAL
    fprintf(reportFile, "intervalHang[%d] %d,%d  interval "F_U64","F_U64" -- ",
            interval,
            leftIntervalHang[interval], rightIntervalHang[interval],
            IL.lo(interval), IL.hi(interval));
#endif

    //  Is this overlap region the largest?
    if (IL.hi(interval) - IL.lo(interval) > intervalMax) {
      intervalBeg = IL.lo(interval);
      intervalEnd = IL.hi(interval);
      intervalMax = intervalEnd - intervalBeg;
#ifdef DEBUG_INTERVAL
      fprintf(reportFile, "overlapregion "F_U32","F_U32" -- ", intervalBeg, intervalEnd);
#endif
    }

    //  If there is a left hanging spur, reset the begin point to the start of the current overlap
    //  region.  If not, leave it as set, including the gap before us and maybe an overlap region
    //  too.
    //
    if (leftIntervalHang[interval] == true)
      currentBeg = IL.lo(interval);

    //  Is this overlap region + whatever was saved before it, the largest?
    if (IL.hi(interval) - currentBeg > intervalMax) {
      intervalBeg = currentBeg;
      intervalEnd = IL.hi(interval);
      intervalMax = intervalEnd - intervalBeg;
#ifdef DEBUG_INTERVAL
      fprintf(reportFile, "+before "F_U32","F_U32" -- ", intervalBeg, intervalEnd);
#endif
    }

    //  If there is right hanging spur, reset the begin point to the end of the current overlap
    //  region.
    //
    if (rightIntervalHang[interval] == true)
      currentBeg = IL.hi(interval);

    //  If no right hanging spur, extend to include the gap.  This is complicated by possibly being
    //  the last gap, and possibly wanting to ignore the last gap.
    //
    //  The test for right hanging spur is actually important.  Without it, a large gap here could
    //  become the winner, even though it has no overlaps.  Something like:
    //
    //  -----------------------------------------
    //  ---------\                  /------------
    //
    //  Without the test, we'd set currentBeg to the left spur point, compute currentEnd as the
    //  right spur point, and this is the biggest interval.
    //
    if (rightIntervalHang[interval] == false) {
      int32 currentEnd = (interval == IL.numberOfIntervals()-1) ? ora : IL.lo(interval+1);

      if (ignoreLast)
        currentEnd = IL.hi(interval);

      if (currentEnd - currentBeg > intervalMax) {
        intervalBeg = currentBeg;
        intervalEnd = currentEnd;
        intervalMax = intervalEnd - intervalBeg;
#ifdef DEBUG_INTERVAL
        fprintf(reportFile, "+after %d,%d -- ", intervalBeg, intervalEnd);
#endif
      }
    }

#ifdef DEBUG_INTERVAL
    fprintf(reportFile, "%d,%d\n", intervalBeg, intervalEnd);
#endif
  }

  //  If we're not chimeric, we can pick the portion covered by overlaps.  Basically, ignore all the
  //  crud above.
  if (isChimera == false) {
    intervalBeg = minOvl;
    intervalEnd = maxOvl;
  }

  assert(intervalBeg < intervalEnd);

  return(chimeraRes(iid, isSpur, isChimera, ola, ora, intervalBeg, intervalEnd));
}





uint32
intervalOverlap(uint32 b1, uint32 e1, uint32 b2, uint32 e2) {
  uint32 minmax = MIN(e1, e2);
  uint32 maxmin = MAX(b1, b2);

  if (minmax > maxmin)
    return(minmax - maxmin);

  return(0);
}





//  Examine overlaps for a specific pattern indicating a read that flips back on itself.
//  This requires multiple overlaps from the same read, in opposite orientation, to work.
//
//  If this read isn't suspected to contain sub reads, all we can do is mark some
//  overlaps as junk.
//
//  If this read is suspected to contain sub reads, we want to mark the junction region.
//
void
processSubRead(const AS_IID           iid,
               chimeraClear          *clear,
               gkStore               *gkp,
               vector<chimeraOvl>    &olist,
               vector<chimeraBad>    &blist,
               FILE                  *subreadFile,
               bool                   subreadFileVerbose) {

  if (olist.size() == 0)
    //  No overlaps, nothing to do.
    return;

  if (clear[olist[0].Aiid].doCheckSubRead == false)
    //  Not supposed to be checking this read for subreads, nothing to do.
    return;

  //  Find the last index for each overlap.

  map<AS_IID, uint32>  secondIdx;
  map<AS_IID, uint32>  numOlaps;

  bool                 largePalindrome = false;
  intervalList<int32>  BAD;
  intervalList<int32>  BADall;

  for (uint32 ii=0; ii<olist.size(); ii++) {
    //fprintf(stderr, "COUNT OLAPS B %u %u-%u\n", olist[ii].Biid, olist[ii].Bbeg, olist[ii].Bend);
    secondIdx[olist[ii].Biid] = ii;
    numOlaps[olist[ii].Biid]++;
  }

  //  Scan overlaps.  For any pair of b_iid, with overlaps in opposite directions, compute
  //  a 'bad' interval where a suspected flip occurs.

  for (uint32 ii=0; ii<olist.size(); ii++) {

    if (numOlaps[olist[ii].Biid] == 1) {
      //  Only one overlap, can't indicate sub read!
      if ((subreadFile) && (subreadFileVerbose))
        fprintf(subreadFile, "oneOverlap                 %u (%lu-%lu) %u (%lu-%lu) -- can't indicate subreads\n",
                olist[ii].Aiid, olist[ii].Abeg, olist[ii].Aend, olist[ii].Biid, olist[ii].Bbeg, olist[ii].Bend);
      continue;
    }

    //  We should never get more than two overlaps per read pair.
    if (numOlaps[olist[ii].Biid] > 2) {
      fprintf(stderr, "WARNING: more overlaps than expected for pair %u %u.\n",
              olist[ii].Aiid, olist[ii].Biid);
      continue;
    }
    assert(numOlaps[olist[ii].Biid] == 2);

    uint32 jj = secondIdx[olist[ii].Biid];

    if (ii == jj) {
      //  Already did this one!
      if ((subreadFile) && (subreadFileVerbose))
        fprintf(subreadFile, "sameOverlap                %u (%lu-%lu) %u (%lu-%lu)\n",
                olist[ii].Aiid, olist[ii].Abeg, olist[ii].Aend, olist[ii].Biid, olist[ii].Bbeg, olist[ii].Bend);
      continue;
    }

    assert(olist[ii].Aiid == olist[jj].Aiid);
    assert(olist[ii].Biid == olist[jj].Biid);

    if (olist[ii].flipped == olist[jj].flipped) {
      fprintf(stderr, "WARNING: same orient duplicate overlaps for pair %u %u\n",
              olist[ii].Aiid, olist[ii].Biid);
      continue;
    }
    assert(olist[ii].flipped != olist[jj].flipped);


    bool  AcheckSub = (clear[olist[ii].Aiid].doCheckSubRead == true);
    bool  BcheckSub = (clear[olist[ii].Biid].doCheckSubRead == true);

    assert(AcheckSub == true);  //  Otherwise we wouldn't be in this function!

    //  Decide what type of duplicate we have.
    //    Overlap on the A read -> B read is potentially sub read containing -> don't use overlaps
    //    Overlap on the B read -> A read is potentially sub read containing -> split this read

    uint32 Aoverlap = intervalOverlap(olist[ii].Abeg, olist[ii].Aend, olist[jj].Abeg, olist[jj].Aend);
    uint32 Boverlap = intervalOverlap(olist[ii].Bbeg, olist[ii].Bend, olist[jj].Bbeg, olist[jj].Bend);

    //  If there is no overlap anywhere, we're not sure what is going on.  This could be a genomic
    //  repeat.  Leave the overlaps alone.
    //
    if ((Aoverlap == 0) &&
        (Boverlap == 0))
      continue;

    //  Remember if the overlapping ovelap is large - we'll later check if the bad region falls
    //  within here, and if there are enough spanning reads not trim.  We also use this as one more
    //  count of BAD.
    //
    if ((AcheckSub) && (Aoverlap > 1000) &&
        (BcheckSub) && (Boverlap > 1000)) {
      uint32  dist = (olist[ii].Aiid > olist[ii].Biid) ? (olist[ii].Aiid - olist[ii].Biid) : (olist[ii].Biid - olist[ii].Aiid);

      if (subreadFile)
        fprintf(subreadFile, "Palindrom ignore overlaps  %u (%lu-%lu) %u (%lu-%lu)  %u (%lu-%lu) %u (%lu-%lu)%s\n",
                olist[ii].Aiid, olist[ii].Abeg, olist[ii].Aend, olist[ii].Biid, olist[ii].Bbeg, olist[ii].Bend,
                olist[jj].Aiid, olist[jj].Abeg, olist[jj].Aend, olist[jj].Biid, olist[jj].Bbeg, olist[jj].Bend,
                (dist > 5) ? " WARNING--FAR-IID--WARNING" : "");
      largePalindrome  = true;
    }

    //  Otherwise, if the overlaps overlap on both reads by significant chunks, don't believe
    //  either.  These are possibly both chimeric reads, at least PacBio junction reads.
    //
    //  Or an inverted repeat.
    //
    if ((AcheckSub) && (Aoverlap > 50) &&
        (BcheckSub) && (Boverlap > 50)) {
      if (subreadFile)
        fprintf(subreadFile, "BothOv    ignore overlaps  %u (%lu-%lu) %u (%lu-%lu)  %u (%lu-%lu) %u (%lu-%lu)\n",
                olist[ii].Aiid, olist[ii].Abeg, olist[ii].Aend, olist[ii].Biid, olist[ii].Bbeg, olist[ii].Bend,
                olist[jj].Aiid, olist[jj].Abeg, olist[jj].Aend, olist[jj].Biid, olist[jj].Bbeg, olist[jj].Bend);
    }

    //  Stronger overlap in the A reads.  The B read looks like it has subreads, which is perfectly fine
    //  evidence for us.  Unless they span a junction.
    //
    if ((BcheckSub) && (Boverlap < Aoverlap)) {
      if (subreadFile)
        fprintf(subreadFile, "BcheckSub ignore overlaps  %u (%lu-%lu) %u (%lu-%lu)  %u (%lu-%lu) %u (%lu-%lu)\n",
                olist[ii].Aiid, olist[ii].Abeg, olist[ii].Aend, olist[ii].Biid, olist[ii].Bbeg, olist[ii].Bend,
                olist[jj].Aiid, olist[jj].Abeg, olist[jj].Aend, olist[jj].Biid, olist[jj].Bbeg, olist[jj].Bend);
    }

    //  We now want to decide if this pair of overlaps is indicating that the A read contains sub reads.
    //  Of course, if the read isn't suspected to have sub reads, there is nothing to do.
    //
    if (AcheckSub == false)
      continue;

    //  It looks like A has sub reads if the B read has a strong overlap in overlaps, and the A read does not
    //  have a strong overlap.

    if ((Aoverlap > 250) ||
        (Boverlap < 250))
      //  A strong overlap in the A read, there isn't a sub read junction we can identifiy, OR
      //  A weak overlap in the B read, and we expected the B read to align to both of the A subreads.
      continue;

    //  Don't trust these overlaps.  We _want_ to split this read, and these might just be keeping
    //  it together (this is logged later).
    //
    //olist[ii].ignore = true;
    //olist[jj].ignore = true;

    assert(olist[ii].Aiid == olist[jj].Aiid);

    //  Decide on a region in the read that is suspected to contain the chimer junction.
    //
    //  In the true case: ii overlap is first on the read; bad region from the end of this overlap
    //  to the start of the jj overlap.
    //
    //  Note that sometimes overlaps extend through the junction.  This will just flip the region
    //  around.  We're expecting to find non-overlapping overlaps, but if we find overlapping ones,
    //  the bad interval is still between the end points.
    //
    //    -------------->                 ------------>
    //                    <---------  vs             <---------
    //  
    uint32  badbgn = (olist[ii].Abeg < olist[jj].Abeg) ? olist[ii].Aend : olist[jj].Aend;
    uint32  badend = (olist[ii].Abeg < olist[jj].Abeg) ? olist[jj].Abeg : olist[ii].Abeg;

    if (badbgn > badend) {
      uint32  a = badbgn;
      badbgn = badend;
      badend = a;
    }
    assert(badbgn <= badend);

    if (subreadFile)
      fprintf(subreadFile, "AcheckSub ignore overlaps  %u (%lu-%lu) %u (%lu-%lu)  %u (%lu-%lu) %u (%lu-%lu)  BAD %u-%u size %u %s\n",
              olist[ii].Aiid, olist[ii].Abeg, olist[ii].Aend, olist[ii].Biid, olist[ii].Bbeg, olist[ii].Bend,
              olist[jj].Aiid, olist[jj].Abeg, olist[jj].Aend, olist[jj].Biid, olist[jj].Bbeg, olist[jj].Bend,
              badbgn, badend, badend - badbgn,
              (badend - badbgn <= SUBREAD_LOOP_MAX_SIZE) ? "(EVIDENCE)" : "(too far)");


    //  A true subread signature will have a small bad interval (10 bases) and largely agree on the
    //  interval.  False signature will have a large size, and not agree.  We only check for size
    //  though.
    //
    if (badend - badbgn <= SUBREAD_LOOP_MAX_SIZE)
      BAD.add(badbgn, badend - badbgn);

    //  Save all plausible pairs.
    //
    if (badend - badbgn <= SUBREAD_LOOP_EXT_SIZE)
      BADall.add(badbgn, badend - badbgn);
  }

  //
  //  Merge all the 'bad' intervals.  Save the merged intervals for later use.
  //

  BAD.merge();
  BADall.merge();

#if 0
  if (subreadFile) {
    for (uint32 bb=0; bb<BAD.numberOfIntervals(); bb++)
      fprintf(subreadFile, "BAD[%d]  %ld-%ld %d hits\n", bb, BAD.lo(bb), BAD.hi(bb), BAD.ct(bb));
    for (uint32 bb=0; bb<BADall.numberOfIntervals(); bb++)
      fprintf(subreadFile, "BADall[%d]  %ld-%ld %d hits\n", bb, BADall.lo(bb), BADall.hi(bb), BADall.ct(bb));
  }
#endif

  for (uint32 bb=0; bb<BAD.numberOfIntervals(); bb++) {
    uint32  numSpan = 0;
    uint32  allHits = 0;

    //  Find the BADall interval that corresponds to this one.  This BAD interval must be contained
    //  in a BADall (because it contains all bad intervals, while BAD is just the close stuff).
    //  Once we find it, remember the number of reads for later use.

    for (uint32 aa=0; aa<BADall.numberOfIntervals(); aa++)
      if ((BADall.lo(aa) <= BAD.lo(bb)) && (BAD.hi(bb) <= BADall.hi(aa)))
        allHits += BADall.ct(aa);

    assert(allHits != 0);

    //  Count the number of reads that span this region.  If the spanning read is not from a library
    //  that might contain subreads, give it more weight.

    for (uint32 ii=0; ii<olist.size(); ii++)
      if ((olist[ii].Abeg + 100 < BAD.lo(bb)) && (BAD.hi(bb) + 100 < olist[ii].Aend))
        numSpan += (clear[olist[ii].Aiid].doCheckSubRead) ? 1 : 2;

    if (subreadFile)
      fprintf(subreadFile, "AcheckSub region %u ("F_S32"-"F_S32") with %u hits %u bighits - span %u largePalindrome %s\n",
              olist[0].Aiid, BAD.lo(bb), BAD.hi(bb), BAD.ct(bb), allHits,
              numSpan, largePalindrome ? "true" : "false");

    if (numSpan > 9)
      //  If there are 10 or more spanning read (equivalents) this is not a subread junction.  There
      //  is plenty of evidence it is true.
      continue;

    if (BAD.ct(bb) + allHits / 4 + largePalindrome < 3)
      //  If 2 or fewer reads claim this is a sub read junction, skip it.  Evidence is weak.
      continue;

    if (subreadFile)
      fprintf(subreadFile, "CONFIRMED BAD REGION %d-%d\n", BAD.lo(bb), BAD.hi(bb));

    blist.push_back(chimeraBad(olist[0].Aiid, BAD.lo(bb), BAD.hi(bb)));
  }
}




int
main(int argc, char **argv) {

  argc = AS_configure(argc, argv);

  gkStore           *gkp                = 0L;
  OverlapStore      *ovlPrimary         = 0L;
  OverlapStore      *ovlSecondary       = 0L;

  double             errorRate          = AS_OVL_ERROR_RATE;
  double             errorLimit         = 2.5;

  bool               doUpdate           = true;
  bool               doSubreadLogging   = false;

  uint32             iidMin = 1;
  uint32             iidMax = UINT32_MAX;

  char              *outputPrefix = NULL;
  char               outputName[FILENAME_MAX];

  FILE              *summaryFile  = NULL;
  FILE              *reportFile   = NULL;
  FILE              *subreadFile  = NULL;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-G", 2) == 0) {
      gkp = new gkStore(argv[++arg], FALSE, doUpdate);
      gkp->gkStore_metadataCaching(true);

    } else if (strncmp(argv[arg], "-O", 2) == 0) {
      if (ovlPrimary == NULL)
        ovlPrimary = AS_OVS_openOverlapStore(argv[++arg]);
      else if (ovlSecondary == NULL)
        ovlSecondary = AS_OVS_openOverlapStore(argv[++arg]);
      else {
        fprintf(stderr, "Only two obtStores allowed.\n");
        err++;
      }

    } else if (strcmp(argv[arg], "-e") == 0) {
      errorRate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      errorLimit = atof(argv[++arg]);

    } else if (strncmp(argv[arg], "-mininniepair", 5) == 0) {
      minInniePair = atoi(argv[++arg]);

    } else if (strncmp(argv[arg], "-minoverhanging", 5) == 0) {
      minOverhang = atoi(argv[++arg]);

    } else if (strncmp(argv[arg], "-usemates", 5) == 0) {
      confirmWithMatedReads = true;

    } else if (strncmp(argv[arg], "-o", 2) == 0) {
      outputPrefix = argv[++arg];

    } else if (strncmp(argv[arg], "-n", 2) == 0) {
      doUpdate = false;

    } else if (strcmp(argv[arg], "-t") == 0) {
      AS_UTL_decodeRange(argv[++arg], iidMin, iidMax);

    } else if (strncmp(argv[arg], "-subreadlog", 2) == 0) {
      doSubreadLogging = true;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }

  if (errorRate < 0.0)
    err++;
  if (errorRate > AS_MAX_ERROR_RATE)
    err++;

  if ((gkp == 0L) || (ovlPrimary == 0L) || (outputPrefix == NULL) || (err)) {
    fprintf(stderr, "usage: %s -G <gkpStore> -O <obtStore> [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e erate           allow 'erate' percent error (default is AS_OVL_ERROR_RATE environment variable)\n");
    fprintf(stderr, "  -E elimit          allow 'elimit' errors\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -mininniepair n    trim if at least n pairs of innie frags detect chimer\n");
    fprintf(stderr, "  -minoverhanging n  trim if at least n frags detect chimer\n");
    fprintf(stderr, "  -usemates          trim if the read is not spanned by overlaps, and not spanned by mates\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o P               write a logging and a summary of fixes to files with prefix P\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -n                 compute and log, but don't update the store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t bgn-end         limit processing to only reads from bgn to end (inclusive)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -subreadlog        write (large) subread logging file\n");
    fprintf(stderr, "\n");

    if (errorRate < 0.0)
      fprintf(stderr, "ERROR: Error rate (-e) value %f too small; must be 'fraction error' and above 0.0\n", errorRate);

    if (errorRate > AS_MAX_ERROR_RATE)
      fprintf(stderr, "ERROR:  Error rate (-e) value %f too large; must be 'fraction error' and below %f\n", errorRate, AS_MAX_ERROR_RATE);

    exit(1);
  }

  chimeraClear    *clear = readClearRanges(gkp);



  libMean   = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  libStdDev = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  libOrient = new uint32 [gkp->gkStore_getNumLibraries() + 1];

  for (uint32 i=0; i<gkp->gkStore_getNumLibraries() + 1; i++) {
    libMean[i]   = gkp->gkStore_getLibrary(i)->mean;
    libStdDev[i] = gkp->gkStore_getLibrary(i)->stddev;
    libOrient[i] = gkp->gkStore_getLibrary(i)->orientation;
  }


  sprintf(outputName, "%s.log",         outputPrefix);
  errno = 0;
  reportFile  = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", outputName, strerror(errno)), exit(1);

  if (doSubreadLogging) {
    sprintf(outputName, "%s.subread.log", outputPrefix);
    errno = 0;
    subreadFile = fopen(outputName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", outputName, strerror(errno)), exit(1);
  }

  uint32      ovlLen = 0;
  uint32      ovlMax = 64 * 1024;
  OVSoverlap *ovl    = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * ovlMax);

  memset(ovl, 0, sizeof(OVSoverlap) * ovlMax);

  vector<chimeraOvl>  olist;
  vector<chimeraBad>  blist;

  if (iidMin < 1)
    iidMin = 1;
  if (iidMax > gkp->gkStore_getNumFragments())
    iidMax = gkp->gkStore_getNumFragments();

  fprintf(stderr, "Processing from IID "F_U32" to "F_U32" out of "F_U32" reads, using errorRate = %.2f, gkpStore will%s be updated.\n",
          iidMin,
          iidMax,
          gkp->gkStore_getNumFragments(),
          errorRate,
          doUpdate ? "" : " NOT");

  for (uint32 iid=iidMin; iid<=iidMax; iid++) {

    if (clear[iid].doFixChimera == false)
      continue;

    //  NOTE!  We DO get multiple overlaps for the same pair of fragments in the partial overlap
    //  output.  We used to pick one of the overlaps (the first seen) and ignore the rest.  We do not
    //  do that anymore.  The multiple overlaps _should_ be opposite orientation.  In PacBio, these
    //  can indicate subreads.

    loadOverlaps(iid, ovl, ovlLen, ovlMax, ovlPrimary, ovlSecondary);

    //  If there are no overlaps for this read, do nothing.

    if ((iid < ovl[0].a_iid) ||
        (ovlLen == 0))
      continue;

    //  Otherwise, trim!

    adjust(ovl, ovlLen, clear, olist, errorRate, errorLimit, reportFile);

    processSubRead(ovl[0].a_iid, clear, gkp, olist, blist, subreadFile, true);

    chimeraRes   res = processChimera(ovl[0].a_iid, clear, gkp, olist, reportFile);

    //  If after chimer trimming a read had a bad interval in the clear, just delete the read.
    //  Evidence said it was both good and bad.
    //
    //  But first, if the bad interval just touches the clear, trim it out.

    if (blist.size() > 0) {
      intervalList<int32>  goodRegions;

      for (uint32 bb=0; bb<blist.size(); bb++) {
        if ((blist[bb].end   <= res.intervalBeg) ||
            (res.intervalEnd <= blist[bb].bgn)) {
          if (subreadFile)
            fprintf(subreadFile, "BAD iid %u trim %u %u region %u %u GOOD_TRIM\n",
                    res.iid, res.intervalBeg, res.intervalEnd, blist[bb].bgn, blist[bb].end);
        }

        if ((blist[bb].bgn <= res.intervalBeg) && (res.intervalBeg <= blist[bb].end)) {
          //  Trim bad interval from the start.
          if (subreadFile)
            fprintf(subreadFile, "BAD iid %u trim %u %u region %u %u TRIM_5'\n",
                    res.iid, res.intervalBeg, res.intervalEnd, blist[bb].bgn, blist[bb].end);
          res.isBadTrim   = true;
          res.intervalBeg = MIN(blist[bb].end, res.intervalEnd);
        }

        if ((blist[bb].bgn <= res.intervalEnd) && (res.intervalEnd <= blist[bb].end)) {
          //  Trim bad interval from the end.
          if (subreadFile)
            fprintf(subreadFile, "BAD iid %u trim %u %u region %u %u TRIM_3'\n",
                    res.iid, res.intervalBeg, res.intervalEnd, blist[bb].bgn, blist[bb].end);
          res.isBadTrim   = true;
          res.intervalEnd = MAX(blist[bb].bgn, res.intervalBeg);
        }

        if ((res.intervalBeg <= blist[bb].bgn) && (blist[bb].end <= res.intervalEnd)) {
          //  Bad interval completely within the clear range.
          if (subreadFile)
            fprintf(subreadFile, "BAD iid %u trim %u %u region %u %u SUBREAD_JUNCTION\n",
                    res.iid, res.intervalBeg, res.intervalEnd, blist[bb].bgn, blist[bb].end);
          
          res.isBadTrim = true;
          res.deleteMe  = false;  //doDeleteAggressive;

#warning more than one bad interval is logged incorrectly.
          res.badBeg    = blist[bb].bgn;
          res.badEnd    = blist[bb].end;

          goodRegions.add(blist[bb].bgn, blist[bb].end - blist[bb].bgn);
        }
      }

      if (goodRegions.numberOfIntervals() > 1)
        if (subreadFile)
          fprintf(subreadFile, "WARNING: read %u has %u potential subreads; logging is inaccurate.\n", res.iid, goodRegions.numberOfIntervals() + 1);

      if (goodRegions.numberOfIntervals() > 0) {
        goodRegions.invert(res.intervalBeg, res.intervalEnd);

        uint32  goodLen = 0;
        uint32  goodBeg = 0;
        uint32  goodEnd = 0;

        for (uint32 ii=0; ii<goodRegions.numberOfIntervals(); ii++) {
          uint32  len = goodRegions.hi(ii) - goodRegions.lo(ii);

          if (subreadFile)
            fprintf(subreadFile, "BAD iid %u trim %u %u region %d %d SUBREAD_REGION len %u\n",
                    res.iid, res.intervalBeg, res.intervalEnd, goodRegions.lo(ii), goodRegions.hi(ii), len);

          if (goodLen < len) {
            goodLen = len;
            goodBeg = goodRegions.lo(ii);
            goodEnd = goodRegions.hi(ii);
          }
        }

        if (subreadFile)
          fprintf(subreadFile, "BAD iid %u trim %u %u region %u %u SUBREAD_FINAL\n",
                  res.iid, res.intervalBeg, res.intervalEnd, goodBeg, goodEnd);

        res.intervalBeg = goodBeg;
        res.intervalEnd = goodEnd;
      }
    }


    //  Decide on a solution, generate some logs and update the store.

    char    typ[256]  = {0};
    char    msg[1024] = {0};
    uint32  len       = res.intervalEnd - res.intervalBeg;

    sprintf(typ, "UNTRIMMED");

    if       (res.isSpur && res.isChimera) {
      sprintf(typ, "BOTH");

      if (len < AS_READ_MIN_LEN)
        bothDeletedSmall++;
      else
        bothFixed++;
    }

    else if  (res.isSpur) {
      sprintf(typ, "SPUR");

      if (len < AS_READ_MIN_LEN)
        spurDeletedSmall++;
      else
        spurFixed++;
    }

    else if  (res.isChimera) {
      sprintf(typ, "CHIMERA");

      if (len < AS_READ_MIN_LEN)
        chimeraDeletedSmall++;
      else
        chimeraFixed++;
    }

    //  Log any bad intervals.

    if (res.isBad) {
      sprintf(msg, "same read overlaps mark "F_U32"-"F_U32" as junction in %s", res.badBeg, res.badEnd, typ);
      sprintf(typ, "SUBREAD");

      badOvlDeleted++;
    }

    if (len < AS_READ_MIN_LEN) {
      sprintf(msg, "New length too small, fragment deleted");

      res.deleteMe = true;
    }

    //  Do the update.  If nothing changed, isGood = true.

    if ((res.isGood    == false) ||
        (res.isBadTrim == true)  ||
        (res.isBad     == true)) {
      fprintf(reportFile, F_IID" %s Trimmed from "F_U32W(4)" "F_U32W(4)" to "F_U32W(4)" "F_U32W(4)".  %s.\n",
              res.iid, typ,
              res.origBeg, res.origEnd,
              res.intervalBeg, res.intervalEnd,
              (msg[0])   ? msg       : "Length OK");

      if (doUpdate) {
        if (res.deleteMe == true) {
          gkp->gkStore_delFragment(res.iid);
        } else {
          gkFragment fr;
          gkp->gkStore_getFragment(res.iid, &fr, GKFRAGMENT_INF);
          fr.gkFragment_setClearRegion(res.intervalBeg, res.intervalEnd, AS_READ_CLEAR_OBTCHIMERA);
          gkp->gkStore_setFragment(&fr);
        }
      }
    }

    olist.clear();
    blist.clear();
  }

  delete gkp;

  //  Close log files

  if (reportFile)
    fclose(reportFile);

  if (subreadFile)
    fclose(subreadFile);

  //  Write the summary

  sprintf(outputName, "%s.summary",     outputPrefix);
  errno = 0;
  summaryFile = fopen(outputName, "w");
  if (errno) {
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", outputName, strerror(errno));
    summaryFile = stdout;
  }

  fprintf(summaryFile, "READS (= ACEEPTED + TRIMMED + DELETED)\n");
  fprintf(summaryFile, "  total processed       "F_U32"\n", readsProcessed);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "ACCEPTED\n");
  fprintf(summaryFile, "  no coverage           "F_U32"\n", noCoverage);
  fprintf(summaryFile, "  full coverage         "F_U32"\n", fullCoverage);
  fprintf(summaryFile, "  no signal, no gaps    "F_U32"\n", noSignalNoGap);
  fprintf(summaryFile, "  no signal, gaps       "F_U32"\n", noSignalButGap);
  fprintf(summaryFile, "  gap spanned by mate   "F_U32"\n", gapConfirmedMate);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "TRIMMED\n");
  fprintf(summaryFile, "  both                  "F_U32"\n", bothFixed);
  fprintf(summaryFile, "  chimera               "F_U32"\n", chimeraFixed);
  fprintf(summaryFile, "  spur                  "F_U32"\n", spurFixed);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "DELETED\n");
  fprintf(summaryFile, "  both                  "F_U32"\n", bothDeletedSmall);
  fprintf(summaryFile, "  chimera               "F_U32"\n", chimeraDeletedSmall);
  fprintf(summaryFile, "  spur                  "F_U32"\n", spurDeletedSmall);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "SPUR TYPES (= TRIMMED/DELETED spur + both)\n");
  fprintf(summaryFile, "  normal                "F_U32"\n", spurDetectedNormal);
  fprintf(summaryFile, "  linker                "F_U32"\n", spurDetectedLinker);
  fprintf(summaryFile, "\n");
  fprintf(summaryFile, "CHIMERA TYPES (= TRIMMED/DELETED chimera + both)\n");
  fprintf(summaryFile, "  innie pair            "F_U32"\n", chimeraDetectedInnie);
  fprintf(summaryFile, "  overhang              "F_U32"\n", chimeraDetectedOverhang);
  fprintf(summaryFile, "  gap                   "F_U32"\n", chimeraDetectedGap);
  fprintf(summaryFile, "  gap (no mate)         "F_U32"\n", chimeraDetectedGapNoMate);
  fprintf(summaryFile, "  linker                "F_U32"\n", chimeraDetectedLinker);

  if (summaryFile != stdout)
    fclose(summaryFile);

  delete [] clear;
  delete [] libMean;
  delete [] libStdDev;
  delete [] libOrient;

  exit(0);
}
