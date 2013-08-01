
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
#include "AS_OVS_overlapStore.H"
#include "AS_UTL_intervalList.H"

#include <vector>
#include <map>

using namespace std;


#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)


//  We trim each overlap end back by this amount.
#define OVLTRIM  ((AS_OVERLAP_MIN_LEN / 2) - 1)

//  But then use only overlaps larger than this for some of the more questionable overlaps.
#define MIN_INTERVAL_OVERLAP  60


//  WITH_REPORT_FULL will ALL overlap evidence.
//  REPORT_OVERLAPS  will print the incoming overlaps in the log.
//
#undef WITH_REPORT_FULL
#undef REPORT_OVERLAPS

#undef DEBUG_ISLINKER
#undef DEBUG_INTERVAL
#undef DEBUG_FILTER
#undef DEBUG_MATES


FILE   *summaryFile = NULL;
FILE   *reportFile  = NULL;

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
  uint64   deleted     : 1;
  uint64   doFix       : 1;
  uint64   needsFix    : 1;

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
    flipped = false;
    deleted = false;
    style   = 0;
    Alhang  = 0;
    Abeg    = 0;
    Aend    = 0;
    Arhang  = 0;
    Blhang  = 0;
    Bbeg    = 0;
    Bend    = 0;
    Brhang  = 0;
    Aiid    = 0;
    Biid    = 0;
  };

  chimeraOvl(uint64 Aiid_, uint32 Alhang_, uint32 Abeg_, uint32 Aend_, uint32 Arhang_,
             uint64 Biid_, uint32 Blhang_, uint32 Bbeg_, uint32 Bend_, uint32 Brhang_, char ori_) {

    flipped = (ori_ == 'r');
    deleted = false;
    style   = 0;
    Alhang  = Alhang_;
    Abeg    = Abeg_;
    Aend    = Aend_;
    Arhang  = Arhang_;
    Blhang  = Blhang_;
    Bbeg    = Bbeg_;
    Bend    = Bend_;
    Brhang  = Brhang_;
    Aiid    = Aiid_;
    Biid    = Biid_;

    assert(Alhang  == Alhang_);
    assert(Abeg    == Abeg_);
    assert(Aend    == Aend_);
    assert(Arhang  == Arhang_);
    assert(Blhang  == Blhang_);
    assert(Bbeg    == Bbeg_);
    assert(Bend    == Bend_);
    assert(Brhang  == Brhang_);
    assert(Aiid    == Aiid_);
    assert(Biid    == Biid_);

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
  uint64   flipped  : 1;
  uint64   deleted  : 1;
  uint64   style    : 4;
  uint64   Alhang   : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Abeg     : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Aend     : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Arhang   : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Blhang   : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Bbeg     : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Bend     : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Brhang   : AS_READ_MAX_NORMAL_LEN_BITS;
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

    clear[iid].deleted       = fr.gkFragment_getIsDeleted() ? 1 : 0;
    clear[iid].doFix         = lr->doRemoveChimericReads;
    clear[iid].needsFix      = false;

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
adjust(OVSoverlap *ovl, uint32 ovlLen, const chimeraClear *clear, vector<chimeraOvl> &olist) {

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
    fprintf(reportFile, F_S32"\t"F_S32"\t%c\t"F_S32"\t"F_S32"\t"F_S32"\t"F_S32"\t"F_S32"\t"F_S32"\t%5.3f\n",
            idA, idB, ori, leftA, righA, lenA, leftB, righB, lenB, error);
#endif

    if ((leftA < 0) ||
        (leftB < 0) ||
        (righA < 0) ||
        (righB < 0) ||
        (leftA + AS_OVERLAP_MIN_LEN > righA) ||
        (leftB + AS_OVERLAP_MIN_LEN > righB))
      //  Trimmed out.
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

    double maxErrHQ = 0.333 * AS_OVL_ERROR_RATE;
    double maxErrMQ = 0.750 * AS_OVL_ERROR_RATE;
    double maxErrLQ = 1.000 * AS_OVL_ERROR_RATE;

    bool isNotCrap = (((righA - leftA >= minLenHQ) && (error <= maxErrHQ)) ||
                      ((righA - leftA >= minLenMQ) && (error <= maxErrMQ)) ||
                      ((righA - leftA >= minLenLQ) && (error <= maxErrLQ)));

#ifdef DEBUG_FILTER
    fprintf(reportFile, "OVERLAP "F_IID" %3d %3d-%3d %3d vs "F_IID" %3d %3d-%3d %3d ori %c len/err %d/%f minLen/maxErr HQ:%d/%f MQ:%d/%f LQ:%d/%f%s\n",
            idA, lhangA, leftA, righA, rhangA,
            idB, lhangB, leftB, righB, rhangB, ori,
            righA - leftA, error,
            minLenHQ, maxErrHQ,
            minLenMQ, maxErrMQ,
            minLenLQ, maxErrLQ,
            (isNotCrap) ? "" : " CRAP");
#endif

    if (isNotCrap == false)
      continue;

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
                               idB, lhangB, leftB, righB, rhangB, ori));
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
                   intervalList          &IL) {

  //  Build a list of mated reads internal to this fragment.  We have choice of algorithm:
  //    N log N -- build a set<> of the IDs for lookup
  //    N^2 / 2 -- compare all-vs-all
  //    log N   -- sort, sort, compare
  //
  //  The log N solution is to build a list of <mateID, readID, postiion>, sort that by mate id.  A
  //  second list of <readID, mateID, position>, also sorted.  The two lists are scanned, and
  //  whenever we find mateID (list 1) == readID (list 2) we have a mate internal to this fragment.

  intervalList   CC;

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
               bool                   doUpdate,
               vector<chimeraOvl>    &olist) {
  int32   ola = clear[iid].mergL;
  int32   ora = clear[iid].mergR;

  if (olist.size() == 0)
    return(chimeraRes(iid, ola, ora));

  if (clear[iid].doFix == false)
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

      if (ovl->deleted == true)
        continue;

      if (ovl->style == 0)
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

        if (ovl->deleted == true)
          continue;

        if (ovl->style == 0)
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



  intervalList     IL;
  uint32           hasPotentialChimera = 0;
  uint32           hasInniePair        = 0;
  uint32           hasOverhang         = 0;

  //  These types are described in Bri's VI Notebook #2 page 33.

  for (uint32 i=0; i<olist.size(); i++) {
    chimeraOvl  *ovl   = &olist[i];
    bool         valid = false;
    int32        bgn   = 0;
    int32        end   = 0;
    bool         ipc   = false;

    if (ovl->deleted == true)
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

        if (ovl->deleted == true)
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

    if (ovl->deleted == true)
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
void
processSubRead(const AS_IID           iid,
               chimeraClear          *clear,
               gkStore               *gkp,
               bool                   doUpdate,
               vector<chimeraOvl>    &olist,
               vector<chimeraBad>    &blist) {

  //  Find the last index for each overlap.

  map<AS_IID, uint32>  secondIdx;
  map<AS_IID, uint32>  numOlaps;

  intervalList         BAD;

  for (uint32 ii=0; ii<olist.size(); ii++) {
    secondIdx[olist[ii].Biid] = ii;
    numOlaps[olist[ii].Biid]++;
  }

  //  Scan overlaps.  For any pair of b_iid, with overlaps in opposite directions, compute
  //  a 'bad' interval where a suspected flip occurs.

  for (uint32 ii=0; ii<olist.size(); ii++) {

    if (numOlaps[olist[ii].Biid] == 1)
      //  Only one overlap, can't indicate sub read!
      continue;

    assert(numOlaps[olist[ii].Biid] == 2);

    uint32 jj = secondIdx[olist[ii].Biid];

    if (ii == jj)
      //  Already did this one!
      continue;

    assert(olist[ii].Aiid    == olist[jj].Aiid);
    assert(olist[ii].Biid    == olist[jj].Biid);

    if (olist[ii].flipped == olist[jj].flipped) {
      fprintf(stderr, "ii %u %u %u flipped %lu\n", ii, olist[ii].Aiid, olist[ii].Biid, olist[ii].flipped);
      fprintf(stderr, "jj %u %u %u flipped %lu\n", jj, olist[jj].Aiid, olist[jj].Biid, olist[jj].flipped);
    }
    assert(olist[ii].flipped != olist[jj].flipped);

    //  Decide what type of duplicate we have.
    //    Overlap on the A read -> B read is bad, don't use overlaps
    //    Overlap on the B read -> A read needs to be split

    uint32 Aoverlap = 0;
    uint32 Boverlap = 0;

    Aoverlap = intervalOverlap(olist[ii].Abeg, olist[ii].Aend, olist[jj].Abeg, olist[jj].Aend);
    Boverlap = intervalOverlap(olist[ii].Bbeg, olist[ii].Bend, olist[jj].Bbeg, olist[jj].Bend);

    //  Stronger overlap in the A reads.  The B read looks like it has subreads and we shouldn't use
    //  this for evidence.
    //
    if (Boverlap < Aoverlap) {
      olist[ii].deleted = true;
      olist[jj].deleted = true;
      continue;
    }

    //  If the overlaps overlap on both reads by significant chunks, don't believe either.
    //  These are possibly both chimeric reads, at least PacBio junction reads.
    //
    if ((Aoverlap > 50) &&
        (Boverlap > 50)) {
      olist[ii].deleted = true;
      olist[jj].deleted = true;
      //continue;
    }

    //  If there is no overlap, we're not sure what is going on.  This could be a genomic repeat.
    //  Leave the overlaps alone.
    //
    //  Tested on e.coli, removing these overlaps results in missing some true chimeric junctions.
    //
    if (Aoverlap == 0) {
      //olist[ii].deleted = true;
      //olist[jj].deleted = true;
      continue;
    }

    //  Otherwise, the overlap is stronger on the B reads.  This is indicating that this we might
    //  have subreads in this read.

    olist[ii].deleted = true;  //  First, also ignore these overlaps.  Why?  They're probably valid.
    olist[jj].deleted = true;

    assert(olist[ii].Aiid == olist[jj].Aiid);

    clear[olist[ii].Aiid].needsFix = true;

    //  Decide on a region in the read that is suspected to contain the chimer junction.

    uint32  badbgn;
    uint32  badend;

    if (olist[ii].Abeg < olist[jj].Abeg) {
      badbgn = olist[ii].Aend;  //  ii overlap is first on the read; bad region from the end of this overlap
      badend = olist[jj].Abeg;  //  to the start of the jj overlap.
    } else {
      badbgn = olist[jj].Aend;
      badend = olist[ii].Abeg;
    }

    //  Sometimes, overlaps can extend through the junction.  This will just flip the region around.
    //
    //  We're expecting to find non-overlapping overlaps, but if we find overlapping ones, the bad interval
    //  is still between the end points.
    //
    //    -------------->                 ------------>
    //                    <---------  vs             <---------
    //  
    if (badbgn < badend)
      BAD.add(badbgn, badend - badbgn);
    else
      BAD.add(badend, badbgn - badend);

#ifdef DEBUG_BAD_REGION
    fprintf(stderr, "debug  BAD read %u %u-%u EVIDENCE  %u %u\n",
            olist[ii].Aiid, badbgn, badend, Aoverlap, Boverlap);

    fprintf(stderr, "debug  %9u %9u %c "F_U64W(5)" "F_U64W(5)"-"F_U64W(5)" "F_U64W(5)"  "F_U64W(5)" "F_U64W(5)"-"F_U64W(5)" "F_U64W(5)"\n",
            olist[ii].Aiid,
            olist[ii].Biid,
            olist[ii].flipped ? 'I' : 'N',
            olist[ii].Alhang, olist[ii].Abeg, olist[ii].Aend, olist[ii].Arhang,
            olist[ii].Blhang, olist[ii].Bbeg, olist[ii].Bend, olist[ii].Brhang);

    fprintf(stderr, "debug  %9u %9u %c "F_U64W(5)" "F_U64W(5)"-"F_U64W(5)" "F_U64W(5)"  "F_U64W(5)" "F_U64W(5)"-"F_U64W(5)" "F_U64W(5)"\n",
            olist[jj].Aiid,
            olist[jj].Biid,
            olist[jj].flipped ? 'I' : 'N',
            olist[jj].Alhang, olist[jj].Abeg, olist[jj].Aend, olist[jj].Arhang,
            olist[jj].Blhang, olist[jj].Bbeg, olist[jj].Bend, olist[jj].Brhang);
#endif
  }

  //  Merge all the 'bad' intervals.  Save the merged intervals for later use.

  BAD.merge();

  for (uint32 bb=0; bb<BAD.numberOfIntervals(); bb++) {
    //fprintf(stderr, "ADD BAD read %u FROM "F_U64"-"F_U64"\n",
    //        olist[0].Aiid, BAD.lo(bb), BAD.hi(bb));
    blist.push_back(chimeraBad(olist[0].Aiid, BAD.lo(bb), BAD.hi(bb)));
  }
}




int
main(int argc, char **argv) {
  bool               doUpdate           = true;
  char              *summaryName        = 0L;
  char              *reportName         = 0L;

  gkStore           *gkp                = 0L;
  OverlapStore      *ovsprimary         = 0L;
  OverlapStore      *ovssecondary       = 0L;

  bool               doDeleteAggressive = false;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-gkp", 2) == 0) {
      gkp = new gkStore(argv[++arg], FALSE, doUpdate);
      gkp->gkStore_metadataCaching(true);

    } else if (strncmp(argv[arg], "-ovs", 2) == 0) {
      if (ovsprimary == NULL)
        ovsprimary = AS_OVS_openOverlapStore(argv[++arg]);
      else if (ovssecondary == NULL)
        ovssecondary = AS_OVS_openOverlapStore(argv[++arg]);
      else {
        fprintf(stderr, "Only two obtStores allowed.\n");
        err++;
      }

    } else if (strncmp(argv[arg], "-mininniepair", 5) == 0) {
      minInniePair = atoi(argv[++arg]);

    } else if (strncmp(argv[arg], "-minoverhanging", 5) == 0) {
      minOverhang = atoi(argv[++arg]);

    } else if (strncmp(argv[arg], "-usemates", 5) == 0) {
      confirmWithMatedReads = true;

    } else if (strncmp(argv[arg], "-summary", 2) == 0) {
      summaryName = argv[++arg];

    } else if (strncmp(argv[arg], "-report", 2) == 0) {
      reportName = argv[++arg];

    } else if (strncmp(argv[arg], "-n", 2) == 0) {
      doUpdate = false;

    } else if (strncmp(argv[arg], "-d", 2) == 0) {
      doDeleteAggressive = true;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkp == 0L) || (ovsprimary == 0L) || (err)) {
    fprintf(stderr, "usage: %s -gkp <gkpStore> -ovs <ovsStore> [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -mininniepair n    trim if at least n pairs of innie frags detect chimer\n");
    fprintf(stderr, "  -minoverhanging n  trim if at least n frags detect chimer\n");
    fprintf(stderr, "  -usemates          trim if the read is not spanned by overlaps, and not spanned by mates\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -summary S         write a summary of the fixes to S\n");
    fprintf(stderr, "  -report R          write a detailed report of the fixes to R\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -n                 compute and log, but don't update the store\n");
    fprintf(stderr, "\n");
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


  if (summaryName) {
    errno = 0;
    summaryFile = fopen(summaryName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", summaryName, strerror(errno)), exit(1);
  }
  if (reportName) {
    errno = 0;
    reportFile  = fopen(reportName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", reportName, strerror(errno)), exit(1);
  }

  uint32      ovlMax = MAX_OVERLAPS_PER_FRAG;
  uint32      ovlLen = 0;
  OVSoverlap *ovl    = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * ovlMax);

  AS_IID      iid    = 1;
  AS_IID      max    = gkp->gkStore_getNumFragments() + 1;

  while ((clear[iid].doFix == false) && (iid < max))
    iid++;

  if (iid > 1) {
    fprintf(stderr, "Reset overlap store range to "F_IID" to "F_IID"\n", iid, max);

    AS_OVS_setRangeOverlapStore(ovsprimary,   iid, max);
    AS_OVS_setRangeOverlapStore(ovssecondary, iid, max);
  }

  ovlLen  = 0;
  ovlLen += AS_OVS_readOverlapsFromStore(ovsprimary,   ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
  ovlLen += AS_OVS_readOverlapsFromStore(ovssecondary, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);

  //  NOTE!  We DO get multiple overlaps for the same pair of fragments in the partial overlap
  //  output.  We used to pick one of the overlaps (the first seen) and ignore the rest.  We do not
  //  do that anymore.  The multiple overlaps _should_ be opposite orientation.  In PacBio, these
  //  can indicate subreads.

  vector<chimeraOvl>  olist;
  vector<chimeraBad>  blist;

  while (ovlLen > 0) {
    adjust(ovl, ovlLen, clear, olist);

    processSubRead(ovl[0].a_iid, clear, gkp, doUpdate, olist, blist);

    chimeraRes   res = processChimera(ovl[0].a_iid, clear, gkp, doUpdate, olist);

    //  If after chimer trimming a read had a bad interval in the clear, just delete the read.
    //  Evidence said it was both good and bad.
    //
    //  But first, if the bad interval just touches the clear, trim it out.

    for (uint32 bb=0; bb<blist.size(); bb++) {

      if ((blist[bb].bgn <= res.intervalBeg) && (res.intervalBeg <= blist[bb].end)) {
        //  Trim bad interval from the start.
        fprintf(stderr, "BAD iid %u trim %u %u region %u %u TRIM_5'\n",
                res.iid, res.intervalBeg, res.intervalEnd, blist[bb].bgn, blist[bb].end);
        res.isBadTrim   = true;
        res.intervalBeg = blist[bb].end;
      }

      if ((blist[bb].bgn <= res.intervalEnd) && (res.intervalEnd <= blist[bb].end)) {
        //  Trim bad interval from the end.
        fprintf(stderr, "BAD iid %u trim %u %u region %u %u TRIM_3'\n",
                res.iid, res.intervalBeg, res.intervalEnd, blist[bb].bgn, blist[bb].end);
        res.isBadTrim   = true;
        res.intervalEnd = blist[bb].bgn;
      }

      if ((res.intervalBeg <= blist[bb].bgn) && (blist[bb].end <= res.intervalEnd)) {
        //  Bad interval completely within the clear range.
        fprintf(stderr, "BAD iid %u trim %u %u region %u %u SUBREAD\n",
                res.iid, res.intervalBeg, res.intervalEnd, blist[bb].bgn, blist[bb].end);

        res.isBad     = true;
        res.deleteMe  = doDeleteAggressive;

        res.badBeg    = blist[bb].bgn;
        res.badEnd    = blist[bb].end;
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
      sprintf(msg, "same read overlaps mark "F_U32"-"F_U32" as bad in %s", res.badBeg, res.badEnd, typ);
      sprintf(typ, "BADOVL");

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
      fprintf(reportFile, F_IID" %s Trimmed from "F_U32W(4)" "F_U32W(4)" to "F_U32W(4)" "F_U32W(4)".  %s, gatekeeper store %s.\n",
              res.iid, typ,
              res.origBeg, res.origEnd,
              res.intervalBeg, res.intervalEnd,
              (msg[0])   ? msg       : "Length OK",
              (doUpdate) ? "updated" : "not updated");

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

    //  Skip past any reads that we aren't going to process

    iid = ovl[0].a_iid + 1;

    while ((clear[iid].doFix == false) && (iid < max))
      iid++;

    if (iid > ovl[0].a_iid + 1) {
      fprintf(stderr, "Reset overlap store range to "F_IID" to "F_IID"\n", iid, max);

      AS_OVS_setRangeOverlapStore(ovsprimary,   iid, max);
      AS_OVS_setRangeOverlapStore(ovssecondary, iid, max);
    }

    //  Read more overlaps.

    ovlLen  = 0;
    ovlLen += AS_OVS_readOverlapsFromStore(ovsprimary,   ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
    ovlLen += AS_OVS_readOverlapsFromStore(ovssecondary, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
  }

  delete gkp;

  if (summaryFile) {
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
  }

  delete [] clear;
  delete [] libMean;
  delete [] libStdDev;
  delete [] libOrient;

  exit(0);
}
