
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

const char *mainid = "$Id: chimera.C,v 1.39 2009-11-30 17:32:54 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "util++.H"
#include "readOverlap.H"

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"


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
uint32   minInniePair = 1;
uint32   minOverhang  = 1000000000;


//  Stats for the summary file.

uint32   readsProcessed      = 0;

uint32   chimeraFixed        = 0;
uint32   chimeraDeletedSmall = 0;
uint32   spurFixed           = 0;
uint32   spurDeletedSmall    = 0;

uint32   chimeraDetectedInnie    = 0;   //  Detected chimera
uint32   chimeraDetectedOverhang = 0;
uint32   chimeraDetectedGap      = 0;
uint32   chimeraDetectedLinker   = 0;

uint32   fullCoverage  = 0;  //  No gap in overlap coverage
uint32   gapNotChimera = 0;  //  Has a gap, but not labeled chimeric
uint32   noChimericOvl = 0;  //  Gap in coverage, but no chimeric looking overlaps



#define F_U32W(X)  "%" #X F_U32P
#define F_U64W(X)  "%" #X F_U64P


class clear_t {
public:
  uint64   deleted  :1;
  uint64   doNotOBT :1;

  uint64   length   :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   initL    :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   initR    :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   mergL    :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   mergR    :AS_READ_MAX_NORMAL_LEN_BITS;

  uint64   tntBeg   :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   tntEnd   :AS_READ_MAX_NORMAL_LEN_BITS;

  AS_UID   uid;
};


class overlap2_t {
public:
  uint64   style   :4;
  uint64   Alhang  :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Abeg    :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Aend    :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Arhang  :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Blhang  :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Bbeg    :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Bend    :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Brhang  :AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   Biid    :32;  //  Debugging only
};




clear_t *
readClearRanges(gkStore *gkp) {
  gkStream       *fs    = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);
  gkFragment      fr;
  clear_t        *clear = new clear_t [gkp->gkStore_getNumFragments() + 1];

  //  Enable OBTCHIMERA, this will copy the current clear range when we access it.
  //
  //  Do NOT enable TNT or the others; we expect those to be invalid unless they're set.  By
  //  enabling them, we would instead get the current clear range.
  //
  gkp->gkStore_enableClearRange(AS_READ_CLEAR_OBTCHIMERA);

  while (fs->next(&fr)) {
    AS_IID       iid  = fr.gkFragment_getReadIID();
    gkLibrary   *lr   = NULL;
    if (fr.gkFragment_getLibraryIID() != 0) {
       lr   = gkp->gkStore_getLibrary(fr.gkFragment_getLibraryIID());
    }

    clear[iid].deleted       = fr.gkFragment_getIsDeleted() ? 1 : 0;
    clear[iid].doNotOBT      = ((lr) && (lr->doNotOverlapTrim)) ? 1 : 0;

    clear[iid].length        = fr.gkFragment_getSequenceLength();
    clear[iid].mergL         = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTMERGE);
    clear[iid].mergR         = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_OBTMERGE);
    clear[iid].initL         = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTINITIAL);
    clear[iid].initR         = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_OBTINITIAL);

    clear[iid].tntBeg        = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT);
    clear[iid].tntEnd        = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT);

    clear[iid].uid           = fr.gkFragment_getReadUID();
  }

  delete fs;

  return(clear);
}







class overlapList {
public:
  overlapList() {
    _iid    = ~(uint64)0;
    _ovlMax = 16;
    _ovlLen = 0;
    _ovl    = new overlap2_t [_ovlMax];
  };
  ~overlapList() {
    delete [] _ovl;
  };

  void          add(uint64 Aiid, uint32 Alhang, uint32 Abeg, uint32 Aend, uint32 Arhang,
                    uint64 Biid, uint32 Blhang, uint32 Bbeg, uint32 Bend, uint32 Brhang, char ori) {
    if (_ovlLen >= _ovlMax) {
      _ovlMax *= 2;
      overlap2_t *O = new overlap2_t [_ovlMax];
      memcpy(O, _ovl, sizeof(overlap2_t) * _ovlLen);
      delete [] _ovl;
      _ovl = O;
    }

    if (_iid == ~(uint64)0)
      _iid = Aiid;
    if (_iid != Aiid)
      fprintf(stderr, "ERROR: adding "F_U64" to overlapList with iid="F_U64"\n", Aiid, _iid), exit(1);
      
    uint32   style = 0;
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
        if ((ori == 'f') && (((Abeg >= Bbeg) && (Abeg - Bbeg < 30)) ||
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
        fprintf(stderr, "UNCLASSIFIED OVERLAP TYPE "F_U32"\n", style);
        exit(1);
        break;
    }

    _ovl[_ovlLen].Biid   = Biid;
    _ovl[_ovlLen].style  = style;
    _ovl[_ovlLen].Alhang = Alhang;
    _ovl[_ovlLen].Abeg   = Abeg;
    _ovl[_ovlLen].Aend   = Aend;
    _ovl[_ovlLen].Arhang = Arhang;
    _ovl[_ovlLen].Blhang = Blhang;
    _ovl[_ovlLen].Bbeg   = Bbeg;
    _ovl[_ovlLen].Bend   = Bend;
    _ovl[_ovlLen].Brhang = Brhang;
    _ovlLen++;
  };

  void          print(FILE *out, uint32 i) const {
    bool  mark = false;

    switch (_ovl[i].style) {
      case 5:
      case 7:
      case 10:
      case 11:
      case 13:
      case 14:
        mark = true;
        break;
    }

    fprintf(out, F_U64W(6)" "F_U64W(6)" "F_U64W(2)" "F_U64W(4)" "F_U64W(4)"-"F_U64W(4)" "F_U64W(4)"  "F_U64W(4)" "F_U64W(4)"-"F_U64W(4)" "F_U64W(4)"%s\n",
            _iid, _ovl[i].Biid, _ovl[i].style,
            _ovl[i].Alhang, _ovl[i].Abeg, _ovl[i].Aend, _ovl[i].Arhang,
            _ovl[i].Blhang, _ovl[i].Bbeg, _ovl[i].Bend, _ovl[i].Brhang,
            (mark) ? " ****" : "");
  };

  uint32        length(void) const {
    return(_ovlLen);
  };
  overlap2_t    *get(uint32 i) const {
    if (i < _ovlLen)
      return(_ovl + i);
    else
      return(0L);
  };

private:
  uint64        _iid;
  overlap2_t   *_ovl;
  uint32        _ovlLen;
  uint32        _ovlMax;
};








void
printReport(char          *type,
            AS_UID         uid,
            AS_IID         iid,
            intervalList  &IL,
            uint32         intervalBeg,
            uint32         intervalEnd,
            uint32         hasPotentialChimera,
            const overlapList  *overlap) {

#ifdef WITH_REPORT_FULL
  fprintf(reportFile, "%s,"F_IID" %s!  "F_U32" intervals ("F_U32","F_U32").  "F_U32" potential chimeric overlaps (%5.2f%%).\n",
          AS_UID_toString(uid), iid, type,
          IL.numberOfIntervals(), intervalBeg, intervalEnd,
          hasPotentialChimera, (double)hasPotentialChimera / (double)overlap->length() * 100);

  for (uint32 i=0; i<overlap->length(); i++)
    overlap->print(reportFile, i);
#endif
}


void
printLogMessage(AS_UID        uid,
                AS_IID        iid,
                uint32        obtBgn,
                uint32        obtEnd,
                uint32        intervalBeg,
                uint32        intervalEnd,
                bool          doUpdate,
                char const   *type,
                char const   *message) {

  fprintf(reportFile, "%s,"F_IID" %s Trimmed from "F_U32W(4)" "F_U32W(4)" to "F_U32W(4)" "F_U32W(4)".  %s, gatekeeper store %s.\n",
          AS_UID_toString(uid), iid, type,
          obtBgn, obtEnd,
          intervalBeg, intervalEnd,
          message,
          doUpdate ? "updated" : "not updated");
}






void
process(const AS_IID           iid,
        const clear_t         *clear,
              gkStore *gkp,
              bool             doUpdate,
        const overlapList     *overlap) {

  if (overlap->length() <= 0)
    return;

  readsProcessed++;

  if ((doUpdate) && (clear[iid].doNotOBT))
    doUpdate = false;

  //fprintf(reportFile, "process %s,"F_IID"\n", AS_UID_toString(clear[iid].uid), iid);

  uint32           loLinker = clear[iid].tntBeg;
  uint32           hiLinker = clear[iid].tntEnd;
  bool             isLinker = false;


  //  If this read has left over linker from gatekeeper, we need to
  //  decide, now, if that region is supported by overlaps.  If it
  //  isn't we need to remove overlaps from here so that it is
  //  properly detected as chimeric.
  //
  if (loLinker < hiLinker) {
    uint32  isectbefore = 0;
    uint32  isect       = 0;
    uint32  isectafter  = 0;

    //  Count the number of overlaps intersecting this region, compare
    //  to the number of overlaps in the surrounding areas.
    //
    //               ---this---
    //               ----------        not isect; doesn't span more than the region
    //     --------------------------- isect
    //        ------                   isectbefore
    //   ------                        not before, too far away
    //                          -----  isectafter
    //        -------------            **
    //
    //  ** - This is trouble.  If this region is genomic dna and not
    //  linker, we expect to have more than enough true isect to
    //  notice it.  Ideally, we'd be using the overlap types (as
    //  above) to notice that this overlap has more sequence spurring
    //  off because it diagrees with the linker.
    //
    for (uint32 i=0; i<overlap->length(); i++) {
      overlap2_t  *ovl   = overlap->get(i);
      uint32       ovllo = ovl->Abeg;
      uint32       ovlhi = ovl->Aend;

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
      fprintf(reportFile, "frag %s,"F_IID" region "F_U32"-"F_U32" isectbefore "F_U32" isect "F_U32" isectafter "F_U32"\n",
              AS_UID_toString(clear[iid].uid), iid,
              loLinker, hiLinker, isectbefore, isect, isectafter);
#endif

      for (uint32 i=0; i<overlap->length(); i++) {
        overlap2_t  *ovl = overlap->get(i);
        uint32       ovllo = ovl->Abeg;
        uint32       ovlhi = ovl->Aend;

        if (ovl->style == 0)
          continue;

        //  Overlap spans the region
        if ((ovllo <= loLinker) && (hiLinker <= ovlhi)) {
          ovl->style = 16;  //  Invalid style, will be ignored
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
            ovl->style = 16;
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
            ovl->style = 16;
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

  for (uint32 i=0; i<overlap->length(); i++) {
    overlap2_t  *ovl   = overlap->get(i);
    bool         valid = false;
    int32        bgn   = 0;
    int32        end   = 0;
    bool         ipc   = false;

    switch (ovl->style) {
      case 5:
      case 7:
        //hasPotentialChimera++;
        //IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg - OVLTRIM);
        valid = true;
        bgn   = ovl->Abeg;
        end   = ovl->Aend - OVLTRIM;
        ipc   = true;
        break;

      case 13:
        if ((ovl->Aend - ovl->Abeg) >= MIN_INTERVAL_OVERLAP) {
          //hasPotentialChimera++;
          //IL.add(ovl->Abeg + OVLTRIM, ovl->Aend - ovl->Abeg - 2*OVLTRIM);
          valid = true;
          bgn   = ovl->Abeg + OVLTRIM;
          end   = ovl->Aend - OVLTRIM;
          ipc   = true;
        }
        break;

      case 10:
      case 11:
        //hasPotentialChimera++;
        //IL.add(ovl->Abeg + OVLTRIM, ovl->Aend - ovl->Abeg - OVLTRIM);
        valid = true;
        bgn   = ovl->Abeg + OVLTRIM;
        end   = ovl->Aend;
        ipc   = true;
        break;

      case 14:
        if ((ovl->Aend - ovl->Abeg) >= MIN_INTERVAL_OVERLAP) {
          //hasPotentialChimera++;
          //IL.add(ovl->Abeg + OVLTRIM, ovl->Aend - ovl->Abeg - 2*OVLTRIM);
          valid = true;
          bgn   = ovl->Abeg + OVLTRIM;
          end   = ovl->Aend - OVLTRIM;
          ipc   = true;
        }
        break;

      case 6:
        //IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg - OVLTRIM);
        valid = true;
        bgn   = ovl->Abeg;
        end   = ovl->Aend - OVLTRIM;
        break;
      case 9:
        //IL.add(ovl->Abeg + OVLTRIM, ovl->Aend - ovl->Abeg - OVLTRIM);
        valid = true;
        bgn   = ovl->Abeg + OVLTRIM;
        end   = ovl->Aend;
        break;

      case 1:
      case 2:
      case 3:
        //IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg);
        valid = true;
        bgn   = ovl->Abeg;
        end   = ovl->Aend;
        break;

      case 4:
        //IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg - OVLTRIM);
        valid = true;
        bgn   = ovl->Abeg;
        end   = ovl->Aend - OVLTRIM;
        break;
      case 8:
        //IL.add(ovl->Abeg + OVLTRIM, ovl->Aend - ovl->Abeg - OVLTRIM);
        valid = true;
        bgn   = ovl->Abeg + OVLTRIM;
        end   = ovl->Aend;
        break;

      case 12:
        //IL.add(ovl->Abeg + OVLTRIM, ovl->Aend - ovl->Abeg - 2*OVLTRIM);
        valid = true;
        bgn   = ovl->Abeg + OVLTRIM;
        end   = ovl->Aend - OVLTRIM;
        break;

      case 0:
        break;

      case 15:
        if ((ovl->Aend - ovl->Abeg) >= MIN_INTERVAL_OVERLAP) {
          //IL.add(ovl->Abeg + OVLTRIM, ovl->Aend - ovl->Abeg - 2*OVLTRIM);
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
      fprintf(reportFile, "%6d "F_U64W(6)" "F_U64W(2)" "F_U64W(4)" "F_U64W(4)"-"F_U64W(4)" "F_U64W(4)"  "F_U64W(4)" "F_U64W(4)"-"F_U64W(4)" "F_U64W(4)" interval "F_U32W(4)"-"F_U32W(4)"%s\n",
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
    fullCoverage++;
    return;
  }

  //  Run through the overlaps again, counting the number of innie pairs across each gap in the
  //  intervals.  Also mark hangs at the ends of intervals.
  // 
  //  The two intervalHang arrays show, for each interval, if there are bases from the other
  //  reads hanging into the gap on that side of the interval.
  //
  bool           leftIntervalHang[1025];
  bool           rightIntervalHang[1025];

  uint32         ola = clear[iid].mergL;
  uint32         ora = clear[iid].mergR;

  for (uint32 interval=0; interval<=IL.numberOfIntervals(); interval++) {
    uint32  begGap = (interval == 0)                      ? ola : IL.hi(interval-1);
    uint32  endGap = (interval == IL.numberOfIntervals()) ? ora : IL.lo(interval);

    uint32  l = 0;
    uint32  r = 0;

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
      for (uint32 i=0; i<overlap->length(); i++) {
        overlap2_t  *ovl = overlap->get(i);

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



  bool   isChimera = false;

  //  Chimera induced by having linker in the middle.
  if ((IL.numberOfIntervals() > 1) && (isLinker)) {
    chimeraDetectedLinker++;
    isChimera = true;
  }

  //  The classic chimera pattern.
  else if ((IL.numberOfIntervals() > 1) &&
           (hasPotentialChimera > 0) &&
           (hasInniePair >= minInniePair)) {
    chimeraDetectedInnie++;
    isChimera = true;
  }

  //  The aggressive chimera pattern.
  else if ((IL.numberOfIntervals() > 1) &&
           (hasPotentialChimera > 0) &&
           (hasOverhang >= minOverhang)) {
    chimeraDetectedOverhang++;
    isChimera = true;
  }

  //  The super aggressive 'any gap is chimeric' pattern.
  else if ((IL.numberOfIntervals() > 1) &&
           (minInniePair == 0) &&
           (minOverhang  == 0)) {
    chimeraDetectedGap++;
    isChimera = true;
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

  uint32  minOvl = AS_READ_MAX_NORMAL_LEN + 1;
  uint32  maxOvl = 0;
  bool    isLeftSpur  = false;
  bool    isRightSpur = false;
  bool    isSpur      = false;

  for (uint32 i=0; i<overlap->length(); i++) {
    overlap2_t  *ovl = overlap->get(i);

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

  isSpur = isLeftSpur || isRightSpur;




  if ((isChimera == false) &&
      (isSpur    == false)) {
    noChimericOvl++;
    return;
  }



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



  if (isSpur) {
    intervalBeg = minOvl;
    intervalEnd = maxOvl;
  }

  assert(intervalBeg < intervalEnd);





  if (isSpur) {
    if (intervalMax < AS_READ_MIN_LEN) {
      spurDeletedSmall++;
      printLogMessage(clear[iid].uid, iid, ola, ora, intervalBeg, intervalEnd, doUpdate, "SPUR", "New length too small, fragment deleted");

      if (doUpdate)
        gkp->gkStore_delFragment(iid);
    } else {
      spurFixed++;
      printLogMessage(clear[iid].uid, iid, ola, ora, intervalBeg, intervalEnd, doUpdate, "SPUR", "Length OK");

      if (doUpdate) {
        gkFragment fr;
        gkp->gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);
        fr.gkFragment_setClearRegion(intervalBeg, intervalEnd, AS_READ_CLEAR_OBTCHIMERA);
        gkp->gkStore_setFragment(&fr);
      }
    }
    printReport("SPUR", clear[iid].uid, iid, IL, intervalBeg, intervalEnd, hasPotentialChimera, overlap);

  } else if (isChimera) {
    if (intervalMax < AS_READ_MIN_LEN) {
      chimeraDeletedSmall++;
      printLogMessage(clear[iid].uid, iid, ola, ora, intervalBeg, intervalEnd, doUpdate, "CHIMERA", "New length too small, fragment deleted");

      if (doUpdate)
        gkp->gkStore_delFragment(iid);
    } else {
      chimeraFixed++;
      printLogMessage(clear[iid].uid, iid, ola, ora, intervalBeg, intervalEnd, doUpdate, "CHIMERA", "Length OK");

      if (doUpdate) {
        gkFragment fr;
        gkp->gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);
        fr.gkFragment_setClearRegion(intervalBeg, intervalEnd, AS_READ_CLEAR_OBTCHIMERA);
        gkp->gkStore_setFragment(&fr);
      }
    }
    printReport("CHIMERA", clear[iid].uid, iid, IL, intervalBeg, intervalEnd, hasPotentialChimera, overlap);

  } else {
    gapNotChimera++;
    printReport("NOT CHIMERA", clear[iid].uid, iid, IL, intervalBeg, intervalEnd, hasPotentialChimera, overlap);
  }
}



int
main(int argc, char **argv) {
  bool    doUpdate          = true;
  char   *summaryName       = 0L;
  char   *reportName        = 0L;

  gkStore           *gkp          = 0L;
  OverlapStore      *ovsprimary   = 0L;
  OverlapStore      *ovssecondary = 0L;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-gkp", 2) == 0) {
      gkp = new gkStore(argv[++arg], FALSE, doUpdate);

      //  The cache is not enabled, as we don't expect many changes to the store.
      //gkp->gkStore_metadataCaching(false);

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

    } else if (strncmp(argv[arg], "-summary", 2) == 0) {
      summaryName = argv[++arg];

    } else if (strncmp(argv[arg], "-report", 2) == 0) {
      reportName = argv[++arg];

    } else if (strncmp(argv[arg], "-test", 2) == 0) {
      doUpdate = false;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkp == 0L) || (ovsprimary == 0L) || (err)) {
    fprintf(stderr, "usage: %s [-1] -gkp <gkpStore> -ovs <ovsStore> [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -mininniepair n    trim if at least n pairs of innie frags detect chimer\n");
    fprintf(stderr, "  -minoverhanging n  trim if at least n frags detect chimer\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -summary S         write a summary of the fixes to S\n");
    fprintf(stderr, "  -report R          write a detailed report of the fixes to R\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  AS_IID           idA=0, idB=0, idAlast=0;

  clear_t         *clear = readClearRanges(gkp);
  overlapList     *overlap = new overlapList;

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

  OVSoverlap   *ovl = readOverlap(ovsprimary, ovssecondary);
  while (ovl) {
    idA    = ovl->a_iid;
    idB    = ovl->b_iid;

    //  Process any last read.

    if (idA != idAlast) {
      process(idAlast, clear, gkp, doUpdate, overlap);
      delete overlap;
      overlap = new overlapList;
    }
    idAlast = ovl->a_iid;

    //  Grab the clear ranges for each read.

    int32  ori    = (ovl->dat.obt.fwd) ? 'f' : 'r';

    int32 initLa = clear[idA].initL;
    int32 initRa = clear[idA].initR;
    int32 mergLa = clear[idA].mergL;
    int32 mergRa = clear[idA].mergR;
    int32 dela = clear[idA].deleted;

    int32 initLb = clear[idB].initL;
    int32 initRb = clear[idB].initR;
    int32 mergLb = clear[idB].mergL;
    int32 mergRb = clear[idB].mergR;
    int32 delb = clear[idA].deleted;

    //  OBT overlaps are relative to the OBTINITIAL clear range.  We convert them to be relative to
    //  the start of the fragment.

    int32 leftA = initLa + ovl->dat.obt.a_beg;
    int32 righA = initLa + ovl->dat.obt.a_end;
    int32 lenA  = clear[idA].length;

    int32 leftB = initLb + ovl->dat.obt.b_beg;
    int32 righB = initLb + ((ovl->dat.obt.b_end_hi << 9) | (ovl->dat.obt.b_end_lo));
    int32 lenB  = clear[idB].length;

    if (ori == 'r') {
      int32 t = leftB;
      leftB   = righB;
      righB   = t;
    }

    assert(leftA < righA);
    assert(leftB < righB);

    double error  = AS_OVS_decodeQuality(ovl->dat.obt.erate);

#ifdef REPORT_OVERLAPS
    if (reportFile)
      fprintf(reportFile, F_U32"\t"F_U32"\t%c\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t%5.3f -- ",
              idA, idB, ori, leftA, righA, lenA, leftB, righB, lenB, error);
#endif

    //  Decide how much we already trimmed off.  This could be described as taking the intersection
    //  of the two OBTMERGE clear ranges for each read, with the overlap.
    //
    int32  trimLa = (leftA  < mergLa) ? (mergLa - leftA)  : 0;
    int32  trimRa = (mergRa < righA)  ? (righA  - mergRa) : 0;
    int32  trimLb = (leftB  < mergLb) ? (mergLb - leftB)  : 0;
    int32  trimRb = (mergRb < righB)  ? (righB  - mergRb) : 0;

    assert(trimLa >= 0);
    assert(trimRa >= 0);
    assert(trimLb >= 0);
    assert(trimRb >= 0);

    int32  trimL = MAX(trimLa, trimLb);
    int32  trimR = MAX(trimRa, trimRb);

#ifdef REPORT_OVERLAPS
    fprintf(reportFile, "TRIM %d %d %d %d -- ", trimLa, trimRa, trimLb, trimRb);
#endif

    leftA += trimL;
    righA -= trimR;

    leftB += trimL;
    righB -= trimR;

    //  Until told otherwise, the overlap is tiny and bad.

    bool isNotCrap = false;


    if ((0 <= leftA) && (0 <= righA) &&
        (0 <= leftB) && (0 <= righB) &&
        (leftA + 40 <= righA) &&
        (leftB + 40 <= righB)) {

      //  We used to use hard and fast cutoffs here.  The original said use the overlap if its
      //  length (righA-leftA) was >=35 and error <= 0.02, or if length >= 70.  This is completely
      //  unfair to 454 reads.
      //
      //  Second idea was:
      //    length >=  40, error <= 1/3 OVL_ERATE
      //    length >= 100, error <= 1/2 OVL_ERATE
      //    length >= 200, error <= 1/1 OVL_ERATE
      //  which is even more unfair.
      //
      //  lenA, the untrimmed length of the A read.
      //
      if (lenA < 200) {
        //  454 mate read -- The smallest overlap we should be seeing from overlapper is 40bp.
        //  That's nearly half of a 454 mated read, so we allow much higher error overlaps at that
        //  length.
        //
        isNotCrap = (((righA - leftA >=  40) && (error <= 0.333 * AS_OVL_ERROR_RATE)) ||
                     ((righA - leftA >=  40) && (error <= 0.750 * AS_OVL_ERROR_RATE)) ||
                     ((righA - leftA >=  70) && (error <= 1.000 * AS_OVL_ERROR_RATE)));
      } else if (lenA < 400) {
        //  454 unmated read
        isNotCrap = (((righA - leftA >=  40) && (error <= 0.333 * AS_OVL_ERROR_RATE)) ||
                     ((righA - leftA >=  70) && (error <= 0.750 * AS_OVL_ERROR_RATE)) ||
                     ((righA - leftA >= 100) && (error <= 1.000 * AS_OVL_ERROR_RATE)));
      } else {
        //  Sanger read
        isNotCrap = (((righA - leftA >=  40) && (error <= 0.333 * AS_OVL_ERROR_RATE)) ||
                     ((righA - leftA >= 100) && (error <= 0.750 * AS_OVL_ERROR_RATE)) ||
                     ((righA - leftA >= 200) && (error <= 1.000 * AS_OVL_ERROR_RATE)));
      }
    }

    if ((isNotCrap) && (dela == false) && (delb == false)) {
      int32 lhangA = leftA  - mergLa;
      int32 rhangA = mergRa - righA;
      int32 lhangB = leftB  - mergLb;
      int32 rhangB = mergRb - righB;

      if (ori == 'r') {
        int32 t = lhangB;
        lhangB  = rhangB;
        rhangB  = t;
      }

      if (lhangA < 15)  lhangA = 0;
      if (rhangA < 15)  rhangA = 0;
      if (lhangB < 15)  lhangB = 0;
      if (rhangB < 15)  rhangB = 0;

      //  Add to the list.

#ifdef REPORT_OVERLAPS
      fprintf(reportFile, "ADD %3d %3d-%3d %3d vs %3d %3d-%3d %3d ori %c\n",
              lhangA, leftA, righA, rhangA,
              lhangB, leftB, righB, rhangB, ori);
#endif

      overlap->add(idA, lhangA, leftA, righA, rhangA,
                   idB, lhangB, leftB, righB, rhangB, ori);
    } else {
#ifdef REPORT_OVERLAPS
      fprintf(reportFile, "SKIP\n");
#endif
    }

    ovl = readOverlap(ovsprimary, ovssecondary);
  }

  process(idAlast, clear, gkp, doUpdate, overlap);
  delete overlap;

  delete gkp;

  if (summaryFile) {
    fprintf(summaryFile, "EXCEPT for the number of things fixed or deleted, these stats\n");
    fprintf(summaryFile, "are not meaningful, and likely wrong.\n");
    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "READS\n");
    fprintf(summaryFile, "total processed       "F_U32"\n", readsProcessed);
    fprintf(summaryFile, "full coverage         "F_U32"\n", fullCoverage);
    fprintf(summaryFile, "gap not chimera       "F_U32"\n", gapNotChimera);
    fprintf(summaryFile, "no chimeric ovl       "F_U32"\n", noChimericOvl);
    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "chimera fixed:        "F_U32"\n", chimeraFixed);
    fprintf(summaryFile, "chimera deleted small "F_U32"\n", chimeraDeletedSmall);
    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "spur fixed            "F_U32"\n", spurFixed);
    fprintf(summaryFile, "spur deleted small    "F_U32"\n", spurDeletedSmall);
    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "CHIMERA\n");
    fprintf(summaryFile, "from innie pair       "F_U32"\n", chimeraDetectedInnie);
    fprintf(summaryFile, "from overhang         "F_U32"\n", chimeraDetectedOverhang);
    fprintf(summaryFile, "from gap              "F_U32"\n", chimeraDetectedGap);
    fprintf(summaryFile, "from linker           "F_U32"\n", chimeraDetectedLinker);
  }

  exit(0);
}
