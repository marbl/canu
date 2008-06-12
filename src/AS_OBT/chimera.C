
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "util++.H"
#include "readOverlap.H"

extern "C" {
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"
}

//  WITH_REPORT_FULL will report unmodified fragments.
//  REPORT_OVERLAPS  will print the incoming overlaps in the log.
//
#undef  WITH_REPORT_FULL
#undef  REPORT_OVERLAPS

FILE   *summaryFile = stdout;
FILE   *reportFile  = stdout;

#define F_U32W(X)  "%" #X F_U32P
#define F_U64W(X)  "%" #X F_U64P


class clear_t {
public:
  uint64   length:11;
  uint64   origL:11;
  uint64   origR:11;
  uint64   ovlpL:11;
  uint64   ovlpR:11;
  uint64   deleted:1;
};


class overlap2_t {
public:
  uint64   style:4;
  uint64   Alhang:11;
  uint64   Abeg:11;
  uint64   Aend:11;
  uint64   Arhang:11;
  uint64   Biid:25;
  uint64   Blhang:11;
  uint64   Bbeg:11;
  uint64   Bend:11;
  uint64   Brhang:11;
};


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
      fprintf(stderr, "ERROR: adding "F_U64" to overlapList with iid="F_U64"\n", Aiid, _iid);

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
        if ((ori == 'f') && (((Abeg >= Bbeg) && (Abeg - Bbeg < 30))
                             || ((Abeg < Bbeg) && (Bbeg - Abeg < 30)))) {
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

  void          print(FILE *out, uint32 i) {
    fprintf(out,
            F_U64"\t"F_U64"\t"F_U64"\t"F_U64"\t"F_U64"\t"F_U64"\t"F_U64"\t--\t"F_U64"\t"F_U64"\t"F_U64"\t"F_U64,
            _iid, _ovl[i].Biid, _ovl[i].style,
            _ovl[i].Alhang, _ovl[i].Abeg, _ovl[i].Aend, _ovl[i].Arhang,
            _ovl[i].Blhang, _ovl[i].Bbeg, _ovl[i].Bend, _ovl[i].Brhang);
    switch (_ovl[i].style) {
      case 5:
      case 7:
      case 10:
      case 11:
      case 13:
      case 14:
        fprintf(out, "\t*****");
    }
    fprintf(out, "\n");
  };

  uint32        length(void) {
    return(_ovlLen);
  };
  overlap2_t    *get(uint32 i) {
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




uint32   chimeraDetected     = 0;
uint32   chimeraDeletedSmall = 0;
uint32   spurDetected        = 0;
uint32   spurDeletedSmall    = 0;
uint32   noInniePair         = 0;
uint32   noChimericOvl       = 0;
uint32   fullCoverage        = 0;





//  Reads both clear ranges, packs into a 64-bit integer.
//
clear_t *
readClearRanges(GateKeeperStore *gkp) {
  FragStream       *fs    = openFragStream(gkp, FRAG_S_INF);
  fragRecord        fr;
  clear_t          *clear = new clear_t [getLastElemFragStore(gkp) + 1];

  while (nextFragStream(fs, &fr)) {
    AS_IID         iid = getFragRecordIID(&fr);
    clear[iid].length  = getFragRecordSequenceLength(&fr);
    clear[iid].ovlpL   = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_OBT);
    clear[iid].ovlpR   = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_OBT);
    clear[iid].origL   = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_OBTINI);
    clear[iid].origR   = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_OBTINI);
    clear[iid].deleted = getFragRecordIsDeleted(&fr) ? 1 : 0;
  }

  closeFragStream(fs);

  return(clear);
}




void
printReport(FILE          *reportFile,
            char          *type,
            AS_UID         uid,
            AS_IID         iid,
            intervalList  &IL,
            uint32         intervalBeg,
            uint32         intervalEnd,
            uint32         hasPotentialChimera,
            overlapList  *overlap) {

  if (reportFile) {
    fprintf(reportFile, "%s,"F_IID" %s!  "F_U32" intervals ("F_U32","F_U32").  "F_U32" potential chimeric overlaps (%5.2f%%).\n",
            AS_UID_toString(uid), iid, type,
            IL.numberOfIntervals(), intervalBeg, intervalEnd,
            hasPotentialChimera, (double)hasPotentialChimera / (double)overlap->length() * 100);

    for (uint32 i=0; i<overlap->length(); i++)
      overlap->print(reportFile, i);
  }
}


void
printLogMessage(FILE         *reportFile,
                fragRecord   *fr,
                AS_UID        uid,
                AS_IID        iid,
                uint32        intervalBeg,
                uint32        intervalEnd,
                bool          doUpdate,
                char         *type,
                char         *message) {

  if (reportFile)
    fprintf(reportFile, "%s,"F_IID" %s Trimmed from "F_U32W(4)" "F_U32W(4)" to "F_U32W(4)" "F_U32W(4)".  %s, gatekeeper store %s.\n",
            AS_UID_toString(uid),
            iid,
            type,
            getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_OBT),
            getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_OBT),
            intervalBeg,
            intervalEnd,
            message,
            doUpdate ? "updated" : "not updated");
}





//    Sort each overlapList?
//    Build an interval list using all overlaps.
//    Squash the list, see if it intersects.
//    Look for chimeric and spur patterns.
//
void
process(uint32           iid,
        GateKeeperStore *gkp,
        bool             doUpdate,
        overlapList     *overlap,
        uint32           ola,
        uint32           ora) {

  if (overlap->length() <= 0)
    return;

  int    slopSm = 20;  //  A little slop

  intervalList   IL;
  uint32         hasPotentialChimera = 0;
  uint32         hasInniePair        = 0;

  for (uint32 i=0; i<overlap->length(); i++) {
    overlap2_t  *ovl = overlap->get(i);

    switch (ovl->style) {
      case 5:
      case 7:
        hasPotentialChimera++;
        IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg - slopSm);
        break;

      case 13:
        if ((ovl->Aend - ovl->Abeg) > 75) {
          hasPotentialChimera++;
          IL.add(ovl->Abeg + slopSm, ovl->Aend - ovl->Abeg - 2*slopSm);
        }
        break;
          
      case 10:
      case 11:
        hasPotentialChimera++;
        IL.add(ovl->Abeg + slopSm, ovl->Aend - ovl->Abeg - slopSm);
        break;

      case 14:
        if ((ovl->Aend - ovl->Abeg) > 75) {
          hasPotentialChimera++;
          IL.add(ovl->Abeg + slopSm, ovl->Aend - ovl->Abeg - 2*slopSm);
        }
        break;

      case 6:
        IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg - slopSm);
        break;
      case 9:
        IL.add(ovl->Abeg + slopSm, ovl->Aend - ovl->Abeg - slopSm);
        break;

      case 1:
      case 2:
      case 3:
        IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg);
        break;          

      case 4:
        IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg - slopSm);
        break;
      case 8:
        IL.add(ovl->Abeg + slopSm, ovl->Aend - ovl->Abeg - slopSm);
        break;

      case 12:
        IL.add(ovl->Abeg + slopSm, ovl->Aend - ovl->Abeg - 2*slopSm);
        break;          

      case 0:
        break;

      case 15:
        if ((ovl->Aend - ovl->Abeg) > 75)
          IL.add(ovl->Abeg + slopSm, ovl->Aend - ovl->Abeg - 2*slopSm);
        break;

      default:
        fprintf(stderr, "UNCLASSIFIED OVERLAP TYPE "F_U64"\n", ovl->style);
        break;
    }
  }

  IL.merge();


  //  Having a single style 0 overlap (no hangs on either side on
  //  either fragment) can do this.
  //
  if (IL.numberOfIntervals() == 0)
    return;


  bool           leftIntervalHang[1025];
  bool           rightIntervalHang[1025];
  bool           isSpur = false;


  //  Run through the overlaps again, counting the number of innie
  //  pairs across each gap in the intervals.
  //  Also mark hangs at the ends of intervals.
  //
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
            if (((ovl->Aend - slopSm) < endGap) && (ovl->Aend >= begGap)) {
              l++;
              assert(interval > 0);
              rightIntervalHang[interval-1] = true;
            }
            break;
          
          case 13:
            if ((ovl->Aend - ovl->Abeg) > 75) {
              //  These should be to the left of the endGap to count.
              if (((ovl->Aend - slopSm) < endGap) && (ovl->Aend >= begGap)) {
                l++;
                assert(interval > 0);
                rightIntervalHang[interval-1] = true;
              }
            }
            break;
          
          case 10:
          case 11:
            //  These should be to the right of the begGap to count.
            if (((ovl->Abeg + slopSm) > begGap) && (ovl->Abeg <= endGap)) {
              r++;
              assert(interval < IL.numberOfIntervals());
              leftIntervalHang[interval] = true;
            }
            break;
	  
          case 14:
            if ((ovl->Aend - ovl->Abeg) > 75) {
              //  These should be to the right of the begGap to count.
              if (((ovl->Abeg + slopSm) > begGap) && (ovl->Abeg <= endGap)) {
                r++;
                assert(interval < IL.numberOfIntervals());
                leftIntervalHang[interval] = true;
              }
            }
            break;
	  
          case 15:
            //  Repeats.
            if ((ovl->Aend - ovl->Abeg) > 75) {
              //  These should be to the left of the endGap to count.
              if (((ovl->Aend - slopSm) < endGap) && (ovl->Aend >= begGap)) {
                l++;
                assert(interval  > 0);
                rightIntervalHang[interval-1] = true;
              }
              //  These should be to the right of the begGap to count.
              if (((ovl->Abeg + slopSm) > begGap) && (ovl->Abeg <= endGap)) {
                r++;
                assert(interval < IL.numberOfIntervals());
                leftIntervalHang[interval] = true;
              }
            }
            break;
        }
      }

      if ((l > 0) && (r > 0))
        hasInniePair++;
    }
  }

  {
    uint32  minOvl = AS_FRAG_MAX_LEN + 1;
    uint32  maxOvl = 0;
    bool    isLeftSpur = false;
    bool    isRightSpur = false;

    //  Check if the overlaps on the left or on the right are spurs.

    for (uint32 i=0; i<overlap->length(); i++) {
      overlap2_t  *ovl = overlap->get(i);

      switch (ovl->style) {
        case 5:
        case 7:
        case 13:
          //  These should be at the end (right) to count as spurs.
          if (ovl->Aend > maxOvl) {
            maxOvl = ovl->Aend;
            isRightSpur = true;
          }
          if (ovl->Abeg <= minOvl) {
            minOvl = ovl->Abeg;
            isLeftSpur = false;
          }
          break;
          
        case 10:
        case 11:
        case 14:
          //  These should be at the beginning (left) to count as spurs.
          if (ovl->Abeg < minOvl) {
            minOvl = ovl->Abeg;
            isLeftSpur = true;
          }
          if (ovl->Aend >= maxOvl) {
            maxOvl = ovl->Aend;
            isRightSpur = false;
          }
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
          if (ovl->Abeg <= minOvl) {
            minOvl = ovl->Abeg;
            isLeftSpur = false;
          }
          if (ovl->Aend >= maxOvl) {
            maxOvl = ovl->Aend;
            isRightSpur = false;
          }
          break;          
      
        case 0:
          //  Duplicate read
          break;
      
        case 15:
          //  Repeat overlap
          if (ovl->Abeg < minOvl) {
            minOvl = ovl->Abeg;
            isLeftSpur = true;
          }
          if (ovl->Aend > maxOvl) {
            maxOvl = ovl->Aend;
            isRightSpur = true;
          }
          break;
      
        default:
          fprintf(stderr, "UNCLASSIFIED OVERLAP TYPE "F_U64"\n", ovl->style);
          break;
      }
    }

    isSpur = isLeftSpur || isRightSpur;
  }


  //  And another pass, this time to find the largest block.  If the
  //  fragment is chimeric or a spur, we reset the clear to this
  //  interval, rather than killing the fragment completely.
  //
  //  In pictures, using ^ to indicate the IntervalHang[]s, and #'s to
  //  indicate the regions that could be saved.
  //
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
  //

  uint32  intervalBeg = 0;
  uint32  intervalEnd = 0;
  uint32  intervalMax = 0;
  {
    uint32  currentBeg = ola;
    bool    ignoreLast = false;

    if (IL.numberOfIntervals() > 1) {
      //  Bad 454 reads -- usually they're too long to begin with --
      //  tend to have only a few overlaps to the second half of the
      //  chimer, or just spurious overlaps on the end, and since the
      //  end is too long to begin with, we can't simply pick the
      //  longest interval.  Instead, we compute the average number of
      //  overlaps per interval (excluding the last interval), and
      //  ignore the last interval if it is way below average.
      //
      //  In practice this had basically no impact.  It did fix a
      //  handful of bad choices.

      uint32  aveOlaps = 0;
      for (uint32 interval=0; interval<IL.numberOfIntervals()-1; interval++)
        aveOlaps += IL.ct(interval);
      aveOlaps /= (IL.numberOfIntervals() - 1);

      if (IL.ct(IL.numberOfIntervals()-1) < 0.33 * aveOlaps)
        ignoreLast = true;
    }


    for (uint32 interval=0; interval<IL.numberOfIntervals(); interval++) {

      if (leftIntervalHang[interval]) {
        //  If this interval (currentBeg <-> IL.lo() - slopSm) is the biggest, save it.
        if (IL.lo(interval) > intervalMax + slopSm + currentBeg) {
	  intervalBeg = currentBeg;
	  intervalEnd = IL.lo(interval) - slopSm;
	  intervalMax = intervalEnd - intervalBeg;
	}
        //  The next interval can begin where this one stops.
	currentBeg = IL.lo(interval);
      }

      if (rightIntervalHang[interval]) {
        //  If this interval (currentBeg <-> IL.hi()) is the biggest, save it.
        if (IL.hi(interval) > intervalMax + currentBeg) {
	  intervalBeg = currentBeg;
	  intervalEnd = IL.hi(interval);
	  intervalMax = intervalEnd - intervalBeg;
	}
        //  The next interval can begin where this one stops.
	currentBeg = IL.hi(interval) + slopSm;
      }
    }

    //  Check the last interval.  Why?  'currentBeg' could start
    //  before the last interval.  But it only does this if the last
    //  interval has a rightIntervalHang[].  We explicitly disallow
    //  the "right spur" case from getting here.
    //
    //  Finally, we only count the overlapped length -- not any
    //  unmapped crud at the end.  (used to be ora).  Again, this had
    //  a minor positive effect.
    //
    if ((rightIntervalHang[IL.numberOfIntervals()-1] == false) &&
        (ignoreLast == false) &&
        (IL.hi(IL.numberOfIntervals()-1) - currentBeg > intervalMax)) {
      intervalBeg = currentBeg;
      intervalEnd = ora;
      intervalMax = intervalEnd - intervalBeg;
    }
  }


  fragRecord        fr;

  getFrag(gkp, iid, &fr, FRAG_S_INF);

  AS_UID uid = getFragRecordUID(&fr);

  GateKeeperLibraryRecord  *gklr = getGateKeeperLibrary(gkp, getFragRecordLibraryIID(&fr));

  if ((doUpdate) && (gklr) && (gklr->doNotOverlapTrim))
    doUpdate = false;

  if (isSpur) {
    spurDetected++;

    if (intervalMax < AS_READ_MIN_LEN) {
      spurDeletedSmall++;

      printLogMessage(reportFile, &fr, uid, iid, intervalBeg, intervalEnd, doUpdate, "SPUR", "New length too small, fragment deleted");

      if (doUpdate)
        delFrag(gkp, iid);
    } else {
      printLogMessage(reportFile, &fr, uid, iid, intervalBeg, intervalEnd, doUpdate, "SPUR", "Length OK");

      if (doUpdate) {
        setFragRecordClearRegion(&fr, intervalBeg, intervalEnd, AS_READ_CLEAR_OBT);
        setFrag(gkp, iid, &fr);
      }
    }

#ifdef WITH_FULL_REPORT
    printReport(reportFile, "SPUR", uid, iid, IL, intervalBeg, intervalEnd, hasPotentialChimera, overlap);
#endif
  } else if ((IL.numberOfIntervals() > 1) &&
             (hasPotentialChimera > 0) &&
             (hasInniePair > 0)) {
    chimeraDetected++;

    if (intervalMax < AS_READ_MIN_LEN) {
      chimeraDeletedSmall++;

      printLogMessage(reportFile, &fr, uid, iid, intervalBeg, intervalEnd, doUpdate, "CHIMERA", "New length too small, fragment deleted");

      if (doUpdate)
        delFrag(gkp, iid);
    } else {
      printLogMessage(reportFile, &fr, uid, iid, intervalBeg, intervalEnd, doUpdate, "CHIMERA", "Length OK");

      if (doUpdate) {
        setFragRecordClearRegion(&fr, intervalBeg, intervalEnd, AS_READ_CLEAR_OBT);
        setFrag(gkp, iid, &fr);
      }
    }

#ifdef WITH_REPORT_FULL
    printReport(reportFile, "CHIMERA", uid, iid, IL, intervalBeg, intervalEnd, hasPotentialChimera, overlap);
#endif
  } else if (IL.numberOfIntervals() == 1) {
    fullCoverage++;

#ifdef WITH_REPORT_FULL
    printReport(reportFile, "FULL COVERAGE", uid, iid, IL, intervalBeg, intervalEnd, hasPotentialChimera, overlap);
#endif
  } else if (hasPotentialChimera == 0) {
    noChimericOvl++;

#ifdef WITH_REPORT_FULL
    printReport(reportFile, "NO CHIMERIC OVERLAPS", uid, iid, IL, intervalBeg, intervalEnd, hasPotentialChimera, overlap);
#endif
  } else if (hasInniePair == 0) {
    noInniePair++;

#ifdef WITH_REPORT_FULL
    printReport(reportFile, "NO INNIE PAIR", uid, iid, IL, intervalBeg, intervalEnd, hasPotentialChimera, overlap);
#endif
  } else {
#ifdef WITH_REPORT_FULL
    printReport(reportFile, "NOT CHIMERA", uid, iid, IL, intervalBeg, intervalEnd, hasPotentialChimera, overlap);
#endif
  }

}



int
main(int argc, char **argv) {
  bool    doUpdate          = true;  //  set to false for testing
  char   *summaryName       = 0L;
  char   *reportName        = 0L;

  GateKeeperStore   *gkp          = 0L;
  OverlapStore      *ovsprimary   = 0L;
  OverlapStore      *ovssecondary = 0L;

  uint32  overflow = 0;
  uint32  notclear = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-gkp", 2) == 0) {
      gkp      = openGateKeeperStore(argv[++arg], doUpdate);
      gkp->frg = convertStoreToMemoryStore(gkp->frg);
    } else if (strncmp(argv[arg], "-ovs", 2) == 0) {
      if (ovsprimary == NULL)
        ovsprimary = AS_OVS_openOverlapStore(argv[++arg]);
      else if (ovssecondary == NULL)
        ovssecondary = AS_OVS_openOverlapStore(argv[++arg]);
      else {
        fprintf(stderr, "Only two obtStores allowed.\n");
        err++;
      }
    } else if (strncmp(argv[arg], "-summary", 2) == 0) {
      summaryName = argv[++arg];
    } else if (strncmp(argv[arg], "-report", 2) == 0) {
      reportName = argv[++arg];
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkp == 0L) || (ovsprimary == 0L) || (err)) {
    fprintf(stderr, "usage: %s [-1] -gkp <gkpStore> -ovs <ovsStore> [opts]\n", argv[0]);
    fprintf(stderr, "  -delete         instead of fixing chimera, delete them\n");
    fprintf(stderr, "  -summary S      write a summary of the fixes to S\n");
    fprintf(stderr, "  -report R       write a detailed report of the fixes to R\n");
    exit(1);
  }

  uint32           idAlast, olalast, oralast;

  uint32           idA, idB;
  char             ori;
  uint32           leftA, rightA, lenA;
  uint32           leftB, rightB, lenB;
  double           error;
  clear_t         *clear = readClearRanges(gkp);
  uint64           maxIID      = 65536;
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
    ori    = (ovl->dat.obt.fwd) ? 'f' : 'r';

    uint32   cra = clear[idA].origR;
    uint32   cla = clear[idA].origL;
    uint32   ora = clear[idA].ovlpR;
    uint32   ola = clear[idA].ovlpL;
    uint32   dla = clear[idA].deleted;

    uint32   crb = clear[idB].origR;
    uint32   clb = clear[idB].origL;
    uint32   orb = clear[idB].ovlpR;
    uint32   olb = clear[idB].ovlpL;
    uint32   dlb = clear[idA].deleted;

    leftA  = ovl->dat.obt.a_beg + cla;
    rightA = ovl->dat.obt.a_end + cla;
    lenA   = clear[idA].length;
    leftB  = ovl->dat.obt.b_beg + clb;
    rightB = ovl->dat.obt.b_end + clb;
    lenB   = clear[idB].length;
    error  = AS_OVS_decodeQuality(ovl->dat.obt.erate);

    if (idA != idAlast) {
      process(idAlast, gkp, doUpdate, overlap, olalast, oralast);
      delete overlap;
      overlap = new overlapList;
    }
    idAlast = ovl->a_iid;
    olalast = ola;
    oralast = ora;

#ifdef REPORT_OVERLAPS
    if (reportFile)
      fprintf(reportFile, F_U32"\t"F_U32"\t%c\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t%5.3f\n",
              idA, idB, ori, leftA, rightA, lenA, leftB, rightB, lenB, error);
#endif

    //  Make sure that the overlap at least intersects both of the
    //  final clear ranges.  If it doesn't, discard it.
    //
    bool   intersectA  = (ola < rightA) && (leftA < ora);
    bool   intersectBf = (olb < rightB) && (leftB < orb) && (ori == 'f');
    bool   intersectBr = (orb > rightB) && (leftB > olb) && (ori == 'r');
    bool   failedTrim  = false;
    uint32 diffA       = 0;
    uint32 diffB       = 0;

    //  failedTrim == true can occur frequently on high error
    //  overlaps.  One case BPW has seen was an overlap with 12.6%
    //  error:
    //
    //  a frag overlap: 247 -> 792 = 545 bp covered
    //  b frag overlap:  17 -> 502 = 485 bp covered
    //
    //  The a frag overlap extended 504bp past the OBT clear range,
    //  and so we trimmed both overlaps back by 504 bp.  This thus
    //  resulted in the b frag overlap ending at 502 - 504 = -2, an
    //  error.
    //
    //  failedTrim == true can also occur if we trim both the left and
    //  the right sides of the overlap.  What happens: fragment A says
    //  to trim X bases off of the right end, and fragment B says to
    //  trim Y bases off the left end.  If that's bigger than our
    //  original overlap, then we have nothing left...and we might be
    //  negative.

    if (intersectA && (intersectBf || intersectBr)) {

      //  Rebuild the overlap, so that it doesn't extend past the
      //  overlap clear range.

      uint32 ltrim = 0;
      uint32 rtrim = 0;

      switch (ori) {
        case 'f':
          if (leftA < ola)
            ltrim = ola - leftA;
          if ((leftB < olb) && (ltrim < olb - leftB))
            ltrim = olb - leftB;

          if (ora < rightA)
            rtrim = rightA - ora;
          if ((orb < rightB ) && (rtrim < rightB - orb))
            rtrim = rightB - orb;

          if ((rtrim > rightA) || (rtrim > rightB))
            failedTrim = true;

          leftA  += ltrim;
          rightA -= rtrim;

          leftB  += ltrim;
          rightB -= rtrim;

          if ((leftA > rightA) || (leftB > rightB))
            failedTrim = true;

          diffA = rightA - leftA;
          diffB = rightB - leftB;

          break;
        case 'r':
          if (leftA < ola)
            ltrim = ola - leftA;
          if ((orb < leftB) && (ltrim < leftB - orb))
            ltrim = leftB - orb;

          if (ora < rightA)
            rtrim = rightA - ora;
          if ((rightB < olb ) && (rtrim < olb - rightB))
            rtrim = olb - rightB;

          if ((rtrim > rightA) || (ltrim > leftB))
            failedTrim = true;

          leftA  += ltrim;
          rightA -= rtrim;

          leftB  -= ltrim;
          rightB += rtrim;

          if ((leftA > rightA) || (rightB > leftB))
            failedTrim = true;

          diffA = rightA - leftA;
          diffB = leftB - rightB;

          break;
        default:
          fprintf(stderr, "Unknown ori '%c'\n", ori);
          break;
      }


      //  Just a simple sanity check.  If we're not a failure already,
      //  make sure the end points are plausible.
      //
      if ((failedTrim == false) &&
          ((leftA > 2048) || (rightA > 2048) || (leftB > 2048) || (rightB > 2048))) {
        overflow++;
        fprintf(stderr, "\n");
        fprintf(stderr, "Overflow!  YIKES!\n");
        fprintf(stderr, "A:\torig:\t"F_U32"\t"F_U32"\tovlp:"F_U32"\t"F_U32"\n", cla, cra, ola, ora);
        fprintf(stderr, "B:\torig:\t"F_U32"\t"F_U32"\tovlp:"F_U32"\t"F_U32"\n", clb, crb, olb, orb);
        fprintf(stderr, F_U32"\t"F_U32"\t%c\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t%5.3f\n",
                idA, idB, ori, leftA, rightA, lenA, leftB, rightB, lenB, error);
        failedTrim = true;
      }

      //  Assume it _IS_ crap.
      bool  isNotCrap = false;

      //  Only chance in being not crap is if the trim succeeded.
      if (failedTrim == false) {
        //  We used to use hard and fast cutoffs here.  The original
        //  said use the overlap if its length (rightA-leftA) was >=35
        //  and error <= 0.02, or if length >= 70.  This is completely
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
          //  454 mate read
          //
          //  The smallest overlap we should be seeing from overlapper
          //  is 40bp.  That's nearly half of a 454 mated read, so we
          //  allow much higher error overlaps at that length.
          //
          isNotCrap = (((rightA - leftA >=  40) && (error <= 0.333 * AS_OVL_ERROR_RATE)) ||
                       ((rightA - leftA >=  40) && (error <= 0.750 * AS_OVL_ERROR_RATE)) ||
                       ((rightA - leftA >=  70) && (error <= 1.000 * AS_OVL_ERROR_RATE)));
        } else if (lenA < 400) {
          //  454 unmated read
          isNotCrap = (((rightA - leftA >=  40) && (error <= 0.333 * AS_OVL_ERROR_RATE)) ||
                       ((rightA - leftA >=  70) && (error <= 0.750 * AS_OVL_ERROR_RATE)) ||
                       ((rightA - leftA >= 100) && (error <= 1.000 * AS_OVL_ERROR_RATE)));
        } else {
          //  Sanger read
          isNotCrap = (((rightA - leftA >=  40) && (error <= 0.333 * AS_OVL_ERROR_RATE)) ||
                       ((rightA - leftA >= 100) && (error <= 0.750 * AS_OVL_ERROR_RATE)) ||
                       ((rightA - leftA >= 200) && (error <= 1.000 * AS_OVL_ERROR_RATE)));
        }
      }

      if (isNotCrap) {
        uint32  lhangA=0, rhangA=0, lhangB=0, rhangB=0;

        //  Are we in the middle, or on an end?
        //
        switch (ori) {
          case 'f':
            lhangA = leftA - ola;
            rhangA = ora   - rightA;
            lhangB = leftB - olb;
            rhangB = orb   - rightB;
            break;
          case 'r':
            lhangA = leftA  - ola;
            rhangA = ora    - rightA;
            lhangB = orb    - leftB;
            rhangB = rightB - olb;
            break;
          default:
            fprintf(stderr, "Unknown ori '%c'\n", ori);
            break;
        }

        if (lhangA < 15)  lhangA = 0;
        if (rhangA < 15)  rhangA = 0;
        if (lhangB < 15)  lhangB = 0;
        if (rhangB < 15)  rhangB = 0;

        //  Add to the list.

        if ((dla == 0) && (dlb == 0))
          overlap->add(idA, lhangA, leftA, rightA, rhangA,
                       idB, lhangB, leftB, rightB, rhangB, ori);

     } else {
        notclear++;
#ifdef DEBUG_NOT_IN_CLEAR
        fprintf(stderr, "Overlap not in clear!  Removed!\n");
        fprintf(stderr, "A:\torig:\t"F_U32"\t"F_U32"\tovlp:"F_U32"\t"F_U32"\n", cla, cra, ola, ora);
        fprintf(stderr, "B:\torig:\t"F_U32"\t"F_U32"\tovlp:"F_U32"\t"F_U32"\n", clb, crb, olb, orb);
        fprintf(stderr, F_U32"\t"F_U32"\t%c\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32" %5.3f\n",
                idA, idB, ori, leftA, rightA, lenA, leftB, rightB, lenB, error);
#endif
      }
    }  //  End of intersection test

    ovl = readOverlap(ovsprimary, ovssecondary);
  }

  process(idAlast, gkp, doUpdate, overlap, olalast, oralast);
  delete overlap;

  closeGateKeeperStore(gkp);

  if (summaryFile) {
    fprintf(summaryFile, "fullCoverage:        "F_U32"\n", fullCoverage);
    fprintf(summaryFile, "noInniePair:         "F_U32"\n", noInniePair);
    fprintf(summaryFile, "noChimericOvl:       "F_U32"\n", noChimericOvl);
    fprintf(summaryFile, "chimeraDetected:     "F_U32"\n", chimeraDetected);
    fprintf(summaryFile, "chimeraDeletedSmall: "F_U32"\n", chimeraDeletedSmall);
    fprintf(summaryFile, "spurDetected:        "F_U32"\n", spurDetected);
    fprintf(summaryFile, "spurDeletedSmall:    "F_U32"\n", spurDeletedSmall);
    fprintf(summaryFile, "overlapsNotInClear:  "F_U32" (overlaps ignored)\n", notclear);
  }

  if (overflow)
    fprintf(stderr, "WARNING!  "F_U32" overflows!\n", overflow);

  exit(0);
}
