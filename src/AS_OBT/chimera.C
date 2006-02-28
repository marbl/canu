#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
}

#include "util++.H"
#include "overlap.H"
#include "maps.H"

//  Reads the output of overlap used to find trim points.  Overlaps
//  were computed using the READSTRUCT_ORIGINAL trimming.  Massages
//  the overlaps to be compliant with the READSTRUCT_OVL trimming.

#define POSITION_BITS     11
#define POSITION_MASK     u64bitMASK(11);
#define LENGTH_MASK       u64bitMASK(11);

#define POSITION_ORIG_R   (0 * POSITION_BITS)
#define POSITION_ORIG_L   (1 * POSITION_BITS)
#define POSITION_OVLP_R   (2 * POSITION_BITS)
#define POSITION_OVLP_L   (3 * POSITION_BITS)
#define POSITION_LENGTH   (4 * POSITION_BITS)

//  WITH_REPORT will report spurs and chimera.  WITH_REPORT_FULL will
//  report unmodified fragments.
//
#define WITH_REPORT
//#define WITH_REPORT_FULL

FILE   *summaryFile = 0L;
FILE   *reportFile  = 0L;

//  Define this to delete chimeric and spur fragments, instead of
//  fixing them.
bool fixChimera = true;

//  General debuggin
//#define DEBUG

//  Yuck!  Needs a better name.  overlapHang_t?
class overlap2_t {
public:
  u64bit   style:4;
  u64bit   Alhang:11;
  u64bit   Abeg:11;
  u64bit   Aend:11;
  u64bit   Arhang:11;
  u64bit   Biid:25;
  u64bit   Blhang:11;
  u64bit   Bbeg:11;
  u64bit   Bend:11;
  u64bit   Brhang:11;
};

class overlapList {
public:
  overlapList() {
    _iid    = ~u64bitZERO;
    _ovlMax = 16;
    _ovlLen = 0;
    _ovl    = new overlap2_t [_ovlMax];
  };
  ~overlapList() {
    delete [] _ovl;
  };

  void          add(u64bit Aiid, u32bit Alhang, u32bit Abeg, u32bit Aend, u32bit Arhang,
                    u64bit Biid, u32bit Blhang, u32bit Bbeg, u32bit Bend, u32bit Brhang, char ori) {
    if (_ovlLen >= _ovlMax) {
      _ovlMax *= 2;
      overlap2_t *O = new overlap2_t [_ovlMax];
      memcpy(O, _ovl, sizeof(overlap2_t) * _ovlLen);
      delete [] _ovl;
      _ovl = O;
    }

    if (_iid == ~u64bitZERO)
      _iid = Aiid;
    if (_iid != Aiid)
      fprintf(stderr, "ERROR: adding "u64bitFMT" to overlapList with iid="u64bitFMT"\n", Aiid, _iid);

    u32bit   style = 0;
    if (Alhang > 0)  style |= (u32bitONE << 3);
    if (Arhang > 0)  style |= (u32bitONE << 2);
    if (Blhang > 0)  style |= (u32bitONE << 1);
    if (Brhang > 0)  style |= (u32bitONE);

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
        fprintf(stderr, "UNCLASSIFIED OVERLAP TYPE "u64bitFMT"\n", style);
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

  void          print(FILE *out, u32bit i) {
    fprintf(out,
            u64bitFMTW(7)" "u64bitFMTW(7)" "u64bitFMTW(02)" "
            u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" -- "
            u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" "u64bitFMTW(4)" %5.3f",
            _iid, _ovl[i].Biid, _ovl[i].style,
            _ovl[i].Alhang, _ovl[i].Abeg, _ovl[i].Aend, _ovl[i].Arhang,
            _ovl[i].Blhang, _ovl[i].Bbeg, _ovl[i].Bend, _ovl[i].Brhang, 0.0);
    switch (_ovl[i].style) {
      case 5:
      case 7:
      case 10:
      case 11:
      case 13:
      case 14:
        fprintf(out, " *****");
    }
    fprintf(out, "\n");
  };

  u32bit        length(void) {
    return(_ovlLen);
  };
  overlap2_t    *get(u32bit i) {
    if (i < _ovlLen)
      return(_ovl + i);
    else
      return(0L);
  };

private:
  u64bit       _iid;
  overlap2_t   *_ovl;
  u32bit       _ovlLen;
  u32bit       _ovlMax;
};




u32bit   chimeraDetected     = 0;
u32bit   chimeraDeletedSmall = 0;
u32bit   spurDetected        = 0;
u32bit   spurDeletedSmall    = 0;
u32bit   noInniePair         = 0;
u32bit   noChimericOvl       = 0;
u32bit   fullCoverage        = 0;





//  Reads both clear ranges, packs into a 64-bit integer.
//
u64bit*
readClearRanges(char *frgStore) {

  FragStoreHandle fs = openFragStore(frgStore, "r");
  if (fs == NULLSTOREHANDLE)
    fprintf(stderr, "Failed to open fragStore %s!\n", frgStore), exit(1);

  u32bit            fe = getFirstElemFragStore(fs);
  u32bit            le  = getLastElemFragStore(fs) + 1;
  ReadStructp       rd = new_ReadStruct();
  u64bit           *clear = new u64bit [le];

  fprintf(stderr, "Reading clear ranges for frags "u32bitFMT" to "u32bitFMT".\n", fe, le-1);

  for (u32bit iid=0; iid<le; iid++)
    clear[iid] = ~u64bitZERO;

  for (u32bit iid=fe; iid<le; iid++) {
    unsigned int  origl=0, origr=0, ovlpl=0, ovlpr=0;

    getFragStore(fs, iid, FRAG_S_ALL, rd);
    getClearRegion_ReadStruct(rd, &origl, &origr, READSTRUCT_ORIGINAL);
    getClearRegion_ReadStruct(rd, &ovlpl, &ovlpr, READSTRUCT_OVL);

    clear[iid]   = getSequence_ReadStruct(rd, 0L, 0L, 0);
    clear[iid] <<= POSITION_BITS;
    clear[iid]  |= ovlpl;
    clear[iid] <<= POSITION_BITS;
    clear[iid]  |= ovlpr;
    clear[iid] <<= POSITION_BITS;
    clear[iid]  |= origl;
    clear[iid] <<= POSITION_BITS;
    clear[iid]  |= origr;
  }

  delete_ReadStruct(rd);
  closeFragStore(fs);

  return(clear);
}








//    Sort each overlapList?
//    Build an interval list using all overlaps.
//    Squash the list, see if it intersects.
//    Look for chimeric and spur patterns.
//
void
process(u32bit           iid,
        FragStoreHandle  fs,
        bool             doUpdate,
        vectorMap       &m,
        overlapList     *overlap,
        u32bit           ola,
        u32bit           ora) {

  if (overlap->length() > 0) {
    intervalList   IL;
    bool           leftIntervalHang[1025], rightIntervalHang[1025];
    u32bit         hasPotentialChimera = 0;
    u32bit         hasInniePair    = 0;
    bool           isSpur = false;

    for (u32bit i=0; i<overlap->length(); i++) {
      overlap2_t  *ovl = overlap->get(i);

      switch (ovl->style) {
        case 5:
        case 7:
          //  A is anchored on the left.
          hasPotentialChimera++;
          IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg - 10);
          break;

        case 13:
	  if ((ovl->Aend - ovl->Abeg) > 75) {
	    hasPotentialChimera++;
	    IL.add(ovl->Abeg + 10, ovl->Aend - ovl->Abeg - 20);
	  }
          break;
          
        case 10:
        case 11:
          //  A is anchored on the right.
          hasPotentialChimera++;
          IL.add(ovl->Abeg + 10, ovl->Aend - ovl->Abeg - 10);
          break;

        case 14:
	  if ((ovl->Aend - ovl->Abeg) > 75) {
	    hasPotentialChimera++;
	    IL.add(ovl->Abeg + 10, ovl->Aend - ovl->Abeg - 20);
	  }
          break;

        case 6:
          //  Trust normal overlaps except the last 10 bp
          IL.add(ovl->Abeg + 10, ovl->Aend - ovl->Abeg - 10);
          break;
        case 9:
          //  Trust normal overlaps except the last 10 bp
          IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg - 10);
          break;

        case 1:
          //  Trust containment, B contains A right
          IL.add(ovl->Abeg + 10, ovl->Aend - ovl->Abeg - 10);
          break;          

        case 2:
          //  Trust containment, B contains A left
          IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg - 10);
          break;          

        case 3:
          //  Trust containment, B contains A both
          IL.add(ovl->Abeg, ovl->Aend - ovl->Abeg);
          break;          

        case 4:
        case 8:
        case 12:
          //  Trust containment, A contains B
          IL.add(ovl->Abeg + 10, ovl->Aend - ovl->Abeg - 20);
          break;          

        case 0:
          //  Don't trust things that look like duplicate reads.
          break;

        case 15:
          //  Repeats.
	  if ((ovl->Aend - ovl->Abeg) > 75) {
	    IL.add(ovl->Abeg + 10, ovl->Aend - ovl->Abeg - 20);
	  }
          break;

        default:
          fprintf(stderr, "UNCLASSIFIED OVERLAP TYPE "u64bitFMT"\n", ovl->style);
          break;
      }
    }

    IL.merge();


    //  Run through the overlaps again, counting the number of innie
    //  pairs across each gap in the intervals.
    //  Also mark hangs at the ends of intervals.
    //
    for (u32bit interval=0; interval<=IL.numberOfIntervals(); interval++) {
      u32bit  begGap = (interval == 0) ? ola : IL.hi(interval-1);
      u32bit  endGap = (interval == IL.numberOfIntervals()) ? ora : IL.lo(interval);

      u32bit  l = 0;
      u32bit  r = 0;

      //  initialize interval hang marks
      //
      assert(interval < 1025);
      leftIntervalHang[interval] = false;
      rightIntervalHang[interval] = false;

      //  Count the number of overlaps with hangs that are
      //  on the correct side to be chimeric.
      //
      for (u32bit i=0; i<overlap->length(); i++) {
        overlap2_t  *ovl = overlap->get(i);

        switch (ovl->style) {
	case 5:
	case 7:
	  //  These should be to the left of the endGap to count.
	  if (((ovl->Aend - 10) < endGap) && (ovl->Aend >= begGap)) {
	    l++;
            assert(interval > 0);
            rightIntervalHang[interval-1] = true;
	  }
	  break;
          
	case 13:
	  if ((ovl->Aend - ovl->Abeg) > 75) {
	    //  These should be to the left of the endGap to count.
	    if (((ovl->Aend - 10) < endGap) && (ovl->Aend >= begGap)) {
	      l++;
	      assert(interval > 0);
	      rightIntervalHang[interval-1] = true;
	    }
	  }
	  break;
          
	case 10:
	case 11:
	  //  These should be to the right of the begGap to count.
	  if (((ovl->Abeg + 10) > begGap) && (ovl->Abeg <= endGap)) {
	    r++;
            assert(interval < IL.numberOfIntervals());
	    leftIntervalHang[interval] = true;
	  }
	  break;
	  
	case 14:
	  if ((ovl->Aend - ovl->Abeg) > 75) {
	    //  These should be to the right of the begGap to count.
	    if (((ovl->Abeg + 10) > begGap) && (ovl->Abeg <= endGap)) {
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
	    if (((ovl->Aend - 10) < endGap) && (ovl->Aend >= begGap)) {
	      l++;
              assert(interval  > 0);
	      rightIntervalHang[interval-1] = true;
	    }
	    //  These should be to the right of the begGap to count.
	    if (((ovl->Abeg + 10) > begGap) && (ovl->Abeg <= endGap)) {
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

    {
      u32bit  minOvl = 65536;
      u32bit  maxOvl = 0;
      bool    isLeftSpur = false;
      bool    isRightSpur = false;

      //  Check if the last overlaps on the left or
      //  on the right are spurs.
      //
      for (u32bit i=0; i<overlap->length(); i++) {
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
	    fprintf(stderr, "UNCLASSIFIED OVERLAP TYPE "u64bitFMT"\n", ovl->style);
	    break;
        }
      }

      isSpur = isLeftSpur || isRightSpur;
    }


    //  And another pass, this time to find the largest interval or
    //  non-interval.  If the fragment is chimeric or a spur, we reset
    //  the clear to this interval, rather than killing the fragment
    //  completely.
    //
    u32bit  intervalBeg = 0;
    u32bit  intervalEnd = 0;
    u32bit  intervalMax = 0;
    u32bit  currentBeg = ola;
    u32bit  currentEnd = 0;
    for (u32bit interval=0; interval<=IL.numberOfIntervals(); interval++) {

      if (interval == IL.numberOfIntervals()) {
	currentEnd = ora;
	if ((currentEnd - currentBeg) > intervalMax) {
	  intervalBeg = currentBeg;
	  intervalEnd = currentEnd;
	  intervalMax = intervalEnd - intervalBeg;
	}
	break;
      }
      if (leftIntervalHang[interval]) {
        //  UL.lo(interval) - currentBeg - 10 > intervalMax + 10
        if (IL.lo(interval) > intervalMax + 10 + currentBeg) {
	  intervalBeg = currentBeg;
	  intervalEnd = IL.lo(interval) - 10;
	  intervalMax = intervalEnd - intervalBeg;
	}
	currentBeg = IL.lo(interval);
      }
      if (rightIntervalHang[interval]) {
        //  IL.hi(interval) - currentBeg > intervalMax
        if (IL.hi(interval) > intervalMax + currentBeg) {
	  intervalBeg = currentBeg;
	  intervalEnd = IL.hi(interval);
	  intervalMax = intervalEnd - intervalBeg;
	}
	currentBeg = IL.hi(interval) + 10;
      }
    }

    ReadStructp       rd    = new_ReadStruct();
    u64bit            uid   = 0;

    if (doUpdate) {
      getFragStore(fs, iid, FRAG_S_ALL, rd);
      getAccID_ReadStruct(rd, &uid);

      //  If we don't have the uid in our map, we are allowed to
      //  modify it.  We just reset the doUpdate flag to disallow
      //  modification.  The flag is passed by value, so the change
      //  is local to here.
      //
      if (m.exists(uid))
        doUpdate = false;
    }

#if defined(WITH_REPORT) || defined(WITH_REPORT_FULL)
    //  If we are with reporting, we need the uid regardless of the
    //  doUpdate status.
    //
    if (uid == 0) {
      getFragStore(fs, iid, FRAG_S_ALL, rd);
      getAccID_ReadStruct(rd, &uid);
    }
#endif


    if (isSpur) {
      spurDetected++;

      if (fixChimera) {
        if (doUpdate) {
          setClearRegion_ReadStruct(rd, intervalBeg, intervalEnd, READSTRUCT_OVL);
          if (setFragStore(fs, iid, rd)) {
            fprintf(stderr, "setFragStore() failed.\n");
            exit(1);
          }
        }

        if (intervalMax < 40) {
          spurDeletedSmall++;
          if (doUpdate)
            deleteFragStore(fs, iid);
        }
      } else {
        if (doUpdate)
          deleteFragStore(fs, iid);
      }

#ifdef WITH_REPORT
      if (reportFile) {
        fprintf(reportFile, "----------------------------------------SPUR!("u32bitFMT","u32bitFMT") update=%d\n",
                intervalBeg, intervalEnd, doUpdate);
        fprintf(reportFile, u64bitFMT","u32bitFMT" has "u32bitFMT" intervals and "u32bitFMT" potential chimeric overlaps (%5.2f%%).\n",
                uid, iid, IL.numberOfIntervals(), hasPotentialChimera,
                (double)hasPotentialChimera / (double)overlap->length() * 100);
        for (u32bit i=0; i<overlap->length(); i++)
          overlap->print(reportFile, i);
      }
#endif
    } else if ((IL.numberOfIntervals() > 1) &&
               (hasPotentialChimera > 0) &&
               (hasInniePair > 0)) {
      chimeraDetected++;

      if (fixChimera) {
        if (doUpdate) {
          setClearRegion_ReadStruct(rd, intervalBeg, intervalEnd, READSTRUCT_OVL);
          if (setFragStore(fs, iid, rd)) {
            fprintf(stderr, "setFragStore() failed.\n");
            exit(1);
          }
        }

        if (intervalMax < 40) {
          chimeraDeletedSmall++;
          if (doUpdate)
            deleteFragStore(fs, iid);
        }
      } else {
        if (doUpdate)
          deleteFragStore(fs, iid);
      }

#ifdef WITH_REPORT
      if (reportFile) {
        fprintf(reportFile, "----------------------------------------CHIMERA!("u32bitFMT","u32bitFMT") update=%d\n",
                intervalBeg, intervalEnd, doUpdate);
        fprintf(reportFile, u64bitFMT","u32bitFMT" has "u32bitFMT" intervals and "u32bitFMT" potential chimeric overlaps (%5.2f%%).\n",
                uid, iid, IL.numberOfIntervals(), hasPotentialChimera,
                (double)hasPotentialChimera / (double)overlap->length() * 100);
        for (u32bit i=0; i<overlap->length(); i++)
          overlap->print(reportFile, i);
      }
#endif
    } else if (IL.numberOfIntervals() == 1) {
      fullCoverage++;

#ifdef WITH_REPORT_FULL
      if (reportFile) {
        fprintf(reportFile, "----------------------------------------FULL COVERAGE\n");
        fprintf(reportFile, u64bitFMT","u32bitFMT" has "u32bitFMT" intervals and "u32bitFMT" potential chimeric overlaps (%5.2f%%).\n",
                uid, iid, IL.numberOfIntervals(), hasPotentialChimera,
                (double)hasPotentialChimera / (double)overlap->length() * 100);
      }
#endif
    } else if (hasPotentialChimera == 0) {
      noChimericOvl++;

#ifdef WITH_REPORT_FULL
      if (reportFile) {
        fprintf(reportFile, "----------------------------------------NO CHIMERIC OVERLAPS\n");
        fprintf(reportFile, u64bitFMT","u32bitFMT" has "u32bitFMT" intervals and "u32bitFMT" potential chimeric overlaps (%5.2f%%).\n",
                uid, iid, IL.numberOfIntervals(), hasPotentialChimera,
                (double)hasPotentialChimera / (double)overlap->length() * 100);
        for (u32bit i=0; i<overlap->length(); i++)
          overlap->print(reportFile, i);
      }
#endif
    } else if (hasInniePair == 0) {
      noInniePair++;

#ifdef WITH_REPORT_FULL
      if (reportFile) {
        fprintf(reportFile, "----------------------------------------NO INNIE PAIR (innie="u32bitFMT")\n",
                hasInniePair);
        fprintf(reportFile, u64bitFMT","u32bitFMT" has "u32bitFMT" intervals and "u32bitFMT" potential chimeric overlaps (%5.2f%%).\n",
                uid, iid, IL.numberOfIntervals(), hasPotentialChimera,
                (double)hasPotentialChimera / (double)overlap->length() * 100);
        for (u32bit i=0; i<overlap->length(); i++)
          overlap->print(reportFile, i);
      }
#endif
    } else {
#ifdef WITH_REPORT_FULL
      if (reportFile) {
        fprintf(reportFile, "----------------------------------------NOT CHIMERA, don't know why\n");
        fprintf(reportFile, u64bitFMT","u32bitFMT" has "u32bitFMT" intervals and "u32bitFMT" potential chimeric overlaps (%5.2f%%).\n",
                uid, iid, IL.numberOfIntervals(), hasPotentialChimera,
                (double)hasPotentialChimera / (double)overlap->length() * 100);
        for (u32bit i=0; i<overlap->length(); i++)
          overlap->print(reportFile, i);
      }
#endif
    }

    //  Clean up our read struct.
    delete_ReadStruct(rd);
  }
}






bool
readOverlap(FILE *file, overlap_t &ovl) {
  static char line[1024];

  if (feof(file))
    return(false);

#ifdef ASCII_OVERLAPS
  fgets(line, 1024, stdin);
  ovl.decode(line, false);
#else
  ovl.load(file);
#endif

  if (feof(file))
    return(false);

  return(true);
}





int
main(int argc, char **argv) {
  char   *frgStore          = 0L;
  bool    doUpdate          = true;  //  set to false for testing
  char   *summaryName       = 0L;
  char   *reportName        = 0L;
  char   *immutableFileName = 0L;

  u32bit  overflow = 0;
  u32bit  notclear = 0;

  int arg=1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-frg", 2) == 0) {
      frgStore = argv[++arg];
    } else if (strncmp(argv[arg], "-immutable", 2) == 0) {
      immutableFileName = argv[++arg];
    } else if (strncmp(argv[arg], "-delete", 2) == 0) {
      fixChimera = false;
    } else if (strncmp(argv[arg], "-summary", 2) == 0) {
      summaryName = argv[++arg];
    } else if (strncmp(argv[arg], "-report", 2) == 0) {
      reportName = argv[++arg];
    }
    arg++;
  }

  if (frgStore == 0L) {
    fprintf(stderr, "usage: %s [-1] -f <fragStore> [opts] < overlap-trim-results\n", argv[0]);
    fprintf(stderr, "  -immutable I    don't modify the uids listed in I\n");
    fprintf(stderr, "  -update         do update the frag store\n");
    fprintf(stderr, "  -delete         instead of fixing chimera, delete them\n");
    fprintf(stderr, "  -summary S      write a summary of the fixes to S\n");
    fprintf(stderr, "  -report R       write a detailed report of the fixes to R\n");
    exit(1);
  }

  u32bit           idAlast, olalast, oralast;

  u32bit           idA, idB;
  char             ori;
  u32bit           leftA, rightA, lenA;
  u32bit           leftB, rightB, lenB;
  double           error;
  u64bit          *clear = readClearRanges(frgStore);
  u64bit           maxIID      = 65536;
  overlapList     *overlap = new overlapList;

  //  Our vectorMap is used exclusively for testing if the fragment is
  //  mutable.  If the uid has a key int the map, then it is NOT
  //  mutable.  Be sure to update tests for this if more stuff is
  //  inserted into the map (search for exists()).
  //
  vectorMap  m;
  m.readImmutableMap(immutableFileName);

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

  FragStoreHandle fs = openFragStore(frgStore, doUpdate ? "r+" : "r");
  if (fs == NULLSTOREHANDLE)
    fprintf(stderr, "Failed to open fragStore %s!\n", frgStore), exit(1);

#ifdef SPEEDCOUNTER_H
  speedCounter  *C = new speedCounter("%7.2f Moverlaps -- %5.2f Moverlaps/second\r",
                                      1000000.0, 0x7fff, true);
  C->enableLiner();
#endif

  overlap_t ovl;
  while (readOverlap(stdin, ovl)) {

    idA    = ovl.Aiid;
    idB    = ovl.Biid;
    ori    = (ovl.ori) ? 'f' : 'r';

    u32bit   cra = (clear[idA] >> POSITION_ORIG_R) & POSITION_MASK;
    u32bit   cla = (clear[idA] >> POSITION_ORIG_L) & POSITION_MASK;
    u32bit   ora = (clear[idA] >> POSITION_OVLP_R) & POSITION_MASK;
    u32bit   ola = (clear[idA] >> POSITION_OVLP_L) & POSITION_MASK;

    u32bit   crb = (clear[idB] >> POSITION_ORIG_R) & POSITION_MASK;
    u32bit   clb = (clear[idB] >> POSITION_ORIG_L) & POSITION_MASK;
    u32bit   orb = (clear[idB] >> POSITION_OVLP_R) & POSITION_MASK;
    u32bit   olb = (clear[idB] >> POSITION_OVLP_L) & POSITION_MASK;

    leftA  = ovl.Abeg + cla;
    rightA = ovl.Aend + cla;
    lenA   = (clear[idA] >> POSITION_LENGTH) & LENGTH_MASK;
    leftB  = ovl.Bbeg + clb;
    rightB = ovl.Bend + clb;
    lenB   = (clear[idB] >> POSITION_LENGTH) & LENGTH_MASK;
    error  = ovl.erate;

    if (idA != idAlast) {
      process(idAlast, fs, doUpdate, m, overlap, olalast, oralast);
      delete overlap;
      overlap = new overlapList;
    }
    idAlast = idA;
    olalast = ola;
    oralast = ora;

#ifdef DEBUG
    if (reportFile) {
      fprintf(reportFile, "----------------------------------------\n");
      fprintf(reportFile, u32bitFMTW(7)" "u32bitFMTW(7)"  %c "u32bitFMTW(4)" "u32bitFMTW(4)" "u32bitFMTW(4)"  "u32bitFMTW(4)" "u32bitFMTW(4)" "u32bitFMTW(4)" %5.3f\n",
              idA, idB, ori, leftA, rightA, lenA, leftB, rightB, lenB, error);
    }
#endif

    //  Make sure that the overlap at least intersects both of the
    //  final clear ranges.  If it doesn't, discard it.
    //
    bool  intersectA  = (ola < rightA) && (leftA < ora);
    bool  intersectBf = (olb < rightB) && (leftB < orb) && (ori == 'f');
    bool  intersectBr = (orb > rightB) && (leftB > olb) && (ori == 'r');
    bool  failed      = false;

    if (intersectA && (intersectBf || intersectBr)) {

      //  Rebuild the overlap, so that it doesn't extend past the
      //  overlap clear range.

      u32bit ltrim = 0;
      u32bit rtrim = 0;

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

          leftA  += ltrim;
          rightA -= rtrim;

          leftB  += ltrim;
          rightB -= rtrim;

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

          leftA  += ltrim;
          rightA -= rtrim;

          leftB  -= ltrim;
          rightB += rtrim;

          break;
        default:
          fprintf(stderr, "Unknown ori '%c'\n", ori);
          break;
      }


#ifdef DEBUG
      if (reportFile) {
        fprintf(reportFile, "A: orig: "u32bitFMT" "u32bitFMT"  ovlp:"u32bitFMT" "u32bitFMT"\n", cla, cra, ola, ora);
        fprintf(reportFile, "B: orig: "u32bitFMT" "u32bitFMT"  ovlp:"u32bitFMT" "u32bitFMT"\n", clb, crb, olb, orb);
        fprintf(reportFile, u32bitFMTW(7)" "u32bitFMTW(7)"  %c "u32bitFMTW(4)" "u32bitFMTW(4)" "u32bitFMTW(4)"  "u32bitFMTW(4)" "u32bitFMTW(4)" "u32bitFMTW(4)" %5.3f\n",
                idA, idB, ori, leftA, rightA, lenA, leftB, rightB, lenB, error);
      }
#endif


      //  We've seen a few cases where the resulting overlap is VERY
      //  small, and near the start of one fragment.  For reasons that
      //  may or may not be OK, the clear then becomes negative.  We don't use
      //  overlaps less than 35 bases long anyway, so just stop thinking now
      //  if we do get a short overlap.
      //
      //  Another case where it would have been useful to not be unsigned.
      //
      u32bit diffA = rightA - leftA;
      if (rightA < leftA)  diffA = leftA - rightA;

      u32bit diffB = rightB - leftB;
      if (rightB < leftB)  diffB = leftB - rightB;


      //  Another simple test -- if the A or B overlap range has
      //  'flipped' then we don't have any overlap left over after
      //  applying the new clear range.  What happens: fragment A says
      //  take X bases off of the right end, and fragment B says take
      //  Y bases off the left end.  If that's bigger than our
      //  original overlap, then we have nothing left.  This is
      //  indicated by the A overlap 'flipping' so that the left is
      //  bigger than the right.
      //
      //  BUT, if it flips AND is near the start of the overlap, we
      //  can see negative coordinates.  As unsigned, these come out
      //  huge.  So, we also check that both ranges are above 34 bp
      //  long.
      //
      //  And, so, we have rather complicated set of conditions to
      //  skip this overlap if it's crap.
      //
      //  The comments:
      //    1 - Both ranges are not small.  If we have a negative, one
      //        range will be REAL big, the other will be small.
      //    2 - Not flipped.  If we had signed values, this would be
      //        all we need.  But we're unsigned, so we also need #1.
      //    3 - Small overlap, but great quality.
      //    4 - Normal, large overlap.
      // 
      if ( (diffA > 34) && (diffB > 34) &&                //  1
           (leftA < rightA) &&                            //  2
           (((rightA - leftA > 34) && (error <= 2.0)) ||  //  3
            (rightA - leftA > 69))) {                     //  4

        //  Just a simple sanity check.
        //
        if ((leftA > 2048) || (rightA > 2048) || (leftB > 2048) || (rightB > 2048)) {
          overflow++;
          fprintf(stderr, "\n");
          fprintf(stderr, "Overflow!  YIKES!\n");
          fprintf(stderr, "A: orig: "u32bitFMT" "u32bitFMT"  ovlp:"u32bitFMT" "u32bitFMT"\n", cla, cra, ola, ora);
          fprintf(stderr, "B: orig: "u32bitFMT" "u32bitFMT"  ovlp:"u32bitFMT" "u32bitFMT"\n", clb, crb, olb, orb);
          fprintf(stderr, u32bitFMTW(7)" "u32bitFMTW(7)"  %c "u32bitFMTW(4)" "u32bitFMTW(4)" "u32bitFMTW(4)"  "u32bitFMTW(4)" "u32bitFMTW(4)" "u32bitFMTW(4)" %5.3f\n",
                  idA, idB, ori, leftA, rightA, lenA, leftB, rightB, lenB, error);
        }

        //  Do something with the overlap
        //
        u32bit  lhangA=0, rhangA=0, lhangB=0, rhangB=0;

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

        overlap->add(idA, lhangA, leftA, rightA, rhangA,
                     idB, lhangB, leftB, rightB, rhangB, ori);

      } else {
        notclear++;
#if 0
        fprintf(stderr, "Overlap not in clear!  Removed!\n");
        fprintf(stderr, "A: orig: "u32bitFMT" "u32bitFMT"  ovlp:"u32bitFMT" "u32bitFMT"\n", cla, cra, ola, ora);
        fprintf(stderr, "B: orig: "u32bitFMT" "u32bitFMT"  ovlp:"u32bitFMT" "u32bitFMT"\n", clb, crb, olb, orb);
        fprintf(stderr, u32bitFMTW(7)" "u32bitFMTW(7)" N%c "u32bitFMTW(4)" "u32bitFMTW(4)" "u32bitFMTW(4)"  "u32bitFMTW(4)" "u32bitFMTW(4)" "u32bitFMTW(4)" %5.3f\n",
                idA, idB, ori, leftA, rightA, lenA, leftB, rightB, lenB, error);
#endif
      }
    }  //  End of intersection test

#ifdef SPEEDCOUNTER_H
    C->tick();
#endif
  }

  process(idAlast, fs, doUpdate, m, overlap, olalast, oralast);
  delete overlap;

#ifdef SPEEDCOUNTER_H
  delete C;
#endif

  if (notclear)
    fprintf(stderr, u32bitFMT" overlaps not in the clear and were discarded.\n", notclear);
  if (overflow)
    fprintf(stderr, u32bitFMT" overflows (that's bad).\n", overflow), exit(1);

  closeFragStore(fs);


  if (summaryFile) {
    fprintf(summaryFile, "fullCoverage:        "u32bitFMT"\n", fullCoverage);
    fprintf(summaryFile, "noInniePair:         "u32bitFMT"\n", noInniePair);
    fprintf(summaryFile, "noChimericOvl:       "u32bitFMT"\n", noChimericOvl);
    fprintf(summaryFile, "chimeraDetected:     "u32bitFMT"\n", chimeraDetected);
    fprintf(summaryFile, "chimeraDeletedSmall: "u32bitFMT"\n", chimeraDeletedSmall);
    fprintf(summaryFile, "spurDetected:        "u32bitFMT"\n", spurDetected);
    fprintf(summaryFile, "spurDeletedSmall:    "u32bitFMT"\n", spurDeletedSmall);
  }

  exit(0);
}
