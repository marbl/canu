
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2011, J. Craig Venter Institute.
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

static const char *rcsid = "$Id: finalTrim-evidenceBased.C,v 1.1 2011-06-03 17:34:19 brianwalenz Exp $";

#include "finalTrim.H"
#include "finalTrim-consolidate.H"


//  If after all that, we have an invalid range, throw it out.
//
//  What has happened (once) is for a fragment:
//
//                OVL                QUAL
//      ------|---------|---------|---------|---------
//              --------
//               -------
//                ------
//
//  (that is, lots of 5' alignment in the overlap, none in 3', no
//  intersection between OVL and QUAL)
//
//  We picked the 5'mode (the right end of OVL) and the left
//  quality (left end of QUAL), which isn't a valid region.
//
//  XXX: We should check if the quality of the OVL region is
//  decent, and then use that if it is also large.


//  Find the mode of the 5'mode.  Read the whole overlap-trim file,
//  counting the 5'mode
//
//mode5 *modes = findModeOfFivePrimeMode(gkp, ovlFile);


bool
evidenceBased(OVSoverlap  *ovl,
              uint32       ovlLen,
              gkFragment  &fr,
              uint32       ibgn,
              uint32       iend,
              uint32      &fbgn,
              uint32      &fend,
              char        *logMsg,
              uint32       errorRate,
              double       errorLimit,
              bool         qvTrimAllowed) {

  obtConsolidate  co(ovl, ovlLen, ibgn);

  uint32  jbgn;  //  If 'ibgn' is the initial quality trim, then 'jbgn' is the
  uint32  jend;  //  second, more restrictive, quality trim.

  assert(fr.gkFragment_getReadIID() == ovl[0].a_iid);

  //  Adjust the mode and max/min(m) counts and values if the
  //  min/max or min/max(m) values are within OBT_MODE_WIGGLE
  //  (currently, 5 bp) of the mode

  if ((co._min5 != co._mode5) && ((co._min5 + OBT_MODE_WIGGLE) >= co._mode5)){
    if (co._minm5 < co._mode5) {
      if (co._min5 < co._minm5) {
        co._mode5c++;
      }
      co._mode5c += co._minm5c;
      co._minm5c = co._mode5c;
    } else if (co._minm5 == co._mode5) {
      co._mode5c++;
      co._minm5c++;
      assert(co._minm5c == co._mode5c);
    } else {
      assert(0);
    }
    co._min5 = co._mode5;
    co._minm5 = co._mode5;
  } else if ((co._minm5 != co._mode5) && ((co._minm5 + OBT_MODE_WIGGLE) >= co._mode5)){
    co._mode5c += co._minm5c;
    co._minm5c = co._mode5c;
    co._minm5 = co._mode5;
  }

  if ((co._max3 != co._mode3) && ((co._max3 - OBT_MODE_WIGGLE) <= co._mode3)){
    if (co._maxm3 > co._mode3) {
      if (co._max3 > co._maxm3) {
        co._mode3c++;
      }
      co._mode3c += co._maxm3c;
      co._maxm3c = co._mode3c;
    } else if (co._maxm3 == co._mode3) {
      co._mode3c++;
      co._maxm3c++;
      assert(co._maxm3c == co._mode3c);
    } else {
      assert(0);
    }
    co._max3 = co._mode3;
    co._maxm3 = co._mode3;
  } else if ((co._maxm3 != co._mode3) && ((co._maxm3 - OBT_MODE_WIGGLE) <= co._mode3)){
    co._mode3c += co._maxm3c;
    co._maxm3c = co._mode3c;
    co._maxm3 = co._mode3;
  }

  double  minQuality = qual.lookupNumber(20);

  //  If we allowed initial trimming based on quality, do another quality trim using
  //  a more stringent parameter.  If not, pick some questionably reasonable values
  //  based on overlaps.
  //
  if (qvTrimAllowed == false) {
    ibgn = 0;
    iend = fr.gkFragment_getSequenceLength();

    jbgn   = co._min5;
    if (co._mode5c > 4) jbgn = co._mode5;
    if (co._minm5c > 8) jbgn = co._minm5;

    jend   = co._max3;
    if (co._mode3c > 4) jend = co._mode3;
    if (co._maxm3c > 8) jend = co._maxm3;
  } else {
    doTrim(&fr, minQuality, jbgn, jend);
  }

  //  This is a flaw in OBT.  If qvTrimAllowed==false, we can occasionally pick jbgn > jend.  Here is one such example.
  //  The 5'mode is picked using the overlaps on the right, and the 3'mode is picked using the overlaps
  //  on the left -- backwards to what we should be doing.
  //
  //  DUMPING PICTURE for ID 1300 in store 0-overlaptrim/test.obtStore (gkp test.gkpStore clear OBTINITIAL)
  //      1300  A:    0  260                                             --------------------------------------------------------------------------------------------------->
  //    141111  A:   10  116 ( 260)  B:   76  179 ( 179)   2.91%        +76 ---------------------------------------->
  //     42961  A:   10  117 ( 260)  B:  143  247 ( 248)   2.88%       +143 -----------------------------------------> +1
  //      5060  A:   10  119 ( 260)  B:  129   23 ( 240)   2.83%       +111 <----------------------------------------- +23
  //     47583  A:   10  119 ( 260)  B:  118   12 ( 252)   2.83%       +134 <----------------------------------------- +12
  //     76173  A:   10  119 ( 260)  B:  122  228 ( 254)   2.83%       +122 -----------------------------------------> +26
  //     93981  A:   10  119 ( 260)  B:  134  240 ( 251)   2.83%       +134 -----------------------------------------> +11
  //     98810  A:   10  119 ( 260)  B:  114  220 ( 241)   2.83%       +114 -----------------------------------------> +21
  //     16771  A:   10  200 ( 260)  B:  186    4 ( 247)   4.40%        +61 <------------------------------------------------------------------------ +4
  //     77751  A:   28  119 ( 260)  B:   60  150 ( 254)   3.33%               +60 ----------------------------------> +104
  //     98695  A:   43  119 ( 260)  B:  110   35 ( 256)   1.33%                    +146 <---------------------------- +35
  //     72005  A:  144  200 ( 260)  B:    6   61 ( 238)   1.82%                                                             +6 --------------------> +177
  //     90602  A:  144  200 ( 260)  B:   55  110 ( 276)   1.82%                                                            +55 --------------------> +166
  //    111027  A:  144  200 ( 260)  B:   67   12 ( 246)   1.82%                                                           +179 <-------------------- +12
  //    159401  A:  144  200 ( 260)  B:    6   61 ( 111)   1.82%                                                             +6 --------------------> +50
  //    179340  A:  144  200 ( 260)  B:   14   69 (  81)   1.82%                                                            +14 --------------------> +12
  //    149520  A:  144  231 ( 260)  B:   84    0 ( 142)   3.57%                                                            +58 <--------------------------------
  //     41765  A:  144  255 ( 260)  B:   58  165 ( 238)   3.74%                                                            +58 ------------------------------------------> +73
  //     83378  A:  144  255 ( 260)  B:  137   30 ( 249)   3.74%                                                           +112 <------------------------------------------ +30
  //    137940  A:  144  255 ( 260)  B:   21  128 ( 171)   3.74%                                                            +21 ------------------------------------------> +43
  //    208013  A:  144  255 ( 260)  B:   63  170 ( 183)   3.74%                                                            +63 ------------------------------------------> +13
  //    223379  A:  144  255 ( 260)  B:  150   43 ( 208)   3.74%                                                            +58 <------------------------------------------ +43
  //
  //  5'mode is 144
  //  3'mode is 119
  //
  //  _min5 = 10,  _minm5 = 10,  _minm5c = 8, _mode5 = 144, _mode5c = 11
  //  _max3 = 255, _maxm3 = 255, _maxm3c = 5, _mode3 = 119, _mode3c = 7
  //
  if (jbgn >= jend) {
    fbgn = jend;
    fend = jbgn;
    strcpy(logMsg, "\t5'mode larger than 3'mode; likely chimeric read");
    return(false);
  }


  //  Use the intersection of the two QV trims.  If there is no intersection, abort.
  //
  if ((iend < jbgn) ||
      (jend < ibgn)) {
    fbgn = jbgn;
    fend = jend;
    strcpy(logMsg, "no overlap between loose QV and strict QV trim");
    return(false);
  }

  //  Guard against errors in trimming or parameters.  Do not allow the
  //  stringent quality trim to exceed the loose quality trim.
  jbgn = MAX(ibgn, jbgn);
  jend = MIN(iend, jend);

#if 0
  //  0) if the left quality trim is less than the mode-of-the-mode,
  //     reset it.
  //
  //  This is exceedingly expensive to compute.  modes[] is storing the mode of all the mode5's for
  //  each library.  To compute it, we must first compute the mode5 for every fragment -- which
  //  needs one entire pass through overlaps.  That is now hugely expensive with large, deep
  //  assemblies.
  //
  if (jbgn < modes[lib].get())
    jbgn = modes[lib].get();
#endif

  assert(jbgn < jend);

  //  1) if the quality-range is < 100, or
  //  2) if the quality-range is < 200 and the intersection with the
  //     overlap-range is < 100
  //  be more conservative
  //
  if ((jbgn + OBT_CQ_LENGTH > jend) ||
      ((jbgn + OBT_CQO_LENGTH > jend) &&
       ((co._min5 + OBT_CQO_OVERLAP > jend) ||
        (jbgn + OBT_CQO_OVERLAP > co._max3)))) {
    //  SHORT QUALITY
    if ((jbgn  + OBT_CQ_SHORT > jend)) {
      //  VERY SHORT QUALITY, or SHORT IN COMMON
      fbgn = jbgn;
      fend = jend;
      strcpy(logMsg, "\tvery short quality, or too short in common");
      return(false);
    }

    //  SHORT QUALITY; use overlap modes instead

    //  This is a valid trim if mode5 < mode3.  For whatever bizarre reasons, sometimes this is not
    //  true.

    if (co._mode5 < co._mode3) {
      fbgn = co._mode5;
      fend = co._mode3;
      return(true);
    } else {
      fbgn = co._mode3;
      fend = co._mode5;
      strcpy(logMsg, "\tshort quality, invalid mode\n");
      return(false);
    }
  }

  //  USE min/max MODE

  if ((co._minm5 < jbgn + OBT_QLT_CLOSE_5) && (co._minm5c > 1)) {
    //  USE MIN>1 5'
    fbgn = co._minm5;
  } else if ((co._mode5 < jbgn + OBT_QLT_CLOSE_5) && (co._mode5c > 0)) {
    //  USE MODE 5'
    fbgn = co._mode5;
  } else if ((co._mode5c > 0) && ((co._min5 <= jbgn) || (co._min5 < jbgn + OBT_QLT_FAR_5))) {
    //  USE MIN 5'
    fbgn = co._min5;
  } else {
    //  USE QUALITY 5;
    fbgn = jbgn;
  }

  if ((co._maxm3 >= jend) && (co._maxm3c > 1)) {
    //  USE MAX>1 3'
    fend = co._maxm3;
  } else if ((co._mode3 == co._max3) && (co._mode3 == co._maxm3) && (co._mode3c > 1) &&
             ((co._mode3 >= jend) || (jend < co._mode3 + OBT_QLT_MODE3))) {
    //  USE MODE 3'
    fend = co._mode3;
  } else if ((co._maxm3c > 1) && (co._maxm3 < jend) && (co._max3 > jend) && (co._max3 < co._maxm3 + OBT_QLT_CLOSE_MAXM3)) {
    //  USE MAX>1 if also CLOSE TO MAX 3'
    fend = co._maxm3;
  } else if ((co._mode3c > 0) && ((co._max3 >= jend) || (jend < co._max3 + OBT_QLT_CLOSE_MAX3))) {
    //  USE MAX 3'
    fend = co._max3;
  } else {
    //  USE QUALITY 3'
    fend = jend;
  }

  assert(fbgn <= fend);
  return(true);
}

