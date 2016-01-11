
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_OVL/AS_OVL_overlap_common.h
 *    src/AS_OVM/overlapInCore-Output.C
 *
 *  Modifications by:
 *
 *    Michael Schatz on 2004-SEP-23
 *      are Copyright 2004 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2005-MAR-22
 *      are Copyright 2005 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUN-16 to 2013-AUG-01
 *      are Copyright 2005-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2005-JUL-15 to 2007-NOV-20
 *      are Copyright 2005,2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2006-MAR-27 to 2006-AUG-21
 *      are Copyright 2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Art Delcher on 2007-FEB-13
 *      are Copyright 2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2007-AUG-27 to 2009-JAN-16
 *      are Copyright 2007,2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2011-MAR-08 to 2011-JUN-17
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUN-16
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "overlapInCore.H"
#include <pthread.h>

//  Output the overlap between strings  S_ID  and  T_ID  which
//  have lengths  S_Len  and  T_Len , respectively.
//  The overlap information is in  (* olap) .
//  S_Dir  indicates the orientation of  S .
//  T is always forward.

void
Output_Overlap(uint32 S_ID, int S_Len, Direction_t S_Dir,
               uint32 T_ID, int T_Len, Olap_Info_t *olap,
               Work_Area_t *WA) {

  ovOverlap  *ovs = WA->overlaps + WA->overlapsLen++;

  //  Overlap is good for UTG only.

  ovs->dat.ovl.forUTG = true;
  ovs->dat.ovl.forOBT = false;
  ovs->dat.ovl.forDUP = false;

  //  Set the span of the alignment.

  ovs->dat.ovl.span  = (olap->s_hi - olap->s_lo);
  ovs->dat.ovl.span += (olap->t_hi - olap->t_lo);
  ovs->dat.ovl.span += (olap->delta_ct);
  assert((ovs->dat.ovl.span % 2) == 0);
  ovs->dat.ovl.span /= 2;

  assert (S_ID < T_ID);

  int32  S_Right_Hang = S_Len - olap->s_hi - 1;
  int32  T_Right_Hang = T_Len - olap->t_hi - 1;

  bool Sleft;

  char    orient = 0;
  int32   ahg    = 0;
  int32   bhg    = 0;

  //fprintf(stderr, "S %d %d = %d len %d T %d %d = %d len %d delta_ct %d delta %d\n",
  //        olap->s_lo, olap->s_hi, olap->s_hi - olap->s_lo, S_Len,
  //        olap->t_lo, olap->t_hi, olap->t_hi - olap->t_lo, T_Len,
  //        olap->delta_ct, olap->delta[0]);

  assert(olap->s_lo < olap->s_hi);
  assert(olap->t_lo < olap->t_hi);

  //  Ensure this is a dovetail or contained overlap.
  assert((olap->s_lo == 0)         || (olap->t_lo == 0));
  assert((olap->s_hi == S_Len - 1) || (olap->t_hi == T_Len - 1));


  if ((olap->s_lo > olap->t_lo) ||
      ((olap->s_lo == olap->t_lo) && (S_Right_Hang > T_Right_Hang))) {
    assert(olap->t_lo == 0);
    assert((olap->s_hi == S_Len - 1) || (olap->t_hi == T_Len - 1));
    Sleft = true;
  } else {
    assert(olap->s_lo == 0);
    assert((olap->s_hi == S_Len - 1) || (olap->t_hi == T_Len - 1));
    Sleft = false;
  }

  if (Sleft) {
    ovs->a_iid = S_ID;
    ovs->b_iid = T_ID;
  } else {
    ovs->a_iid = T_ID;
    ovs->b_iid = S_ID;
  }

  //if (Sleft) {
  //  if (S_Right_Hang >= T_Right_Hang)
  //    overlap_type = AS_CONTAINMENT;
  //  else
  //    overlap_type =  AS_DOVETAIL;
  //} else {
  //  if (T_Right_Hang >= S_Right_Hang)
  //    overlap_type = AS_CONTAINMENT;
  //  else
  //    overlap_type =  AS_DOVETAIL;
  //}

  if (Sleft) {
    orient  = (S_Dir == FORWARD) ? 'N' : 'O';
    ahg     = olap->s_lo;
    bhg     = T_Right_Hang - S_Right_Hang;

  } else {
    orient  = (S_Dir == FORWARD) ? 'N' : 'I';
    ahg     = olap->t_lo;
    bhg     = S_Right_Hang - T_Right_Hang;
  }

  //  CMM: Regularize the reverse orientated containment overlaps to a common orientation.
  //
  //  This catches the case where a reverse orient S (T is always fowrard) is placed
  //  in the A position; we flip the overlap to make S be forward and T be reverse.
  //
  if ((orient == 'O') && (S_Right_Hang >= T_Right_Hang)) {
    orient  = 'I';
    ahg     = -(T_Right_Hang - S_Right_Hang);
    bhg     = -(olap->s_lo);
  }

  ovs->erate(olap->quality);

  switch (orient) {
    case 'N':
      ovs->a_hang(ahg);
      ovs->b_hang(bhg);
      ovs->dat.ovl.flipped  = false;
      break;

    case 'I':
      ovs->a_hang(ahg);
      ovs->b_hang(bhg);
      ovs->dat.ovl.flipped  = true;
      break;

    case 'O':
      ovs->a_hang(-bhg);
      ovs->b_hang(-ahg);
      ovs->dat.ovl.flipped  = true;
      break;

    case 'A':
      ovs->a_hang(-bhg);
      ovs->b_hang(-ahg);
      ovs->dat.ovl.flipped  = false;
      break;
  }

#if  OUTPUT_OVERLAP_DELTAS
  signed char deltas[2 * AS_READ_MAX_NORMAL_LEN];
  signed char *deltaCursor = deltas;

  if (Sleft == false)
    for (i = 0;  i < olap->delta_ct;  i ++)
      olap->delta [i] *= -1;

  for (int i = 0;  i < olap->delta_ct;  i ++) {
    for (int j = abs (olap->delta [i]);  j > 0;  j -= AS_LONGEST_DELTA) {
      if (j > AS_LONGEST_DELTA)
        *deltaCursor++ = AS_LONG_DELTA_CODE;
      else
        *deltaCursor++ = j * Sign (olap->delta [i]);
    }
  }

  *deltaCursor = AS_ENDOF_DELTA_CODE;
#endif


  WA->Total_Overlaps ++;

  if (bhg <= 0)
    WA->Contained_Overlap_Ct ++;
  else
    WA->Dovetail_Overlap_Ct ++;



  //  We also flush the file at the end of a thread

  if (WA->overlapsLen >= WA->overlapsMax) {
    pthread_mutex_lock (& Write_Proto_Mutex);

    for (int32 zz=0; zz<WA->overlapsLen; zz++)
      Out_BOF->writeOverlap(WA->overlaps + zz);
    WA->overlapsLen = 0;

    pthread_mutex_unlock (& Write_Proto_Mutex);
  }
}



void
Output_Partial_Overlap(uint32 s_id,
                       uint32 t_id,
                       Direction_t dir,
                       const Olap_Info_t *olap,
                       int s_len,
                       int t_len,
                       Work_Area_t  *WA) {

  Total_Overlaps++;

  ovOverlap  *ovl = WA->overlaps + WA->overlapsLen++;

  assert(s_id < t_id);

  ovl->a_iid = s_id;
  ovl->b_iid = t_id;

  //  Overlap is good for OBT or DUP.  It will be refined more during the store build.

  ovl->dat.ovl.forUTG = false;
  ovl->dat.ovl.forOBT = true;
  ovl->dat.ovl.forDUP = true;

  //  Set the span of the alignment.

  ovl->dat.ovl.span  = (olap->s_hi - olap->s_lo);
  ovl->dat.ovl.span += (olap->t_hi - olap->t_lo);
  ovl->dat.ovl.span += (olap->delta_ct);
  assert((ovl->dat.ovl.span % 2) == 0);
  ovl->dat.ovl.span /= 2;

  // Convert to canonical form with s forward and use space-based
  // coordinates

  if (dir == FORWARD) {
    ovl->dat.ovl.ahg5 =         (olap->s_lo);
    ovl->dat.ovl.ahg3 = s_len - (olap->s_hi + 1);
    ovl->dat.ovl.bhg5 =         (olap->t_lo);
    ovl->dat.ovl.bhg3 = t_len - (olap->t_hi + 1);

    ovl->dat.ovl.flipped = false;

    //assert(a < b);
    //assert(c < d);

  } else {
    ovl->dat.ovl.ahg5 = s_len - (olap->s_hi + 1);
    ovl->dat.ovl.ahg3 =         (olap->s_lo);
    ovl->dat.ovl.bhg5 = t_len - (olap->t_hi + 1);
    ovl->dat.ovl.bhg3 =         (olap->t_lo);

    ovl->dat.ovl.flipped = true;

    //assert(a < b);
    //assert(c > d);  //  Reverse!
  }

  ovl->erate(olap->quality);

  //  We also flush the file at the end of a thread

  if (WA->overlapsLen >= WA->overlapsMax) {
    pthread_mutex_lock(&Write_Proto_Mutex);

    for (int32 zz=0; zz<WA->overlapsLen; zz++)
      Out_BOF->writeOverlap(WA->overlaps + zz);

    WA->overlapsLen = 0;

    pthread_mutex_unlock(&Write_Proto_Mutex);
  }
}

