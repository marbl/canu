
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

static const char *rcsid = "$Id$";

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

  ovsOverlap  *ovs = WA->overlaps + WA->overlapsLen++;

  //  Overlap is good for UTG only.

  ovs->dat.ovl.forUTG = true;
  ovs->dat.ovl.forOBT = false;
  ovs->dat.ovl.forDUP = false;

  //ovs->clear(olap->quality);

  assert (S_ID < T_ID);

  int32  S_Right_Hang = S_Len - olap->s_hi - 1;
  int32  T_Right_Hang = T_Len - olap->t_hi - 1;

  bool Sleft;

  char    orient = 0;
  int32   ahg    = 0;
  int32   bhg    = 0;

  //fprintf(stderr, "S %d %d %d T %d %d %d\n",
  //        olap->s_lo, olap->s_hi, S_Len,
  //        olap->t_lo, olap->t_hi, T_Len);

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
                       const Olap_Info_t *p,
                       int s_len,
                       int t_len,
                       Work_Area_t  *WA) {

  Total_Overlaps++;

  ovsOverlap  *ovl = WA->overlaps + WA->overlapsLen++;

  assert(s_id < t_id);

  ovl->a_iid = s_id;
  ovl->b_iid = t_id;

  //  Overlap is good for OBT or DUP.  It will be refined more during the store build.

  ovl->dat.ovl.forUTG = false;
  ovl->dat.ovl.forOBT = true;
  ovl->dat.ovl.forDUP = true;


  // Convert to canonical form with s forward and use space-based
  // coordinates

  ovl->dat.ovl.span = 0;

  fprintf(stdout, "S: %6u %6d-%6d  T: %6u %6d-%6d  dir %d\n",
          s_id, p->s_lo, p->s_hi,
          t_id, p->t_lo, p->t_hi, dir);

  if (dir == FORWARD) {
    ovl->dat.ovl.ahg5 =         (p->s_lo);
    ovl->dat.ovl.ahg3 = s_len - (p->s_hi + 1);
    ovl->dat.ovl.bhg5 =         (p->t_lo);
    ovl->dat.ovl.bhg3 = t_len - (p->t_hi + 1);

    ovl->dat.ovl.flipped = false;

    //assert(a < b);
    //assert(c < d);

  } else {
    ovl->dat.ovl.ahg5 = s_len - (p->s_hi + 1);
    ovl->dat.ovl.ahg3 =         (p->s_lo);
    ovl->dat.ovl.bhg5 = t_len - (p->t_hi + 1);
    ovl->dat.ovl.bhg3 =         (p->t_lo);

    ovl->dat.ovl.flipped = true;

    //assert(a < b);
    //assert(c > d);  //  Reverse!
  }

  ovl->erate(p->quality);

  //  We also flush the file at the end of a thread

  if (WA->overlapsLen >= WA->overlapsMax) {
    pthread_mutex_lock(&Write_Proto_Mutex);

    for (int32 zz=0; zz<WA->overlapsLen; zz++)
      Out_BOF->writeOverlap(WA->overlaps + zz);

    WA->overlapsLen = 0;

    pthread_mutex_unlock(&Write_Proto_Mutex);
  }
}

