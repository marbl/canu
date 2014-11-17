
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

//  Output the overlap between strings  S_ID  and  T_ID  which
//  have lengths  S_Len  and  T_Len , respectively.
//  The overlap information is in  (* olap) .
//  S_Dir  indicates the orientation of  S .
//  T is always forward.

void
Output_Overlap(AS_IID S_ID, int S_Len, Direction_t S_Dir,
               AS_IID T_ID, int T_Len, Olap_Info_t *olap,
               Work_Area_t *WA) {

  OVSoverlap  *ovs = WA->overlaps + WA->overlapsLen++;

  ovs->a_iid              = 0;
  ovs->b_iid              = 0;
  ovs->dat.ovl.a_hang     = 0;
  ovs->dat.ovl.b_hang     = 0;
  ovs->dat.ovl.flipped    = FALSE;
  ovs->dat.ovl.orig_erate = AS_OVS_encodeQuality(olap->quality);
  ovs->dat.ovl.corr_erate = ovs->dat.ovl.orig_erate;
  ovs->dat.ovl.type       = AS_OVS_TYPE_OVL;

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

  //  The asserts test that we're storing all the bits in a bit-packed field.

  switch (orient) {
    case 'N':
      ovs->dat.ovl.a_hang   = ahg;
      ovs->dat.ovl.b_hang   = bhg;
      ovs->dat.ovl.flipped  = FALSE;

      assert(ovs->dat.ovl.a_hang == ahg);
      assert(ovs->dat.ovl.b_hang == bhg);
      break;

    case 'I':
      ovs->dat.ovl.a_hang   = ahg;
      ovs->dat.ovl.b_hang   = bhg;
      ovs->dat.ovl.flipped  = TRUE;

      assert(ovs->dat.ovl.a_hang == ahg);
      assert(ovs->dat.ovl.b_hang == bhg);
      break;

    case 'O':
      ovs->dat.ovl.a_hang   = -bhg;
      ovs->dat.ovl.b_hang   = -ahg;
      ovs->dat.ovl.flipped  = TRUE;

      assert(ovs->dat.ovl.a_hang == -bhg);
      assert(ovs->dat.ovl.b_hang == -ahg);
      break;

    case 'A':
      ovs->dat.ovl.a_hang   = -bhg;
      ovs->dat.ovl.b_hang   = -ahg;
      ovs->dat.ovl.flipped  = FALSE;

      assert(ovs->dat.ovl.a_hang == -bhg);
      assert(ovs->dat.ovl.b_hang == -ahg);
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
      AS_OVS_writeOverlap(Out_BOF, WA->overlaps + zz);
    WA->overlapsLen = 0;

    pthread_mutex_unlock (& Write_Proto_Mutex);
  }
}



void
Output_Partial_Overlap(AS_IID s_id,
                       AS_IID t_id,
                       Direction_t dir,
                       const Olap_Info_t *p,
                       int s_len,
                       int t_len,
                       Work_Area_t  *WA) {

  int  a, b, c, d;
  char  dir_ch;

  Total_Overlaps ++;

  // Convert to canonical form with s forward and use space-based
  // coordinates
  if (dir == FORWARD) {
    a = p->s_lo;
    b = p->s_hi + 1;
    c = p->t_lo;
    d = p->t_hi + 1;
    dir_ch = 'f';
  } else {
    a = s_len - p->s_hi - 1;
    b = s_len - p->s_lo;
    c = p->t_hi + 1;
    d = p->t_lo;
    dir_ch = 'r';
  }

  OVSoverlap  *ovl = WA->overlaps + WA->overlapsLen++;
  ovl->a_iid            = s_id;
  ovl->b_iid            = t_id;

  ovl->dat.dat[0]       = 0;
  ovl->dat.dat[1]       = 0;
#if AS_OVS_NWORDS > 2
  ovl->dat.dat[2]       = 0;
#endif

  ovl->dat.obt.fwd      = (dir == FORWARD);
  ovl->dat.obt.a_beg    = a;
  ovl->dat.obt.a_end    = b;
  ovl->dat.obt.b_beg    = c;
  ovl->dat.obt.b_end_hi = d >> 9;
  ovl->dat.obt.b_end_lo = d & 0x1ff;
  ovl->dat.obt.erate    = AS_OVS_encodeQuality(p->quality);
  ovl->dat.obt.type     = AS_OVS_TYPE_OBT;

  //  We also flush the file at the end of a thread

  if (WA->overlapsLen >= WA->overlapsMax) {
    int zz;

    pthread_mutex_lock (& Write_Proto_Mutex);

    for (zz=0; zz<WA->overlapsLen; zz++)
      AS_OVS_writeOverlap(Out_BOF, WA->overlaps + zz);
    WA->overlapsLen = 0;

    pthread_mutex_unlock (& Write_Proto_Mutex);
  }

  return;
}

