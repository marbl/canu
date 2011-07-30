
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

static const char *rcsid = "$Id: overlapInCore-Output.C,v 1.1 2011-07-30 01:16:05 brianwalenz Exp $";

#include "overlapInCore.H"

//  Output the overlap between strings  S_ID  and  T_ID  which
//  have lengths  S_Len  and  T_Len , respectively.
//  The overlap information is in  (* olap) .
//  S_Dir  indicates the orientation of  S .

void
Output_Overlap(AS_IID S_ID, int S_Len, Direction_t S_Dir,
               AS_IID T_ID, int T_Len, Olap_Info_t * olap,
               Work_Area_t *WA) {


#ifdef COMPARE
  fprintf(stderr, "OUTPUT %d/%d %d/%d\n", S_ID, S_Len, T_ID, T_Len);
#endif

  int  S_Right_Hang, T_Right_Hang;
  int  this_diag;
  OverlapMesg ovMesg;
  GenericMesg outputMesg;
  signed char deltas[2 * AS_READ_MAX_NORMAL_LEN];
  signed char *deltaCursor = deltas;
  outputMesg.m = &ovMesg;
  outputMesg.t = MESG_OVL;


#if  OUTPUT_OVERLAP_DELTAS
  ovMesg.alignment_delta = deltas;
  *deltas = '\0';
#endif
  assert (S_ID < T_ID);

  S_Right_Hang = S_Len - olap -> s_hi - 1;
  T_Right_Hang = T_Len - olap -> t_hi - 1;
  this_diag = olap -> t_lo - olap -> s_lo;

  if  (olap -> s_lo > olap -> t_lo
       || (olap -> s_lo == olap -> t_lo
           && S_Right_Hang > T_Right_Hang))
    {       // S is on the left
      ovMesg.aifrag = (AS_IID) S_ID;
      ovMesg.bifrag = (AS_IID) T_ID;

      if  (S_Dir == FORWARD)
        ovMesg.orientation.setIsNormal();
      else
        ovMesg.orientation.setIsOuttie();

      if  (S_Right_Hang >= T_Right_Hang)
        ovMesg.overlap_type = AS_CONTAINMENT;
      else
        ovMesg.overlap_type =  AS_DOVETAIL;

      ovMesg.ahg = olap -> s_lo;
      ovMesg.bhg = T_Right_Hang - S_Right_Hang;
      ovMesg.quality = olap -> quality;
      ovMesg.min_offset = olap -> s_lo - (olap -> max_diag - this_diag);
      ovMesg.max_offset = olap -> s_lo + (this_diag - olap -> min_diag);
    }
  else
    {
#if  OUTPUT_OVERLAP_DELTAS
      for  (i = 0;  i < olap -> delta_ct;  i ++)
        olap -> delta [i] *= -1;
#endif

      ovMesg.bifrag = (AS_IID) S_ID;
      ovMesg.aifrag = (AS_IID) T_ID;
      if  (S_Dir == FORWARD)
        ovMesg.orientation.setIsNormal();
      else
        ovMesg.orientation.setIsInnie();

      if  (T_Right_Hang >= S_Right_Hang)
        ovMesg.overlap_type = AS_CONTAINMENT;
      else
        ovMesg.overlap_type =  AS_DOVETAIL;
      ovMesg.ahg = olap -> t_lo;
      ovMesg.bhg = S_Right_Hang - T_Right_Hang;
      ovMesg.quality =  olap -> quality;
      ovMesg.min_offset = olap -> t_lo - (this_diag - olap -> min_diag);
      ovMesg.max_offset = olap -> t_lo + (olap -> max_diag - this_diag);
    }

  assert (ovMesg.min_offset <= ovMesg.ahg && ovMesg.ahg <= ovMesg.max_offset);
  ovMesg.polymorph_ct = 0;

#if  OUTPUT_OVERLAP_DELTAS
  for  (i = 0;  i < olap -> delta_ct;  i ++)
    {
      int  j;

      for  (j = abs (olap -> delta [i]);  j > 0;  j -= AS_LONGEST_DELTA)
        {
          if  (j > AS_LONGEST_DELTA)
            *deltaCursor++ = AS_LONG_DELTA_CODE;
          else{
            *deltaCursor++ = j * Sign (olap -> delta [i]);
          }
        }
    }
#endif
  *deltaCursor = AS_ENDOF_DELTA_CODE;


  if((ovMesg.overlap_type == AS_CONTAINMENT) &&
     (ovMesg.orientation.isOuttie())) {
    // CMM: Regularize the reverse orientated containment overlaps to
    // a common orientation.
    const int ahg = ovMesg.ahg;
    const int bhg = ovMesg.bhg;
    const int min_delta = ovMesg.min_offset - ahg;
    const int max_delta = ovMesg.max_offset - ahg;

    ovMesg.orientation.setIsInnie();
    ovMesg.ahg = -bhg;
    ovMesg.bhg = -ahg;
    ovMesg.min_offset = ovMesg.ahg - max_delta;
    ovMesg.max_offset = ovMesg.ahg - min_delta;
  }

  WA->Total_Overlaps ++;
  if  (ovMesg . bhg <= 0)
    WA->Contained_Overlap_Ct ++;
  else
    WA->Dovetail_Overlap_Ct ++;

  AS_OVS_convertOverlapMesgToOVSoverlap(&ovMesg, WA->overlaps + WA->overlapsLen++);

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



void
Output_Partial_Overlap(AS_IID s_id,
                       AS_IID t_id,
                       Direction_t dir,
                       const Olap_Info_t * p,
                       int s_len,
                       int t_len,
                       Work_Area_t  *WA) {

  int  a, b, c, d;
  char  dir_ch;

  Total_Overlaps ++;

  // Convert to canonical form with s forward and use space-based
  // coordinates
  if  (dir == FORWARD)
    {
      a = p -> s_lo;
      b = p -> s_hi + 1;
      c = p -> t_lo;
      d = p -> t_hi + 1;
      dir_ch = 'f';
    }
  else
    {
      a = s_len - p -> s_hi - 1;
      b = s_len - p -> s_lo;
      c = p -> t_hi + 1;
      d = p -> t_lo;
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

