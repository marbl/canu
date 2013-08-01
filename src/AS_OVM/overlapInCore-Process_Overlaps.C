
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


//  Find and output all overlaps between strings in  stream
//  and those in the global hash table.   (* WA)  has the
//  data structures that used to be global.

void
Process_Overlaps (gkStream *stream, Work_Area_t * WA){
  char  Frag [AS_READ_MAX_NORMAL_LEN + 1];
  char  quality [AS_READ_MAX_NORMAL_LEN + 1];
  AS_IID  Curr_String_Num;
  uint32  last_old_frag_read;
  int  frag_status;
  int  Len;

  WA->overlapsLen                = 0;

  WA->Total_Overlaps             = 0;
  WA->Contained_Overlap_Ct       = 0;
  WA->Dovetail_Overlap_Ct        = 0;

  WA->Kmer_Hits_Without_Olap_Ct  = 0;
  WA->Kmer_Hits_With_Olap_Ct     = 0;
  WA->Multi_Overlap_Ct           = 0;

  while  ((frag_status = Read_Next_Frag (Frag,
                                         quality, stream, &WA -> myRead,
                                         & last_old_frag_read,
                                         minLibToRef,
                                         maxLibToRef))) {

    if  (frag_status == DELETED_FRAG)
      continue;

    Curr_String_Num = WA->myRead.gkFragment_getReadIID ();

    Len = strlen (Frag);

    if  (Len < Min_Olap_Len)
      continue;

    Find_Overlaps (Frag, Len, quality, Curr_String_Num, FORWARD, WA);

    reverseComplement (Frag, quality, Len);

    Find_Overlaps (Frag, Len, quality, Curr_String_Num, REVERSE, WA);
  }


  pthread_mutex_lock (& Write_Proto_Mutex);

  for (int zz=0; zz<WA->overlapsLen; zz++)
    AS_OVS_writeOverlap(Out_BOF, WA->overlaps + zz);

  WA->overlapsLen = 0;

  Total_Overlaps            += WA->Total_Overlaps;
  Contained_Overlap_Ct      += WA->Contained_Overlap_Ct;
  Dovetail_Overlap_Ct       += WA->Dovetail_Overlap_Ct;

  Kmer_Hits_Without_Olap_Ct += WA->Kmer_Hits_Without_Olap_Ct;
  Kmer_Hits_With_Olap_Ct    += WA->Kmer_Hits_With_Olap_Ct;
  Multi_Overlap_Ct          += WA->Multi_Overlap_Ct;

  pthread_mutex_unlock (& Write_Proto_Mutex);
}


