
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
#include "AS_UTL_reverseComplement.H"
#include <pthread.h>

//  Find and output all overlaps between strings in store and those in the global hash table.
//  This is the entry point for each compute thread.

void *
Process_Overlaps(void *ptr){
  Work_Area_t  *WA = (Work_Area_t *)ptr;

  WA->overlapsLen                = 0;

  WA->Total_Overlaps             = 0;
  WA->Contained_Overlap_Ct       = 0;
  WA->Dovetail_Overlap_Ct        = 0;

  WA->Kmer_Hits_Without_Olap_Ct  = 0;
  WA->Kmer_Hits_With_Olap_Ct     = 0;
  WA->Multi_Overlap_Ct           = 0;

#warning gkReadData wants to be persistent, but isnt implemented as such.
  gkReadData    readData;

  char         *bases = new char [AS_MAX_READLEN + 1];
  char         *quals = new char [AS_MAX_READLEN + 1];

  for (uint32 fi=WA->frag_segment_lo; fi<=WA->frag_segment_hi; fi++) {

    //  Load sequence/quality data
    //  Duplicated in Build_Hash_Index()

    gkRead   *read = WA->gkpStore->gkStore_getRead(fi);

    if (read->gkRead_isDeleted())
      continue;

    if ((read->gkRead_libraryID() < G.minLibToRef) ||
        (read->gkRead_libraryID() > G.maxLibToRef))
      continue;

    uint32 bgn = read->gkRead_clearRegionBegin();
    uint32 end = read->gkRead_clearRegionEnd();
    uint32 len = read->gkRead_clearRegionLength();

    if (G.Ignore_Clear_Range == true) {
      bgn = 0;
      end = read->gkRead_sequenceLength();
      len = end;
    }

    if (len < G.Min_Olap_Len)
      continue;

    WA->gkpStore->gkStore_loadReadData(read, &readData);

    char   *seqptr   = readData.gkReadData_getSequence()  + bgn;
    char   *qltptr   = readData.gkReadData_getQualities() + bgn;

    for (uint32 i=0; i<len; i++) {
      bases[i] = tolower(seqptr[i]);
      quals[i] = qltptr[i] - QUALITY_BASE_CHAR;
    }

    bases[len] = 0;
    quals[len] = 0;

    //  Generate overlaps.

    Find_Overlaps(bases, len, quals, read->gkRead_readID(), FORWARD, WA);

    reverseComplement(bases, quals, len);

    Find_Overlaps(bases, len, quals, read->gkRead_readID(), REVERSE, WA);
  }

  delete [] bases;
  delete [] quals;

  pthread_mutex_lock(& Write_Proto_Mutex);

  for (int zz=0; zz<WA->overlapsLen; zz++)
    Out_BOF->writeOverlap(WA->overlaps + zz);

  WA->overlapsLen = 0;

  Total_Overlaps            += WA->Total_Overlaps;
  Contained_Overlap_Ct      += WA->Contained_Overlap_Ct;
  Dovetail_Overlap_Ct       += WA->Dovetail_Overlap_Ct;

  Kmer_Hits_Without_Olap_Ct += WA->Kmer_Hits_Without_Olap_Ct;
  Kmer_Hits_With_Olap_Ct    += WA->Kmer_Hits_With_Olap_Ct;
  Multi_Overlap_Ct          += WA->Multi_Overlap_Ct;

  pthread_mutex_unlock(& Write_Proto_Mutex);

  return(ptr);
}


