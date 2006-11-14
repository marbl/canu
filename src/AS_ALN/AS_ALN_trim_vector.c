
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
/*********************************************************************/
// headers
/*********************************************************************/
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// project headers
#include "AS_global.h"
#include "AS_ALN_aligners.h"

#define MAX_END_GAP   1

void TrimFivePrimeThreePrime( InternalFragMesg * vector,
                              InternalFragMesg * frag,
                              cds_float32 variation,
                              cds_int32 min_length )
{
  OverlapMesg * ovl;
  int where;

  ovl =
    DP_Compare_AS( vector, frag,
                   -(frag->clear_rng.end - frag->clear_rng.bgn),
                   vector->clear_rng.end,
                   0, 1e-1, (double) variation, min_length,
                   AS_FIND_OVERLAP, &where );

  /* ovl must be non-NULL &
     if dovetail, 5' of vector overlaps 3' of fragment
       frag   = ovl->aifrag:   5'------------>3'
       vector = ovl->bifrag:          5'----------->3'
     else if containment,
       if vector contains fragment - trim to nothing
       or if fragment contains vector near 3' end of fragment, trim
  */
  if( ovl )
  {
    switch( ovl->overlap_type )
    {
      case AS_DOVETAIL:
        if( ovl->aifrag == frag->iaccession &&
            ovl->orientation == AS_NORMAL )
          frag->clear_rng.end -=
            (frag->clear_rng.end - frag->clear_rng.bgn) - ovl->ahg;
        break;
      case AS_CONTAINMENT:
        if( ovl->aifrag == vector->iaccession )
          frag->clear_rng.end = frag->clear_rng.bgn;
        else if( ovl->bhg <= -min_length )
          frag->clear_rng.end += ovl->bhg;
        return;
        break;
      default:
        assert(0);
        break;
    }
  }
}

void TrimThreePrimeFivePrime( InternalFragMesg * vector,
                              InternalFragMesg * frag,
                              cds_float32 variation,
                              cds_int32 min_length )
{
  OverlapMesg * ovl;
  int where;

  ovl =
    DP_Compare_AS( vector, frag,
                   -(frag->clear_rng.end - frag->clear_rng.bgn),
                   vector->clear_rng.end,
                   0, 1e-1, (double) variation, min_length,
                   AS_FIND_OVERLAP, &where );

  /* ovl must be non-NULL &
     if dovetail, 3' of vector overlaps 5' of fragment
       vector = ovl->aifrag:   5'------------>3'
       frag   = ovl->bifrag:          5'----------->3'
     else if containment,
       if vector contains fragment - trim to nothing
       or if fragment contains vector near 5' end of fragment, trim
  */
  if( ovl )
  {
    switch( ovl->overlap_type )
    {
      case AS_DOVETAIL:
        if( ovl->bifrag == frag->iaccession &&
            ovl->orientation == AS_NORMAL )
          frag->clear_rng.bgn +=
            (frag->clear_rng.end - frag->clear_rng.bgn) - ovl->bhg;
        break;
      case AS_CONTAINMENT:
        if( ovl->aifrag == vector->iaccession )
          frag->clear_rng.end = frag->clear_rng.bgn;
        else if( ovl->ahg <= min_length )
          frag->clear_rng.bgn += ovl->ahg;
        return;
        break;
      default:
        assert(0);
        break;
    }
  }
}

void TrimThreePrimeComplementThreePrime( InternalFragMesg * vector,
                                         InternalFragMesg * frag,
                                         cds_float32 variation,
                                         cds_int32 min_length )
{
  OverlapMesg * ovl;
  int where;

  ovl =
    DP_Compare_AS( vector, frag,
                   -(frag->clear_rng.end - frag->clear_rng.bgn),
                   vector->clear_rng.end,
                   1, 1e-1, (double) variation, min_length,
                   AS_FIND_OVERLAP, &where );

  /* ovl must be non-NULL &
     if dovetail, 3' of vector overlaps 3' of fragment
       frag   = ovl->aifrag:   5'------------>3'
       vector = ovl->bifrag:         3'<-----------5'
       OR
       vector = ovl->aifrag:   5'------------>3'
       frag   = ovl->bifrag:         3'<-----------5'
     else if containment,
       if vector contains fragment - trim to nothing
       or if fragment contains vector near 3' end of fragment, trim
  */
  if( ovl )
  {
    switch( ovl->overlap_type )
    {
      case AS_DOVETAIL:
        if( ovl->orientation == AS_INNIE )
        {
          if( ovl->aifrag == frag->iaccession )
            frag->clear_rng.end -=
              (frag->clear_rng.end - frag->clear_rng.bgn) - ovl->ahg;
          else
            frag->clear_rng.end -=
              (frag->clear_rng.end - frag->clear_rng.bgn) - ovl->bhg;
        }
        break;
      case AS_CONTAINMENT:
        if( ovl->aifrag == vector->iaccession )
          frag->clear_rng.end = frag->clear_rng.bgn;
        else if( ovl->bhg <= -min_length )
          frag->clear_rng.end += ovl->bhg;
        return;
        break;
      default:
        assert(0);
        break;
    }
  }
}


void TrimFivePrimeComplementFivePrime( InternalFragMesg * vector,
                                       InternalFragMesg * frag,
                                       cds_float32 variation,
                                       cds_int32 min_length )
{
  OverlapMesg * ovl;
  int where;

  ovl =
    DP_Compare_AS( vector, frag,
                   -(frag->clear_rng.end - frag->clear_rng.bgn),
                   vector->clear_rng.end,
                   1, 1e-1, (double) variation, min_length,
                   AS_FIND_OVERLAP, &where );

  /* ovl must be non-NULL &
     if dovetail, 5' of vector overlaps 5' of fragment
       frag   = ovl->aifrag: 3'<-----------5'
       vector = ovl->bifrag:        5'------------>3'
       OR
       vector = ovl->aifrag: 3'<-----------5'
       frag   = ovl->bifrag:        5'------------>3'
     else if containment,
       if vector contains fragment - trim to nothing
       or if fragment contains vector near 3' end of fragment, trim
  */
  if( ovl )
  {
    switch( ovl->overlap_type )
    {
      case AS_DOVETAIL:
        if( ovl->orientation == AS_OUTTIE )
        {
          if( ovl->aifrag == frag->iaccession )
            frag->clear_rng.bgn +=
              (frag->clear_rng.end - frag->clear_rng.bgn) - ovl->ahg;
          else
            frag->clear_rng.bgn +=
              (frag->clear_rng.end - frag->clear_rng.bgn) - ovl->bhg;
        }
        break;
      case AS_CONTAINMENT:
        if( ovl->aifrag == vector->iaccession )
          frag->clear_rng.end = frag->clear_rng.bgn;
        else if( ovl->bhg <= -min_length )
          frag->clear_rng.bgn -= ovl->bhg;
        return;
        break;
      default:
        assert(0);
        break;
    }
  }
}


int DP_Trim_Vector_AS( InternalScreenItemMesg * isn,
                       InternalFragMesg       * ifg,
                       cds_float32              variation,
                       cds_int32                min_length,
                       VectorType               vector_type )
{
  // NOTE: DP_Compare_AS ignores clear range of ifg message
  InternalFragMesg   vector;
  InternalFragMesg   frag;
  cds_int32          temp_bgn_end;
  char               clr_seq[AS_READ_MAX_LEN];

  // make an ifg message for the isn message
  vector.iaccession = isn->iaccession;
  vector.sequence = isn->sequence;
  vector.quality = NULL;
  vector.clear_rng.bgn = 0;
  vector.clear_rng.end = strlen( isn->sequence );

  // DP_Compare_AS works with entire sequence, ignoring clear range
  // make a clear range trimmed sequence for the ifg
  // keep the clear_rng values the same
  memcpy( &frag, ifg, sizeof( InternalFragMesg ) );
  strncpy( clr_seq, &(ifg->sequence[ifg->clear_rng.bgn]),
           ifg->clear_rng.end - ifg->clear_rng.bgn );
  clr_seq[ifg->clear_rng.end - ifg->clear_rng.bgn] = '\0';
  frag.sequence = clr_seq;
  frag.quality = NULL;

  // Check for different vector/fragment relationships
  switch( vector_type )
  {
    case AS_SEQUENCING_VECTOR:
      // NOTE: not tested yet
      return 1;
      // 3' end of vector into 5' end of sequence
      TrimThreePrimeFivePrime( &vector, &frag, variation, min_length );
      ifg->clear_rng.bgn = frag.clear_rng.bgn;
      break;
    case AS_CLONING_VECTOR:
      // NOTE: not tested yet
      return 1;
      // 3' end of vector into 5' end of sequence
      TrimThreePrimeFivePrime( &vector, &frag, variation, min_length );
      temp_bgn_end = frag.clear_rng.bgn;
      frag.clear_rng.bgn = ifg->clear_rng.bgn;
      
      // complement of 5' end of vector into 5' end of sequence
      TrimFivePrimeComplementFivePrime( &vector, &frag,
                                        variation, min_length );
      ifg->clear_rng.bgn = MAX( frag.clear_rng.bgn, temp_bgn_end );
      break;
    case AS_INSERT_VECTOR:
      // 5' end of vector into 3' end of sequence
      TrimFivePrimeThreePrime( &vector, &frag, variation, min_length );
      temp_bgn_end = frag.clear_rng.end;
      frag.clear_rng.end = ifg->clear_rng.end;
      
      // complement of 3' end of vector into 3' end of sequence
      TrimThreePrimeComplementThreePrime( &vector, &frag,
                                           variation, min_length );
      ifg->clear_rng.end = min( frag.clear_rng.end, temp_bgn_end );
      break;
    default:
      // unknown vector type - bug out
      return 1;
      break;
  }
  
  return 0;
}
