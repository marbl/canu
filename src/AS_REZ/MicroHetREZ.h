
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

#ifndef AS_REZ_MICROHETREZ_H
#define AS_REZ_MICROHETREZ_H

#include "AS_global.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_fragStorePartition.h"

typedef struct{
  int *part;
  int len;
  int groups;
} Partition_t;

typedef struct{
  int *set;
  int len;
} Marker_t;

// The order of the enum is important:
// UNITIG_IS_UNKOWN is greater than UNITIG_IS_SIMPLE
//
typedef enum {
  UNITIG_IS_SHALLOW,
  UNITIG_IS_SIMPLE, 
  UNITIG_IS_UNKNOWN, 
  UNITIG_IS_REPETITIVE
} UnitigStatus_t;

// The alignment structure. 
// The array in countA can be changed using col_count 
// with a specific marker
//
typedef struct{
  char** ali;  // the actual alignment array
  int* countA; // an array containing the number of 'A' in each column
  int* countC; // an array containing the number of 'C' in each column
  int* countG; // an array containing the number of 'G' in each column
  int* countT; // an array containing the number of 'T' in each column
  int* countDash; // an array containing the number of '-' in each column
  int* countBlank; // an array containing the number of ' ' in each column
  int cols;
  int rows;
  int noOfSegs;
  double **seqErrArray; 
  int hasQuality;
} Alignment_t;





/* AS_REZ_MP_MicroHet_prob() 

   RESULT: The function returns a (double) pvalue (probability) of an
   input unitig being SIMPLE -- meaning, having mismatches due to randomly 
   distributed sequencing errors.

   If the returned value is sufficiently small, the unitig should be treated
   as a likely repeat.

   A return value of 1.0 may indicate that the unitig was not deep enough for
   a meaningful test.

   Some false positives may be induced by polymorphisms; however, the 
   calculation should not be drastically misled by multibase indel 
   polymorphisms.

   INPUT:

   bqarray : an array of size [depth*2]*len of bases and quality values
             in alternative rows, giving a multialignment
   idarray : an array of size depth*len giving the fragment iid of each base
             in the multialignment
   handle  : the fragStore from which locale information for each fragment iid
             will be obtained  (-1 (NULLFRAGSTOREHANDLE) if paritioned store is used.)
   phandle  : the partitioned fragStore from which locale information for each fragment iid
             will be obtained (NULL if a traditional unpartitioned store is used);
   len     : number of columns in the multialignment
   depth   : number of rows in the multialignment
*/

double AS_REZ_MP_MicroHet_prob(char **bqarray,
                               int **idarray,
                               FragStoreHandle handle,
                               tFragStorePartition *phandle,
                               int len,
                               int depth);



//  Used in MicroHetIUM.c, not for general consumption.
//
void compress_shreds_and_null_indels(int c,
                                     int r,
                                     FragStoreHandle frag_store, 
                                     tFragStorePartition *pfrag_store,
                                     char **array,
                                     int **id_array,
                                     int verbose);
Alignment_t *AS_REZ_convert_array_to_alignment(char **ar, int c, int r);
Marker_t *AS_REZ_allocate_marker(int l);


////////////////////////////////////////////////////////////////////////////////

#endif
