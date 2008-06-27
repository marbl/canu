
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
#include "AS_PER_gkpStore.h"

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



double       AS_REZ_MP_MicroHet_prob(char **bqarray,
                                     int  **idarray,
                                     GateKeeperStore  *handle,
                                     int len,
                                     int depth);



//  Used in MicroHetIUM.c and AS_CNS/colCorr_CNS.c, not for general consumption.
//
void            AS_REZ_count_columns(Alignment_t* a, Marker_t* m);
UnitigStatus_t  AS_REZ_test_MPsimple(Alignment_t *ali, double thresh, Marker_t* m,
                                     int start, int end,double *pval);
void            AS_REZ_compress_shreds_and_null_indels(int c,
                                                       int r,
                                                       GateKeeperStore *gkpstore,
                                                       char **array,
                                                       int **id_array,
                                                       int verbose);
Alignment_t    *AS_REZ_convert_array_to_alignment(char **ar, int c, int r);

Marker_t       *AS_REZ_allocate_marker(int l);
void            AS_REZ_free_marker(Marker_t *m);

void            AS_REZ_print_alignment(Alignment_t *a,  int w);
void            AS_REZ_free_alignment(Alignment_t* a);




////////////////////////////////////////////////////////////////////////////////

#endif
