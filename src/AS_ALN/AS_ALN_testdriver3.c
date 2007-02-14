
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
/* MODULE FOR READING FASTA SEQUENCES, COMPUTING OVERLAPS, AND THE COLUMN
   SETS FOR THE CORRELATED-DIFFERENCES DETECTOR.
*/

#undef INPUT

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AS_global.h"
#include "AS_ALN_aligners.h"

//void *safe_malloc(int size)	/* Guarded malloc utility */
//{ void *new;
//
//  new = (void *) malloc(size);
//  if (new == NULL)
//    { fprintf(stderr,"Out of memory\n");
//      exit (1);
//    }
//  return (new);
//}

#define SEQLEN 500

void GenPair(char **a, char **b)
{ static char A[SEQLEN+1], B[SEQLEN+1];
  static DNA[] = { 'a', 'c', 'g', 't' };
  int i;

  for (i = 0; i < SEQLEN; i++)
    A[i] = DNA[(int) (4.*drand48())];
  A[SEQLEN] = '\0';
  for (i = 0; i < SEQLEN; i++)
    B[i] = DNA[(int) (4.*drand48())];
  B[SEQLEN] = '\0';
  *a = A;
  *b = B;
}

main(int argc, char *argv[])
{ InternalFragMesg  A, B;
  OverlapMesg  *O;
  char qlty1[AS_READ_MAX_LEN+1], qlty2[AS_READ_MAX_LEN+1];
  int j, where;

  srand(getpid());
  A.quality = qlty1;
  B.quality = qlty2;
  for (j = 0; j < 1000; j++)
    { GenPair(&(A.sequence),&(B.sequence));
      A.iaccession = A.eaccession = j+1;
      B.iaccession = B.eaccession = j+2;
      printf("%4d\n",j);
      O = DP_Compare_AS(&A,&B,-strlen(B.sequence),strlen(A.sequence),
                        0,.13,1e-6,40,AS_FIND_ALIGN,&where);
    }

  exit (0);
}
