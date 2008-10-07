
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
#ifndef CA_ALN_SCAFCOMP
#define CA_ALN_SCAFCOMP

typedef struct {
  int gap_length;  /* Length of the inter-contig gap */
  float gap_var;   /* Variance of inter-contig gap */
} Scaffold_Gap;

typedef struct {
  int length;      /* Length of contig */
  int lft_end;     /* Position of 5p end, relative to scaffold origin */
  int insert_pnt;  /* Location of contig seq in compacted sequence */
} Scaffold_Tig;

typedef struct {
  int            length;      /* Length of scaffold, including mean gap len */
  char          *packed_seq;  /* 0-terminated compacted scaffold sequence */
  int            num_gaps;    /* Number of inter-contig gaps in scaffold */
  Scaffold_Gap  *gaps;        /* gaps[i] for i in [0,num_gaps-1] is the  */
                              /* description of the (i+1)'st inter-contig gap */
  Scaffold_Tig  *ctgs;        /* ctgs[i] for i in [0,num_gaps] is the  */
                              /* description of the (i+1)'st contig */
} Scaffold;


typedef struct _Segment_tag {
  struct _Segment_tag *next;     /* NULL-terminated linked-list of said */
  int                  a_contig; /* A contig containing overlap */
  int                  b_contig; /* B contig containing overlap */
  int                  alow;     /* A[alow,ahgh] is involved in the overlap */
  int                  ahgh;
  int                  blow;     /* B[blow,bhgh] is involved in the overlap */
  int                  bhgh;
  Local_Overlap       *overlap;  /* Overlap between the two contigs */
} Segment;


/* find overlap alignment between two scaffolds at given sigma */

Segment *Align_Scaffold(Segment *seglist, int numsegs, int varwin,
			Scaffold *AF, Scaffold *BF, int *best,
			int bandbeg, int bandend);

#endif





