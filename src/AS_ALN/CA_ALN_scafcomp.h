
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

/* SCAFFOLD DATA ABSTRACTION */

#if defined(__cplusplus)
extern "C"
{
#endif


#define AHANG_BAND_TEST
  /* whether to restrict Align_Scaffold to a given ahang band */


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

typedef struct {
  char *bigstr;
  int  *parts;
  int  *scafs;
  int   slen, plen;
} AggregateString;

void Complement_Scaffold(Scaffold *S);

/* In place, complement the sequence of scaffold S and adjust its gap
   records accordingly.  A second call returns S to its original state. */

void Free_Scaffold(Scaffold *S);

/* Free the storage for scaffold S */

/* LIST OF SCAFFOLDS DATA ABSTRACTION */

typedef struct _ScafList_tag {
  struct _ScafList_tag *next;     /* NULL-terminated linked list of scafs */
  Scaffold             *scaffold; /* Scaffold for this list element */
} Scaffold_List;

Scaffold_List *Read_Multi_Fasta(char *fname);

/* Read a multi-fasta file of the set of scaffolds in an assembly.
   It is assumed that inter-contig gaps are signaled by lower-case
   'n's (upper-case 'N's denote ambiguous bases) and that the header
   line contains an appropriately formated list of the standard
   deviations of the inter-contig gaps.                               */

void Free_Scaffold_List(Scaffold_List *SL);

/* Free list of scaffolds */

/* LIST OF SCAFFOLD OVERLAPS ABSTRACTION */

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

typedef struct CO_tag {
   Segment       *seg;
   int            best;
   struct CO_tag *trace;
   struct CO_tag *Alink;
   struct CO_tag *Blink;
} COvlps;

typedef struct _ScafLap_tag {
  struct _ScafLap_tag *next;    /* NULL-terminated linked-list of said */
  Scaffold            *ascaf;   /* A scaffold */
  Scaffold            *bscaf;   /* B scaffold */
  int                  asnum;   /* A scaffold number in assembly */
  int                  bsnum;   /* B scaffold number in assembly */
  int                  abase;   /* #(in assembly) of 1st contig in A scaffold */
  int                  bbase;   /* #(in assembly) of 1st contig in B scaffold */
  Overlap             *overlap; /* Overlap between scaffolds (sequence and
                                     trace fields are NULL!)              */
  int                  score;   /* Sum of lengths contig overlaps involved
                                     in this scaffold overlap.         */
  float                erate;   /* Error rate of matching parts. */
  int                  firm;    /* Alignment satisfies scaffold constraints */
  int                  D_delta; /* Width of band containing all contig
                                     overlaps of this scaffold overlap */
  int                  D_var;   /* Sum of all gap variation between
                                     contigs involved in the overlap */
  int                  A_delta; /* Sum of all A-scaffold segments that should
                                     have been but were not aligned     */
  int                  B_delta; /* Sum of all B-scaffold segments that should
                                     have been but were not aligned     */
  Segment             *seglist; /* List of contig overlaps constituting the
                                     scaffold overlap                      */
} Scaffold_Overlap;


typedef struct {
  int ascaf, bscaf;
  int acntg, bcntg;
} Local_Address;

typedef struct {
  int num_locals;
  Local_Segment *locals;
  Local_Address *address;
} Local_Pool;




Scaffold_Overlap *Compare_Scaffolds(Scaffold_List *A, Scaffold_List *B);

/* Produce a list of all scaffold pairs that overlap, with a complete
   description of each overlap (as embodied in the Scaffold_Overlap struct). */

void Print_Scaffold_Overlaps(Scaffold_Overlap *SO);

/* Produce a print out of the overlaps in SO on stdout. */

void Print_Anno_Info(Scaffold_Overlap *SO);

/* Produce a print out of the anno info needed to
   highlight overlaps in SO on stdout.               */

void Free_Scaffold_Overlaps(Scaffold_Overlap *SO);

/* Free a list of scaffold overlaps */

char *ScfCmp_Next_Line(FILE *ifile);

/* Store next line of file in returned string */

void Free_Segments_ScafComp(Segment *seglist);

/* Free a list of segments */

int Link_Horizontal(Scaffold *, Scaffold *, int,
                           int, int, int, int, COvlps *);

/* traverse Pillar (or Beam?) */

int Link_Vertical  (Scaffold *, Scaffold *, int,
                           int, int, int, int, COvlps *);

/* traverse Beam (or Pillar?) */


Segment *Find_All_Overlaps(Scaffold *AF, Scaffold *BF, Local_Pool *pool,
                                  int as, int bs, int ac, int bc,
                                  int *fing, int *segs, int comp);

/* find all overlaps between two lists of scaffolds */


int MSORT(const void *l, const void *r);

/* sort function for matches */


#ifdef AHANG_BAND_TEST
Segment *Align_Scaffold(Segment *seglist, int numsegs, int varwin,
			Scaffold *AF, Scaffold *BF, int *best,
			int bandbeg, int bandend);
#else
Segment *Align_Scaffold(Segment *seglist, int numsegs, int varwin,
			Scaffold *AF, Scaffold *BF, int *best);
#endif

#ifdef AHANG_BAND_TEST
Segment *Align_Scaffold_ala_Aaron(Segment *seglist, int numsegs, int varwin,
			Scaffold *AF, Scaffold *BF, int *best,
			int bandbeg, int bandend);
#else
Segment *Align_Scaffold_ala_Aaron(Segment *seglist, int numsegs, int varwin,
			Scaffold *AF, Scaffold *BF, int *best);
#endif

/* find overlap alignment between two scaffolds at given sigma */

Segment *Compare_Scaffold(Segment *seglist, int numsegs,
			  Scaffold *AF, Scaffold *BF, int *best);

/* look for poor quality overlap between two scaffolds */

Scaffold_Overlap *Analyze_Overlap(Scaffold *AF, Scaffold *BF,
				  Segment *seglist, int score, int comp);

/* compute statistics on an overlap between two scaffolds */

#if defined(__cplusplus)
}
#endif

#endif





