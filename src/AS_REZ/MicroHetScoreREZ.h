
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
/*********************************************************************
   CVS_ID:  $Id: MicroHetScoreREZ.h,v 1.3 2005-03-22 19:07:57 jason_miller Exp $
 *********************************************************************/
#ifndef MICROHETSCOREREZ_H
#define MICROHETSCOREREZ_H

#include "MicroHetPartitionsREZ.h"
#include "MicroHetREZ.h"

/* seq_err returns the probability of a sequencing
   error turning character 'given' at position (c,r) 
   into character a->ali[c][r].
   The sum of the conditional probabilities should be se.
   The probability of not changing the character is (1-s_ij) */
double seq_err(Alignment_t* a, char given, int c, int r);

/* this function inspects an alignment and assuming all changes are sequencing
   error returns the average of the most likely explanation for each column */
double guess_seqErr(Alignment_t *a, Marker_t* m, int s, int e);

/* functions for allocating and freeing a Alignment_t */
Alignment_t *allocate_alignment(int c, int r);
void         free_alignment(Alignment_t *ali);

/* heuristic score functions to determine the actual partition */
double value(Partition_t *partition, Alignment_t *a, int ct, Marker_t* m);
double value_two_fixed(Partition_t *partition, Alignment_t *a, int r1, int r2, int start, int end, Marker_t* m);
double col_value(Partition_t *partition, Alignment_t *a, int col, Marker_t* m);


/* this function returns the 1-p, where p
   is the probability that there are at least two
   different characters in a column that occur 
   at least x times, given the sequencing error is SeqErr */
double column_prob(int x, int r, double seqErr);

/* this function computes the probability that for two fixed rows
   in a column the two characters are the same, and that it is different
   from the majority character of that column */
double column_prob_two_fixed(int r, double seqErr);

/* this function returns the expected number of "steps" which can be "saved"
   given the optimal partitioning of a column of data under the null hypothesis
   of a simple alignment */
double expected_savedSteps(int r, double seqErr);
//Define maxfive for use in expected_savedSteps
#define maxfive(v,w,x,y,z) max(v,max(max(w,x),max(y,z)))

int no_col_contributing(Alignment_t *a, int t, int s, int e);
int col_contributing(Alignment_t* a, int col, int thresh);

int col_contributing_two_fixed(Alignment_t *a, int c, int r1, int r2);
int no_col_contributing_two_fixed(Alignment_t *a, int r1, int r2, int s, int e);

/* this function returns the value k, for which the probability
   density becomes greater than 1.0-thresh, p is the success probability
   of the binomial distribution and n the number of trials */
int critical_binomial_value(double p, double thresh, int n);

//void set_xContributing(Alignment_t *a, int s, int e);
//void set_tfContributing(Alignment_t *a, int r1, int r2, int s, int e);

#endif


double binomial_prob(int k,int n, double r);
double binomial(int c1,int n, double q);
double binomial_tail_prob(int n,int k,double p);


//Function to give probability of alignment under null (simple) model
//based only on pairwise mismatch counts
double Pairwise_Only_Prob_Simple(Alignment_t *aln,int begincol,int endcol,Marker_t* m);

//Give p-value associated with a chi-square statistic 
//with df degress of freedom
double chi_to_p(double chisq,int df);
