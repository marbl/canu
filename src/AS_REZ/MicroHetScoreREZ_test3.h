
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
   CVS_ID:  $Id: MicroHetScoreREZ_test3.h,v 1.5 2006-11-14 17:52:18 eliv Exp $
 *********************************************************************/
#ifndef AS_REZ_MICROHETSCOREREZ_H
#define AS_REZ_MICROHETSCOREREZ_H

#include "MicroHetPartitionsREZ_test3.h"
#include "MicroHetREZ_test3.h"

/* seq_err returns the probability of a sequencing
   error turning character 'given' at position (c,r) 
   into character a->ali[c][r].
   The sum of the conditional probabilities should be se.
   The probability of not changing the character is (1-s_ij) */
double seq_err(Alignment_t* a, char given, int c, int r);

/* this function inspects an alignment and assuming all changes are sequencing
   error returns the average of the most likely explanation for each column */
double AS_REZ_guess_seqErr(Alignment_t *a, Marker_t* m, int s, int e);

/* functions for allocating and freeing a Alignment_t */
Alignment_t *AS_REZ_allocate_alignment(int c, int r);
void         AS_REZ_free_alignment(Alignment_t *ali);

/* this function returns the expected number of "steps" which can be "saved"
   given the optimal partitioning of a column of data under the null hypothesis
   of a simple alignment */
double AS_REZ_expected_savedSteps(int r, double seqErr);
//Define maxfive for use in expected_savedSteps
#define maxfive(v,w,x,y,z) MAX(v,MAX(MAX(w,x),MAX(y,z)))

double AS_REZ_binomial_prob(int k,int n, double r);
double AS_REZ_binomial(int c1,int n, double q);
double AS_REZ_binomial_tail_prob(int n,int k,double p);

double AS_REZ_fournomial(double seqErr, int c1, int c2, int c3, int c4, int n);
#endif


