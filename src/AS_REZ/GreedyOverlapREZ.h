
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
   CVS_ID:  $Id: GreedyOverlapREZ.h,v 1.1.1.1 2004-04-14 13:53:18 catmandew Exp $
 *********************************************************************/
#ifndef GREEDY_OVERLAP_H
#define GREEDY_OVERLAP_H

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "ScaffoldGraph_CGW.h"

#define GREEDYDEBUG 0

// A scaling function to have a meaningfull integer QV score
#define SCALE_QV 20


// this function computes the quality for all CI overlaps and
// and tests whether the best is a true overlap or repetitive.
void checkEdgeQuality(GraphCGW_T *ScaffoldGraph, int MinLength, int qt, char* fp);


// The lookup table (onedimensional) so if you want
// the score of a match with quality values 30,40 look it up in MatchTableREZ[30,40]
extern int* MatchTableREZ;
extern int* MismatchTableREZ;

// the probability of a match respectively a mismatch given the quality values
int prob_match(int qa, int qb);
int prob_mismatch(int qa, int qb);

// These functions fill the tables
int *fill_QV_Match_table(int (*nscore)(int qa, int qb));
int *fill_QV_Mismatch_table(int (*nscore)(int qa, int qb));

// tests, whether the overlap between two fragments is normal or repetitive
int repetitive_overlap_sim(InternalFragMesg* IFG1,InternalFragMesg* IFG2, OverlapMesg* olap);

/* This function fills into the preallocated arrays S1,Q1,S2,Q2
   the matching regions and the respective quality values
   of the match between seq1,seq2 with quality values q1,q2 that is given in olap
   id1 and id2 are the iaccession numbers of the first resp. second fragment
*/
int get_match_region(char* S1,char* Q1,char* S2,char* Q2,
		     char* seq1, char* seq2, char* q1, char* q2,
		     int id1 ,int id2,
		     OverlapMesg* olap);

/* computes the average probability with which a mismatch occurs
   given to null-terminated quality value strings */
double get_average_mismatch_prob(char* S1, char* S2, char* Q1, char* Q2, int qt, int *cc, int wdash);

/* simply counts the number of mismatches in two null-terminated strings */
int    get_no_mismatches(char* S1, char* S2, char* Q1, char* Q2, int qt, int wdash);

/* The main quality function. It assigns a value between 0 (good) and 1 (bad)
   to an overlap between two fragments. NOTE that it is also used to compute
   the overlap between ANY two sequences with quality values. It returns the 
   number of matches */
int compute_bayesian_quality(InternalFragMesg* IFG1, InternalFragMesg* IFG2, OverlapMesg* olap, int qt, int* l, FILE* fp);

#endif








