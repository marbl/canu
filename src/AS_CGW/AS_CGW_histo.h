
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
#ifndef AS_CGW_HISTO_H
#define AS_CGW_HISTO_H

#include <assert.h>

typedef struct {
  int nsamples;
  int sum_bases;
  int min_bases;
  int max_bases;
  int sum_span;
  int min_span;
  int max_span;
  int sum_cedges;
  int min_cedges;
  int max_cedges;
  int sum_wcedges; /* weight of confirmed edges */
  int min_wcedges;
  int max_wcedges;
  int sum_uedges;
  int min_uedges;
  int max_uedges;
  int sum_oedges;
  int min_oedges;
  int max_oedges;
  float sum_ratio;
  float min_ratio;
  float max_ratio;
  int sum_nScaffolds;
  int min_nScaffolds;
  int max_nScaffolds;
  int sum_nScaffoldsConfirmed;
  int min_nScaffoldsConfirmed;
  int max_nScaffoldsConfirmed;
  int sum_nScaffoldsProblem;
  int min_nScaffoldsProblem;
  int max_nScaffoldsProblem;
} ChunkAggregate;

static void SetChunkAggregate(ChunkAggregate *ca, int bases, int span, int uEdges, 
			      int cEdges, int oEdges, int wcEdges,
			      int nScaffolds, int nScaffoldsConfirmed, int nScaffoldsProblem){
  float ratio;
  AssertPtr(ca);
  assert(bases != 0 && span != 0);
  ratio = (float)bases/(float)span;
  ca->nsamples = 1;
  ca->sum_bases = ca->min_bases = ca->max_bases = bases;
  ca->sum_span = ca->min_span = ca->max_span = span;
  ca->sum_uedges = ca->min_uedges = ca->max_uedges = uEdges;
  ca->sum_cedges = ca->min_cedges = ca->max_cedges = cEdges;
  ca->sum_oedges = ca->min_oedges = ca->max_oedges = oEdges;
  ca->sum_wcedges = ca->min_wcedges = ca->max_wcedges = wcEdges;
  ca->sum_ratio = ca->min_ratio = ca->max_ratio = ratio;
  ca->sum_nScaffoldsProblem = ca->min_nScaffoldsProblem = ca->max_nScaffoldsProblem = nScaffoldsProblem;
  ca->sum_nScaffoldsConfirmed = ca->min_nScaffoldsConfirmed = ca->max_nScaffoldsConfirmed = nScaffoldsConfirmed;
  ca->sum_nScaffolds = ca->min_nScaffolds = ca->max_nScaffolds = nScaffolds;
}

extern void printChunks(FILE *fout,
		 DataType *d,
		 DataType *s,
		 DataType *a);

extern void aggregateChunks(DataType *aa,int i,DataType *bb) ;

extern void printChunkAggregate(FILE *fout, DataType *aa);

#endif
