
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
static char CM_ID[] = "$Id: AS_CGW_histo.c,v 1.6 2006-11-14 17:52:14 eliv Exp $";
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "AS_global.h"

#include "AS_UTL_histo.h"
#include "AS_CGW_histo.h"

#define AGGREGATE_FIELD_VALUE(field,value) \
  a[i].sum_##field += value ;\
  a[i].min_##field  = min(a[i].min_##field,value);\
  a[i].max_##field  = MAX(a[i].max_##field,value);

#define AGGREGATE_FIELD(field) \
  a[i].sum_##field += b->sum_##field ;\
  a[i].min_##field  = min(a[i].min_##field,b->min_##field);\
  a[i].max_##field  = MAX(a[i].max_##field,b->max_##field);

void printChunkAggregate(FILE *fout, DataType *aa){
  ChunkAggregate * a = (ChunkAggregate *) aa;
  assert (a==aa);

  fprintf(fout,"* nsamples = %d sum_bases = %d sum_span = %d sum_cedges = %d, min_ratio %g max_ratio %g\n",
	  a->nsamples, a->sum_bases, a->sum_span, a->sum_cedges, a->min_ratio, a->max_ratio);
  fprintf(fout,"* nscaff = %d nscaffC = %d nscaffP = %d\n",
	  a->sum_nScaffolds, a->sum_nScaffoldsConfirmed, a->sum_nScaffoldsProblem);
}

void aggregateChunks(DataType *aa,int i,DataType *bb) {
  ChunkAggregate * a = (ChunkAggregate *) aa;
  ChunkAggregate * b = (ChunkAggregate *) bb;
  double ratio;
  assert (a==aa);
  assert (b==bb);

  if(b->nsamples > 0){
    if(b->min_bases && b->min_span)
      ratio = b->min_span / b->min_bases;
    else{
      fprintf(stderr,"* Ratio 0 : ");
      ratio = 0.0;
    }
#if 0
    fprintf(stderr,"* aggChunks:\n");
    printChunkAggregate(stderr, b);
    printChunkAggregate(stderr, a + i);
#endif
    if(a[i].nsamples){
      a[i].nsamples  += b->nsamples;
      AGGREGATE_FIELD(bases)
        AGGREGATE_FIELD(span)
        AGGREGATE_FIELD(oedges)
        AGGREGATE_FIELD(cedges)
        AGGREGATE_FIELD(uedges)
        AGGREGATE_FIELD(nScaffolds)
        AGGREGATE_FIELD(nScaffoldsConfirmed)
        AGGREGATE_FIELD(nScaffoldsProblem)
        if(b->sum_cedges > 0){
          if(a[i].sum_cedges > 0){
            AGGREGATE_FIELD(wcedges)
              }else{
                a[i].sum_wcedges = a[i].min_wcedges = a[i].max_wcedges  = b->sum_wcedges;
              }
        }
      AGGREGATE_FIELD(ratio)
        }else{
          a[i] = *b;
        }
#if 0
    printChunkAggregate(stderr, a + i);
#endif
  }
}

#define PRINT_FIELD(field) \
	  data->sum_##field, \
	  scan_data->sum_##field,  \
	  (aggr_data->sum_##field == 0?0:((float)scan_data->sum_##field)/(((float)aggr_data->sum_##field))), \
	  data->min_##field , \
	  (data->nsamples > 0 ? data->sum_##field / data->nsamples : 0), \
	  data->max_##field

#define PRINT_FIELD_ALT(field) \
	  data->sum_##field, \
	  data->min_##field , \
	  (data->nsamples > 0 ? data->sum_##field / data->nsamples : 0), \
	  data->max_##field


void printChunks(FILE *fout,
		 DataType *d,
		 DataType *s,
		 DataType *a)
{
  ChunkAggregate * data = (ChunkAggregate *) d;
  ChunkAggregate * scan_data = (ChunkAggregate *) s;
  ChunkAggregate * aggr_data = (ChunkAggregate *) a;
  assert (data==d);
  assert (scan_data==s);
  assert (aggr_data==a);

  /* Assert that the correct number of items are in the bin? */
  fprintf(fout,"\n"
	  "  bases    % 10d % 10d % 10.5f   % 10d % 10d % 10d\n"
	  "  span     % 10d % 10d % 10.5f   % 10d % 10d % 10d\n"
	  "  nscaff   % 10d % 10d % 10.5f   % 10d % 10d % 10d\n"
	  "  nscaffC  % 10d % 10d % 10.5f   % 10d % 10d % 10d\n"
	  "  nscaffP  % 10d % 10d % 10.5f   % 10d % 10d % 10d\n"
	  "  cedges   % 10d % 10d % 10.5f   % 10d % 10d % 10d\n"
	  "  wcedges  % 10d                          % 10d % 10d % 10d\n"
	  "  uedges   % 10d % 10d % 10.5f   % 10d % 10d % 10d\n"
	  "  oedges   % 10d % 10d % 10.5f   % 10d % 10d % 10d\n"
	  "  avgRatio % 10.4f                         % 10.4f % 10.4f % 10.4f\n"
	  "  ratioAvg                                    % 10.4f % 10.4f % 10.4f\n",
	  PRINT_FIELD(bases),
	  PRINT_FIELD(span),
	  PRINT_FIELD(nScaffolds),
	  PRINT_FIELD(nScaffoldsConfirmed),
	  PRINT_FIELD(nScaffoldsProblem),
	  PRINT_FIELD(cedges),
	  PRINT_FIELD_ALT(wcedges),
	  PRINT_FIELD(uedges),
	  PRINT_FIELD(oedges),
	  PRINT_FIELD_ALT(ratio),
          (float)data->min_bases/(float)data->min_span,
          (float)data->sum_bases/(float)data->sum_span,
          (float)data->max_bases/(float)data->max_span);
}
