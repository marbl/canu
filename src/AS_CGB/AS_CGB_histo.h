
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

#ifndef AS_CGB_HISTO_INCLUDE
#define AS_CGB_HISTO_INCLUDE

static const char *rcsid_AS_CGB_HISTO_INCLUDE = "$Id: AS_CGB_histo.h,v 1.7 2008-10-08 22:02:54 brianwalenz Exp $";

#include "AS_CGB_all.h"

typedef void HistoDataType;
typedef void Histogram_t;

//  Create an initially empty histogram where buckets will be selected
//  so that "sync" will be the low end of some bucket.

Histogram_t *create_histogram(int nsample,
                              int nbucket,
                              //int sync0,
                              int filled,
                              int logarithmic);


void extend_histogram(Histogram_t *h,
		      size_t nbytes,
		      HistoDataType * (*indexdata)(HistoDataType *b,int ib),
		      void (*setdata)(HistoDataType *a,int ib,HistoDataType *b),
		      void (*aggregate)(HistoDataType *a,int ib,HistoDataType *b),
		      void (*printdata)(FILE *fout,
					HistoDataType *,
					HistoDataType *,
					HistoDataType *));

void free_histogram(Histogram_t *);
void add_to_histogram(Histogram_t *, int thescore, HistoDataType *);
void print_histogram(FILE *,Histogram_t *, int rez, int indent);
double histogram_avg(Histogram_t *);
double histogram_stdev(Histogram_t *);



typedef struct {
  int nsamples;
  int min_frags;
  int sum_frags;
  int max_frags;
  int min_rs_frags; // randomly sampled fragments
  int sum_rs_frags;
  int max_rs_frags;
  int min_nr_frags; // non-randomly sampled fragments
  int sum_nr_frags;
  int max_nr_frags;
  int64  min_bp;
  int64  sum_bp;
  int64  max_bp;
  int64  min_rho;
  int64  sum_rho;
  int64  max_rho;
  int min_arrival;
  int sum_arrival;
  int max_arrival;
  int min_discr;
  int sum_discr;
  int max_discr;
} MyHistoDataType;


static HistoDataType *myindexdata (HistoDataType *aa,int i) {
  MyHistoDataType *a = (MyHistoDataType *)aa;
  return (HistoDataType *)&(a[i]);
}

static void mysetdata (HistoDataType *aa,int i,HistoDataType *bb) {
  MyHistoDataType *a = (MyHistoDataType *)aa;
  MyHistoDataType *b = (MyHistoDataType *)bb;
  a[i] = *b;
}

static void myaggregate (HistoDataType *aa,int i,HistoDataType *bb) {
  MyHistoDataType *a = (MyHistoDataType *)aa;
  MyHistoDataType *b = (MyHistoDataType *)bb;
  a[i].nsamples  += b->nsamples;
  a[i].sum_frags += b->sum_frags;
  a[i].min_frags  = MIN(a[i].min_frags,b->min_frags);
  a[i].max_frags  = MAX(a[i].max_frags,b->max_frags);
  // Just the randomly sampled fragments:
  a[i].sum_rs_frags += b->sum_rs_frags;
  a[i].min_rs_frags  = MIN(a[i].min_rs_frags,b->min_rs_frags);
  a[i].max_rs_frags  = MAX(a[i].max_rs_frags,b->max_rs_frags);
  // Just the non-randomly sampled fragments:
  a[i].sum_nr_frags += b->sum_nr_frags;
  a[i].min_nr_frags  = MIN(a[i].min_nr_frags,b->min_nr_frags);
  a[i].max_nr_frags  = MAX(a[i].max_nr_frags,b->max_nr_frags);
  a[i].sum_bp    += b->sum_bp;
  a[i].min_bp     = MIN(a[i].min_bp,b->min_bp);
  a[i].max_bp     = MAX(a[i].max_bp,b->max_bp);
  a[i].sum_rho   += b->sum_rho;
  a[i].min_rho    = MIN(a[i].min_rho,b->min_rho);
  a[i].max_rho    = MAX(a[i].max_rho,b->max_rho);
  a[i].sum_arrival += b->sum_arrival;
  a[i].min_arrival  = MIN(a[i].min_arrival,b->min_arrival);
  a[i].max_arrival  = MAX(a[i].max_arrival,b->max_arrival);
  a[i].sum_discr   += b->sum_discr;
  a[i].min_discr    = MIN(a[i].min_discr,b->min_discr);
  a[i].max_discr    = MAX(a[i].max_discr,b->max_discr);
}

static void myprintdata (FILE *fout,
                         HistoDataType *d,
                         HistoDataType *s,
                         HistoDataType *a) {
  MyHistoDataType *data      = (MyHistoDataType *) d;
  MyHistoDataType *scan_data = (MyHistoDataType *) s;
  MyHistoDataType *aggr_data = (MyHistoDataType *) a;

  fprintf(fout,"\n"
	  "  frags   %10d %10d %10.5f   %6d %6d %6d\n"
	  "  rs frag %10d %10d %10.5f   %6d %6d %6d\n"
	  "  nr frag %10d %10d %10.5f   %6d %6d %6d\n"
	  "  bases   %10"F_S64P " %10"F_S64P " %10.5f   %6"F_S64P" %6"F_S64P" %6"F_S64P"\n"
	  "  rho     %10"F_S64P " %10"F_S64P " %10.5f   %6"F_S64P" %6"F_S64P" %6"F_S64P"\n"
	  "  arrival %10d %10d %10.5f   %6d %6d %6d\n"
	  "  discr   %10d %10d %10.5f   %6d %6d %6d\n",
	  /* data->nsamples, */

	  data->sum_frags,
	  scan_data->sum_frags,
	  ( aggr_data->sum_frags > 0 ? ((float)scan_data->sum_frags)/ ((float)aggr_data->sum_frags) : 0.),
	  data->min_frags,
	  (data->nsamples > 0 ? data->sum_frags / data->nsamples : 0),
	  data->max_frags,

	  data->sum_rs_frags,
	  scan_data->sum_rs_frags,
	  ( aggr_data->sum_rs_frags > 0 ? ((float)scan_data->sum_rs_frags)/ ((float)aggr_data->sum_rs_frags) : 0.),
	  data->min_rs_frags,
	  (data->nsamples > 0 ? data->sum_rs_frags / data->nsamples : 0),
	  data->max_rs_frags,

	  data->sum_nr_frags,
	  scan_data->sum_nr_frags,
	  ( aggr_data->sum_nr_frags > 0 ? ((float)scan_data->sum_nr_frags)/ ((float)aggr_data->sum_nr_frags) : 0.),
	  data->min_nr_frags,
	  (data->nsamples > 0 ? data->sum_nr_frags / data->nsamples : 0),
	  data->max_nr_frags,

	  data->sum_bp,
	  scan_data->sum_bp,
	  ((float)scan_data->sum_bp)/((float)aggr_data->sum_bp),
	  data->min_bp,
	  (data->nsamples > 0 ? data->sum_bp / data->nsamples : 0),
	  data->max_bp,

	  data->sum_rho,
	  scan_data->sum_rho,
	  ( aggr_data->sum_rho > 0 ? ((float)scan_data->sum_rho)/ ((float)aggr_data->sum_rho) : 0.),
	  data->min_rho,
	  (data->nsamples > 0 ? data->sum_rho/data->nsamples : 0),
	  data->max_rho,

	  data->sum_arrival,
	  scan_data->sum_arrival,
	  ( aggr_data->sum_arrival > 0 ? ((float)scan_data->sum_arrival)/ ((float)aggr_data->sum_arrival) : 0.),
	  data->min_arrival,
	  (data->nsamples > 0 ? data->sum_arrival / data->nsamples : 0),
	  data->max_arrival,

	  data->sum_discr,
	  scan_data->sum_discr,
	  ( aggr_data->sum_discr > 0 ? ((float)scan_data->sum_discr)/ ((float)aggr_data->sum_discr) : 0.),
	  data->min_discr,
	  (data->nsamples > 0 ? data->sum_discr / data->nsamples : 0),
	  data->max_discr);
}


#endif /*AS_CGB_HISTO_INCLUDE*/
