
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
#ifndef AS_UTL_HISTO_H
#define AS_UTL_HISTO_H

typedef void DataType;
typedef void (*AggregateFn)(DataType *a,int ib,DataType *b);
typedef void (*PrintFn)(FILE *fout,
			DataType *,
			DataType *,
			DataType *);
typedef void (*PrintAg)(FILE *fout,
			DataType *);


typedef struct {
  int        cnt;     /* # of data points */
  int        nsample; /* size of sample array */
  /* Number of samples collected before bucketing,
     and number of buckets thereafter */
  int   *thescore;       /* score array */
  int    filled;  /* Have NSAMPLE samples been collected?
		     Are we in bucket mode? */
  int   sync;    /* value to sync lower bucket bounds to */
  int   low;     /* Lowest bucket score, an integer */
  int   hgh;     /* Highest bucket score */
  int   min;
  int   max;
  int   bucket_width;    /* number of integral score values per bucket */
  int        nbucket;  /* size of bucket array */
  int        *bucket_cnt;  /* bucket array */
  int        *bucket_min;
  int        *bucket_max;
  /* */
#if 0
  int        extended;
#endif
  int        logarithmic;
  int        datasize;
  DataType   *temp_data;
  DataType   *scan_data;
  DataType   *aggr_data;
  DataType   *sample_data; /* sample_data array */
  DataType   *bucket_data;
#ifdef USE_CUSTOM_INDEXING
  DataType   *(*indexdata)(DataType *b,int ib);
  void       (*setdata)(DataType *a,int ib,DataType *b);
#endif
  AggregateFn aggregate;
  PrintFn     printdata;
  PrintAg     printAg;
} HISTOGRAM;


/* Create an initially empty histogram where buckets will be
   selected so that "sync" will be the low end of some bucket.  */

extern HISTOGRAM *create_histogram(
			    int nsample,
			    int nbucket,
			    int sync,
			    int logarithmic);


extern HISTOGRAM *create_extended_histogram(
			    int nsample,
			    int nbucket,
			    int sync,
			    int logarithmic,
			    int datasize,
			    AggregateFn aggregate,
			    PrintFn printdata,
			    PrintAg printAg);


extern void extend_histogram(
		      HISTOGRAM *h,
		      int dataSize,
#ifdef USE_CUSTOM_INDEXING
		      DataType * (*indexdata)(DataType *b,int ib),
		      void (*setdata)(DataType *a,int ib,DataType *b),
#endif
		      void (*aggregate)(DataType *a,int ib,DataType *b),
		      void (*printdata)(FILE *fout,
					DataType *,
					DataType *,
					DataType *),
			    PrintAg printAg
	);

/* Add data point "data" to histogram "h" */
extern void add_to_histogram(HISTOGRAM *, int thescore, DataType *);

/* Print a histogram of "h" with approximately "rez" buckets and with
   the  ":" of the output in column "indent".  If "rez" = 0 then histogram
   is printed in the finest bucket resolution possible.  (N.B. A call
   to print histogram has the effect of shifting it to bucket mode
   regardless of how many samples have been added thus far.)  The
   function returns the field width used to print the frequency (in
   case one would like to add some aligned columns after the output). */

extern void print_histogram(FILE *,HISTOGRAM *, int rez, int indent);

/* Return the average of the data points in histogram "h" */

extern double histogram_avg(HISTOGRAM *);

/* Return the standard deviation of the data points in histogram "h" */

extern double histogram_stdev(HISTOGRAM *);

/* Free the data structure for histogram "h" */

extern void free_histogram(HISTOGRAM *);

#endif
