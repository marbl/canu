
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
static char CM_ID[] = "$Id: AS_UTL_histo.c,v 1.2 2004-09-23 20:25:29 mcschatz Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "AS_UTL_histo.h"

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)
#define FALSE 0
#define DataType void

#define DEBUG
#undef DEBUG


/* Histogram library */



#ifndef USE_CUSTOM_INDEXING
static DataType *getData(DataType *data, int elementSize, int indx){
  return ((DataType *)((char *)data + (indx * elementSize)));
}

static void      setData(DataType *data, int elementSize, int indx, DataType *d){
  DataType *ref = getData(data, elementSize, indx);
  //  fprintf(stderr,"*setData %d ==> address 0x%x size %d\n",
  //  indx, ref, elementSize);
  memcpy(ref, d, elementSize);
}
#endif


HISTOGRAM *create_extended_histogram(
			    int nsample,
			    int nbucket,
			    int sync,
			    int logarithmic,
			    int datasize,
			    AggregateFn aggregate,
			    PrintFn printdata,
			    PrintAg printag){
  HISTOGRAM *h = create_histogram(nsample, nbucket, sync, logarithmic);
  assert(datasize > 0);
  extend_histogram(h, datasize, aggregate, printdata, printag);

  return h;
}


HISTOGRAM *create_histogram(
			    int nsample,
			    int nbucket,
			    int sync,
			    int logarithmic)
/* Create an initially empty histogram where buckets will be
   selected so that "sync" will be the low end of some bucket.  */
{ 
  HISTOGRAM *h;
  int *bucket_cnt,*bucket_min,*bucket_max;
  int *thescore;
  h = (HISTOGRAM *) malloc(sizeof(HISTOGRAM));
  thescore   = (int *)calloc(nsample,sizeof(int));
  bucket_cnt = (int *)calloc(nbucket,sizeof(int));
  h->cnt         = 0;
  h->filled      = 0;
  h->sync        = sync;
  h->nsample     = nsample;
  h->nbucket     = nbucket;
  h->thescore    = thescore;
  h->bucket_cnt  = bucket_cnt; 
  h->logarithmic = logarithmic;
#ifdef DEBUG
  fprintf(stderr,"* nsamplle %d nbucket %d bucket_cnt %d log %d\n",
	  h->nsample, h->nbucket, h->bucket_cnt, h->logarithmic);
#endif
  /* */
  h->datasize    = 0;
  h->temp_data   = NULL;
  h->scan_data   = NULL;
  h->aggr_data   = NULL;
  h->sample_data = NULL;
  h->bucket_data = NULL;
#ifdef USE_CUSTOM_INDEXING
  h->indexdata   = NULL;
  h->setdata     = NULL;
#endif
  h->aggregate   = NULL;
  h->printdata   = NULL;
  /* */
  /* #ifdef LOGARITHMIC */
  bucket_min = (int *)calloc(nbucket+10,sizeof(int));
  bucket_max = (int *)calloc(nbucket,sizeof(int));
  h->bucket_min  = bucket_min;
  h->bucket_max  = bucket_max;
  /* #endif */
  return h;
}

void extend_histogram(  
		      HISTOGRAM *h,
		      int datasize,

#ifdef USE_CUSTOM_INDEXING
		      DataType * (*indexdata)(DataType *b,int ib),
		      void (*setdata)(DataType *a,int ib,DataType *b),
#endif
		      void (*aggregate)(DataType *a,int ib,DataType *b),
		      void (*printdata)(FILE *fout,
					DataType *data,
					DataType *scan_data,
					DataType *aggr_data),
		      void (*printAg) (FILE *fout,
				       DataType *data)
	)
{
  DataType *sample_data;
  DataType *bucket_data;
  DataType *temp_data, *scan_data, *aggr_data;

  assert(h != NULL);
  assert(datasize > 0 );
  h->datasize = datasize;
  sample_data = (DataType *) calloc((h->nsample),datasize);
  bucket_data = (DataType *) calloc((h->nbucket + 10),datasize); /* Code is broken */
  temp_data   = (DataType *) calloc(datasize,1);
  scan_data   = (DataType *) calloc(datasize,1);
  aggr_data   = (DataType *) calloc(datasize,1);
  assert(temp_data != NULL);
#ifdef USE_CUSTOM_INDEXING
  assert(indexdata != NULL);
  assert(setdata   != NULL);
#endif
  assert(aggregate != NULL);
  assert(printdata != NULL);

  h->temp_data   = temp_data;
  h->scan_data   = scan_data;
  h->aggr_data   = aggr_data;
  h->sample_data = sample_data;
  h->bucket_data = bucket_data;
#ifdef USE_CUSTOM_INDEXING
  h->indexdata   = indexdata;
  h->setdata     = setdata;
#endif

#if 0
  h->extended    = 1;
#endif
  h->aggregate   = aggregate;
  h->printdata   = printdata;
  h->printAg = printAg;
}

void free_histogram(HISTOGRAM *h)
/* Free the data structure for histogram "h" */
{ 
  assert(h != NULL);
  free(h->thescore);
  free(h->bucket_cnt);
  free(h->sample_data);
  free(h->bucket_data);
  free(h->temp_data);
  free(h->scan_data);
  free(h->aggr_data);
  free(h->bucket_min);
  free(h->bucket_max);
  free((char *) h); 
}

static int bucket_from_score(HISTOGRAM *h, int thescore) {
  int ib;
  h->max = max(thescore,h->max);
  h->min = min(thescore,h->min);
  ib = thescore;
  if(h->logarithmic == 0){
    ib = (thescore - h->low - 1)/h->bucket_width + 1;
  } else  {
    int thesign, decade, divisor, theoffset;
    int middlebucket,ii;
    middlebucket = h->nbucket/2;
    thesign = (thescore < 0 ? -1 : 1 );
    /* decade  = floor(log10(1.*thesign*thescore)); *** IMPORTANT TO AVOID ROUNDING PROBLEMS ***/
    decade=0; for(ii=thesign*thescore; ii>=10; ii /=10 ) { decade++;}
    /* divisor = (int) pow(10.,(double)decade); *** IMPORTANT TO AVOID ROUNDING PROBLEMS ***/
    divisor = 1; for(ii=0;ii<decade;ii++) { divisor *= 10;}
    theoffset = thesign*thescore/divisor;
    ib = thesign*(theoffset + decade*10) + middlebucket;
#ifdef NEVER
    printf("thescore,thesign,decade,theoffset,ib = %5d %3d %5d %5d %5d\n",
	   thescore,thesign,decade,theoffset,ib);
#endif
  }
  /* Place a floor and ceiling to the histogram. */
  ib = max(ib,0);
  ib = min(ib,h->nbucket-1);
  return ib;
}

static void fill(HISTOGRAM *h)
{ 
  /* Converts from sample-mode to bucket-mode */
  int is, ib;
  int low, hgh;

  assert(h != NULL);
  /* At this point h->sample[.] is a data value. */
  if( h->cnt > 0 ) { 
    low = hgh = h->thescore[0];
  } else {
    low = hgh = 0;
  }
  for (is = 1; is < h->cnt; is++) {
    int thescore;
    thescore = h->thescore[is];
    hgh = max(hgh,thescore);
    low = min(low,thescore);
  }
  /* hgh and low are the bounds on the sampled data values. */
  h->hgh = hgh;
  h->low = low;
  h->max = hgh;
  h->min = low;
  h->bucket_width = (h->hgh - h->low - 1)/(h->nbucket) + 1;
#ifdef DEBUG
  printf("in fill: h->low,h->hgh=%d,%d\n",
	 h->low,h->hgh);
#endif

  for (ib = 0; ib < h->nbucket; ib++) {
    h->bucket_cnt[ib] = 0;
  }
  for (is = 0 ; is < h->cnt; is++) {
    int thescore;
    thescore =  h->thescore[is];
    ib = bucket_from_score(h,thescore);
#ifdef DEBUG
    printf("in fill is,score,ib = %d,%d,%d\n",is,thescore,ib);
#endif
    h->bucket_cnt[ib] += 1;
    if( h->bucket_cnt[ib] == 1) {
      h->bucket_min[ib] = thescore;
      h->bucket_max[ib] = thescore;
    } else {
      h->bucket_min[ib] = min(h->bucket_min[ib],thescore);
      h->bucket_max[ib] = max(h->bucket_max[ib],thescore);
    }
    if(h->datasize) {
      DataType *data;
#ifdef DEBUG
      printf("heeee %d %d\n",is,ib);
#endif
#ifdef USE_CUSTOM_INDEXING
      data = (h->indexdata)(h->sample_data,is);
#else
      assert(is >= 0 && is < h->nsample);
      data = getData(h->sample_data,h->datasize, is);

#endif
      if( h->bucket_cnt[ib] == 1) {
#ifdef USE_CUSTOM_INDEXING
	(h->setdata)(h->bucket_data,ib,data);
#else
	setData(h->bucket_data, h->datasize, ib, data);
#endif
      } else {
	(h->aggregate)(h->bucket_data,ib,data);
      }
    }
  }
  h->filled = 1;
}

/* Add data point "data" to histogram "h" */


void add_to_histogram(HISTOGRAM *h,int thescore,DataType *data)
{ 
  assert(h != NULL);
  assert( !(h->datasize) || (data != NULL) );

  if(h->datasize) {
    if( h->cnt == 0) {
#ifdef USE_CUSTOM_INDEXING
      (h->setdata)(h->aggr_data,0,data);
#else
	setData(h->aggr_data, h->datasize, 0, data);
#endif
    } else {
      (h->aggregate)(h->aggr_data,0,data);
    }
  }

  if(( h->cnt == h->nsample) && (!(h->filled)))
    { 
      fprintf(stderr,"CONVERTING HISTOGRAM TO BUCKET MODE\n");
      fill(h);
    }
  if(! (h->filled) ) {
    h->thescore[h->cnt] = thescore;

    if(h->datasize) {
#ifdef USE_CUSTOM_INDEXING
      (h->setdata)(h->sample_data,h->cnt,data);
#else
	setData(h->sample_data, h->datasize, h->cnt, data);
#endif
    }
  } else {
    int ib;
    ib = bucket_from_score(h,thescore);
    h->bucket_cnt[ib] += 1;
    if( h->bucket_cnt[ib] == 1) {
      h->bucket_min[ib] = thescore;
      h->bucket_max[ib] = thescore;
    } else {
      h->bucket_min[ib] = min(h->bucket_min[ib],thescore);
      h->bucket_max[ib] = max(h->bucket_max[ib],thescore);
    }
    if(h->datasize) {
      if( h->bucket_cnt[ib] == 1) {
#ifdef USE_CUSTOM_INDEXING
	(h->setdata)(h->bucket_data,ib,data);
#else
	setData(h->bucket_data, h->datasize, ib, data);
#endif
      } else {
	(h->aggregate)(h->bucket_data,ib,data);
      }
    }
  }
  h->cnt ++;
}

static int compute_cm(HISTOGRAM *h, int rez, int numb) {
  /* 
     Compute the number of buckets per printed bucket.
   */
  int cm; 
  if (rez == 0)
    { cm = 1;}
  else
    { 
      int nbuck_printed,bck;
      nbuck_printed = (numb-1)/rez+1;
      bck = nbuck_printed*(h->bucket_width);
#ifdef NEVER
      { 
	double nrm;
	nrm = pow(10.,floor(log10(1.*bck)));
	if (bck/nrm < 1.001)
	  bck = nrm;
	else if (bck/nrm < 2.)
	  bck = 2.*nrm;
	else if (bck/nrm < 5.)
	  bck = 5.*nrm;
	else
        bck = 10.*nrm;
	if (bck % h->bucket_width != 0) bck *= 2;
      }
#endif
      cm = bck/h->bucket_width;
    }
  //  fprintf(stderr,"* compute_cm returning %d rez = %d num = %d\n",
  //  cm, rez, numb);

  return cm;
}

static int compute_precision(HISTOGRAM *h,int min,int max) 
{
  int i, prec;
  /* "prec" is the number of digits necessary to print the integer. */
  i = (h->low) + (max+1)*(h->bucket_width) - 1;
  if (i == 0)
    { prec = 0;}
  else
    { prec = (int) log10(1.*i);}
  i = h->low + min*h->bucket_width;
  if ((i < 0) && (prec < 1. + log10(-1.*i)))
    prec = 1 + (int) log10(-1.*i);
  prec += 1;
  return prec;
}

static int compute_proc(HISTOGRAM *h,int cs,int cm)
{
  int max,ib,proc;
  /* Compute the maximum number of samples in a consolidated bucket.
     Call it "max". */
  max = 0;
  for (ib = cs; ib+cm >= 1; ib -= cm) {
    int j;
    int s;
    s = 0;
    for (j = 0; j < cm; j++)
      if( (ib+j >= 0) && (ib+j < h->nsample))
	{
	  s += h->thescore[ib+j];
	}
    if (s > max) max = s;
  }
  /* max is now the maximum value in a consolidated bucket */
  
  if (max == 0)
    proc = 1;
  else
    proc = (int) log10(1.*max) + 1;
  
  return proc;
}

/* Print a histogram of "h" with approximately "rez" buckets and with
   the  ":" of the output in column "indent".  If "rez" = 0 then histogram
   is printed in the finest bucket resolution possible.  (N.B. A call
   to print histogram has the effect of shifting it to bucket mode
   regardless of how many samples have been added thus far.)  The
   function returns the field width used to print the frequency (in
   case one would like to add some aligned columns after the output). */

void print_histogram(FILE *fout, HISTOGRAM *h, int rez, int indent)
{ 
  int numb; /* The number of bins that have data. */
  int minb=0; /* The lowest index of a bin with data. */
  int maxb=0; /* The highest index of a bin with data. */
  register int cm, cs;
  int prec = 6; /* the field width used to print the bin size  */
  int proc = 6; /* the field width used to print the frequency */
  int single;
  register int ib;

#ifdef DEBUG
  printf("h->cnt = %d\n",h->cnt);
  printf("h->filled = %d\n",h->filled);
#endif

  if (h->cnt == 0) 
    { return ;}

#ifdef DEBUG3
  {
    int is;
    for (is = 0; is < min((h->cnt),(h->nsample)); is++) {
      DataType *data;
      fprintf(fout," %d : %d ",is,h->thescore[is]);
      if(h->datasize) {
#ifdef USE_CUSTOM_INDEXING 
	data = (h->indexdata)(h->sample_data,is);
#else
	assert(is >= 0 && is < h->nbucket);
	data = getData(h->sample_data,h->datasize,is);
#endif
	(h->printdata)(fout,data);
      }
      fprintf(fout,"\n");
    }
  }
#endif

  if (! h->filled)
    { fill(h);}
  
  numb = 0;
  for (ib = 0; ib < h->nbucket; ib++) {
    int icnt;
    icnt   =  h->bucket_cnt[ib];
    if(icnt > 0)
      { if (numb == 0) minb = ib;
        numb += 1;
        maxb = ib;
#ifdef DEBUG3
	fprintf(fout," %d : %d ",ib, icnt);
	if(h->datasize) {
	  DataType *data;
#if 0
	  data = (h->indexdata)(h->bucket_data,ib);
#else
	  assert(ib >= 0 && ib < h->nbucket);
	  data = getData(h->bucket_data,h->datasize,ib);
#endif
	  (h->printdata)(fout,data);
	}
	fprintf(fout,"\n");
#endif
      }
  }
#ifdef DEBUG
  h->bucket_width = (h->hgh-h->low-1)/(h->nbucket) + 1;
  printf("numb,minb,maxb = %d,%d,%d\n",numb,minb,maxb);
  printf("h->sync,h->low,h->bucket_width = %d,%d,%d\n",
	 h->sync,h->low,h->bucket_width);
#endif
  cm = compute_cm(h, rez, numb);
  cs = (h->sync - h->low)/h->bucket_width;
  //  printf("1* cm,cs = %d,%d\n",cm,cs);
  if (cs > 0) {
    cs = cs%cm - cm;
  } else {
    cs = - ( (-cs)%cm );
  }
  //  printf("2* cm,cs = %d,%d\n",cm,cs);
  cs += ((h->nbucket - cs)/cm+1)*cm;
  //  printf("3* cm,cs = %d,%d\n",cm,cs);

#ifdef DEBUG
  printf("cm,cs = %d,%d\n",cm,cs);
#endif

#ifdef NEVER
  prec = compute_precision(h,minb,maxb);
  proc = compute_proc(h,cs,cm);
#endif

  single = (cm*h->bucket_width == 1);
  single = FALSE;
  if (single)
    { if (prec >= indent) indent = prec+1; }
  else
    { if (2*prec + 3 >= indent) indent = 2*prec+4; }

  if(h->datasize) {
    assert(cm == 1);
  }
  for (ib = cs; ib+cm >= 1; ib -= cm) {
    int j;
    int sum_of_cnt = 0;
    int min_score, max_score;

#ifdef NEXT
    (h->cleardata)(h->temp_data);

    min_score = (h->low) + ib*(h->bucket_width);
    max_score = (h->low) + (ib+cm)*(h->bucket_width) - 1;
#endif

    min_score = h->bucket_min[ib];
    max_score = h->bucket_max[ib];
    min_score = max(min_score,h->min);
    max_score = min(max_score,h->max);
    for (j = 0; j < cm; j++) {
      if( (ib+j >= 0) && (ib+j < h->nbucket) )
	{ 
	  sum_of_cnt += h->bucket_cnt[ib+j];
	  min_score = min(min_score,h->bucket_min[ib+j]);
	  max_score = max(max_score,h->bucket_max[ib+j]);
	  if(h->datasize) {
	    DataType *data;
#ifdef USE_CUSTOM_INDEXING
	    data = (h->indexdata)(h->bucket_data,ib+j);
#else
	    assert((ib+j) >= 0 && (ib+j) < h->nbucket);
	    data = getData(h->bucket_data,h->datasize,ib+j);
#endif
	    //	    fprintf(stderr,"* j = %d ib + j = %d data = 0x%x min %d max %d\n", j,ib + j, 
	    //    data, min_score, max_score);
	    if( j == 0 ) {
#ifdef USE_CUSTOM_INDEXING
	      (h->setdata)(h->temp_data,0,data);
#else
	      memcpy(h->temp_data, data, h->datasize);
	      //setData(h->temp_data, h->datasize, 0, data);
#endif
	    } else {
	      (h->aggregate)(h->temp_data,0,data);
	    }
	  }
	}
      if(h->datasize) {
	if( ib == cs ) {
#ifdef USE_CUSTOM_INDEXING
	  (h->setdata)(h->scan_data,0,h->temp_data);
#else
	  memcpy(h->scan_data, h->temp_data, h->datasize);
	  // setData(h->scan_data, h->datasize, 0, h->temp_data);
#endif	
	} else {
	  (h->aggregate)(h->scan_data,0,h->temp_data);
	}
      }
    
    if (sum_of_cnt > 0)
      {
	//	fprintf(stderr,"* single %d min %d max %d sum_of_cnt = %d\n",
	//single, min_score, max_score, sum_of_cnt);
	if(single || (min_score == max_score)) {
	  fprintf(fout, "%*s%*d: %*d ",indent-prec,
		  "",prec,min_score,proc,sum_of_cnt);
	  if(h->datasize) {
	    (h->printdata)(fout,h->temp_data,h->scan_data,h->aggr_data);
	  }
	  fprintf(fout,"\n");
	} else {
	  fprintf(fout,"%*s%*d - %*d: %*d ",
		  indent-2*prec-3,"",
		  prec,min_score,prec,max_score,proc,sum_of_cnt);
	  if(h->datasize) {
	    (h->printdata)(fout,h->temp_data,h->scan_data,h->aggr_data);
	  }
	  fprintf(fout,"\n");
	}
      }
  }
  }
}
#if 0
#ifdef NEVER
 h->cnt += 1;
 h->sum += data;
 h->mom += data*data;
 if( (h->cnt == 1) || (h->min > data))
   { h->min = data;}
 if( (h->cnt == 1) || (h->max < data))
   { h->max = data;}

/* Return the average of the data points in histogram "h" */

double histogram_avg(HISTOGRAM *h)
{ if (h->cnt == 0)
    return (0.);
  else
    return ((1.*h->sum)/h->cnt);
}

/* Return the standard deviation of the data points in histogram "h" */

double histogram_stdev(HISTOGRAM *h)
{ double fmom, fsum, fcnt;

  fmom = h->mom;
  fsum = h->sum;
  fcnt = h->cnt;
  if (h->cnt == 0)
    return (0.);
  else
    return (sqrt((fcnt*fmom - fsum*fsum)/(fcnt*fcnt)));
}

#endif



#ifdef NEVER
  {
    /* Pretty up something */
    int bck,idiff;
    double nrm;
    bck = (hgh-low)/(h->nbucket-20) + 1;
    nrm = pow(10.,floor(log10(1.*bck)));
    if (bck/nrm < 1.001)
      bck = nrm;
    else if (bck/nrm < 2.)
      bck = 2.*nrm;
    else if (bck/nrm < 5.)
      bck = 5.*nrm;
    else
      bck = 10.*nrm;
    h->bucket_width = bck;
    low = low - 10*bck;
    
    idiff = h->sync - low;
    if (idiff < 0) {
      h->low = low - (-idiff)%(h->bucket_width);
    } else {
      h->low = low + ( idiff)%(h->bucket_width);
    }
  }

#endif /* NEVER */
#endif

