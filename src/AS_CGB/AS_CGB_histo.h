
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
 * $Id: AS_CGB_histo.h,v 1.2 2004-09-23 20:25:01 mcschatz Exp $
 *
 * Module: AS_CGB_histo.h
 * Description: The header file for the hitogramming routines.
 * Assumptions:
 *********************************************************************/

#ifndef AS_CGB_HISTO_INCLUDE
#define AS_CGB_HISTO_INCLUDE

typedef void HistoDataType;
typedef void Histogram_t;

/* Create an initially empty histogram where buckets will be
   selected so that "sync" will be the low end of some bucket.  */

extern Histogram_t *create_histogram(
			    int nsample,
			    int nbucket,
			    //int sync0,
			    int filled,
			    int logarithmic);


extern void extend_histogram(  
		      Histogram_t *h,
		      size_t nbytes,
		      HistoDataType * (*indexdata)(HistoDataType *b,int ib),
		      void (*setdata)(HistoDataType *a,int ib,HistoDataType *b),
		      void (*aggregate)(HistoDataType *a,int ib,HistoDataType *b),
		      void (*printdata)(FILE *fout,
					HistoDataType *,
					HistoDataType *,
					HistoDataType *)
	);

/* Free the data structure for histogram "h" */

extern void free_histogram(Histogram_t *);

/* Add data point "data" to histogram "h" */

extern void add_to_histogram(Histogram_t *, int thescore, HistoDataType *);

/* Print a histogram of "h" with approximately "rez" buckets and with
   the  ":" of the output in column "indent".  If "rez" = 0 then histogram
   is printed in the finest bucket resolution possible.  (N.B. A call
   to print histogram has the effect of shifting it to bucket mode
   regardless of how many samples have been added thus far.)  The
   function returns the field width used to print the frequency (in
   case one would like to add some aligned columns after the output). */

extern void print_histogram(FILE *,Histogram_t *, int rez, int indent);

/* Return the average of the data points in histogram "h" */

extern double histogram_avg(Histogram_t *);

/* Return the standard deviation of the data points in histogram "h" */

extern double histogram_stdev(Histogram_t *);

#endif /*AS_CGB_HISTO_INCLUDE*/
