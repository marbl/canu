
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
   Module:  AS_OVL
   Description:  Assembly Overlap Module.  Computes overlaps between
      pairs of DNA strings.
   Assumptions:  Input meets specifications in the ProtoIO documents
 *********************************************************************/

/* RCS info
 * $Id: AS_OVL_olapstats.h,v 1.7 2008-10-08 22:02:58 brianwalenz Exp $
 * $Revision: 1.7 $
*/


#ifndef AS_OVL_OLAPSTATS_H
#define AS_OVL_OLAPSTATS_H

static const char *rcsid_AS_OVL_OLAPSTATS_H = "$Id: AS_OVL_olapstats.h,v 1.7 2008-10-08 22:02:58 brianwalenz Exp $";

//
// Component:
//   AS_OVL_overlap
//
//   Art. Delcher
//   Last revised:  18 Jan 99
//
//  Declarations used to collect counts of various events in
//  the overlap module.

#include <float.h>

#define  EDIT_DIST_MULTIPLE  10

int64  Kmer_Hits_Ct;
int64  Branch_Point_Ct = 0;
int64  Duplicate_Olap_Ct = 0;
int64  Regular_Olap_Ct = 0;
int64  Save_Olap_Ct;
int64  Too_Short_Ct = 0;

int  Is_Duplicate_Olap;


typedef  struct Distribution
  {
   int  N;
   int  * Ct;
   double  * Thold, Max, Sum;
  }  Distrib_t;

//  Holds  N+1  counts.  For  0 <= i <= N ,  Ct [i]  is the
//  number of values that were  <= Thold [i] .


void  Init_Distrib  (Distrib_t * D, int N)

//  Constructor to allocate memory for a new  D  with  D -> N = N .

  {
   int  i;

   D -> Max = - DBL_MAX;
   D -> Sum = 0.0;
   D -> N = N;
   D -> Ct = (int *) safe_calloc (N + 1, sizeof (int));
   D -> Thold = (double *) safe_calloc (N + 1, sizeof (double));

   return;
  }


void  Incr_Distrib  (Distrib_t * D, double Val)

//  Increment the count in  D  that corresponds to  Val .

  {
   int  i;

   if  (Val > D -> Max)
       D -> Max = Val;
   D -> Sum += Val;
   for  (i = 0;  i < D -> N;  i ++)
     if  (Val <= D -> Thold [i])
         {
          D -> Ct [i] ++;
          return;
         }

   D -> Ct [D -> N] ++;

   return;
  }


void  Print_Distrib  (Distrib_t D, char * Heading)

//  Display the values in  D  on the standard error stream.

  {
   int64  Total = 0;
   int  i;

   fprintf (Stat_File, "\n%s\n", Heading);
   fprintf (Stat_File, "     <=     Count\n");
   for  (i = 0;  i < D . N;  i ++)
     {
      fprintf (Stat_File, "%7.0f  %8d\n", D . Thold [i], D . Ct [i]);
      Total += D . Ct [i];
     }
   fprintf (Stat_File, "  Above   %7d\n", D . Ct [D . N]);
   Total += D . Ct [D . N];
   fprintf (Stat_File, "  Total %9ld\n", Total);
   if  (Total > 0)
       {
        fprintf (Stat_File, "  Max   %9.0f\n", D . Max);
        fprintf (Stat_File, "  Sum  %10.0f\n", D . Sum);
        fprintf (Stat_File, "  Avg   %11.1f\n", D . Sum / Total);
       }

   return;
  }


#endif

