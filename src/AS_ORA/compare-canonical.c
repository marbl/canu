
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
/**********************************************************************
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/compare-canonical.c,v $
$Revision: 1.1.1.1 $
**********************************************************************/

/**********************************************************************
Module: compare-canonical

Description: Reads two canonical-format overlap files and outputs
             discrepancies.  Both files must be sorted as if by
             Unix  "sort -n"  command.

**********************************************************************/


/*********************************************************************/
// headers
// standard headers
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <ctype.h>
#include <unistd.h>
#include  <assert.h>

// project headers
#include  "AS_global.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_fragStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_ORA_fragments.h"
#include  "AS_ORA_overlaps.h"
#include  "AS_ORA_statistics.h"
#include  "AS_ORA_inlines.h"
/*********************************************************************/


/*********************************************************************/
// defines
#ifndef MAX_SEQUENCE_LENGTH
#define MAX_SEQUENCE_LENGTH 2048
#endif

#ifndef MAX_SOURCE_LENGTH
#define MAX_SOURCE_LENGTH 512
#endif

#define  HANG_TOLERANCE   6        // Max allowed difference between
                                   // reported hangs to be considered
                                   // a match
#define  MAX_ERROR_RATE   6.00
#define  MAX_LINE         1000

/*********************************************************************/

typedef  enum {NEED_BOTH, NEED_A, NEED_B, DONE}
    State_t;

/*********************************************************************/
// structures
/*********************************************************************/


/*********************************************************************/
// function prototypeS


/*********************************************************************/


/*********************************************************************/

int  main
    (int argc, char * argv [] )
  {
   FILE  * a_file, * b_file;
   State_t  state = NEED_BOTH;
   char  a_line [MAX_LINE];
   char  b_line [MAX_LINE];
   char  a_orient [10], b_orient [10];
   int  a_lo, a_hi, a_hang, b_lo, b_hi, b_hang;
   int  olap_len_1, olap_len_2;
   double  a_error, b_error;
   int  a_ct = 0, b_ct = 0;
   int  extra_a_ct = 0, extra_b_ct = 0, diff_ct = 0;

   if  (argc < 3)
       {
        fprintf (stderr, "USAGE:  %s <afile> <bfile>\n", argv [0]);
        exit (EXIT_FAILURE);
       }

   a_file = fopen (argv [1], "r");
   if  (a_file == NULL)
       {
        fprintf (stderr, "Could not open file \"%s\"\n", argv [1]);
        exit (EXIT_FAILURE);
       }

   b_file = fopen (argv [2], "r");
   if  (b_file == NULL)
       {
        fprintf (stderr, "Could not open file \"%s\"\n", argv [2]);
        exit (EXIT_FAILURE);
       }


   while  (state != DONE)
     {
      switch  (state)
        {
         case  NEED_BOTH :
           if  (fgets (a_line, MAX_LINE, a_file) == NULL)
               a_lo = INT_MAX;
             else
               {
                sscanf (a_line, "%d %d %s %d %lf",
                        & a_lo, & a_hi, a_orient, & a_hang, & a_error);
                a_ct ++;
               }
           if  (fgets (b_line, MAX_LINE, b_file) == NULL)
               b_lo = INT_MAX;
             else
               {
                sscanf (b_line, "%d %d %s %d %lf %d %d",
                        & b_lo, & b_hi, b_orient, & b_hang, & b_error,
                        & olap_len_1, & olap_len_2);
                b_ct ++;
               }
           break;
         case  NEED_A :
           if  (fgets (a_line, MAX_LINE, a_file) == NULL)
               a_lo = INT_MAX;
             else
               {
                sscanf (a_line, "%d %d %s %d %lf",
                        & a_lo, & a_hi, a_orient, & a_hang, & a_error);
                a_ct ++;
               }
           break;
         case  NEED_B :
           if  (fgets (b_line, MAX_LINE, b_file) == NULL)
               b_lo = INT_MAX;
             else
               {
                sscanf (b_line, "%d %d %s %d %lf %d %d",
                        & b_lo, & b_hi, b_orient, & b_hang, & b_error,
                        & olap_len_1, & olap_len_2);
                b_ct ++;
               }
           break;
         default :
           assert (FALSE);
        }

      if  (a_lo == INT_MAX && b_lo == INT_MAX)
          break;

      if  (a_lo == INT_MAX
             || b_lo < a_lo
             || (a_lo == b_lo && b_hi < a_hi)
             || (a_lo == b_lo && a_hi == b_hi && b_orient [0] < a_orient [0]))
          {
           double  max_error;

           state = NEED_B;
           if  (olap_len_1 < olap_len_2)
               max_error = (int) (b_error * (olap_len_1 + olap_len_2) + 0.5)
                                    / (2.0 * olap_len_1);
             else
               max_error = (int) (b_error * (olap_len_1 + olap_len_2) + 0.5)
                                    / (2.0 * olap_len_2);
           if  (max_error <= MAX_ERROR_RATE)
               {
                printf ("---\n");
                printf ("b> %s", b_line);
                extra_b_ct ++;
               }
          }
      else if  (b_lo == INT_MAX
             || a_lo < b_lo
             || (a_lo == b_lo && a_hi < b_hi)
             || (a_lo == b_lo && a_hi == b_hi && a_orient [0] < b_orient [0]))
          {
           state = NEED_A;
           printf ("---\n");
           printf ("a> %s", a_line);
           extra_a_ct ++;
          }
        else 
          {
           state = NEED_BOTH;
           if  (abs (a_hang - b_hang) > HANG_TOLERANCE)
               {
                printf ("---\n");
                printf ("a> %s", a_line);
                printf ("b> %s", b_line);
                diff_ct ++;
               }
          }
     }

   fclose (a_file);
   fclose (b_file);

   if  (extra_a_ct + extra_b_ct + diff_ct > 0)
       printf ("---\n");

   printf ("\na lines = %d\n", a_ct);
   printf ("b lines = %d\n", b_ct);
   printf ("   unmatched a olaps = %d\n", extra_a_ct);
   printf ("   unmatched b olaps = %d\n", extra_b_ct);
   printf ("mismatched a/b olaps = %d\n", diff_ct);
   printf ("   total differences = %d\n", extra_a_ct + extra_b_ct + diff_ct);

   return 0;
  }



