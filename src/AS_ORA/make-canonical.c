
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/make-canonical.c,v $
$Revision: 1.3 $
**********************************************************************/

/**********************************************************************
Module: make-canonical

Description: Reads output from  get-olaps  and converts it to
             "canonical format" used by  overlap_regressor
             so the two can be compared.

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

#define  MAX_LINE         1000

/*********************************************************************/


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
   char  line [MAX_LINE];
   int  a_id, b_id, a_hang, b_hang;
   int  lo_id, hi_id, output_hang;
   char  output_comp;
   char  orient [MAX_LINE];
   double  error_rate;

   while  (fgets (line, MAX_LINE, stdin) != NULL)
     {
      sscanf (line, "%d %d %d %d %s %lf",
              & a_id, & b_id, & a_hang, & b_hang, orient, & error_rate);

      assert (a_id != b_id);
      assert (a_hang >= 0);

      if  (a_id < b_id)
          {
           lo_id = a_id;
           hi_id = b_id;

           switch  (toupper (orient [0]))
             {
              case  'N' :
                output_hang = a_hang;
                output_comp = FALSE;
                break;
              case  'I' :
                output_hang = a_hang;
                output_comp = TRUE;
                break;
              case  'O' :
                output_hang = - b_hang;
                output_comp = TRUE;
                break;
              default :
                fprintf (stderr, "ERROR:  Bad orientation.  Line is:\n%s",
                         line);
                exit (EXIT_FAILURE);
             }
          }
        else
          {
           lo_id = b_id;
           hi_id = a_id;

           switch  (toupper (orient [0]))
             {
              case  'N' :
                output_hang = - a_hang;
                output_comp = FALSE;
                break;
              case  'I' :
                output_hang = b_hang;
                output_comp = TRUE;
                break;
              case  'O' :
                output_hang = - a_hang;
                output_comp = TRUE;
                break;
              default :
                fprintf (stderr, "ERROR:  Bad orientation.  Line is:\n%s",
                         line);
                exit (EXIT_FAILURE);
             }
          }

      printf ("%8u %8u %c %4d %5.2f\n",
              lo_id, hi_id, output_comp ? 'r' : 'f', output_hang,
              error_rate);
     }

   fprintf (stderr, "Done.\n");

   return 0;
  }



