
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
/* ************************************************
   Module:
   Programmer:  A. Delcher
   Written:  3 November 1998
   File:  ~adelcher/Assembly/Overlaps/gendat.c
   Description:
   Generate sample fragment data with known overlaps to test
   overlap  program.
   Assumptions:
*************************************************/

/* RCS info
 * $Author: jason_miller $
 * $Locker:  $
 * $Date: 2005-03-22 19:06:40 $
 * $Id: AS_OVL_getdat.c,v 1.3 2005-03-22 19:06:40 jason_miller Exp $
 * $Revision: 1.3 $
 * $State: Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.2  2004/06/17 13:25:42  mpop
 * A few bug fixes dealing with memory allocation, specifically pointers into
 * arrays become invalid when the array is realloc'd.  The current code also
 * compiles on Linux, OSX, OSF1 without changes.
 *
 * Revision 1.1  2004/02/10 14:15:07  dewim
 * Initial release0 of the assembler source code.
 * Compiles on aix, tru64/osf, and linux.
 * Runs a006 on aix but not on tru64/osf or linux.
 *
 * Revision 1.2  1998/12/03 17:35:47  cmobarry
 *  Added some CDS style requirements
 *
 */

#include  "delcher.h"


#define  AVG_FRAG_LEN     500
#define  GENOME_LEN       700
#define  INCR              20
#define  UNIF_DBL_RAND  (((double) rand ()) / RAND_MAX)

#define  SUBSTITUTION_ERROR_PROB  0.005
#define  INSERT_ERROR_PROB        0.0025
#define  DELETE_ERROR_PROB        0.0025


char  Pattern [] = "acgttgcaacgttgcaa";


char  Complement
          (char Ch);
void  Output_With_Error
    (char Ch);


int  main  (int argc, char * argv [])

  {
   char  Genome [GENOME_LEN];
   int  i, j, Ct = 0;

   srand (time (NULL));

   for  (i = 0;  i < GENOME_LEN;  i ++)
     Genome [i] = Pattern [rand () % 17];

   for  (i = 0;  i + AVG_FRAG_LEN - 1 < GENOME_LEN;  i += INCR)
     if  (rand () % 2)
         {
fprintf (stderr, "%3d FORWARD\n", ++ Ct);          
          for  (j = 0;  j < AVG_FRAG_LEN;  j ++)
            Output_With_Error (Genome [i + j]);
          putchar ('\n');
         }
       else
         {
fprintf (stderr, "%3d REVERSE\n", ++ Ct);          
          for  (j = 0;  j < AVG_FRAG_LEN;  j ++)
            Output_With_Error (Complement (Genome [i + AVG_FRAG_LEN - 1 - j]));
          putchar ('\n');
         }

   return  0;
  }



char  Complement
    (char Ch)

  {
   switch  (Ch)
     {
      case  'a' :
        return  't';
      case  'c' :
        return  'g';
      case  'g' :
        return  'c';
      case  't' :
        return  'a';
     }

   exit (-1);
  }



void  Output_With_Error
    (char Ch)

  {
   double  X;
   char  Other_Ch;

   X = UNIF_DBL_RAND;

   if  (X < DELETE_ERROR_PROB)
       return;

   if  (X < DELETE_ERROR_PROB + INSERT_ERROR_PROB)
       {
        putchar (Pattern [rand () % 17]);
        putchar (Ch);
        return;
       }

   if  (X < DELETE_ERROR_PROB + INSERT_ERROR_PROB + SUBSTITUTION_ERROR_PROB)
       {
        while ((Other_Ch = Pattern [rand () % 17]) == Ch)
          ;
        putchar (Other_Ch);
        return;
       }

   putchar (Ch);
   return;
  }
