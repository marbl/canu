
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
//  delcher.c
//
//  Definitions of functions declared in  delcher.h


#include  "delcher.h"



int  All_White_Space
    (const char * S)

/* Return  TRUE  iff string  S  contains no characters other than
*  whitespace characters as defined by the builtin function
*  isspace () . */

  {
   while  (* S != '\0' && isspace (* S))
     S ++;

   return  (* S == '\0');
  }



FILE *  File_Open
    (const char * Filename, const char * Mode)

/* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FILE  *  fp;

   fp = fopen (Filename, Mode);
   if  (fp == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
        exit (EXIT_FAILURE);
       }

   return  fp;
  }



void *  Safe_calloc
    (size_t N, size_t Len)

/* Allocate and return a pointer to an array of  N  elements of
*   Len  bytes each.  All set to 0.  Exit if fail. */

  {
   void  * P;

   P = calloc (N, Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  calloc failed\n");
        exit (EXIT_FAILURE);
       }

   return  P;
  }



void *  Safe_malloc
    (size_t Len)

/* Allocate and return a pointer to  Len  bytes of memory.
*  Exit if fail. */

  {
   void  * P;

   P = malloc (Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  malloc failed\n");
        exit (EXIT_FAILURE);
       }

   return  P;
  }



void *  Safe_realloc
    (void * Q, size_t Len)

/* Reallocate memory for  Q  to  Len  bytes and return a
*  pointer to the new memory.  Exit if fail. */

  {
   void  * P;

   P = realloc (Q, Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  realloc failed\n");
        exit (EXIT_FAILURE);
       }

   return  P;
  }



