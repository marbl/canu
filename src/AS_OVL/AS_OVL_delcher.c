
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
      Definitions of functions declared in  delcher.h
   Assumptions:  Input meets specifications in the ProtoIO documents
 *********************************************************************/
static char fileID[] = "$Id: AS_OVL_delcher.c,v 1.3 2005-03-22 19:06:40 jason_miller Exp $";


/* RCS info
 * $Id: AS_OVL_delcher.c,v 1.3 2005-03-22 19:06:40 jason_miller Exp $
 * $Revision: 1.3 $
*/



#include  "AS_OVL_delcher.h"



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
   int  retry;

   fp = fopen (Filename, Mode);
   for  (retry = 0;  fp == NULL && retry < 3;  retry ++)
     {
      sleep (10);
      fp = fopen (Filename, Mode);
     }
   if  (fp == NULL)
       {
#if 0
        char  buff [1000];
#endif

        fprintf (stderr, "ERROR %d:  Could not open file  %s \n",
                 errno, Filename);
        perror (strerror (errno));
#if  0
        sprintf (buff, "mail Randall.Bolanos \n"
                 "Overlap error %d  file %s\n"
                 "%s\n%c\n",
                 errno, Filename, strerror (errno), '\04');
        system (buff);
#endif
        exit (FILE_OPEN_FAILURE);
       }

   return  fp;
  }



int File_Exists (const char * Filename){

  /* Test for filename existence */

   FILE  *  fp;

   fp = fopen (Filename, "r+");
   if  (fp)
       {
	 fclose(fp);
	 return 1;
       }

   return  0;

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



size_t  Safe_fread
    (void * ptr, size_t size, size_t ct, FILE * fp)

//  Call  fread (ptr, size, ct, fp)  and check for errors

  {
   size_t  result;

   result = fread (ptr, size, ct, fp);
   if  (result != ct)
       {
        fprintf (stderr, "ERROR:  fread failed\n");
        exit (EXIT_FAILURE);
       }

   return  result;
  }



size_t  Safe_fwrite
    (const void * ptr, size_t size, size_t ct, FILE * fp)

//  Call  fwrite (ptr, size, ct, fp)  and check for errors

  {
   size_t  result;

   result = fwrite (ptr, size, ct, fp);
   if  (result != ct)
       {
        fprintf (stderr, "ERROR:  fwrite failed\n");
        exit (EXIT_FAILURE);
       }

   return  result;
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

   if  (Q == NULL)
       P = malloc (Len);
     else
       P = realloc (Q, Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  realloc failed\n");
        exit (EXIT_FAILURE);
       }

   return  P;
  }



int  Safe_remove
    (char * filename)

//  Delete file named  filename .  Print error message and die if
//  fail.

  {
   int  result;

   result = remove (filename);

   if  (result != 0)
       {
        fprintf (stderr, "ERROR:  Could not delete file \"%s\"\n",
                 filename);
        exit (EXIT_FAILURE);
       }

   return  result;
  }



int  Safe_rename
    (char * oldname, char * newname)

//  Rename file  oldname  to  newname .  Print error message and die if
//  fail.

  {
   int  result;

   result = rename (oldname, newname);

   if  (result != 0)
       {
        fprintf (stderr, "ERROR:  Could not rename file \"%s\" to \"%s\"\n",
                 oldname, newname);
        exit (EXIT_FAILURE);
       }

   return  result;
  }

