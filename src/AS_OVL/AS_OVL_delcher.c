
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
static char fileID[] = "$Id: AS_OVL_delcher.c,v 1.8 2008-06-27 06:29:17 brianwalenz Exp $";


/* RCS info
 * $Id: AS_OVL_delcher.c,v 1.8 2008-06-27 06:29:17 brianwalenz Exp $
 * $Revision: 1.8 $
*/



#include  "AS_OVL_delcher.h"


int  Global_Debug_Flag = FALSE;
  // Flag for debugging
int  Verbose_Level = 0;
  // Determines amount of diagnostic printout



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
        fprintf (stderr, "ERROR %d:  Could not open file  %s \n",
                 errno, Filename);
        perror (strerror (errno));
        exit (1);
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

