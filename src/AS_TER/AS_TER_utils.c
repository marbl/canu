
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
 Module: AS_TER
 Description: Contains small utility functions not directly related
with the Terminator. Could be used by anyone.
 Assumptions: none
**********************************************************************/

static char CM_ID[] = "$Id: AS_TER_utils.c,v 1.2 2004-09-23 20:25:29 mcschatz Exp $";

#include "AS_TER_utils.h"

/*--------------------------------------------------------------------*/
/* Error Handling */
/*--------------------------------------------------------------------*/

static char ErrorString[AS_TER_NumOfErrors][30] = 
{
  "FILE ERROR",
  "MEMORY ERROR",
  "PRECONDITION ERROR",
  "UID SERVER ERROR"
};


void error(ErrorCode err, const char* errorMesg,  ExitStatus ex, 
	   const char* file, int line)
{
  /************************************************************************/ 
  /* The function checks whether the error code err is valid. If not it   */
  /* issues an error message. Otherwise it issues the error message given */
  /* by errorMesg and exits with status ex.                               */
  /************************************************************************/

  int intErr = (int)err;
  if( intErr < 0 || intErr >= AS_TER_NumOfErrors )
    {
    fprintf(stderr,"ERROR in error function : Invalid ErrorCode %d ",intErr);
    exit(AS_TER_EXIT_FAILURE);
  }
  fprintf(stderr,"%s in file %s, line %d : %s\n",
	  ErrorString[intErr],file,line,errorMesg);
  
  exit(ex);
}



/*--------------------------------------------------------------------*/
/*  File Handling routines */
/*--------------------------------------------------------------------*/


FILE*  file_open(const char* fileName, const char* mode)
{
  /*****************************************************************/
  /* Open  Filename  in  Mode  and return a pointer to its control */
  /* block.  If fail, print a message and exit.                    */
  /*****************************************************************/
   FILE* fp;

   fp = fopen (fileName,mode);
   if(fp == NULL)
     {
        char dummy[40];
	sprintf(dummy,"Could not open file  %s \n",fileName);
        error(AS_TER_FILE_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
   }
   return  fp;
}


FileStatus file_exists (const char * fileName)
{
  /*****************************************************************/
  /* Test for filename existence.                                  */
  /*****************************************************************/
  FILE* fp;

  fp = fopen (fileName, "r+");
  if(fp)
    {
      fclose(fp);
      return AS_TER_FILE_EXISTS;
    }
  
  return AS_TER_FILE_EXISTS_NOT;
}


/*--------------------------------------------------------------------*/
/*  Memory Managment routines */
/*--------------------------------------------------------------------*/


void *safe_calloc(size_t num, size_t len)
{
  /*****************************************************************/
  /* Allocate and return a pointer to an array of  num  elements of*/
  /* len  bytes each.  All are set to 0.  Exit if fai.             */
  /*****************************************************************/
  void  *p;

   p = calloc (num, len);
   if  (p == NULL)
     {
       char dummy[40];
       sprintf(dummy,"Could not calloc memory (%d number * %d bytes) \n",(int) num,(int) len);
       error(AS_TER_MEMORY_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
     }
   return  p;
}



void *safe_malloc(size_t len)
{
  /*****************************************************************/
  /* Allocate and return a pointer to len bytes of memory.         */
  /* Len  bytes each.  Exit if fail.                               */
  /*****************************************************************/
  void  *p;

  p = malloc (len);
  if(p == NULL)
    {
      char dummy[40];
      sprintf(dummy,"Could not malloc memory (%d bytes) \n",(int) len);
      error(AS_TER_MEMORY_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }
  
  return  p;
}



void *safe_realloc(void *q, size_t len)
{
  /*****************************************************************/
  /* Reallocate memory for q to len  bytes and return a pointer    */
  /* to the new memory.  Exit if fail.                             */     
  /*****************************************************************/
  void  *p;

  p = realloc (q, len);
  if(p == NULL)
    {
      char dummy[40];
      sprintf(dummy,"Could not malloc memory (%d bytes) \n",(int) len);
      error(AS_TER_MEMORY_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }
  
  return  p;
}


