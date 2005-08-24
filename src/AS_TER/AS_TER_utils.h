
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
   CVS_ID:  $Id: AS_TER_utils.h,v 1.5 2005-08-24 07:47:15 brianwalenz Exp $
 *********************************************************************/
#ifndef AS_TER_UTILS_H
#define AS_TER_UTILS_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "AS_global.h"


/* this typedef gives a list of possible error codes
   that index into a table used by the error Handler 
   Edit this table and accordingly update the number
   of error messages as well as the actual error messages
   in the array ErrorString in the file AS_TER_utils.c   */


/*---------------------------------------------------------*/
/* #defines */
/*---------------------------------------------------------*/

#define AS_TER_NumOfErrors 4

/*---------------------------------------------------------*/
/* typedefs */
/*---------------------------------------------------------*/

typedef enum {
  AS_TER_FILE_ERROR,
  AS_TER_MEMORY_ERROR,
  AS_TER_PRECONDITION_ERROR,
  AS_TER_UIDSERVER_ERROR
} ErrorCode; 

typedef enum {
  AS_TER_EXIT_SUCCESS = 0,
  AS_TER_EXIT_FAILURE = 1
} ExitStatus;


typedef enum {
  AS_TER_FILE_EXISTS     = 0,
  AS_TER_FILE_EXISTS_NOT = 1
} FileStatus;



/*---------------------------------------------------------*/
/* functions */
/*---------------------------------------------------------*/


void error(ErrorCode err, const char* errorMesg,  ExitStatus ex, 
	   const char* file, int line);
/************************************************************************/ 
/* The function checks whether the error code err is valid. If not it   */
/* issues an error message. Otherwise it issues the error message given */
/* by errorMesg and exits with status ex.                               */
/************************************************************************/


FileStatus file_exists (const char* filename);
/*****************************************************************/
/* Test for filename existence. Returns AS_TER_FILE_EXISTS if    */
/* Filename exists and AS_TER_FILE_EXISTS_NOT o.w.               */
/*****************************************************************/


FILE* file_open(const char* filename, const char* mode);
/*****************************************************************/
/* Open  Filename  in  Mode  and return a pointer to its control */
/* block.  If it fails it issues an error message and exits.     */
/*****************************************************************/




/**********************************************************************/
/* Initial size of the variable arrays */
/**********************************************************************/


#define ARRAYSIZE  2048



#endif








