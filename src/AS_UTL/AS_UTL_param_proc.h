
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
/* 	$Id: AS_UTL_param_proc.h,v 1.2 2004-09-23 20:25:29 mcschatz Exp $	 */
#ifndef AS_UTL_PARAM_PROC
#define AS_UTL_PARAM_PROC
#endif
/*************************************************************************
 Module:  AS_UTL_param_proc
 Description:
     This module allows for the reading of a set of parameters from a 
 file.  The parameters should be in the form "module.param value", one to
 a line.

 Assumptions:
      None.
 Document:
      TBD

 *************************************************************************/

#define INPUTMAX 256
#define MAXRETURNBUFFERLENGTH 2 * 1024
#define PARAM_PROC_SUCCESS 1
#define PARAM_PROC_FAILURE 0

typedef struct paramEntry
{
	  char paramModule[INPUTMAX];
	  char paramName[INPUTMAX];
	  char paramValue[INPUTMAX];
	  struct paramEntry *next;
} paramEntryT;

int loadParams(const char * const filename);
// Returns PARAM_PROC_{SUCCESS, FAILURE}.


char* getParam(const char * const paramNameIn);
// The user is responsible for releasing the memory in the returned
// string using free().  The function returns NULL apon failure.

int getAllParams(const char * const moduleNameIn, char ** returnBuffer);
// The user is responsible for releasing the memory in the returned
// buffer using free().  

