
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
      Common include files and a few simple functions.
   Assumptions:  Input meets specifications in the ProtoIO documents
 *********************************************************************/


/* RCS info
 * $Id: AS_OVL_delcher.h,v 1.2 2004-09-23 20:25:25 mcschatz Exp $
 * $Revision: 1.2 $
*/


#ifndef  __DELCHER_H_INCLUDED
#define  __DELCHER_H_INCLUDED


#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>
#include  <errno.h>
#include  <unistd.h>

#include  "AS_global.h"

#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif

#define  FILE_OPEN_FAILURE   134

int  All_White_Space
    (const char *);
FILE *  File_Open
    (const char *, const char *);
int  File_Exists
    (const char *);
void *  Safe_calloc
    (size_t, size_t);
size_t  Safe_fread
    (void * ptr, size_t size, size_t ct, FILE * fp);
size_t  Safe_fwrite
    (const void * ptr, size_t size, size_t ct, FILE * fp);
void *  Safe_malloc
    (size_t);
void *  Safe_realloc
    (void *, size_t);
int  Safe_remove
    (char * filename);
int  Safe_rename
    (char * oldname, char * newname);


#endif

