/**************************************************************************
 * This file is part of A2Amapper.
 * Copyright (c) 2004 Applera Corporation
 * Author: Clark Mobarry
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
**************************************************************************/


#ifndef _ALIGN__HALIGN_H_
#define _ALIGN__HALIGN_H_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

typedef 
struct H_Alignment_t {
        int    offset1, offset2;
        int    len1, len2;
        int    score;
        int    first;
  // CMM 2004 Feb 05: The variable "first" does not really belong in
  // the alignment object.  We should put it into a separate iterator
  // over alignment struct/class.  Currently, this is called by Python
  // so I would need to learn how to make C++ glue code for Python
  // objects.  Currently I just make glue code for Python modules.
        int  * scriptAsArray;
} H_Alignment_t;

void halignStart
( char const * const string1,
  char const * const string2,
  int const offset1,
  int const offset2,
  int const len1,
  int const len2,
  H_Alignment_t ** alignment
  );

void Free_align(H_Alignment_t* aln_ptr);
/* To free the alignment returned by HALIGN. */
  
void printUngappedAlign(H_Alignment_t* aln);
void printUngappedAlignSharpEnds(H_Alignment_t* aln);
void printUngappedAlignSharpEndsOnConsole
( H_Alignment_t* aln, char const * const seq1, char const * const seq2, int const opt);

int iterateUngappedAlignSharpEnds
( H_Alignment_t* aln,
  int& bgn1, int& bgn2, int& len1, int& len2, int& nmat);

#endif

