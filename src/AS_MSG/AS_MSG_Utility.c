
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
 * Celera Whole Genome Shotgun Assembler.
 * AS_MSG - Directory for "Proto-I/O", the file-based data pipeline.
 * Utility - Utility functions for the I/O. 
 * Jason Miller, 7/2001.
 *********************************************************************/

static char CM_ID[]= "$Id: AS_MSG_Utility.c,v 1.4 2005-03-22 19:48:58 jason_miller Exp $";

#include <assert.h>
#include "AS_MSG_Utility.h"

/* ---------------------------------------------------
 * Convert values of type: char
 * to the enumerated type: FragType
 *
 * This violates the abstraction provided by enumerated types.
 * But our storage requirements demand this byte-saving step.
 *
 * Param 'strict' controls behavior on bad input.
 * If strict, function dies when input is not valid.
 * Else, function returns AS_EXTR when input is not valid.
   ----------------------------------------------------*/
FragType AS_MSG_SafeConvert_charToFragType (const char input, bool strict) {
  FragType output;
  // Note on efficiency:
  // The vast majority of data will satisfy the first test.
  if (input=='R') {
    output=AS_READ;
  }
  else {
    if (input=='X') {output=AS_EXTR;}
    else if (input=='T') {output=AS_TRNR;}
    else if (input=='E') {output=AS_EBAC;}
    else if (input=='L') {output=AS_LBAC;}
    else if (input=='U') {output=AS_UBAC;}
    else if (input=='F') {output=AS_FBAC;}
    else if (input=='S') {output=AS_STS;}
    else if (input=='u') {output=AS_UNITIG;}
    else if (input=='c') {output=AS_CONTIG;}
    else if (input=='B') {output=AS_BACTIG;}
    else if (input=='C') {output=AS_FULLBAC;}
    else if (input=='G') {output=AS_B_READ;}
    else {
      if (strict) {
	assert(FALSE);
      }
      output=AS_EXTR; // arbitrary choice .. there is no AS_UNDEFINED
    }
  }
  return output;
}
