
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
/*************************************************
* Module:  AS_OVL_overlap.c
* Description:
*   Find matches between pairs of DNA strings and output overlaps
*   (matches that go to an end of both strings involved) and branch points
*   (matches that do not go to the end of either string).
*   Use hashing to select candidates and simple edit distance
*   to confirm matches.
*
*    Programmer:  A. Delcher
*       Written:  29 May 98
*  Last Revised:  30 Nov 99
*
*
* Assumptions:
*  argv[1] is <filename>[.<ext>]
*  Input is from file  argv [1]
*  Output is to file  <filename> . OUTPUT_FILENAME_EXTENSION
*
* Notes:
*   - Removed  ceil  function everywhere.  Overlaps are all <= ERROR_RATE
*     except where indels change length of olap region.
*   - Use hash table for index of strings that match
*   - Do memory allocation by self.
*   - "Band" edit distance calculations
*   - Avoid "hopeless" edit distance calculations
*   - Find branch points in overlaps.
*
*************************************************/

/* RCS info
 * $Id: AS_OVL_overlap.c,v 1.7 2008-06-27 06:29:18 brianwalenz Exp $
 * $Revision: 1.7 $
*/

static const char CM_ID[] = "$Id: AS_OVL_overlap.c,v 1.7 2008-06-27 06:29:18 brianwalenz Exp $";

/****************************************
 * Version of the main program
 ****************************************/

/* this is the frag version binary */
#undef CONTIG_OVERLAPPER_VERSION

/* this is for the contig version binary */
//#ifndef CONTIG_OVERLAPPER_VERSION
//#define CONTIG_OVERLAPPER_VERSION
//#endif

#include "AS_OVL_overlap_common.h"

