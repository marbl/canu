
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
/* $Id: agrep.h,v 1.2 2004-09-23 20:25:30 mcschatz Exp $ */

/* REGULAR EXPRESSION DATA TYPE:
     An e-nfa and work storage to do approximate matching.  Occupies
     at most 60N bytes for a regular expression of length N.
*/

typedef void *regexp;

/* RE_PARSE:
     Normally returns regexp object used to search for instances of a
     regular expression.  The parameter i_opt is a boolean input
     parameter and should be set if the search is to be case-insensitive.
     Returns NULL and sets most recent error if the regular expression
     is not syntactically valid or if there is insufficient memory.
*/

regexp re_parse(int i_opt, char *pat, char **mesg, int *pos);

/* RE_MATCH:
     Return 0 or 1 according to whether string is recognized by machine
     regexp with at most thresh differences.  For algorithm details
     see Myers & Miller (1989) Bulletin of Mathematical Biology.
*/

int re_match(regexp rexpr, char *string, int thresh);

/*
 RE_FREE:
     Free regexp data type.
*/

void re_free(regexp rexpr);
