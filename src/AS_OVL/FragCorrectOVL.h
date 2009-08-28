
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

#ifndef FRAGCORRECTOVL_H
#define FRAGCORRECTOVL_H

static const char *rcsid_FRAGCORRECTOVL_H = "$Id: FragCorrectOVL.h,v 1.7 2009-08-28 03:43:44 brianwalenz Exp $";

typedef  enum {
   DELETE,
   A_SUBST,
   C_SUBST,
   G_SUBST,
   T_SUBST,
   A_INSERT,
   C_INSERT,
   G_INSERT,
   T_INSERT,
   NO_VOTE,
   EXTENSION
}  Vote_Value_t;

typedef  struct {
  uint64  is_ID       : 1;
  uint64  keep_left   : 1;    // set true if left overlap degree is low
  uint64  keep_right  : 1;    // set true if right overlap degree is low
  uint64  iid         : 32;
} Frag_ID_t;

typedef  struct {
  uint64  is_ID : 1;
  uint64  pos   : AS_READ_MAX_LONG_LEN_BITS;    // position in fragment
  uint64  type  : 11;
} Correction_t;

typedef  union {
  Frag_ID_t     frag;
  Correction_t  corr;
} Correction_Output_t;


#endif
