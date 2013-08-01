
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

#ifndef OVERLAPSTORE_H
#define OVERLAPSTORE_H

static const char *rcsid_OVERLAPSTORE_H = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "AS_global.H"

typedef  enum {
  PLC_NONE,
  PLC_ALL,
  PLC_INTERNAL
} Ovl_Skip_Type_t;


void
updateErates(char *storeName, char *eratesName);

void
dumpStore(char *storeName, uint32 dumpBinary, double dumpERate, uint32 dumpType, uint32 bgnIID, uint32 endIID, uint32 qryIID);

void
dumpPicture(char *ovlName, char *gkpName, uint32 clearRegion, double dumpERate, uint32 dumpType, uint32 qryIID);

void
estimateGenomeLength(char    *ovlName,
                     char    *gkpName,
                     uint32   overlapLimit,  //  K
                     uint32   bgnIID,
                     uint32   endIID,
                     uint32   into,
                     uint32   windowSize,
                     uint32   minOvl);

int
OVSoverlap_sort(const void *a, const void *b);


#define OP_NONE           0
#define OP_DUMP           3
#define OP_DUMP_PICTURE   4
#define OP_GENOME_LENGTH  5
#define OP_UPDATE_ERATES  6

#define DUMP_5p         1
#define DUMP_3p         2
#define DUMP_CONTAINED  4
#define DUMP_CONTAINS   8

#endif  //  OVERLAPSTORE_H
