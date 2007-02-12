
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

/* $Id: AS_PER_ReadStruct.h,v 1.7 2007-02-12 22:16:58 brianwalenz Exp $ */

/*************************************************************************
 Module:  AS_PER_ReadStruct
 Description:
     This module defines the interface and implementation of the 
 opaque datatype used by the Fragment Store.

     On all set operations, data is COPIED BY VALUE, so no references to
 client-owned memory remain in the ReadStruct (e.g., setSequence).

     On all get operations, data is COPIED BY VALUE into client-owned
 memory.  Client never acquires references to internal ReadStruct data.

 Assumptions:
      
 Document:
      FragStore.rtf

 *************************************************************************/

#ifndef AS_PER_READSTRUCT_H
#define AS_PER_READSTRUCT_H

#include "AS_global.h"
#include "AS_PER_gkpStore.h"


// Clients must specify which clear range to use.
// Prior to Oct 2001, the original clear range was the only one stored.
// Starting Oct 2001, Scaffolder modifies the CGW clear range on some frags.
// The default 'get' function returns the latest clear range.
// Latest is defined as the highest number here.
//
#define READSTRUCT_LATEST   0
#define READSTRUCT_ORIGINAL 1
#define READSTRUCT_OVL      2
#define READSTRUCT_CNS      3
#define READSTRUCT_CGW      4

ReadStructp new_ReadStruct(void);
void        delete_ReadStruct(ReadStructp r);
void        clear_ReadStruct(ReadStructp r);


int setClearRegion_ReadStruct(ReadStructp rs, 
                              uint32 start,
                              uint32 end,
                              uint32 flags);


int getAccID_ReadStruct(ReadStructp rs, uint64 *accID);
int getReadIndex_ReadStruct(ReadStructp rs, uint32 *readIndex);
int getLocalIndex_ReadStruct(ReadStructp rs, uint32 *localIndex);
int getIsDeleted_ReadStruct(ReadStructp rs, uint32 *isDeleted);

int getClearRegion_ReadStruct(ReadStructp rs, uint32 *start, uint32 *end, uint32 flags);

int getSequence_ReadStruct(ReadStructp rs, char *sequence, char *quality, int length);
int getSource_ReadStruct(ReadStructp rs, char *source, int length);

int dump_ReadStruct(ReadStructp rs, FILE *fp, int clearRangeOnly);

#endif
