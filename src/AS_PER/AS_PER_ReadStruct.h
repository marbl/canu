
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

/* RCS Info
 * $Date: 2007-01-28 21:52:25 $
 * $Id: AS_PER_ReadStruct.h,v 1.5 2007-01-28 21:52:25 brianwalenz Exp $
 * $Revision: 1.5 $
 *
 */
#ifndef AS_PER_READSTRUCT_H
#define AS_PER_READSTRUCT_H

#include "AS_global.h"

/* ***************************************************************** */

typedef void *ReadStructp;   /* The handle returned by open/create operations */

/*
 Clients must specify which clear range to use.
 Prior to Oct 2001, the original clear range was the only one stored.
 Starting Oct 2001, Scaffolder modifies the CGW clear range on some frags.
 The default 'get' function returns the latest clear range.
 Latest is defined as the highest number here.
 Client code should specify the NAME because the NUMBERING may have to change.
*/
#define READSTRUCT_LATEST   0
#define READSTRUCT_ORIGINAL 1
#define READSTRUCT_OVL      2
#define READSTRUCT_CNS      3
#define READSTRUCT_CGW      4



/*****************************************************************************
 * Function: newReadStruct
 * Description:
 *     Allocate a readStruct and return an opaque pointer
 *
 * Return Value:
 *     NULL if failure
 ****************************************************************************/
ReadStructp new_ReadStruct(void);


/*****************************************************************************
 * Function: delete_ReadStruct
 * Description:
 *     Delete a ReadStruct
 * Input:
 *       r   The readStruct to delete
 *
 * Return Value:
 *     NULL if failure
 ****************************************************************************/
void        delete_ReadStruct(ReadStructp r);


/*****************************************************************************
 * Function: clear_ReadStruct
 * Description:
 *     Reset the fields of a ReadStruct
 * Input:
 *       r   The readStruct to delete
 *
 ****************************************************************************/
void        clear_ReadStruct(ReadStructp r);


/* Mutators */

/* set<field>_ReadStruct sets field with the passed values */
int setAccID_ReadStruct(ReadStructp rs, uint64 accID);
int setReadIndex_ReadStruct(ReadStructp rs, uint32 readID);
int setReadType_ReadStruct(ReadStructp rs, FragType r);
int setSource_ReadStruct(ReadStructp rs, const char *src);
int setEntryTime_ReadStruct(ReadStructp rs, time_t entryTime);
int setLocID_ReadStruct(ReadStructp rs, uint64 locID);
int setLocalePos_ReadStruct(ReadStructp rs, uint32 start, uint32 end);

// Valid flags are: READSTRUCT_ORIGINAL, 
// READSTRUCT_OVL, READSTRUCT_CGW, READSTRUCT_CNS.
// Note that READSTRUCT_LATEST is not valid.
int setClearRegion_ReadStruct(ReadStructp rs, uint32 start, uint32 end, uint32 flags);

// Removed by Jason, Oct 2001
//int setHasModifiedClearRegion_ReadStruct(ReadStructp rs, uint32 hasModifiedClearRegion);

/*****************************************************************************
 * Function: setSequence
 * Description:
 *     Set the sequence and quality data of a sequence.  This hides any encoding of
 * sequence + quality in the store.
 * Input:
 *       sequence           A string of sequence values
 *       quality            A string of quality values (should be bytes?)
 * Return Value:
 *     Zero if OK
 ****************************************************************************/
int setSequence_ReadStruct(ReadStructp rs, char *sequence, char *quality);


/* Accessors */

/* get<field>_ReadStruct gets the field */

int getAccID_ReadStruct(ReadStructp rs, uint64 *accID);
int getReadIndex_ReadStruct(ReadStructp rs, uint32 *readIndex);
int getLocalIndex_ReadStruct(ReadStructp rs, uint32 *localIndex);
int getReadType_ReadStruct(ReadStructp rs, FragType *r);
int getEntryTime_ReadStruct(ReadStructp rs, time_t *entryTime);
int getLocID_ReadStruct(ReadStructp rs, uint64 *locID);
int getLocalePos_ReadStruct(ReadStructp rs, uint32 *start, uint32 *end);
int getIsDeleted_ReadStruct(ReadStructp rs, uint32 *isDeleted);

// Valid flags are: READSTRUCT_LATEST, READSTRUCT_ORIGINAL, 
// READSTRUCT_OVL, READSTRUCT_CGW, READSTRUCT_CNS.
int getClearRegion_ReadStruct(ReadStructp rs, uint32 *start, uint32 *end, uint32 flags);

/*****************************************************************************
 * Function: getSequence_ReadStruct
 * Description:
 *     Get the sequence and quality data of a sequence.  This hides any encoding of
 * sequence + quality in the store.
 * Input
 *       length             The maximum length of the buffer allocated for quality/sequence
 * Output
 *       sequence           A string of sequence values
 *       quality            A string of quality values (should be bytes?)
 * Return Value:
 *     Zero if OK
 *     Length of required buffer, if buffer is too short
 ****************************************************************************/
int getSequence_ReadStruct(ReadStructp rs, char *sequence, char *quality, int length);


/*****************************************************************************
 * Function: getSource_ReadStruct
 * Description:
 *     Get the source data.
 * Input
 *       length             The length of the buffer allocated for source.
 * Output
 *       source            A string of source.
 * Return Value:
 *     Zero if OK
 ****************************************************************************/
int getSource_ReadStruct(ReadStructp rs, char *source, int length);


/*****************************************************************************
 * Function: dump_ReadStruct
 * Description:
 *     Dump readStruct to FILE
 * Input
 *       rs                ReadStruct
 *       fp                File for output
 * Return Value:
 *     Zero if OK
 ****************************************************************************/
int dump_ReadStruct(ReadStructp rs, FILE *fp, int clearRangeOnly);

#endif
