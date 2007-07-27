
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

#ifndef AS_SDB_SEQUENCEDB_PARTITION_H
#define AS_SDB_SEQUENCEDB_PARTITION_H

#include "AS_UTL_Hash.h"
#include "AS_SDB_SequenceDB.h"

typedef struct {
  int elemID;
  int partitionID;
}tPartitionElement;

VA_DEF(tPartitionElement)

typedef struct {
  char *sequenceDBPath;      
  VA_TYPE(tMARecord) *multiAligns; // the fileID in the tMARecord is used for the multi-alignID
  FILE *datafp;                    // single FILE* for the data, all is loaded into memory
  HashTable_AS *index;
} tSequenceDBPartition;

tSequenceDBPartition *openSequenceDBPartition(char *path, int32 revision, int32 partition);

// Returns a reference to a multiAlignT from this partition
//
MultiAlignT *loadFromSequenceDBPartition(tSequenceDBPartition *partition, int32 indx);


#endif
