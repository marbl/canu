
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

/* $Id: AS_PER_fragStorePartition.h,v 1.5 2007-02-09 20:59:56 brianwalenz Exp $ */

#ifndef AS_PER_FRAGSTORE_PARTITION_H
#define AS_PER_FRAGSTORE_PARTITION_H

#include "AS_PER_genericStore.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_fragStore_private.h"
#include "AS_UTL_Hash.h"

typedef struct {
  int32         partition;
  char         *fragStorePath;
  StoreHandle   partitionStore; // a collection of ShortFragRecords corresponding to the records in this partition
  StoreHandle   sequenceStore;
  StoreHandle   sourceStore;
  StoreHandle   globalFrgStore; // a collection of ShortFragRecords corresponding to the records in this partition
  HashTable_AS *index;
}tFragStorePartition;

int32                testOpenFragStorePartition(char *fragStorePath,
                                                int32 partition);

tFragStorePartition *openFragStorePartition(char *fragStorePath,
                                            int32 partition,
                                            int loadData);

void                 closeFragStorePartition(tFragStorePartition *partition);

int                  getFragStorePartition(tFragStorePartition *partition,
                                           int32 indx,
                                           int32 getFlags,
                                           ReadStructp rs);

int                  isMemberFragStorePartition(tFragStorePartition *partition,
                                                int32 indx);



#endif
