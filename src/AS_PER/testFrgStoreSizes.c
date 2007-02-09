
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

#include <stdlib.h> 
#include <stdio.h> 
#include <assert.h>

#include <unistd.h>
#include "AS_UTL_PHash.h"
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_fragStore_private.h" 
#include "AS_PER_gkpStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"
#include "AS_UTL_Var.h"
#include "AS_SDB_SequenceDB.h"

/* the following copied from AS_UTL_Var.c so that somethign not of
   general use doesn't have to be put in AS_UTL_Var.h */
#define PORTABLE_FILE_IO
#ifndef PORTABLE_FILE_IO
typedef VarArrayType FileVarArrayType;
#else // PORTABLE_FILE_IO
typedef struct {
  uint64 Elements; /* The Data pointer. Must be cast to the appropriate type */
  uint64 sizeofElement;   /* The size in bytes of the appropriate type. */
  uint64 numElements;     
  uint64 allocatedElements;
  char typeofElement[VA_TYPENAMELEN]; /* The name of the data type of 
					  each element. */
} FileVarArrayType;
#endif // PORTABLE_FILE_IO



int main(int argc, char *argv[]){
  printf("ShortFragRecord size is %ld\n",sizeof(ShortFragRecord));
  printf("FragStore size is %ld\n",sizeof(FragStore));
  printf("VLSTRING_SIZE_T size is %ld\n",sizeof(VLSTRING_SIZE_T));
  printf("GateKeeperBatchRecord size is %ld\n",sizeof(GateKeeperBatchRecord));
  printf("GateKeeperFragmentRecord size is %ld\n",sizeof(GateKeeperFragmentRecord));
  printf("GateKeeperDistanceRecord size is %ld\n",sizeof(GateKeeperDistanceRecord));
  printf("GateKeeperLinkRecord size is %ld\n",sizeof(GateKeeperLinkRecord));
  printf("PHashTable_AS size is %ld\n",sizeof(PHashTable_AS));
  printf("tSequenceDB size is %ld\n",sizeof(tSequenceDB));
  printf("tMARrecord size is %ld\n",sizeof(tMARecord));
  printf("IntMultiPos size is %ld\n",sizeof(IntMultiPos));
  printf("VarArrayType size is %ld\n",sizeof(VarArrayType));
  printf("FileVarArrayType size is %ld\n",sizeof(FileVarArrayType));
  printf("StoreStat size is %ld\n",sizeof(StoreStat));
  printf("CDS_IID_t size is %ld\n",sizeof(CDS_IID_t));
  printf("CDS_UID_t size is %ld\n",sizeof(CDS_UID_t));
  printf("time_t size is %ld\n",sizeof(time_t));
  printf("size_t size is %ld\n",sizeof(size_t));
  printf("off_t size is %ld\n",sizeof(off_t));
  return (0);
}

