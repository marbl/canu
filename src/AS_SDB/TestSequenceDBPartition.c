
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

#include "AS_global.h" 
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"
#include "AS_SDB_SequenceDB.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignStore_CNS.h"
#include "AS_SDB_SequenceDBPartition.h"



int main(int argc, char **argv){
  char *storeName;
  int storeVersion;
  int partitionNum;
  tSequenceDBPartition *SDBPartition;
  MultiAlignT *ma;
  int i;
  VA_TYPE(int32) *members;

  storeName = argv[1];
  storeVersion = atoi(argv[2]);
  partitionNum = atoi(argv[3]);

  fprintf(stderr,"* Opening store %s version %d partition %d\n",
	  storeName, storeVersion, partitionNum);

  SDBPartition = openSequenceDBPartition(storeName, storeVersion, partitionNum);
  members = GetContentSequenceDBPartition(SDBPartition);

  for(i = 0; i < GetNumint32s(members); i++){
    CDS_CID_t id = *Getint32(members,i);
    ma = loadFromSequenceDBPartition(SDBPartition,id);
    
    fprintf(stderr,"* id " F_CID " has " F_SIZE_T " fragments\n",
	    id, GetNumIntMultiPoss(ma->f_list));
  }

  return 0;
}
