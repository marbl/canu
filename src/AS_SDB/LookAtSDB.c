
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

FILE *myerr=NULL;
int main(int argc, char *argv[]){
  
// args: id1 bgn end id2 bgn end
  VA_TYPE(IntElementPos) *pos_va=CreateVA_IntElementPos(2);
  IntElementPos pos;
  FragStoreHandle fragStore;
  tSequenceDB *sequenceDB = OpenSequenceDB(argv[1], FALSE, atoi(argv[2]));	
  MultiAlignT *ma1 =  LoadMultiAlignTFromSequenceDB(sequenceDB, atoi(argv[4]), FALSE);
  MultiAlignT *ma2 =  LoadMultiAlignTFromSequenceDB(sequenceDB, atoi(argv[7]), FALSE);
  MultiAlignT *newma;
  int bgn1=atoi(argv[5]);
  int end1=atoi(argv[6]);
  int bgn2=atoi(argv[8]);
  int end2=atoi(argv[9]);
  int min_range=(bgn1<end1)?bgn1:end1;
  min_range=(min_range<bgn2)?min_range:bgn2;
  min_range=(min_range<end2)?min_range:end2;
  bgn1-=min_range;
  end1-=min_range;
  bgn2-=min_range;
  end2-=min_range;
  fragStore = openFragStore(argv[3], "rb");
  pos.type = AS_UNITIG;
  pos.ident = atoi(argv[4]);
  pos.position.bgn=bgn1;
  pos.position.end=end1;
  AppendVA_IntElementPos(pos_va,&pos);
  pos.ident = atoi(argv[7]);
  pos.position.bgn=bgn2;
  pos.position.end=end2;
  AppendVA_IntElementPos(pos_va,&pos);
  
  myerr=stderr;
  fprintf(stderr,"Accomplished the loading of the multialigns\n");
  newma = MergeMultiAlignsFast_new(sequenceDB, fragStore, pos_va, 0, 1, DP_Compare);
  
  fprintf(stderr,"Returned from MergeMultiAlignsFast_new\n");
  
  return 0;
}
