
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
int i;
uint32 contig_id=atoi(argv[4]),unitig_id;
FragStoreHandle fragStore;
tSequenceDB *sequenceDB = OpenSequenceDB(argv[1], FALSE, atoi(argv[2]));	
MultiAlignT *ma1 =  LoadMultiAlignTFromSequenceDB(sequenceDB, contig_id, FALSE);
MultiAlignT *ma2;
MultiAlignT *newma;
fragStore = openFragStore(argv[3], "rb");
unitig_id=GetIntUnitigPos(ma1->u_list,0)->ident;
ma2 =  LoadMultiAlignTFromSequenceDB(sequenceDB, unitig_id, TRUE);
PrintMultiAlignT(stderr,ma1,fragStore,NULL,NULLFRAGSTOREHANDLE,0,0,READSTRUCT_CNS);

myerr=stderr;
fprintf(stderr,"Accomplished the loading of the multialigns\n");
newma = ReplaceEndUnitigInContig(sequenceDB, fragStore, contig_id, unitig_id, 1, DP_Compare, NULL);

PrintMultiAlignT(stderr,newma,fragStore,NULL,NULLFRAGSTOREHANDLE,0,0,READSTRUCT_CNS);

fprintf(stderr,"Returned from ReplaceEndUnitigInContig\n");

}
