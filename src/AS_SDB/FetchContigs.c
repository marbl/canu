
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
#include "MultiAlignStore_CNS.h"

/* Fetch two contigs from an SDB an output them to stdout in
   a format suitable for readying by test-laligner_ca
   
*/

int main(int argc, char *argv[]){
  
// args: seqStore revision_no  contig_id1 contig_id2 orientation
#if 0
  VA_TYPE(IntElementPos) *positions = CreateVA_IntElementPos(500);
  FragStoreHandle global_fragStore = openFragStore(argv[3], "rb");
#endif
  tSequenceDB *sequenceDB = OpenSequenceDB(argv[1], FALSE, atoi(argv[2]));	
  int id1 = atoi(argv[3]);
  int id2 = atoi(argv[4]);
  MultiAlignT *ma1;
  MultiAlignT *ma2;
  VA_TYPE(char) *seq1;
  VA_TYPE(char) *seq2;
  VA_TYPE(char) *qua;
  
  fprintf(stderr,"* Fetching %d\n", id1); fflush(NULL);
  ma1 =  LoadMultiAlignTFromSequenceDB(sequenceDB, id1, FALSE);
  fprintf(stderr,"* Fetching %d\n", id2); fflush(NULL);
  ma2 =  LoadMultiAlignTFromSequenceDB(sequenceDB, id2, FALSE);
  seq1 = CreateVA_char(10000);
  seq2 = CreateVA_char(10000);
  qua = CreateVA_char(10000);
  fprintf(stderr,"* Getting Sequence %d\n", id1); fflush(NULL);
  GetMultiAlignUngappedConsensus(ma1, seq1, qua); 
  fprintf(stdout,"%d > Sequence %d\n%s\njunk\n",
          id1, id1,Getchar(seq1,0));
  fprintf(stderr,"* Getting Sequence %d\n", id2); fflush(NULL);
  GetMultiAlignUngappedConsensus(ma2, seq2, qua); 
  fprintf(stdout,"%d > Sequence %d\n%s\njunk\n",
          id2, id2,Getchar(seq2,0));
  exit(0);
  
}
