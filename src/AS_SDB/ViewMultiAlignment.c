
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
#include <unistd.h> /* man 3 getopt */

#include "AS_global.h" 
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"
#include "AS_SDB_SequenceDB.h"
#include "Array_CNS.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignStore_CNS.h"

FILE *myerr=NULL;
#if 0
int PrintMultiAlignT(FILE *out,MultiAlignT *ma,FragStoreHandle frag_store) {
  int depth;
  int rc,i;
  int window,length;
  // char **multia; 
  int **ids;
  char *consensus = Getchar(ma->consensus,0);
  char *quality = Getchar(ma->quality,0);
  char *nonblank;
  length = strlen(consensus);
  rc = MultiAlignT2Array(ma, frag_store, NULL, NULLFRAGSTOREHANDLE, &depth, &multia, &ids);
  if (rc) {
    for (window=0;window<=length;) {
      fprintf(out,"\n%d\r",window);
      fprintf(out,"\n");
      fprintf(out,"%-100.100s <-- cns\n",consensus+window);
      fprintf(out,"%-100.100s\n",quality+window);
      for (i=0;i<depth;i++) {
        if (multia[2*i] == NULL) continue;
        nonblank = strpbrk(multia[2*i]+window,"ACGT");
        if ( nonblank == NULL || nonblank-(multia[2*i]+window) > 100 ) continue;
        fprintf(out,"%-100.100s <- %10d\n",multia[2*i]+window,ids[i]+window);
        fprintf(out,"%-100.100s <- %10d\n",multia[2*i+1]+window,ids[i]+window);
      }
      window+=100;
    }
  } else {
    fprintf(stderr,"Error returned from IMP2Array.\n");
  }
  if (multia) {
    for (i=0;i<2*depth;i++) {
      free((char *)multia[i]);
    }
  }
  return 1;
}
#endif
int main(int argc, char *argv[]){
  
// args: seqStore revision_no frag_store contig_id
  int i;
  // int **id_array;
  // int depth;
  int show_qv=0;
  int dots=0;
  int num_tigs;
  int sequence_version=-1;
  int multialign_id=0;
  char SDBFilename[1024];
  char frgFilename[1024];
  int is_unitig=0;
  int ch;
  
  optarg = NULL;
  while ( ((ch = getopt(argc, argv, "h?dqcuS:v:F:i:")) != EOF)) {
    switch(ch) {
      case 'd':
        dots = 1;
        break;
      case 'q':
        show_qv=1;
        break;
      case 'c':
        is_unitig=0;
        break;
      case 'u':
        is_unitig=1;
        break;
      case 'S':
        strcpy(SDBFilename, optarg);
        break;
      case 'v':
        sequence_version=atoi(optarg);
        break;
      case 'i':
        multialign_id=atoi(optarg);
        break;
      case 'F':
        strcpy(frgFilename, optarg);
        break;
      case 'h':
      case '?':
        fprintf(stderr,"  Usage:\n\n");
        fprintf(stderr," %s [-c] [-u] [-q] [-d] -i MultiAlignID -S SDBStore -v SDB_version -F fragStore\n",argv[0]);
        fprintf(stderr,"\t\t-c\t(opt)\tmultialign_id refers to a contig (default)\n");
        fprintf(stderr,"\t\t-u\t(opt)\tmultialign_id refers to a unitig\n");
        fprintf(stderr,"\t\t-q\t(opt)\tdisplay fragment quality values\n");
        fprintf(stderr,"\t\t-d\t(opt)\tshow only differences from consensus\n");
        fprintf(stderr,"\t\t-i\t(req)\tmultialign_id\n");
        fprintf(stderr,"\t\t-S\t(req)\tSDBStore\n");
        fprintf(stderr,"\t\t-v\t(req)\tSDB_version number\n");
        fprintf(stderr,"\t\t-F\t(req)\tfragStore\n");
        exit(1);
      default: 
        fprintf(stderr,"  Incorrect usage. Please try %s -h\n",argv[0]);
        exit(1);
    }
  }
  if ( sequence_version == -1 ) { 
    fprintf(stderr," -v (SDB_version) flag required. See -h for usage.\n");
    exit(1);
  }
  if ( multialign_id == 0 ) { 
    fprintf(stderr," -i (multialign_id) flag required. See -h for usage.\n");
    exit(1);
  }
  
  { 
    tSequenceDB *sequenceDB = OpenSequenceDB(SDBFilename, FALSE, sequence_version);
    VA_TYPE(IntElementPos) *positions = CreateVA_IntElementPos(500);
    FragStoreHandle global_fragStore = openFragStore(frgFilename, "rb");
    MultiAlignT *ma =  LoadMultiAlignTFromSequenceDB(sequenceDB, multialign_id, is_unitig);
    MultiAlignT *newma;
    myerr=stderr;
    fprintf(stderr,"Accomplished the loading of the multialign from SeqStore\n");
    num_tigs = GetNumIntUnitigPoss(ma->u_list);
    for (i=0;i<num_tigs;i++) {
      IntUnitigPos *upos = GetIntUnitigPos(ma->u_list,i);
      IntElementPos pos;
      pos.type = AS_UNITIG;
      pos.ident = upos->ident;
      pos.position.bgn = upos->position.bgn;
      pos.position.end = upos->position.end;
      SetVA_IntElementPos(positions,i,&pos);
    }
    
    newma = MergeMultiAlignsFast_new(sequenceDB, NULLFRAGSTOREHANDLE, positions, 0, 1, NULL);
    //MultiAlignT2Array(ma, global_fragStore, NULLFRAGSTOREHANDLE, &depth, &multia, &id_array);
    //fprintf(stderr,"Converted into a character array\n");
    PrintMultiAlignT(myerr,newma,global_fragStore,NULL, NULLFRAGSTOREHANDLE,show_qv,dots,READSTRUCT_LATEST);
  } 
  fprintf(stderr,"You can stop the debugger here if helpful.\n");
  
  return 0;
}
