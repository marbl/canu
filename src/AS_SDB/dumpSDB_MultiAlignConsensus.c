
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
#include "cds.h" 
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"
#include "AS_SDB_SequenceDB.h"
#include "Array_CNS.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignStore_CNS.h"

FILE *myerr=NULL;

int main(int argc, char *argv[]){

// args: seqStore revision_no frag_store contig_id
  int i;
  int sequence_version=-1;
  int multialign_id=NULLINDEX;
  char SDBFilename[1024];
  int is_unitig=0;
  int reverse=0; 
  int len;
  char ch;

  optarg = NULL;
  while ( ((ch = getopt(argc, argv, "h?cuS:v:i:r")) != EOF)) {
    switch(ch) {
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
    case 'r':
      reverse=1;
      break;
    case 'i':
      multialign_id=atoi(optarg);
      break;
    case 'h':
    case '?':
      fprintf(stderr,"  Usage:\n\n");
      fprintf(stderr," %s [-c] [-u] [-q] [-d] -i MultiAlignID -S SDBStore -v SDB_version -F fragStore\n",argv[0]);
      fprintf(stderr,"\t\t-c\t(opt)\tmultialign_id refers to a contig (default)\n");
      fprintf(stderr,"\t\t-u\t(opt)\tmultialign_id refers to a unitig\n");
      fprintf(stderr,"\t\t-r\t(opt)\treverse-complement multialign_id\n");
      fprintf(stderr,"\t\t-i\t(req)\tmultialign_id\n");
      fprintf(stderr,"\t\t-S\t(req)\tSDBStore\n");
      fprintf(stderr,"\t\t-v\t(req)\tSDB_version number\n");
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
  if ( multialign_id == NULLINDEX ) { 
    fprintf(stderr," -i (multialign_id) flag required. See -h for usage.\n");
    exit(1);
  }

  { 
    tSequenceDB *sequenceDB = OpenSequenceDB(SDBFilename, FALSE, sequence_version);
    VA_TYPE(char) *ungappedCNS = CreateVA_char(256*1024);
    VA_TYPE(char) *ungappedQLT = CreateVA_char(256*1024);
    MultiAlignT *ma =  LoadMultiAlignTFromSequenceDB(sequenceDB, multialign_id, is_unitig);
    char *sequence;
    GetMultiAlignUngappedConsensus(ma, ungappedCNS, ungappedQLT);
    sequence=Getchar(ungappedCNS,0);
    if(reverse){
      SequenceComplement( sequence, NULL);
    }

    len = strlen(sequence);
    fprintf(stdout,">%s %d len %d\n",(is_unitig ? "unitig" : "contig"),multialign_id,len);
    for(i=0;i<len;i+=60){
      fprintf(stdout,"%.60s\n",sequence+i);
    }    

  } 
  fprintf(stderr,"You can stop the debugger here if helpful.\n");

  return 0;
}
