
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
#include "AS_global.h"
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"


int main(int argc, char *argv[]){

  char fileName1[1024];
  char fileName2[1024];

  FragStoreHandle source, target;
  FragStreamHandle jane;
  ReadStructp myRead;


  if(argc  == 2 ){
    sprintf(fileName1,"1%s",argv[1]);
    fprintf(stderr,"* Creating memory based target\n");
    source = createFragStore(fileName1, "Proto0.0", 1);
    target = createFragStore(NULL, "Proto0.1", 1);
  }else if(argc == 3){
    sprintf(fileName1,"1%s",argv[1]);
    sprintf(fileName2,"2%s",argv[1]);
    source = createFragStore(fileName1, "Proto0.1", 1);
    fprintf(stderr,"* Creating file-based target %s\n", fileName2);
    target = createFragStore(fileName2, "Proto0.1", 1);
  }else{
    fprintf(stderr,"Usage: %s <StorePath1> [StorePath2]\nExiting...\n",
	    argv[0]);
    exit(1);
  }

  myRead =  new_ReadStruct();

  setAccID_ReadStruct(myRead, 3);
  setReadType_ReadStruct(myRead, AS_READ);
  setReadIndex_ReadStruct(myRead, 3);
  setSequence_ReadStruct(myRead, "aaccggttaaccggtt", "0000000000000000");
  setSource_ReadStruct(myRead, "Go Home");

  appendFragStore(source,myRead);

  setAccID_ReadStruct(myRead, 5);
  setReadType_ReadStruct(myRead, AS_EXTR);
  setReadIndex_ReadStruct(myRead, 5);
  setSequence_ReadStruct(myRead, "acgtacgt", "00000000");
  setSource_ReadStruct(myRead, "external read");
  setClearRegion_ReadStruct(myRead,20,30,READSTRUCT_CGW);
  setClearRegion_ReadStruct(myRead,10,20,READSTRUCT_OVL); // ovl overwrites cgw
  appendFragStore(source,myRead);

  setAccID_ReadStruct(myRead, 6);
  setReadType_ReadStruct(myRead, AS_TRNR);
  setReadIndex_ReadStruct(myRead, 6);
  setSequence_ReadStruct(myRead, "acgtacgtggggg", "0000000000000");
  setSource_ReadStruct(myRead, "transposon library read");
  appendFragStore(source,myRead);

  commitFragStore(source);
#if 0
  closeFragStore(source);
  source = openFragStore(fileName1,"rb");
#endif
  
#if 1
  concatFragStore(target, source);

  fprintf(stderr,"* After concatting source\n");
  fflush(stderr);

  if(argc < 2){
    jane = openFragStream(target,NULL,0);

    while(nextFragStream(jane, myRead, FRAG_S_ALL)){
      dump_ReadStruct(myRead, stderr, FALSE);
    }

    closeFragStream(jane);
    fprintf(stderr,"* Closing jane, source, target\n");
    fprintf(stderr,"* Bye Bye\n");
    exit(0);
  }

  fprintf(stderr,"***GORK***\n");
  fflush(stderr);

  closeFragStore(target);
  closeFragStore(source);

  fprintf(stderr,"*** Closed source and target \n");
  fflush(stderr);
  fprintf(stderr,"* Opening target %s\n", fileName2);
  target = openFragStore(fileName2, "rb");
#else
  fprintf(stderr,"* Opening target %s\n", fileName1);
  target = openFragStore(fileName1, "rb");
#endif
 
  fprintf(stderr,"\n\n* Contents of %s\n", fileName2);
  /* FragStore_get(source, 1, FRAG_S_ALL, myRead); */
  jane = openFragStream(target,NULL,0);

  while(nextFragStream(jane, myRead, FRAG_S_ALL)){
    dump_ReadStruct(myRead, stderr, FALSE);
  }

  closeFragStore(target);
  closeFragStream(jane);
  return 0;
}

