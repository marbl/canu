
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
  IntScreenMatch match[5];;


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

  match[0].next = NULL;
  match[0].where.bgn = 1;
  match[0].where.end = 100;
  match[0].iwhat = 3;
  match[0].repeat_id = 12;
  match[0].relevance = 9;
  match[0].portion_of.bgn = 1;
  match[0].portion_of.end = 100;
  match[0].direction = AS_FORWARD;
  match[0].next = match + 1;
  match[1].where.bgn = 100;
  match[1].where.end = 1000;
  match[1].iwhat = 5;
  match[1].repeat_id = 12;
  match[1].relevance = 9;
  match[1].portion_of.bgn = 10;
  match[1].portion_of.end = 1000;
  match[1].direction = AS_REVERSE;
  match[1].next = NULL;

  myRead =  new_ReadStruct();

  setAccID_ReadStruct(myRead, 1);
  setReadIndex_ReadStruct(myRead, 1);
  setReadType_ReadStruct(myRead, AS_UBAC);
  setSequence_ReadStruct(myRead, "aaccaaccaaccaacc", NULL);
  setSource_ReadStruct(myRead, "Just testing...sucker");
  setScreenMatches_ReadStruct(myRead, 1, match);
  setEntryTime_ReadStruct(myRead,time(0));
  setClearRegion_ReadStruct(myRead,1,14,READSTRUCT_ORIGINAL);
  setLocalePos_ReadStruct(myRead,99,100);
  setLocID_ReadStruct(myRead,0x0ffffffffffffffeLL);
  appendFragStore(source,myRead);

  setAccID_ReadStruct(myRead, 2);
  setReadType_ReadStruct(myRead, AS_EBAC);
  setReadIndex_ReadStruct(myRead, 2);
  setSequence_ReadStruct(myRead, "aaccggttaaccggtt", "0000000000000000");
  setSource_ReadStruct(myRead, "Go Home");
  setLocalePos_ReadStruct(myRead,100,101);
  setLocID_ReadStruct(myRead,0x0fffffffffffffffLL);
  setScreenMatches_ReadStruct(myRead, 2, match);
  appendFragStore(source,myRead);

  setAccID_ReadStruct(myRead, 3);
  setReadType_ReadStruct(myRead, AS_READ);
  setReadIndex_ReadStruct(myRead, 3);
  setSequence_ReadStruct(myRead, "aaccggttaaccggtt", "0000000000000000");
  setSource_ReadStruct(myRead, "Go Home");
  setLocalePos_ReadStruct(myRead,102,103);
  setScreenMatches_ReadStruct(myRead, 1, &match[1]);
  setLocID_ReadStruct(myRead,0x0ffffffffffffffaLL);

  appendFragStore(source,myRead);

  setAccID_ReadStruct(myRead, 4);
  setReadIndex_ReadStruct(myRead, 4);
  setReadType_ReadStruct(myRead, AS_UBAC);
  setSequence_ReadStruct(myRead, "aaccaaccaaccaacc", NULL);
  setSource_ReadStruct(myRead, "Just testing...sucker");
  setEntryTime_ReadStruct(myRead,time(0));
  setClearRegion_ReadStruct(myRead,1,14,READSTRUCT_ORIGINAL);
  setLocalePos_ReadStruct(myRead,99,100);
  setLocID_ReadStruct(myRead,0x0ffffffffffffffeLL);
  setScreenMatches_ReadStruct(myRead, 2, &match[0]);
  // ADDED BY JASON - OCT 2001
  setClearRegion_ReadStruct(myRead,10,20,READSTRUCT_OVL);
  setClearRegion_ReadStruct(myRead,20,30,READSTRUCT_CGW); // both get saved
  appendFragStore(source,myRead);

  setAccID_ReadStruct(myRead, 5);
  setReadType_ReadStruct(myRead, AS_EXTR);
  setReadIndex_ReadStruct(myRead, 5);
  setSequence_ReadStruct(myRead, "acgtacgt", "00000000");
  setSource_ReadStruct(myRead, "external read");
  setLocalePos_ReadStruct(myRead,102,103);
  setScreenMatches_ReadStruct(myRead, 1, &match[1]);
  setLocID_ReadStruct(myRead,0x0ffffffffffffffaLL);
  // ADDED BY JASON - OCT 2001
  setClearRegion_ReadStruct(myRead,20,30,READSTRUCT_CGW);
  setClearRegion_ReadStruct(myRead,10,20,READSTRUCT_OVL); // ovl overwrites cgw
  appendFragStore(source,myRead);

  setAccID_ReadStruct(myRead, 6);
  setReadType_ReadStruct(myRead, AS_TRNR);
  setReadIndex_ReadStruct(myRead, 6);
  setSequence_ReadStruct(myRead, "acgtacgtggggg", "0000000000000");
  setSource_ReadStruct(myRead, "transposon library read");
  setLocalePos_ReadStruct(myRead,102,103);
  setScreenMatches_ReadStruct(myRead, 1, &match[0]);
  setLocID_ReadStruct(myRead,0x0ffffffffffffffaLL);
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

