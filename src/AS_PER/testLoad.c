
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


  FragStoreHandle target;
  FragStreamHandle jane;
  ReadStructp myRead =  new_ReadStruct();
#if 0
  FragStoreHandle source;
  char fileName1[1024];
#endif


#if 0
  if(argc  >= 2 ){
    sprintf(fileName1,"%s",argv[1]);
    fprintf(stderr,"* Creating memory based target from %s\n", fileName1);
    source = openFragStore(fileName1,"r");
    target = createFragStore(NULL, "Proto0.1", 1);
  }

  concatFragStore(target, source);
#else
  target = loadFragStore(argv[1]);
#endif

  fprintf(stderr,"* After concatting source\n");
  fflush(stderr);

    jane = openFragStream(target,NULL,0);

    while(nextFragStream(jane, myRead, FRAG_S_ALL)){
      dump_ReadStruct(myRead, stderr, FALSE);
    }

    closeFragStream(jane);
    fprintf(stderr,"* Closing jane, source, target\n");
    fprintf(stderr,"* Bye Bye\n");
    exit(0);


  fprintf(stderr,"***GORK***\n");
  fflush(stderr);

  closeFragStore(target);
  //  closeFragStore(source);

  return 0;
}

