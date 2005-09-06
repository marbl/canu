
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
#include "AS_PER_gkpStore.h"

int usage(char *pgmname){
  fprintf(stderr,"Usage: %s [-b firstIID] [-e lastIID] [-U] [-D]<FrgStorePath> <GkpStorePath>\n"
	  "\t-U\tcauses UIDs to be printed\n"
	  "\t-D\tcauses info for deleted fragments to be printed\n" 
	  ,
	    pgmname);
    return 0;
}

int main(int argc, char *argv[]){

  FragStoreHandle source;
  ReadStructp myRead;
  int32 begin = -1, end = -1;
  int i;
  int ch;
  int printUIDs=0;
  uint32 clr_bgn,clr_end;
  GateKeeperStore gkpStore;
  GateKeeperFragmentRecord gkpFrag;
  int includeDeleted=0;


  while ((ch = getopt(argc, argv, "b:e:UD")) != EOF){
    switch(ch) {
    case 'D':
      includeDeleted=1;
      break;
    case 'e':
      end = atoi(optarg);
      fprintf(stderr,"* end = %d\n", end);
      break;
    case 'b':
      begin = atoi(optarg);
      fprintf(stderr,"* begin = %d\n", begin);
      break;
    case 'U':
      printUIDs=1;
      break;
    default:
      fprintf(stderr,"* Unknown option %s\n", optarg);
      usage(argv[0]);
      exit(-1);
      break;
    }
  }

  if(argc-optind!=2){
    usage(argv[0]);
    exit(1);
  }

  fprintf(stderr,"* Opening FragStore %s\n", argv[optind]);
  source = openFragStore(argv[optind++],"r");
  if(source == NULLSTOREHANDLE){
    exit(1);
  }

  myRead =  new_ReadStruct();

  InitGateKeeperStore(&gkpStore,argv[optind]);
  OpenReadOnlyGateKeeperStore(&gkpStore);
  //  assert(TestOpenGateKeeperStore(&gkpStore) == TRUE);

  fprintf(stderr,"* Dumping clear range lengths from fragStore %s (%d,%d) of (" F_S64 "," F_S64 ")\n",
	  argv[optind],begin,end,
	  getFirstElemFragStore(source), getLastElemFragStore(source));

  if(begin < getFirstElemFragStore(source) || begin > getLastElemFragStore(source)){
    begin = getFirstElemFragStore(source);
  }
  if(end > getLastElemFragStore(source) || end < getFirstElemFragStore(source)){
    end = getLastElemFragStore(source);
  }
  for(i = begin; i <= end; i++){
    if(getGateKeeperFragmentStore(gkpStore.frgStore,i,&gkpFrag)!=0)
      assert(0);
    if(gkpFrag.deleted&&!includeDeleted)continue;
    getFragStore(source, i, FRAG_S_ALL, myRead);
    getClearRegion_ReadStruct(myRead, &clr_bgn,&clr_end, READSTRUCT_LATEST);
    if(printUIDs){
      if(includeDeleted){
	printf(F_UID "\t%d\t%d\n",gkpFrag.readUID,clr_end-clr_bgn,gkpFrag.deleted);
      } else {
	printf(F_UID "\t%d\n",gkpFrag.readUID,clr_end-clr_bgn);
      }
    } else {
      if(includeDeleted){
	printf("%d\t%d\t%d\n",i,clr_end-clr_bgn,gkpFrag.deleted);
      } else {
	printf("%d\t%d\n",i,clr_end-clr_bgn);
      }
    }
  }

  fprintf(stderr,"* Bye Bye\n");

  exit(0);
}
