
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
static char CM_ID[] = "$Id: computeGCcontent.c,v 1.1 2004-09-23 20:32:56 mcschatz Exp $";


#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>

#include "AS_global.h"
#include "MicroHetREZ.h"
#include "MicroHetScoreREZ.h"
#include "MicroHetPartitionsREZ.h"
#include "MicroHetInterfaceREZ.h"

#include "UtilsREZ.h"
#include "AS_UTL_skiplist.h"
#include "AS_UTL_Var.h"

#include "PrimitiveVA_MSG.h"



#define DEBUG -1



#define USAGETROUBLE 0 /* 0 used to cause assert to fail */
int main (int argc, char *argv[]) {
  FragStoreHandle storeHandle = 0;
  GateKeeperStore gkpStore;
  GateKeeperFragmentRecord gkpFrag;
  ReadStructp fsread=new_ReadStruct();
  char           *storeName;
  char           *gkpName;
  char           *fileName;
  char *qul,*seq,*toprint,*src;
  int alloclen=1000;
  int alloclen2=1000;
  CDS_IID_t iid,lastIID;
  uint32 clr_bgn,clr_end;
  FILE* input;
  CDS_UID_t gc=0;
  CDS_UID_t L=0;

  if((seq=(char*)malloc(sizeof(char)*alloclen)) == NULL ||
     (qul=(char*)malloc(sizeof(char)*alloclen)) == NULL ||
     (toprint=(char*)malloc(sizeof(char)*alloclen)) == NULL ||
     (src=(char*)malloc(sizeof(char)*alloclen2)) == NULL)
    assert(0);

  if(argc!=4){
    fprintf(stderr,"Usage: %s frgStore gkpStore list_of_fragIIDs\n",argv[0]);
    exit(-1);
  }

  storeName = strdup(argv[1]);
  gkpName = strdup(argv[2]);

  lastIID = atoi(argv[3]);


  InitGateKeeperStore(&gkpStore,gkpName);
  assert(existsFragStore(storeName) == TRUE);
  //  assert(TestOpenGateKeeperStore(&gkpStore) == TRUE);
  storeHandle = openFragStore(storeName,"r");
  OpenReadOnlyGateKeeperStore(&gkpStore);

  L=0;
  gc=0;
  for(iid=1;iid<=lastIID;iid++){
    if(getFragStore(storeHandle,iid,FRAG_S_ALL,fsread)!=0){
      fprintf(stderr,"Couldn't get fragment from frgStore for iid " F_IID "\n",
              iid);
      continue;
    }

    if(getGateKeeperFragmentStore(gkpStore.frgStore,iid,&gkpFrag)!=0)
      assert(0);

    if(gkpFrag.deleted)continue;

    getClearRegion_ReadStruct(fsread, &clr_bgn,&clr_end, READSTRUCT_LATEST);
    while(getSequence_ReadStruct(fsread,seq,qul,alloclen)!=0){
      alloclen*=2;
      seq=(char*)realloc(seq,alloclen*sizeof(char));
      qul=(char*)realloc(qul,alloclen*sizeof(char));
      toprint=(char*)realloc(toprint,alloclen*sizeof(char));
    }

    strcpy(toprint,seq+clr_bgn);
    toprint[clr_end-clr_bgn]='\0';
    { 
      int i,n;
      n=clr_end-clr_bgn;
      L+=n;
      for(i=0;i<n;i++){
	if(toprint[i] == 'c' || toprint[i] == 'C' ||
	   toprint[i] == 'g' || toprint[i] == 'G')
	  gc++;
      }
    }
  }
  printf("GC content %f ( " F_U64 " / " F_U64 " )\n",(double)gc/(double)L,gc,L);
  closeFragStore(storeHandle);
  free(seq);
  free(qul);
  free(toprint);
  exit(0);
}
