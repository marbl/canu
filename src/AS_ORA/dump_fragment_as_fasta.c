
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
static char CM_ID[] = "$Id: dump_fragment_as_fasta.c,v 1.2 2004-09-23 20:25:25 mcschatz Exp $";


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
#undef BEG_END
#ifdef BEG_END
  int beg,end;
#endif
  CDS_IID_t iid;
  uint32 clr_bgn,clr_end;
  FILE* input;

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
  fileName  = strdup(argv[3]);

  InitGateKeeperStore(&gkpStore,gkpName);
  assert(existsFragStore(storeName) == TRUE);
  //  assert(TestOpenGateKeeperStore(&gkpStore) == TRUE);
  storeHandle = openFragStore(storeName,"r");
  OpenReadOnlyGateKeeperStore(&gkpStore);
  input = fopen (fileName,"r");
  assert(input != NULL);

#ifdef BEG_END
  while(fscanf(input,F_IID " %d %d",&iid,&beg,&end)!=EOF){
#else
  while(fscanf(input,F_IID,&iid)!=EOF){
#endif
    //    fprintf(stderr,"Trying to handle " F_IID "\n",iid);
    if(getFragStore(storeHandle,iid,FRAG_S_ALL,fsread)!=0){
      fprintf(stderr,"Couldn't get fragment from frgStore for iid " F_IID "\n",
              iid);
      continue;
    }
    getClearRegion_ReadStruct(fsread, &clr_bgn,&clr_end, READSTRUCT_LATEST);
    fprintf(stderr,
            "Obtained fragment for iid " F_IID "; clear range = %u -- %u\n",
            iid,clr_bgn,clr_end);
    while(getSequence_ReadStruct(fsread,seq,qul,alloclen)!=0){
      alloclen*=2;
      seq=(char*)realloc(seq,alloclen*sizeof(char));
      qul=(char*)realloc(qul,alloclen*sizeof(char));
      toprint=(char*)realloc(toprint,alloclen*sizeof(char));
    }

    /* DON'T USE FrgStore SOURCE INFORMATION!!
    while(getSource_ReadStruct(fsread,src,alloclen2)!=0){
      alloclen2*=2;
      src=realloc(seq,alloclen2*sizeof(char));
    }*/

    //    sscanf(src,"??\n[%d,%d]\n",&beg,&end);
    //    fprintf(stderr,"Src: %s\n",src);
    //    fprintf(stderr,"Parsed Src: beg=%d end=%d\n",beg,end);
    strcpy(toprint,seq+clr_bgn);
    toprint[clr_end-clr_bgn]='\0';
    { 
      if(getGateKeeperFragmentStore(gkpStore.frgStore,iid,&gkpFrag)!=0)
        assert(0);
#ifdef BEG_END
      printf(">" F_IID "_%d_%d " F_UID "\n%s\n",
             iid,beg,end,(gkpFrag).readUID,toprint);
#else
      printf(">" F_IID " " F_UID "\n%s\n",iid,(gkpFrag).readUID,toprint);
#endif
    }
  }

  fclose(input);
  closeFragStore(storeHandle);
  free(seq);
  free(qul);
  free(toprint);
  exit(0);
}
