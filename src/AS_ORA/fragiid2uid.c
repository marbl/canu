
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
static char CM_ID[] = "$Id: fragiid2uid.c,v 1.2 2004-09-23 20:25:25 mcschatz Exp $";


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


int main (int argc, char *argv[]) {
  FragStoreHandle storeHandle = 0;
  GateKeeperStore gkpStore;
  GateKeeperFragmentRecord gkpFrag;
  ReadStructp fsread=new_ReadStruct();
  char           *storeName;
  char           *gkpName;
  char           *fileName;
  CDS_IID_t iid;
  FILE* input;

  if(argc!=4){
    fprintf(stderr,"Usage: %s frgStore gkpStore iidlist\n",argv[0]);
    exit(-1);
  }

  { /* Parse the argument list using "man 3 getopt". */
    storeName = strdup(argv[1]);
    gkpName = strdup(argv[2]);
    fileName  = strdup(argv[3]);
  }

  InitGateKeeperStore(&gkpStore,gkpName);
  storeHandle = openFragStore(storeName,"r");
  OpenReadOnlyGateKeeperStore(&gkpStore);
  input = fopen (fileName,"r");
  if(input==NULL||!existsFragStore(storeName)){
    fprintf(stderr,"Usage: %s frgStore gkpStore iidlist\n",argv[0]);
    exit(-1);
  }

  assert(input != NULL);
  assert(existsFragStore(storeName) == TRUE);
  assert(TestOpenGateKeeperStore(&gkpStore) == TRUE);

  fprintf(stderr,"%s output format:\niid,uid,status(1=deleted)\n",argv[0]);

  while(fscanf(input,F_IID,&iid)!=EOF){
    //    fprintf(stderr,"Trying to handle %d\n",iid);
    if(getFragStore(storeHandle,iid,FRAG_S_ALL,fsread)!=0){
      fprintf(stderr,"Couldn't get fragment from frgStore for iid %d\n",iid);
      continue;
    }
    if(getGateKeeperFragmentStore(gkpStore.frgStore,iid,&gkpFrag)!=0)
      assert(0);
    printf(F_IID " " F_UID " %d\n",iid,(gkpFrag).readUID,(gkpFrag).deleted);
  }
  fclose(input);
  closeFragStore(storeHandle);
  exit(0);
}
