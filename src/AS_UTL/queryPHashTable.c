
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
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "AS_global.h"
#include "AS_UTL_HashCommon.h"
#include "AS_UTL_PHash.h"

int main(int argc, char **argv){
  PHashTable_AS *hashtable = OpenPHashTable_AS(argv[1]); 
  CDS_UID_t uid;

  while(fscanf(stderr, F_UID, &uid) != EOF){
    PHashValue_AS value;
    int ret = LookupInPHashTable_AS(hashtable,AS_UID_NAMESPACE, uid, &value);

    if(ret == HASH_SUCCESS){
      fprintf(stderr,"*+* Found uid " F_UID "\n", uid);
      fprintf(stderr,
              "* key:" F_UID " IID:" F_IID " del:%d type:%d refCount:%d \n",
	      uid, value.IID, value.deleted, value.type, value.refCount);
      fflush(stderr);
    }else{

      fprintf(stderr,"*?* Couldn't find uid " F_UID "\n", uid);
      fflush(stderr);
    }
  }
  return 0;
}
