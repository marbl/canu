
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

#include "AS_global.h" 
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"
#include "AS_SDB_SequenceDB.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignStore_CNS.h"
#include "AS_UTL_rand.h"

FILE *myerr=NULL;
int main(int argc, char *argv[]){
  
// args: id1 bgn end id2 bgn end
  tSequenceDB *sequenceDB;
  MultiAlignT *ma =  CreateEmptyMultiAlignT();
  char *sdbName = argv[1];
  int sdbRev = atoi(argv[2]);
  int32 count = 0;
  int i;
  int firstID = -1;
  int lastID = -1;
  int iter; 
  fprintf(stderr,"* Opening %s rev %d\n", sdbName, sdbRev);
  
  InitRandom_AS(0);
  sequenceDB = OpenSequenceDB(sdbName, FALSE, sdbRev);
  
  firstID = 1;
  lastID = GetNumtMARecords(sequenceDB->Unitigs)-1;
  
  
  for(iter = 0; iter < 10; iter++){
    time_t t = time(0);
    fprintf(stderr,"* Iteration %d %s\n",iter, ctime(&t));
    t = time(0);
    
    fprintf(stderr,"* Sequential access for ids [%d,%d]\n",
            firstID, lastID);
    fflush(stderr);
    
    for(i = firstID; i < lastID; i++){
      if(i%100000 == 0){
        fprintf(stderr,".");
        if(i%1000000 == 0){
          t = time(0);
          fprintf(stderr,"\n%s %d ",ctime(&t),i);
        }
        fflush(stderr);
      }
      ReLoadMultiAlignTFromSequenceDB(sequenceDB, ma, i, TRUE);    // This will load the ma from the disk file
    }
    
    fprintf(stderr,"\n* Random access for ids [%d,%d] %s\n",
            firstID, lastID, ctime(&t));
    fflush(stderr);
    
    for(i = firstID; i < lastID; i++){
      int32 id = GetRand_AS(firstID,lastID,TRUE);
      count++;
      if(count%10000 == 0){
        fprintf(stderr," %d ", id);
        if(count%1000000 == 0){
          t = time(0);
          fprintf(stderr,"\n%s %d ",ctime(&t), count);
        }
        fflush(stderr);
      }
      ReLoadMultiAlignTFromSequenceDB(sequenceDB, ma, id, TRUE);    // This will load the ma from the disk file
    }
    
    
  }
  
  fprintf(stderr,"* Successful completion ... bye\n");
  return 0;
}
