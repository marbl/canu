
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
static char CM_ID[] = "$Id: BuildSequenceDB.c,v 1.5 2005-09-15 15:20:16 eliv Exp $";

//#define DEBUG 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_PER_SafeIO.h"
#include "AS_SDB_SequenceDB.h"
#include "AS_MSG_pmesg.h"


void usage(void){
  fprintf(stderr,"* usage: BuildSequenceDB <ium file> <storeName> \n");
  exit(1);
}

void ReadInput(tSequenceDB *sequenceDB, FILE *file, int check){
  GenericMesg   *pmesg;
  MesgReader reader = (MesgReader)InputFileType_AS(file);
  int32 totalFrags = 1;

  while(  (EOF != (reader)(file, &pmesg))){

    switch(pmesg->t){

    case MESG_IDT:
      break;
    case MESG_UOM:  // UnitigOverlapMesg
      break;
    case MESG_IUM:  // IntUnitigMesg
      {
	IntUnitigMesg *ium_mesg = (IntUnitigMesg *)pmesg->m;
	MultiAlignT *ma = CreateMultiAlignTFromIUM(ium_mesg, totalFrags, FALSE);
	totalFrags += ium_mesg->num_frags;
	if(check){
	  MultiAlignT *loadma = LoadMultiAlignTFromSequenceDB(sequenceDB, ium_mesg->iaccession,TRUE);
	  if(CompareMultiAlignT(ma, loadma)){
	    fprintf(stderr,"* MultiAlignment %d differ\n",
		    ium_mesg->iaccession);
	  }
	  DeleteMultiAlignT(loadma);
	  DeleteMultiAlignT(ma);
	}else{
	  InsertMultiAlignTInSequenceDB(sequenceDB, ium_mesg->iaccession, TRUE, ma, FALSE);
	}
      }
      break;


    case MESG_ADT:
      break;
    case MESG_ILK:

    default:
      break;
    }
    
  }

}

int main(int argc, char **argv){
  char *inputName, *storeName;
  FILE *inputFile;
  tSequenceDB *sequenceDB;
  if(argc < 3)
    usage();


  inputName = argv[1];
  storeName = argv[2];
  
  inputFile = fopen(inputName, "r");
  if(!inputFile){
    fprintf(stderr,"* Couldn't open file %s ...bye\n", inputName);
    exit(1);
  }

  sequenceDB = CreateSequenceDB(storeName, 0, FALSE);

  fprintf(stderr,"* Calling ReadInput *\n");
  ReadInput(sequenceDB, inputFile, FALSE);

  fprintf(stderr,"* Calling SaveSequenceDB *\n");
  SaveSequenceDB(sequenceDB);
  SaveSequenceDB(sequenceDB);

  fprintf(stderr,"* Calling DeleteSequenceDB *\n");
  DeleteSequenceDB(sequenceDB);

  rewind(inputFile);

  sequenceDB = OpenSequenceDB(storeName,FALSE,1);
  ReadInput(sequenceDB, inputFile, TRUE);

  rewind(inputFile);
  ReadInput(sequenceDB, inputFile, TRUE);

  rewind(inputFile);
  ClearCacheSequenceDB(sequenceDB,TRUE);
  ReadInput(sequenceDB, inputFile, TRUE);



  fprintf(stderr,"* Calling DeleteSequenceDB *\n");
  DeleteSequenceDB(sequenceDB);


  exit(0);

}
