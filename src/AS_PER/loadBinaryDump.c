
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
#include "AS_PER_genericStore.h"
#include "AS_PER_fragStore.h" 
#include "AS_PER_fragStore_private.h" 
#include "AS_PER_distStore.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA.h"

int main(int argc, char *argv[]){


 FragStoreHandle target;
 ReadStructp myRead;
 char buffer[2048];
 FILE *infp;
 VA_TYPE(char) *sequence = CreateVA_char(2048);
 VA_TYPE(char) *source = CreateVA_char(2048);
 ShortFragRecord fr;

 if(argc  != 2 ){
  fprintf(stderr,"Usage: %s <DumpFileName> \n",
	  argv[0]);
  exit(1);
 }

   fprintf(stdout,"* Opening dumpfile %s\n", argv[1]);
   infp = fopen(argv[1],"r");

   AssertPtr(infp);


 sprintf(buffer,"%s.Store", argv[1]);

 target = createFragStore(buffer,argv[1],1);
 myRead =  new_ReadStruct();

 
 fprintf(stderr,"* Starting load\n");
 fflush(NULL);

 while(loadDumpFragRecord(infp, &fr, sequence, source)){
   appendDumpToFragStore(target, &fr, sequence, source);
 }


 fclose(infp);
 closeFragStore(target);
   exit(0);

}
