
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

int main(int argc, char *argv[]){


 FragStoreHandle source;
 ReadStructp myRead;
 int clearRangeOnly = FALSE;
 int i;
 if(argc  < 2 ){
  fprintf(stderr,"Usage: %s <StorePath1> [anything]\nIf [anything], then only clear ranges of fragments will be displayed\n",
	  argv[0]);
  exit(1);
 }

 if(argc > 2){
   fprintf(stderr,"* Only clear ranges will be displayed!\n");
   clearRangeOnly = TRUE;
 }

 fprintf(stdout,"* Opening FragStore %s\n", argv[1]);
 source = openFragStore(argv[1],"r");

 if(source == NULLSTOREHANDLE){
   exit(1);
 }

 myRead =  new_ReadStruct();

 fprintf(stdout,"* Dumping fragStore %s (" F_S64 "," F_S64 ")\n",
	 argv[1],getFirstElemFragStore(source), getLastElemFragStore(source));

 while(scanf("%d",&i) == 1){
     getFragStore(source, i, FRAG_S_ALL, myRead);
     dump_ReadStruct(myRead, stdout, clearRangeOnly);
 }


   fprintf(stdout,"* Bye Bye\n");

   exit(0);
}
