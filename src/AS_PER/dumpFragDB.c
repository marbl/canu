
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
#include "AS_PER_fragDB.h"


int main(int argc, char *argv[]){


 ReadStructp myRead;
 int i;

 if(argc  < 2 ){
  fprintf(stderr,"Usage: %s <StorePath1> \n",
	  argv[0]);
  exit(1);
 }

 myRead =  new_ReadStruct();

 fprintf(stdout,"* Dumping fragDB\n");


 for(i = 1; i < 1000; i++){
     getFragDB(i, FRAG_S_ALL, myRead, 1);
     dump_ReadStruct(myRead, stdout, FALSE);
 }

   exit(0);
}
