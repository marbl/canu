
/**************************************************************************
 * This file is based on 
 * ... the Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * 
 * The Celera Assembler is free software; you can redistribute it and/or modify
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

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _OSF_SOURCE
#include <sys/mode.h>
#endif
#include <unistd.h>
#include <dirent.h>
#include "assert.h"
#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_Hash.h"
#include "AS_UTL_ID_store.h"
#include "PrimitiveVA.h"
#include "PrimitiveVA_MSG.h"
#include "MultiAlignStore_CNS.h"

#define FASTA_SEQ_LINE_LENGTH 60

int main(int argc, char *argv[])
{ GenericMesg *pmesg;
  SnapConConMesg *contig;
  MesgReader   reader;

  if ( argc> 1 ) {
     fprintf(stderr,"Usage: %s < asmfile > contigs_fasta_file\n",argv[0]);
     exit(1);
  }
  
  reader = InputFileType_AS( stdin );

  while (reader(stdin,&pmesg) != EOF){
    if (pmesg->t ==MESG_CCO)  {
      int intoline=0;
      int i;
      contig = pmesg->m;
      printf(">" F_S64 " contig\n",contig->eaccession);
      assert ( strlen(contig->consensus) == contig->length);
      for(i=0;i<contig->length;i++){
	if(contig->consensus[i] != '-'){
	  if(intoline==FASTA_SEQ_LINE_LENGTH){
	    intoline=0;
	    printf("\n");
	  }
	  printf("%c",contig->consensus[i]);
	  intoline++;
	}
      }
      printf("\n");
    }
 }
 exit (0);
}
