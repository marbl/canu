
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

/**********************************************************************
 Module:
 Description:
 Assumptions:
**********************************************************************/

static char CM_ID[] = "$Id: MicroHetIUMdiffsPrint.c,v 1.5 2006-10-08 08:47:40 brianwalenz Exp $";


#include <unistd.h> /* man 3 getopt */
#include "MicroHetREZ_test3.h"
//#include "MicroHetScoreREZ_test3.h"
//#include "MicroHetPartitionsREZ_test3.h"
#include "MicroHetInterfaceREZ_test3.h"

//#include "UtilsREZ.h"
//#include "AS_UTL_skiplist.h"
//#include "AS_UTL_Var.h"

//#include "PrimitiveVA_MSG.h"

#include <assert.h>
//#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h> /* man 3 getopt */
//#include <string.h>
//#include <time.h>
//#include <ctype.h>
//#include <math.h>


typedef enum {
  PRINT_DOTS,
  PRINT_SPLITS
} PrintTypes;

void  AS_REZ_print_informatives_alignment(Alignment_t *,int);

void AS_REZ_print_IUM_diffs(IntUnitigMesg* ium, FragStoreHandle handle,int printwhat)
{
  Alignment_t *ali;  
  ali = AS_REZ_convert_IUM_to_alignment(ium,handle,(tFragStorePartition *)NULL,FALSE);
  switch(printwhat){
  case PRINT_DOTS:
    AS_REZ_print_informatives_alignment(ali,ium->num_frags);
    break;
  case PRINT_SPLITS:
    AS_REZ_print_informative_splits(ali,ium->num_frags);
    break;
  default:
    fprintf(stderr,"AS_REZ_print_IUM_diffs doesn't know what you want to print\n");
    AS_REZ_free_alignment(ali);
    exit(-1);
  }
  AS_REZ_free_alignment(ali);
}

void Usage(char *pn){
  fprintf(stderr,
	  "USAGE: %s [-a a-stat-cutoff] [-S] frgstoreName utgFile\n"
	  "\t-a a-stat-cutoff: process only unitigs with a-stat < this value (default = 0)\n"
	  "\t-S: print splits as unresolved trees (default = dot-diff. data matrix)\n",
	  pn);
}

int main (int argc, char *argv[]) {
  GenericMesg   *pmesg;  
  IntUnitigMesg *iunitig;
  MessageType    imesgtype;
  MesgReader     reader;
  FragStoreHandle storeHandle = 0;
  char           *storeName;
  char           *fileName;
  FILE* input;
  PrintTypes printwhat = PRINT_DOTS;

  float cthresh=0;
  int ch,errflg=0;

  while (!errflg && ((ch = getopt(argc, argv,
		       "a:S")) != EOF)){
      switch(ch) {
      case 'S':
	printwhat = PRINT_SPLITS;
	break;
      case 'a':
	cthresh    = atof(optarg);
	optind++;
	break;
      default:
	Usage(argv[0]);
	exit(-1);
      }
  }	
  if(optind != argc-2){
    Usage(argv[0]);
    exit(-1);
  }
  storeName = strdup(argv[optind++]);
  fileName  = strdup(argv[optind++]);
  fprintf(stderr,"storeName = %s\tfileName = %s\n",storeName,fileName);

  assert(existsFragStore(storeName) == TRUE);
  storeHandle = openFragStore(storeName,"rb");
  input = fopen (fileName,"r");
  assert(input != NULL);


  reader = (MesgReader)InputFileType_AS(input);
  while( reader(input,&pmesg) != EOF ) 
    {
      CGB_Type type;
      imesgtype = pmesg->t;

      switch(imesgtype){
      case MESG_IUM:
	{
	  Alignment_t *ali3;
	  int rc;
	  int rows;
	  int i;
	  sl_item it;
	  int fpf3=FALSE;
	  double pval3=1;
	  iunitig = (IntUnitigMesg*) pmesg->m;

	  if( iunitig->num_frags > 3 && (argc<4 ||iunitig->coverage_stat < cthresh)){
	    printf("\nInspecting Unitig %d\n",iunitig->iaccession);
	    printf("Number of frags = %d\n",iunitig->num_frags);	
	    printf("Length          = %d\n",iunitig->length);	
	    printf("Source          = %s\n",iunitig->source);	

	    AS_REZ_print_IUM_diffs(iunitig,storeHandle,printwhat);

	  }
	}
	break;
      default:
	{
	  //       Swallowing all other messages 
	}
      }
    }

  fclose(input);
  closeFragStore(storeHandle);
  free(storeName);
  free(fileName);
  return(0);
}
