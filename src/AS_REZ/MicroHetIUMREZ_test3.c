
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
static char CM_ID[] = "$Id: MicroHetIUMREZ_test3.c,v 1.5 2005-09-15 15:20:16 eliv Exp $";


#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* man 3 getopt */
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>

#include "AS_global.h"

#include "MicroHetREZ_test3.h"
#include "MicroHetScoreREZ_test3.h"
#include "MicroHetPartitionsREZ_test3.h"
#include "MicroHetInterfaceREZ_test3.h"

#include "UtilsREZ.h"
#include "AS_UTL_skiplist.h"
#include "AS_UTL_Var.h"

#include "PrimitiveVA_MSG.h"


#define DEBUG -1

/* GETDECISION toggles whether to use old interface and return a decision, or to use new interface and return only a pvalue  (#ifdef and #ifndef respectively)*/
#undef GETDECISION 
//#define GETDECISION 


int main (int argc, char *argv[]) {
  GenericMesg   *pmesg;  
  int numunitig=0;
  IntUnitigMesg *iunitig;
  MessageType    imesgtype;
  MesgReader     reader;
  FragStoreHandle storeHandle = 0;
  char           *storeName;
  char           *fileName;
  int            noOfIUMs=0;
  FILE* input;
  int fp3=0;
  int fn3=0;
  int cor3=0;
  int tn3=0;
  int nt3=0;
#if DEBUG > -1
  int rep_aln_printed=0;
#endif

  float thresh  = 0.001;
  float cthresh = 5.0;

  { /* Parse the argument list using "man 3 getopt". */
    storeName = strdup(argv[1]);
    fileName  = strdup(argv[2]);
    if( argc > 3)
      thresh    = atof(argv[3]);
    if( argc > 4) 
      cthresh   = atof(argv[4]); 
  }

  assert(existsFragStore(storeName) == TRUE);
  storeHandle = openFragStore(storeName,"rb");
  input = fopen (fileName,"r");
  assert(input != NULL);


  reader = (MesgReader)InputFileType_AS(input);
  while( reader(input,&pmesg) != EOF ) 
    {
      CGB_Type type;
      imesgtype = pmesg->t;

      //      if(noOfIUMs>100)break;
      //      printf("Have handled %d unitigs<<<<<<<<<<<<<<<<<<,\n",noOfIUMs);

      switch(imesgtype){
      case MESG_IUM:
	{
#ifdef GETDECISION
	  Alignment_t *ali3;
	  int fpf3=FALSE;
#endif
	  double pval3=1;
	  iunitig = (IntUnitigMesg*) pmesg->m;

#ifdef TRUSTCOVSTAT	  
	  if( iunitig->num_frags > 3 && iunitig->coverage_stat < cthresh){
#else
	  if( iunitig->num_frags > 3 ){
#endif
#ifdef GETDECISION
	    int simple3;
	    int repetitiv = FALSE;
#endif

#if DEBUG > -1
	    printf("\nInspecting Unitig " F_IID "\n",iunitig->iaccession);
	    printf("Number of frags = %d\n",iunitig->num_frags);	
	    printf("Length          = " F_COORD "\n",iunitig->length);	
	    printf("Source          = %s\n",iunitig->source);	
#endif

#if DEBUG > 0
	    printf("\nTEST3 (Aaron whole ali) ======= \n\n");
#endif

	    type = AS_REZ_get_simulator_type(iunitig);

#ifdef GETDECISION
	    simple3 = AS_REZ_is_IUM_MPsimple(iunitig,storeHandle,NULL,&ali3,thresh,2,&pval3);

	    if( type == RU_CGBTYPE || type == RR_CGBTYPE ){
#if DEBUG > 0
	      printf("Unitig is repetitive , A-stat=%f \n",
                     iunitig->coverage_stat);
#endif
#if DEBUG > -1
	      if(rep_aln_printed++<100){
		AS_REZ_print_alignment(ali3,100);
	      }
#endif
	      repetitiv = TRUE;
	    }

#if DEBUG > 0
	    else
	      printf("Unitig is not repetitive , A-stat=%f \n",
                     iunitig->coverage_stat);  
#endif

#if DEBUG > 0
	    printf("\nRESULTS ======= \n\n");
	    printf("Total Number of tested Unitigs = %d\n",noOfIUMs);

#else
	    printf(F_IID " %d %f %d %e\n",
                   iunitig->iaccession,
                   (int)type,
                   iunitig->coverage_stat,
                   simple3,pval3);
#endif


	    fpf3 = 0;

	    switch(simple3){
	    case UNITIG_IS_SIMPLE:
#if DEBUG > 0
	      printf("Aaron's whole test : Assume unitig is simple\n");
#endif
	      if( repetitiv ){
#if DEBUG > 0
		printf("Aaron's whole test : FALSE NEGATIV\n");
#endif
		fn3++; 
	      }else
		tn3++;

	      break;
	    case UNITIG_IS_UNKNOWN:
	      nt3++;
#if DEBUG > 0
	      printf("Aaron's whole test : Assume unitig is unknown\n");
#endif
	      break;
	    case UNITIG_IS_SHALLOW:
	      nt3++;
#if DEBUG > 0
	      printf("Aaron's whole test : Unitig is too shallow\n");
#endif
	      break;
	    case UNITIG_IS_REPETITIVE:
#if DEBUG > 0
	      printf("Aaron's whole test : Assume unitig is repetitive\n");
#endif
	      if( ! repetitiv ){
		fp3++;
		if( fp3 < 100 )
		  fpf3 = TRUE;
	      }
	      else
		cor3++;
	      break;
	    }
	    
	    

#if DEBUG > -1	    
	    if( fpf3 )
	      AS_REZ_print_alignment(ali3,100);
#endif


	    AS_REZ_free_alignment(ali3);

#if DEBUG > 0
	    printf("Threshold = %f\n",thresh);
	    printf("Unitig is %s\n",(repetitiv ? " repetitiv\n" : " not repetitive\n"));
	    printf("Test3: %0.4d correct, %0.4d fp, %0.4d fn decisions\n",cor3,fp3,fn3);
#endif

#else
	    pval3 = AS_REZ_prob_IUM_MPsimple(iunitig,storeHandle,NULL);
	    printf(F_IID " %d %f %e\n",
                   iunitig->iaccession,(int)type,iunitig->coverage_stat,pval3);
#endif

	    noOfIUMs++;
	  } else {
	    type = AS_REZ_get_simulator_type(iunitig);
	    printf(F_IID " %d %f NOTEST\n",
                   iunitig->iaccession,type,iunitig->coverage_stat);
	  }
        }

	  numunitig++;
	break;
      default:
	{
	  //       Swallowing all other messages 
	}
	}
   }

	    printf("Examined %d IUMs, could possibly test %d (i.e. no test on %d)\n\n",
                   numunitig,noOfIUMs,numunitig-noOfIUMs);
	    printf("TEST3 (Aaron whole ali) ======= \n");
	    printf("Threshold = %f\n",thresh);
	    printf("\nTest #tp #fp #fn #tn #nt\n\n");
	    printf("3 %04d %04d %04d %04d %04d\n",cor3,fp3,fn3,tn3,nt3);

  fclose(input);
  closeFragStore(storeHandle);
  return(0);
}
