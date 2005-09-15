
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
static char CM_ID[] = "$Id: MicroHetIUMREZ.c,v 1.5 2005-09-15 15:20:16 eliv Exp $";


#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h> /* man 3 getopt */
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


#define DEBUG -1



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
  int guessedUnknown = 0;
  // int tooShallow = 0;
  FILE* input;
  int fp1=0;
  int fn1=0;
  int cor1=0;
  int fp2=0;
  int fn2=0;
  int cor2=0;
  int fp3=0;
  int fn3=0;
  int cor3=0;
  int fp4=0;
  int fn4=0;
  int cor4=0;
  int fp5=0;
  int fn5=0;
  int cor5=0;
  int fp6=0;
  int fn6=0;
  int cor6=0;
  int fp7=0;
  int fn7=0;
  int cor7=0;
  int fp8=0;
  int fn8=0;
  int cor8=0;
  int tn1=0,tn2=0,tn3=0,tn4=0,tn5=0,tn6=0,tn7=0,tn8=0;
  int nt1=0,nt2=0,nt3=0,nt4=0,nt5=0,nt6=0,nt7=0,nt8=0;
  int rep_aln_printed=0;

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

      //      if(noOfIUMs>10)break;
      //      printf("Have handled %d unitigs<<<<<<<<<<<<<<<<<<,\n",noOfIUMs);

      switch(imesgtype){
      case MESG_IUM:
	{
	  Alignment_t *ali1,*ali2,*ali3,*ali4,*ali5,*ali6,*ali7,*ali8;
	  int fpf1=FALSE;
	  int fpf2=FALSE;
	  int fpf3=FALSE;
	  int fpf4=FALSE;
	  int fpf5=FALSE;
	  int fpf6=FALSE;
	  int fpf7=FALSE;
	  int fpf8=FALSE;
	  double pval1=1,pval2=1,pval3=1,pval4=1,pval5=1,pval6=1,pval7=1,pval8=1;
	  iunitig = (IntUnitigMesg*) pmesg->m;

#ifdef TRUSTCOVSTAT	  
	  if( iunitig->num_frags > 3 && iunitig->coverage_stat < cthresh){
#else
	  if( iunitig->num_frags > 3 ){
#endif
	    int simple1, simple2, simple3, simple4, simple5,simple6,
	      simple7,simple8;
	    int repetitiv = FALSE;

#if DEBUG > -1
	    printf("\nInspecting Unitig %d\n",iunitig->iaccession);
	    printf("Number of frags = %d\n",iunitig->num_frags);	
	    printf("Length          = %d\n",iunitig->length);	
	    printf("Source          = %s\n",iunitig->source);	
#endif

#if DEBUG > 0
	    printf("\nTEST1 (Knut segments )  ======= \n\n");
	    printf("\nTEST2 (Aaron segments)  ======= \n\n");
	    printf("\nTEST3 (Aaron whole ali) ======= \n\n");
	    printf("\nTEST4 (Knut window)     ======= \n\n");
	    printf("\nTEST5 (Aaron window)    ======= \n\n");
	    printf("\nTEST6 (Aaron Pairwise segments)    ======= \n\n");
	    printf("\nTEST7 (Aaron Pairwise whole ali)    ======= \n\n");
	    printf("\nTEST8 (Aaron Pairwise window)    ======= \n\n");
#endif

#ifdef RETURNPVALS
	    simple1 = is_IUM_simple(iunitig,storeHandle,&ali1,thresh,0,&pval1);
	    simple2 = is_IUM_MPsimple(iunitig,storeHandle,&ali2,thresh,0,&pval2);
	    simple3 = is_IUM_MPsimple(iunitig,storeHandle,&ali3,thresh,2,&pval3);
	    simple4 = is_IUM_simple(iunitig,storeHandle,&ali4,thresh,1,&pval4);
	    simple5 = is_IUM_MPsimple(iunitig,storeHandle,&ali5,thresh,1,&pval5);
	    simple6 = is_IUM_PWsimple(iunitig,storeHandle,&ali6,thresh,0,&pval6);
	    simple7 = is_IUM_PWsimple(iunitig,storeHandle,&ali7,thresh,2,&pval7);
	    simple8 = is_IUM_PWsimple(iunitig,storeHandle,&ali8,thresh,1,&pval8);
#else
	    simple1 = is_IUM_simple(iunitig,storeHandle,&ali1,thresh,0);
	    simple2 = is_IUM_MPsimple(iunitig,storeHandle,&ali2,thresh,0);
	    simple3 = is_IUM_MPsimple(iunitig,storeHandle,&ali3,thresh,2);
	    simple4 = is_IUM_simple(iunitig,storeHandle,&ali4,thresh,1);
	    simple5 = is_IUM_MPsimple(iunitig,storeHandle,&ali5,thresh,1);
	    simple6 = is_IUM_PWsimple(iunitig,storeHandle,&ali6,thresh,0);
	    simple7 = is_IUM_PWsimple(iunitig,storeHandle,&ali7,thresh,2);
	    simple8 = is_IUM_PWsimple(iunitig,storeHandle,&ali8,thresh,1);
#endif
	    type = get_simulator_type(iunitig);

	    if( type == RU_CGBTYPE || type == RR_CGBTYPE ){
#if DEBUG > 0
	      printf("Unitig is repetitive , A-stat=%f \n",iunitig->coverage_stat);
#endif
	      if(rep_aln_printed++<50){
		print_alignment(ali1,100);
	      }
	      repetitiv = TRUE;
	    }

#if DEBUG > 0
	    else
	      printf("Unitig is not repetitive , A-stat=%f \n",iunitig->coverage_stat);  
#endif

#if DEBUG > 0
	    printf("\nRESULTS ======= \n\n");
	    printf("Total Number of tested Unitigs = %d\n",noOfIUMs);

#else
	    printf(F_IID " %d %f %d %e %d %e %d %e %d %e %d %e %d %e %d %e %d %e\n",
                   iunitig->iaccession,
                   (int)type,
                   iunitig->coverage_stat,
                   simple1,pval1,
                   simple2,pval2,
                   simple3,pval3,
                   simple4,pval4,
                   simple5,pval5,
                   simple6,pval6,
                   simple7,pval7,
                   simple8,pval8);
#endif


	    fpf1 = fpf2 = fpf3 = fpf4 =fpf5 = fpf6 = fpf7 =fpf8 = FALSE;
	    switch(simple1){
	    case UNITIG_IS_SIMPLE:
#if DEBUG > 0
	      printf("Knut's test  : Assume unitig is simple\n");
#endif
	      //	      print_alignment(ali1,100);
	      if( repetitiv )
		fn1++;
	      else
		tn1++;
      
	      break;
	    case UNITIG_IS_UNKNOWN:
#if DEBUG > 0
	      printf("Knut's test  : Assume unitig is unknown (repetitive?)\n");
	      //	      print_alignment(ali1,100);
#endif
	      guessedUnknown++;
	      if( repetitiv ){
		fp1++;
		if( fp1 < 40 )
		  fpf1 = TRUE;
	      } else
		tn1++;
	      break;
	    case UNITIG_IS_SHALLOW:
	      nt1++;
#if DEBUG > 0
	      printf("Knut's test  : Unitig is too shallow\n");
#endif
	      // print_alignment(ali1,100);
	      break;
	    case UNITIG_IS_REPETITIVE:
#if DEBUG > 0
	      printf("Knut's test  : Assume unitig is repetitive\n");
#endif
	      //	      print_alignment(ali1,100);  
       	      if( ! repetitiv ){
		fp1++;
		if( fp1 < 40 )
		  fpf1 = TRUE;
	      }
	      else
		cor1++;
	      break;
	    }
	    
	    
	    switch(simple2){
	    case UNITIG_IS_SIMPLE:
#if DEBUG > 0
	      printf("Aaron's test : Assume unitig is simple\n");
#endif
	      //		print_alignment(ali2,100);
	      if( repetitiv )
		fn2++; 
	      else
		tn2++;

	      break;
	    case UNITIG_IS_UNKNOWN:
	      nt2++;
#if DEBUG > 0
	      printf("Aaron's test : Assume unitig is unknown\n");
#endif
	      //    print_alignment(ali2,100); 
	      break;
	    case UNITIG_IS_SHALLOW:
	      nt2++;
	      // print_alignment(ali2,100);
#if DEBUG > 0
	      printf("Aaron's test : Unitig is too shallow\n");
#endif
	      //	      tooShallow++;
	      break;
	    case UNITIG_IS_REPETITIVE:
#if DEBUG > 0
	      printf("Aaron's test : Assume unitig is repetitive\n");
#endif
	      //  print_alignment(ali2,100);
	      if( ! repetitiv ){
		fp2++;
		if( fp2 < 40 )
		  fpf2 = TRUE;	
	      }
	      else
		cor2++;
	      break;
	    }

	    switch(simple3){
	    case UNITIG_IS_SIMPLE:
#if DEBUG > 0
	      printf("Aaron's whole test : Assume unitig is simple\n");
#endif
	      //		print_alignment(ali3,100);
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
	      //    print_alignment(ali3,100); 
	      break;
	    case UNITIG_IS_SHALLOW:
	      nt3++;
	      // print_alignment(ali3,100);
#if DEBUG > 0
	      printf("Aaron's whole test : Unitig is too shallow\n");
#endif
	      //	      tooShallow++;
	      break;
	    case UNITIG_IS_REPETITIVE:
#if DEBUG > 0
	      printf("Aaron's whole test : Assume unitig is repetitive\n");
#endif
	      //  print_alignment(ali3,100);
	      if( ! repetitiv ){
		fp3++;
		if( fp3 < 40 )
		  fpf3 = TRUE;
	      }
	      else
		cor3++;
	      break;
	    }
	    
	    
	    switch(simple4){
	    case UNITIG_IS_SIMPLE:
#if DEBUG > 0
	      printf("Knuts's window : Assume unitig is simple\n");
#endif
	      //		print_alignment(ali4,100);
	      if( repetitiv ){
#if DEBUG > 0
		printf("Knut's whole test : FALSE NEGATIV\n");
#endif
		fn4++; 
	      }else
		tn4++;

	      break;
	    case UNITIG_IS_UNKNOWN:

#if DEBUG > 0
	      printf("Knut's window test : Assume unitig is unknown\n");
#endif
	      //    print_alignment(ali4,100); 
	      if( ! repetitiv )
		fp4++;
	      else
		cor4++;
	      break;
	    case UNITIG_IS_SHALLOW:
	      nt4++;
	      // print_alignment(ali4,100);
#if DEBUG > 0
	      printf("Knut's whole test : Unitig is too shallow\n");
#endif
	      //	      tooShallow++;
	      break;
	    case UNITIG_IS_REPETITIVE:
#if DEBUG > 0
	      printf("Knut's whole test : Assume unitig is repetitive\n");
#endif
	      //  print_alignment(ali4,100);
	      if( ! repetitiv ){
		fp4++;
		if( fp4 < 40 )
		  fpf4 = TRUE;
	      }
	      else
		cor4++;
	      break;
	    }
	    
	     
	    switch(simple5){
	    case UNITIG_IS_SIMPLE:
#if DEBUG > 0
	      printf("Aaron's window test : Assume unitig is simple\n");
#endif
	      //		print_alignment(ali5,100);
	      if( repetitiv ){
#if DEBUG > 0
		printf("Aaron's window test : FALSE NEGATIV\n");
#endif
		fn5++; 
	      }else
		tn5++;

	      break;
	    case UNITIG_IS_UNKNOWN:
	      nt5++;
#if DEBUG > 0
	      printf("Aaron's whole test : Assume unitig is unknown\n");
#endif
	      //    print_alignment(ali5,100); 
	      break;
	    case UNITIG_IS_SHALLOW:
	      nt5++;
	      // print_alignment(ali5,100);
#if DEBUG > 0
	      printf("Aaron's whole test : Unitig is too shallow\n");
#endif
	      //	      tooShallow++;
	      break;
	    case UNITIG_IS_REPETITIVE:
#if DEBUG > 0
	      printf("Aaron's whole test : Assume unitig is repetitive\n");
#endif
	      //  print_alignment(ali5,100);
	      if( ! repetitiv ){
		fp5++;
		if( fp5 < 40 )
		  fpf5 = TRUE;	
	      }
	      else
		cor5++;
	      break;
	    }



	    switch(simple6){
	    case UNITIG_IS_SIMPLE:
#if DEBUG > 0
	      printf("Aaron's Pairwise segment test : Assume unitig is simple\n");
#endif
	      //		print_alignment(ali6,100);
	      if( repetitiv ){
#if DEBUG > 0
		printf("Aaron's Pairwise segment test : FALSE NEGATIV\n");
#endif
		fn6++; 
	      }else
		tn6++;

	      break;
	    case UNITIG_IS_UNKNOWN:
	      nt6++;
#if DEBUG > 0
	      printf("Aaron's Pairwise segment test : Assume unitig is unknown\n");
#endif
	      //    print_alignment(ali6,100); 
	      break;
	    case UNITIG_IS_SHALLOW:
	      nt6++;
	      // print_alignment(ali6,100);
#if DEBUG > 0
	      printf("Aaron's Pairwise segment test : Unitig is too shallow\n");
#endif
	      //	      tooShallow++;
	      break;
	    case UNITIG_IS_REPETITIVE:
#if DEBUG > 0
	      printf("Aaron's Pairwise segment test : Assume unitig is repetitive\n");
#endif
	      //  print_alignment(ali6,100);
	      if( ! repetitiv ){
		fp6++;
		if( fp6 < 40 )
		  fpf6 = TRUE;	
	      }
	      else
		cor6++;
	      break;
	    }




	    switch(simple7){
	    case UNITIG_IS_SIMPLE:
#if DEBUG > 0
	      printf("Aaron's Pairwise whole test : Assume unitig is simple\n");
#endif
	      //		print_alignment(ali7,100);
	      if( repetitiv ){
#if DEBUG > 0
		printf("Aaron's Pairwise whole test : FALSE NEGATIV\n");
#endif
		fn7++; 
	      }else
		tn7++;

	      break;
	    case UNITIG_IS_UNKNOWN:
	      nt7++;
#if DEBUG > 0
	      printf("Aaron's Pairwise whole test : Assume unitig is unknown\n");
#endif
	      //    print_alignment(ali7,100); 
	      break;
	    case UNITIG_IS_SHALLOW:
	      nt7++;
	      // print_alignment(ali7,100);
#if DEBUG > 0
	      printf("Aaron's Pairwise whole test : Unitig is too shallow\n");
#endif
	      //	      tooShallow++;
	      break;
	    case UNITIG_IS_REPETITIVE:
#if DEBUG > 0
	      printf("Aaron's Pairwise whole test : Assume unitig is repetitive\n");
#endif
	      //  print_alignment(ali7,100);
	      if( ! repetitiv ){
		fp7++;
		if( fp7 < 40 )
		  fpf7 = TRUE;	
	      }
	      else
		cor7++;
	      break;
	    }




	    switch(simple8){
	    case UNITIG_IS_SIMPLE:
#if DEBUG > 0
	      printf("Aaron's Pairwise window test : Assume unitig is simple\n");
#endif
	      //		print_alignment(ali8,100);
	      if( repetitiv ){
#if DEBUG > 0
		printf("Aaron's Pairwise window test : FALSE NEGATIV\n");
#endif
		fn8++; 
	      }else
		tn8++;
	      break;
	    case UNITIG_IS_UNKNOWN:
	      nt8++;
#if DEBUG > 0
	      printf("Aaron's Pairwise window test : Assume unitig is unknown\n");
#endif
	      //    print_alignment(ali8,100); 
	      break;
	    case UNITIG_IS_SHALLOW:
	      nt8++;
	      // print_alignment(ali8,100);
#if DEBUG > 0
	      printf("Aaron's Pairwise window test : Unitig is too shallow\n");
#endif
	      //	      tooShallow++;
	      break;
	    case UNITIG_IS_REPETITIVE:
#if DEBUG > 0
	      printf("Aaron's Pairwise window test : Assume unitig is repetitive\n");
#endif
	      //  print_alignment(ali8,100);
	      if( ! repetitiv ){
		fp8++;
		if( fp8 < 40 )
		  fpf8 = TRUE;	
	      }
	      else
		cor8++;
	      break;
	    }
	    
	    if( fpf1 || fpf2 || fpf3 || fpf4 || fpf5 || fpf6 || fpf7 || fpf8)
	      print_alignment(ali1,100);


	
	    free_alignment(ali1);
	    free_alignment(ali2);
	    free_alignment(ali3);
	    free_alignment(ali4);
	    free_alignment(ali5);
	    free_alignment(ali6);
	    free_alignment(ali7);
	    free_alignment(ali8);

	    /*
	    if( simple == UNITIG_IS_REPETITIVE )
	      {
		Partition_t *p   = allocate_partition(ali->rows);
		Marker_t    *all = allocate_marker(ali->rows);
		bipartition(ali,all,p,0,ali->cols,0.001,1);
		print_part(p);
		free_marker(all);
		free_partition(p);
	      } 

	    */

#if DEBUG > 0
	    printf("Threshold = %f\n",thresh);
	    printf("Unitig is %s\n",(repetitiv ? " repetitiv\n" : " not repetitive\n"));
	    printf("Test1: %0.4d correct, %0.4d fp, %0.4d fn decisions\n",cor1,fp1,fn1);
	    printf("Test2: %0.4d correct, %0.4d fp, %0.4d fn decisions\n",cor2,fp2,fn2);
	    printf("Test3: %0.4d correct, %0.4d fp, %0.4d fn decisions\n",cor3,fp3,fn3);
	    printf("Test4: %0.4d correct, %0.4d fp, %0.4d fn decisions\n",cor4,fp4,fn4);
	    printf("Test5: %0.4d correct, %0.4d fp, %0.4d fn decisions\n",cor5,fp5,fn5);
	    printf("Test6: %0.4d correct, %0.4d fp, %0.4d fn decisions\n",cor6,fp6,fn6);
	    printf("Test7: %0.4d correct, %0.4d fp, %0.4d fn decisions\n",cor7,fp7,fn7);
	    printf("Test8: %0.4d correct, %0.4d fp, %0.4d fn decisions\n",cor8,fp8,fn8);
#endif
	    noOfIUMs++;
	  } else {
	    type = get_simulator_type(iunitig);
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

	    printf("Examined %d IUMs, could possibly test %d (i.e. no test on %d)\n\n",numunitig,noOfIUMs,numunitig-noOfIUMs);
	    printf("TEST1 (Knut segments )  ======= \n");
	    printf("TEST2 (Aaron segments)  ======= \n");
	    printf("TEST3 (Aaron whole ali) ======= \n");
	    printf("TEST4 (Knut window)     ======= \n");
	    printf("TEST5 (Aaron window)    ======= \n");
	    printf("TEST6 (Aaron Pairwise segments)    ======= \n");
	    printf("TEST7 (Aaron Pairwise whole ali)    ======= \n");
	    printf("TEST8 (Aaron Pairwise window)    ======= \n");
	    printf("Threshold = %f\n",thresh);
	    printf("\nTest #tp #fp #fn #tn #nt\n\n");
	    printf("1 %04d %04d %04d %04d %04d\n",cor1,fp1,fn1,tn1,nt1);
	    printf("2 %04d %04d %04d %04d %04d\n",cor2,fp2,fn2,tn2,nt2);
	    printf("3 %04d %04d %04d %04d %04d\n",cor3,fp3,fn3,tn3,nt3);
	    printf("4 %04d %04d %04d %04d %04d\n",cor4,fp4,fn4,tn4,nt4);
	    printf("5 %04d %04d %04d %04d %04d\n",cor5,fp5,fn5,tn5,nt5);
	    printf("6 %04d %04d %04d %04d %04d\n",cor6,fp6,fn6,tn6,nt6);
	    printf("7 %04d %04d %04d %04d %04d\n",cor7,fp7,fn7,tn7,nt7);
	    printf("8 %04d %04d %04d %04d %04d\n",cor8,fp8,fn8,tn8,nt8);


  fclose(input);
  closeFragStore(storeHandle);
  exit(0);
}
