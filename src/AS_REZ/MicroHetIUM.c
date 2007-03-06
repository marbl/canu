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

#include <stdio.h>
#include <assert.h>
#include <errno.h>

#include "AS_global.h"
#include "GraphCGW_T.h"  //  for CGB_Type
#include "PrimitiveVA_MSG.h"
#include "Array_CNS.h"
#include "MicroHetREZ.h"


void AS_REZ_print_informatives_alignment(Alignment_t *a,int nfrag)
{
  int i,j,l;
  int c = a->cols;
  int r = a->rows;
  int iter = 0;
  char *consensus,*printArray;
  int count[7]; // A C G T Dash N and total
  int *row2frag,lastfrag,infcols,currcol;
    
  consensus = (char*)safe_malloc (sizeof(char)*(a->cols));
  row2frag = (int*)safe_malloc(sizeof(int)*r);

  infcols=0;
  for(i=0; i<c; i++)
    {
      int submax,max,k;
      for(k=0; k<7; k++)
	count[k] = 0;
      for(j=0; j<r; j++)
	{
	  switch(a->ali[i][j])
	    {
              case 'A' :
                count[0]++;
                count[6]++;
                break;	
              case 'C' :
                count[1]++;
                count[6]++;
                break;
              case 'G' :
                count[2]++;
                count[6]++;
                break;
              case 'T' :
                count[3]++;
                count[6]++;
                break;
              case '-' :
                count[4]++;
                count[6]++;
                break;
	    }	
	}	

      max = 0;
      for(k=0; k<6; k++){
	if( count[max] < count[k] ){
	  max = k;
	}
      }
	
      submax = (max==0 ? 1 : 0);
      for(k=0;k<5;k++){
	if(k!=max)
	  if(count[submax] <= count[k])
	    submax=k;
      }
      if(count[submax]>1){
	switch( max )
	  {
            case 0 :
              consensus[i] = 'A';
              infcols++;
              break;	
            case 1 :
              consensus[i] = 'C';
              infcols++;
              break;
            case 2 :
              consensus[i] = 'G';
              infcols++;
              break;
            case 3 :
              consensus[i] = 'T';
              infcols++;
              break;
            case 4 :
              consensus[i] = '-';
              infcols++;
              break;
	  }
      } else {
	consensus[i] = '.';
      }

    }


  printArray = (char *) safe_malloc(sizeof(char)*(infcols+1)*nfrag);
  assert(printArray!=NULL);
  for(i=0;i<nfrag;i++){
    for(j=0;j<infcols;j++)
      (printArray+i*(infcols+1))[j]=' ';
    (printArray+i*(infcols+1))[j]='\0';
  }

  for(i=0;i<r;i++){
    row2frag[i] = -1;
  }
  lastfrag=-1;
  currcol=-1;
  for(j=0;j<c;j++){
    if(consensus[j]!='.')
      currcol++;
    for(i=0;i<r;i++){
      if(a->ali[j][i] != ' '){
	if(j==0||a->ali[j-1][i]==' '){
	  row2frag[i] = ++lastfrag;
	}
	if(consensus[j]!='.')
	  (printArray+row2frag[i]*(infcols+1))[currcol] =
	    (consensus[j] != a->ali[j][i] ? a->ali[j][i] : '.' );
      }
    }
  }

  assert(lastfrag+1==nfrag);
  for(i=0;i<nfrag;i++)
    printf("%s\t%d\n",printArray+i*(infcols+1),i);

  safe_free(consensus);
  safe_free(printArray);
  safe_free(row2frag);
}



void AS_REZ_print_informative_splits(Alignment_t *a,int nfrag)
{
  int i,j,l;
  int c = a->cols;
  int r = a->rows;
  int iter = 0;
  char *consensus,*printArray;
  int count[7]; // A C G T Dash N and total
  int *row2frag,lastfrag,infcols,currcol;
  int submax,max,k,subct,currct,clades;
    
  consensus = (char*)safe_malloc (sizeof(char)*(a->cols));
  row2frag = (int*)safe_malloc(sizeof(int)*r);

  infcols=0;
  for(i=0; i<c; i++)
    {
      for(k=0; k<7; k++)
	count[k] = 0;
      for(j=0; j<r; j++)
	{
	  switch(a->ali[i][j])
	    {
              case 'A' :
                count[0]++;
                count[6]++;
                break;	
              case 'C' :
                count[1]++;
                count[6]++;
                break;
              case 'G' :
                count[2]++;
                count[6]++;
                break;
              case 'T' :
                count[3]++;
                count[6]++;
                break;
              case '-' :
                count[4]++;
                count[6]++;
                break;
	    }	
	}	

      max = 0;
      for(k=0; k<6; k++){
	if( count[max] < count[k] ){
	  max = k;
	}
      }
	
      submax = (max==0 ? 1 : 0);
      for(k=0;k<5;k++){
	if(k!=max)
	  if(count[submax] <= count[k])
	    submax=k;
      }
      if(count[submax]>1){
	switch( max )
	  {
            case 0 :
              consensus[i] = 'A';
              infcols++;
              break;	
            case 1 :
              consensus[i] = 'C';
              infcols++;
              break;
            case 2 :
              consensus[i] = 'G';
              infcols++;
              break;
            case 3 :
              consensus[i] = 'T';
              infcols++;
              break;
            case 4 :
              consensus[i] = '-';
              infcols++;
              break;
	  }
      } else {
	consensus[i] = '.';
      }

    }


  printArray = (char *) safe_malloc(sizeof(char)*(infcols+1)*nfrag);
  for(i=0;i<nfrag;i++){
    for(j=0;j<infcols;j++)
      (printArray+i*(infcols+1))[j]=' ';
    (printArray+i*(infcols+1))[j]='\0';
  }

  for(i=0;i<r;i++){
    row2frag[i] = -1;
  }
  lastfrag=-1;
  currcol=-1;
  for(j=0;j<c;j++){
    // determine what fragment is in each row at this column
    for(i=0;i<r;i++)
      if(a->ali[j][i] != ' ')
	if(j==0||a->ali[j-1][i]==' ')
	  row2frag[i] = ++lastfrag;

    // if the column is informative ...
    if(consensus[j]!='.'){
      printf("(");

      clades=0; // number of non-empty character classes already seen

      // HANDLE FRAGMENTS WITH 'A' AT THIS COLUMN

      // count number of As
      currct=0;
      for(i=0;i<r;i++){
	if(a->ali[j][i]=='A')currct++;
      } 

      // if non-zero count, we will print for this state
      if(currct>0){

	// if not the first char state, need a separating comma
	if(clades>0)	  printf(",");

	// if count is 2 or more, we need an open bracket
	if(currct>1)	  printf("(");

	// find and print the fragments for this state
	subct=0;
	for(i=0;i<r;i++){
	  if(a->ali[j][i]=='A'){
	    // if not the first, separate from previous with comma
	    if(subct>0)	      printf(",");
	    // print the identifier
	    printf("F%d",row2frag[i]);
	    subct++;
	  }
	}
	assert(subct==currct);

	// if two or more, need close bracket
	if(currct>1) printf(")");

	// increment character states (classes) already output
	clades++;
      }


      // HANDLE FRAGMENTS WITH 'C' AT THIS COLUMN

      // count number of Cs
      currct=0;
      for(i=0;i<r;i++){
	if(a->ali[j][i]=='C')currct++;
      } 

      // if non-zero count, we will print for this state
      if(currct>0){

	// if not the first char state, need a separating comma
	if(clades>0)	  printf(",");

	// if count is 2 or more, we need an open bracket
	if(currct>1)	  printf("(");

	// find and print the fragments for this state
	subct=0;
	for(i=0;i<r;i++){
	  if(a->ali[j][i]=='C'){
	    // if not the first, separate from previous with comma
	    if(subct>0)	      printf(",");
	    // print the identifier
	    printf("F%d",row2frag[i]);
	    subct++;
	  }
	}
	assert(subct==currct);

	// if two or more, need close bracket
	if(currct>1) printf(")");

	// increment character states (classes) already output
	clades++;
      }


      // HANDLE FRAGMENTS WITH 'G' AT THIS COLUMN

      // count number of Cs
      currct=0;
      for(i=0;i<r;i++){
	if(a->ali[j][i]=='G')currct++;
      } 

      // if non-zero count, we will print for this state
      if(currct>0){

	// if not the first char state, need a separating comma
	if(clades>0)	  printf(",");

	// if count is 2 or more, we need an open bracket
	if(currct>1)	  printf("(");

	// find and print the fragments for this state
	subct=0;
	for(i=0;i<r;i++){
	  if(a->ali[j][i]=='G'){
	    // if not the first, separate from previous with comma
	    if(subct>0)	      printf(",");
	    // print the identifier
	    printf("F%d",row2frag[i]);
	    subct++;
	  }
	}
	assert(subct==currct);

	// if two or more, need close bracket
	if(currct>1) printf(")");

	// increment character states (classes) already output
	clades++;
      }


      // HANDLE FRAGMENTS WITH 'T' AT THIS COLUMN

      // count number of Cs
      currct=0;
      for(i=0;i<r;i++){
	if(a->ali[j][i]=='T')currct++;
      } 

      // if non-zero count, we will print for this state
      if(currct>0){

	// if not the first char state, need a separating comma
	if(clades>0)	  printf(",");

	// if count is 2 or more, we need an open bracket
	if(currct>1)	  printf("(");

	// find and print the fragments for this state
	subct=0;
	for(i=0;i<r;i++){
	  if(a->ali[j][i]=='T'){
	    // if not the first, separate from previous with comma
	    if(subct>0)	      printf(",");
	    // print the identifier
	    printf("F%d",row2frag[i]);
	    subct++;
	  }
	}
	assert(subct==currct);

	// if two or more, need close bracket
	if(currct>1) printf(")");

	// increment character states (classes) already output
	clades++;
      }


      // HANDLE FRAGMENTS WITH '-' AT THIS COLUMN

      // count number of Cs
      currct=0;
      for(i=0;i<r;i++){
	if(a->ali[j][i]=='-')currct++;
      } 

      // if non-zero count, we will print for this state
      if(currct>0){

	// if not the first char state, need a separating comma
	if(clades>0)	  printf(",");

	// if count is 2 or more, we need an open bracket
	if(currct>1)	  printf("(");

	// find and print the fragments for this state
	subct=0;
	for(i=0;i<r;i++){
	  if(a->ali[j][i]=='-'){
	    // if not the first, separate from previous with comma
	    if(subct>0)	      printf(",");
	    // print the identifier
	    printf("F%d",row2frag[i]);
	    subct++;
	  }
	}
	assert(subct==currct);

	// if two or more, need close bracket
	if(currct>1) printf(")");

	// increment character states (classes) already output
	clades++;
      }


      // CLOSE PARTITION
      printf(");\n");


    }
  }

  safe_free(consensus);
  safe_free(printArray);
  safe_free(row2frag);
}



////////////////////////////////////////////////////////////////////////////////




static int doPrintMPobs=0;

#define DEBUG -1

/* GETDECISION toggles whether to use old interface and return a
 *  decision, or to use new interface and return only a pvalue
 */
#undef GETDECISION 


CGB_Type get_simulator_type(IntUnitigMesg* ium_mesg){
  char *type;
  CGB_Type t;
  // See if this is a repeat, or we can pin it down to an interval
  type = strstr(ium_mesg->source,"gen> ");
  t = (unsigned int)XX_CGBTYPE;
  if(type){
    type += 5;
    if(!strncmp(type,"uu",2)){
      t = (unsigned int)UU_CGBTYPE;
    }else if(!strncmp(type,"ru",2)){
      t = (unsigned int)RU_CGBTYPE;
    }else if(!strncmp(type,"rr",2)){
      t = (unsigned int)RR_CGBTYPE;
    }else if(!strncmp(type,"ur",2)){
      t = (unsigned int)UR_CGBTYPE;
    }
  }
  return t;
}



/* this is the main test function for a unitig. 
   It returns a pvalue (roughly, a probability) that the unitig is simple
   (mismatches are random errors).
   A return value of 1.0 may indicate that the unitig was not deep enough for meaningful test.
*/
double AS_REZ_prob_IUM_MPsimple(IntUnitigMesg* ium, GateKeeperStore *handle)
{
  int i;
  int rows;
  double ret;

  char **bqarray;
  int **idarray;
  int **oriarray;
  IMP2Array(ium->f_list,
            ium->num_frags,
            ium->length,
	    handle,
            &rows,
            &bqarray,
            &idarray,
            &oriarray,
            0,
            AS_READ_CLEAR_LATEST);

  ret = AS_REZ_MP_MicroHet_prob(bqarray,idarray,handle, ium->length,rows);
  /* free the space that is allocated by IMP2Array */
  for(i=0; i<2*rows; i++)
    safe_free(bqarray[i]);
  safe_free(bqarray);

  for(i=0; i<rows; i++)
    safe_free(idarray[i]);
  safe_free(idarray);

  for(i=0; i<rows; i++)
    safe_free(oriarray[i]);
  safe_free(oriarray);

  return(ret);
}








Alignment_t* AS_REZ_convert_IUM_to_alignment(IntUnitigMesg* ium,
                                             GateKeeperStore *handle,
					     int compress)
{
  int i;
  int rows;
  char **bqarray;
  int **idarray;
  int **oriarray;
  Alignment_t *ali;

  IMP2Array(ium->f_list,
            ium->num_frags,
            ium->length,
	    handle,
            &rows,
            &bqarray,
            &idarray,
            &oriarray,
            0,
            AS_READ_CLEAR_LATEST);

  ali = AS_REZ_convert_array_to_alignment(bqarray,ium->length,rows);

  if(compress){
    //printf("BEFORE COMPRESSION\n");
    //AS_REZ_print_alignment(ali,90);

    AS_REZ_compress_shreds_and_null_indels(ium->length,
                                           rows,
                                           handle,
                                           ali->ali,
                                           idarray,
                                           0);  //  verbose

    //printf("AFTER COMPRESSION\n");
    //AS_REZ_print_alignment(ali,90);
  }

  /* free the space that is allocated by IMP2Array */
  for(i=0; i<2*rows; i++)
    safe_free(bqarray[i]);
  safe_free(bqarray);

  for(i=0; i<rows; i++)
    safe_free(idarray[i]);
  safe_free(idarray);

  for(i=0; i<rows; i++)
    safe_free(oriarray[i]);
  safe_free(oriarray);

  return ali;
}






/* this is the main test function for a unitig. 
   It returns 
   UNITIG_IS_SIMPLE     if the unitig is simple
   UNITIG_IS_UNKNOWN    if we cannot determine the status
   since there is no big enough segment
   or the critical value is below the
   minimum critical value.
   UNITIG_IS_REPETITIVE if the test could determine that
   the unitig is repetitive.
   UNITIG_IS_SHALLOW if the number of rows is smaller than 4
   The Alignment_t ali contains the alignment in which the list of segments
   can be inspected in order to find repetitive segments
*/
UnitigStatus_t AS_REZ_is_IUM_MPsimple(IntUnitigMesg* ium, GateKeeperStore *handle,
                                      Alignment_t **ali, double thresh, int variant, 
                                      double *pval)
{
  Marker_t *m;
  UnitigStatus_t ret = UNITIG_IS_SIMPLE;
  
  *pval=1.0; //Sets up default return for cases when no test can be made.

  *ali = AS_REZ_convert_IUM_to_alignment(ium,handle,TRUE);

  /* for this test we allocate a marker that is by default TRUE
     for all columns */

  if( (*ali)->rows < 4) {
    ret = UNITIG_IS_SHALLOW;
  } else {
    m = AS_REZ_allocate_marker((*ali)->rows);
    AS_REZ_count_columns(*ali,m);

    //doPrintMPobs=1;
    ret = AS_REZ_test_MPsimple(*ali,thresh,m,0,(*ali)->cols-1,pval);
    //doPrintMPobs=0;

    /* clean up dynamically allocated stuff */
    AS_REZ_free_marker(m);

  }

  return ret;
}








////////////////////////////////////////////////////////////////////////////////

#define PRINT_DOTS   0
#define PRINT_SPLITS 1

int
main(int argc, char **argv) {

  char   *storeName  = 0L;
  char   *fileName   = 0L;

  int     doWhat     = 0;
  int     printWhat  = PRINT_DOTS;
  double  thresh     = 0.0;
  double  cthresh    = 0.0;

  GateKeeperStore *storeHandle;
  FILE            *input;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      cthresh = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-t") == 0) {
      thresh = atof(argv[++arg]);
      doWhat = 1;
    } else if (strcmp(argv[arg], "-S") == 0) {
      printWhat = PRINT_SPLITS;
    } else if (strcmp(argv[arg], "-f") == 0) {
      storeName = argv[++arg];
    } else if (strcmp(argv[arg], "-i") == 0) {
      fileName = argv[++arg];
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }

    arg++;
  }

  if ((storeName == 0L) || (fileName == 0L) || (err)) {
    fprintf(stderr, "usage: %s [-a thresh] [-t thresh] [-S] [-f frgStore] [-i input]\n", argv[0]);
    exit(1);
  }

  storeHandle = openGateKeeperStore(storeName, FALSE);
  input       = fopen(fileName,"r");


  if (doWhat == 0) {
    GenericMesg  *pmesg;

    while( ReadProtoMesg_AS(input,&pmesg) != EOF ) {
      if (pmesg->t == MESG_IUM) {
        Alignment_t   *ali;  
        IntUnitigMesg *iunitig = (IntUnitigMesg*) pmesg->m;

        if ((iunitig->num_frags > 3) &&
            (iunitig->coverage_stat < cthresh)) {
          printf("\nInspecting Unitig %d\n",iunitig->iaccession);
          printf("Number of frags = %d\n",iunitig->num_frags);	
          printf("Length          = %d\n",iunitig->length);	
          printf("Source          = %s\n",iunitig->source);	

          ali = AS_REZ_convert_IUM_to_alignment(iunitig,storeHandle,FALSE);

          switch(printWhat){
            case PRINT_DOTS:
              AS_REZ_print_informatives_alignment(ali,iunitig->num_frags);
              break;
            case PRINT_SPLITS:
              AS_REZ_print_informative_splits(ali,iunitig->num_frags);
              break;
            default:
              fprintf(stderr,"AS_REZ_print_IUM_diffs doesn't know what you want to print\n");
              AS_REZ_free_alignment(ali);
              exit(-1);
          }

          AS_REZ_free_alignment(ali);
        }
      }
    }
  }








  if (doWhat == 1) {
    GenericMesg  *pmesg;

    int fp3=0;
    int fn3=0;
    int cor3=0;
    int tn3=0;
    int nt3=0;
#if DEBUG > -1
    int rep_aln_printed=0;
#endif
    int numunitig=0;


    while( ReadProtoMesg_AS(input,&pmesg) != EOF ) {
      CGB_Type type;

      if (pmesg->t == MESG_IUM) {
#ifdef GETDECISION
        Alignment_t *ali3;
        int fpf3=FALSE;
#endif
        double          pval3=1;
        IntUnitigMesg  *iunitig = (IntUnitigMesg*) pmesg->m;


        if ((iunitig->num_frags > 3)
#ifdef TRUSTCOVSTAT
            && (iunitig->coverage_stat < cthresh)
#endif
            ) {

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

          type = get_simulator_type(iunitig);

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
          pval3 = AS_REZ_prob_IUM_MPsimple(iunitig,storeHandle);
          printf(F_IID " %d %f %e\n",
                 iunitig->iaccession,(int)type,iunitig->coverage_stat,pval3);
#endif

          numunitig++;
        } else {
          type = get_simulator_type(iunitig);
          printf(F_IID " %d %f NOTEST\n",
                 iunitig->iaccession,type,iunitig->coverage_stat);
        }

        numunitig++;
      }
    }

    printf("Examined %d IUMs\n\n", numunitig);
    printf("TEST3 (Aaron whole ali) ======= \n");
    printf("Threshold = %f\n",thresh);
    printf("\nTest #tp #fp #fn #tn #nt\n\n");
    printf("3 %04d %04d %04d %04d %04d\n",cor3,fp3,fn3,tn3,nt3);
  }
}
