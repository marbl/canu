
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
static char CM_ID[] = "$Id: colCorr_CNS.c,v 1.14 2007-04-28 08:46:22 brianwalenz Exp $";

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"

#include "MicroHetREZ.h"
#include "Array_CNS.h"
#include "colCorr_CNS.h"

#define ALLOW_NULL_IN_BQARRAY		

static int col_contributing(Alignment_t *a, int c, int t)

{
  int observed=0;
 
  if( a->countA[c] >= t)
    observed++;
  if( a->countC[c] >= t)
    observed++;
  if( a->countG[c] >= t)
    observed++;
  if( a->countT[c] >= t)
    observed++;
  if( a->countDash[c] >= t)
    observed++;

  if( observed >= 2 )
    return TRUE;
  else
    return FALSE;
}



static void count_columns_copy(Alignment_t* a)
{
  register int i,j;
  for(i=0; i<a->cols; i++)
    {
      a->countA[i] = 0;
      a->countC[i] = 0;
      a->countG[i] = 0;
      a->countT[i] = 0;
      a->countDash[i]  = 0;
      a->countBlank[i] = 0;

      for(j=0; j<a->rows; j++)
	  switch(a->ali[i][j])
	    {
	    case 'A' :
	    case 'a' :
	      a->countA[i]++;
	      break;	
	    case 'C' :
	    case 'c' :
	      a->countC[i]++;
	      break;
	    case 'G' :
	    case 'g' :
	      a->countG[i]++;
	      break;
	    case 'T' :
	    case 't' :
	      a->countT[i]++;
	      break;
	    case '-' :
	      a->countDash[i]++;
	      break;
	    case ' ' :
	    case 'N' : // NOTE we count Ns as blanks
	    case 'n' : // NOTE we count Ns as blanks
	      a->countBlank[i]++;
	      break;
	    }	
    }
}



static Alignment_t* convert_MultiAlignT_to_alignment(MultiAlignT* inAlign, GateKeeperStore *handle, int ***idArray)
{
  int i;
  int rows;
  char **bqarray;
  int **idarray;
  int **oriarray;
  Alignment_t *ali;
  int length;
  int num_frags;
  int rc;

  length=GetMultiAlignLength(inAlign);
  num_frags=GetNumIntMultiPoss(inAlign->f_list);
#if DEBUG > 0
  printf("\nInspecting MultiAlignT %d\n",inAlign->id);
  printf("Number of frags = %d\n",inAlign->num_frags);	
  printf("Length          = %d\n",length);
  //  printf("Source          = %s\n",inAlign->source);	
#endif 

  rc = IMP2Array(GetIntMultiPos(inAlign->f_list,0),num_frags,length,
                 handle,&rows,&bqarray,&idarray,&oriarray,0,AS_READ_CLEAR_LATEST);

  // need idarray outside
  *idArray=idarray;


  ali = AS_REZ_convert_array_to_alignment(bqarray,length,rows);
  
  //  AS_REZ_print_alignment(ali,90);
  if( inAlign->id == -1 ){
    printf("BEFORE COMPRESSION\n");
    AS_REZ_print_alignment(ali,90);
    AS_REZ_compress_shreds_and_null_indels(length,rows,handle,
                                                ali->ali, idarray,1);
    printf("AFTER COMPRESSION\n");
    AS_REZ_print_alignment(ali,90);
  }
  else
    AS_REZ_compress_shreds_and_null_indels(length,rows,handle,
			    ali->ali, idarray,0);
  //AS_REZ_print_alignment(ali,90);

  /* free the space that is allocated by IMP2Array */
  for(i=0; i<2*rows; i++)
    safe_free(bqarray[i]);
  safe_free(bqarray);

  for(i=0; i<rows; i++)
    safe_free(oriarray[i]);
  safe_free(oriarray);

  count_columns_copy(ali);

  return ali;
}


static char 	max_count_char(Alignment_t *a, int col){
  char ret='A';
  int count;
  count=a->countA[col];
  if(a->countC[col]>count){
    count=a->countC[col];
    ret='C';
  }
  if(a->countG[col]>count){
    count=a->countG[col];
    ret='G';
  }
  if(a->countT[col]>count){
    count=a->countT[col];
    ret='T';
  }
  if(a->countDash[col]>count){
    count=a->countDash[col];
    ret='-';
  }
  return ret;
}



// test whether two "columns" in the alignment are consistent

static int is_consistent_partition(Alignment_t *a,int prevmm,int thismm, int rows, int **idArray, int ctgId){

  int j;
  int idprev,idthis,countprev,countthis;
  char majprev, majthis, charprev, charthis;
  int maxRows=100;
  int prevPart[100];
  int thisPart[100];
  char prevPartChar[100];
  char thisPartChar[100];
  char prevChar[100];
  char thisChar[100];
  char prevUsed[100];
  char thisUsed[100];
  char *partSpell,origSpell[3]={'.','0','1'};
  int n=0, i;
  int prevNonC=0,thisNonC=0;
  int ret;
  
  partSpell=&(origSpell[0])+1;

  // find the dominant character for each column
  // (or lexicographically smallest, in case of tie)
  majprev=max_count_char(a,prevmm);
  majthis=max_count_char(a,thismm);

  // over all values belonging to the alignment column
  for(j=0;j<rows;j++){

    // check the fragment ids at [prevmm][j] and [thismm][j] -- if not the
    // same, this means different fragment occupies the prev and this column in the
    // "row" in question.
    idprev=idArray[j][prevmm];
    idthis=idArray[j][thismm];
    if(idprev!=idthis||idprev==0){
      if(idprev==0){
	prevChar[j]=' ';
      } else {
	prevChar[j]=':';
      }
      if(idthis==0){
	thisChar[j]=' ';
      } else {
	thisChar[j]=':';
      }
      continue;
    }

    // figure out how many times the character at [prevmm][j] in prevmm
    charprev=a->ali[prevmm][j];
    switch(charprev){
    case 'A':
    case 'a':
      countprev=a->countA[prevmm];
      break;
    case 'C':
    case 'c':
      countprev=a->countC[prevmm];
      break;
    case 'G':
    case 'g':
      countprev=a->countG[prevmm];
      break;
    case 'T':
    case 't':
      countprev=a->countT[prevmm];
      break;
    case '-':
      countprev=a->countDash[prevmm];
      break;
    case ' ':
    case 'N':
    case 'n':
#ifdef ALLOW_NULL_IN_BQARRAY		
    case '\0':
#endif
      countprev=a->countBlank[prevmm];
      break;
    default:
      fprintf(stderr,"Trouble with unexpected character _%c_ (ascii %d)\n",a->ali[prevmm][j],(int)(a->ali[prevmm][j]));
      assert(0);
    }

    // figure out how many times the character at [thismm][j] in thismm
    charthis=a->ali[thismm][j];
    switch(charthis){
    case 'A':
    case 'a':
      countthis=a->countA[thismm];
      break;
    case 'C':
    case 'c':
      countthis=a->countC[thismm];
      break;
    case 'G':
    case 'g':
      countthis=a->countG[thismm];
      break;
    case 'T':
    case 't':
      countthis=a->countT[thismm];
      break;
    case '-':
      countthis=a->countDash[thismm];
      break;
    case ' ':
    case 'N':
    case 'n':
#ifdef ALLOW_NULL_IN_BQARRAY		
    case '\0':
#endif
      countthis=a->countBlank[thismm];
      break;
    default:
      fprintf(stderr,"Trouble with unexpected character _%c_ (ascii %d)\n",a->ali[thismm][j],(int)(a->ali[thismm][j]));
      assert(0);
    }
    
    // should find non-zero counts
    assert(countprev>0&&countthis>0);


    prevChar[j]=charprev;
    thisChar[j]=charthis;

    // exclude unconfirmed differences from the comparison of partitions
    if(countprev==1||countthis==1){
      continue;
    }

    prevNonC += 
      prevPart[n] = (a->ali[prevmm][j]==majprev) ? 0 : 1;

    thisNonC += 
      thisPart[n] = (a->ali[thismm][j]==majthis) ? 0 : 1;

    prevPartChar[n] = partSpell[prevPart[n]];
    thisPartChar[n] = partSpell[thisPart[n]];

    prevUsed[n] = a->ali[prevmm][j];
    thisUsed[n] = a->ali[thismm][j];

    n++;
    assert(n<maxRows-1);
  }


  if(n<=1||prevNonC==0||prevNonC==n||thisNonC==0||thisNonC==n){
    return -1;
  }




  ret=1;

  // check whether partition labels are all the same
  for(i=0;i<n;i++){
    if(prevPart[i]!=thisPart[i]){
      ret=0;
      break;
    }
  }

  if(!ret){

    //alternatively, completely opposite partition labels are also compatible
    ret=1;
    for(i=0;i<n;i++){
      if(prevPart[i]==thisPart[i]){
	ret=0;
	break;
      }
    }
  }

  return ret;
}

ColumnCorrelationT *test_correlated_columns(MultiAlignT* ma, 
					    GateKeeperStore *handle) {
  Alignment_t *ali;
  int i;
  int rows,length;
  int **idArray;
  static ColumnCorrelationT *colcorr=NULL;
  static int sizeColCorr=100;
  ColumnCorrelationT * ret;

  if(colcorr==NULL)
    colcorr=(ColumnCorrelationT*)safe_malloc(sizeof(ColumnCorrelationT)*sizeColCorr);
  
  ali = convert_MultiAlignT_to_alignment(ma,handle,&idArray);

  rows=ali->rows;
  length=GetMultiAlignLength(ma);

  if(rows < 4)
    ret = NULL;
  else {

      int i;
      int confirmCols=0;
      int prevmm=-1;

      //add_stop_to_colCorr(confirmCols,colcorr);
      colcorr[confirmCols].col  = -1;
      colcorr[confirmCols].corr = -1;
      
      for(i=0;i<length;i++){

	// if most freq. and 2nd most freq. counts both > 1, confirmed "SNP", so process:

	if(col_contributing(ali,i,2)){


	  if(prevmm>=0){
	    colcorr[confirmCols].corr=is_consistent_partition(ali,prevmm,i,rows,idArray,ma->id);
	  }else {
	    colcorr[confirmCols].corr=-1;
	  }

	  colcorr[confirmCols].col=i;

	  prevmm=i;
	  confirmCols++;

	  if(confirmCols>=sizeColCorr){
	    sizeColCorr*=2;
	    colcorr=(ColumnCorrelationT*)safe_realloc(colcorr,sizeof(ColumnCorrelationT)*sizeColCorr);
	  }

	}

      }
	
      //add_stop_to_colCorr(confirmCols,colcorr);
      colcorr[confirmCols].col  = -1;
      colcorr[confirmCols].corr = -1;

      ret=colcorr;
  }

  for(i=0; i<rows; i++)
    safe_free(idArray[i]);
  safe_free(idArray);

  AS_REZ_free_alignment(ali);

  return ret;
}
