
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
static char CM_ID[] = "$Id: colCorr_CNS.c,v 1.7 2007-02-04 09:30:45 brianwalenz Exp $";

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA_MSG.h"

//#include "UtilsREZ.h"
#include "MicroHetREZ.h"

#include "colCorr_CNS.h"
#include "Array_CNS.h"

#undef DEBUG_CORR
//#define DEBUG_CORR 1
#define ALLOW_NULL_IN_BQARRAY		

/* this function actually computes whether a ncolumn in an alignment is xContributing 
 precondition is that the counts are set according to the desired marker */
static int col_contributing_copy(Alignment_t *a, int c, int t)

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



/* memory allocation and deallocation for Alignment_t */
static Alignment_t *allocate_alignment_copy (int c, int r)
{
  int i;
  Alignment_t *a = (Alignment_t*) safe_malloc(sizeof(Alignment_t));

  if( c > 0 && r > 0){
    a->cols = c;
    a->rows = r; 
    a->ali        = (char**) safe_calloc(sizeof(char*),c);
    a->countA     = (int*) safe_calloc(sizeof(int),c);
    a->countC     = (int*) safe_calloc(sizeof(int),c);
    a->countG     = (int*) safe_calloc(sizeof(int),c);
    a->countT     = (int*) safe_calloc(sizeof(int),c);
    a->countDash  = (int*) safe_calloc(sizeof(int),c);
    a->countBlank = (int*) safe_calloc(sizeof(int),c);
    for(i=0; i<c; i++)
      a->ali[i] = (char*) safe_calloc(sizeof(char),r);
    a->seqErrArray  = (double**) safe_calloc(sizeof(double*),c);
    for(i=0; i<c; i++)
      a->seqErrArray[i] = (double*) safe_calloc(sizeof(double),r);
    a->hasQuality = FALSE;
  }
  else
    {
      a->cols = 0;
      a->rows = 0;
      a->ali        = NULL;
      a->countA     = NULL;
      a->countC     = NULL;
      a->countG     = NULL;
      a->countT     = NULL;
      a->countDash  = NULL;
      a->countBlank = NULL;
      a->seqErrArray  = NULL;
    }
  return a;
}

static void free_alignment_copy(Alignment_t* a)
{
  int i;

  if(a->ali != NULL){
    for(i=0; i<a->cols; i++){
      assert(a->ali[i] != NULL);  
    free(a->ali[i]);
    }
    free(a->ali);

    for(i=0; i<a->cols; i++){
      assert(a->seqErrArray[i] != NULL);
      free(a->seqErrArray[i]);
    }
    free(a->seqErrArray);
    free(a->countA);
    free(a->countC);
    free(a->countG);
    free(a->countT);
    free(a->countDash);
    free(a->countBlank);
  }
  free(a);
}


/* a print function for alignments */

static void print_alignment_copy(Alignment_t *a,  int w)
{
  int i,j,l;
  int c = a->cols;
  int r = a->rows;
  int iter = 0;
  char *consensus;
  int count[6]; // A C G T Dash N
    
  consensus=(char*)safe_malloc(sizeof(char)*a->cols);
  assert(consensus!=NULL);

  for(i=0; i<a->cols; i++)
    {
      int max,k;
      for(k=0; k<6; k++)
	count[k] = 0;
      for(j=0; j<a->rows; j++)
	{
	  switch(a->ali[i][j])
	    {
	    case 'A' :
	      count[0]++;
	      break;	
	    case 'C' :
	      count[1]++;
	      break;
	    case 'G' :
	      count[2]++;
	      break;
	    case 'T' :
	      count[3]++;
	      break;
	    case '-' :
	      count[4]++;
	      break;
	    case 'N' :
	      count[5]++;
	      break;
	    }	
	}	
      max = 0;
      for(k=0; k<6; k++){
	if( count[max] < count[k] )
	  max = k;
      
	switch( max )
	  {
	  case 0 :
	    consensus[i] = 'A';
	    break;	
	  case 1 :
	    consensus[i] = 'C';
	    break;
	  case 2 :
	    consensus[i] = 'G';
	    break;
	  case 3 :
	    consensus[i] = 'T';
 	    break;
	  case 4 :
	    consensus[i] = '-';
	    break;
	  case 5 :
	    consensus[i] = 'N';
	    break;
	  }	
      }
    }
  do{
    int p = w*iter+w;
    l = ( p < c ? p : c );
    
    for(i=iter*w; i<l; i++)
      printf("%c",consensus[i]);
    printf(" %d\n",l);
    
    for(i=0; i<r; i++){
      for(j=iter*w; j<l; j++)
        if( a->ali[j][i] != consensus[j] )
          printf("%c",a->ali[j][i]);
        else
          printf("%c",'.');
      printf(" %d\n",l);
    }

    iter++;
    printf("\n");
  }while(iter*w < c);
  
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
	      a->countBlank[i]++;
	    case 'N' : // NOTE we count Ns as blanks
	    case 'n' : // NOTE we count Ns as blanks
	      a->countBlank[i]++;
	      break;
	    }	
    }
}


/* this function converts an array computed by consensus into
   an Alignment_t. The function that converts an IMP message 
   into the array is called IMP2Array. The array contains in the odd
  rows the sequence and in the even rows the quality values. */

#define USE_QUALITY
#define DEBUG 0

static Alignment_t *convert_array_to_alignment_copy(char **ar, int c, int r){
  int i,j;

  Alignment_t *a = allocate_alignment_copy(c,r);
  for(i=0; i<c; i++)
    for(j=0; j<r; j++)
      {
	/*
	  assert( ar[2*j][i] == 'A' || ar[2*j][i] == 'C' || 
		ar[2*j][i] == 'G' || ar[2*j][i] == 'T' ||
		ar[2*j][i] == '-' || ar[2*j][i] == ' ' ||
		ar[2*j][i] == 'N'
#ifdef ALLOW_NULL_IN_BQARRAY		
		|| ar[2*j][i] == '\0'
#endif
		);
	*/

	a->ali[i][j] = ar[2*j][i]; 
	
#ifdef USE_QUALITY
	a->hasQuality = TRUE;
	if( ar[2*j][i] != ' '){
	  a->seqErrArray[i][j] = 1.0/pow(10,((double)ar[2*j+1][i]-'0')/10.0);
	  //	  printf("Seq = %c, Quality = %c, Seq err %f\n",ar[2*j][i],ar[2*j+1][i],a->seqErrArray[i][j]);
	}
	else
	  a->seqErrArray[i][j] = 0.0;
#else
	a->seqErrArray[i][j] = 0.0;
#endif

      } 

  return a;
}


static int32 *fetch_int32(VA_TYPE(int32) *map, uint32 iid){
  int32 *ret;
  ret = Getint32(map,iid);
  if( ret == NULL )
    return NULL;
  else
    if( *ret == 0 )
      return NULL;
  else
    return ret;
}                   

static uint32 *fetch_uint32(VA_TYPE(uint32) *map, uint32 iid){
  uint32 *ret;
  ret = Getuint32(map,iid);
  if( ret == NULL )
    return NULL;
  else
    if( *ret == 0 )
      return NULL;
  else
    return ret;
} 

static uint64 *fetch_uint64(VA_TYPE(uint64) *map, uint32 iid){
  uint64 *ret;
  ret = Getuint64(map,iid);
  if( ret == NULL )
    return NULL;
  else
    if( *ret == 0 )
      return NULL;
  else
    return ret;
}          

static uint64 *get_uint64(VA_TYPE(uint64) *map, uint32 iid){
  uint64 *ret;
  ret = Getuint64(map,iid);
  if( ret == NULL )
    return NULL;
  else
    return ret;
}            

       
static uint32 *get_uint32(VA_TYPE(uint32) *map, uint32 iid){
  uint32 *ret;
  ret = Getuint32(map,iid);
  if( ret == NULL )
    return NULL;
  else
    return ret;
}            

       
static void my_get_info(uint32 iid, FragStoreHandle frag_store, 
			tFragStorePartition *pfrag_store,
		     uint64 *locale, uint32 *beg, uint32 *end, FragType *type,
		     VA_TYPE(uint64) *locales, VA_TYPE(uint32) *fragtype,
		     VA_TYPE(uint32) *locbeg, VA_TYPE(uint32) *locend,
		     ReadStructp input){
  if( iid != 0 ){ // blank positions have IID 0
    if(  fetch_uint64(locales,iid) == NULL ){ // we have not seen that iid before 
      if(getFragStore(frag_store,iid,FRAG_S_FIXED,input) != 0)
        assert(0);
      getReadType_ReadStruct(input,type); 
      getLocID_ReadStruct(input,locale);
      getLocalePos_ReadStruct(input,beg,end);
      
      // we store the locales for later lookup
      Setuint64(locales,iid,locale);
      Setuint32(locbeg,iid,beg);
      Setuint32(locend,iid,end);
      Setuint32(fragtype,iid,(uint32*)type); // type is actually an "enum declared without a tag" which is signed.
    }
    else{
      *locale = *fetch_uint64(locales,iid);
      *beg    = *get_uint32(locbeg,iid);	
      *end    = *get_uint32(locend,iid);	
      *type   = *fetch_uint32(fragtype,iid);
    }
  }
  else{
    *beg = *end = *locale = 0;
    *type = 'X';
  }
}




// These statics accumulate the information for each
// fragment seen in this run of the Microhetrez tool
// they will grow to the length of the number of frags
// in the entire genome (not the number in this partition). SAK

static  VA_TYPE(uint64)  *locales  = NULL;
static  VA_TYPE(uint32)  *locbeg   = NULL;
static  VA_TYPE(uint32)  *locend   = NULL;
static  VA_TYPE(uint32)  *fragtype = NULL;
static  ReadStructp input;

/* This function compresses shredded fragments from the same location 
   into basically a 1x coverage, such that there are no aritfical microhets;
   it also nulls out all but one position of a multibase gap, to reduce the
   effect of multibase indel polymorphisms on microhet detection. */
static void compress_shreds_and_null_indels_copy(int c, int r,  FragStoreHandle frag_store, 
                               tFragStorePartition *pfrag_store,
			       char **array, int **id_array, int verbose){
  int i,j,k;
  CDS_UID_t l1,l2;
  uint32 b1,b2;
  uint32 e1,e2;
  CDS_IID_t iid1,iid2;
  FragType t1,t2;

  if(locales == NULL){
    locales = CreateVA_uint64(10000);
    locbeg = CreateVA_uint32(10000);
    locend = CreateVA_uint32(10000);
    fragtype = CreateVA_uint32(10000);
    input = new_ReadStruct();
  }

  if(verbose > 0)
    printf("Alignment has %d columns and %d rows\n",c,r);

  if( verbose > 0 )
    for(i=0; i<c; i++)
      for(j=0; j<r; j++)
	printf("first loop (%d,%d) iid1=%d\n",i,j,id_array[j][i]);


  for(i=0; i<c; i++){
    if(i==c-1){

      for(j=0; j<r-1; j++){
	iid1 = id_array[j][i];
	if( verbose > 0 )
	  printf("loop (%d,%d) iid1=%d\n",i,j,iid1);
	if( iid1 != 0 ){ // blank positions have IID 0
	  my_get_info(iid1,frag_store,pfrag_store,&l1,&b1,&e1,&t1,locales,fragtype,locbeg,locend,input);
	  if( verbose > 0 )
	    printf("(id,loc,beg,end,type) = (" F_IID "," F_UID ",%u,%u,%c)\n",
                   iid1,l1,b1,e1,(char) t1);
	  
	  if(AS_FA_SHREDDED(t1)){ 
	    for(k=j+1; k<r; k++){ // now eliminate the duplicates
	      if( verbose > 0 )
		printf("pos (%d,%d)\n",i,k);
	      iid2 = id_array[k][i];
	      my_get_info(iid2,frag_store,pfrag_store,&l2,&b2,&e2,&t2,locales,fragtype,locbeg,locend,input);
	      if( verbose > 0)
		printf("(id2,loc2,beg2,end2,type2) = (" F_IID "," F_UID ",%u,%u,%c)\n",
                       iid2,l2,b2,e2,(char)t2);
	      
	      if( t1 == t2  && l2 == l1 ){
		uint32 min2,max1;
		if( verbose > 0 )
		  printf("changed position (%d,%d) from %c to blank, types (%c,%c), locales (" F_UID "," F_UID ")\n",
                         i,k,array[i][k],
                         (char) t1,(char) t2,l1,l2);
		min2 = ( b2 < e2 ? b2 : e2);
		max1 = ( b1 > e1 ? b1 : e1);
		
		if(max1 < min2)
		  ;//		 fprintf(stdout,"WARNING begin of frag2 is less than end of second !!\n");
		else{
		  id_array[k][i] = 0;
		  array[i][k] = ' ';
		}
	      }
	    }
	  }
	}
      }
    } else {
      for(j=0; j<r; j++){

	if(array[i][j]=='-'&&array[i+1][j]=='-'){
	  //	  fprintf(stdout,"Multibase gap!\n");
	  array[i][j]=' ';
	}	

	if(j<r-1){

	  iid1 = id_array[j][i];
	  if( verbose > 0 )
	    printf("loop (%d,%d) iid1=" F_IID "\n",i,j,iid1);
	  if( iid1 != 0 ){ // blank positions have IID 0
	    my_get_info(iid1,frag_store,pfrag_store,&l1,&b1,&e1,&t1,locales,fragtype,locbeg,locend,input);
	    if( verbose > 0 )
	      printf("(id,loc,beg,end,type) = (" F_IID "," F_UID ",%u,%u,%c)\n",
                     iid1,l1,b1,e1,(char) t1);
	    
	    if(AS_FA_SHREDDED(t1)){ 
	      for(k=j+1; k<r; k++){ // now eliminate the duplicates
		if( verbose > 0 )
		  printf("pos (%d,%d)\n",i,k);
		iid2 = id_array[k][i];
		my_get_info(iid2,frag_store,pfrag_store,&l2,&b2,&e2,&t2,locales,fragtype,locbeg,locend,input);
		if( verbose > 0)
		  printf("(id2,loc2,beg2,end2,type2) = (" F_IID "," F_UID ",%u,%u,%c)\n",
                         iid2,l2,b2,e2,(char) t2);
		
		if( t1 == t2  && l2 == l1 ){
		  uint32 min2,max1;
		  if( verbose > 0 )
		    printf("changed position (%d,%d) from %c to blank, types (%c,%c), locales (" F_UID "," F_UID ")\n",
                           i,k,array[i][k],(char) t1,(char) t2,l1,l2);
		  min2 = ( b2 < e2 ? b2 : e2);
		  max1 = ( b1 > e1 ? b1 : e1);
		  
		  if(max1 < min2)
		    ;//		 fprintf(stdout,"WARNING begin of frag2 is less than end of second !!\n");
		  else{
		    id_array[k][i] = 0;
		    array[i][k] = ' ';
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  // cleanup (none, we made 'em static)
}


static Alignment_t* convert_MultiAlignT_to_alignment_copy(MultiAlignT* inAlign, FragStoreHandle handle, tFragStorePartition *pfraghandle, int ***idArray)
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



  rc=IMP2Array(GetIntMultiPos(inAlign->f_list,0),num_frags,length,
	    handle,pfraghandle,NULLFRAGSTOREHANDLE,&rows,&bqarray,&idarray,&oriarray,0,READSTRUCT_LATEST);

  // need idarray outside
  *idArray=idarray;


  ali = convert_array_to_alignment_copy(bqarray,length,rows);
  
  //  print_alignment_copy(ali,90);
  if( inAlign->id == -1 ){
    printf("BEFORE COMPRESSION\n");
    print_alignment_copy(ali,90);
    compress_shreds_and_null_indels_copy(length,rows,handle,pfraghandle,
			  ali->ali, idarray,1);
    printf("AFTER COMPRESSION\n");
    print_alignment_copy(ali,90);
  }
  else
    compress_shreds_and_null_indels_copy(length,rows,handle,pfraghandle,
			    ali->ali, idarray,0);
  // print_alignment_copy(ali,90);

  /* free the space that is allocated by IMP2Array */
  for(i=0; i<2*rows; i++)
    free(bqarray[i]);
  free(bqarray);

  for(i=0; i<rows; i++)
    free(oriarray[i]);
  free(oriarray);

  count_columns_copy(ali);

  return ali;
}

static void add_stop_to_colCorr(int col, ColumnCorrelationT *colcorr){
  colcorr[col].col=-1;
  colcorr[col].corr=-1;
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

#ifdef DEBUG_CORR 
  prevChar[rows]='\0';
  thisChar[rows]='\0';
  prevUsed[n]='\0';
  thisUsed[n]='\0';
  prevPartChar[n]='\0';
  thisPartChar[n]='\0';

  fprintf(stderr,
	  "ctg %-8d col %-8d\n"
	  "\t%s\t%s\n"
	  "\t%s\n"
	  "\t%s\n"
	  "\t%s\t%s\n"
	  "ctg %-8d col %-8d\n"
	  "corr: %d\n",
	  ctgId, prevmm,
	  prevUsed,prevChar,
	  prevPartChar,
	  thisPartChar,
	  thisUsed,thisChar,
	  ctgId,thismm,
	  ret);
#endif

  return ret;
}

ColumnCorrelationT *test_correlated_columns(MultiAlignT* ma, 
					    FragStoreHandle handle,
					    tFragStorePartition *pfraghandle){
  Alignment_t *ali;
  int i;
  int rows,length;
  int **idArray;
  static ColumnCorrelationT *colcorr=NULL;
  static int sizeColCorr=100;
  ColumnCorrelationT * ret;

  if(colcorr==NULL){
    colcorr=(ColumnCorrelationT*)safe_malloc(sizeof(ColumnCorrelationT)*sizeColCorr);
    assert(colcorr!=NULL);
  }
  
  ali = convert_MultiAlignT_to_alignment_copy(ma,handle,pfraghandle,&idArray);

  rows=ali->rows;
  length=GetMultiAlignLength(ma);

  if(rows < 4)
    ret = NULL;
  else {

      int i;
      int confirmCols=0;
      int prevmm=-1;

      add_stop_to_colCorr(confirmCols,colcorr);
      
      for(i=0;i<length;i++){

	// if most freq. and 2nd most freq. counts both > 1, confirmed "SNP", so process:

	if(col_contributing_copy(ali,i,2)){


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
	    colcorr=(ColumnCorrelationT*)realloc(colcorr,sizeof(ColumnCorrelationT)*sizeColCorr);
	    assert(colcorr!=NULL);
	  }

	}

      }
	
      add_stop_to_colCorr(confirmCols,colcorr);
      ret=colcorr;
  }

  for(i=0; i<rows; i++)
    free(idArray[i]);
  free(idArray);

  free_alignment_copy(ali);

  return ret;
}
