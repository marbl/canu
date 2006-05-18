
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
 Module: AS_REZ
 Description: This file contains interface functions that convert
              assembler data structures into an Alignment_t which is 
	      data structure the microhet is working with.
	      It contains at the moment the following functions that
	      all return a pointer to a dynamically allocated
	      Alignment_t.
 Assumptions: 
**********************************************************************/

static char CM_ID[] = "$Id: MicroHetInterfaceREZ_test3.c,v 1.5 2006-05-18 18:30:31 vrainish Exp $";

#include <math.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"
#include "ScaffoldGraph_CGW.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "PrimitiveVA_MSG.h"
#include "GraphCGW_T.h"

#include "Array_CNS.h"
#include "MicroHetREZ_test3.h"


/* this function converts an array computed by consensus into
   an Alignment_t. The function that converts an IMP message 
   into the array is called IMP2Array. The array contains in the odd
  rows the sequence and in the even rows the quality values. */

#define USE_QUALITY
#define DEBUG 0

Alignment_t *AS_REZ_convert_array_to_alignment(char **ar, int c, int r){
  int i,j;

  //  fprintf(stderr,"Working on an array of size rXc %d X %d\n",r,c);
  Alignment_t *a = AS_REZ_allocate_alignment(c,r);
  for(i=0; i<c; i++)
    for(j=0; j<r; j++)
      {
	/*
	if(!( ar[2*j][i] == 'A' || ar[2*j][i] == 'C' || 
		ar[2*j][i] == 'G' || ar[2*j][i] == 'T' ||
		ar[2*j][i] == '-' || ar[2*j][i] == ' ' ||
		ar[2*j][i] == 'N'	))
	  fprintf(stderr,"Unexpected char %c at row %d col %d\n",
		  ar[2*j][i],j,i);
	*/
	assert( ar[2*j][i] == 'A' || ar[2*j][i] == 'C' || 
		ar[2*j][i] == 'G' || ar[2*j][i] == 'T' ||
		ar[2*j][i] == '-' || ar[2*j][i] == ' ' ||
		ar[2*j][i] == 'N'	);

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

static int32 *AS_REZ_fetch_int32(VA_TYPE(int32) *map, CDS_IID_t iid){
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

static uint32 *AS_REZ_fetch_uint32(VA_TYPE(uint32) *map, CDS_IID_t iid){
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

static uint64 *AS_REZ_fetch_uint64(VA_TYPE(uint64) *map, CDS_IID_t iid){
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

static uint64 *AS_REZ_get_uint64(VA_TYPE(uint64) *map, CDS_IID_t iid){
  uint64 *ret;
  ret = Getuint64(map,iid);
  if( ret == NULL )
    return NULL;
  else
    return ret;
}            

       
static uint32 *AS_REZ_get_uint32(VA_TYPE(uint32) *map, CDS_IID_t iid){
  uint32 *ret;
  ret = Getuint32(map,iid);
  if( ret == NULL )
    return NULL;
  else
    return ret;
}            

       



static void AS_REZ_get_info(CDS_IID_t iid,
                            FragStoreHandle frag_store,
                            tFragStorePartition *pfrag_store,
                            CDS_UID_t *locale,
                            uint32 *beg, uint32 *end,
                            FragType *type,
                            VA_TYPE(CDS_UID_t) *locales,
                            VA_TYPE(uint32) *fragtype,
                            VA_TYPE(uint32) *locbeg,
                            VA_TYPE(uint32) *locend,
                            ReadStructp input){
  if( iid != 0 ){ // blank positions have IID 0
    uint32 *mytype = AS_REZ_fetch_uint32(fragtype,iid);
    if(!mytype ){ // we have not seen that iid before 
      if ( frag_store != NULLFRAGSTOREHANDLE ) {
        if(getFragStore(frag_store,iid,FRAG_S_FIXED,input) != 0)
          assert(0);
      } else {
        if(getFragStorePartition(pfrag_store,iid,FRAG_S_FIXED,input) != 0)
          assert(0);
      }
      getReadType_ReadStruct(input,type); 
      getLocID_ReadStruct(input,locale);
      getLocalePos_ReadStruct(input,beg,end);
      
      // we store the locales for later lookup
      SetCDS_UID_t(locales,iid,locale);
      Setuint32(locbeg,iid,beg);
      Setuint32(locend,iid,end);
      Setuint32(fragtype,iid,(uint32 *)type);
    }
    else{
      // We know these entries have been initialized, so we can use them directly
      // without going through the fetch routines that behave special for uninitialized entries.
      *locale = *GetCDS_UID_t(locales,iid);
      *beg    = *Getuint32(locbeg,iid);	
      *end    = *Getuint32(locend,iid);	
      *type   = *Getuint32(fragtype,iid);
    }
  }
  else{
    *beg = *end = *locale = 0;
    *type = 'X';
  }
}



CGB_Type AS_REZ_get_simulator_type(IntUnitigMesg* ium_mesg){
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


// These statics accumulate the information for each
// fragment seen in this run of the Microhetrez tool
// they will grow to the length of the number of frags
// in the entire genome (not the number in this partition). SAK
static  VA_TYPE(CDS_UID_t)  *locales  = NULL;
static  VA_TYPE(uint32)  *locbeg   = NULL;
static  VA_TYPE(uint32)  *locend   = NULL;
static  VA_TYPE(uint32)  *fragtype = NULL;
static  ReadStructp input;

/* This function compresses shredded fragments from the same location 
   into basically a 1x coverage, such that there are no aritfical microhets;
   it also nulls out all but one position of a multibase gap, to reduce the
   effect of multibase indel polymorphisms on microhet detection. */
void compress_shreds_and_null_indels(int c, int r,  FragStoreHandle frag_store, 
                               tFragStorePartition *pfrag_store,
			       char **array, int **id_array, int verbose){
  int i,j,k;
  CDS_UID_t l1,l2;
  uint32 b1,b2;
  uint32 e1,e2;
  CDS_IID_t iid1,iid2;
  FragType t1,t2;

  if(locales == NULL){
    locales = CreateVA_CDS_UID_t(10000);
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
	  printf("loop (%d,%d) iid1=" F_IID "\n",i,j,iid1);
	if( iid1 != 0 ){ // blank positions have IID 0
	  AS_REZ_get_info(iid1,frag_store,pfrag_store,&l1,&b1,&e1,&t1,locales,fragtype,locbeg,locend,input);
	  if( verbose > 0 )
	    printf("(id,loc,beg,end,type) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                   iid1,l1,b1,e1,t1);
	  
	  if(AS_FA_SHREDDED(t1)){ 
	    for(k=j+1; k<r; k++){ // now eliminate the duplicates
	      if( verbose > 0 )
		printf("pos (%d,%d)\n",i,k);
	      iid2 = id_array[k][i];
	      AS_REZ_get_info(iid2,frag_store,pfrag_store,&l2,&b2,&e2,&t2,locales,fragtype,locbeg,locend,input);
	      if( verbose > 0)
		printf("(id2,loc2,beg2,end2,type2) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                       iid2,l2,b2,e2,t2);
	      
	      if( t1 == t2  && l2 == l1 ){
		uint32 min2,max1;
		if( verbose > 0 )
		  printf("changed position (%d,%d) from %c to blank, types (%c,%c), locales (" F_UID "," F_UID ")\n",i,k,array[i][k],t1,t2,l1,l2);
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
	    AS_REZ_get_info(iid1,frag_store,pfrag_store,&l1,&b1,&e1,&t1,locales,fragtype,locbeg,locend,input);
	    if( verbose > 0 )
	      printf("(id,loc,beg,end,type) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                     iid1,l1,b1,e1,t1);
	    
	    if(AS_FA_SHREDDED(t1)){ 
	      for(k=j+1; k<r; k++){ // now eliminate the duplicates
		if( verbose > 0 )
		  printf("pos (%d,%d)\n",i,k);
		iid2 = id_array[k][i];
		AS_REZ_get_info(iid2,frag_store,pfrag_store,&l2,&b2,&e2,&t2,locales,fragtype,locbeg,locend,input);
		if( verbose > 0)
		  printf("(id2,loc2,beg2,end2,type2) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                         iid2,l2,b2,e2,t2);
		
		if( t1 == t2  && l2 == l1 ){
		  uint32 min2,max1;
		  if( verbose > 0 )
		    printf("changed position (%d,%d) from %c to blank, types (%c,%c), locales (" F_UID "," F_UID ")\n",
                           i,k,array[i][k],t1,t2,l1,l2);
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




/* this function compresses shredded fragments from the same location into basically
   a 1x coverage, such that there are no aritfical microhets */
static void compress_shredded_frags(int c, int r,  FragStoreHandle frag_store, 
                               tFragStorePartition *pfrag_store,
			       char **array, int **id_array, int verbose){
  int i,j,k;
  CDS_UID_t l1,l2;
  uint32 b1,b2;
  uint32 e1,e2;
  CDS_IID_t iid1,iid2;
  FragType t1,t2;

  if(locales == NULL){
    locales = CreateVA_CDS_UID_t(10000);
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
    for(j=0; j<r-1; j++){
      iid1 = id_array[j][i];
      if( verbose > 0 )
	printf("loop (%d,%d) iid1=" F_IID "\n",i,j,iid1);
      if( iid1 != 0 ){ // blank positions have IID 0
	AS_REZ_get_info(iid1,frag_store,pfrag_store,&l1,&b1,&e1,&t1,locales,fragtype,locbeg,locend,input);
	if( verbose > 0 )
	  printf("(id,loc,beg,end,type) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                 iid1,l1,b1,e1,t1);
	
        if(AS_FA_SHREDDED(t1)){ 
	  for(k=j+1; k<r; k++){ // now eliminate the duplicates
	    if( verbose > 0 )
	      printf("pos (%d,%d)\n",i,k);
	    iid2 = id_array[k][i];
	    AS_REZ_get_info(iid2,frag_store,pfrag_store,&l2,&b2,&e2,&t2,locales,fragtype,locbeg,locend,input);
	    if( verbose > 0)
	      printf("(id2,loc2,beg2,end2,type2) = (" F_IID "," F_UID ",%d,%d,%c)\n",
                     iid2,l2,b2,e2,t2);
	    
	    if( t1 == t2  && l2 == l1 ){
	      uint32 min2,max1;
	      if( verbose > 0 )
		printf("changed position (%d,%d) from %c to blank, types (%c,%c), locales (" F_UID "," F_UID ")\n",
                       i,k,array[i][k],t1,t2,l1,l2);
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
  
  // cleanup (none, we made 'em static)
}


/* this function converts runs of gaps into a single gap followed by spaces; 
   the effect should be to reduce/eliminate the impact of multibase indels on
   the microhet detector, with slight loss of sensitivity */
static void MicroHet_discount_multibase_gaps(int c, int r,  char **array){

  int i,j;

  //  fprintf(stdout,"Testing for block gaps...\n");

  for(i=0; i<c-1; i++)
    for(j=0; j<r; j++)
      if(array[i][j]=='-'&&array[i+1][j]=='-'){
	//        fprintf(stdout,"Multibase gap!\n");
	array[i][j]=' ';
      }	
}


Alignment_t* AS_REZ_convert_IUM_to_alignment(IntUnitigMesg* ium,
                                             FragStoreHandle handle,
                                             tFragStorePartition *phandle,
					     int compress)
{
  int i;
  int rows;
  char **bqarray;
  int **idarray;
  int **oriarray;
  Alignment_t *ali;
#if DEBUG > 0
  printf("\nInspecting Unitig " F_IID "\n",ium->iaccession);
  printf("Number of frags = %d\n",ium->num_frags);	
  printf("Length          = " F_COORD "\n",ium->length);	
  printf("Source          = %s\n",ium->source);	
#endif 

  IMP2Array(ium->f_list,ium->num_frags,ium->length,
	    handle,phandle,NULLFRAGSTOREHANDLE,
            &rows,&bqarray,
            &idarray,&oriarray,
            0,READSTRUCT_LATEST);

  ali = AS_REZ_convert_array_to_alignment(bqarray,ium->length,rows);

  if(compress){
#define DEBUG_COMPRESS_AND_NULL 0
#if DEBUG_COMPRESS_AND_NULL  == 1
  printf("BEFORE COMPRESSION\n");
  AS_REZ_print_alignment(ali,90);
#endif
  compress_shreds_and_null_indels(ium->length,rows,handle,  phandle,
    		  ali->ali, idarray,DEBUG_COMPRESS_AND_NULL);
#if DEBUG_COMPRESS_AND_NULL == 1
  printf("AFTER COMPRESSION\n");
  AS_REZ_print_alignment(ali,90);
#endif
  }

  /* free the space that is allocated by IMP2Array */
  for(i=0; i<2*rows; i++)
    free(bqarray[i]);
  free(bqarray);

  for(i=0; i<rows; i++)
    free(idarray[i]);
  free(idarray);

  for(i=0; i<rows; i++)
    free(oriarray[i]);
  free(oriarray);

  return ali;
}


