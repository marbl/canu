
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

static char CM_ID[] = "$Id: MicroHetInterfaceREZ.c,v 1.1.1.1 2004-04-14 13:53:25 catmandew Exp $";

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
#include "MicroHetREZ.h"


/* this function converts an array computed by consensus into
   an Alignment_t. The function that converts an IMP message 
   into the array is called IMP2Array. The array contains in the odd
  rows the sequence and in the even rows the quality values. */

#define USE_QUALITY
#define DEBUG 0

Alignment_t *convert_array_to_alignment(char **ar, int c, int r){
  int i,j,k,count=0;

  Alignment_t *a = allocate_alignment(c,r);
  for(i=0; i<c; i++)
    for(j=0; j<r; j++)
      {
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

       


static void get_info(uint32 iid, FragStoreHandle frag_store, 
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
      Setuint32(fragtype,iid,type); // type is actually an "enum declared without a tag" which is signed.
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



CGB_Type get_simulator_type(IntUnitigMesg* ium_mesg){
  char *type;
  char *labels;
  char *end;
  int result;
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


/* this function compresses shredded fragments from the same location into basically
   a 1x coverage, such that there are no aritfical microhets */
static void compress_shredded_frags(int c, int r,  FragStoreHandle frag_store, 
			       char **array, int **id_array, int verbose){
  VA_TYPE(uint64)  *locales  = CreateVA_uint64(10000);
  VA_TYPE(uint32)  *locbeg   = CreateVA_uint32(10000);        
  VA_TYPE(uint32)  *locend   = CreateVA_uint32(10000);
  VA_TYPE(uint32)  *fragtype = CreateVA_uint32(10000);

  ReadStructp input = new_ReadStruct();  

  int i,j,k;
  uint64 l1,l2;
  uint32 b1,b2;
  uint32 e1,e2;
  uint32 iid1,iid2;
  FragType t1,t2;

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
	printf("loop (%d,%d) iid1=%d\n",i,j,iid1);
      if( iid1 != 0 ){ // blank positions have IID 0
	get_info(iid1,frag_store,&l1,&b1,&e1,&t1,locales,fragtype,locbeg,locend,input);
	if( verbose > 0 )
	  printf("(id,loc,beg,end,type) = (%d,%lu,%d,%d,%c)\n",iid1,l1,b1,e1,t1);
	
	if( t1 == AS_UBAC || t1 == AS_FBAC ){
	  for(k=j+1; k<r; k++){ // now eliminate the duplicates
	    if( verbose > 0 )
	      printf("pos (%d,%d)\n",i,k);
	    iid2 = id_array[k][i];
	    get_info(iid2,frag_store,&l2,&b2,&e2,&t2,locales,fragtype,locbeg,locend,input);
	    if( verbose > 0)
	      printf("(id2,loc2,beg2,end2,type2) = (%d,%lu,%d,%d,%c)\n",iid2,l2,b2,e2,t2);
	    
	    if( t1 == t2  && l2 == l1 ){
	      uint32 min2,max1;
	      if( verbose > 0 )
		printf("changed position (%d,%d) from %c to blank, types (%c,%c), locales (%lu,%lu)\n",i,k,array[i][k],t1,t2,l1,l2);
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
  
  // cleanup
  DeleteVA_uint64(locales);
  DeleteVA_uint32(locbeg);
  DeleteVA_uint32(locend);
  DeleteVA_uint32(fragtype);

  delete_ReadStruct(input);       
}




Alignment_t* convert_IUM_to_alignment(IntUnitigMesg* ium, FragStoreHandle handle)
{
  int i;
  int rows;
  UnitigStatus_t simple;
  char **bqarray;
  int **idarray;
  int **oriarray;
  Alignment_t *ali;
#if DEBUG > 0
  printf("\nInspecting Unitig %d\n",ium->iaccession);
  printf("Number of frags = %d\n",ium->num_frags);	
  printf("Length          = %d\n",ium->length);	
  printf("Source          = %s\n",ium->source);	
#endif 

  IMP2Array(ium->f_list,ium->num_frags,ium->length,
	    handle,NULLFRAGSTOREHANDLE,&rows,&bqarray,&idarray,&oriarray,0,READSTRUCT_LATEST);

  ali = convert_array_to_alignment(bqarray,ium->length,rows);
  
  //  print_alignment(ali,90);
  if( ium->iaccession == -1 ){
    printf("BEFORE COMPRESSION\n");
    print_alignment(ali,90);
    compress_shredded_frags(ium->length,rows,handle, 
			  ali->ali, idarray,1);
    printf("AFTER COMPRESSION\n");
    print_alignment(ali,90);
  }
  else
    compress_shredded_frags(ium->length,rows,handle, 
			    ali->ali, idarray,0);
  // print_alignment(ali,90);

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




/* WAIT for KARIN to write the routine CI2Array */

Alignment_t* convert_CI_to_alignment(ChunkInstanceT* ci, FragStoreHandle handle)
{
  int i;
  int rows;
  UnitigStatus_t simple;
  char **bqarray;
  char **oriarray;
  Alignment_t *ali;
  /* 
 //#if DEBUG > 0
  printf("\nInspecting Unitig %d\n",ium->iaccession);
  printf("Number of frags = %d\n",ium->num_frags);	
  printf("Length          = %d\n",ium->length);	
  //#endif

  IMP2Array(ium->f_list,ium->num_frags,ium->length,
	    handle,NULL,&rows,&bqarray,&oriarray,0,READSTRUCT_LATEST);

  ali = convert_array_to_alignment(bqarray,ium->length,rows);

  for(i=0; i<2*rows; i++)
    free(bqarray[i]);
  free(bqarray);

  */

  return ali;
}



