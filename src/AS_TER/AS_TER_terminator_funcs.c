
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
 Module:      AS_TER
 Description: Functions used by the Terminator module to output 
              snapshot messages.
 Assumptions: There is no UID 0
**********************************************************************/

static char CM_ID[] = "$Id: AS_TER_terminator_funcs.c,v 1.3 2005-03-22 19:08:36 jason_miller Exp $";



#include "AS_TER_terminator.h"
#include "AS_TER_terminator_funcs.h"
#include "AS_TER_utils.h"
#include "AS_TER_alloc.h"

#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_distStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_version.h"
#include "PrimitiveVA.h"
#include "PrimitiveVA_MSG.h"

#include <SYS_UIDcommon.h>
#include <SYS_UIDclient.h>

#define DEBUG 0

/*--------------------------------------------------------------------*/
/*  Output Routine for all Messages in the genome snapshot. The static*/
/*  functions allocate a new message and make sure that the IIDs are  */
/*  correctly replaced by the UIDs                                    */
/*--------------------------------------------------------------------*/


/* allocate dynamic arrays for each message type that has his */
/* accession number replaced. These array are allocated and
   freed by output snapshot */


// artifical counter to keep track of assigned UIDs for scaffolds
// and output them to the map file

     //static ISFcount = 0;
static int32 simulator = FALSE;

static VA_TYPE(short)  *FRGpresent; // NEW
static VA_TYPE(short)  *BTGpresent; // NEW
static VA_TYPE(int)    *ICM1inISFmap;
static VA_TYPE(int)    *ICM2inISFmap;
static VA_TYPE(CDS_UID_t) *IUMmap;
static VA_TYPE(CDS_UID_t) *ICMmap;
static VA_TYPE(CDS_UID_t) *ISFmap;
static VA_TYPE(CDS_UID_t) *FRGmap;
static VA_TYPE(CDS_UID_t) *BTGmap; // NEW
static VA_TYPE(CDS_UID_t) *DSCmap; // NEW
static VA_TYPE(CDS_UID_t) *DSTmap;
static VA_TYPE(CDS_UID_t) *SCNmap;
static VA_TYPE(CDS_UID_t) *RPTmap;
static VA_TYPE(uint32) *ClearStartMap;
static VA_TYPE(uint32) *ClearEndMap;
static VA_TYPE(uint32) *ClearBactigStartMap;
static VA_TYPE(uint32) *ClearBactigEndMap;
static VA_TYPE(PtrT)   *SMAmap;
static VA_TYPE(uint32) *SMASizeMap;

static int PipeIn=FALSE;
static int PipeOut=FALSE;
static char *OutFileNameREZ = NULL;
static char *IID2UIDFileNameREZ = NULL;
static FILE *IFPREZ = NULL;
static FILE *OFPREZ = NULL;
static FILE *IID2UIDFPREZ = NULL;

/* the store handles are global such that we can open them
   in output_snapshot and use them in read_stores and the various fetch functions */
static FragStoreHandle FSHandle,BSHandle;
static GateKeeperStore GKPStore;

extern CDS_UID_t AS_TER_uidStart;

static void remove_output(void)
{
  char command[100];
  if( IFPREZ != NULL )
    fclose(IFPREZ);
  if( OFPREZ != NULL )
    fclose(OFPREZ);
  if( IID2UIDFPREZ != NULL )
    fclose(IID2UIDFPREZ);
  if( ! PipeOut ){ 
    sprintf(command,"rm -f %s ",OutFileNameREZ);
    system(command);
    sprintf(command,"rm -f %s ",IID2UIDFileNameREZ);
    system(command);
  }
}



void DumpIID2UIDmap(VA_TYPE(CDS_UID_t) *map,FILE *file){
  int i;
  for(i=0; i<GetNumCDS_UID_ts(map); i++){
      CDS_UID_t *di;
      di = GetCDS_UID_t(map,i);
      if( *di != 0 )
	fprintf(file,"%d " F_U64 "\n",i,*di);
    }
}


static short test_frg_present(CDS_IID_t iid){
  short* present;
  present = Getshort(FRGpresent,iid);
  if( present == NULL ){
    return FALSE;
  }
  else{
    if( *present == 0 )
      return FALSE;
    else
      return TRUE;
  }
}


static short test_btg_present(CDS_IID_t iid){
  short* present;
  present = Getshort(BTGpresent,iid);
  if( present == NULL ){
    return FALSE;
  }
  else{
    if( *present == 0 )
      return FALSE;
    else
      return TRUE;
  }
}


/* the function returns NULL if either the VA is not defined or
   holds only the default value 0 */
static uint32 *fetch_range(VA_TYPE(uint32) *map, CDS_IID_t iid){
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


/* the function returns NULL if the VA is not defined */
static uint32 *fetch_range_with_null(VA_TYPE(uint32) *map, CDS_IID_t iid){
  uint32 *ret;
  ret = Getuint32(map,iid);
  if( ret == NULL )
    return NULL;
  else
    return ret;
}


/* the function returns NULL if either the VA is not defined or
   holds only the default value 0 */

static CDS_UID_t *fetch_UID(VA_TYPE(CDS_UID_t) *map, CDS_IID_t iid){
  CDS_UID_t *ret;
  ret = GetCDS_UID_t(map,iid);
  if( ret == NULL )
    return NULL;
  else
    if( *ret == 0 )
      return NULL;
  else
    return ret;
}


/* this function first looks up the map to see whether we fetched it already from
   the fragStore. If not the values are retrieved from the store and stored in 
   the respective maps */
static uint32 *fetch_screenMatches_from_fragStore(VA_TYPE(uint32) *map, CDS_IID_t iid){
  uint32 *ret;
  assert(test_frg_present(iid));
  ret = fetch_range_with_null(map,iid);
  return ret;
}


static CDS_UID_t *fetch_UID_from_bactigStore(CDS_IID_t iid){
  if( test_btg_present(iid) ){
    CDS_UID_t* ret = fetch_UID(BTGmap,iid);
    assert( ret != NULL);
    return ret;
  }
  else{ 
    short truedummy=TRUE;
    ReadStructp input;
    assert( ! test_btg_present(iid) );

#if DEBUG > 1
    fprintf(stderr,"Setting BTGpresent for iid = %d\n",iid);
#endif
    Setshort(BTGpresent,iid,&truedummy);
    
    input = new_ReadStruct();
    if( 0 == getFragStore(BSHandle,iid,FRAG_S_FIXED,input) ){
      CDS_UID_t uid;
      CDS_UID_t *di;
      getAccID_ReadStruct(input,&uid);
      di = GetCDS_UID_t(BTGmap,iid);
      if( di != NULL )
	if( *di != 0 ){
	  char dummy[40];
	  sprintf(dummy,"Internal fragment ID %d occurred twice (BTGmap)",iid);
	  error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}
#if DEBUG >1
      fprintf(stderr,"Setting BTGmap for iid = %d\n",iid);
#endif
      SetCDS_UID_t(BTGmap,iid,&uid);
      delete_ReadStruct(input); 
    }
    else{
      char dummy[40];
      sprintf(dummy,"Internal fragment ID %d is not present in the bactig store",iid);
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
    input = new_ReadStruct();
    if( 0 == getFragStore(BSHandle,iid,FRAG_S_FIXED,input) ){
      uint32 cStart,cEnd;
      
      /* get the clear region from the frag store */
      getClearRegion_ReadStruct(input,&cStart,&cEnd,READSTRUCT_LATEST);
      Setuint32(ClearBactigStartMap,iid,&cStart);
      Setuint32(ClearBactigEndMap,iid,&cEnd);
    }
    delete_ReadStruct(input); 
    return GetCDS_UID_t(BTGmap,iid);   
  }
}




static CDS_UID_t *fetch_UID_from_fragStore(CDS_IID_t iid){
  if( test_frg_present(iid) ){
    CDS_UID_t* ret = fetch_UID(FRGmap,iid);
    assert( ret != NULL);
    return ret;
  }
  else{ 
    short truedummy=TRUE;
    ReadStructp input;
    assert( ! test_frg_present(iid) );
    Setshort(FRGpresent,iid,&truedummy);
    
    input = new_ReadStruct();
    if( 0 == getFragStore(FSHandle,iid,FRAG_S_FIXED,input) ){
      CDS_UID_t uid;
      CDS_UID_t *di;
      getAccID_ReadStruct(input,&uid);
      delete_ReadStruct(input);
      di = GetCDS_UID_t(FRGmap,iid);
      if( di != NULL )
	if( *di != 0 ){
	  char dummy[40];
	  sprintf(dummy,"Internal fragment ID %d occurred twice (FRGmap)",iid);
	  error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}
      SetCDS_UID_t(FRGmap,iid,&uid);
    }
    else{
      char dummy[40];
      sprintf(dummy,"Internal fragment ID %d is not present in the fragment store",iid);
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
    { 
      IntScreenMatch *matches;
      ReadStructp input =  new_ReadStruct();
      /* record the Screen Matches of each fragment */
      /* first get the number of screen matches and allocate the buffer */
      uint32 noOfMatches = getScreenMatches_ReadStruct(input,NULL,0);
      /* record the number of matches. 
       This is used in the IFG_to_AFG conversion
       to test whether the list is empty */
      Setuint32(SMASizeMap,iid,&noOfMatches);
      if(noOfMatches > 0){
	matches = (IntScreenMatch*) safe_calloc(noOfMatches,sizeof(IntScreenMatch));
	getScreenMatches_ReadStruct(input,matches,noOfMatches);
	SetPtrT(SMAmap,iid,(PtrT)&matches);
      }
      delete_ReadStruct(input);
    }
    {  
      ReadStructp input =  new_ReadStruct();
      if( 0 == getFragStore(FSHandle,iid,FRAG_S_FIXED,input) ){
	uint32 cStart,cEnd;
	
	/* get the clear region from the frag store */
	getClearRegion_ReadStruct(input,&cStart,&cEnd,READSTRUCT_LATEST);
	Setuint32(ClearStartMap,iid,&cStart);
	Setuint32(ClearEndMap,iid,&cEnd);
      }
      delete_ReadStruct(input); 
    }
    return GetCDS_UID_t(FRGmap,iid);   
  }
}
  

static uint32 *fetch_clearRange_from_fragStore(VA_TYPE(uint32) *map, CDS_IID_t iid){
  uint32* ret;
  assert( test_frg_present(iid) );
  ret = fetch_range_with_null(map,iid);
  assert( ret != NULL);
  return ret;
}


static uint32 *fetch_clearRange_from_bactigStore(VA_TYPE(uint32) *map, CDS_IID_t iid){
  uint32* ret;
  assert( test_btg_present(iid) );
  ret = fetch_range_with_null(map,iid);
  assert( ret != NULL);
  return ret;
}



static CDS_UID_t *fetch_UID_from_distStore(VA_TYPE(CDS_UID_t) *map, CDS_IID_t iid){
  CDS_UID_t* ret = fetch_UID(map,iid);
  if(ret == NULL)
  {
    GateKeeperDistanceRecord gkpd;
    if( 0 == getGateKeeperDistanceStore(GKPStore.dstStore,iid,&gkpd) ){
      CDS_UID_t *di;
      di = GetCDS_UID_t(map,iid);
      if( di != NULL )
	if( *di != 0 ){
	  char dummy[40];
	  sprintf(dummy,"Internal DST ID %d occurred twice",iid);
	  error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}
      SetCDS_UID_t(map,iid,&gkpd.UID);
      return GetCDS_UID_t(map,iid);
    }
    else{
      char dummy[40];
      sprintf(dummy,"Internal dist ID %d is not present in the dist store",iid);
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  }
  return ret;
}


static CDS_UID_t *fetch_UID_from_repeatStore(VA_TYPE(CDS_UID_t) *map, CDS_IID_t iid){
  CDS_UID_t* ret = fetch_UID(map,iid);
  // BEWARE. I return NULL in the case the repeat and screen library 
  // are not submitted (normally in the simulator)
  if( simulator )
    return NULL;
  if( ret != NULL )
    return ret;
  else{
    GateKeeperRepeatRecord gkpr;
    if( 0 == getGateKeeperRepeatStore(GKPStore.rptStore,iid,&gkpr) ){
      CDS_UID_t *di;
      di = GetCDS_UID_t(map,iid);
      if( di != NULL )
	if( *di != 0 ){
	  char dummy[40];
	  sprintf(dummy,"Internal RPT ID %d occurred twice",iid);
	  error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}
      SetCDS_UID_t(map,iid,&gkpr.UID);
      return GetCDS_UID_t(map,iid);
    }
    else{
      char dummy[40];
      // BEWARE. I return NULL in the case the repeat and screen library 
      // are not submitted (normally in the simulator)
      if( simulator )
	return NULL;

      sprintf(dummy,"Internal repeat ID %d is not present in the dist store",iid);
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  }
}



static CDS_UID_t *fetch_UID_from_screenStore(VA_TYPE(CDS_UID_t) *map, CDS_IID_t iid){
  CDS_UID_t* ret = fetch_UID(map,iid);
  // BEWARE. I return NULL in the case the repeat and screen library are not
  // submitted (normally in the simulator)
  if( simulator )
    return NULL;
  if( ret != NULL )
    return ret;
  else{
    GateKeeperScreenRecord gkps;
    if( 0 == getGateKeeperScreenStore(GKPStore.scnStore,iid,&gkps) ){
      CDS_UID_t *di;
      di = GetCDS_UID_t(map,iid);
      if( di != NULL )
	if( *di != 0 ){
	  char dummy[40];
	  sprintf(dummy,"Internal SCN ID %d occurred twice",iid);
	  error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}
      SetCDS_UID_t(map,iid,&gkps.UID);
      return GetCDS_UID_t(map,iid);
    }
    else{
      char dummy[40];
      // BEWARE. I return NULL in the case the repeat and screen library are not
      // submitted (normally in the simulator)
     if( simulator )
       return NULL;

      sprintf(dummy,"Internal SCN ID %d is not present in the dist store",iid);
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  }
}




/***********************/
/* conversion routines */
/***********************/


static SnapUnitigMesg* convert_IUM_to_UTG(IntUnitigMesg* iumMesg, int64 blockSize, int32 real)
     /*
       converts an IntUnitigMesg to a SnapUnitigMessage.
       What happend ?
       - A UID is assigned and recorded in IUMmap 
     */
{
  int i;
  CDS_UID_t uid;
  int32  uidStatus;
  CDS_UID_t interval_UID[4];

  CDS_UID_t *di;
  SnapUnitigMesg *utgMesg = 
    (SnapUnitigMesg*) safe_malloc(sizeof(SnapUnitigMesg));

#if DEBUG > 1
  fprintf(stderr,"IUM internal acc %u\n",iumMesg->iaccession);
#endif
  di = fetch_UID(IUMmap,iumMesg->iaccession);
  if( (di != NULL) && (*di != 0))
    {
      char dummy[40];
      sprintf(dummy,"Spotted IUM internal ID %d second time \n",
	      iumMesg->iaccession);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }

  /* This is a new message. We assign it a UID */
  uidStatus = get_next_uid(&uid,real);
  if( uidStatus != UID_CODE_OK )
    {
      uidStatus = get_uids(blockSize,interval_UID,real);
      get_next_uid(&uid,real);
    }	  
  if( UID_CODE_OK != uidStatus )
    { 
      char dummy[40];
      sprintf(dummy,"Could not get UID \n");
      remove_output();
      error(AS_TER_UIDSERVER_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  SetCDS_UID_t(IUMmap,iumMesg->iaccession,&uid);
  
  /* Set all toplevel fields */
  
  utgMesg->iaccession      = iumMesg->iaccession;
  utgMesg->eaccession      = uid;
#ifdef AS_ENABLE_SOURCE
  utgMesg->source         = strdup(iumMesg->source);
#endif
  utgMesg->coverage_stat  = iumMesg->coverage_stat;
  utgMesg->status         = iumMesg->status;
  utgMesg->a_branch_point = iumMesg->a_branch_point;
  utgMesg->b_branch_point = iumMesg->b_branch_point;
  utgMesg->length         = iumMesg->length;
  utgMesg->consensus      = strdup(iumMesg->consensus);
  utgMesg->quality        = strdup(iumMesg->quality);
  utgMesg->forced         = iumMesg->forced;
  utgMesg->num_frags      = iumMesg->num_frags;

  if( iumMesg->num_frags > 0 ){
    utgMesg->f_list = (SnapMultiPos*) safe_malloc(iumMesg->num_frags*sizeof(SnapMultiPos));

    for(i=0; i<iumMesg->num_frags; i++){
      utgMesg->f_list[i].type = iumMesg->f_list[i].type;
#ifdef AS_ENABLE_SOURCE
      utgMesg->f_list[i].source = strdup(iumMesg->f_list[i].source);
#endif
      if( iumMesg->f_list[i].type == AS_BACTIG )
	di = fetch_UID_from_bactigStore(iumMesg->f_list[i].ident);
      else
	di = fetch_UID_from_fragStore(iumMesg->f_list[i].ident);
      
      if( di == NULL ){
	char dummy[40];
	sprintf(dummy,"Reference before definition for fragment/bactig ID %d\n",
		iumMesg->f_list[i].ident);
	remove_output();
	error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
      }
      
      utgMesg->f_list[i].eident       = *di;
      utgMesg->f_list[i].delta_length = iumMesg->f_list[i].delta_length;
      utgMesg->f_list[i].position     = iumMesg->f_list[i].position;
      
      if( utgMesg->f_list[i].delta_length > 0 ){
	utgMesg->f_list[i].delta = iumMesg->f_list[i].delta;
      }
      else
	utgMesg->f_list[i].delta = NULL; 
    }
  } 
  return utgMesg;
}


/************************************************/


static SnapUnitigLinkMesg* convert_IUL_to_ULK(IntUnitigLinkMesg* iulMesg, int64 blockSize, int32 real)
     /*
       converts an IntUnitigLinkMesg to a SnapUnitigLinkMessage.
       What happend ?
       - The two internal Chunk IDs are replaced by UIDs
       - The jump_list is traversed and the internal fragment IDs in
       IntMatePairs are replaced by the external IDs
     */
{
  CDS_UID_t *di;
  int32  jumplistLength;
  int i;

  SnapUnitigLinkMesg *ulkMesg = (SnapUnitigLinkMesg*) safe_malloc(sizeof(SnapUnitigLinkMesg));
  
  /* Set all toplevel fields */
  /* look up the external IDs of the unitigs */
  di = GetCDS_UID_t(IUMmap,iulMesg->unitig1);
  if( di == NULL )
    {
      char dummy[40];
      sprintf(dummy,"IUL reference before definition error for unitig ID %d",iulMesg->unitig1);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  ulkMesg->eunitig1 = *di;
  di = GetCDS_UID_t(IUMmap,iulMesg->unitig2);
  if( di == NULL )
    {
      char dummy[40];
      sprintf(dummy,"IUL reference before definition error for unitig ID %d",iulMesg->unitig2);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  ulkMesg->eunitig2        = *di;

  ulkMesg->orientation    = iulMesg->orientation;
  ulkMesg->overlap_type   = iulMesg->overlap_type;
  ulkMesg->is_possible_chimera = iulMesg->is_possible_chimera;
  ulkMesg->includes_guide = iulMesg->includes_guide;
  ulkMesg->mean_distance  = iulMesg->mean_distance;
  ulkMesg->std_deviation  = iulMesg->std_deviation;
  ulkMesg->num_contributing = iulMesg->num_contributing;
  ulkMesg->status         = iulMesg->status;
    
  /* a case distinction to find out the number of elements in jump_list */
  if( iulMesg->overlap_type == AS_NO_OVERLAP )
    jumplistLength = ulkMesg->num_contributing;
  else
    jumplistLength = ulkMesg->num_contributing-1;

  ulkMesg->jump_list = (SnapMate_Pairs*) safe_malloc(jumplistLength*sizeof(SnapMate_Pairs));
  
  /* traverse the jump list and replace in the MatePairs the internal IDs
     by the external IDs */

  for(i=0; i<jumplistLength; i++)
    {
      /* NEW : nothing should change. 
	 Potentially we might have some links in between Bactigs? */
      
      di = fetch_UID_from_fragStore(iulMesg->jump_list[i].in1);
      if( di == NULL )
	{
	  char dummy[40];
	  sprintf(dummy,"Internal Fragment ID %d does not exist in Fragstore",iulMesg->jump_list[i].in1);
	  remove_output();
	  error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}  
      ulkMesg->jump_list[i].in1 = *di;
      di = fetch_UID_from_fragStore(iulMesg->jump_list[i].in2);
      if( di == NULL )
	{
	  char dummy[40];
	  sprintf(dummy,"Internal Fragment ID %d does not exist in Fragstore",iulMesg->jump_list[i].in2);
	  remove_output();
	  error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}  
      ulkMesg->jump_list[i].in2 = *di; 
      // the next line is due to SAULS change in proto I/O
      ulkMesg->jump_list[i].type = iulMesg->jump_list[i].type;
    }
  
  return ulkMesg;
}


/************************************************/


static SnapConConMesg* convert_ICM_to_CCO(IntConConMesg* icmMesg, int64 blockSize, int32 real)
     /* converts an IntConConMesg to a SnapConConMessage.
	What happend ?
	- A new UID is assigned to the ConConMesg
	- We initialize the ICMinISFmap to FALSE
     */
{
  int i;
  CDS_UID_t uid;
  int32  uidStatus;
  CDS_UID_t interval_UID[4];
  CDS_UID_t *di;
  int false = FALSE;

  SnapConConMesg *ccoMesg = (SnapConConMesg*) safe_malloc(sizeof(SnapConConMesg));
  
  /* Set all toplevel fields */
  /* This is a new message. Hence we get a new UID */
  uidStatus = get_next_uid(&uid,real);
  if( uidStatus != UID_CODE_OK )
    {
      uidStatus = get_uids(blockSize,interval_UID,real);
      get_next_uid(&uid,real);
    }	  
  if( UID_CODE_OK != uidStatus )
    { 
      char dummy[40];
      sprintf(dummy,"Could not get UID \n");
      remove_output();
      error(AS_TER_UIDSERVER_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }

  /* we assume that the numbers are in ascending order */
  di = GetCDS_UID_t(ICMmap,icmMesg->iaccession);
  if( (di != NULL) && (*di != 0)){
    char dummy[40];
    sprintf(dummy,"ICM internal contig number %d spotted second time \n",
	      icmMesg->iaccession);
    remove_output();
    error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }
  SetCDS_UID_t(ICMmap,icmMesg->iaccession,&uid);

  // record that we have used that icm so far in no scaffold
  Setint(ICM1inISFmap,icmMesg->iaccession,&false);
  Setint(ICM2inISFmap,icmMesg->iaccession,&false);

  ccoMesg->eaccession = uid;
  ccoMesg->iaccession = icmMesg->iaccession;
  ccoMesg->placed     = icmMesg->placed;
  ccoMesg->length     = icmMesg->length;
  ccoMesg->consensus  = strdup(icmMesg->consensus);
  ccoMesg->quality    = strdup(icmMesg->quality);
  ccoMesg->forced     = icmMesg->forced;
  ccoMesg->num_pieces = icmMesg->num_pieces;
  ccoMesg->num_unitigs= icmMesg->num_unitigs;

  if( ccoMesg->num_pieces > 0 ){ 
    ccoMesg->pieces = (SnapMultiPos*) safe_malloc(icmMesg->num_pieces*sizeof(SnapMultiPos));
    for(i=0; i<icmMesg->num_pieces; i++){// i loop
      ccoMesg->pieces[i].type = icmMesg->pieces[i].type;
#ifdef AS_ENABLE_SOURCE
      ccoMesg->pieces[i].source = strdup(icmMesg->pieces[i].source);
#endif

#if DEBUG > 1
      fprintf(stderr,"ICM = %d, no of pieces = %d, piece = %d, piece id = %d \n",icmMesg->iaccession,icmMesg->num_pieces,i,icmMesg->pieces[i].ident);
#endif

      if( icmMesg->pieces[i].type == AS_BACTIG ){
#if DEBUG > 1
	fprintf(stderr,"Fetching from bactig Store\n");
#endif
	di = fetch_UID_from_bactigStore(icmMesg->pieces[i].ident);
      }
      else{	
#if DEBUG > 1
	fprintf(stderr,"Fetching from frag Store\n");
#endif
	di = fetch_UID_from_fragStore(icmMesg->pieces[i].ident);
      }
      if( di == NULL ){
	char dummy[40];
	sprintf(dummy,"Reference before definition for fragment ID %d\n",
		icmMesg->pieces[i].ident);
	remove_output();
	error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
      }
      ccoMesg->pieces[i].eident       = *di;
      ccoMesg->pieces[i].delta_length = icmMesg->pieces[i].delta_length;
      ccoMesg->pieces[i].position     = icmMesg->pieces[i].position;
      
      if( ccoMesg->pieces[i].delta_length > 0 ){
	ccoMesg->pieces[i].delta = icmMesg->pieces[i].delta; /*** COPY BY REFERENCE ***/
      }
      else
	ccoMesg->pieces[i].delta = NULL;
      
    }  
  }
  
  if( ccoMesg->num_unitigs > 0 ){
    ccoMesg->unitigs = (UnitigPos*) safe_malloc(icmMesg->num_unitigs*sizeof(UnitigPos));
    for(i=0; i<icmMesg->num_unitigs; i++){
      ccoMesg->unitigs[i].type  = icmMesg->unitigs[i].type;
      di = GetCDS_UID_t(IUMmap,icmMesg->unitigs[i].ident);
      if( di == NULL ){
	char dummy[40];
	sprintf(dummy,"Reference before definition for unitig ID %d\n",
		icmMesg->pieces[i].ident);
	remove_output();
	error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
      }
      ccoMesg->unitigs[i].eident = *di;
      ccoMesg->unitigs[i].position = icmMesg->unitigs[i].position;
      ccoMesg->unitigs[i].delta = icmMesg->unitigs[i].delta; /*** COPY BY REFERENCE ***/
      ccoMesg->unitigs[i].delta_length = icmMesg->unitigs[i].delta_length;
    }
  }
  return ccoMesg;
}

/************************************************/

static SnapContigLinkMesg* convert_ICL_to_CLK(IntContigLinkMesg* iclMesg, int64 blockSize, int32 real)
     /* converts an IntContigLinkMesg to a SnapContigLinkMessage.
	What happend ? 
	- The two internal contig IDs are replaced by UIDs
	- The jump_list is traversed and the internal fragment IDs in
	IntMatePairs are replaced by the external IDs
     */
{
  CDS_UID_t *di;
  int32  jumplistLength;
  int i;

  SnapContigLinkMesg *clkMesg = (SnapContigLinkMesg*) safe_malloc(sizeof(SnapContigLinkMesg));
  
  /* Set all toplevel fields */
  /* replace the internal IDs by the external IDs */
  di = GetCDS_UID_t(ICMmap,iclMesg->contig1);
  if( di == NULL )
    {
      char dummy[40];
      sprintf(dummy,"ICL reference before definition error for contig ID %d",iclMesg->contig1);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  clkMesg->econtig1 = *di;
  di = GetCDS_UID_t(ICMmap,iclMesg->contig2);
  if( di == NULL )
    {
      char dummy[40];
      sprintf(dummy,"ICL reference before definition error for contig ID %d",iclMesg->contig2);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  clkMesg->econtig2 = *di;

  clkMesg->orientation    = iclMesg->orientation;
  clkMesg->overlap_type   = iclMesg->overlap_type;
  clkMesg->is_possible_chimera = iclMesg->is_possible_chimera;
  clkMesg->includes_guide = iclMesg->includes_guide;
  clkMesg->mean_distance  = iclMesg->mean_distance;
  clkMesg->std_deviation  = iclMesg->std_deviation;
  clkMesg->num_contributing = iclMesg->num_contributing;
  clkMesg->status         = iclMesg->status;


  /* make a case distinction to determine the length of the jump_list */
  if( iclMesg->overlap_type == AS_NO_OVERLAP )
    jumplistLength = clkMesg->num_contributing;
  else
    jumplistLength = clkMesg->num_contributing-1;

  if( jumplistLength > 0 )
    clkMesg->jump_list = (SnapMate_Pairs*) safe_malloc(jumplistLength*sizeof(SnapMate_Pairs));
  else
    clkMesg->jump_list = NULL;
 
  /* traverse the jump_list and and replace the internal fragment IDs
     by external fragment IDs */

  for(i=0; i<jumplistLength; i++){
    di = fetch_UID_from_fragStore(iclMesg->jump_list[i].in1);
    
    if( di == NULL ){
      char dummy[40];
      sprintf(dummy,"Internal Fragment ID %d does not exist in Fragstore",iclMesg->jump_list[i].in1);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }  
    clkMesg->jump_list[i].in1 = *di;
    
    di = fetch_UID_from_fragStore(iclMesg->jump_list[i].in2);
    if( di == NULL ){
      char dummy[40];
      sprintf(dummy,"Internal Fragment ID %d does not exist in Fragstore",iclMesg->jump_list[i].in2);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }  
    clkMesg->jump_list[i].in2 = *di; 
    // the next line is due to SAULS change in proto I/O
    clkMesg->jump_list[i].type = iclMesg->jump_list[i].type;
  }
  
  return clkMesg;
}



static SnapScaffoldLinkMesg* convert_ISL_to_SLK(InternalScaffoldLinkMesg* islMesg, int64 blockSize, int32 real)
     /* converts an InternalScaffoldLinkMesg to a SnapScaffoldLinkMessage.
	What happend ?
	- The two internal scaffold IDs are replaced by UIDs
	- The jump_list is traversed and the internal fragment IDs in
	IntMatePairs are replaced by the external IDs
     */
{
  CDS_UID_t *di;
  int32  jumplistLength;
  int i;

  SnapScaffoldLinkMesg *slkMesg = (SnapScaffoldLinkMesg*) safe_malloc(sizeof(SnapScaffoldLinkMesg));
  
  /* Set all toplevel fields */
  /* replace the internal IDs by the external IDs */
  di = GetCDS_UID_t(ISFmap,islMesg->iscaffold1);
  if( di == NULL )
    {
      char dummy[40];
      sprintf(dummy,"ISL reference before definition error for scaffold ID %d",islMesg->iscaffold1);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  slkMesg->escaffold1 = *di;
  di = GetCDS_UID_t(ISFmap,islMesg->iscaffold2);
  if( di == NULL )
    {
      char dummy[40];
      sprintf(dummy,"ISL reference before definition error for scaffold ID %d",islMesg->iscaffold2);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  slkMesg->escaffold2 = *di;


  slkMesg->orientation    = islMesg->orientation;
  slkMesg->includes_guide = islMesg->includes_guide;
  slkMesg->mean_distance  = islMesg->mean_distance;
  slkMesg->std_deviation  = islMesg->std_deviation;
  slkMesg->num_contributing = islMesg->num_contributing;

  jumplistLength = slkMesg->num_contributing;

  if( jumplistLength > 0 )
    slkMesg->jump_list = (SnapMate_Pairs*) safe_malloc(jumplistLength*sizeof(SnapMate_Pairs));
  else
    slkMesg->jump_list = NULL;
 
  /* traverse the jump_list and and replace the internal fragment IDs
     by external fragment IDs */

  for(i=0; i<jumplistLength; i++){
    di = fetch_UID_from_fragStore(islMesg->jump_list[i].in1);

    if( di == NULL ){
      char dummy[40];
      sprintf(dummy,"Internal Fragment ID %d does not exist in Fragstore",islMesg->jump_list[i].in1);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }  
    slkMesg->jump_list[i].in1 = *di;
    
    di = fetch_UID_from_fragStore(islMesg->jump_list[i].in2);
    if( di == NULL ){
      char dummy[40];
      sprintf(dummy,"Internal Fragment ID %d does not exist in Fragstore",islMesg->jump_list[i].in2);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }  
    slkMesg->jump_list[i].in2 = *di; 
    // the next line is due to SAULS change in proto I/O
    slkMesg->jump_list[i].type = islMesg->jump_list[i].type;
  }
  
  return slkMesg;
}




/************************************************/


static SnapScaffoldMesg* convert_ISF_to_SCF(IntScaffoldMesg* isfMesg, int64 blockSize, int32 real)
     /* converts an IntScaffoldMesg to a SnapScaffoldMessage.
	What happend ?
	- The contigs list is traversed and for each IntContigPair
	a SnapContigPair is created with the internal contig IDs replaced
	by external IDs
	- (07/01/00) in addition we store now the assigned UIDs in a VA
	in order to be able to dump them in the map file.
     */
{
  CDS_UID_t uid;
  int32  uidStatus;
  CDS_UID_t interval_UID[4];
  CDS_UID_t *di;
  int *used;
  int i;
  int true  = TRUE;

  SnapScaffoldMesg *scfMesg = (SnapScaffoldMesg*) safe_malloc(sizeof(SnapScaffoldMesg));

  /* NEW : the ISF messages HAVE now an internal ID
     so we can now check for the IID to UID mapping */

  di = GetCDS_UID_t(ISFmap,isfMesg->iaccession);
  if( di != NULL )
    {
      char dummy[40];
      sprintf(dummy,"ISF internal id %d spotted second time \n",
	      isfMesg->iaccession);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  uidStatus = get_next_uid(&uid,real);
  if( uidStatus != UID_CODE_OK )
    {
      uidStatus = get_uids(blockSize,interval_UID,real);
      get_next_uid(&uid,real);
    }	  
  if( UID_CODE_OK != uidStatus )
    { 
      char dummy[40];
      sprintf(dummy,"Could not get UID \n");
      remove_output();
      error(AS_TER_UIDSERVER_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  SetCDS_UID_t(ISFmap,isfMesg->iaccession,&uid);


  scfMesg->iaccession = isfMesg->iaccession;
  scfMesg->eaccession = uid;

  scfMesg->num_contig_pairs = isfMesg->num_contig_pairs;
  scfMesg->contig_pairs = NULL;

  /* check whether there are any contig pairs */
  if( scfMesg->num_contig_pairs > 0 )
    {
      scfMesg->contig_pairs = (SnapContigPairs*) safe_malloc(scfMesg->num_contig_pairs*sizeof(SnapContigPairs));
  

      for(i=0; i<isfMesg->num_contig_pairs; i++)
	{
	  int32 con1 = isfMesg->contig_pairs[i].contig1;
	  int32 con2 = isfMesg->contig_pairs[i].contig2;
	  
	  // try to find out whether we used the con1 before
	  used = Getint(ICM1inISFmap,con1);
	  if( used == NULL )
	    {
	      char dummy[40];
	      sprintf(dummy,"ICM1inISFmap is not initialized for contig ID %d",con1);
	      remove_output();
	      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  if( *used == TRUE )
	    {
	      char dummy[40];
	      sprintf(dummy,"Contig ID %d was already used as ct1 in a scaffold",con1);
	      remove_output();
	      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  else
	    Setint(ICM1inISFmap,con1,&true);

	  // try to find out whether we used the con2 before
	  used = Getint(ICM2inISFmap,con2);
	  if( used == NULL )
	    {
	      char dummy[40];
	      sprintf(dummy,"ICM2inISFmap is not initialized for contig ID %d",con2);
	      remove_output();
	      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  if( *used == TRUE )
	    {
	      char dummy[40];
	      sprintf(dummy,"Contig ID %d was already used as ct2 in a scaffold",con2);
	      remove_output();
	      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  else
	    Setint(ICM2inISFmap,con2,&true);



	  di = GetCDS_UID_t(ICMmap,con1);
	  if( di == NULL )
	    {
	      char dummy[40];
	      sprintf(dummy,"ISF reference before definition for contig ID %d",con1);
	      remove_output();
	      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  scfMesg->contig_pairs[i].econtig1 = *di;
	  
	  di = GetCDS_UID_t(ICMmap,con2);
	  if( di == NULL )
	    {
	      char dummy[40];
	      sprintf(dummy,"ISF reference before definition for contig ID %d",con2);
	      remove_output();
	      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  scfMesg->contig_pairs[i].econtig2 = *di; 
      
	  scfMesg->contig_pairs[i].mean   = isfMesg->contig_pairs[i].mean;
	  scfMesg->contig_pairs[i].stddev = isfMesg->contig_pairs[i].stddev;
	  scfMesg->contig_pairs[i].orient = isfMesg->contig_pairs[i].orient;
	}
    }
  else // special case if there are zero contig pairs there is ONE
       // with the second contig id -1
    {
      scfMesg->contig_pairs = (SnapContigPairs*) safe_malloc(sizeof(SnapContigPairs));
      {
	int32 con1 = isfMesg->contig_pairs[0].contig1;
	// try to find out whether we used the con1 before
	used = Getint(ICM1inISFmap,con1);
	if( used == NULL )
	  {
	    char dummy[40];
	    sprintf(dummy,"ICM1inISFmap is not initialized for contig ID %d",con1);
	    remove_output();
	    error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	  }  
	if( *used == TRUE )
	  {
	    char dummy[40];
	    sprintf(dummy,"Contig ID %d was already used as ct1 in a scaffold",con1);
	    remove_output();
	    error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	  }  
	else
	  Setint(ICM1inISFmap,con1,&true);	

	di = GetCDS_UID_t(ICMmap,con1);
	if( di == NULL )
	  {
	    char dummy[40];
	    sprintf(dummy,"ISF reference before definition for contig ID %d",con1);
	    remove_output();
	    error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	  }  	

	scfMesg->contig_pairs[0].econtig1 = *di;
      	scfMesg->contig_pairs[0].econtig2 = *di;

	scfMesg->contig_pairs[0].mean   = isfMesg->contig_pairs[0].mean;
	scfMesg->contig_pairs[0].stddev = isfMesg->contig_pairs[0].stddev;
	scfMesg->contig_pairs[0].orient = isfMesg->contig_pairs[0].orient;
      }
    }

  return scfMesg;
}


/************************************************/


static SnapMateDistMesg* convert_IMD_to_MDI(IntMateDistMesg* imdMesg, int64 blockSize, int32 real)
     /* converts an IntMateDistMesg to a SnapMateDistMessage.
	What happend ?
	- The intDistance ID is replaced by the corresponding UIDs
     */
{
  CDS_UID_t *di;

  SnapMateDistMesg *mdiMesg = (SnapMateDistMesg*) safe_malloc(sizeof(SnapMateDistMesg));
  
  /* Set all toplevel fields */
  di = fetch_UID_from_distStore(DSTmap,imdMesg->refines);
  if( di == NULL )
    {
      char dummy[40];
      sprintf(dummy,"IMD reference before definition error for ID %d",imdMesg->refines);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  mdiMesg->erefines = *di;
  mdiMesg->irefines = imdMesg->refines;
  mdiMesg->min     = imdMesg->min;
  mdiMesg->max     = imdMesg->max;
  mdiMesg->mean    = imdMesg->mean;
  mdiMesg->stddev  = imdMesg->stddev;
  mdiMesg->num_buckets = imdMesg->num_buckets;
  mdiMesg->histogram = NULL;

  /* the histogram s not changed. We can memcpy it */  
  if(  mdiMesg->num_buckets > 0 ){
    mdiMesg->histogram = (int32*) safe_malloc(mdiMesg->num_buckets*sizeof(int32));
    memcpy(mdiMesg->histogram,imdMesg->histogram,mdiMesg->num_buckets*sizeof(int32));
  }
  return mdiMesg;
}


/************************************************/


static SnapDegenerateScaffoldMesg* convert_IDS_to_DSC(IntDegenerateScaffoldMesg* idsMesg, int64 blockSize, int32 real)
     /* converts an  IntDegenerateScaffoldMesg to a  SnapDegenerateScaffoldMesg
	What happend ?
	- The icontig is replaced by the corresponding UID
	- the degenerate scaffold is assigned a new UID
     */
{
  CDS_UID_t uid;
  int32  uidStatus;
  CDS_UID_t interval_UID[4];
  CDS_UID_t *di;

  SnapDegenerateScaffoldMesg* dscMesg = (SnapDegenerateScaffoldMesg*) safe_malloc(sizeof(SnapDegenerateScaffoldMesg));
  
  /* check whether the internal id is defined */
  di = GetCDS_UID_t(ICMmap,idsMesg->icontig);
  if( di == NULL )
    {
      char dummy[40];
      sprintf(dummy,"Could not find Contig iid %d in ICMmap \n",
	      idsMesg->icontig);
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }

  dscMesg->econtig = *di;

  /* This is a new message. Hence we get a new UID */
  uidStatus = get_next_uid(&uid,real);
  if( uidStatus != UID_CODE_OK )
    {
      uidStatus = get_uids(blockSize,interval_UID,real);
      get_next_uid(&uid,real);
    }	  
  if( UID_CODE_OK != uidStatus )
    { 
      char dummy[40];
      sprintf(dummy,"Could not get UID \n");
      remove_output();
      error(AS_TER_UIDSERVER_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  dscMesg->eaccession = uid;
  SetCDS_UID_t(DSCmap,idsMesg->icontig,&uid);

  return dscMesg;
}




/************************************************/


static AugFragMesg* convert_IAF_to_AFG(IntAugFragMesg* iafMesg, int64 blockSize, int32 real)
     /* converts an IntAugFragMesg to a AugFragMessage.
	What happend ?
	- The internal frgament ID is replaced with the external UID
	- If the clear range in the IAF was corrected it is taken
	otherwise it is taken from the ClearStartMap and ClearEndMap
	- We look up in SMAmap the screen match list for the fragment.
	Then we walk through it, replacing internal screen item IDs
	with external screen item IDs.
     */
{
  CDS_UID_t *di;
  uint32 *sdi;
  uint32 *size;
  uint32 dummyZero=0;

  IntScreenMatch *pi;
  AugFragMesg *afgMesg = (AugFragMesg*) safe_malloc(sizeof(AugFragMesg));
  
  /* Set all toplevel fields */
  if( iafMesg->type == AS_BACTIG )
    di = fetch_UID_from_bactigStore(iafMesg->iaccession);
  else
    di = fetch_UID_from_fragStore(iafMesg->iaccession);

  if( di == NULL ){
    char dummy[40];
    sprintf(dummy,"No fragment with ID %d in the frag/bactigstore",iafMesg->iaccession);
    remove_output();
    error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }

  /*  fprintf(stderr,"*Fragment with iid %d has uid %ld\n",
      iafMesg->iaccession, *di); */

  afgMesg->eaccession = *di;
  afgMesg->iaccession = iafMesg->iaccession;
  afgMesg->chimeric    = iafMesg->chimeric;
  afgMesg->mate_status = iafMesg->mate_status;
  afgMesg->screened = NULL;
  afgMesg->chaff = iafMesg->chaff;

  /* check whether the clear range was corrected */
  /* the range in the iaf are different from -1 */  
  if( iafMesg->clear_rng.bgn != -1 ) {
    afgMesg->clear_rng = iafMesg->clear_rng;

  // This assertion was added by Jason, Oct 2001,
  // while adding modified clear range fields to the frag store.
  // If we reach this point in the code, then 
  // some program upstream of Terminator has written
  // a modified clear range in the IAF message.
  // That violates my assumption that all
  // clear range modifications get written to the frag store.
    assert(0);
  }
  else
    {
      if( iafMesg->type == AS_BACTIG )
	sdi = fetch_clearRange_from_bactigStore(ClearBactigStartMap,iafMesg->iaccession);
      else
	sdi = fetch_clearRange_from_fragStore(ClearStartMap,iafMesg->iaccession);
      if( sdi == NULL )
	{
	  char dummy[40];
	  sprintf(dummy,"No SeqInterval associated with ID %d in the frag/bactig store",iafMesg->iaccession);
	  remove_output();
	  error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}
      afgMesg->clear_rng.bgn = *sdi;

      if( iafMesg->type == AS_BACTIG )
	sdi = fetch_clearRange_from_bactigStore(ClearBactigEndMap,iafMesg->iaccession);
      else
	sdi = fetch_clearRange_from_fragStore(ClearEndMap,iafMesg->iaccession);
      if( sdi == NULL )
	{
	  char dummy[40];
	  sprintf(dummy,"No SeqInterval associated with ID %d in the frag store",iafMesg->iaccession);
	  remove_output();
	  error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}
      afgMesg->clear_rng.end = *sdi;
    }
    

  /* we now walk through the null terminated list of IntScreenMatches */
  if( iafMesg->type == AS_BACTIG ){
    size = &dummyZero;
  }
  else
    size = fetch_screenMatches_from_fragStore(SMASizeMap,iafMesg->iaccession);

  if( *size > 0 )
    {
      IntScreenMatch *ismMesg;
      ScreenMatch *smaMesg = NULL;
      ScreenMatch *last    = NULL;
      ScreenMatch *first   = NULL;
    
      pi = *(IntScreenMatch **)GetPtrT(SMAmap,iafMesg->iaccession);
      ismMesg = pi;
      
      while( ismMesg != NULL )
	{
	  smaMesg = (ScreenMatch*) safe_malloc(sizeof(ScreenMatch));
	  if( first == NULL )
	    first = smaMesg;

	  smaMesg->where     = ismMesg->where;
	  smaMesg->relevance = ismMesg->relevance;
	  smaMesg->portion_of= ismMesg->portion_of;
	  smaMesg->direction = ismMesg->direction;
	  /* we look up the external id of the screen item and the repeat item */
	  di = fetch_UID_from_screenStore(SCNmap,ismMesg->iwhat);
	  if( di == NULL )
	    {
	      char dummy[40];
	      sprintf(dummy,"No Screen match with ID %d in the frag store",ismMesg->iwhat);
	      remove_output();
	      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }	  
	  // BEWARE. I return NULL in the case the repeat and screen library are not
	  // submitted (normally in the simulator)
	  if( simulator ){
	    if( di == NULL )
	      smaMesg->what = ismMesg->iwhat;
	  }
	  else
	    smaMesg->what = *di;

	  di = fetch_UID_from_repeatStore(RPTmap,ismMesg->repeat_id);
	  if( di == NULL )
	    {
	      char dummy[40];
	      sprintf(dummy,"No Repeat item with ID %d in the frag store",ismMesg->repeat_id);
	      remove_output();
	      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }
	  // BEWARE. I return NULL in the case the repeat and screen library are not
	  // submitted (normally in the simulator)

	  if( simulator ){
	    if( di == NULL )
	      smaMesg->repeat_id = ismMesg->repeat_id;
	  }
	  else
	    smaMesg->repeat_id = *di;
	  
	  if( last != NULL )
	    {
	      last->next = smaMesg;
	      last = smaMesg;
	    }
	  else
	    last = smaMesg;
	  
	  ismMesg = ismMesg->next;
	}
      if( smaMesg != NULL )
	smaMesg->next = NULL;
      afgMesg->screened = first;    
    }

  return afgMesg;
}






/*****************/
/* free routines */
/*****************/

static void free_DSC(SnapDegenerateScaffoldMesg* dscMesg){
  free(dscMesg);
}

static void free_UTG(SnapUnitigMesg* utgMesg){
  /*** DO NOT FREE DELTA ARRAYS...THESE WERE COPIED BY REFERENCE */
  if( utgMesg->f_list != NULL )
    free(utgMesg->f_list);
  if( utgMesg->consensus != NULL )
    free(utgMesg->consensus);
  if( utgMesg->quality != NULL )
    free(utgMesg->quality);
  free(utgMesg);
}

static void free_ULK(SnapUnitigLinkMesg* ulkMesg){
  if( ulkMesg->jump_list != NULL )
    free(ulkMesg->jump_list);
  free(ulkMesg);
}

static void free_CCO(SnapConConMesg* ccoMesg){
  /*** DO NOT FREE DELTA ARRAYS...THESE WERE COPIED BY REFERENCE */
  if( ccoMesg->num_pieces > 0)
    free(ccoMesg->pieces);
  if( ccoMesg->num_unitigs > 0)
    free(ccoMesg->unitigs);
  if( ccoMesg->consensus != NULL )
    free(ccoMesg->consensus);
  if( ccoMesg->quality != NULL )
    free(ccoMesg->quality);
  free(ccoMesg);
}

static void free_CLK(SnapContigLinkMesg* clkMesg){
  if( clkMesg->jump_list != NULL )
    free(clkMesg->jump_list);
  free(clkMesg);
}

static void free_SLK(SnapScaffoldLinkMesg* slkMesg){
  if( slkMesg->jump_list != NULL )
    free(slkMesg->jump_list);
  free(slkMesg);
}

static void free_SCF(SnapScaffoldMesg* scfMesg){
  if( scfMesg->contig_pairs != NULL )
    free(scfMesg->contig_pairs);
  free(scfMesg);
}

static void free_MDI(SnapMateDistMesg* mdiMesg){
  if( mdiMesg->histogram != NULL )
    free(mdiMesg->histogram);
  free(mdiMesg);
}

static void free_AFG(AugFragMesg* afgMesg){
  ScreenMatch *match = afgMesg->screened;
  ScreenMatch *next  = NULL;
  while( match != NULL)
    {
      next = match->next;
      free(match);
      match = next;
    }
  free(afgMesg);
}



/********************/
/* read the stores  */
/********************/

static void read_stores(char* fragStoreName, char* bactigStoreName, char* gkpStoreName);

/************************/
/* main output routine  */
/************************/


void output_snapshot(char* fragStoreName, char* bactigStoreName, 
		     char* gkpStoreName, char** inputFileList, int32 numInputFiles,
		     char* outputFileName, char* mapFileName, OutputType output, 
		     int32 blockSize, int32 real, int32 quiet,
		     int32 random, CDS_UID_t uidStart, 
		     int argc, char *argv[])
{
  GenericMesg *pmesg       = NULL; 
  FILE        *fileInput   = NULL;
  FILE        *fileOutput  = NULL;
  FILE        *fileIID2UID = NULL;
  char        *IID2UIDFileName = NULL;
  MesgReader   readerFn   = NULL;
  MesgWriter   writerFn   = NULL;
  CDS_UID_t       interval_UID[4];
  int32        ifile;
  char        *inputFileName;
  int numIAF = 0;
  int numIUM = 0;
  int numICM = 0;

  /* Test preconditions */
  if( outputFileName == NULL )
    {
      char dummy[40];
      sprintf(dummy,"Argument outputFileName is NULL \n");
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,
	    AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }

  if( strlen(outputFileName) == 0 )
    {
      char dummy[40];
      sprintf(dummy,"Argument outputFileName is empty \n");
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,
	    AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }
  /* Test preconditions */
  if( mapFileName == NULL )
    {
      char dummy[40];
      sprintf(dummy,"Argument mapFileName is NULL \n");
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,
	    AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }

  if( strlen(mapFileName) == 0 )
    {
      char dummy[40];
      sprintf(dummy,"Argument mapFileName is empty \n");
      remove_output();
      error(AS_TER_PRECONDITION_ERROR,dummy,
	    AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }


  /* allocate dynamic arrays as maps for IID->UID */
  FRGpresent   = CreateVA_short(ARRAYSIZE);
  BTGpresent   = CreateVA_short(ARRAYSIZE);
  ICM1inISFmap = CreateVA_int(ARRAYSIZE);
  ICM2inISFmap = CreateVA_int(ARRAYSIZE);
  IUMmap = CreateVA_CDS_UID_t(ARRAYSIZE);
  ICMmap = CreateVA_CDS_UID_t(ARRAYSIZE);
  ISFmap = CreateVA_CDS_UID_t(ARRAYSIZE);
  FRGmap = CreateVA_CDS_UID_t(ARRAYSIZE);
  BTGmap = CreateVA_CDS_UID_t(ARRAYSIZE);
  DSCmap = CreateVA_CDS_UID_t(ARRAYSIZE);
  DSTmap = CreateVA_CDS_UID_t(ARRAYSIZE);
  SCNmap = CreateVA_CDS_UID_t(ARRAYSIZE);
  RPTmap = CreateVA_CDS_UID_t(ARRAYSIZE);
  ClearStartMap = CreateVA_uint32(ARRAYSIZE);
  ClearEndMap   = CreateVA_uint32(ARRAYSIZE);  
  ClearBactigStartMap = CreateVA_uint32(ARRAYSIZE);
  ClearBactigEndMap   = CreateVA_uint32(ARRAYSIZE);  

  /* allocate a dynamic parray for keeping track of the ISM list
     of each fragment in the frag store */
  SMAmap       = CreateVA_PtrT(ARRAYSIZE);
  SMASizeMap   = CreateVA_uint32(ARRAYSIZE);

  /* set the UID start to the number given by uidStart */
  AS_TER_uidStart = uidStart;

  /* set the simulator bool to the value of quiet (-Q flag) */
  simulator = quiet;

  /* if no bactig store is given we are in grande mode */
  if( bactigStoreName == NULL )
    fprintf(stderr,"*** Terminator : Assembly Grande ***\n");
  else
    fprintf(stderr,"*** Terminator : Overlay Assembly ***\n");

  /* if no gatekeeper store is given we assume it to reside in the frag store */
  if( gkpStoreName == NULL )
    gkpStoreName = fragStoreName;

//  fprintf(stderr,"%s\n",gkpStoreName);
//  fprintf(stderr,"%s\n",fragStoreName);
 
  /* Read the The Frag store, the gatekeeper Store 
     and maybe the Bactig store */  
  // first initialize the global static variables
  if( existsFragStore(fragStoreName) == FALSE ){
    char dummy[60];
    sprintf(dummy,"Frag Store %s does not exist \n",fragStoreName);
    error(AS_TER_PRECONDITION_ERROR,dummy,
          AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
  }
  FSHandle     = openFragStore(fragStoreName,"r");
  if( bactigStoreName != NULL )
    BSHandle     = openFragStore(bactigStoreName,"r");
  InitGateKeeperStore(&GKPStore,gkpStoreName);
  OpenReadOnlyGateKeeperStore(&GKPStore);
  
  if( random == FALSE )
    read_stores(fragStoreName,bactigStoreName,gkpStoreName);
  

  IID2UIDFileName    = mapFileName;
  OutFileNameREZ     = outputFileName;
  IID2UIDFileNameREZ = IID2UIDFileName;




  if( strcmp(outputFileName,"-") == 0 ){
    PipeOut = TRUE;
    fprintf(stderr,"*** Terminator : Print output to stdout ***\n");
  }
    
  if( ! PipeOut )
    fileOutput  = file_open(outputFileName,"w");  
  else
    fileOutput = stdout;

  fileIID2UID = file_open(IID2UIDFileName,"w");  


  /* set static global variables for emergency exit */
  OFPREZ = fileOutput;

  /* Allocate the buffer of UIDS. If real is TRUE, the server is queried
     Otherwise a dummy contiguous number is assigned */
  get_uids(blockSize,interval_UID,real);


  for(ifile = 0; ifile < numInputFiles; ifile++){
    inputFileName = inputFileList[ifile];


    /* Test preconditions */
    if( inputFileName == NULL )
      {
	char dummy[40];
	sprintf(dummy,"Argument inputFileName is NULL \n");
	remove_output();
	error(AS_TER_PRECONDITION_ERROR,dummy,
	      AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
      }

    if( strlen(inputFileName) == 0 )
      {
	char dummy[40];
	sprintf(dummy,"Argument inputFileName is empty \n");
	remove_output();
	error(AS_TER_PRECONDITION_ERROR,dummy,
	      AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
      }

    if( strcmp(inputFileName,"-") == 0 ){
      PipeIn = TRUE;
      fprintf(stderr,"*** Terminator : Read input from stdin ***\n");
    }

    if( ! PipeIn )
      fileInput = file_open(inputFileName,"r");  
    else
      fileInput = stdin;

    /* set static global variables for emergency exit */
    IFPREZ = fileInput;

    /* detect whether reader is ASCII or binary and set
       writer to provided parameter */
    readerFn  = InputFileType_AS(fileInput);
    writerFn  = OutputFileType_AS(output);


    /* main loop  */
    fprintf(stderr,"*** Terminator : Reading cgw file %d %s ***\n", ifile, inputFileName);

    while(readerFn(fileInput,&pmesg) != EOF){
#if DEBUG > 0
      //    writerFn(stderr,pmesg);    
#endif
      switch(pmesg->t){
      case MESG_ADT : {
	AuditMesg *adt_mesg;
	adt_mesg = pmesg->m;
	pmesg->t = MESG_ADT;
      
	VersionStampADT(adt_mesg,argc,argv);
	writerFn(fileOutput,pmesg);    
	writerFn(stderr,pmesg);    
      }
      break; 
      case MESG_IAF :{
	IntAugFragMesg *iafMesg = (IntAugFragMesg*) pmesg->m;
	AugFragMesg *afgMesg = convert_IAF_to_AFG(iafMesg,blockSize,real);
	pmesg->m = afgMesg;
	pmesg->t = MESG_AFG;
	writerFn(fileOutput,pmesg);
	free_AFG(afgMesg);
	if((++numIAF % 1000000) == 0){
	  fprintf(stderr,"* Processed %d frags\n", numIAF);
	  fflush(stderr);
	}
      }
      break;
      case MESG_IUM :{
	IntUnitigMesg  *iumMesg = (IntUnitigMesg*) pmesg->m;
	SnapUnitigMesg *utgMesg = convert_IUM_to_UTG(iumMesg,blockSize,real);
	pmesg->m = utgMesg;
	pmesg->t = MESG_UTG;
	writerFn(fileOutput,pmesg);
	free_UTG(utgMesg);
	if((++numIUM % 1000000) == 0){
	  fprintf(stderr,"* Processed %d unitigs\n", numIUM);
	  fflush(stderr);
	}
      }
      break;
      case MESG_IUL : {
	IntUnitigLinkMesg  *iulMesg = (IntUnitigLinkMesg*) pmesg->m;
	SnapUnitigLinkMesg *ulkMesg = convert_IUL_to_ULK(iulMesg,blockSize,real);
	pmesg->m = ulkMesg;
	pmesg->t = MESG_ULK;
	writerFn(fileOutput,pmesg);
	free_ULK(ulkMesg);   
      }
      break;
      case MESG_ICM : {
	IntConConMesg *icmMesg = (IntConConMesg*) pmesg->m;
	SnapConConMesg *ccoMesg = convert_ICM_to_CCO(icmMesg,blockSize,real);
	pmesg->m = ccoMesg;
	pmesg->t = MESG_CCO;
	writerFn(fileOutput,pmesg);
	free_CCO(ccoMesg);    
	if((++numICM % 1000000) == 0){
	  fprintf(stderr,"* Processed %d contigs\n", numICM);
	  fflush(stderr);
	}
      }
      break;
      case MESG_ICL : {
	IntContigLinkMesg *iclMesg = (IntContigLinkMesg*) pmesg->m;
	SnapContigLinkMesg *clkMesg = convert_ICL_to_CLK(iclMesg,blockSize,real);
	pmesg->m = clkMesg;
	pmesg->t = MESG_CLK;
	writerFn(fileOutput,pmesg);
	free_CLK(clkMesg);    
      }
      break;
      case MESG_ISL : {
	InternalScaffoldLinkMesg *islMesg = (InternalScaffoldLinkMesg*) pmesg->m;
	SnapScaffoldLinkMesg *slkMesg = convert_ISL_to_SLK(islMesg,blockSize,real);
	pmesg->m = slkMesg;
	pmesg->t = MESG_SLK;
	writerFn(fileOutput,pmesg);
	free_SLK(slkMesg);    
      }
      break;
      case MESG_ISF : {
	IntScaffoldMesg *isfMesg = (IntScaffoldMesg*) pmesg->m;
	SnapScaffoldMesg *scfMesg = convert_ISF_to_SCF(isfMesg,blockSize,real);
	pmesg->m = scfMesg;
	pmesg->t = MESG_SCF;
	writerFn(fileOutput,pmesg);
	free_SCF(scfMesg);    
      }
      break;
      case MESG_IMD : {
	IntMateDistMesg *imdMesg = (IntMateDistMesg*) pmesg->m;
	SnapMateDistMesg *mdiMesg = convert_IMD_to_MDI(imdMesg,blockSize,real);
	pmesg->m = mdiMesg;
	pmesg->t = MESG_MDI;
	writerFn(fileOutput,pmesg);
	free_MDI(mdiMesg);    
      }
      break;
      case MESG_IDS : {
	IntDegenerateScaffoldMesg *idsMesg = (IntDegenerateScaffoldMesg*) pmesg->m;
	SnapDegenerateScaffoldMesg *dscMesg = convert_IDS_to_DSC(idsMesg,blockSize,real);
	pmesg->m = dscMesg;
	pmesg->t = MESG_DSC;
	writerFn(fileOutput,pmesg);
	free_DSC(dscMesg);    
      }
      break;
      default:
	{
#if DEBUG > 0
	  fprintf(stderr,"Swallowing message of type %d\n",pmesg->t); 
	  writerFn(stderr,pmesg);
#endif
	}
      break;
      }
    }
  }

  fprintf(stderr,"*** Terminator : Successfully generated genome snapshot ***\n");
  fflush(NULL);

  fclose(fileInput);
  fclose(fileOutput);
  
  fprintf(fileIID2UID,"Fragment IID2UID map\n");
  DumpIID2UIDmap(FRGmap,fileIID2UID);
  fprintf(fileIID2UID,"Fragment IID2UID map\n");
  DumpIID2UIDmap(FRGmap,fileIID2UID);
  fprintf(fileIID2UID,"Unitig IID2UID map\n");
  DumpIID2UIDmap(IUMmap,fileIID2UID);
  fprintf(fileIID2UID,"Contig IID2UID map\n");
  DumpIID2UIDmap(ICMmap,fileIID2UID);
  fprintf(fileIID2UID,"Scaffold IID2UID map\n");
  DumpIID2UIDmap(ISFmap,fileIID2UID);
  fprintf(fileIID2UID,"Bactig IID2UID map\n");
  DumpIID2UIDmap(BTGmap,fileIID2UID);
  fprintf(fileIID2UID,"Distrib IID2UID map\n");
  DumpIID2UIDmap(DSTmap,fileIID2UID);
  fprintf(fileIID2UID,"Screen Item IID2UID map\n");
  DumpIID2UIDmap(SCNmap,fileIID2UID);
  fprintf(fileIID2UID,"Repeat Item IID2UID map\n");
  DumpIID2UIDmap(RPTmap,fileIID2UID);
  fprintf(fileIID2UID,"Degenerate Contig IID to Scaffold UID map\n");
  DumpIID2UIDmap(DSCmap,fileIID2UID);

  fflush(NULL);
  fclose(fileIID2UID);

  closeFragStore(FSHandle);
  

  fprintf(stderr,"*** Terminator : Successfully generated IID to UID map ***\n");

  /* deallocate dynamic arrays as maps for IID->UID */
  DeleteVA_short(FRGpresent);
  DeleteVA_short(BTGpresent);
  DeleteVA_int(ICM1inISFmap);
  DeleteVA_int(ICM2inISFmap);
  DeleteVA_CDS_UID_t(IUMmap);
  DeleteVA_CDS_UID_t(ICMmap);
  DeleteVA_CDS_UID_t(ISFmap);
  DeleteVA_CDS_UID_t(FRGmap);
  DeleteVA_CDS_UID_t(BTGmap);
  DeleteVA_CDS_UID_t(DSCmap);
  DeleteVA_CDS_UID_t(DSTmap);
  DeleteVA_CDS_UID_t(SCNmap);
  DeleteVA_CDS_UID_t(RPTmap);
  DeleteVA_uint32(ClearStartMap);
  DeleteVA_uint32(ClearEndMap);
  DeleteVA_uint32(ClearBactigStartMap);
  DeleteVA_uint32(ClearBactigEndMap);
  DeleteVA_uint32(SMASizeMap);
  DeleteVA_PtrT(SMAmap);
}






/* read information from the fragment, gatekeeper and bactig store */


void read_stores(char* fragStoreName, char* bactigStoreName, char* gkpStoreName){
  ReadStructp  input;
  //  FragStoreHandle storeHandle = 0;
  FragStoreHandle streamHandle = 0;
  CDS_UID_t *di;
  int32 first = TRUE;

  fprintf(stderr,"*** Terminator : Reading Fragment Store ***\n");

  /* Test preconditions */
  if( fragStoreName == NULL ){
    char dummy[40];
    sprintf(dummy,"Argument fragStoreName is NULL \n");
    error(AS_TER_PRECONDITION_ERROR,dummy,
	  AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
  }

  if( strlen(fragStoreName) == 0 ){
    char dummy[40];
    sprintf(dummy,"Argument fragStoreName is empty \n");
    error(AS_TER_PRECONDITION_ERROR,dummy,
	  AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
  }

    
  if( existsFragStore(fragStoreName) == FALSE ){
    char dummy[60];
    sprintf(dummy,"Frag Store %s does not exist \n",fragStoreName);
    error(AS_TER_PRECONDITION_ERROR,dummy,
	  AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
  }

  streamHandle = openFragStream(FSHandle,NULL,0);

  input =  new_ReadStruct();
  while( nextFragStream(streamHandle,input,FRAG_S_FIXED) ){
    CDS_IID_t iid;
    CDS_UID_t uid;
    uint32 cStart,cEnd;
    // CDS_UID_t *di;
    IntScreenMatch *pi;
    IntScreenMatch *matches;
    uint32         noOfMatches;
    short truedummy  = TRUE;

    /* read the iid and uid pair */
    getReadIndex_ReadStruct(input,&iid);
    getAccID_ReadStruct(input,&uid);

    if(((iid+1) % 1000000) == 0){
      fprintf(stderr,"* Read frag %d from fragStore\n", iid);
      fflush(stderr);
    }
    /* test whether the iid is not already present */
    //di = GetCDS_UID_t(FRGmap,iid);
    //    if( di != NULL ){
    if( test_frg_present(iid) ){
      char dummy[40];
      sprintf(dummy,"Internal fragment ID %d occurred twice (FRGmap)",iid);
      error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
    SetCDS_UID_t(FRGmap,iid,&uid);
    Setshort(FRGpresent,iid,&truedummy);

    if( first ){
      AS_TER_uidStart = uid+1;
      first = FALSE;
    }
    else
      if(AS_TER_uidStart <= uid)
	AS_TER_uidStart = uid+1;
    
#if DEBUG > 1
    fprintf(stderr,"AS_TER_uidStart after reading frags = %lu\n",AS_TER_uidStart);
#endif

    /* get the clear region from the frag store */
    getClearRegion_ReadStruct(input,&cStart,&cEnd,READSTRUCT_LATEST);
    Setuint32(ClearStartMap,iid,&cStart);
    Setuint32(ClearEndMap,iid,&cEnd);

    /* record the Screen Matches of each fragment */
    /* first get the number of screen matches and allocate the buffer */
    noOfMatches = getScreenMatches_ReadStruct(input,NULL,0);

#if DEBUG > 1
    fprintf(stderr,"Found %u screen matches for fragment %u\n",noOfMatches,iid);
#endif

    /* record the number of matches. This is used in the IFG_to_AFG conversion
       to test whether the list is empty */
    Setuint32(SMASizeMap,iid,&noOfMatches);
    if(noOfMatches > 0){
      matches = (IntScreenMatch*) safe_calloc(noOfMatches,sizeof(IntScreenMatch));
      getScreenMatches_ReadStruct(input,matches,noOfMatches);
      
      pi = (IntScreenMatch *)GetPtrT(SMAmap,iid);
      if( pi != NULL ){
	char dummy[40];
	sprintf(dummy,"Two screen match lists for internal ID %d in the frag store",iid);
	error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
      }
      SetPtrT(SMAmap,iid,(PtrT)&matches);
    }
    
  }
  

  closeFragStream(streamHandle);
  /* open the gatekeeper store */
  fprintf(stderr,"*** Terminator : Reading Gatekeeper Store ***\n");
  {
    /* reading distributions */
    {
     GateKeeperDistanceRecord gkpd;
     StoreStat stat;
     int i ;
     statsStore(GKPStore.dstStore, &stat);
#if DEBUG > 0
     fprintf(stderr,"* Stats for Dist Store are first:%lu last :%lu\n",
	     stat.firstElem, stat.lastElem);
#endif
     for(i = stat.firstElem; i <= stat.lastElem; i++){
       getGateKeeperDistanceStore(GKPStore.dstStore,i,&gkpd);
#if DEBUG > 1
       fprintf(stderr,"* Dist %d UID:%lu del:%d red:%d mean:%f std:%f batch(%d,%d) prevID:%d prevInstanceID:%d\n",
	       i,gkpd.UID, gkpd.deleted, gkpd.redefined, gkpd.mean, gkpd.stddev,
	       gkpd.birthBatch, gkpd.deathBatch, gkpd.prevID, gkpd.prevInstanceID);
#endif

       di = GetCDS_UID_t(DSTmap,i);
       if( di != NULL ){
	 char dummy[40];
	 sprintf(dummy,"Internal DST ID %d occurred twice",i);
	 error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
       }
       SetCDS_UID_t(DSTmap,i,&gkpd.UID);
       if( AS_TER_uidStart <= gkpd.UID )
	 AS_TER_uidStart = gkpd.UID+1;

#if DEBUG > 1
       fprintf(stderr,"AS_TER_uidStart after reading distribs %lu\n",AS_TER_uidStart);
#endif
     }
    }

    /* reading RPT IID to UID mapping */
    {
      GateKeeperRepeatRecord gkpr;
      StoreStat stat;
      int i ;
      statsStore(GKPStore.rptStore, &stat);
#if DEBUG > 0
      fprintf(stderr,"* Stats for Repeat Store are first:%lu last :%lu\n",
	      stat.firstElem, stat.lastElem);
#endif
     for(i = stat.firstElem; i <= stat.lastElem; i++){
       getGateKeeperRepeatStore(GKPStore.rptStore,i,&gkpr);
#if DEBUG > 1
       fprintf(stderr,"* Repeat %d UID:%lu which:%s \n",
	       i, gkpr.UID, gkpr.which);
#endif
       di = GetCDS_UID_t(RPTmap,i);
       if( di != NULL ){
	 char dummy[40];
	 sprintf(dummy,"Internal RPT ID %d occurred twice",i);
	 error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
       }
       SetCDS_UID_t(RPTmap,i,&gkpr.UID);
       if(AS_TER_uidStart <= gkpr.UID)
	 AS_TER_uidStart = gkpr.UID+1;
     }
    }
    /* reading SCN IID to UID mapping */
    {
      GateKeeperScreenRecord gkps;
      StoreStat stat;
      int i ;
      statsStore(GKPStore.scnStore, &stat);
#if DEBUG > 0
      fprintf(stderr,"* Stats for Screen Store are first:%lu last :%lu\n",
	     stat.firstElem, stat.lastElem);
#endif
     for(i = stat.firstElem; i <= stat.lastElem; i++){
       getGateKeeperScreenStore(GKPStore.scnStore,i,&gkps);
#if DEBUG > 1
       fprintf(stderr,"* Screen %d UID:%lu repeatID:%d batch:(%d,%d)\n",
	       i, gkps.UID, gkps.repeatID, gkps.birthBatch, gkps.deathBatch);
#endif
       
       di = GetCDS_UID_t(SCNmap,i);
       if( di != NULL ){
	 char dummy[40];
	 sprintf(dummy,"Internal SCN ID %d occurred twice",i);
	 error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
       }
       SetCDS_UID_t(SCNmap,i,&gkps.UID);
       if(AS_TER_uidStart <= gkps.UID)
	 AS_TER_uidStart = gkps.UID+1;
     }
    } 
#if DEBUG > 1
       fprintf(stderr,"AS_TER_uidStart after reading repeats %lu\n",AS_TER_uidStart);
#endif
 

  }

  // from the bactig store we only want to know
  // its uid and the clear range
  
  if( bactigStoreName != NULL ){
    fprintf(stderr,"*** Terminator : Reading Bactig Store ***\n");
    if( strlen(bactigStoreName) == 0 ){
      char dummy[40];
      sprintf(dummy,"Argument bactigStoreName is empty \n");
      error(AS_TER_PRECONDITION_ERROR,dummy,
	    AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }
    
    if( existsFragStore(bactigStoreName) == FALSE ){
      char dummy[60];
      sprintf(dummy,"Bactig Store %s does not exist \n",fragStoreName);
      error(AS_TER_PRECONDITION_ERROR,dummy,
	    AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }
    
    streamHandle = openFragStream(BSHandle,NULL,0);
    
    input = new_ReadStruct();
    while( nextFragStream(streamHandle,input,FRAG_S_ALL) ){
      short truedummy  = TRUE;
      CDS_IID_t iid;
      CDS_UID_t uid;
      uint32 cStart,cEnd;
      // CDS_UID_t *di;
      
      /* read the iid and uid pair */
      getReadIndex_ReadStruct(input,&iid);
      getAccID_ReadStruct(input,&uid);
      
      /* test whether the iid is not already present */
      //di = GetCDS_UID_t(BTGmap,iid);
      //      if( di != NULL )
      if( test_btg_present(iid) ){
	char dummy[40];
	sprintf(dummy,"Internal bactig ID %d occurred twice (BTGmap)",iid);
	error(AS_TER_PRECONDITION_ERROR,dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
      }
      SetCDS_UID_t(BTGmap,iid,&uid);
      Setshort(BTGpresent,iid,&truedummy);
      
      if(AS_TER_uidStart <= uid)
	AS_TER_uidStart = uid+1;

      /* get the clear region from the bactig store */
      getClearRegion_ReadStruct(input,&cStart,&cEnd,READSTRUCT_LATEST);
      Setuint32(ClearBactigStartMap,iid,&cStart);
      Setuint32(ClearBactigEndMap,iid,&cEnd);
    }
  }
}


