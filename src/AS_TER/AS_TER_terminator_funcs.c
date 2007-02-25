
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

static char CM_ID[] = "$Id: AS_TER_terminator_funcs.c,v 1.31 2007-02-25 08:13:38 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

#include "AS_UTL_Var.h"
#include "AS_UTL_version.h"

#include "AS_TER_terminator_funcs.h"

#include "SYS_UIDcommon.h"
#include "SYS_UIDclient.h"

VA_DEF(short);
VA_DEF(int);
VA_DEF(uint32);
VA_DEF(CDS_UID_t);

#define DEBUG       0
#define DEBUG_UID   0

/*--------------------------------------------------------------------*/
/*  Output Routine for all Messages in the genome snapshot. The static*/
/*  functions allocate a new message and make sure that the IIDs are  */
/*  correctly replaced by the UIDs                                    */
/*--------------------------------------------------------------------*/

static int32 simulator = FALSE;

static VA_TYPE(short)  *FRGpresent;
static VA_TYPE(short)  *BTGpresent;
static VA_TYPE(int)    *ICM1inISFmap;
static VA_TYPE(int)    *ICM2inISFmap;
static VA_TYPE(CDS_UID_t) *IUMmap;
static VA_TYPE(CDS_UID_t) *ICMmap;
static VA_TYPE(CDS_UID_t) *ISFmap;
static VA_TYPE(CDS_UID_t) *FRGmap;
static VA_TYPE(CDS_UID_t) *BTGmap;
static VA_TYPE(CDS_UID_t) *DSCmap;
static VA_TYPE(CDS_UID_t) *DSTmap;
static VA_TYPE(uint32) *ClearStartMap;
static VA_TYPE(uint32) *ClearEndMap;

static int PipeIn=FALSE;
static int PipeOut=FALSE;

// the store handles are global such that we can open them in
// output_snapshot and use them in read_stores and the various fetch
// functions
//
static GateKeeperStore *FSHandle;
static GateKeeperStore *BSHandle;



//  Initial size of our arrays
//
#define ARRAYSIZE  2048



// For reporting errors
//
static char errorreport[1024];

#define   AS_TER_EXIT_SUCCESS 0
#define   AS_TER_EXIT_FAILURE 1

void
error(const char* errorMesg,
      int ex, 
      const char* file,
      int line) {

  fprintf(stderr,"%s:%d %s\n", file, line, errorMesg);
  exit(ex);
}

FILE*
file_open(const char* fileName, const char* mode) {
   FILE* fp;

   fp = fopen (fileName,mode);
   if(fp == NULL) {
     char dummy[40];
     sprintf(dummy,"Could not open file '%s'",fileName);
     error(dummy,AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
   }
   return  fp;
}




void DumpIID2UIDmap(VA_TYPE(CDS_UID_t) *map,FILE *file){
  int i;
  for(i=0; i<GetNumCDS_UID_ts(map); i++){
      CDS_UID_t *di;
      di = GetCDS_UID_t(map,i);
      if (*di != 0)
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



static CDS_UID_t *fetch_UID_from_fragStore(CDS_IID_t iid){
  short truedummy=TRUE;
  fragRecord *input;
  CDS_UID_t uid;
  CDS_UID_t *di;
  uint32 cStart, cEnd;

  if( test_frg_present(iid) ){
    CDS_UID_t* ret = fetch_UID(FRGmap,iid);
    assert( ret != NULL);
    return ret;
  }

  assert( ! test_frg_present(iid) );
  Setshort(FRGpresent,iid,&truedummy);
    
  input = new_fragRecord();
  getFrag(FSHandle,iid,input,FRAG_S_INF);

  uid = getFragRecordUID(input);

  di = GetCDS_UID_t(FRGmap,iid);
  if ((di != NULL) && (*di != 0)) {
    sprintf(errorreport,"Internal fragment ID %d occurred twice (FRGmap)",iid);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }

  SetCDS_UID_t(FRGmap,iid,&uid);
    
  cStart = getFragRecordClearRegionBegin(input, AS_READ_CLEAR_LATEST);
  cEnd   = getFragRecordClearRegionEnd  (input, AS_READ_CLEAR_LATEST);

  Setuint32(ClearStartMap,iid,&cStart);
  Setuint32(ClearEndMap,iid,&cEnd);

  del_fragRecord(input); 

  return GetCDS_UID_t(FRGmap,iid);   
}
  

static uint32 *fetch_clearRange_from_fragStore(VA_TYPE(uint32) *map, CDS_IID_t iid){
  uint32* ret;
  assert( test_frg_present(iid) );
  ret = fetch_range_with_null(map,iid);
  assert( ret != NULL);
  return ret;
}



static CDS_UID_t *fetch_UID_from_distStore(VA_TYPE(CDS_UID_t) *map, CDS_IID_t iid){
  CDS_UID_t* ret = fetch_UID(map,iid);
  if(ret == NULL)
  {
    GateKeeperLibraryRecord gkpl;
    if( 0 == getGateKeeperLibraryStore(FSHandle->lib,iid,&gkpl) ){
      CDS_UID_t *di;
      di = GetCDS_UID_t(map,iid);
      if ((di != NULL) && (*di != 0)) {
        sprintf(errorreport,"Internal DST ID %d occurred twice",iid);
        error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
      }
      SetCDS_UID_t(map,iid,&gkpl.UID);
      return GetCDS_UID_t(map,iid);
    }
    else{
      sprintf(errorreport,"Internal dist ID %d is not present in the dist store",iid);
      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
  }
  return ret;
}



// This is an inelegant way to achieve large block sizes at VI
// and small block sizes at TIGR. Later, define an EUID service
// configuration file, to be parsed at run time, containing
// primary and backup URLs plus the block size.
// Also, take the opportunity to clean up the dozens of places
// in the code where block size is set explicitly, usually to 300.
// -- Jason Miller, 8/22/2006.
#ifdef USE_SOAP_UID
#define UID_BLOCK_SIZE 100
#else
#define UID_BLOCK_SIZE 1024
#endif

static
void
initUID(int32 real) {
  CDS_UID_t  interval_UID[4];

  //  Allocate the buffer of UIDS. If real is TRUE, the server is
  //  queried, otherwise a dummy contiguous number is assigned
  //
  get_uids(UID_BLOCK_SIZE, interval_UID, real);
}


static
CDS_UID_t
getUID(int32 real) {
  CDS_UID_t  uid = 0;
  CDS_UID_t  interval_UID[4];
  int32      uidStatus = 0;

  uidStatus = get_next_uid(&uid, real);
  if (uidStatus != UID_CODE_OK) {
    get_uids(UID_BLOCK_SIZE, interval_UID, real);
    uidStatus = get_next_uid(&uid, real);
  }	  
  if (uidStatus != UID_CODE_OK) {
    error("Could not get UID.", AS_TER_EXIT_FAILURE, __FILE__, __LINE__); 
  }
  return(uid);
}



/***********************/
/* conversion routines */
/***********************/


static SnapUnitigMesg* convert_IUM_to_UTG(IntUnitigMesg* iumMesg, int32 real)
     /*
       converts an IntUnitigMesg to a SnapUnitigMessage.
       What happend ?
       - A UID is assigned and recorded in IUMmap 
     */
{
  int i;
  CDS_UID_t uid;

  CDS_UID_t *di;
  SnapUnitigMesg *utgMesg = 
    (SnapUnitigMesg*) safe_malloc(sizeof(SnapUnitigMesg));

#if DEBUG > 1
  fprintf(stderr,"IUM internal acc %u\n",iumMesg->iaccession);
#endif
  di = fetch_UID(IUMmap,iumMesg->iaccession);
  if( (di != NULL) && (*di != 0)) {
    sprintf(errorreport,"Spotted IUM internal ID %d second time",
            iumMesg->iaccession);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }

  /* This is a new message. We assign it a UID */
  uid = getUID(real);
  SetCDS_UID_t(IUMmap,iumMesg->iaccession,&uid);
  
  /* Set all toplevel fields */
  
  utgMesg->eaccession      = uid;
  utgMesg->iaccession      = iumMesg->iaccession;
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
  utgMesg->num_vars       = 0;
  utgMesg->f_list         = NULL;
  utgMesg->v_list         = NULL;

  if( iumMesg->num_frags > 0 ){
    utgMesg->f_list = (SnapMultiPos*) safe_malloc(iumMesg->num_frags*sizeof(SnapMultiPos));

    for(i=0; i<iumMesg->num_frags; i++){
      utgMesg->f_list[i].type = iumMesg->f_list[i].type;
#ifdef AS_ENABLE_SOURCE
      utgMesg->f_list[i].source = NULL;
#endif
      di = fetch_UID_from_fragStore(iumMesg->f_list[i].ident);
      
      if( di == NULL ){
	sprintf(errorreport,"Reference before definition for fragment ID %d",
		iumMesg->f_list[i].ident);
	error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
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


static SnapUnitigLinkMesg* convert_IUL_to_ULK(IntUnitigLinkMesg* iulMesg, int32 real)
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
  if ((di == NULL) || (*di == 0)) {
    sprintf(errorreport,"IUL reference before definition error for unitig ID %d",iulMesg->unitig1);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }
  ulkMesg->eunitig1 = *di;
  di = GetCDS_UID_t(IUMmap,iulMesg->unitig2);
  if ((di == NULL) || (*di == 0)) {
    sprintf(errorreport,"IUL reference before definition error for unitig ID %d",iulMesg->unitig2);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
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
     
      di = fetch_UID_from_fragStore(iulMesg->jump_list[i].in1);
      if( di == NULL )
	{
	  sprintf(errorreport,"Internal Fragment ID %d does not exist in Fragstore",iulMesg->jump_list[i].in1);
	  error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}  
      ulkMesg->jump_list[i].in1 = *di;
      di = fetch_UID_from_fragStore(iulMesg->jump_list[i].in2);
      if( di == NULL )
	{
	  sprintf(errorreport,"Internal Fragment ID %d does not exist in Fragstore",iulMesg->jump_list[i].in2);
	  error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}  
      ulkMesg->jump_list[i].in2 = *di; 
      // the next line is due to SAULS change in proto I/O
      ulkMesg->jump_list[i].type = iulMesg->jump_list[i].type;
    }
  
  return ulkMesg;
}


/************************************************/


static SnapConConMesg* convert_ICM_to_CCO(IntConConMesg* icmMesg, int32 real)
     /* converts an IntConConMesg to a SnapConConMessage.
	What happend ?
	- A new UID is assigned to the ConConMesg
	- We initialize the ICMinISFmap to FALSE
     */
{
  int i;
  CDS_UID_t uid;
  CDS_UID_t *di;
  int false = FALSE;

  SnapConConMesg *ccoMesg = (SnapConConMesg*) safe_malloc(sizeof(SnapConConMesg));
  
  /* Set all toplevel fields */

  /* we assume that the numbers are in ascending order */
  di = GetCDS_UID_t(ICMmap,icmMesg->iaccession);
  if ((di != NULL) && (*di != 0)) {
    sprintf(errorreport,"ICM internal contig number %d spotted second time",
	      icmMesg->iaccession);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }

  uid = getUID(real);
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
  ccoMesg->num_vars   = icmMesg->num_vars;  // affects .asm/CCO
  ccoMesg->pieces     = NULL;
  ccoMesg->vars       = NULL;
  ccoMesg->unitigs    = NULL;

  if (ccoMesg->num_vars > 0) {
     ccoMesg->vars = (IntMultiVar*) safe_malloc(icmMesg->num_vars*sizeof(IntMultiVar));
     for(i=0; i<icmMesg->num_vars; i++) // i loop
     {
        ccoMesg->vars[i].position         = icmMesg->v_list[i].position;
        ccoMesg->vars[i].num_reads        = icmMesg->v_list[i].num_reads; 
        ccoMesg->vars[i].num_conf_alleles = icmMesg->v_list[i].num_conf_alleles;
        ccoMesg->vars[i].anchor_size      = icmMesg->v_list[i].anchor_size;
        ccoMesg->vars[i].var_length       = icmMesg->v_list[i].var_length ;
        ccoMesg->vars[i].nr_conf_alleles  = strdup(icmMesg->v_list[i].nr_conf_alleles);
        ccoMesg->vars[i].weights          = strdup(icmMesg->v_list[i].weights);
        ccoMesg->vars[i].var_seq          = strdup(icmMesg->v_list[i].var_seq);  
     }
      
  }

  if( ccoMesg->num_pieces > 0 ){ 
    ccoMesg->pieces = (SnapMultiPos*) safe_malloc(icmMesg->num_pieces*sizeof(SnapMultiPos));
    for(i=0; i<icmMesg->num_pieces; i++){// i loop
      ccoMesg->pieces[i].type = icmMesg->pieces[i].type;
#ifdef AS_ENABLE_SOURCE
      ccoMesg->pieces[i].source = NULL;
#endif

#if DEBUG > 1
      fprintf(stderr,"ICM = %d, no of pieces = %d, piece = %d, piece id = %d \n",icmMesg->iaccession,icmMesg->num_pieces,i,icmMesg->pieces[i].ident);
#endif

      di = fetch_UID_from_fragStore(icmMesg->pieces[i].ident);
      if( di == NULL ){
	sprintf(errorreport,"Reference before definition for fragment ID %d",
		icmMesg->pieces[i].ident);
	error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
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
      if ((di == NULL) || (*di == 0)) {
	sprintf(errorreport,"Reference before definition for unitig ID %d",
		icmMesg->unitigs[i].ident);
	error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
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

static SnapContigLinkMesg* convert_ICL_to_CLK(IntContigLinkMesg* iclMesg, int32 real)
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
  if ((di == NULL) || (*di == 0)) {
    sprintf(errorreport,"ICL reference before definition error for contig ID %d",iclMesg->contig1);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }
  clkMesg->econtig1 = *di;
  di = GetCDS_UID_t(ICMmap,iclMesg->contig2);
  if ((di == NULL) || (*di == 0)) {
    sprintf(errorreport,"ICL reference before definition error for contig ID %d",iclMesg->contig2);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
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
      sprintf(errorreport,"Internal Fragment ID %d does not exist in Fragstore",iclMesg->jump_list[i].in1);
      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }  
    clkMesg->jump_list[i].in1 = *di;
    
    di = fetch_UID_from_fragStore(iclMesg->jump_list[i].in2);
    if( di == NULL ){
      sprintf(errorreport,"Internal Fragment ID %d does not exist in Fragstore",iclMesg->jump_list[i].in2);
      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }  
    clkMesg->jump_list[i].in2 = *di; 
    // the next line is due to SAULS change in proto I/O
    clkMesg->jump_list[i].type = iclMesg->jump_list[i].type;
  }
  
  return clkMesg;
}



static SnapScaffoldLinkMesg* convert_ISL_to_SLK(InternalScaffoldLinkMesg* islMesg, int32 real)
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
  if ((di == NULL) || (*di == 0)) {
    sprintf(errorreport,"ISL reference before definition error for scaffold ID %d",islMesg->iscaffold1);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }
  slkMesg->escaffold1 = *di;
  di = GetCDS_UID_t(ISFmap,islMesg->iscaffold2);
  if ((di == NULL) || (*di == 0)) {
    sprintf(errorreport,"ISL reference before definition error for scaffold ID %d",islMesg->iscaffold2);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
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
      sprintf(errorreport,"Internal Fragment ID %d does not exist in Fragstore",islMesg->jump_list[i].in1);
      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }  
    slkMesg->jump_list[i].in1 = *di;
    
    di = fetch_UID_from_fragStore(islMesg->jump_list[i].in2);
    if( di == NULL ){
      sprintf(errorreport,"Internal Fragment ID %d does not exist in Fragstore",islMesg->jump_list[i].in2);
      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }  
    slkMesg->jump_list[i].in2 = *di; 
    // the next line is due to SAULS change in proto I/O
    slkMesg->jump_list[i].type = islMesg->jump_list[i].type;
  }
  
  return slkMesg;
}




/************************************************/


/* converts an IntScaffoldMesg to a SnapScaffoldMessage.
   What happend ?
   - The contigs list is traversed and for each IntContigPair
   a SnapContigPair is created with the internal contig IDs replaced
   by external IDs
   - (07/01/00) in addition we store now the assigned UIDs in a VA
   in order to be able to dump them in the map file.
*/
static SnapScaffoldMesg* convert_ISF_to_SCF(IntScaffoldMesg* isfMesg, int32 real)
{
  CDS_UID_t uid;
  CDS_UID_t *di;
  int *used;
  int i;
  int true  = TRUE;

  SnapScaffoldMesg *scfMesg = (SnapScaffoldMesg*) safe_malloc(sizeof(SnapScaffoldMesg));

  /* NEW : the ISF messages HAVE now an internal ID
     so we can now check for the IID to UID mapping */

  di = GetCDS_UID_t(ISFmap,isfMesg->iaccession);
  if ((di != NULL) && (*di != 0)) {
    sprintf(errorreport,"ISF internal id %d spotted second time", isfMesg->iaccession);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }

  uid = getUID(real);
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
	      sprintf(errorreport,"ICM1inISFmap is not initialized for contig ID %d",con1);
	      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  if( *used == TRUE )
	    {
	      sprintf(errorreport,"Contig ID %d was already used as ct1 in a scaffold",con1);
	      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  else
	    Setint(ICM1inISFmap,con1,&true);

	  // try to find out whether we used the con2 before
	  used = Getint(ICM2inISFmap,con2);
	  if( used == NULL )
	    {
	      sprintf(errorreport,"ICM2inISFmap is not initialized for contig ID %d",con2);
	      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  if( *used == TRUE )
	    {
	      sprintf(errorreport,"Contig ID %d was already used as ct2 in a scaffold",con2);
	      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	    }  
	  else
	    Setint(ICM2inISFmap,con2,&true);



	  di = GetCDS_UID_t(ICMmap,con1);
	  if ((di == NULL) || (*di == 0)) {
            sprintf(errorreport,"ISF reference before definition for contig ID %d",con1);
            error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
          }  
	  scfMesg->contig_pairs[i].econtig1 = *di;
	  
	  di = GetCDS_UID_t(ICMmap,con2);
	  if ((di == NULL) || (*di == 0)) {
            sprintf(errorreport,"ISF reference before definition for contig ID %d",con2);
            error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
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
	    sprintf(errorreport,"ICM1inISFmap is not initialized for contig ID %d",con1);
	    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	  }  
	if( *used == TRUE )
	  {
	    sprintf(errorreport,"Contig ID %d was already used as ct1 in a scaffold",con1);
	    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	  }  
	else
	  Setint(ICM1inISFmap,con1,&true);	

	di = GetCDS_UID_t(ICMmap,con1);
	if ((di == NULL) || (*di == 0)) {
          sprintf(errorreport,"ISF reference before definition for contig ID %d",con1);
          error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
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


static SnapMateDistMesg* convert_IMD_to_MDI(IntMateDistMesg* imdMesg, int32 real)
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
      sprintf(errorreport,"IMD reference before definition error for ID %d",imdMesg->refines);
      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
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


static SnapDegenerateScaffoldMesg* convert_IDS_to_DSC(IntDegenerateScaffoldMesg* idsMesg, int32 real)
     /* converts an  IntDegenerateScaffoldMesg to a  SnapDegenerateScaffoldMesg
	What happend ?
	- The icontig is replaced by the corresponding UID
	- the degenerate scaffold is assigned a new UID
     */
{
  CDS_UID_t uid;
  CDS_UID_t *di;

  SnapDegenerateScaffoldMesg* dscMesg = (SnapDegenerateScaffoldMesg*) safe_malloc(sizeof(SnapDegenerateScaffoldMesg));
  
  /* check whether the internal id is defined */
  di = GetCDS_UID_t(ICMmap,idsMesg->icontig);
  if ((di == NULL) || (*di == 0)) {
    sprintf(errorreport,"Could not find Contig iid %d in ICMmap", idsMesg->icontig);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }

  dscMesg->econtig = *di;

  /* This is a new message. Hence we get a new UID */

  uid = getUID(real);
  dscMesg->eaccession = uid;
  SetCDS_UID_t(DSCmap,idsMesg->icontig,&uid);

  return dscMesg;
}




/************************************************/


static AugFragMesg* convert_IAF_to_AFG(IntAugFragMesg* iafMesg, int32 real)
     /* converts an IntAugFragMesg to a AugFragMessage.
	What happend ?
	- The internal frgament ID is replaced with the external UID
	- If the clear range in the IAF was corrected it is taken
	otherwise it is taken from the ClearStartMap and ClearEndMap
     */
{
  CDS_UID_t *di;
  uint32 *sdi;
  uint32 *size;
  uint32 dummyZero=0;

  AugFragMesg *afgMesg = (AugFragMesg*) safe_malloc(sizeof(AugFragMesg));
  
  di = fetch_UID_from_fragStore(iafMesg->iaccession);
  if( di == NULL ){
    sprintf(errorreport,"No fragment with ID %d in the frag store",iafMesg->iaccession);
    error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
  }

  afgMesg->eaccession = *di;
  afgMesg->iaccession = iafMesg->iaccession;
  afgMesg->chimeric    = iafMesg->chimeric;
  afgMesg->mate_status = iafMesg->mate_status;
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
      sdi = fetch_clearRange_from_fragStore(ClearStartMap,iafMesg->iaccession);
      if( sdi == NULL )
	{
	  sprintf(errorreport,"No SeqInterval associated with ID %d in the frag store",iafMesg->iaccession);
	  error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}
      afgMesg->clear_rng.bgn = *sdi;

      sdi = fetch_clearRange_from_fragStore(ClearEndMap,iafMesg->iaccession);
      if( sdi == NULL )
	{
	  sprintf(errorreport,"No SeqInterval associated with ID %d in the frag store",iafMesg->iaccession);
	  error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
	}
      afgMesg->clear_rng.end = *sdi;
    }
    
  return afgMesg;
}






/*****************/
/* free routines */
/*****************/


static void free_DSC(SnapDegenerateScaffoldMesg* dscMesg){
  safe_free(dscMesg);
}

static void free_UTG(SnapUnitigMesg* utgMesg){
  /*** DO NOT FREE DELTA ARRAYS...THESE WERE COPIED BY REFERENCE */
  if( utgMesg->f_list != NULL )
    safe_free(utgMesg->f_list);
  if( utgMesg->consensus != NULL )
    safe_free(utgMesg->consensus);
  if( utgMesg->quality != NULL )
    safe_free(utgMesg->quality);
  safe_free(utgMesg);
}

static void free_ULK(SnapUnitigLinkMesg* ulkMesg){
  if( ulkMesg->jump_list != NULL )
    safe_free(ulkMesg->jump_list);
  safe_free(ulkMesg);
}

static void free_CCO(SnapConConMesg* ccoMesg){
  /*** DO NOT FREE DELTA ARRAYS...THESE WERE COPIED BY REFERENCE */
  if( ccoMesg->num_vars > 0)
    safe_free(ccoMesg->vars);
  if( ccoMesg->num_pieces > 0)
    safe_free(ccoMesg->pieces);
  if( ccoMesg->num_unitigs > 0)
    safe_free(ccoMesg->unitigs);
  if( ccoMesg->consensus != NULL )
    safe_free(ccoMesg->consensus);
  if( ccoMesg->quality != NULL )
    safe_free(ccoMesg->quality);
  safe_free(ccoMesg);
}

static void free_CLK(SnapContigLinkMesg* clkMesg){
  if( clkMesg->jump_list != NULL )
    safe_free(clkMesg->jump_list);
  safe_free(clkMesg);
}

static void free_SLK(SnapScaffoldLinkMesg* slkMesg){
  if( slkMesg->jump_list != NULL )
    safe_free(slkMesg->jump_list);
  safe_free(slkMesg);
}

static void free_SCF(SnapScaffoldMesg* scfMesg){
  if( scfMesg->contig_pairs != NULL )
    safe_free(scfMesg->contig_pairs);
  safe_free(scfMesg);
}

static void free_MDI(SnapMateDistMesg* mdiMesg){
  if( mdiMesg->histogram != NULL )
    safe_free(mdiMesg->histogram);
  safe_free(mdiMesg);
}

static void free_AFG(AugFragMesg* afgMesg){
  safe_free(afgMesg);
}



/********************/
/* read the stores  */
/********************/

static void read_stores(char* fragStoreName);

/************************/
/* main output routine  */
/************************/


void output_snapshot(char* fragStoreName,
		     char** inputFileList, int32 numInputFiles,
		     char* outputFileName, char* mapFileName,
		     int32 real, int32 quiet,
		     int32 random, CDS_UID_t uidStart, 
		     int argc, char *argv[])
{
  GenericMesg *pmesg       = NULL; 
  FILE        *fileInput   = NULL;
  FILE        *fileOutput  = NULL;
  int32        ifile;
  char        *inputFileName;
  int numIAF = 0;
  int numIUM = 0;
  int numICM = 0;

  /* Test preconditions */
  if( outputFileName == NULL )
    {
      sprintf(errorreport,"Argument outputFileName is NULL");
      error(errorreport, AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }

  if( strlen(outputFileName) == 0 )
    {
      sprintf(errorreport,"Argument outputFileName is empty");
      error(errorreport, AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }
  /* Test preconditions */
  if( mapFileName == NULL )
    {
      sprintf(errorreport,"Argument mapFileName is NULL");
      error(errorreport, AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
    }

  if( strlen(mapFileName) == 0 )
    {
      sprintf(errorreport,"Argument mapFileName is empty");
      error(errorreport, AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
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
  ClearStartMap = CreateVA_uint32(ARRAYSIZE);
  ClearEndMap   = CreateVA_uint32(ARRAYSIZE);  

  /* set the UID start to the number given by uidStart */
  set_start_uid(uidStart);

  /* set the simulator bool to the value of quiet (-Q flag) */
  simulator = quiet;

  FSHandle     = openGateKeeperStore(fragStoreName, FALSE);

  if( random == FALSE )
    read_stores(fragStoreName);
  

  if( strcmp(outputFileName,"-") == 0 ){
    PipeOut = TRUE;
    fprintf(stderr,"*** Terminator : Print output to stdout ***\n");
  }
    
  if( ! PipeOut )
    fileOutput  = file_open(outputFileName,"w");  
  else
    fileOutput = stdout;


  initUID(real);


  for(ifile = 0; ifile < numInputFiles; ifile++){
    inputFileName = inputFileList[ifile];


    /* Test preconditions */
    if( inputFileName == NULL )
      {
	sprintf(errorreport,"Argument inputFileName is NULL");
	error(errorreport, AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
      }

    if( strlen(inputFileName) == 0 )
      {
	sprintf(errorreport,"Argument inputFileName is empty");
	error(errorreport, AS_TER_EXIT_FAILURE,__FILE__,__LINE__);
      }

    if( strcmp(inputFileName,"-") == 0 ){
      PipeIn = TRUE;
      fprintf(stderr,"*** Terminator : Read input from stdin ***\n");
    }

    if( ! PipeIn )
      fileInput = file_open(inputFileName,"r");  
    else
      fileInput = stdin;


    /* main loop  */
    fprintf(stderr,"*** Terminator : Reading cgw file %d %s ***\n", ifile, inputFileName);

    while(ReadProtoMesg_AS(fileInput,&pmesg) != EOF){
      switch(pmesg->t){
      case MESG_ADT : {
	AuditMesg *adt_mesg;
	adt_mesg = pmesg->m;
	pmesg->t = MESG_ADT;
      
	VersionStampADT(adt_mesg,argc,argv);
	WriteProtoMesg_AS(fileOutput,pmesg);    
	WriteProtoMesg_AS(stderr,pmesg);    
      }
      break; 
      case MESG_IAF :{
	IntAugFragMesg *iafMesg = (IntAugFragMesg*) pmesg->m;
	AugFragMesg *afgMesg = convert_IAF_to_AFG(iafMesg,real);
	pmesg->m = afgMesg;
	pmesg->t = MESG_AFG;
	WriteProtoMesg_AS(fileOutput,pmesg);
	free_AFG(afgMesg);
	if((++numIAF % 1000000) == 0){
	  fprintf(stderr,"* Processed %d frags\n", numIAF);
	  fflush(stderr);
	}
      }
      break;
      case MESG_IUM :{
	IntUnitigMesg  *iumMesg = (IntUnitigMesg*) pmesg->m;
	SnapUnitigMesg *utgMesg = convert_IUM_to_UTG(iumMesg,real);
	pmesg->m = utgMesg;
	pmesg->t = MESG_UTG;
	WriteProtoMesg_AS(fileOutput,pmesg);
	free_UTG(utgMesg);
	if((++numIUM % 1000000) == 0){
	  fprintf(stderr,"* Processed %d unitigs\n", numIUM);
	  fflush(stderr);
	}
      }
      break;
      case MESG_IUL : {
	IntUnitigLinkMesg  *iulMesg = (IntUnitigLinkMesg*) pmesg->m;
	SnapUnitigLinkMesg *ulkMesg = convert_IUL_to_ULK(iulMesg,real);
	pmesg->m = ulkMesg;
	pmesg->t = MESG_ULK;
	WriteProtoMesg_AS(fileOutput,pmesg);
	free_ULK(ulkMesg);   
      }
      break;
      case MESG_ICM : {
	IntConConMesg *icmMesg = (IntConConMesg*) pmesg->m;
	SnapConConMesg *ccoMesg = convert_ICM_to_CCO(icmMesg,real);
	pmesg->m = ccoMesg;
	pmesg->t = MESG_CCO;
	WriteProtoMesg_AS(fileOutput,pmesg);
	free_CCO(ccoMesg);    
	if((++numICM % 1000000) == 0){
	  fprintf(stderr,"* Processed %d contigs\n", numICM);
	  fflush(stderr);
	}
      }
      break;
      case MESG_ICL : {
	IntContigLinkMesg *iclMesg = (IntContigLinkMesg*) pmesg->m;
	SnapContigLinkMesg *clkMesg = convert_ICL_to_CLK(iclMesg,real);
	pmesg->m = clkMesg;
	pmesg->t = MESG_CLK;
	WriteProtoMesg_AS(fileOutput,pmesg);
	free_CLK(clkMesg);    
      }
      break;
      case MESG_ISL : {
	InternalScaffoldLinkMesg *islMesg = (InternalScaffoldLinkMesg*) pmesg->m;
	SnapScaffoldLinkMesg *slkMesg = convert_ISL_to_SLK(islMesg,real);
	pmesg->m = slkMesg;
	pmesg->t = MESG_SLK;
	WriteProtoMesg_AS(fileOutput,pmesg);
	free_SLK(slkMesg);    
      }
      break;
      case MESG_ISF : {
	IntScaffoldMesg *isfMesg = (IntScaffoldMesg*) pmesg->m;
	SnapScaffoldMesg *scfMesg = convert_ISF_to_SCF(isfMesg,real);
	pmesg->m = scfMesg;
	pmesg->t = MESG_SCF;
	WriteProtoMesg_AS(fileOutput,pmesg);
	free_SCF(scfMesg);    
      }
      break;
      case MESG_IMD : {
	IntMateDistMesg *imdMesg = (IntMateDistMesg*) pmesg->m;
	SnapMateDistMesg *mdiMesg = convert_IMD_to_MDI(imdMesg,real);
	pmesg->m = mdiMesg;
	pmesg->t = MESG_MDI;
	WriteProtoMesg_AS(fileOutput,pmesg);
	free_MDI(mdiMesg);    
      }
      break;
      case MESG_IDS : {
	IntDegenerateScaffoldMesg *idsMesg = (IntDegenerateScaffoldMesg*) pmesg->m;
	SnapDegenerateScaffoldMesg *dscMesg = convert_IDS_to_DSC(idsMesg,real);
	pmesg->m = dscMesg;
	pmesg->t = MESG_DSC;
	WriteProtoMesg_AS(fileOutput,pmesg);
	free_DSC(dscMesg);    
      }
      break;
      default:
	{
#if DEBUG > 0
	  fprintf(stderr,"Swallowing message of type %d\n",pmesg->t); 
	  WriteProtoMesg_AS(stderr,pmesg);
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


  {
    FILE    *F = NULL;
    char     N[1024];

    sprintf(N, "%s.fragment.iidtouid", mapFileName);
    F = file_open(N, "w");  
    fprintf(F,"Fragment IID2UID map\n");
    DumpIID2UIDmap(FRGmap,F);
    fclose(F);

    sprintf(N, "%s.unitig.iidtouid", mapFileName);
    F = file_open(N, "w");  
    fprintf(F,"Unitig IID2UID map\n");
    DumpIID2UIDmap(IUMmap,F);
    fclose(F);

    sprintf(N, "%s.contig.iidtouid", mapFileName);
    F = file_open(N, "w");  
    fprintf(F,"Contig IID2UID map\n");
    DumpIID2UIDmap(ICMmap,F);
    fclose(F);

    sprintf(N, "%s.scaffold.iidtouid", mapFileName);
    F = file_open(N, "w");  
    fprintf(F,"Scaffold IID2UID map\n");
    DumpIID2UIDmap(ISFmap,F);
    fclose(F);

    sprintf(N, "%s.distrib.iidtouid", mapFileName);
    F = file_open(N, "w");  
    fprintf(F,"Distrib IID2UID map\n");
    DumpIID2UIDmap(DSTmap,F);
    fclose(F);

    sprintf(N, "%s.degeneratecontig.iidtouid", mapFileName);
    F = file_open(N, "w");  
    fprintf(F,"Degenerate Contig IID to Scaffold UID map\n");
    DumpIID2UIDmap(DSCmap,F);
    fclose(F);
  }


  closeGateKeeperStore(FSHandle);
  

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
  DeleteVA_uint32(ClearStartMap);
  DeleteVA_uint32(ClearEndMap);
}






/* read information from the fragment, gatekeeper store */


void read_stores(char* fragStoreName)
{
  fragRecord *input;
  FragStream *streamHandle = 0;
  CDS_UID_t *di;

  fprintf(stderr,"*** Terminator : Reading Fragment Store ***\n");

  streamHandle = openFragStream(FSHandle, FRAG_S_INF);

  input =  new_fragRecord();
  while( nextFragStream(streamHandle,input) ){
    CDS_IID_t iid;
    CDS_UID_t uid;
    uint32 cStart,cEnd;
    uint32         noOfMatches;
    short truedummy  = TRUE;

    /* read the iid and uid pair */
    iid = getFragRecordIID(input);
    uid = getFragRecordUID(input);

    if(((iid+1) % 1000000) == 0){
      fprintf(stderr,"* Read frag %d from fragStore\n", iid);
      fflush(stderr);
    }
    /* test whether the iid is not already present */
    if (test_frg_present(iid)) {
      sprintf(errorreport,"Internal fragment ID %d occurred twice (FRGmap)",iid);
      error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
    }
    SetCDS_UID_t(FRGmap,iid,&uid);
    Setshort(FRGpresent,iid,&truedummy);

    //  BPW - before 21nov2005, we used to always set_start_uid(uid+1)
    //  on the first time through here.  What this did is to just
    //  reset the start uid specified on the command line to be the
    //  highest uid in the .asm + 1.  So, we don't do that anymore.
    //
    if(get_start_uid() <= uid)
      set_start_uid(uid+1);
    
#if DEBUG_UID
    fprintf(stderr,"SYS_UID_uidStart after reading frags = %lu\n",get_start_uid());
#endif

    /* get the clear region from the frag store */

    cStart = getFragRecordClearRegionBegin(input, AS_READ_CLEAR_LATEST);
    cEnd   = getFragRecordClearRegionEnd  (input, AS_READ_CLEAR_LATEST);

    Setuint32(ClearStartMap,iid,&cStart);
    Setuint32(ClearEndMap,iid,&cEnd);
  }
  

  closeFragStream(streamHandle);

  /* open the gatekeeper store */
  fprintf(stderr,"*** Terminator : Reading Gatekeeper Store ***\n");
  {
    /* reading distributions */
    {
     GateKeeperLibraryRecord gkpl;
     StoreStat stat;
     int i ;
     statsStore(FSHandle->lib, &stat);
#if DEBUG > 0
     fprintf(stderr,"* Stats for Dist Store are first:%lu last :%lu\n",
	     stat.firstElem, stat.lastElem);
#endif
     for(i = stat.firstElem; i <= stat.lastElem; i++){
       getGateKeeperLibraryStore(FSHandle->lib,i,&gkpl);
#if DEBUG > 1
       fprintf(stderr,"* Dist %d UID:%lu del:%d red:%d mean:%f std:%f batch(%d,%d) prevID:%d prevInstanceID:%d\n",
	       i,gkpl.UID, gkpl.deleted, gkpl.redefined, gkpl.mean, gkpl.stddev,
	       gkpl.birthBatch, gkpl.deathBatch, gkpl.prevID, gkpl.prevInstanceID);
#endif

       di = GetCDS_UID_t(DSTmap,i);
       if ((di != NULL) && (*di != 0)) {
	 sprintf(errorreport,"Internal DST ID %d occurred twice",i);
	 error(errorreport,AS_TER_EXIT_FAILURE,__FILE__,__LINE__); 
       }
       SetCDS_UID_t(DSTmap,i,&gkpl.UID);
       if( get_start_uid() <= gkpl.UID )
         set_start_uid(gkpl.UID+1);

#if DEBUG_UID
       fprintf(stderr,"SYS_UID_uidStart after reading distribs %lu\n",get_start_uid());
#endif
     }
    }
  }
}


