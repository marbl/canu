
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
/*********************************************************************
 * $Id: AS_CGB_miniunitigger.c,v 1.3 2005-03-22 19:02:16 jason_miller Exp $
 * Module:  AS_CGB_miniunitigger
 * Description: 
 * Assumptions:
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "AS_CGB_unitigger_globals.h"
#include "AS_CGB_miniunitigger.h"

#if 0
#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)
#define FALSE 0
#undef DEBUGGING1
#undef DEBUGGING2
#endif

/*************************************************************************/
static char CM_ID[] = "$Id: AS_CGB_miniunitigger.c,v 1.3 2005-03-22 19:02:16 jason_miller Exp $";
/*************************************************************************/

extern int REAPER_VALIDATION;

typedef struct { 

  UnitiggerGlobals unitiggerGlobals;
  TStateGlobals    stateGlobals;
  THeapGlobals     heapGlobals;
  
#if 0  /* Dynamic binding */

  int             yourDataTypeExtended;
  YourDataType   *temp_data;
  YourDataType   *scan_data;
  YourDataType   *aggr_data;
  YourDataType   *sample_data; /* sample_data array */
  YourDataType   *bucket_data;
  YourataType   *(*indexdata)(YourDataType *b,int ib);
  void       (*setdata)(YourDataType *a,int ib,YourDataType *b);
  void       (*aggregate)(YourDataType *a,int ib,YourDataType *b);
  void       (*printdata)(FILE *fout,YourDataType *,YourDataType *,YourDataType *);
#endif /* Dynamic binding */
} MiniUnitiggerStruct;


// set default or uninitialized values for important variables
static void InitializeGlobals
(
  UnitiggerGlobals * rg,
  char * program_name
  )
{
  memset( rg, 0, sizeof( UnitiggerGlobals ) );
#ifndef USE_STATIC_CMD_STRINGS
  rg->program_name = program_name;
#else // USE_STATIC_CMD_TRINGS
  strcpy(rg->program_name,program_name);
#endif //USE_STATIC_CMD_STRINGS
  rg->input_store_flag = FALSE;
  rg->output_store_flag = FALSE;
  rg->create_store_flag = FALSE;
  rg->append_store_flag = FALSE;
  rg->clobber_store_flag = FALSE;
  rg->developer_mode_flag = FALSE;
  rg->num_threads = 1;
  rg->work_limit_per_candidate_edge = 1000;
  rg->remove_blizzard_overlaps = FALSE;
  rg->use_all_overlaps_in_reaper_pass = TRUE;
  rg->as_cgb_max_frag_iid = 100000 * CGB_MULTIPLIER;
  rg->dvt_double_sided_threshold_fragment_end_degree = 100;
  rg->con_double_sided_threshold_fragment_end_degree = 100;
  rg->cutoff_fragment_end_degree = 1000000;
#ifdef REPAIR_BREAKERS
  rg->breaker_fix = 3;
#endif // REPAIR_BREAKERS
  rg->dechord_the_graph = TRUE;
  rg->check_point_level = 0;
  rg->compress_the_graph = TRUE;
  rg->cgb_unique_cutoff = CGB_UNIQUE_CUTOFF;
  rg->walk_depth = 100;
  rg->iv_start = 0;
  rg->overlap_error_threshold = PerMil_to_CGB_ERATE_TYPE(1000);
  rg->work_limit_placing_contained_fragments = 20;
  rg->walk_depth=100;
  rg->iv_start=0;
  rg->num_threads = 4;
  rg->output_iterations_flag = TRUE;
  rg->aggressive_spur_fragment_marking = TRUE;
}

MiniUnitiggerObject *createMiniUnitigger
(
 int maxfrags,
 int maxedges,
 int maxtext
)
{ 
  MiniUnitiggerStruct * self 
    = (MiniUnitiggerStruct *) malloc(sizeof(MiniUnitiggerStruct));
  assert(NULL != self);
  memset( self, 0, sizeof(MiniUnitiggerStruct));

  {  
    UnitiggerGlobals    * rg     = &(self->unitiggerGlobals);
    TStateGlobals       * gstate = &(self->stateGlobals);
    THeapGlobals        * heapva = &(self->heapGlobals);

    InitializeGlobals( rg, "miniunitigger");
    
    // Fill in default values for MiniUnitiggerStruct here:
    rg->Input_Graph_Store = NULL;
    rg->maxfrags = maxfrags;
    rg->maxedges = maxedges;
    rg->maxtext  = maxtext;
    rg->dvt_double_sided_threshold_fragment_end_degree = 10;
    rg->con_double_sided_threshold_fragment_end_degree = 10;
    rg->cutoff_fragment_end_degree = 1000000;
    
    open_fgb_store( gstate, heapva );
    
    /* Create globals */ {
      // read the fgb store to initialize the heap.
      time_t tp1 = time(NULL);
      time_t tp2;
      
      fprintf(stderr,"Begin reading the fragment graph store.\n");
      system_date();
      read_fgb_store(
		     rg->Input_Graph_Store,
		     gstate, heapva,
		     rg->maxfrags, 
		     rg->maxedges, 
		     rg->maxtext);
      tp2 = time(NULL); 
      fprintf(stderr,"%10" F_TIME_TP " sec: Finished reading the fragment graph store.\n",
	      (tp2-tp1));
      system_date();
    }
  }
  // Allocate memory for work space here.

  // NULL the extension functions here.
  
  return (MiniUnitiggerObject *)self;
}

int destroyMiniUnitigger
(MiniUnitiggerObject * _self)
{
  int error = 0; // the succesful return value.
  MiniUnitiggerStruct * self = (MiniUnitiggerStruct *) _self;
  assert(NULL != self);
  {
    // UnitiggerGlobals    * rg     = &(self->unitiggerGlobals);
    TStateGlobals       * gstate = &(self->stateGlobals);
    THeapGlobals        * heapva = &(self->heapGlobals);
    
    // Free the allocated data.
    close_fgb_store( gstate, heapva );
  }
  free(self);
  return error;
}

int set_cgb_unique_cutoff_MiniUnitigger
( MiniUnitiggerObject * _self, float cgb_unique_cutoff)
  // A typical value used for the Human genome is five.
{
  int error = 0;
  MiniUnitiggerStruct * self = (MiniUnitiggerStruct *) _self;
  self->unitiggerGlobals.cgb_unique_cutoff = cgb_unique_cutoff;
  return error; 
}

int set_overlap_error_threshold_MiniUnitigger
( MiniUnitiggerObject * _self, float overlap_error_threshold)
// An overlap_error_threshold of 0.6 means six percent.  This value is
// used to filter the overlaps on input.  To turn it off, set the
// value to one.
{ int error = 0;
  MiniUnitiggerStruct * self = (MiniUnitiggerStruct *) _self;
  self->unitiggerGlobals.overlap_error_threshold = overlap_error_threshold;
 return error; 
}

int set_as_cgb_max_frag_iid_MiniUnitigger
( MiniUnitiggerObject * _self, int as_cgb_max_frag_iid)
// Unfortunately I do not have a hash function that maps the fragment
// IIDs in the OFG messages to a dense interval of integers. Instead I
// just allocate an array of with range [0 .. as_cgb_max_frag_iid] to
// store the mapping.
{
  int error = 0;
  MiniUnitiggerStruct * self = (MiniUnitiggerStruct *) _self;
  self->unitiggerGlobals.as_cgb_max_frag_iid = as_cgb_max_frag_iid;
  return error;
} 

#if 0 // Dynamic binding
extern int extend_MiniUnitigger_with_YourDataType
(  
 MiniUnitiggerObject *self,
 size_t sizeof_YourDataType,
 YourDataType * (*indexdata)(YourDataType *b,int ib),
 void (*setdata)(YourDataType *a,int ib,YourDataType *b),
 void (*aggregate)(YourDataType *a,int ib,YourDataType *b),
 void (*printdata)(FILE *fout,
                   YourDataType *,
                   YourDataType *,
                   YourDataType *)
 )
// If later we want to dynamicly bind a function inside the Mini
// Unitigger, then I propose to use a binding method similar to that
// used by the AS_CGB histogramming utility.
{
  int error = 0;
  return error;
} 
#endif  // Dynamic binding

// typedef struct { } RunMiniUnitiggerParms  is in AS_CGB_miniunitigger.
// typedef struct { } RunMiniUnitiggerResults  is in AS_CGB_miniunitigger.

int run_MiniUnitigger
(
 MiniUnitiggerObject     * _self,
 RunMiniUnitiggerParams  * params,
 RunMiniUnitiggerResults * results
 ) {
  int error = 0;
  int status = 0;
  int argc = 1;
  char * program_name = "runMiniUnitigger";
  char *argv[1] = { program_name };
  
  MiniUnitiggerStruct * self = (MiniUnitiggerStruct *) _self;
  UnitiggerGlobals    * rg   = &(self->unitiggerGlobals);
  TStateGlobals       * gstate = &(self->stateGlobals);
  THeapGlobals        * heapva = &(self->heapGlobals);

  assert(NULL == heapva->the_ofg_messages);
  assert(NULL == heapva->the_ovl_messages);
  assert(NULL == heapva->the_ofg_source);
  assert(NULL == heapva->the_ovl_source);
  assert(NULL == heapva->the_imp_messages);
  assert(NULL == heapva->the_ium_messages); 
  assert(NULL == heapva->the_imp_source);
  assert(NULL == heapva->the_ium_source);
  //assert(NULL == heapva->frag_annotations);
  //assert(NULL == heapva->chunksrcs);

  heapva->the_ofg_messages = params->the_ofg_messages;
  heapva->the_ovl_messages = params->the_ovl_messages;
  heapva->the_ofg_source   = params->the_ofg_source;
  heapva->the_ovl_source   = params->the_ovl_source;

  heapva->the_imp_messages = results->the_imp_messages;
  heapva->the_ium_messages = results->the_ium_messages;
  heapva->the_imp_source   = results->the_imp_source;
  heapva->the_ium_source   = results->the_ium_source;
  //heapva->frag_annotations = results->the_imp_source;
  //heapva->chunksrcs = results->the_ium_source;

  assert(NULL != heapva->the_ofg_messages);
  assert(NULL != heapva->the_ovl_messages);
  assert(NULL != heapva->the_ofg_source);
  assert(NULL != heapva->the_ovl_source);
  assert(NULL != heapva->the_imp_messages);
  assert(NULL != heapva->the_ium_messages);
  assert(NULL != heapva->the_imp_source);
  assert(NULL != heapva->the_ium_source);
  
  clear_fgb_store( gstate, heapva, TRUE, TRUE, TRUE, TRUE);
  
  status = main_fgb( argc, argv, gstate, heapva, rg);
  fprintf(stderr,"main_fgb status = %d\n", status);

  status = main_cgb( argc, argv, gstate, heapva, rg);
  fprintf(stderr,"main_cgb status = %d\n", status);

#if 0
  view_fgb_chkpnt( "deleteme", heapva->frags, heapva->edges);
#endif
  
  convert_the_chunks_to_IUM
    (/* Input Only*/
     heapva->frags,
     heapva->edges,
     heapva->frag_annotations,
     heapva->chunkfrags,
     heapva->thechunks,
     heapva->chunkseqs,
     heapva->chunkquas,
     heapva->chunksrcs,
     rg->analysis_level,
     gstate->global_fragment_arrival_rate,
     
     /* Output Only */
     results->the_imp_messages,
     results->the_ium_messages,
     results->the_imp_source,
     results->the_ium_source
     );

  // Reset to run again!  We assume that the client owns the VAs in
  // the params and results structures.
  
  heapva->the_ofg_messages = NULL;
  heapva->the_ovl_messages = NULL;
  heapva->the_ofg_source   = NULL;
  heapva->the_ovl_source   = NULL;
  heapva->the_imp_messages = NULL;
  heapva->the_ium_messages = NULL;
  heapva->the_imp_source   = NULL;
  heapva->the_ium_source   = NULL;

  assert(NULL == heapva->the_ofg_messages);
  assert(NULL == heapva->the_ovl_messages);
  assert(NULL == heapva->the_ofg_source);
  assert(NULL == heapva->the_ovl_source);
  assert(NULL == heapva->the_imp_messages);
  assert(NULL == heapva->the_ium_messages);
  assert(NULL == heapva->the_imp_source);
  assert(NULL == heapva->the_ium_source);

  return error;
}



