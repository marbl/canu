
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
 *
 * Module: AS_CGB_store.c
 * Description: Reads and writes the check point data.
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
//#include <dirent.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_CGB_all.h"

//  Define this to disable reading/writing the store from/to disk.
//
#define DONT_USE_STORE

/*************************************************************************/
static char CM_ID[] 
= "$Id: AS_CGB_store.c,v 1.6 2006-11-14 17:52:14 eliv Exp $";
/*************************************************************************/

void open_fgb_store
(
 TStateGlobals    * gstate,
 THeapGlobals     * heapva
 //UnitiggerGlobals * rg
 )
{
  assert(NULL != heapva);

  assert(NULL == heapva->frags);
  assert(NULL == heapva->edges);
  assert(NULL == heapva->next_edge_obj);
  assert(NULL == heapva->chunkfrags);
  assert(NULL == heapva->thechunks);
  assert(NULL == heapva->chunkseqs);
  assert(NULL == heapva->chunkquas); 
  assert(NULL == heapva->frag_annotations);
  assert(NULL == heapva->chunksrcs);
  
  assert(NULL == heapva->the_ofg_messages);
  assert(NULL == heapva->the_ovl_messages);
  assert(NULL == heapva->the_ofg_source);
  assert(NULL == heapva->the_ovl_source);
  assert(NULL == heapva->the_imp_messages);
  assert(NULL == heapva->the_ium_messages);
  assert(NULL == heapva->the_imp_source);
  assert(NULL == heapva->the_ium_source);

}

#define SAFE_Delete_VA(x) Delete_VA(x); x = NULL;

void close_fgb_store
(
 TStateGlobals    * gstate,
 THeapGlobals     * heapva
 )
{
  assert(NULL != heapva);

  // ierr = closeFragStore(TheFragStore); assert(ierr == 0);
  
  fprintf(stderr,"Unitigger memory clean up\n");
  
  SAFE_Delete_VA(heapva->frags);
  SAFE_Delete_VA(heapva->edges);
  SAFE_Delete_VA(heapva->next_edge_obj);
  SAFE_Delete_VA(heapva->chunkfrags);
  SAFE_Delete_VA(heapva->thechunks);
  SAFE_Delete_VA(heapva->chunkseqs);
  SAFE_Delete_VA(heapva->chunkquas); 
  SAFE_Delete_VA(heapva->frag_annotations);
  SAFE_Delete_VA(heapva->chunksrcs);
  
  assert(NULL == heapva->the_ofg_messages);
  assert(NULL == heapva->the_ovl_messages);
  assert(NULL == heapva->the_ofg_source);
  assert(NULL == heapva->the_ovl_source);
  assert(NULL == heapva->the_imp_messages);
  assert(NULL == heapva->the_ium_messages);
  assert(NULL == heapva->the_imp_source);
  assert(NULL == heapva->the_ium_source);
}

void read_fgb_store
( const char * const theStorePath,
  TStateGlobals * gstate,
  THeapGlobals  * heapva,
  const size_t new_additional_number_of_frags,
  const size_t new_additional_number_of_edges,
  const size_t new_additional_amount_of_text)
{
  FILE *fp = NULL;
  int my_createmode = (NULL == theStorePath) || ('\0' == theStorePath[0]);
  int ierr;

  assert(NULL != gstate);
  assert(NULL != heapva);
  fprintf(stderr,"* Reading the store\n");
  
  //fprintf(stderr, "* theStorePath = %p\n", theStorePath);
  //fprintf(stderr, "* theStorePath [0] = %d\n", (int) (theStorePath [0]));

#ifdef DONT_USE_STORE
  if(! my_createmode){
    fprintf(stderr, "ERROR: Someone wanted to read the store '%s', but DONT_USE_STORE is defined.\n",
            theStorePath);
    assert(0);
  }
#endif

  if(! my_createmode){
    char buffer[FILENAME_MAX];
    fprintf(stderr, 
	    "* Reading the fragment graph store %s\n", theStorePath);
    sprintf(buffer,"%s/fgb.ckp", theStorePath);
    assert(NULL == fp);
    fp = fopen(buffer,"r");
    if(fp == NULL) {
      fprintf
	(stderr,"* The graph store check point file <%s> was not found.\n",
	 buffer);
      exit(1);
    }
  }

  if(! my_createmode) {
    size_t nitems = 0;
    fprintf(stderr,"Read globals...\n");
    nitems = fread(gstate,sizeof(TStateGlobals),1,fp);
    if(nitems != 1) {
      fprintf
	(stderr,
	 "* The check point file <%s> appears corrupted.\n",
	 theStorePath);
      exit(1);
    }
    assert(nitems == 1);
    
    assert(NULL == heapva->frags);
    heapva->frags = CreateFromFileVA_Afragment
      (fp,new_additional_number_of_frags);
    assert(NULL != heapva->frags);
    fprintf(stderr,"Finished reading frags\n");

    assert(NULL == heapva->edges);
    heapva->edges = CreateFromFileVA_Aedge
      (fp,new_additional_number_of_edges);
    assert(NULL != heapva->edges);
    fprintf(stderr,"Finished read edges\n");

    assert(NULL == heapva->frag_annotations);
    heapva->frag_annotations = CreateFromFileVA_char
      (fp,new_additional_amount_of_text);
    assert(NULL != heapva->frag_annotations);
    fprintf(stderr,"Finished reading annotations\n");

    assert(NULL == heapva->next_edge_obj);
    heapva->next_edge_obj = CreateFromFileVA_IntEdge_ID
      (fp,new_additional_number_of_edges);
    assert(NULL != heapva->next_edge_obj);
    fprintf(stderr,"Finished reading next_edge_obj...\n");

#if 1
    assert(NULL == heapva->chunkfrags);
    heapva->chunkfrags = CreateFromFileVA_AChunkFrag(fp,0);
    assert(NULL != heapva->chunkfrags);
    fprintf(stderr,"Finished reading chunkfrags...\n");

    assert(NULL == heapva->thechunks);
    heapva->thechunks = CreateFromFileVA_AChunkMesg(fp,0);
    assert(NULL != heapva->thechunks);
    fprintf(stderr,"Finished reading thechunks...\n");
#else
    assert(NULL == heapva->chunkfrags);
    heapva->chunkfrags = CreateVA_AChunkFrag(0);
    assert(NULL == heapva->thechunks);
    heapva->thechunks = CreateVA_AChunkMesg(0);
#endif
    
    ierr = fclose(fp); fp = NULL;
    assert(ierr == 0);

  } else {
    /* createmode */

    const size_t initial_number_of_frags =
      MAX(1,new_additional_number_of_frags);
    const size_t initial_number_of_edges =
      MAX(1,new_additional_number_of_edges);
    const size_t initial_text_length =
      MAX(1,new_additional_amount_of_text);
    
    /* Initialize the fragment and edge data structures. */
    gstate->store_version = 1;
    gstate->state_of_the_store = 1; // Just initialized.
    gstate->min_frag_iid = 0;
    gstate->max_frag_iid = 0;
    gstate->nbase_in_genome = 0;

    // consult "AS_UTL_Var.h"
    assert(NULL == heapva->frags);
    heapva->frags = CreateVA_Afragment(initial_number_of_frags);
    assert(NULL == heapva->edges);
    heapva->edges = CreateVA_Aedge(initial_number_of_edges); 
    assert(NULL == heapva->frag_annotations);
    heapva->frag_annotations = CreateVA_char(initial_text_length);
    assert(NULL == heapva->next_edge_obj);
    heapva->next_edge_obj = CreateVA_IntEdge_ID(initial_number_of_edges); 

    assert(NULL == heapva->chunkfrags);
    heapva->chunkfrags = CreateVA_AChunkFrag(0);
    assert(NULL == heapva->thechunks);
    heapva->thechunks = CreateVA_AChunkMesg(0);
  }

  assert(NULL == heapva->chunkseqs);
  heapva->chunkseqs = CreateVA_char(0);
  assert(NULL == heapva->chunkquas);
  heapva->chunkquas = CreateVA_char(0);
  assert(NULL == heapva->chunksrcs);
  heapva->chunksrcs = CreateVA_char(0);

  assert(NULL == heapva->the_ofg_messages);
  assert(NULL == heapva->the_ovl_messages);
  assert(NULL == heapva->the_ofg_source);
  assert(NULL == heapva->the_ovl_source);
  assert(NULL == heapva->the_imp_messages);
  assert(NULL == heapva->the_ium_messages);
  assert(NULL == heapva->the_imp_source);
  assert(NULL == heapva->the_ium_source);
}


void write_fgb_store
( const char * const theStorePath,
  const TStateGlobals * const gstate,
  const THeapGlobals  * const heapva
  )
{
  int ierr=0;

#ifdef DONT_USE_STORE
  return;
#endif

  assert( NULL != theStorePath );
  assert( '\0' != theStorePath[0] );

  fprintf(stderr,"* Overwriting the store %s\n", theStorePath);
  //fprintf(stderr,"* theStorePath = %p\n", theStorePath);
  //fprintf(stderr,"* theStorePath[0] = %d\n", (int) (theStorePath[0]));
  //fprintf(stderr,"* Writing the store %s\n", theStorePath);

  // sprintf(buffer,"%s/fgb.db", theStorePath);

  /* Make a structure suitable for regression testing. */
  { 
    size_t nitems = 0;
    int ierr = fprintf(stderr,"theStorePath=<%s>\n", theStorePath);
    FILE * fp = fopen(theStorePath,"w");
    assert(fp != NULL);

    nitems = fwrite(gstate,sizeof(TStateGlobals),1,fp);
    assert(nitems == 1);
    
    CopyToFileVA_Afragment((heapva->frags),fp);
    CopyToFileVA_Aedge((heapva->edges),fp);
    CopyToFileVA_char((heapva->frag_annotations),fp);
    CopyToFileVA_IntEdge_ID((heapva->next_edge_obj),fp);
    
    CopyToFileVA_AChunkFrag(heapva->chunkfrags,fp);
    CopyToFileVA_AChunkMesg(heapva->thechunks,fp);
    //CopyToFileVA_char(heapva->chunkseqs,fp);
    //CopyToFileVA_char(heapva->chunkquas,fp);
    //CopyToFileVA_char(heapva->chunksrcs,fp);
  
    ierr = fflush(fp);
    assert(ierr == 0);
    ierr = fclose(fp); fp = NULL;
    assert(ierr == 0);
  }

#if 1
  { 
    // Now check the store for correctness!
    FILE * fp = fopen(theStorePath,"r");
    size_t nitems = 0;
    TStateGlobals ctmp = {0};
    assert(fp != NULL);
    fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS\n");
    nitems = fread(&ctmp,sizeof(TStateGlobals),1,fp);
    assert(nitems == 1);
    assert(ctmp.min_frag_iid == gstate->min_frag_iid);
    assert(ctmp.max_frag_iid == gstate->max_frag_iid);
    assert(ctmp.nbase_in_genome == gstate->nbase_in_genome);
    
    fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS frags\n");
    CheckFile_VA((heapva->frags),fp);
    fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS edges\n");
    CheckFile_VA((heapva->edges),fp);
    fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS annos\n");
    CheckFile_VA((heapva->frag_annotations),fp);
    fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS next_edge_obj\n");
    CheckFile_VA((heapva->next_edge_obj),fp);

    fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS chunkfrags\n");
    CheckFile_VA(heapva->chunkfrags,fp);
    fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS thechunks\n");
    CheckFile_VA(heapva->thechunks,fp);

    ierr = fclose(fp); fp = NULL;
    assert(ierr == 0);
    fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS done\n");
  }
#endif
  fprintf(stderr,"wrote the store\n");
}

void check_fgb_store
( const char * const theStorePath,
  const TStateGlobals * const gstate,
  const THeapGlobals  * const heapva
)
{
  TStateGlobals gtmp = {0};
  int ierr=0;
  size_t nitems = 0;
  FILE *fp;

#ifdef DONT_USE_STORE
  return;
#endif

  fp = fopen(theStorePath,"r");
  assert(fp != NULL);

  fprintf(stderr,"Now check the store for correctness!\n");

  fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS\n");
  nitems = fread(&gtmp,sizeof(TStateGlobals),1,fp);
  assert(nitems == 1);
  assert(gtmp.min_frag_iid == gstate->min_frag_iid);
  assert(gtmp.max_frag_iid == gstate->max_frag_iid);
  assert(gtmp.nbase_in_genome == gstate->nbase_in_genome);
  
  fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS frags\n");
  CheckFile_VA((heapva->frags),fp);
  fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS edges\n");
  CheckFile_VA((heapva->edges),fp);
  fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS annos\n");
  CheckFile_VA((heapva->frag_annotations),fp);
  fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS next_edge_obj\n");
  CheckFile_VA((heapva->next_edge_obj),fp);
  
  fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS chunkfrags\n");
  CheckFile_VA((heapva->chunkfrags),fp);
  fprintf(stderr,"CHECK CHECK-POINTING FOR CORRECTNESS thechunks\n");
  CheckFile_VA((heapva->thechunks),fp);

  ierr = fclose(fp); fp = NULL;
  assert(ierr == 0);

  fprintf(stderr,"checked the store\n");
  return;
}

void clear_fgb_store
(
  TStateGlobals * gstate,
  THeapGlobals  * heapva,
  int clear_the_frags,
  int clear_the_edges,
  int clear_the_frag_annotations,
  int clear_the_next_edge)
{

  fprintf(stderr,"* Clearing the FGB store\n");
  fprintf(stderr,"** clear_the_frags=%d\n", clear_the_frags);
  fprintf(stderr,"** clear_the_edges=%d\n", clear_the_edges);
  fprintf(stderr,"** clear_the_frag_annotations=%d\n", clear_the_frag_annotations);
  fprintf(stderr,"** clear_the_next_edge=%d\n", clear_the_next_edge);
  
  assert(gstate != NULL);
  assert(heapva != NULL);
  
  if(clear_the_frags) {
    clear_the_edges |= TRUE;
    clear_the_frag_annotations |= TRUE;

    /* Initialize the fragment and edge data structures. */
    gstate->store_version = 1;
    gstate->state_of_the_store = 1; // Just initialized.
    
    gstate->min_frag_iid = 0;
    gstate->max_frag_iid = 0;
    gstate->nbase_in_genome = 0;
    ResetToRange_VA(heapva->frags,0);
  }

  if(clear_the_frag_annotations) {
    ResetToRange_VA(heapva->frag_annotations,0);
  }
    
  if(clear_the_edges) {
    clear_the_next_edge |= TRUE;
#if 0
    clear_non_champion_edges |= TRUE;
#endif
    
    // Clear_VA(globals->edges); // Frees the memory
    ResetToRange_VA(heapva->edges,0);
    
    {
      // Need to re-set adjacency data for Reaper:
      const IntFragment_ID nfrag = GetNumFragments(heapva->frags);
      //const IntEdge_ID nedge = GetNumEdges(heapva->edges);
      IntFragment_ID iv;
      for(iv=0;iv<nfrag;iv++) {
        set_lab_fragment(heapva->frags,iv,AS_CGB_UNLABELED_FRAG);
        set_seglen_vertex(heapva->frags,iv,0,0);
        set_seglen_vertex(heapva->frags,iv,1,0);
      }
    }
  }

#if 0
  if(clear_non_champion_edges) {
    clear_the_next_edge |= TRUE;
    
    // Clear_VA(globals->edges); // Frees the memory
    // ResetToRange_VA(heapva->edges,0);
    
    {
      // Need to re-set adjacency data for Reaper:
      // const IntFragment_ID nfrag = GetNumFragments(heapva->frags);
      const IntEdge_ID nedge = GetNumEdges(heapva->edges);
      //IntFragment_ID iv;
      IntEdge_ID ie;
      for(ie=0;ie<nege;ie++) {

        set_lab_fragment(heapva->frags,iv,AS_CGB_UNLABELED_FRAG);
        set_seglen_vertex(heapva->frags,iv,0,0);
        set_seglen_vertex(heapva->frags,iv,1,0);
      }
    }
  }
#endif
  
  if(clear_the_next_edge) {

    // Clear_VA(heapva->next_edge_obj); // Frees the memory
    ResetToRange_VA(heapva->next_edge_obj,0);

    {
      // Need to re-set adjacency data for Reaper:
      const IntFragment_ID nfrag = GetNumFragments(heapva->frags);
      //const IntEdge_ID nedge = GetNumEdges(heapva->edges);
      IntFragment_ID iv;
      for(iv=0;iv<nfrag;iv++) {
        set_seglen_dvt_vertex(heapva->frags,iv,FALSE,0);
        set_seglen_dvt_vertex(heapva->frags,iv,TRUE,0);
        set_seglen_frc_vertex(heapva->frags,iv,FALSE,0);
        set_seglen_frc_vertex(heapva->frags,iv,TRUE,0);

        set_segstart_vertex(heapva->frags,iv,FALSE,0);
        set_segstart_vertex(heapva->frags,iv,TRUE,0);
        set_segend_vertex(heapva->frags,iv,FALSE,0);
        set_segend_vertex(heapva->frags,iv,TRUE,0);

        set_raw_dvt_count_vertex(heapva->frags,iv,FALSE,0);
        set_raw_dvt_count_vertex(heapva->frags,iv,TRUE,0);
        set_raw_frc_count_fragment(heapva->frags,iv,0);
        set_raw_toc_count_fragment(heapva->frags,iv,0);
      }
    }
  }

  assert(heapva->edges != NULL);
  assert(heapva->frags != NULL);
  assert(heapva->frag_annotations != NULL);
  assert(heapva->next_edge_obj != NULL);
  
  fprintf(stderr,"* Partially cleared the FGB store\n");
  fprintf(stderr,"** cleared the frags=%d\n", clear_the_frags);
  fprintf(stderr,"** cleared the edges=%d\n", clear_the_edges);
  fprintf(stderr,"** cleared the frag annotations=%d\n", clear_the_frag_annotations);
  fprintf(stderr,"** cleared the next edge=%d\n", clear_the_next_edge);

}
