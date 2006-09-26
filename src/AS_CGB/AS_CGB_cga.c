
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
static char CM_ID[] 
= "$Id: AS_CGB_cga.c,v 1.7 2006-09-26 22:21:13 brianwalenz Exp $";
/*********************************************************************
 *
 * Module: AS_CGB_cga.c
 * 
 * Description: A chunk graph analyzer. This functional unit computes
 * graph statistics, and writes the chunk graph in the term
 * representation which is suitable for "DaVinci".
 *
 * Assumptions: 
 * 
 * Dependencies: Fragment overlap store and the chunk graph builder
 * output.
 *
 * Author: Clark Mobarry
 ********************************************************************/

/*************************************************************************/
/* System include files */

/*************************************************************************/
/* Local include files */
#include "AS_CGB_all.h"
#include "AS_CGB_myhisto.h"

/*************************************************************************/
/* File Scope Global Variables */
static int TIMINGS = TRUE;

/*************************************************************************/
/* Conditional compilation */

#define DEBUGGING
#undef DEBUGGING

#define DEBUG_VISUAL
#undef DEBUG_VISUAL

#ifdef DEBUGGING
#define ORDERING
#define DEBUG01
#define DEBUG8
#define DEBUG78
#undef DEBUG14
#undef DEBUG15
#define DEBUG16
#define DEBUG17
#define DEBUG77
#undef DEBUG_SHORTCUT
#endif // DEBUGGING

/*************************************************************************/
/* Global Defines */

#define GRAPHID(ifrag) ifrag

#define NLETTERS 26

typedef struct {
  BPTYPE atip,btip;    /* Coordinates from 0..(length of genome). */
  int abforward;
  // BPTYPE gmin,gmax;    /* Coordinates from 0..(length of genome). */
  char essential_type; /* 'u'="true unique",'r'="repeat", '@'="unknown" */
  char contained_type; /* 'u'="true unique",'r'="repeat", '@'="unknown" */
} ChunkAnnotation;


typedef enum {
  CGA_SINGLETON_CONTAINED,
  // multiply contained and orphaned contained fragments that form
  // singleton unitigs.
  CGA_SINGLETON_NONCONTAINED,
  // The rest of the singleton unitigs.
  CGA_NONSINGLETON_SPANNED,
  // The non-singleton unitigs that are spanned by one fragment.
  CGA_NONSINGLETON_NONSPANNED
} ChunkLabel;

#define MAX_NUM_CHUNK_LABELS 4

char * ChunkLabelDesc[MAX_NUM_CHUNK_LABELS] = {
  "singleton contained ",  
  // Contained fragments that form singleton unitigs.
  "singleton non-contained",
  // The rest of the singleton unitigs.
  "non-singleton spanned",
  // The non-singleton unitigs that are spanned by one fragment.
  "non-singleton non-spanned"
  // The rest of the non-singleton unitigs.
};



/*************************************************************************/
/*************************************************************************/

//#define ABS(a) (((a)>0)?(a):(-(a)))

typedef struct {
  int32 key;
  int32 index;
} Arec;

static int compare1(const Arec * const a, const Arec * const b) 
{
  /* This comparison function is to be used with ANSI qsort() for
     sorting the edges to form contiguous segments. */
  int icom;
  icom = (a->key - b->key);
  return icom ;
}

static void plot_edges_from_vertex
(
 FILE *                   fp, 
 const Tfragment * const  frags,
 const Tedge * const      edges,
 const IntFragment_ID * const fragment_visited,
 const Tlab               fragment_types_to_plot,
 const Tnes               edge_types_to_plot,
 const IntFragment_ID     ifrag, 
 const IntFragment_ID     iend, /* The vertex is a fragment and a suffix flag. */
 const int shortcut)
 
{	
  IntEdge_ID ir;
  const IntEdge_ID ir0  = get_segstart_vertex(frags,ifrag,iend);
  const int       nnode = get_seglen_vertex(frags,ifrag,iend);
  
  for(ir=ir0; (ir<ir0+nnode); ir++) 
    /* Begin plotting the edges from the vertex. */ {
    IntFragment_ID icount=0;
    IntFragment_ID iavx = get_avx_edge(edges,ir);
    IntFragment_ID ibvx = get_bvx_edge(edges,ir);
    int iasx = get_asx_edge(edges,ir);
    int ibsx = get_bsx_edge(edges,ir);
    Tlab inesvb = get_lab_fragment(frags,ibvx);
    Tnes inese = get_nes_edge(edges,ir);

#ifdef DEBUG01
    if((ifrag != iavx) || (iend != iasx)) {
      int ibhg = get_bhg_edge(edges,ir);
      int iahg = get_ahg_edge(edges,ir);
      printf(F_IID "," F_IID ", " F_IID " : " F_IID ",%d,%d, " F_IID ",%d,%d : %d\n",
	     ifrag,iend,ir,iavx,iasx,iahg,ibvx,ibsx,ibhg,inese);
    }
#endif
    assert(ifrag == iavx);
    assert(iend == iasx);

    if( shortcut &&
	((get_lab_fragment(frags,iavx) == AS_CGB_INTERCHUNK_FRAG) ||
	 (get_lab_fragment(frags,iavx) == AS_CGB_HANGING_CHUNK_FRAG) ) &&
	(inese == AS_CGB_INTRACHUNK_EDGE) )
      /* Begin shortcut to the end of a chunk. */ {
      IntFragment_ID iv0;
      int        is0;
      IntEdge_ID ie1,ie0;
      int        ne0;
      FILE *     fout;
      fout = stdout;
#ifdef DEBUG_SHORTCUT
      { 
	const IntFragment_ID iafr = get_iid_fragment(frags,iavx);
	const IntFragment_ID ibfr = get_iid_fragment(frags,ibvx);
	fprintf(stderr,"Start shortcut "
		"iafr=" F_IID ",iavx=" F_IID ",ibvx=" F_IID ",inese=%d=\n",
		iafr,iavx,ibvx,inese);
      }
#endif
      for(iv0=ibvx,is0= !ibsx;
	    (get_lab_fragment(frags,iv0) == AS_CGB_INTRACHUNK_FRAG);)
	/* We process the INTRACHUNK fragments here. */ {
	
	/* Search all vertices adjacent from "iv0" */
	ie0 = get_segstart_vertex(frags,iv0,is0);
	ne0 = get_seglen_vertex(frags,iv0,is0);
#ifdef DEBUG_SHORTCUT
	fprintf(fout,"shortcut: fragment " F_IID " " F_IID " %d \n",iv0,ie0,ne0);
#endif
	for(ie1=ie0; ie1 < ie0+ne0; ie1++) {
	  const IntFragment_ID iv1 = get_bvx_edge(edges,ie1);
	  const int        is1 = get_bsx_edge(edges,ie1);
	  const Tnes       ines1 = get_nes_edge(edges,ie1);
	  if((ines1 == AS_CGB_INTRACHUNK_EDGE)
	     )
	    { /* iv1 is unexplored */
#ifdef DEBUG_SHORTCUT
	      fprintf(fout,"shortcut: iv1=" F_IID ",ines1=%d\n",iv1,ines1);
#endif
	      ibvx = iv1; 
	      ibsx = is1;
	      iv0 = iv1; 
	      is0 = ! is1; /* search from the other side of the fragment. */
	      icount ++;
	      /* break; */
	    }
	}
      }
      inesvb = get_lab_fragment(frags,ibvx);
#ifdef DEBUG_SHORTCUT
      { 
	const IntFragment_ID iafr = get_iid_fragment(frags,iavx);
	const IntFragment_ID ibfr = get_iid_fragment(frags,ibvx);
	fprintf(stderr,"End shortcut "
		"iafr=" F_IID ",iavx=" F_IID ",iasx=%d, ibfr=" F_IID ",ibvx=" F_IID ",ibsx=%d,"
		"inesvb=%d\n",
		iafr,iavx,iasx,get_iid_fragment(frags,ibvx),ibvx,ibsx,
		inesvb);
      }
#endif
      assert(inesvb==AS_CGB_INTERCHUNK_FRAG ||
	     inesvb==AS_CGB_HANGING_CHUNK_FRAG );
    } /* End shortcut to the end of a chunk. */
	
    /* LOGICAL START OF THE EDGE PLOTTING. */

    if( ((inese <= edge_types_to_plot) ||
         (inese == AS_CGB_TOUCHES_CONTAINED_EDGE)) &&
        ((fragment_visited[ibvx] != FRAGMENT_NOT_VISITED) ||
         ((inese == AS_CGB_CONTAINED_EDGE)&&
          (inesvb <= fragment_types_to_plot))) )
      /* Begin do we want to plot this edge? */ {
      char color_string[80];
      char essential_string[80];
      char vs3[80];
	  
#ifdef ORDERING
      // printf(" icount=" F_IID " iavx=" F_IID ",ibvx=" F_IID "\n",icount,iavx,ibvx);
#endif	    
#ifdef DEBUG_SHORTCUT
      fprintf(stderr,"plot edge iavx=" F_IID ",ibvx=" F_IID "\n",iavx,ibvx);
#endif
      /* Format the bvx index to a string. */
      sprintf(vs3,F_IID,GRAPHID(ibvx));
	      
      sprintf(essential_string," ");
      if( inese == AS_CGB_INTRACHUNK_EDGE ) 
	{
	  sprintf(essential_string,"a(\"EDGEPATTERN\",\"thick\"),");
	  sprintf(essential_string,"a(\"EDGEPATTERN\",\"double\"),");
	}
	  

      if(
         inese == AS_CGB_CONTAINED_EDGE
         ) /* Contained overlaps only */
	{
	  sprintf(essential_string,"a(\"EDGEPATTERN\",\"dashed\"),");
	}
	  
      if(iasx && !ibsx ) /* Normal */
	{
	  sprintf(color_string,"a(\"EDGECOLOR\",\"blue\"),");
	}
      if( !iasx && ibsx ) /* Anti-normal */
	{
	  sprintf(color_string,
		  "a(\"EDGECOLOR\",\"cyan\"),a(\"_DIR\",\"inverse\"),");
	}
      if( iasx && ibsx ) /* AS_INNIE */
	{
	  sprintf(color_string,
		  "a(\"EDGECOLOR\",\"green\"),a(\"_DIR\",\"none\"),");
	}
      if( !iasx && !ibsx ) /* AS_OUTTIE */
	{
	  sprintf(color_string,
		  "a(\"EDGECOLOR\",\"red\"),a(\"_DIR\",\"both\"),");
	}
      assert(NULL != fp);
      fprintf(fp,
	      "       e(\"\",[%s %s],r(\"%s\")),\n", 
	      color_string, essential_string,vs3);
    } /* End of do we plot this particular edge ? */
  } /* End plotting the edges from the fragment. */
}

static void exhale_term_rep
(FILE *fp, 
 const Tfragment * const frags,
 const Tedge * const edges,
 const Tlab fragment_types_to_plot,
 const Tnes edge_types_to_plot,
 const int shortcut,
 const Tfraginfo * const fraginfo)
{
  const IntFragment_ID nfrag = GetNumFragments(frags);

  Arec * fragment_mapping = NULL;
  IntFragment_ID * fragment_ranking = NULL;
  IntFragment_ID * fragment_visited = NULL;

  fragment_mapping = safe_malloc(sizeof(Arec) * nfrag);
  fragment_ranking = safe_malloc(sizeof(IntFragment_ID) * nfrag);
  fragment_visited = safe_malloc(sizeof(IntFragment_ID) * nfrag);

  {
    time_t tp1,tp2;
    IntFragment_ID ifrag;
    
#ifdef DEBUGGING
    fprintf(stderr,"Sort the fragments by genome position.\n");
#endif
    if(TIMINGS) {
      time(&tp1); fprintf(stderr,"Begin sort\n");
    }
    for(ifrag=0; ifrag<nfrag; ifrag++) {
      if(fraginfo != NULL) {
	fragment_mapping[ifrag].key = (int32) 
	  get_genbgn_fraginfo(fraginfo,ifrag);
      } else {
	fragment_mapping[ifrag].key = ifrag;
      }
      fragment_mapping[ifrag].index = ifrag;
    }
    qsort(fragment_mapping,nfrag,sizeof(Arec),
	  (int (*)(const void *,const void *))compare1); 
    for(ifrag=0; ifrag<nfrag; ifrag++) {
      fragment_ranking[fragment_mapping[ifrag].index] = ifrag;
    }
    
    /* Check for a valid permutation. */
    for(ifrag=0;ifrag<nfrag;ifrag++) 
      { fragment_visited[ifrag] = 0;}
    for(ifrag=0;ifrag<nfrag;ifrag++) 
      { fragment_visited[fragment_ranking[ifrag]] ++;}
    for(ifrag=0;ifrag<nfrag;ifrag++) 
      { assert(fragment_visited[ifrag] == 1);}
      
    if(TIMINGS) {
      time(&tp2); fprintf(stderr,"%10" F_TIME_TP " sec: Finished sort\n",(tp2-tp1));
    }
  }

  /* Initialize a flag for chunk following. */
  { 
    IntFragment_ID ifrag;
    for(ifrag=0;ifrag<nfrag;ifrag++) 
      { fragment_visited[ifrag] = FRAGMENT_NOT_VISITED;}
  }
  
  assert(NULL != fp);
  { 
    IntFragment_ID ifrag_raw;
    int pass;
    fprintf(fp,"[\n"); /* beginning the graph */

    for(pass=0; pass<2; pass++) {
    for(ifrag_raw=0;ifrag_raw<nfrag;ifrag_raw++)
      /* beginning loop through fragments */ {
      /* fragment strings */
      char vs0[200], vs1[200], color_string[80];
      const IntFragment_ID ifrag = fragment_mapping[ifrag_raw].index;
      const IntFragment_ID iafr = get_iid_fragment(frags,ifrag);
      const Tlab       ilab = get_lab_fragment(frags,ifrag);
      
    if(
      (fragment_visited[ifrag] == FRAGMENT_NOT_VISITED) &&
      ((ilab <= fragment_types_to_plot) ||
       (ilab == AS_CGB_BRANCHMULTICONT_FRAG)) &&
      (
        ((pass==0)&&(
          (!shortcut || ilab != AS_CGB_INTRACHUNK_FRAG) &&
          (FALSE == get_con_fragment(frags,ifrag))
          ))
        ||
        ((pass==1)&&(
          (!shortcut || ilab != AS_CGB_SINGLECONT_FRAG) &&
           (TRUE == get_con_fragment(frags,ifrag))
          ))
        )

	)
      /* begin do we plot this fragment ? */ {

      fragment_visited[ifrag] = ifrag; 
      /* mark this fragment unique to this pass */

      /* Format the avx index into a string. */
      sprintf(vs0,"\"" F_IID "\"",GRAPHID(ifrag));
      /* Format the afr index into a string. */
      if(fraginfo == NULL) {
	sprintf(vs1,"a(\"OBJECT\",\"(" F_IID ":" F_IID ")\"),",
		iafr,ifrag);
      } else {
#ifdef NEVER
	sprintf(vs1,"a(\"OBJECT\",\"[%d,%d] %c %c (" F_IID ":" F_IID ")\"),",
		(int32)get_genbgn_fraginfo(fraginfo,ifrag),
		(int32)get_genend_fraginfo(fraginfo,ifrag),
		get_pre_let_fraginfo(fraginfo,ifrag),
		get_suf_let_fraginfo(fraginfo,ifrag),
		iafr,ifrag);
#else
#ifdef GENINFO
	sprintf(vs1,"a(\"OBJECT\",\"[%d,%d] (" F_IID ":" F_IID ")\"),",
		(int32)get_genbgn_fraginfo(fraginfo,ifrag),
		(int32)get_genend_fraginfo(fraginfo,ifrag),
		iafr,ifrag);
#endif // GENINFO
#endif
      }
#ifdef NEVER
      fprintf(stdout,"[%d,%d] (" F_IID ":" F_IID ") %c,%d,%d,%d %c,%d,%d,%d\n",
	      (int32)get_genbgn_fraginfo(fraginfo,ifrag),
	      (int32)get_genend_fraginfo(fraginfo,ifrag),
	      iafr,ifrag,
	      get_pre_let_fraginfo(fraginfo,ifrag),
	      get_pre_ins_fraginfo(fraginfo,ifrag),
	      get_pre_brp_fraginfo(fraginfo,ifrag),
	      get_pre_end_fraginfo(fraginfo,ifrag),
	      get_suf_let_fraginfo(fraginfo,ifrag),
	      get_suf_ins_fraginfo(fraginfo,ifrag),
	      get_suf_brp_fraginfo(fraginfo,ifrag),
	      get_suf_end_fraginfo(fraginfo,ifrag));
#endif
      
      switch(ilab) {
      case AS_CGB_SOLO_FRAG:
	sprintf(color_string,"a(\"COLOR\",\"magenta\"),"); break;
      case AS_CGB_HANGING_FRAG:
	sprintf(color_string,"a(\"COLOR\",\"brown\"),"); break;
      case AS_CGB_THRU_FRAG:
	sprintf(color_string,"a(\"COLOR\",\"lightblue\"),"); break;

      case AS_CGB_ORPHANEDCONT_FRAG:
	sprintf(color_string, 
		"a(\"COLOR\",\"brown\"),a(\"_GO\",\"ellipse\"),"); break;
      case AS_CGB_BRANCHMULTICONT_FRAG:
	sprintf(color_string,
                "a(\"COLOR\",\"brown\"),a(\"_GO\",\"ellipse\"),"); break;
      case AS_CGB_MULTICONT_FRAG:
	sprintf(color_string, 
		"a(\"COLOR\",\"lightblue\"),a(\"_GO\",\"ellipse\"),"); break;
      case AS_CGB_SINGLECONT_FRAG:
	sprintf(color_string, 
		"a(\"COLOR\",\"lightgreen\"),a(\"_GO\",\"ellipse\"),"); break;

      case AS_CGB_INTERCHUNK_FRAG:
	sprintf(color_string,"a(\"COLOR\",\"green\"),"); break;
      case AS_CGB_INTRACHUNK_FRAG:
        sprintf(color_string,"a(\"COLOR\",\"lightgreen\"),"); break;
      case AS_CGB_HANGING_CHUNK_FRAG:
	sprintf(color_string,"a(\"COLOR\",\"yellow\"),"); break;

      case AS_CGB_HANGING_CRAPPY_FRAG:
	sprintf(color_string,"a(\"COLOR\",\"lightred\"),"); break;
      case AS_CGB_DELETED_FRAG:
	sprintf(color_string,"a(\"COLOR\",\"red\"),"); break;
      default:
	assert(FALSE); break;
      }

      fprintf(fp,
	      "l(%s,n(\"\",[%s %s],\n",
	      vs0,vs1,color_string);
      fprintf(fp,"       [\n");
      
      plot_edges_from_vertex(
			       fp, 
			       frags,
			       edges,
			       fragment_visited,
                               fragment_types_to_plot,
			       edge_types_to_plot,
			       ifrag, FALSE,
			       shortcut);
      
      plot_edges_from_vertex(
			       fp, 
			       frags,
			       edges,
			       fragment_visited,
                               fragment_types_to_plot,
			       edge_types_to_plot,
			       ifrag, TRUE,
			       shortcut);
    
    fprintf(fp,"])),\n");  /* ending old fragment in term representation */
    } /* end do we plot this fragment? */
  } /* ending loop through fragments */
    } /* ending loop through passes */
  fprintf(fp,"]\n");    /* ending the graph */
  }

  safe_free(fragment_ranking);
  safe_free(fragment_mapping);
  safe_free(fragment_visited);
}


static void graph_diagnostics
(
 const char root[],
 const Tfragment * const frags,
 const Tedge * const edges,
 const Tlab fragment_types_to_plot,
 const Tnes edge_types_to_plot,
 const int shortcut,
 const Tfraginfo * const fraginfo
 ) 
{
  char fname[80]={0};
  { 
    FILE  *fpg2 = NULL;
    assert(NULL != root);
    assert('\0' != root[0]);
    sprintf(fname,"%s.daVinci",root);
    fpg2 = fopen(fname,"w");
    assert(fpg2 != NULL);
    exhale_term_rep(fpg2, frags, edges, 
		    fragment_types_to_plot, edge_types_to_plot,shortcut,
		    fraginfo);
    if(fclose(fpg2) == EOF)
      assert(0);
  }
}

static void annotate_the_chunks_with_coordinate_info
(/*Input Only*/
 const IntFragment_ID max_frag_iid,
 const Tfragment frags[],
 const Tedge edges[],
 const Tfraginfo fraginfo[],
 const TChunkFrag chunkfrags[],
 const TChunkMesg thechunks[],
 /*Output Only*/
 ChunkAnnotation chunkinfo[])
{ 
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);


  IntChunk_ID ichunk;

  assert((!0) == 1); /* Needed for the following bitwise XOR operator. */
  for(ichunk=0;ichunk<nchunks;ichunk++)
    /* beginning loop through fragment */ {
    // Process the chunk-end fragments first.
    const IntFragment_ID chunk_avx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_avx;
    const IntFragment_ID chunk_bvx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_bvx;
    const int chunk_asx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_asx;
    const int chunk_bsx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_bsx;
    const IntFragment_ID nfrag_in_chunk
      = GetVA_AChunkMesg(thechunks,ichunk)->num_frags;
    // const BPTYPE nbase_essential_in_chunk 
    // = GetVA_AChunkMesg(thechunks,ichunk)->bp_length;
    const IntFragment_ID irec_start_of_chunk
      = GetVA_AChunkMesg(thechunks,ichunk)->f_list;

#ifdef GENINFO
    BPTYPE atip=0,btip=0;
    BPTYPE lowest_genome_coordinate=0,highest_genome_coordinate=0;
    int ivote_repeat_essential=0;
    int ivote_repeat_contained=0;
    int abforward=-1, direction_error = 0;
#endif
    int iavx=-1,ibvx=-1; 
    int iasx=-1,ibsx=-1;

#ifdef DEBUG07
    fprintf(stderr,
	    "Process ichunk,nchunks,nfrag_in_chunk,nbase_essential_in_chunk=\n"
	   F_IID "," F_IID "," F_IID "," BPFORMAT "\n",
	   ichunk,nchunks,nfrag_in_chunk,nbase_essential_in_chunk);
#endif

    { // Begin processing the A fragment of the chunk.
      const IntFragment_ID vid = chunk_avx;
      // Synthesize the original INTRACHUNK edge.
      iavx = chunk_avx;
      iasx = get_forward_fragment(frags,chunk_avx);
#ifdef GENINFO
      { /* shotgun simulator info */
	const BPTYPE ibgn = get_genbgn_fraginfo(fraginfo,vid);
	const BPTYPE iend = get_genend_fraginfo(fraginfo,vid);
	const BPTYPE ilow = (ibgn < iend ? ibgn : iend);
	const BPTYPE ihgh = (ibgn > iend ? ibgn : iend);
	atip = (iasx ? ibgn : iend);

	abforward = ((iend<ibgn) ^ iasx);
	lowest_genome_coordinate  = ilow;
	highest_genome_coordinate = ihgh;
      }
#endif // GENINFO
    } // Finished processing the A fragment of the chunk.

    { // Begin processing the B fragment of the chunk.
      const IntFragment_ID vid = chunk_bvx;
      ibvx = chunk_bvx;
      ibsx = ! get_forward_fragment(frags,chunk_bvx);
#ifdef GENINFO
      { /* shotgun simulator info */
	const BPTYPE ibgn = get_genbgn_fraginfo(fraginfo,vid);
	const BPTYPE iend = get_genend_fraginfo(fraginfo,vid);
	const BPTYPE ilow = (ibgn < iend ? ibgn : iend);
	const BPTYPE ihgh = (ibgn > iend ? ibgn : iend);
	btip = (ibsx ? ibgn : iend);
	lowest_genome_coordinate = (lowest_genome_coordinate < ilow ?
				    lowest_genome_coordinate : ilow);
	highest_genome_coordinate = (highest_genome_coordinate > ihgh ?
				     highest_genome_coordinate : ihgh);
      }
#endif
    } // Finished processing the B fragment of the chunk.

    if((chunk_avx == chunk_bvx)&&(chunk_asx != chunk_bsx) ) {
      // A spanned chunk.
    } else if ((chunk_avx == chunk_bvx)&&(chunk_asx == chunk_bsx) ) {
      // A circular chunk?
      fprintf(stderr,"A circular chunk???\n");
      assert(FALSE);
    }
    
    // Check if all of the essential edges in the chunk are valid
    // overlaps.  Also check is the contained edges are valid
    // overlaps.
    { IntFragment_ID ifrag;
    for(ifrag=0;ifrag<nfrag_in_chunk;ifrag++){
      const IntFragment_ID ivc = irec_start_of_chunk + ifrag;
      const IntFragment_ID vid = GetVA_AChunkFrag(chunkfrags,ivc)->vid;
      const Tlab ilabel = get_lab_fragment(frags,vid);
#ifdef GENINFO
      /* shotgun simulator info */
      const BPTYPE ibgn = get_genbgn_fraginfo(fraginfo,vid);
      const BPTYPE iend = get_genend_fraginfo(fraginfo,vid);
#endif
      
      // Synthesize the original INTRACHUNK edge.
      ibvx = vid;
      ibsx = ! get_forward_fragment(frags,vid);
#ifdef GENINFO2
      // Find the genome coordinate span of the essential fragments of
      // the chunk.
      if( ilabel != AS_CGB_SINGLECONT_FRAG ) // Work-around for H.Flu.
      {
	BPTYPE ilow,ihgh;
	ilow = (ibgn < iend ? ibgn : iend);
	ihgh = (ibgn < iend ? iend : ibgn);
	if((ibgn == 0) && (iend == 0)) {
	  lowest_genome_coordinate = ilow;
	  highest_genome_coordinate = ihgh;
	} else {
	  lowest_genome_coordinate = (lowest_genome_coordinate < ilow ?
				      lowest_genome_coordinate : ilow);
	  highest_genome_coordinate = (highest_genome_coordinate > ihgh ?
				       highest_genome_coordinate : ihgh);
	}
      }
#endif
      
      switch(ilabel){
      case AS_CGB_SOLO_FRAG:
      case AS_CGB_HANGING_FRAG:
      case AS_CGB_THRU_FRAG:
      case AS_CGB_MULTICONT_FRAG:
      case AS_CGB_BRANCHMULTICONT_FRAG:
      case AS_CGB_HANGING_CHUNK_FRAG:
      case AS_CGB_ORPHANEDCONT_FRAG:
      case AS_CGB_INTERCHUNK_FRAG:
      case AS_CGB_INTRACHUNK_FRAG:
      case AS_CGB_HANGING_CRAPPY_FRAG:
	{
	  int invalid;
	  assert(iavx != -1);
	  assert(iasx != -1);
#ifdef GENINFO
	  assert(abforward != -1);
	  direction_error += (abforward ^ ((iend>ibgn) ^ ibsx));
#if 0
          {
            const IntFragment_ID iid = get_iid_fragment(frags,vid);

            fprintf(stderr,"iid,ifrag,abforward,ibgn,iend,ibsx=" F_IID "," F_IID ",%d," BPFORMAT "," BPFORMAT ",%d\n",
                    iid,ifrag,abforward,ibgn,iend,ibsx);
          }
#endif
	  invalid = check_overlap_with_simulator(nfrag,fraginfo,iavx,ibvx);
	  if(invalid) {
            const int acon = get_con_fragment(frags,iavx);
            const int bcon = get_con_fragment(frags,ibvx);
            const int alab = get_lab_fragment(frags,iavx);
            const int blab = get_lab_fragment(frags,ibvx);
            const IntFragment_ID aiid = get_iid_fragment(frags,iavx);
            const IntFragment_ID biid = get_iid_fragment(frags,ibvx);
            const IntChunk_ID acid = get_cid_fragment(frags,iavx);
            const IntChunk_ID bcid = get_cid_fragment(frags,ibvx);
            
	    fprintf(stdout,"INVALID OVL: "
                    "afr=" F_IID " asx=%d acon=%d alab=%d acid=" F_IID " bfr=" F_IID " bsx=%d bcon=%d blab=%d bcid=" F_IID "\n",
		    aiid,iasx, acon, alab, acid,
		    biid,ibsx, bcon, blab, bcid);
	  }
	  ivote_repeat_essential += invalid;
#endif // GENINFO
	  /* Store this information in the CHK message? */
	}
	break;
      case AS_CGB_SINGLECONT_FRAG:
	{
#ifdef GENINFO
	  assert(abforward != -1);
	  assert(iavx != -1);
	  assert(iasx != -1);
	  {
	    const int invalid
	      = check_overlap_with_simulator(nfrag,fraginfo,iavx,ibvx);
	    ivote_repeat_contained += invalid;

	    // Should the contained fragments contribute to a
	    // unique/repeat discrimination?
#endif // GENINFO
	  }
	}
	break;
      default:
	assert(FALSE);
      }
      if(ilabel != AS_CGB_SINGLECONT_FRAG) { 
	// Set up for the next fragment.
	assert(ibvx != -1);
	assert(ibsx != -1);
	iavx = ibvx; 
	iasx = ! ibsx;
      }
    }}

    chunkinfo[ichunk].atip    = 0;
    chunkinfo[ichunk].btip    = 0;
    chunkinfo[ichunk].abforward = abforward;
    // chunkinfo[ichunk].gmin    = 0;
    // chunkinfo[ichunk].gmax    = 0;
    chunkinfo[ichunk].essential_type = '@';
    chunkinfo[ichunk].contained_type = '@';
#ifdef GENINFO
    {

#ifdef GENINFO3
      const BPTYPE span_of_genome_coordinates
	= highest_genome_coordinate - lowest_genome_coordinate;
      const BPTYPE temp1 = ABS(span_of_genome_coordinates 
			       - GetVA_AChunkMesg(thechunks,ichunk)->bp_length);
      const int unique_in_genome 
	=  temp1 < 0.04*(GetVA_AChunkMesg(thechunks,ichunk)->bp_length);
#endif // GENINFO3

      chunkinfo[ichunk].atip = atip;
      chunkinfo[ichunk].btip = btip;

      if( (ivote_repeat_essential == 0) &&
	  ((min(atip,btip) != lowest_genome_coordinate) ||
	   (max(atip,btip) != highest_genome_coordinate))) {
	fprintf(stderr,"CGA simulator coordinate problem with ichunk=" F_IID "\n",ichunk);
	fprintf(stderr,"    atip=" BPFORMAT ",btip=" BPFORMAT "\n", atip, btip);
	fprintf(stderr,"    lowest_genome_coordinate=" BPFORMAT ",highest_genome_coordinate=" BPFORMAT "\n",
		lowest_genome_coordinate, highest_genome_coordinate);
      }


      chunkinfo[ichunk].essential_type = 
	( ivote_repeat_essential == 0 ? 'u' : 'r' );
      chunkinfo[ichunk].contained_type = 
	( ivote_repeat_contained == 0 ? 'u' : 'r' );

#ifdef DEBUG77
      if( (chunkinfo[ichunk].essential_type == 'r' ) && 
	  (GetVA_AChunkMesg(thechunks,ichunk)->coverage_stat 
	   >= cgb_unique_cutoff) ) {
	fprintf(stderr,
		"ichunk=" F_IID " "
		"is an essential invalid unitig with coverage_stat=%f\n",
		ichunk,
		GetVA_AChunkMesg(thechunks,ichunk)->coverage_stat);
      }
      if( (chunkinfo[ichunk].contained_type == 'r' ) && 
	  (GetVA_AChunkMesg(thechunks,ichunk)->coverage_stat 
	   >= cgb_unique_cutoff) ) {
	fprintf(stderr,
		"ichunk=" F_IID " "
		"is a contained invalid unitig with coverage_stat=%f\n",
		ichunk,
		GetVA_AChunkMesg(thechunks,ichunk)->coverage_stat);
      }
#endif/*DEBUG77*/
#ifdef DEBUG78
      if( ((!(ivote_repeat_essential==0))&&(unique_in_genome)) ||
	  (((ivote_repeat_essential==0))&&(!unique_in_genome)) ) {
	float diff = (float)(span_of_genome_coordinates 
			     - GetVA_AChunkMesg(thechunks,ichunk)->bp_length)
	  /(float)GetVA_AChunkMesg(thechunks,ichunk)->bp_length;
	fprintf(stderr,
		"ichunk=%8" F_IIDP ", unique_in_genome=%d, ivote_repeat_essential=%d, "
		"low=" BPFORMAT15 
		",high=" BPFORMAT15 " diff:%g\n",
		ichunk,
		unique_in_genome,ivote_repeat_essential,
		lowest_genome_coordinate,highest_genome_coordinate, diff);
      }
#endif/*DEBUG78*/
    }
#endif
  }

#ifdef GENINFO
#ifdef DEBUG08
  fprintf(stderr,"direction_error_in_uniques=%d\n",direction_error_in_uniques);
  fprintf(stderr,"direction_error_in_repeats=%d\n",direction_error_in_repeats);
#endif
#endif /*GENINFO*/

}


/*************************************************************************/

static void analyze_the_fragment_overlap_graph
(
 FILE *fout,
 const IntFragment_ID max_frag_iid,
 Tfragment frags[],
 Tedge edges[]
)
{
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);
  assert(NULL != fout);
  fprintf(fout,"FRAGMENT OVERLAP GRAPH INFORMATION\n\n");
  {
    IntFragment_ID 
      ifrag,
      n_as_cgb_solo_frag=0,
      n_as_cgb_hanging_frag=0,
      n_as_cgb_thru_frag=0,
      n_as_cgb_interchunk_frag=0,
      n_as_cgb_intrachunk_frag=0,
      n_as_cgb_unplacedcont_frag=0,
      n_as_cgb_singlecont_frag=0,
      n_as_cgb_multicont_frag=0,
      n_as_cgb_branchmulticont_frag=0,
      n_as_cgb_hanging_chunk_frag=0,
      n_as_cgb_hanging_crappy_frag=0,
      n_as_cgb_orphanedcont_frag=0,
      n_as_cgb_deleted_frag=0;
    
    fprintf(fout,"Fragment types\n");
    for(ifrag=0; ifrag<nfrag; ifrag++){
      const Tlab ilab = get_lab_fragment(frags,ifrag);
      switch(ilab) {
      case AS_CGB_SOLO_FRAG: 
	n_as_cgb_solo_frag++; break;
      case AS_CGB_HANGING_FRAG: 
	n_as_cgb_hanging_frag++; break;
      case AS_CGB_THRU_FRAG: 
	n_as_cgb_thru_frag++; break;
      case AS_CGB_INTERCHUNK_FRAG:
	n_as_cgb_interchunk_frag++; break;
      case AS_CGB_INTRACHUNK_FRAG:
	n_as_cgb_intrachunk_frag++; break;
      case AS_CGB_UNPLACEDCONT_FRAG:
	n_as_cgb_unplacedcont_frag++; break;
      case AS_CGB_SINGLECONT_FRAG:
	n_as_cgb_singlecont_frag++; break;
      case AS_CGB_MULTICONT_FRAG:
	n_as_cgb_multicont_frag++; break;
      case AS_CGB_BRANCHMULTICONT_FRAG:
	n_as_cgb_branchmulticont_frag++; break;
      case AS_CGB_HANGING_CHUNK_FRAG:
	n_as_cgb_hanging_chunk_frag++; break;
      case AS_CGB_HANGING_CRAPPY_FRAG:
	n_as_cgb_hanging_crappy_frag++; break;
      case AS_CGB_ORPHANEDCONT_FRAG:
	n_as_cgb_orphanedcont_frag++; break;
      case AS_CGB_DELETED_FRAG:
	n_as_cgb_deleted_frag++; break;
      default:
	fprintf(stderr,"ifrag=" F_IID ", Unknown fragment type %d\n",ifrag,ilab);
	assert(FALSE);
      }
    }
      fprintf(fout,
	      "%15" F_IIDP " : total number of fragment reads\n"
	      "%15" F_IIDP " :   solo\n"
	      "%15" F_IIDP " :   hanging    alone\n"
	      "%15" F_IIDP " :   hanging    chunk-end\n"
	      "%15" F_IIDP " :   hanging    spur\n"
	      "%15" F_IIDP " :   thru       alone\n"
	      "%15" F_IIDP " :   thru       chunk-end\n"
	      "%15" F_IIDP " :   thru       intrachunk\n"
	      "%15" F_IIDP " :   contained  unplaced\n"
	      "%15" F_IIDP " :   contained  singly\n"
	      "%15" F_IIDP " :   contained  multiply non-branch\n"
	      "%15" F_IIDP " :   contained  multiply branch\n"
	      "%15" F_IIDP " :   contained  orphaned\n"
	      "%15" F_IIDP " :   deleted\n",
	      nfrag,
	      n_as_cgb_solo_frag,
	      n_as_cgb_hanging_frag, 
	      n_as_cgb_hanging_chunk_frag,
	      n_as_cgb_hanging_crappy_frag,
	      n_as_cgb_thru_frag, 
	      n_as_cgb_interchunk_frag,
	      n_as_cgb_intrachunk_frag, 
	      n_as_cgb_unplacedcont_frag,
	      n_as_cgb_singlecont_frag,
	      n_as_cgb_multicont_frag,
	      n_as_cgb_branchmulticont_frag,
	      n_as_cgb_orphanedcont_frag,
	      n_as_cgb_deleted_frag
	      );
    }

    count_fragment_and_edge_labels( frags, edges, "in cga");

    {
      IntEdge_ID
	iedge,

	n_as_cgb_dovetail=0,
	n_as_cgb_thickest=0,
	/* The inter-chunk skeleton overlaps: */
	n_as_cgb_interchunk=0,
	n_as_cgb_containment_interchunk=0,
	n_as_cgb_touches_contained=0,
	n_as_cgb_between_contained=0,
	n_as_cgb_touches_crappy_dvt=0,
	n_as_cgb_between_crappy_dvt=0,
	n_as_cgb_touches_crappy_con=0,
	n_as_cgb_between_crappy_con=0,

	/* The intra-chunk skeleton overlaps: */
	n_as_cgb_intrachunk=0,
	n_as_cgb_containment_intrachunk=0,
	n_as_cgb_touches_singly_contained=0,

	/* The marked overlap family: */
	n_as_cgb_marked_by_branch_dvt=0,
        
	/* The removed overlap family: */
	n_as_cgb_removed_by_transitivity_dvt=0,
	n_as_cgb_removed_by_threshold_dvt=0,
	n_as_cgb_removed_by_transitivity_con=0,
	n_as_cgb_removed_by_threshold_con=0,

	n_as_cgb_removed_by_duplicate_dvt=0,
	n_as_cgb_removed_by_duplicate_con=0;

      fprintf(fout,"Overlap types\n");
      for(iedge=0; iedge<nedge; iedge++){
	Tnes ines = get_nes_edge(edges,iedge);
	switch(ines) {
	case AS_CGB_DOVETAIL_EDGE: 
	  n_as_cgb_dovetail++; break;
	case AS_CGB_THICKEST_EDGE: 
	  n_as_cgb_thickest++; break;
	case AS_CGB_INTERCHUNK_EDGE: 
	  n_as_cgb_interchunk++; break;
	case AS_CGB_INTRACHUNK_EDGE:
	  n_as_cgb_intrachunk++; break;
	case AS_CGB_TOUCHES_CONTAINED_EDGE:
	  {
	    const IntFragment_ID iavx = get_avx_edge(edges,iedge);
	    /* An index into the Tfragment array of the fragment
	       at the proximal vertex of the edge. */
	    const IntFragment_ID ibvx = get_bvx_edge(edges,iedge);
	    /* An index into the Tfragment array of the fragment
	       at the distal vertex of the edge. */

	    assert( get_con_fragment(frags,iavx) ||
		    get_con_fragment(frags,ibvx) );
	    assert( !( get_con_fragment(frags,iavx) &&
		       get_con_fragment(frags,ibvx) ));

	    if (
	     !((AS_CGB_SINGLECONT_FRAG == get_lab_fragment(frags,iavx)) ||
	       (AS_CGB_SINGLECONT_FRAG == get_lab_fragment(frags,ibvx)))
		) {
	      n_as_cgb_touches_contained++;
	    } else {
	      n_as_cgb_touches_singly_contained++;
	    }
	  }
	  break;
	case AS_CGB_BETWEEN_CONTAINED_EDGE:
	  n_as_cgb_between_contained++; break;
        case AS_CGB_TOUCHES_CRAPPY_DVT:
          n_as_cgb_touches_crappy_dvt++; break;
        case AS_CGB_BETWEEN_CRAPPY_DVT:
          n_as_cgb_between_crappy_dvt++; break;
        case AS_CGB_TOUCHES_CRAPPY_CON:
          n_as_cgb_touches_crappy_con++; break;
        case AS_CGB_BETWEEN_CRAPPY_CON:
          n_as_cgb_between_crappy_con++; break;
	case AS_CGB_CONTAINED_EDGE:
	  {
	    const IntFragment_ID iavx = get_avx_edge(edges,iedge);
	    /* An index into the Tfragment array of the fragment
	       at the proximal vertex of the edge. */
	    const IntFragment_ID ibvx = get_bvx_edge(edges,iedge);
	    /* An index into the Tfragment array of the fragment
	       at the distal vertex of the edge. */
	    if(
	       (AS_CGB_SINGLECONT_FRAG == get_lab_fragment(frags,iavx)) ||
	       (AS_CGB_SINGLECONT_FRAG == get_lab_fragment(frags,ibvx)) ) {
	      n_as_cgb_containment_intrachunk++;
	    } else {
	      n_as_cgb_containment_interchunk++; 
	    }
	  }
	  break;

	case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
	  n_as_cgb_removed_by_transitivity_dvt++; break;
	case AS_CGB_REMOVED_BY_THRESHOLD_DVT:
	  n_as_cgb_removed_by_threshold_dvt++; break;
	case AS_CGB_MARKED_BY_BRANCH_DVT:
	  n_as_cgb_marked_by_branch_dvt++; break;
	case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
	  n_as_cgb_removed_by_duplicate_dvt++; break;

	case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
	  n_as_cgb_removed_by_transitivity_con++; break;
	case AS_CGB_REMOVED_BY_THRESHOLD_CON:
	  n_as_cgb_removed_by_threshold_con++; break;
	case AS_CGB_REMOVED_BY_DUPLICATE_CON:
	  n_as_cgb_removed_by_duplicate_con++; break;
	default:
	  assert(FALSE);
	}
      }

      fprintf(fout,
	      "%15" F_IIDP " : total number of fragment overlaps\n"
	      "      The non-chunking overlaps:\n"
	      "%15" F_IIDP " : dovetail non-chunking\n"
	      "%15" F_IIDP " : thickest non-chunking\n"
	      "      The inter-chunk skeleton overlaps:\n"
	      "%15" F_IIDP " : dovetail inter-chunk (between non-contained)\n"
	      "%15" F_IIDP " : dovetail non-contained touches contained\n"
	      "%15" F_IIDP " : dovetail between contained\n"
	      "%15" F_IIDP " : dovetail touches spur and non-spur\n"
	      "%15" F_IIDP " : dovetail between spurs\n"
	      "%15" F_IIDP " : containment to non-singly contained\n"
	      "%15" F_IIDP " : containment touches spur and non-spur\n"
	      "%15" F_IIDP " : containment between spurs\n"
	      "      The intra-chunk skeleton overlaps:\n"
	      "%15" F_IIDP " : dovetail intra-chunk (between non-contained)\n"
	      "%15" F_IIDP " : dovetail non-contained touches singly contained\n"
	      "%15" F_IIDP " : containment to singly contained\n"
	      "      The marked overlaps:\n"
	      "%15" F_IIDP " : dovetail marked by branch points\n"
	      "      The removed overlaps:\n"
	      "%15" F_IIDP " : dovetail removed by transitivity\n"
	      "%15" F_IIDP " : dovetail removed by threshold\n"
	      "%15" F_IIDP " : containment removed by transitivity\n"
	      "%15" F_IIDP " : containment removed by containment\n"
	      "%15" F_IIDP " : duplicate dvt\n"
	      "%15" F_IIDP " : duplicate con\n"
              ,
	      nedge/2,
	      n_as_cgb_dovetail/2,
	      n_as_cgb_thickest/2,
	      /* The inter-chunk skeleton overlaps: */
	      n_as_cgb_interchunk/2,
	      n_as_cgb_touches_contained/2,
	      n_as_cgb_between_contained/2,
              n_as_cgb_touches_crappy_dvt/2,
              n_as_cgb_between_crappy_dvt/2,
	      n_as_cgb_containment_interchunk/2,
              n_as_cgb_touches_crappy_con/2,
              n_as_cgb_between_crappy_con/2,
              
	      /* The intra-chunk skeleton overlaps: */
	      n_as_cgb_intrachunk/2, 
	      n_as_cgb_touches_singly_contained/2,
	      n_as_cgb_containment_intrachunk/2,

              /* The marked overlaps: */
	      n_as_cgb_marked_by_branch_dvt/2,
	      /* The removed overlaps: */
	      n_as_cgb_removed_by_transitivity_dvt/2,
	      n_as_cgb_removed_by_threshold_dvt/2,
	      n_as_cgb_removed_by_transitivity_con/2,
	      n_as_cgb_removed_by_threshold_con/2,
	      n_as_cgb_removed_by_duplicate_dvt/2,
	      n_as_cgb_removed_by_duplicate_con/2
	      );
    }
#ifdef GENINFO
    { // 
      IntEdge_ID
	iedge,
	count_invalid_dovetail_edges=0,
	count_invalid_thickest_edges=0,
	count_invalid_interchunk_edges=0,
	count_invalid_intrachunk_edges=0,
	count_invalid_containment_edges=0,
	count_invalid_touches_contained_edges=0,
	count_invalid_between_contained_edges=0,
	count_invalid_touches_crappy_dvt=0,
	count_invalid_touches_crappy_con=0,
	count_invalid_between_crappy_dvt=0,
	count_invalid_between_crappy_con=0,
	count_invalid_marked_by_branch_dvt_edges=0,
	count_invalid_removed_by_transitivity_dvt_edges=0,
	count_invalid_removed_by_transitivity_con_edges=0,
	count_invalid_removed_by_threshold_dvt_edges=0,
	count_invalid_removed_by_threshold_con_edges=0,
	count_invalid_removed_by_duplicate_dvt_edges=0,
	count_invalid_removed_by_duplicate_con_edges=0;

      
      for(iedge=0;iedge<nedge;iedge++){
	const Tnes ines = get_nes_edge(edges,iedge);
#ifdef STORE_OVERLAP_EXTREMES
	const int iamn = get_amn_edge(edges,iedge);
	const int iamx = get_amx_edge(edges,iedge);
#endif // STORE_OVERLAP_EXTREMES

	const int invalid = get_inv_edge(edges,iedge);

	if( invalid ) {
	  switch(ines) {
	  case AS_CGB_DOVETAIL_EDGE:
	    count_invalid_dovetail_edges++; break;
	  case AS_CGB_THICKEST_EDGE:
	    count_invalid_thickest_edges++; break;
	  case AS_CGB_INTERCHUNK_EDGE:
	    count_invalid_interchunk_edges++; break;
	  case AS_CGB_INTRACHUNK_EDGE:
	    count_invalid_intrachunk_edges++; break;

	  case AS_CGB_TOUCHES_CONTAINED_EDGE:
	    count_invalid_touches_contained_edges++; break;
	  case AS_CGB_BETWEEN_CONTAINED_EDGE:
	    count_invalid_between_contained_edges++; break;
          case AS_CGB_TOUCHES_CRAPPY_DVT:
            count_invalid_touches_crappy_dvt++; break;
          case AS_CGB_BETWEEN_CRAPPY_DVT:
            count_invalid_between_crappy_dvt++; break;
	  case AS_CGB_MARKED_BY_BRANCH_DVT:
	    count_invalid_marked_by_branch_dvt_edges++; break;
          case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
            count_invalid_removed_by_transitivity_dvt_edges++; break;
	  case AS_CGB_REMOVED_BY_THRESHOLD_DVT:
	    count_invalid_removed_by_threshold_dvt_edges++; break;
	  case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
	    count_invalid_removed_by_duplicate_dvt_edges++; break;

	  case AS_CGB_CONTAINED_EDGE:
#ifdef DEBUGGING
          {
            const IntFragment_ID iavx = get_avx_edge(edges,iedge);
            const IntFragment_ID ibvx = get_bvx_edge(edges,iedge);
            const int iasx = get_asx_edge(edges,iedge);
            const int ibsx = get_bsx_edge(edges,iedge);

	    printf("invalid_containment_edges: %d,%d %d,%d\n",
		   get_iid_fragment(frags,iavx),iasx,
		   get_iid_fragment(frags,ibvx),ibsx);
          }
#endif // DEBUGGING
	    count_invalid_containment_edges++; break;

          case AS_CGB_TOUCHES_CRAPPY_CON:
            count_invalid_touches_crappy_con++; break;
          case AS_CGB_BETWEEN_CRAPPY_CON:
            count_invalid_between_crappy_con++; break;
          case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
	    count_invalid_removed_by_transitivity_con_edges++; break;
	  case AS_CGB_REMOVED_BY_THRESHOLD_CON:
	    count_invalid_removed_by_threshold_con_edges++; break;
	  case AS_CGB_REMOVED_BY_DUPLICATE_CON:
	    count_invalid_removed_by_duplicate_con_edges++; break;

	  default:
	    fprintf(stderr,"Unsupported edge label nes=%d\n",ines);
	    assert(FALSE);
	  }
	}
      }
      
      fprintf(fout,"Fragment information from the simulator\n");
      
      fprintf(fout,
	      "%15" F_IIDP " : number of invalid dovetail undistinguished overlaps.\n"
	      "%15" F_IIDP " : number of invalid thickest undistinguished overlaps.\n"
	      "%15" F_IIDP " : number of invalid inter-chunk overlaps.\n"
	      "%15" F_IIDP " : number of invalid intra-chunk overlaps.\n"
	      "%15" F_IIDP " : number of invalid containment overlaps.\n"
	      "%15" F_IIDP " : number of invalid touches contained overlaps.\n"
	      "%15" F_IIDP " : number of invalid between contained overlaps.\n"
	      "%15" F_IIDP " : number of invalid touches spur dvt overlaps.\n"
	      "%15" F_IIDP " : number of invalid between spur dvt overlaps.\n"
	      "%15" F_IIDP " : number of invalid touches spur toc overlaps.\n"
	      "%15" F_IIDP " : number of invalid between spur toc overlaps.\n"
	      "%15" F_IIDP " : number of invalid marked_by_branch point dvt.\n"
	      "%15" F_IIDP " : number of invalid removed_by_transitivity overlaps dvt.\n"
	      "%15" F_IIDP " : number of invalid removed_by_transitivity overlaps con.\n"
	      "%15" F_IIDP " : number of invalid removed_by_threshold overlaps dvt.\n"
	      "%15" F_IIDP " : number of invalid removed_by_threshold overlaps con.\n"
	      "%15" F_IIDP " : number of invalid removed_by_duplicate overlaps dvt.\n"
	      "%15" F_IIDP " : number of invalid removed_by_duplicate overlaps con.\n",
	      count_invalid_dovetail_edges/2,
	      count_invalid_thickest_edges/2,
	      count_invalid_interchunk_edges/2,
	      count_invalid_intrachunk_edges/2,
	      count_invalid_containment_edges/2,
	      count_invalid_touches_contained_edges/2,
	      count_invalid_between_contained_edges/2,
              count_invalid_touches_crappy_dvt/2,
              count_invalid_between_crappy_dvt/2,
              count_invalid_touches_crappy_con/2,
              count_invalid_between_crappy_con/2,
	      count_invalid_marked_by_branch_dvt_edges/2,
	      count_invalid_removed_by_transitivity_dvt_edges/2,
	      count_invalid_removed_by_transitivity_con_edges/2,
	      count_invalid_removed_by_threshold_dvt_edges/2,
	      count_invalid_removed_by_threshold_con_edges/2,
	      count_invalid_removed_by_duplicate_dvt_edges/2,
	      count_invalid_removed_by_duplicate_con_edges/2
	      );
    }
#endif // GENINFO
    
    
    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*solo_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*hanging_alone_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*hanging_chunk_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*hanging_crappy_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*thru_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*interchunk_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*intrachunk_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*orphanedcont_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*multicont_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*branchmulticont_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*singlecont_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE),
	*unplacedcont_histogram 
	= create_histogram(nsample,nbucket,TRUE,FALSE);

      for(ifrag=0;ifrag<nfrag;ifrag++) {
	const FragType type = get_typ_fragment(frags,ifrag);
	switch(get_lab_fragment(frags,ifrag)) {
	case AS_CGB_SOLO_FRAG:
	  add_to_histogram(solo_histogram, (int)type, NULL); break;
	case AS_CGB_HANGING_FRAG:
	  add_to_histogram(hanging_alone_histogram, (int)type, NULL); break;
	case AS_CGB_HANGING_CHUNK_FRAG:
	  add_to_histogram(hanging_chunk_histogram, (int)type, NULL); break;
	case AS_CGB_HANGING_CRAPPY_FRAG:
	  add_to_histogram(hanging_crappy_histogram, (int)type, NULL); break;
	case AS_CGB_THRU_FRAG:
	  add_to_histogram(thru_histogram, (int)type, NULL); break;
	case AS_CGB_INTERCHUNK_FRAG:
	  add_to_histogram(interchunk_histogram, (int)type, NULL); break;
	case AS_CGB_INTRACHUNK_FRAG:
	  add_to_histogram(intrachunk_histogram, (int)type, NULL); break;
	case AS_CGB_ORPHANEDCONT_FRAG:
	  add_to_histogram(orphanedcont_histogram, (int)type, NULL); break;
	case AS_CGB_MULTICONT_FRAG:
	  add_to_histogram(multicont_histogram, (int)type, NULL); break;
	case AS_CGB_BRANCHMULTICONT_FRAG:
	  add_to_histogram(branchmulticont_histogram, (int)type, NULL); break;
	case AS_CGB_SINGLECONT_FRAG:
	  add_to_histogram(singlecont_histogram, (int)type, NULL); break;
	case AS_CGB_UNPLACEDCONT_FRAG:
	  add_to_histogram(unplacedcont_histogram, (int)type, NULL); break;
	case AS_CGB_DELETED_FRAG:
	  break;
	default:
	  assert(FALSE);
	}
      }
      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "solo fragments.\n");
      print_histogram(fout,solo_histogram, 0, 1);
      free_histogram(solo_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "hanging alone fragments.\n");
      print_histogram(fout,hanging_alone_histogram, 0, 1);
      free_histogram(hanging_alone_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "hanging chunk fragments.\n");
      print_histogram(fout,hanging_chunk_histogram, 0, 1);
      free_histogram(hanging_chunk_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "hanging spur fragments.\n");
      print_histogram(fout,hanging_crappy_histogram, 0, 1);
      free_histogram(hanging_crappy_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "thru fragments.\n");
      print_histogram(fout,thru_histogram, 0, 1);
      free_histogram(thru_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "interchunk fragments.\n");
      print_histogram(fout,interchunk_histogram, 0, 1);
      free_histogram(interchunk_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "intrachunk fragments.\n");
      print_histogram(fout,intrachunk_histogram, 0, 1);
      free_histogram(intrachunk_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "orphaned contained fragments.\n");
      print_histogram(fout,orphanedcont_histogram, 0, 1);
      free_histogram(orphanedcont_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "unplaced contained fragments.\n");
      print_histogram(fout,unplacedcont_histogram, 0, 1);
      free_histogram(unplacedcont_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "singly contained fragments.\n");
      print_histogram(fout,singlecont_histogram, 0, 1);
      free_histogram(singlecont_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "multiply contained fragments.\n");
      print_histogram(fout,multicont_histogram, 0, 1);
      free_histogram(multicont_histogram);

      fprintf(fout,"\n\nHistogram of the fragment type of "
	      "branch multiply contained fragments.\n");
      print_histogram(fout,branchmulticont_histogram, 0, 1);
      free_histogram(branchmulticont_histogram);
    }
    
    {
      IntFragment_ID ifrag;
      int isuff;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*raw_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE),
	*work_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE),
	*reduced_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE),
	*dovetail_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE),
	*interchunk_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE),
	*intrachunk_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE),
	*singly_containment_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE),
	*multiply_containment_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE),
	*touches_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE),
	*between_edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE);

      for(ifrag=0;ifrag<nfrag;ifrag++) { 
	for(isuff=0;isuff<2;isuff++) {
	  int count_dovetail_edges=0;
	  int count_thickest_edges=0;
	  int count_interchunk_edges=0;
	  int count_intrachunk_edges=0;
	  int count_touches_contained_edges=0;
	  int count_between_contained_edges=0;
	  int count_singly_containment_edges=0;
	  int count_multiply_containment_edges=0;
	  int count_touches_crappy_dvt=0;
	  int count_touches_crappy_con=0;
	  int count_between_crappy_dvt=0;
	  int count_between_crappy_con=0;
          // int count_marked_by_containment_dvt_edges=0;
          // int count_marked_by_containment_con_edges=0;
	  int count_removed_by_transitivity_dvt_edges=0;
	  int count_removed_by_transitivity_con_edges=0;
	  // int count_removed_by_containment_dvt_edges=0;
	  // int count_removed_by_containment_con_edges=0;
	  int count_removed_by_threshold_dvt_edges=0;
	  int count_removed_by_threshold_con_edges=0;
	  int count_marked_by_branch_dvt_edges=0;
	  int count_removed_by_duplicate_dvt_edges=0;
	  int count_removed_by_duplicate_con_edges=0;
	  int count = 0;
	  { 
	    IntEdge_ID snode = get_segstart_vertex(frags,ifrag,isuff);
	    int nnode = get_seglen_vertex(frags,ifrag,isuff);
	    { IntEdge_ID ie; for(ie=snode;ie<snode+nnode;ie++) {
	      Tnes nes = get_nes_edge(edges,ie);
	      IntFragment_ID avx = get_avx_edge(edges,ie);
	      IntFragment_ID bvx = get_bvx_edge(edges,ie);
	      Tlab alab = get_lab_fragment(frags,avx);
	      Tlab blab = get_lab_fragment(frags,bvx);
	      switch(nes) {
	      case AS_CGB_DOVETAIL_EDGE:
		count_dovetail_edges++; break;
	      case AS_CGB_THICKEST_EDGE:
		count_thickest_edges++; break;
	      case AS_CGB_INTERCHUNK_EDGE:
		count_interchunk_edges++; break;
	      case AS_CGB_INTRACHUNK_EDGE:
		count_intrachunk_edges++; break;

	      case AS_CGB_TOUCHES_CONTAINED_EDGE:
		count_touches_contained_edges++; break;
	      case AS_CGB_BETWEEN_CONTAINED_EDGE:
		count_between_contained_edges++; break;
              case AS_CGB_TOUCHES_CRAPPY_DVT:
                count_touches_crappy_dvt++; break;
              case AS_CGB_BETWEEN_CRAPPY_DVT:
                count_between_crappy_dvt++; break;
	      case AS_CGB_MARKED_BY_BRANCH_DVT:
		count_marked_by_branch_dvt_edges++; break;
       	      case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
		count_removed_by_transitivity_dvt_edges++; break;
	      case AS_CGB_REMOVED_BY_THRESHOLD_DVT:
		count_removed_by_threshold_dvt_edges++; break;
	      case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
		count_removed_by_duplicate_dvt_edges++;


	      case AS_CGB_CONTAINED_EDGE:
		if((AS_CGB_SINGLECONT_FRAG == alab) ||
		   (AS_CGB_SINGLECONT_FRAG == blab) ) {
		  count_singly_containment_edges++; 
		} else {
		  count_multiply_containment_edges++; 
		}
		break;

              case AS_CGB_TOUCHES_CRAPPY_CON:
                count_touches_crappy_con++; break;
              case AS_CGB_BETWEEN_CRAPPY_CON:
                count_between_crappy_con++; break;
	      case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
		count_removed_by_transitivity_con_edges++; break;
	      case AS_CGB_REMOVED_BY_THRESHOLD_CON:
		count_removed_by_threshold_con_edges++; break;
	      case AS_CGB_REMOVED_BY_DUPLICATE_CON:
		count_removed_by_duplicate_con_edges++;
		break;

	      default:
		assert(FALSE);
	      }
	    }}
	  }

	  count =
	    count_dovetail_edges +
	    count_thickest_edges +
	    count_interchunk_edges +
	    count_intrachunk_edges +
	    count_singly_containment_edges +
	    count_multiply_containment_edges +
	    count_touches_contained_edges +
	    count_between_contained_edges +
            count_touches_crappy_dvt +
	    count_marked_by_branch_dvt_edges +
	    count_removed_by_transitivity_dvt_edges +
	    count_removed_by_transitivity_con_edges +
	    count_removed_by_threshold_dvt_edges +
	    count_removed_by_threshold_con_edges;
	  add_to_histogram(raw_edges_per_vertex_histogram, count, NULL);

	  count =
	    count_dovetail_edges +
	    count_thickest_edges +
	    count_interchunk_edges +
	    count_intrachunk_edges +
	    count_singly_containment_edges +
	    count_multiply_containment_edges +
	    count_touches_contained_edges +
	    count_between_contained_edges +
            count_touches_crappy_dvt +
	    count_marked_by_branch_dvt_edges;
	  add_to_histogram(work_edges_per_vertex_histogram, count, NULL);

	  count =
	    count_dovetail_edges +
	    count_thickest_edges +
	    count_interchunk_edges +
	    count_intrachunk_edges +
	    count_singly_containment_edges +
	    count_multiply_containment_edges +
	    count_touches_contained_edges +
	    count_between_contained_edges +
            count_touches_crappy_dvt;

	  add_to_histogram(reduced_edges_per_vertex_histogram, count, NULL);

	  count = count_dovetail_edges;
	  add_to_histogram(dovetail_edges_per_vertex_histogram, count, NULL);

	  count = count_thickest_edges;
	  add_to_histogram(dovetail_edges_per_vertex_histogram, count, NULL);

	  count = count_interchunk_edges;
	  add_to_histogram(interchunk_edges_per_vertex_histogram, count, NULL);

	  count = count_touches_contained_edges;
	  add_to_histogram(touches_edges_per_vertex_histogram, count, NULL);

	  count = count_between_contained_edges;
	  add_to_histogram(between_edges_per_vertex_histogram, count, NULL);

	  count = count_multiply_containment_edges;
	  add_to_histogram(multiply_containment_edges_per_vertex_histogram,
			   count, NULL);

	  count = count_intrachunk_edges;
	  add_to_histogram(intrachunk_edges_per_vertex_histogram, count, NULL);

	  count = count_singly_containment_edges;
	  add_to_histogram(singly_containment_edges_per_vertex_histogram,
			   count, NULL);

	}
      }

      fprintf(fout,"\n\nHistogram of the raw degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,raw_edges_per_vertex_histogram, 0, 1);
      free_histogram(raw_edges_per_vertex_histogram);

      fprintf(fout,"\n\nHistogram of the working degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,work_edges_per_vertex_histogram, 0, 1);
      free_histogram(work_edges_per_vertex_histogram);

      fprintf(fout,"\n\nHistogram of the reduced degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,reduced_edges_per_vertex_histogram, 0, 1);
      free_histogram(reduced_edges_per_vertex_histogram);

      fprintf(fout,"\n\nHistogram of the non-chunking dovetail degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,dovetail_edges_per_vertex_histogram, 0, 1);
      free_histogram(dovetail_edges_per_vertex_histogram);

      fprintf(fout,"\n\nHistogram of the interchunk degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,interchunk_edges_per_vertex_histogram, 0, 1);
      free_histogram(interchunk_edges_per_vertex_histogram);

      fprintf(fout,"\n\nHistogram of the touches-contained degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,touches_edges_per_vertex_histogram, 0, 1);
      free_histogram(touches_edges_per_vertex_histogram);

      fprintf(fout,"\n\nHistogram of the between-contained degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,between_edges_per_vertex_histogram, 0, 1);
      free_histogram(between_edges_per_vertex_histogram);

      fprintf(fout,"\n\nHistogram of the multiply contained degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,multiply_containment_edges_per_vertex_histogram, 0, 1);
      free_histogram(multiply_containment_edges_per_vertex_histogram);


      fprintf(fout,"\n\nHistogram of the intrachunk degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,intrachunk_edges_per_vertex_histogram, 0, 1);
      free_histogram(intrachunk_edges_per_vertex_histogram);

      fprintf(fout,"\n\nHistogram of the singly_contained degree of the fragment-ends"
	      " in the fragment overlap graph\n");
      print_histogram(fout,singly_containment_edges_per_vertex_histogram, 0, 1);
      free_histogram(singly_containment_edges_per_vertex_histogram);

    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram 
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_SOLO_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for solo fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram 
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_HANGING_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for hanging fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram 
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_THRU_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for thru fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_ORPHANEDCONT_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for orphaned fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_MULTICONT_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for multiply contained fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_BRANCHMULTICONT_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for branch multiply contained fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_HANGING_CRAPPY_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for hanging crappy fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_INTERCHUNK_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for interchunk fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_INTRACHUNK_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for intrachunk fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_SINGLECONT_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for uniquely contained fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }

    {
      IntFragment_ID ifrag;
      const int nsample=500;
      const int nbucket=500;
      Histogram_t 
	*edges_per_vertex_histogram
	= create_histogram(nsample,nbucket,0,TRUE);
      for(ifrag=0;ifrag<nfrag;ifrag++) {
	if( AS_CGB_DELETED_FRAG == get_lab_fragment(frags,ifrag)) {
	  int nnode;
	  nnode = get_seglen_vertex(frags,ifrag,FALSE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  nnode = get_seglen_vertex(frags,ifrag,TRUE);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
      }
      fprintf(fout,"\n\nHistogram of the degree of the fragment-ends "
	      "in the fragment overlap graph\n"
	      "for deleted fragments.\n");
      print_histogram(fout,edges_per_vertex_histogram, 0, 1);
      free_histogram(edges_per_vertex_histogram);
    }
}

static void analyze_the_chunks
(
 FILE *fout,
 FILE *fp_unitig_statistics,
 const IntFragment_ID max_frag_iid,
 Tfragment frags[],
 Tedge edges[],
 ChunkAnnotation chunkinfo[],
 TChunkFrag chunkfrags[],
 TChunkMesg thechunks[],
 const BPTYPE nbase_in_genome,
 const float cgb_unique_cutoff,
 const float global_fragment_arrival_rate )
{
  IntChunk_ID ichunk;
  IntFragment_ID num_of_chunks[MAX_NUM_CHUNK_LABELS]={0};
  IntFragment_ID nfrag_in_all_chunks=0;
  IntFragment_ID nfrag_essential_in_all_chunks=0;
  IntFragment_ID nfrag_contained_in_all_chunks=0;
  BPTYPE nbase_essential_in_all_chunks=0;
  BPTYPE rho_in_all_chunks = 0;
  IntFragment_ID n_rs_frag_in_all_chunks = 0;
  IntFragment_ID n_nr_frag_in_all_chunks = 0;

#ifdef GENINFO
  BPTYPE
    rho_in_true_unique=0, rho_in_true_repeat=0,
    rho_in_true_unique_with_cpos=0, rho_in_true_repeat_with_cpos=0,
    rho_in_true_unique_with_cneg=0, rho_in_true_repeat_with_cneg=0;
#endif

  const int nsample=500;
  const int nbucket=500;
  MyHistoDataType zork;

  IntFragment_ID * fragment_visited = NULL;
  int * fragment_timesinchunks = NULL;
  FragmentHashObject  * afr_to_avx = NULL;
  const IntFragment_ID  nfrag   = GetNumFragments(frags);
  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  Histogram_t * length_of_unitigs_histogram
    = create_histogram(nsample,nbucket,0,TRUE);


#ifdef GENINFO 
  Histogram_t * true_uniques = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * true_repeats = create_histogram(nsample,nbucket,0,TRUE);
#endif
  Histogram_t * rho_histogram 
    = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * coverage_histogram 
    = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * nfrag_in_chunk_histogram 
    = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * nfrag_essential_in_chunk_histogram 
    = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * nbase_essential_in_chunk_histogram 
    = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * fragment_timesinchunks_histogram 
    = create_histogram(nsample,nbucket,0,TRUE);

  Histogram_t * labeled_unitig_histogram[MAX_NUM_CHUNK_LABELS];

  {
    int ii;
    for(ii=0;ii<MAX_NUM_CHUNK_LABELS;ii++){
      labeled_unitig_histogram[ii]
	= create_histogram(nsample,nbucket,0,TRUE);
      extend_histogram
	(labeled_unitig_histogram[ii],sizeof(MyHistoDataType),
	 myindexdata,mysetdata,myaggregate,myprintdata);
    }
  }

  extend_histogram(length_of_unitigs_histogram, sizeof(MyHistoDataType),
		   myindexdata,mysetdata,myaggregate,myprintdata);
#ifdef GENINFO
  extend_histogram(true_uniques,sizeof(MyHistoDataType),
		   myindexdata,mysetdata,myaggregate,myprintdata);
  extend_histogram(true_repeats,sizeof(MyHistoDataType),
		   myindexdata,mysetdata,myaggregate,myprintdata);
#endif /*GENINFO*/

  fragment_visited = (IntFragment_ID *)safe_malloc(nfrag*sizeof(IntFragment_ID));
  afr_to_avx = create_FragmentHash((max_frag_iid+1));
  fragment_timesinchunks = (int *)safe_malloc(nfrag*sizeof(int));

  assert(fragment_visited != NULL);
  assert(fragment_timesinchunks != NULL);
  assert(afr_to_avx != NULL);

  /* Initialize a flag for chunk following. */
  {
    IntFragment_ID ifrag;
    for(ifrag=0;ifrag<nfrag;ifrag++) { 
      const IntFragment_ID iid = get_iid_fragment(frags,ifrag);
      fragment_visited[ifrag] = FRAGMENT_NOT_VISITED;
      fragment_timesinchunks[ifrag] = 0;
      set_vid_FragmentHash(afr_to_avx,iid,ifrag);
    }
  }
  
  assert((!0) == 1); /* Needed for the following bitwise XOR operator. */
  for(ichunk=0;ichunk<nchunks;ichunk++) {

    const IntFragment_ID irec_start_of_chunk 
      = GetVA_AChunkMesg(thechunks,ichunk)->f_list;
    const BPTYPE rho 
      = GetVA_AChunkMesg(thechunks,ichunk)->rho;
    const IntFragment_ID nfrag_in_chunk 
      = GetVA_AChunkMesg(thechunks,ichunk)->num_frags;
    const BPTYPE nbase_essential_in_chunk 
      = GetVA_AChunkMesg(thechunks,ichunk)->bp_length;

    const int number_of_randomly_sampled_fragments_in_chunk
      = count_the_randomly_sampled_fragments_in_a_chunk
      ( frags, chunkfrags, thechunks, ichunk);
    const float coverage_statistic = compute_coverage_statistic
      ( rho,
        number_of_randomly_sampled_fragments_in_chunk,
        global_fragment_arrival_rate );
    const int number_of_non_randomly_sampled_fragments_in_chunk = 
      nfrag_in_chunk - number_of_randomly_sampled_fragments_in_chunk;

    /* The coverage statistic should be positive for single coverage,
       negative for multiple coverage, and near zero for indecisive. */
    const float coverage_resolution = 1.f;
    const int coverage_index 
      = (int)(coverage_statistic/coverage_resolution);
    const int arrival_distance =
      (nfrag_in_chunk > 1 ? rho/(nfrag_in_chunk-1) : 0);
    // The fragment arrival distance in base pairs.
    // CMM: Should this be for Celera reads only?

    IntFragment_ID nfrag_essential_in_chunk=0;
    IntFragment_ID nfrag_contained_in_chunk=0;
    BPTYPE nbase_sampled_in_chunk=0;
    BPTYPE nbase_essential_sampled_in_chunk=0;
    BPTYPE nbase_contained_sampled_in_chunk=0;

    ChunkLabel chunk_label;

#ifdef DEBUG07
    fprintf(stderr,
	    "Process ichunk,nchunks,nfrag_in_chunk,nbase_essential_in_chunk=\n"
	   F_IID "," F_IID "," F_IID "," BPFORMAT "\n",
	   ichunk,nchunks,nfrag_in_chunk,nbase_essential_in_chunk);
#endif

    {
      // Process the chunk-end fragments first to label the chunks.
      const IntFragment_ID chunk_avx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_avx;
      const IntFragment_ID chunk_bvx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_bvx;
      const int chunk_asx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_asx;
      const int chunk_bsx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_bsx;

      assert( nfrag_in_chunk > 0);
      if(chunk_avx == chunk_bvx ) {
	// This is a spanned chunk.
	assert(chunk_asx != chunk_bsx);
	if( nfrag_in_chunk == 1) {
	  // This is a singleton chunk.
	  const Tlab lab=get_lab_fragment(frags,chunk_avx);
	  if( (AS_CGB_MULTICONT_FRAG == lab) ||
              (AS_CGB_BRANCHMULTICONT_FRAG == lab) ||
	      (AS_CGB_ORPHANEDCONT_FRAG == lab) ) {
	    chunk_label = CGA_SINGLETON_CONTAINED;
	  } else {
	    chunk_label = CGA_SINGLETON_NONCONTAINED;
	  }
	} else {
	  chunk_label = CGA_NONSINGLETON_SPANNED;
	}
      } else {
	// This is a non-spanned chunk.
	assert(nfrag_in_chunk > 1);
	chunk_label = CGA_NONSINGLETON_NONSPANNED;
      }
    }
    num_of_chunks[chunk_label] ++;

    // Process the fragments of the chunk.
    {
      IntFragment_ID ifrag;
      for(ifrag=0;ifrag<nfrag_in_chunk;ifrag++){
	
      const IntFragment_ID ivc = irec_start_of_chunk + ifrag;
      const IntFragment_ID vid = GetVA_AChunkFrag(chunkfrags,ivc)->vid;

      const IntFragment_ID iid  = get_iid_fragment(frags,vid);
      const Tlab ilabel = get_lab_fragment(frags,vid);
      const int ilen = get_length_fragment(frags,vid);
      const IntFragment_ID ibvx = get_vid_FragmentHash(afr_to_avx,iid);

      fragment_visited[ibvx] = ichunk;
      fragment_timesinchunks[ibvx] ++;
      
      switch(ilabel){
      case AS_CGB_SOLO_FRAG:
      case AS_CGB_HANGING_FRAG:
      case AS_CGB_THRU_FRAG:
      case AS_CGB_MULTICONT_FRAG:
      case AS_CGB_BRANCHMULTICONT_FRAG:
      case AS_CGB_HANGING_CHUNK_FRAG:
      case AS_CGB_HANGING_CRAPPY_FRAG:
      case AS_CGB_ORPHANEDCONT_FRAG:
      case AS_CGB_INTERCHUNK_FRAG:
      case AS_CGB_INTRACHUNK_FRAG:
	{
	  nfrag_essential_in_chunk++;
	  nbase_essential_sampled_in_chunk += ilen;
	}
	break;
      case AS_CGB_SINGLECONT_FRAG:
	{
	  assert(TRUE == get_con_fragment(frags,vid));
	  nfrag_contained_in_chunk++;
	  nbase_contained_sampled_in_chunk += ilen;
	} 
	break;
      default:
	assert(FALSE);
      }
    }
    }
    
    assert(nfrag_in_chunk ==
	   nfrag_essential_in_chunk + nfrag_contained_in_chunk);
    nbase_sampled_in_chunk = nbase_essential_sampled_in_chunk
      + nbase_contained_sampled_in_chunk;
      
    nfrag_in_all_chunks += nfrag_in_chunk;
    n_rs_frag_in_all_chunks += number_of_randomly_sampled_fragments_in_chunk;
    n_nr_frag_in_all_chunks += number_of_non_randomly_sampled_fragments_in_chunk;

    assert(nfrag_in_all_chunks == 
	   n_rs_frag_in_all_chunks + n_nr_frag_in_all_chunks);
    nfrag_essential_in_all_chunks += nfrag_essential_in_chunk;
    nfrag_contained_in_all_chunks += nfrag_contained_in_chunk;
    nbase_essential_in_all_chunks += nbase_essential_in_chunk;
    rho_in_all_chunks += rho;

    zork.nsamples  = 1;
    zork.sum_frags = nfrag_in_chunk;
    zork.min_frags = nfrag_in_chunk;
    zork.max_frags = nfrag_in_chunk;
    zork.sum_rs_frags = number_of_randomly_sampled_fragments_in_chunk;
    zork.min_rs_frags = number_of_randomly_sampled_fragments_in_chunk;
    zork.max_rs_frags = number_of_randomly_sampled_fragments_in_chunk;
    zork.sum_nr_frags = number_of_non_randomly_sampled_fragments_in_chunk;
    zork.min_nr_frags = number_of_non_randomly_sampled_fragments_in_chunk;
    zork.max_nr_frags = number_of_non_randomly_sampled_fragments_in_chunk;
    zork.sum_bp    = nbase_essential_in_chunk;
    zork.min_bp    = nbase_essential_in_chunk;
    zork.max_bp    = nbase_essential_in_chunk;
    zork.sum_rho   = rho;
    zork.min_rho   = rho;
    zork.max_rho   = rho;
    zork.sum_arrival = arrival_distance;
    zork.min_arrival = arrival_distance;
    zork.max_arrival = arrival_distance;
    zork.sum_discr = coverage_index;
    zork.min_discr = coverage_index;
    zork.max_discr = coverage_index;

    add_to_histogram(length_of_unitigs_histogram,
                     nbase_essential_in_chunk,&zork);

#ifdef GENINFO

    if( ! (chunkinfo[ichunk].essential_type == 'r') ) { 
      add_to_histogram(true_uniques, coverage_index, &zork);
      rho_in_true_unique += rho;
      if( coverage_statistic >= cgb_unique_cutoff ) {
	rho_in_true_unique_with_cpos += rho;
      } else {
	rho_in_true_unique_with_cneg += rho;
      }
    }
    if( (chunkinfo[ichunk].essential_type == 'r') ) { 
      add_to_histogram(true_repeats, coverage_index, &zork);
      rho_in_true_repeat += rho;
      if( coverage_statistic >= cgb_unique_cutoff ) {
	rho_in_true_repeat_with_cpos += rho;
      } else {
	rho_in_true_repeat_with_cneg += rho;
      }
    }
#endif /*GENINFO*/
    add_to_histogram(labeled_unitig_histogram[chunk_label],
		     coverage_index, &zork);

    { // For Gene Myer^s Jan 2000 paper:
      int 
	num_as_overlap[2]={0},                   // O
	num_as_touches_contained_overlap[2]={0}, // M
	num_as_between_contained_overlap[2]={0}, // Y
	num_as_1_contains_2_overlap[2]={0},      // C
	num_as_2_contains_1_overlap[2]={0};      // I
      
      // Process the chunk-end fragments first to label the chunks.
      const IntFragment_ID chunk_avx
        = GetVA_AChunkMesg(thechunks,ichunk)->chunk_avx;
      const IntFragment_ID chunk_bvx
        = GetVA_AChunkMesg(thechunks,ichunk)->chunk_bvx;
      const int chunk_asx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_asx;
      const int chunk_bsx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_bsx;
      const int chunk_spanned = (chunk_avx == chunk_bvx);
      int chunk_contained = FALSE;

      assert( nfrag_in_chunk > 0);
      if(chunk_avx == chunk_bvx ) {
	// This is a spanned chunk.
	assert(chunk_asx != chunk_bsx);
	if( nfrag_in_chunk == 1) {
	  // This is a singleton chunk.
	  const Tlab lab=get_lab_fragment(frags,chunk_avx);
	  if( (AS_CGB_MULTICONT_FRAG == lab) ||
              (AS_CGB_BRANCHMULTICONT_FRAG == lab) ||
	      (AS_CGB_ORPHANEDCONT_FRAG == lab) ) {
	    chunk_contained = TRUE;
	  }
	}
      }
      assert( (!chunk_contained)||(chunk_spanned) );

      { 
	int isuffix;
	for(isuffix=0;isuffix<2;isuffix++) {
	  {
	    const AChunkMesg * chunk = GetVA_AChunkMesg( thechunks, ichunk);
	    IntFragment_ID ifrag = (isuffix == 0 ? chunk->chunk_avx : chunk->chunk_bvx);
	    int        isuff = (isuffix == 0 ? chunk->chunk_asx : chunk->chunk_bsx);
	    IntEdge_ID ir;
	    const IntEdge_ID ir0  = get_segstart_vertex(frags,ifrag,isuff);
	    const int       nnode = get_seglen_vertex(frags,ifrag,isuff);
	    for(ir=ir0; (ir<ir0+nnode); ir++) {
	      Tnes nes = get_nes_edge(edges,ir);
	      switch(nes) {
	      case AS_CGB_INTERCHUNK_EDGE:
	      case AS_CGB_MARKED_BY_BRANCH_DVT:
              case AS_CGB_TOUCHES_CRAPPY_DVT:
              case AS_CGB_BETWEEN_CRAPPY_DVT:
		num_as_overlap[isuffix]++; break;

              case AS_CGB_INTRACHUNK_EDGE:
              case AS_CGB_DOVETAIL_EDGE:
              case AS_CGB_THICKEST_EDGE:
#if 0
                fprintf(stderr,"Unexpected overlap label nes=%d interior to chunk\n", nes);
                fprintf(stderr,"ichunk=" F_IID " isuffix=%d ifrag=" F_IID " isuff=%d\n",
                        ichunk, isuffix, ifrag, isuff);
                // Should this edge type be seen here??
		//assert(FALSE);
#endif
                break;
	      case AS_CGB_TOUCHES_CONTAINED_EDGE:
		num_as_touches_contained_overlap[isuffix]++; break;
              case AS_CGB_BETWEEN_CONTAINED_EDGE:
		num_as_between_contained_overlap[isuffix]++; break;
	      case AS_CGB_CONTAINED_EDGE:
	      case AS_CGB_TOUCHES_CRAPPY_CON:
              case AS_CGB_BETWEEN_CRAPPY_CON:
                if(is_a_frc_edge(edges,ir)) {
                  num_as_2_contains_1_overlap[isuffix]++;
                } else {
                  num_as_1_contains_2_overlap[isuffix]++;
                }
		break;

	      default:
                fprintf(stderr,"Unexpected overlap label nes=%d\n", nes);
		assert(FALSE);
		break;
	      }
	    }
	  }
	}
      }

      if(NULL != fp_unitig_statistics) {
      if(ichunk==0) {
	fprintf(fp_unitig_statistics,
		"# For each row:\n"
		"#    chunk contained spanned nfrag_in_chunk \n"
		"#    nbase_essential_in_chunk rho arrival_distance coverage_index\n"
                "#    essential_valid_type(ur@) contained_valid_type(ur@)\n"
		"#    prefix_O prefix_M prefix_C prefix_I prefix_Y\n"
		"#    suffix_O suffix_M suffix_C suffix_I suffix_Y\n"
		);
      }
      fprintf(fp_unitig_statistics,
	      F_IID " %d %d " F_IID " " BPFORMAT " %d %d %d "
              "%c %c "
              "%d %d %d %d %d %d %d %d %d %d\n", 
	      ichunk,            // An identity field
	      (chunk_contained), // A selection field
	      (chunk_spanned),
	      nfrag_in_chunk, // A data field 
	      nbase_essential_in_chunk,
	      (int)rho,
	      arrival_distance,
	      coverage_index, // A selection field

              chunkinfo[ichunk].essential_type,
              chunkinfo[ichunk].contained_type,

	      num_as_overlap[0],                   // O
	      num_as_touches_contained_overlap[0], // M
	      num_as_1_contains_2_overlap[0],      // C
	      num_as_2_contains_1_overlap[0],      // I
	      num_as_between_contained_overlap[0], // Y
	      num_as_overlap[1],                   // O
	      num_as_touches_contained_overlap[1], // M
	      num_as_1_contains_2_overlap[1],      // C
	      num_as_2_contains_1_overlap[1],      // I
	      num_as_between_contained_overlap[1]  // Y
	      );
      }

    add_to_histogram(nfrag_in_chunk_histogram,
		     nfrag_in_chunk, NULL);
    add_to_histogram(nfrag_essential_in_chunk_histogram,
		     nfrag_essential_in_chunk, NULL);
    add_to_histogram(nbase_essential_in_chunk_histogram,
		     nbase_essential_in_chunk, NULL);
    add_to_histogram(rho_histogram,(int)rho, NULL);
    add_to_histogram(coverage_histogram,
		     (int)(nbase_essential_in_chunk/nfrag_in_chunk), NULL);
    
  }
  }

  assert(NULL != fout);
  {
    fprintf(fout,"\n\nUNITIG OVERLAP GRAPH INFORMATION\n\n");
    fprintf(fout,FRAG_FORMAT
	    " : Total number of unitigs\n",nchunks);
    {
      int ii;
      for(ii=0;ii<MAX_NUM_CHUNK_LABELS;ii++){
	switch(ii) {
	case CGA_SINGLETON_CONTAINED:
	  fprintf(fout,FRAG_FORMAT
		  " : Total number of singleton, contained unitigs\n",
		  num_of_chunks[ii]);
	  break;
	case CGA_SINGLETON_NONCONTAINED:
	  fprintf(fout,FRAG_FORMAT
		  " : Total number of singleton, non-contained unitigs\n",
		  num_of_chunks[ii]);
	  break;
	case CGA_NONSINGLETON_SPANNED:
	  fprintf(fout,FRAG_FORMAT
		  " : Total number of non-singleton, spanned unitigs\n",
		  num_of_chunks[ii]);
	  break;
	case CGA_NONSINGLETON_NONSPANNED:
	  fprintf(fout,FRAG_FORMAT
		  " : Total number of non-singleton, non-spanned unitigs\n",
		  num_of_chunks[ii]);
	  break;
	default:
	  break;
	}
      }
    }

    fprintf(fout,FRAG_FORMAT
	    " : Total number of fragments\n",
	    nfrag);
    fprintf(fout,FRAG_FORMAT
	    " : Total number of fragments in all unitigs\n",
	    nfrag_in_all_chunks);
    fprintf(fout,FRAG_FORMAT
	    " : Total number of essential fragments in all unitigs\n",
	    nfrag_essential_in_all_chunks);
    fprintf(fout,FRAG_FORMAT
	    " : Total number of contained fragments in all unitigs\n",
	    nfrag_contained_in_all_chunks);
    fprintf(fout,"% 15.10f"
	    " : Randomly sampled fragment arrival rate per bp\n",
	    global_fragment_arrival_rate);
    fprintf(fout,BPFORMAT15
	    " : The sum of overhangs in all the unitigs\n",
	    rho_in_all_chunks);
    fprintf(fout,BPFORMAT15
	    " : Total number of bases in all unitigs\n",
	    nbase_essential_in_all_chunks);
    fprintf(fout,BPFORMAT15
	    " : Estimated number of base pairs in the genome.\n", 
	    nbase_in_genome);

#ifdef GENINFO
    fprintf(fout,BPFORMAT15
	    " : The sum of overhangs in the truly unique unitigs\n",
	    rho_in_true_unique);
    fprintf(fout,BPFORMAT15
	    " : The sum of overhangs in the truly unique unitigs\n"
	    "                  with the Unique/Repeat discriminator score"
	    " >= %5.1f\n",
	    rho_in_true_unique_with_cpos, cgb_unique_cutoff);
    fprintf(fout,BPFORMAT15
	    " : The sum of overhangs in the truly unique unitigs\n"
	    "                  with the Unique/Repeat discriminator score"
	    " <  %5.1f\n",
	    rho_in_true_unique_with_cneg, cgb_unique_cutoff);
    
    /* fprintf(fout,"\n\n\n");*/
    fprintf(fout,BPFORMAT15
	    " : The sum of overhangs in the truly repeat unitigs\n",
	    rho_in_true_repeat);
    fprintf(fout,BPFORMAT15
	    " : The sum of overhangs in the truly repeat unitigs\n"
	    "                  with the Unique/Repeat discriminator score"
	    " >= %5.1f\n",
	    rho_in_true_repeat_with_cpos, cgb_unique_cutoff);
    fprintf(fout,BPFORMAT15
	    " : The sum of overhangs in the truly repeat unitigs\n"
	    "                  with the Unique/Repeat discriminator score"
	    " <  %5.1f\n",
	    rho_in_true_repeat_with_cneg, cgb_unique_cutoff);
#endif /*GENINFO*/

    {
      IntFragment_ID nfound = 0;
      IntFragment_ID ifrag;
      for(ifrag=0;ifrag<nfrag;ifrag++) /* beginning loop through fragments */ {
	Tlab lab = get_lab_fragment(frags,ifrag);
	if( (fragment_visited[ifrag] == FRAGMENT_NOT_VISITED) &&
	    (AS_CGB_DELETED_FRAG != lab) ) {
	  // #ifdef DEBUG10
	  fprintf(stderr,"Not visited: fragment iid=" F_IID ",ifrag=" F_IID ",ilab=%d\n",
		  get_iid_fragment(frags,ifrag),ifrag,lab);
	  // #endif
	  nfound++;
	}
      }
      fprintf(fout,"%15" F_IIDP " : Total number of contained fragments not connected\n"
	      "                  by containment edges to essential fragments.\n",
	      nfound);
    }

    {
      int unitigger_bpt = 0;
      for(ichunk=0;ichunk<nchunks; ichunk++) {
        int branching, position;
	
        branching = GetVA_AChunkMesg(thechunks,ichunk)->a_branch_type;
        position  = GetVA_AChunkMesg(thechunks,ichunk)->a_branch_point;
        assert(branching != AS_INTO_REPEAT);
        unitigger_bpt += (branching == AS_INTO_UNIQUE);
        
        branching = GetVA_AChunkMesg(thechunks,ichunk)->b_branch_type;
        position  = GetVA_AChunkMesg(thechunks,ichunk)->b_branch_point;
        assert(branching != AS_INTO_REPEAT);
        unitigger_bpt += (branching == AS_INTO_UNIQUE);
        
      }
      fprintf(fout,"%15d : Total number of branch points.\n",
              unitigger_bpt);
    }

    {
      compute_the_global_fragment_arrival_rate
        ( fout, nbase_in_genome, frags, edges,
          global_fragment_arrival_rate, chunkfrags, thechunks );
    }
    
    {
      IntFragment_ID ifrag;
      for(ifrag=0; ifrag<nfrag; ifrag++) {
	const Tlab ilab = get_lab_fragment(frags,ifrag);
	assert(!(ilab==AS_CGB_SOLO_FRAG) 
	       ||(fragment_timesinchunks[ifrag]==1));
	assert(!(ilab==AS_CGB_HANGING_FRAG) 
	       ||(fragment_timesinchunks[ifrag]==1));
	assert(!(ilab==AS_CGB_THRU_FRAG) 
	       ||(fragment_timesinchunks[ifrag]==1));
#if 1
	if(!(!(ilab==AS_CGB_INTERCHUNK_FRAG)
	     ||(fragment_timesinchunks[ifrag]==1))) {
	  fprintf(stderr,"** iid=" F_IID " ifrag=" F_IID " ilab=%d fragment_timesinchunks=%d\n",
		  get_iid_fragment(frags,ifrag),
		  ifrag, ilab,
		  fragment_timesinchunks[ifrag]
		  );
	}
#endif
	assert(!(ilab==AS_CGB_INTERCHUNK_FRAG)
	       ||(fragment_timesinchunks[ifrag]==1));
	if(!(!(ilab==AS_CGB_INTRACHUNK_FRAG)
	     ||(fragment_timesinchunks[ifrag]==1))) {
	  fprintf(stderr,"Not visited? ifrag,ilab,"
		  "fragment_timesinchunks[ifrag]="
		  F_IID ",%d,%d\n",
		  ifrag,ilab,fragment_timesinchunks[ifrag]);
	}
	assert(!(ilab==AS_CGB_INTRACHUNK_FRAG)
	       ||(fragment_timesinchunks[ifrag]==1)
	       ||(fragment_timesinchunks[ifrag]==0));
	if(AS_CGB_SINGLECONT_FRAG == ilab) {
	  assert(TRUE == get_con_fragment(frags,ifrag));
	  add_to_histogram(fragment_timesinchunks_histogram,
			   fragment_timesinchunks[ifrag],NULL);
	}
      }
    }

    fprintf(fout,"\n\nHistogram of the number of "
	    "fragment reads in a chunk.\n");
    print_histogram(fout,nfrag_in_chunk_histogram, 0, 1);

    fprintf(fout,"\n\nHistogram of the number of "
	    "non-contained fragment reads in a chunk.\n");
    print_histogram(fout,nfrag_essential_in_chunk_histogram, 0, 1);

    fprintf(fout,"\n\nHistogram of the number of "
	    "base pairs in a chunk\n");
    print_histogram(fout,nbase_essential_in_chunk_histogram, 0, 1);

    fprintf(fout,"\n\nHistogram of the sum of overhangs for chunks\n");
    print_histogram(fout,rho_histogram, 0, 1);
    
    fprintf(fout,"\n\nHistogram of the average bps per fragment for chunks\n");
    print_histogram(fout,coverage_histogram, 0, 1);
    
    fprintf(fout,"\n\nHistogram of the number of copies of"
	    " a particular contained fragment.\n");
    print_histogram(fout,fragment_timesinchunks_histogram,0,1);



    /////////////////////////////////////////////////////////////


    fprintf(fout,"\n\nUnitig Length Extended Histograms\n");
    fprintf(fout,"  Legend for extended histograms\n"
	    "length\t:  count\n"
	    "frags\tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    "... again for randomly sampled fragments\n"
	    "... again for non-randomly sampled fragments\n"
	    "bases\tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    "rho  \tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    "arrival\tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    "discr\tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    );

    print_histogram(fout,length_of_unitigs_histogram, 0, 1);

    /////////////////////////////////////////////////////////////

    fprintf(fout,"\n\nUnique/Repeat Extended Histograms\n");
    fprintf(fout,"  Legend for extended histograms\n"
	    "score\t:  count\n"
	    "frags\tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    "... again for randomly sampled fragments\n"
	    "... again for non-randomly sampled fragments\n"
	    "bases\tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    "rho  \tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    "arrival\tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    "discr\tsum\tcummulative\tcummulative   min  average  max\n"
	    "     \t   \t sum       \t fraction\n"
	    );


#ifdef GENINFO
    fprintf(fout,"\n\nHistogram of the Unique/Repeat discriminator score"
	    " for valid unitigs\n");
    print_histogram(fout,true_uniques, 0, 1);
    fprintf(fout,"\n\nHistogram of the Unique/Repeat discriminator score"
	    " for invalid unitigs\n");
    print_histogram(fout,true_repeats, 0, 1);

    {
      int ii;
      for(ii=0;ii<MAX_NUM_CHUNK_LABELS;ii++){
	fprintf(fout,"\n\nHistogram of the "
		"Unique/Repeat discriminator score for %s\n",
		ChunkLabelDesc[ii]);
	print_histogram(fout,labeled_unitig_histogram[ii],0,1);
      }
    }

#endif /*GENINFO*/

#ifdef NEVER
    fprintf(fout,"\n\nHistogram of fragment_timesinchunks\n");
    print_histogram(fout,fragment_timesinchunks_histogram, 0, 1);
#endif

#ifdef NEVER
    fprintf(fout,"\n\nHistogram of the degree of the ends of the chunks\n");
    print_histogram(fout,ndegree_of_chunkends_histogram, 0, 1);

    fprintf(fout,"\n\nHistogram of the degree squared of the ends of the chunks\n");
    print_histogram(fout,ndegreesq_of_chunkends_histogram, 0, 1);

    fprintf(fout,"\n\nHistogram of the degree of the ends of the "
	    "discriminator unique chunks\n");
    print_histogram(fout,ndegree_of_disc_unique_chunkends_histogram, 0, 1);

    fprintf(fout,"\n\nHistogram of the unique-unique degree of the ends of the "
	    "discriminator unique chunks\n");
    print_histogram(fout,ndegree_of_disc_unique_to_unique_histogram, 0, 1);

#ifdef GENINFO
    fprintf(fout,"\n\nHistogram of the degree of the ends of the "
	    "true unique chunks\n");
    print_histogram(fout,ndegree_of_gen_unique_chunkends_histogram, 0, 1);
    fprintf(fout,"\n\nHistogram of the unique-unique degree of the ends of the "
	    "true unique chunks\n");
    print_histogram(fout,ndegree_of_true_unique_to_unique_histogram, 0, 1);
#endif
#ifdef SIMINFO2
    {
      int ii;
      for(ii=0;ii<=NLETTERS;ii++){
	fprintf(fout,"\n\nHistogram of the degree of the ends of the chunks\n"
		" with repeat letter %c\n",(char)(ii+64));
	print_histogram(fout,ndegree_of_chunkends_by_letter_histogram[ii], 0, 1);
      }
    }
#endif // NEVER
#endif /*SIMINFO*/
  }
  free_histogram(nfrag_in_chunk_histogram);
  free_histogram(nfrag_essential_in_chunk_histogram);
  free_histogram(nbase_essential_in_chunk_histogram);
  free_histogram(fragment_timesinchunks_histogram);
  free_histogram(rho_histogram);
  free_histogram(coverage_histogram);
  {
    int ii;
    for(ii=0;ii<MAX_NUM_CHUNK_LABELS;ii++){
      free_histogram(labeled_unitig_histogram[ii]);
    }
  }

  free_histogram(length_of_unitigs_histogram);
#ifdef GENINFO
  free_histogram(true_uniques);
  free_histogram(true_repeats);
#endif
#ifdef NEVER
  free_histogram(ndegree_of_chunkends_histogram);
  free_histogram(ndegreesq_of_chunkends_histogram);
  free_histogram(ndegree_of_disc_unique_chunkends_histogram);
  free_histogram(ndegree_of_disc_unique_to_unique_histogram);
  free_histogram(ndegree_of_gen_unique_chunkends_histogram);
  free_histogram(ndegree_of_true_unique_to_unique_histogram);
#ifdef SIMINFO
  {
    int ii;
    for(ii=0;ii<=NLETTERS;ii++){
      free_histogram(ndegree_of_chunkends_by_letter_histogram[ii]);
    }
  }
#endif /*SIMINFO*/
#endif //NEVER
  safe_free(fragment_visited);
  safe_free(fragment_timesinchunks);
  destroy_FragmentHash(afr_to_avx);
} 

void chunk_graph_analysis
(/* Input Only */
 const int        analysis_flag,
 const IntFragment_ID max_frag_iid,
 Tfragment        frags[],     /* The internal representation of
			   the fragment reads. I have one
			   extra for the segstart field. */
 Tedge            edges[],     /* The internal representation of the
			   overlaps. */
 VA_TYPE(char)    frag_annotations[], /* Fragment annotations. */
 const BPTYPE     nbase_in_genome,
 const float      cgb_unique_cutoff,
 const float      global_fragment_arrival_rate,
 const char       bubble_boundaries_filename[],
 TChunkFrag       chunkfrags[],

 /* Modify the annotation,  by setting the "isrc" index into 
	chunksrc array. */
 TChunkMesg       thechunks[],
 /* Output only */
 VA_TYPE(char) chunksrcs[], /* The character array used to store the annotations. */
 FILE    *fcga,
 FILE    *fcam,
 FILE    *fp_unitig_statistics,
 FILE    *fwrn
 ) 
{
  const int ProcessFragmentAnnotationsForSimulatorCoordinates = (analysis_flag > 1);
  const IntFragment_ID nfrag = GetNumFragments(frags);
  //const IntEdge_ID nedge = GetNumEdges(edges);
  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);

  ChunkAnnotation *chunkinfo
    = (ChunkAnnotation *)safe_malloc((nchunks)*sizeof(ChunkAnnotation));

#ifdef DEBUG
  time_t tp1,tp2;
#endif
  Tfraginfo     *fraginfo = NULL;

  assert(frags      != NULL);
  assert(edges      != NULL);

  assert(chunkfrags != NULL);
  assert(thechunks  != NULL);
  assert(chunkinfo !=NULL);

#if 0
  assert(fcga != NULL);
  assert(fp_unitig_statistics != NULL);
  assert(fwrn != NULL);
#else
  warning_AS(fcga != NULL);
  warning_AS(fp_unitig_statistics != NULL);
  warning_AS(fwrn != NULL);
#endif

  if( ProcessFragmentAnnotationsForSimulatorCoordinates ) {
    assert(NULL == fraginfo);
    fraginfo   = CreateVA_Afraginfo(nfrag);
    assert(NULL != fraginfo);
    EnableRangeVA_Afraginfo(fraginfo,nfrag);
    setup_fraginfo( max_frag_iid, frags, frag_annotations, fraginfo);
    //  mark_invalid_edges( nfrag, fraginfo, edges);
    
    annotate_the_chunks_with_coordinate_info
      (/*Input Only*/
       max_frag_iid,
       frags,
       edges,
       fraginfo,
       chunkfrags,
       thechunks,
       /*Output Only*/
       chunkinfo);
    fprintf(stderr,"after annotate_the_chunks_with_coordinate_info\n");
  } else {
    const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
    IntChunk_ID ichunk;
    for(ichunk=0;ichunk<nchunks;ichunk++) {
      chunkinfo[ichunk].atip    = 0;
      chunkinfo[ichunk].btip    = 0;
      chunkinfo[ichunk].abforward = -1; // unknown orientation
      // chunkinfo[ichunk].gmin    = 0;
      // chunkinfo[ichunk].gmax    = 0;
      chunkinfo[ichunk].essential_type = '@';
      chunkinfo[ichunk].contained_type = '@';
    }
  }

#ifdef DEBUG
  time(&tp1); fprintf(stderr,"Begin Statistical Graph Diagnostics\n");
#endif
  
  /* BEGIN CHUNK GRAPH ANALYSIS */ {

    /******************************************************
     * Statistical Graph Diagnostics
     ******************************************************/

    fprintf(stderr,"Begin Statistical Graph Diagnostics " __FILE__ "\n");
    fprintf(stderr,"nbase_in_genome=" BPFORMAT "\n",nbase_in_genome);
    fprintf(stderr,
	    "global_fragment_arrival_rate=%f\n",
	    global_fragment_arrival_rate);

    analyze_the_fragment_overlap_graph
      (fcga,
       max_frag_iid,
       frags,
       edges);

    analyze_the_chunks
      (fcga, fp_unitig_statistics,
       max_frag_iid,
       frags,
       edges,
       chunkinfo,
       chunkfrags,
       thechunks,
       nbase_in_genome,
       cgb_unique_cutoff,
       global_fragment_arrival_rate );

    fprintf(stderr,"End Statistical Graph Diagnostics " __FILE__ "\n");
  } /* END CHUNK GRAPH ANALYSIS */

  if(NULL == fraginfo) {
    IntChunk_ID ichunk;
    for(ichunk=0;ichunk<nchunks; ichunk++) {
      { 
	size_t nchunksrc, nsource;
	char source[1024];
	sprintf(source,	"gen> @@ [0,0]\n");
	nsource = strlen(source);
	nchunksrc = GetNumVA_char(chunksrcs);
	GetVA_AChunkMesg(thechunks,ichunk)->isrc = nchunksrc;
	EnableRangeVA_char(chunksrcs,nchunksrc+nsource+1);
	/* Remember the terminal null character for a char string. */
	strcpy(GetVA_char(chunksrcs,nchunksrc),source);
      }
    }
  }
#if defined(GENINFO) || defined(SIMINFO)
  else /* (NULL != fraginfo) */ {
  /* 
     genome location type:
     '@'=no information (default)
     'u'=true unique chunk (a.k.a. valid unitig)
     'r'=true repeat chunk (a.k.a. invalid unitig)
     
     simulator annotation letter for repeat sequences:
     '@'= no information (default)
     'A'-'Z'=repeat sequence letter
     
     Proposed format for the first two line of the source text:
     "%1[@ur] [%ld,%ld]\n%1[@-Z]%1[A-Z]*\n"
     The first line gives information about the fragmentation of the
     genome. The second line gives information about the annotation 
     of the DNA sequence in the chunk.  
  */
  {
    // FILE *fout=stdout;
    IntChunk_ID ichunk;
#ifdef SIMINFO
    int false_positive_brc = 0, false_negative_brc = 0;
    int simulator_bpts=0, unitigger_bpts=0;
    int unitigger_bpt=FALSE, simulator_bpt=FALSE;
#endif

#define  NUM_COLOURS   16

    enum {
      UNUSED_COLOUR,

      /* proper unitig colors */
      UNIQUE_UNITIG_COLOUR,
      CONSISTENT_UNITIG_COLOUR,
      REPEAT_UNITIG_COLOUR,
      BADUNIQUE_UNITIG_COLOUR,
      CONT_BADUNIQUE_UNITIG_COLOUR,

      /* singleton unitig colors */
      ONE_FRAG_COLOUR,
      SPANNING_FRAG_COLOUR,
      MULTICONT_FRAG_COLOUR,
      BRANCHMULTICONT_FRAG_COLOUR,
      HANGING_CRAPPY_FRAG_COLOUR,
      SOLO_FRAG_COLOUR,
      HANGING_FRAG_COLOUR,
      ORPHAN_FRAG_COLOUR,

    /* branch point colors */
      LEFTBP_COLOUR,
      RIGHTBP_COLOUR
    } CelamyColours;
    
    char  * Colour_String [NUM_COLOURS]
      = {
	"C000000 T2 S  # Unused",

	"CFFFF00 T2 S  # UniqueUnitig",
	"C00FF00 T2 S  # ConsistentUnitig",
	"CFF0000 T2 S  # RepeatUnitig",
	"CFF00FF T2 S  # BadUniqueUnitig",
	"CFF9A11 T2 S  # ContBadUniqueUnitig",


	"C0F0F0F T2 S  # OneFrag",
        "C0000FF T2 S  # SpanningFrag",
        "C00FFFF T2 S  # MultiContFrag",  // Cyan
	"CE000A0 T2 S  # BranchMultiContFrag",
        "CFF9A9A T2 S  # HangingSpurFrag",     
        "C009A00 T2 S  # SoloFrag",
        "CCF9A00 T2 S  # HangingFrag",
        "C9A9A00 T2 S  # OrphanFrag",

	"C00FF00 # LeftBP",
	"CFF0000 # RightBP"
      };

    if(NULL != fcam) {
      int icolour;
      for(icolour=0; icolour<NUM_COLOURS; icolour++) {
        fprintf(fcam,"%d: %s\n",icolour,Colour_String[icolour]);
      }
    }

    if(bubble_boundaries_filename){
      FILE * cgb_bubbles = NULL;
      FragmentHashObject *afr_to_avx = NULL; 

      SAFE_FOPEN(cgb_bubbles, bubble_boundaries_filename,"r");
      afr_to_avx = create_FragmentHash(max_frag_iid+1);

      {
	IntFragment_ID vid;
	for(vid=0;vid<nfrag;vid++) { 
	  const IntFragment_ID iid = get_iid_fragment(frags,vid);
	  set_vid_FragmentHash(afr_to_avx,iid,vid);
	}
      }

      fprintf(fcam, "8BubbleColor: C0A0A0A T7 # PossibleBubble\n");
      if( cgb_bubbles != NULL ) {
        int bubble_count=0;
        IntFragment_ID start_iid;
        IntFragment_ID end_iid;
        int start_sx;
        int end_sx;

        clearerr( cgb_bubbles);
        while ( ! feof(cgb_bubbles) ) {
          int nitems;
	  IntFragment_ID start_vid, end_vid;

          nitems = fscanf(cgb_bubbles,
                          F_IID " %d " F_IID " %d\n", 
                          &start_iid, &start_sx, &end_iid, &end_sx);
          if(nitems != 4 ) break;
#if 0
          fprintf(stdout,
                 "CGABUBBLES: " F_IID " %d " F_IID " %d\n", 
                 start_iid, start_sx, end_iid, end_sx);
#endif
	  start_vid = get_vid_FragmentHash(afr_to_avx,start_iid);
	  end_vid   = get_vid_FragmentHash(afr_to_avx,end_iid);
	  assert(get_iid_fragment(frags,start_vid) == start_iid);
	  assert(get_iid_fragment(frags,end_vid) == end_iid);

          { 
            const BPTYPE start_pos
	      = ( !start_sx ?  // Include the whole start fragment.
		  get_genend_fraginfo(fraginfo,start_vid) :
		  get_genbgn_fraginfo(fraginfo,start_vid));
            const BPTYPE end_pos
	      = ( !end_sx ?    // Include the whole end fragment.
		  get_genend_fraginfo(fraginfo,end_vid) :
		  get_genbgn_fraginfo(fraginfo,end_vid));
            const int forward_bubble = ( start_pos <  end_pos );
            
            const BPTYPE low_pos = (  forward_bubble ? start_pos : end_pos);
            const BPTYPE hgh_pos = ( !forward_bubble ? start_pos : end_pos);
            const IntFragment_ID low_iid = (  forward_bubble ? start_iid : end_iid);
            const IntFragment_ID hgh_iid = ( !forward_bubble ? start_iid : end_iid);
            const int low_sx  = (  forward_bubble ? start_sx : end_sx);
            const int hgh_sx  = ( !forward_bubble ? start_sx : end_sx);
            char comment[100];
            
	    if (hgh_pos - low_pos <= 10000) {

	      sprintf(comment, F_IID ":%d " F_IID ":%d",
		      low_iid, low_sx, hgh_iid, hgh_sx );
	      
	      fprintf(fcam, 
		      "9Bubble%d: "
		      BPFORMAT " A8BubbleColor " BPFORMAT " R0 # %s\n",
		      bubble_count, low_pos, hgh_pos, comment);
	    }
          }
          bubble_count++;
        }
        SAFE_FCLOSE(cgb_bubbles);
      }
      destroy_FragmentHash(afr_to_avx);
    }

    fprintf(stderr,"XXXXB\n");

    for(ichunk=0;ichunk<nchunks; ichunk++) {
      const int abforward = (chunkinfo[ichunk].abforward);
      const BPTYPE gen_low_coord 
	= (abforward ? chunkinfo[ichunk].atip : chunkinfo[ichunk].btip);
      const BPTYPE gen_hgh_coord 
	= (chunkinfo[ichunk].essential_type == 'u' ?
	   (! abforward ? chunkinfo[ichunk].atip : chunkinfo[ichunk].btip) :
	   (gen_low_coord + GetVA_AChunkMesg(thechunks,ichunk)->bp_length)    );
      // Is this after consensus?  If so then this is wrong!!!

      const BPTYPE rho 
	= GetVA_AChunkMesg(thechunks,ichunk)->rho;
      const int number_of_randomly_sampled_fragments_in_chunk
	= count_the_randomly_sampled_fragments_in_a_chunk
	( frags, chunkfrags, thechunks, ichunk);
      const float coverage_statistic 
	= compute_coverage_statistic
	( rho,
	  number_of_randomly_sampled_fragments_in_chunk,
	  global_fragment_arrival_rate );

      { 
	size_t nchunksrc, nsource;
	char source[1024];
	sprintf(source,
		"gen> %c%c [" BPFORMAT "," BPFORMAT "]\n",
		chunkinfo[ichunk].essential_type,
		chunkinfo[ichunk].contained_type,
		(abforward ? gen_low_coord : gen_hgh_coord), 
		(abforward ? gen_hgh_coord : gen_low_coord));
	nsource = strlen(source);
	nchunksrc = GetNumVA_char(chunksrcs);
	GetVA_AChunkMesg(thechunks,ichunk)->isrc = nchunksrc;
	EnableRangeVA_char(chunksrcs,nchunksrc+nsource+1);
	/* Remember the terminal null character for a char string. */
	strcpy(GetVA_char(chunksrcs,nchunksrc),source);
      }

#ifdef SIMINFO
      /* Begin checking the new chunk based branch points. */ {
 
	/* The simulator^s definition of a branch point is as
           follows. It is the first transition from a repeat region to
           another repeat region or the unique region of the genome is
           a branch point.  */

	/* A false positive branch point is defined as a branch point
           more than 10 bps into the unique region of the genome.
           This can be checked by taking the branch points of the
           unitig and finding the fragments that include the branch
           points. If the branch points are in the unique region of
           the genome as determined by the fragments simulator
           information, then it is a false positive branch point.  */
	

	/* A false negative branch point is defined when a simulator
	 branch point is detected in a unitig, but a branch point was
	 not found.  */
	
	const IntFragment_ID chunk_avx
          = GetVA_AChunkMesg(thechunks,ichunk)->chunk_avx;
	const IntFragment_ID chunk_bvx
          = GetVA_AChunkMesg(thechunks,ichunk)->chunk_bvx;

	const int chunk_asx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_asx;
	const int chunk_bsx = GetVA_AChunkMesg(thechunks,ichunk)->chunk_bsx;
	int gen_frag_bplength;
	IntFragment_ID ifrag
	int isuffix, branching, position;
	
	branching = GetVA_AChunkMesg(thechunks,ichunk)->a_branch_type;
	position  = GetVA_AChunkMesg(thechunks,ichunk)->a_branch_point;
	ifrag   = chunk_avx;
	isuffix = chunk_asx;
	gen_frag_bplength = get_genend_fraginfo(fraginfo,ifrag) 
	  - get_genbgn_fraginfo(fraginfo,ifrag);
	gen_frag_bplength = (gen_frag_bplength > 0 ? 
			     gen_frag_bplength : -gen_frag_bplength); 

	assert(branching != AS_INTO_REPEAT);
	unitigger_bpt = (branching == AS_INTO_UNIQUE);
	simulator_bpt = 
	  ((isuffix == TRUE)&&
	   // Is there enough repeated sequence?
	   (get_pre_end_fraginfo(fraginfo,ifrag) -
	    get_pre_brp_fraginfo(fraginfo,ifrag) >= BPT_MIN_PREFIX)
	   // And not contained in the repeat region.
	   &&(get_pre_brp_fraginfo(fraginfo,ifrag) != 0)
#ifdef NEVER
	   &&
	   (get_pre_end_fraginfo(fraginfo,ifrag) -
	    get_pre_brp_fraginfo(fraginfo,ifrag) <= 
	    get_length_fragment(frags,ifrag) - BPT_MIN_SUFFIX) 
#endif
	   )
	  ||
	  ((isuffix == FALSE)&&
	   (get_suf_brp_fraginfo(fraginfo,ifrag) -
	    get_suf_end_fraginfo(fraginfo,ifrag) >= BPT_MIN_PREFIX)
	   // And not contained in the repeat region.
	   &&(get_suf_brp_fraginfo(fraginfo,ifrag) != gen_frag_bplength)
#ifdef NEVER
	   &&(get_suf_brp_fraginfo(fraginfo,ifrag) -
	      get_suf_end_fraginfo(fraginfo,ifrag) <= 
	      get_length_fragment(frags,ifrag) - BPT_MIN_SUFFIX)
#endif
	   );
	
	// Screen against a repeat region wholly contained by the fragment.
	simulator_bpt = simulator_bpt &&
	  ((get_pre_brp_fraginfo(fraginfo,ifrag) >
	    get_suf_brp_fraginfo(fraginfo,ifrag)) &&
	   (get_pre_brp_fraginfo(fraginfo,ifrag) != 0) &&
	   (get_suf_brp_fraginfo(fraginfo,ifrag) != gen_frag_bplength));

	if(unitigger_bpt) {unitigger_bpts++;} 
	if(simulator_bpt) {simulator_bpts++;} 

	if(unitigger_bpt && !simulator_bpt
	   // Associated with an ALU
	   // &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'B') ||
	   //      (get_suf_let_fraginfo(fraginfo,ifrag) == 'B')))
	   // Associated with a short tandem repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'E') ||
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'E')))
	   // Associated with a medium tandem repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'G') ||
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'G')))
#ifdef NEVER
	    // Inside a ALU repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'B') &&
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'B')))
	   // Inside a short tandem repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'E') &&
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'E')))
	   // Inside a medium tandem repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'G') &&
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'G')))
#endif
	   ) {
	  false_positive_brc ++;
#ifdef DEBUG16
	  fprintf(fout,"The simulator missed a BRC\n");
#endif
	}
	if(!unitigger_bpt && simulator_bpt) {
	  false_negative_brc ++;
#ifdef DEBUG17
	  fprintf(fout,"FN BPT: " F_IID ",A,%c,%d " F_IID ",%d, %c,%d,%d:%c,%d,%d\n",
		  ichunk,branching,position,
		  get_iid_fragment(frags,ifrag),
		  isuffix,
		  get_suf_let_fraginfo(fraginfo,ifrag),
		  get_suf_brp_fraginfo(fraginfo,ifrag),
		  get_suf_end_fraginfo(fraginfo,ifrag),
		  get_pre_let_fraginfo(fraginfo,ifrag),
		  get_pre_brp_fraginfo(fraginfo,ifrag),
		  get_pre_end_fraginfo(fraginfo,ifrag)
		  );
#endif
	}

	branching = GetVA_AChunkMesg(thechunks,ichunk)->b_branch_type;
	position  = GetVA_AChunkMesg(thechunks,ichunk)->b_branch_point;
	ifrag   = chunk_bvx;
	isuffix = chunk_bsx;
	gen_frag_bplength = get_genend_fraginfo(fraginfo,ifrag) 
	  - get_genbgn_fraginfo(fraginfo,ifrag);
	gen_frag_bplength = (gen_frag_bplength > 0 ?
			     gen_frag_bplength : -gen_frag_bplength); 

	assert(branching != AS_INTO_REPEAT);
	unitigger_bpt = (branching == AS_INTO_UNIQUE);
	simulator_bpt = 
	  ((isuffix == TRUE)&&
	   (get_pre_end_fraginfo(fraginfo,ifrag) -
	    get_pre_brp_fraginfo(fraginfo,ifrag) >= BPT_MIN_PREFIX)
	   // And not contained in the repeat region.
	   &&(get_pre_brp_fraginfo(fraginfo,ifrag) != 0)  
#ifdef NEVER
	   &&
	   (get_pre_end_fraginfo(fraginfo,ifrag) -
	    get_pre_brp_fraginfo(fraginfo,ifrag) <= 
	    get_length_fragment(frags,ifrag) - BPT_MIN_SUFFIX) 
#endif
	   )
	  ||
	  ((isuffix == FALSE)&&
	   (get_suf_brp_fraginfo(fraginfo,ifrag) -
	    get_suf_end_fraginfo(fraginfo,ifrag) >= BPT_MIN_PREFIX)
	   // And not contained in the repeat region.
	   &&(get_suf_brp_fraginfo(fraginfo,ifrag) != gen_frag_bplength)
#ifdef NEVER
	   &&(get_suf_brp_fraginfo(fraginfo,ifrag) -
	      get_suf_end_fraginfo(fraginfo,ifrag) <= 
	      get_length_fragment(frags,ifrag) - BPT_MIN_SUFFIX)
#endif
	   );

	// Screen against a repeat region wholly contained by the fragment.
	simulator_bpt = simulator_bpt &&
	  ((get_pre_brp_fraginfo(fraginfo,ifrag) >
	    get_suf_brp_fraginfo(fraginfo,ifrag)) &&
	   (get_pre_brp_fraginfo(fraginfo,ifrag) != 0) &&
	   (get_suf_brp_fraginfo(fraginfo,ifrag) != gen_frag_bplength));

	if(unitigger_bpt) {unitigger_bpts++;} 
	if(simulator_bpt) {simulator_bpts++;} 

	if(unitigger_bpt && !simulator_bpt
	   // Associated with an ALU
	   //&&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'B') ||
	   //    (get_suf_let_fraginfo(fraginfo,ifrag) == 'B')))
	   // Associated with a short tandem repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'E') ||
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'E')))
	   // Associated with a medium tandem repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'G') ||
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'G')))
#ifdef NEVER
	    // Inside a ALU repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'B') &&
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'B')))
	   // Inside a short tandem repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'E') &&
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'E')))
	   // Inside a medium tandem repeat
	   &&(!((get_pre_let_fraginfo(fraginfo,ifrag) == 'G') &&
		(get_suf_let_fraginfo(fraginfo,ifrag) == 'G')))
#endif
	   ) {
	  false_positive_brc ++;
#ifdef DEBUG16
	  fprintf(fout,"The simulator missed a BRC\n");
#endif
	}
	if(!unitigger_bpt && simulator_bpt) {
	  false_negative_brc ++;
#ifdef DEBUG17
	  fprintf(fout,"FN BPT: " F_IID ",B,%c,%d " F_IID ",%d, %c,%d,%d:%c,%d,%d\n",
		  ichunk,branching,position,
		  get_iid_fragment(frags,ifrag),
		  isuffix,
		  get_suf_let_fraginfo(fraginfo,ifrag),
		  get_suf_brp_fraginfo(fraginfo,ifrag),
		  get_suf_end_fraginfo(fraginfo,ifrag),
		  get_pre_let_fraginfo(fraginfo,ifrag),
		  get_pre_brp_fraginfo(fraginfo,ifrag),
		  get_pre_end_fraginfo(fraginfo,ifrag)
		  );
#endif
	}
      } /* End checking the new branch points. */
#endif // SIMINFO

      {
	int ia,im_left,im_right;

	int valid=TRUE, evalid=TRUE, cvalid=TRUE,
	  unique=FALSE, singleton=TRUE;
	BPTYPE pos_left_end,pos_right_end,pos_left_bpt,pos_right_bpt;
	char chunk_name[1024] = {0};
	char source[1024] = {0};

	// assert(ichunk >= 0);
	evalid = ( chunkinfo[ichunk].essential_type == 'u' );
	cvalid = ( chunkinfo[ichunk].contained_type == 'u' );
	valid  = evalid && cvalid;
	unique = ( coverage_statistic >= cgb_unique_cutoff);
	singleton = (GetVA_AChunkMesg(thechunks,ichunk)->num_frags == 1);
	  
	/* Set the color of each unitig for Gene^s genome visualizer */
	ia = -1; // The default colour.

	// override for discriminator unique chunks
	ia = ( ( valid) && ( unique) ? UNIQUE_UNITIG_COLOUR : ia);
	ia = ( ( valid) && (!unique) ? CONSISTENT_UNITIG_COLOUR : ia);

	// override for non-discriminator unique chunks
	ia = ( (!valid) && (!unique)
               ? REPEAT_UNITIG_COLOUR : ia);
	ia = ( (!valid) && ( unique) && ( evalid)
               ? CONT_BADUNIQUE_UNITIG_COLOUR : ia);
	ia = ( (!valid) && ( unique) && (!evalid)
               ? BADUNIQUE_UNITIG_COLOUR : ia);

	assert(ia != -1);

	// override for singletons, they were CONSISTENT_UNIQUE_COLOUR.
	if( singleton ) {
	  ia = ( singleton ? ONE_FRAG_COLOUR : ia);

	  // override for various singleton unitigs,
	  // assuming the special frag is first in the list.
	  { 
	    IntFragment_ID ivc = GetVA_AChunkMesg(thechunks,ichunk)->f_list;
	    IntFragment_ID vid = GetVA_AChunkFrag(chunkfrags,ivc)->vid;
	    ia = (get_lab_fragment(frags,vid) == AS_CGB_SOLO_FRAG
		  ? SOLO_FRAG_COLOUR : ia);
	    ia = (get_lab_fragment(frags,vid) == AS_CGB_HANGING_FRAG
		  ? HANGING_FRAG_COLOUR : ia);
	    ia = (get_lab_fragment(frags,vid) == AS_CGB_THRU_FRAG
		  ? SPANNING_FRAG_COLOUR : ia);
	    ia = (get_lab_fragment(frags,vid) == AS_CGB_ORPHANEDCONT_FRAG 
		  ? ORPHAN_FRAG_COLOUR : ia);
	    ia = (get_lab_fragment(frags,vid) == AS_CGB_MULTICONT_FRAG
		  ? MULTICONT_FRAG_COLOUR : ia);
	    ia = (get_lab_fragment(frags,vid) == AS_CGB_BRANCHMULTICONT_FRAG
		  ? BRANCHMULTICONT_FRAG_COLOUR : ia);
	    ia = (get_lab_fragment(frags,vid) == AS_CGB_HANGING_CRAPPY_FRAG
		  ? HANGING_CRAPPY_FRAG_COLOUR : ia);
	  }
	}

	sprintf(chunk_name,F_IID "%c" F_IID "%c",
		get_iid_fragment
		(frags,GetVA_AChunkMesg(thechunks,ichunk)->chunk_avx),
		(GetVA_AChunkMesg(thechunks,ichunk)->chunk_asx ? 'T' : 'F' ),
		get_iid_fragment
		(frags,GetVA_AChunkMesg(thechunks,ichunk)->chunk_bvx),
		(GetVA_AChunkMesg(thechunks,ichunk)->chunk_bsx ? 'T' : 'F' ));
		
	sprintf(source,"unitig "
		F_IID ",%c:"
		F_IID ":%d " F_IID ":%d",
		(ichunk), (abforward ? 'f' : 'r'),
		get_iid_fragment
		(frags,(abforward ? 
			GetVA_AChunkMesg(thechunks,ichunk)->chunk_avx:
			GetVA_AChunkMesg(thechunks,ichunk)->chunk_bvx)),
		 (abforward ? 
		 GetVA_AChunkMesg(thechunks,ichunk)->chunk_asx:
		 GetVA_AChunkMesg(thechunks,ichunk)->chunk_bsx),
		get_iid_fragment
		(frags,(!abforward ? 
		 GetVA_AChunkMesg(thechunks,ichunk)->chunk_avx:
		 GetVA_AChunkMesg(thechunks,ichunk)->chunk_bvx)),
		(!abforward ? 
		 GetVA_AChunkMesg(thechunks,ichunk)->chunk_asx:
		 GetVA_AChunkMesg(thechunks,ichunk)->chunk_bsx)
		);

	{
	  BPTYPE bp_length = GetVA_AChunkMesg(thechunks,ichunk)->bp_length;
	  // If this is after consensus, then this is wrong.

	  BPTYPE a_branch_point 
	    = GetVA_AChunkMesg(thechunks,ichunk)->a_branch_point;
	  BPTYPE b_branch_point 
	    = GetVA_AChunkMesg(thechunks,ichunk)->b_branch_point;
	  assert(bp_length > 0);
	  assert(a_branch_point >= 0);
	  assert(b_branch_point >= 0);
	  a_branch_point = min(a_branch_point,bp_length);
	  b_branch_point = min(b_branch_point,bp_length);
	  pos_left_end  = gen_low_coord;
	  pos_right_end = pos_left_end + bp_length;

#if 0
          {
            // Plot unitigs with sum of overhangs length instead of maximal length.
            BPTYPE slop_a = bp_length - (GetVA_AChunkMesg(thechunks,ichunk)->rho);
            BPTYPE slop_b = max(slop_a,0);
            BPTYPE slop_l = slop_b / 2;
            BPTYPE slop_r = slop_b - slop_l;
            
            pos_left_end  += slop_l;
            pos_right_end -= slop_r;
          }
#endif
          
	  pos_left_bpt = pos_left_end 
	    + (abforward ? a_branch_point : b_branch_point);
	  pos_right_bpt = pos_right_end
	    - (abforward ? b_branch_point : a_branch_point);
	}
	if( pos_left_bpt < pos_right_bpt ) {
	  // The BLACK branch point is associated with the left chunk-end.
	  im_left = LEFTBP_COLOUR;
	  im_right = RIGHTBP_COLOUR;
	} else { 
	  BPTYPE tmp;
	  tmp = pos_left_bpt; pos_left_bpt = pos_right_bpt; pos_right_bpt = tmp;
	  // The BLACK branch point is associated with the left chunk-end.
	  im_left = RIGHTBP_COLOUR;
	  im_right = LEFTBP_COLOUR;
	}
	if( pos_left_end < pos_left_bpt && 
	    pos_right_bpt < pos_right_end ) {
	  fprintf(fcam,
		  "%s: " BPFORMAT " A%d M" BPFORMAT " A%d M" BPFORMAT " A%d " BPFORMAT " #%s\n",
		  chunk_name,
		  // place the chunk in one of its valid locations.
		  pos_left_end,
		  im_left,
		  pos_left_bpt,	
		  im_right,
		  pos_right_bpt,
		  ia,
		  pos_right_end,
		  source
		  );
	} else if( pos_left_end >= pos_left_bpt && 
		   pos_right_bpt < pos_right_end ) {
	  fprintf(fcam,
		  "%s: " BPFORMAT " A%d M" BPFORMAT " A%d " BPFORMAT " #%s\n",
		  chunk_name,
		  // place the chunk in one of its valid locations.
		  pos_left_end,
		  im_right,
		  pos_right_bpt,
		  ia,
		  pos_right_end,
		  source
		  );
	} else if( pos_left_end < pos_left_bpt && 
		   pos_right_bpt >= pos_right_end ) {
	  fprintf(fcam,
		  "%s: " BPFORMAT " A%d M" BPFORMAT " A%d " BPFORMAT " #%s\n",
		  chunk_name,
		  // place the chunk in one of its valid locations.
		  pos_left_end,
		  im_left,
		  pos_left_bpt,
		  ia,
		  pos_right_end,
		  source
		  );
	} else if( pos_left_end >= pos_left_bpt && 
		   pos_right_bpt >= pos_right_end ) {
	  fprintf(fcam,
		  "%s: " BPFORMAT " A%d " BPFORMAT " #%s\n",
		  chunk_name,
		  // place the chunk in one of its valid locations.
		  pos_left_end,
		  ia,
		  pos_right_end,
		  source
		  );
	} else { assert(FALSE);}
      }

    }
#ifdef DEBUGGING
    fprintf(fout,"Chunk end branch points\n"
	    "simulator bpts=%d, unitigger bpts=%d,"
	    "false_positive_brc=%d,false_negative_brc=%d\n",
	    simulator_bpts, unitigger_bpts,
	    false_positive_brc,false_negative_brc);
#endif
  }
    fprintf(stderr,"XXXXC\n");

  }
#endif /* defined(GENINFO) || defined(SIMINFO) */

#ifdef DEBUG
  time(&tp2); fprintf(stderr,"%10" F_TIME_TP " sec: Finished Statistical Graph Diagnostics\n",
		      (tp2-tp1));
#endif

#ifdef DEBUG_VISUAL
  /******************************************************
   * Visual Graph Diagnostics
   ******************************************************/

#ifdef DEBUG
  if(TIMINGS) {
    time(&tp1); fprintf(stderr,"Begin visual diagnostics\n");
  }
#endif
  if(NULL != fraginfo) {
    char File_Prefix[] = "TMP";
    char strtmp[1024]={0};
    
    strcpy(strtmp,File_Prefix);
    strcat(strtmp,"graph0f");
    graph_diagnostics(strtmp,
		      frags, edges, 
		      AS_CGB_SINGLECONT_FRAG,
		      AS_CGB_MARKED_BY_BRANCH_DVT,
		      FALSE,fraginfo);

    strcpy(strtmp,File_Prefix);
    strcat(strtmp,"graph1f");
    graph_diagnostics(strtmp,
		      frags, edges, 
		      AS_CGB_MULTICONT_FRAG,
		      AS_CGB_TO_CONTAINED_EDGE,
		      FALSE,fraginfo);

    strcpy(strtmp,File_Prefix);
    strcat(strtmp,"graph2f");
    graph_diagnostics(strtmp,
		      frags, edges, 
		      AS_CGB_MULTICONT_FRAG,
		      AS_CGB_TO_CONTAINED_EDGE,
		      TRUE,fraginfo);

    strcpy(strtmp,File_Prefix);
    strcat(strtmp,"graph3f");
    graph_diagnostics(strtmp, frags, edges, 
		      AS_CGB_INTERCHUNK_FRAG,
		      AS_CGB_INTRACHUNK_EDGE,
		      TRUE,fraginfo);
  }

#ifdef DEBUG
  if(TIMINGS) {
    time(&tp2); fprintf(stderr,"%10" F_TIME_TP " sec: Finished visual diagnostics\n",
			(tp2-tp1));
  }
#endif
#endif /*DEBUG_VISUAL*/

  safe_free(chunkinfo);
  if( ProcessFragmentAnnotationsForSimulatorCoordinates ) {
    DeleteVA_Afraginfo(fraginfo);
  }
}
