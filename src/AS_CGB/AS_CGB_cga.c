
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

static char CM_ID[] = "$Id: AS_CGB_cga.c,v 1.17 2007-07-20 04:47:37 brianwalenz Exp $";

//  A chunk graph analyzer. This functional unit computes graph
//  statistics, and writes the chunk graph in the term representation
//  which is suitable for "DaVinci".

#include "AS_CGB_all.h"
#include "AS_CGB_histo.h"

typedef enum {
  CGA_SINGLETON_CONTAINED,    // multiply contained and orphaned contained fragments that form singleton unitigs.
  CGA_SINGLETON_NONCONTAINED, // The rest of the singleton unitigs.
  CGA_NONSINGLETON_SPANNED,   // The non-singleton unitigs that are spanned by one fragment.
  CGA_NONSINGLETON_NONSPANNED
} ChunkLabel;

#define MAX_NUM_CHUNK_LABELS 4

char * ChunkLabelDesc[MAX_NUM_CHUNK_LABELS] = {
  "singleton contained ",      // Contained fragments that form singleton unitigs.
  "singleton non-contained",   // The rest of the singleton unitigs.
  "non-singleton spanned",     // The non-singleton unitigs that are spanned by one fragment.
  "non-singleton non-spanned"  // The rest of the non-singleton unitigs.
};


static void analyze_the_fragment_overlap_graph(FILE *fout,
                                               const IntFragment_ID max_frag_iid,
                                               Tfragment frags[],
                                               Tedge edges[]) {

  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

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




static void analyze_the_chunks(FILE *fout,
                               FILE *fp_unitig_statistics,
                               const IntFragment_ID max_frag_iid,
                               Tfragment frags[],
                               Tedge edges[],
                               TChunkFrag chunkfrags[],
                               TChunkMesg thechunks[],
                               const int64  nbase_in_genome,
                               const int recalibrate_global_arrival_rate,
                               const float cgb_unique_cutoff,
                               const float global_fragment_arrival_rate) {

  IntChunk_ID ichunk;
  IntFragment_ID num_of_chunks[MAX_NUM_CHUNK_LABELS]={0};
  IntFragment_ID nfrag_in_all_chunks=0;
  IntFragment_ID nfrag_essential_in_all_chunks=0;
  IntFragment_ID nfrag_contained_in_all_chunks=0;
  int64  nbase_essential_in_all_chunks=0;
  int64  rho_in_all_chunks = 0;
  IntFragment_ID n_rs_frag_in_all_chunks = 0;
  IntFragment_ID n_nr_frag_in_all_chunks = 0;

  const int nsample=500;
  const int nbucket=500;
  MyHistoDataType zork;

  IntFragment_ID *fragment_visited = NULL;
  int            *fragment_timesinchunks = NULL;
  IntFragment_ID *afr_to_avx = NULL;
  IntFragment_ID  nfrag   = GetNumFragments(frags);
  IntChunk_ID     nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  Histogram_t    *length_of_unitigs_histogram = create_histogram(nsample,nbucket,0,TRUE);

  Histogram_t * rho_histogram = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * coverage_histogram = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * nfrag_in_chunk_histogram = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * nfrag_essential_in_chunk_histogram = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * nbase_essential_in_chunk_histogram = create_histogram(nsample,nbucket,0,TRUE);
  Histogram_t * fragment_timesinchunks_histogram = create_histogram(nsample,nbucket,0,TRUE);

  Histogram_t * labeled_unitig_histogram[MAX_NUM_CHUNK_LABELS];

  {
    int ii;
    for(ii=0;ii<MAX_NUM_CHUNK_LABELS;ii++){
      labeled_unitig_histogram[ii] = create_histogram(nsample,nbucket,0,TRUE);
      extend_histogram(labeled_unitig_histogram[ii],sizeof(MyHistoDataType), myindexdata,mysetdata,myaggregate,myprintdata);
    }
  }

  extend_histogram(length_of_unitigs_histogram, sizeof(MyHistoDataType),
		   myindexdata,mysetdata,myaggregate,myprintdata);

  fragment_visited       = (IntFragment_ID *)safe_calloc(nfrag,     sizeof(IntFragment_ID));
  fragment_timesinchunks = (int            *)safe_calloc(nfrag,     sizeof(int));

  /* Initialize a flag for chunk following. */
  {
    IntFragment_ID ifrag;
    for(ifrag=0;ifrag<nfrag;ifrag++) { 
      fragment_visited[ifrag]          = FRAGMENT_NOT_VISITED;
      fragment_timesinchunks[ifrag]    = 0;
    }
  }

  // Re-hash the fragment IID to fragment VID mapping using the
  // fragments in the store.  (duplicated, search for BUILD_AFR_TO_AVX
  {
    IntFragment_ID  iv = 0;
    IntFragment_ID  max_frag_iid = 0;
    IntFragment_ID  nfrag        = GetNumFragments(frags);

    for (iv=0; iv<nfrag; iv++) {
      IntFragment_ID iid = get_iid_fragment(frags,iv);
      max_frag_iid = MAX(max_frag_iid, iid);
    }

    assert(max_frag_iid < AS_CGB_NOT_SEEN_YET);

    afr_to_avx = safe_calloc(max_frag_iid + 1, sizeof(IntFragment_ID));

    for(iv=0; iv<nfrag; iv++)
      afr_to_avx[get_iid_fragment(frags,iv)] = iv;
  }


  
  assert((!0) == 1); /* Needed for the following bitwise XOR operator. */
  for(ichunk=0;ichunk<nchunks;ichunk++) {

    const IntFragment_ID irec_start_of_chunk = GetVA_AChunkMesg(thechunks,ichunk)->f_list;
    const int64  rho = GetVA_AChunkMesg(thechunks,ichunk)->rho;
    const IntFragment_ID nfrag_in_chunk = GetVA_AChunkMesg(thechunks,ichunk)->num_frags;
    const int64  nbase_essential_in_chunk = GetVA_AChunkMesg(thechunks,ichunk)->bp_length;

    const int number_of_randomly_sampled_fragments_in_chunk
      = count_the_randomly_sampled_fragments_in_a_chunk ( frags, chunkfrags, thechunks, ichunk);
    const float coverage_statistic = compute_coverage_statistic ( rho,
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
    int64  nbase_sampled_in_chunk=0;
    int64  nbase_essential_sampled_in_chunk=0;
    int64  nbase_contained_sampled_in_chunk=0;

    ChunkLabel chunk_label;

#ifdef DEBUG07
    fprintf(stderr,
	    "Process ichunk,nchunks,nfrag_in_chunk,nbase_essential_in_chunk=\n"
            F_IID "," F_IID "," F_IID "," F_S64 "\n",
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
        const IntFragment_ID ibvx = afr_to_avx[iid];

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

                //  '@' -- essentia_type unknown -- no simulator
                //  '@' -- contained_type unknown -- no simulator

                fprintf(fp_unitig_statistics, F_IID" %d %d "F_IID" "F_S64" %d %d %d @ @ %d %d %d %d %d %d %d %d %d %d\n", 
                        ichunk,            // An identity field
                        (chunk_contained), // A selection field
                        (chunk_spanned),
                        nfrag_in_chunk, // A data field 
                        nbase_essential_in_chunk,
                        (int)rho,
                        arrival_distance,
                        coverage_index, // A selection field
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
    fprintf(fout,"%15"F_IIDP" : Total number of unitigs\n",nchunks);
    {
      int ii;
      for(ii=0;ii<MAX_NUM_CHUNK_LABELS;ii++){
	switch(ii) {
          case CGA_SINGLETON_CONTAINED:
            fprintf(fout,"%15"F_IIDP " : Total number of singleton, contained unitigs\n",
                    num_of_chunks[ii]);
            break;
          case CGA_SINGLETON_NONCONTAINED:
            fprintf(fout,"%15"F_IIDP" : Total number of singleton, non-contained unitigs\n",
                    num_of_chunks[ii]);
            break;
          case CGA_NONSINGLETON_SPANNED:
            fprintf(fout,"%15"F_IIDP" : Total number of non-singleton, spanned unitigs\n",
                    num_of_chunks[ii]);
            break;
          case CGA_NONSINGLETON_NONSPANNED:
            fprintf(fout,"%15"F_IIDP" : Total number of non-singleton, non-spanned unitigs\n",
                    num_of_chunks[ii]);
            break;
          default:
            break;
	}
      }
    }

    fprintf(fout,"%15"F_IIDP" : Total number of fragments\n", nfrag);
    fprintf(fout,"%15"F_IIDP" : Total number of fragments in all unitigs\n", nfrag_in_all_chunks);
    fprintf(fout,"%15"F_IIDP" : Total number of essential fragments in all unitigs\n", nfrag_essential_in_all_chunks);
    fprintf(fout,"%15"F_IIDP" : Total number of contained fragments in all unitigs\n", nfrag_contained_in_all_chunks);
    fprintf(fout,"%15.10f"" : Randomly sampled fragment arrival rate per bp\n", global_fragment_arrival_rate);
    fprintf(fout,"%15"F_S64P" : The sum of overhangs in all the unitigs\n", rho_in_all_chunks);
    fprintf(fout,"%15"F_S64P" : Total number of bases in all unitigs\n", nbase_essential_in_all_chunks);
    fprintf(fout,"%15"F_S64P" : Estimated number of base pairs in the genome.\n", nbase_in_genome);

    {
      IntFragment_ID nfound = 0;
      IntFragment_ID ifrag;
      for(ifrag=0;ifrag<nfrag;ifrag++) /* beginning loop through fragments */ {
	Tlab lab = get_lab_fragment(frags,ifrag);
	if( (fragment_visited[ifrag] == FRAGMENT_NOT_VISITED) &&
	    (AS_CGB_DELETED_FRAG != lab) ) {
	  fprintf(stderr,"Not visited: fragment iid=" F_IID ",ifrag=" F_IID ",ilab=%d\n",
		  get_iid_fragment(frags,ifrag),ifrag,lab);
	  nfound++;
	}
      }
      fprintf(fout,"%15" F_IIDP " : Total number of contained fragments not connected\n"
	      "                  by containment edges to essential fragments.\n",
	      nfound);
    }

    {
      compute_the_global_fragment_arrival_rate
        ( recalibrate_global_arrival_rate, cgb_unique_cutoff, fout, nbase_in_genome, frags, edges,
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
	if(!(!(ilab==AS_CGB_INTERCHUNK_FRAG)
	     ||(fragment_timesinchunks[ifrag]==1))) {
	  fprintf(stderr,"** iid=" F_IID " ifrag=" F_IID " ilab=%d fragment_timesinchunks=%d\n",
		  get_iid_fragment(frags,ifrag),
		  ifrag, ilab,
		  fragment_timesinchunks[ifrag]
		  );
	}
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
  safe_free(fragment_visited);
  safe_free(fragment_timesinchunks);
  safe_free(afr_to_avx);
} 



void
chunk_graph_analysis(THeapGlobals *heapva,
                     TStateGlobals *gstate,
                     UnitiggerGlobals *rg) {
  char strtmp2[FILENAME_MAX];

  sprintf(strtmp2,"%s.cga.0",rg->Output_Graph_Store_Prefix);
  FILE *fcga = fopen(strtmp2, "w");

  sprintf(strtmp2,"%s.cam.0",rg->Output_Graph_Store_Prefix);
  FILE *fcam = fopen(strtmp2, "w");

  sprintf(strtmp2,"%s.cus.0",rg->Output_Graph_Store_Prefix);
  FILE *fcus = fopen(strtmp2, "w");

  analyze_the_fragment_overlap_graph(fcga,
                                     gstate->max_frag_iid,
                                     heapva->frags,
                                     heapva->edges);

  analyze_the_chunks(fcga,
                     fcus,
                     gstate->max_frag_iid,
                     heapva->frags,
                     heapva->edges,
                     heapva->chunkfrags,
                     heapva->thechunks,
                     gstate->nbase_in_genome,
                     rg->recalibrate_global_arrival_rate,
                     rg->cgb_unique_cutoff,
                     gstate->global_fragment_arrival_rate);

  fclose(fcga);
  fclose(fcam);
  fclose(fcus);
}
