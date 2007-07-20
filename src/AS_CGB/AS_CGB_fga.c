
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

static char CM_ID[] = "$Id: AS_CGB_fga.c,v 1.13 2007-07-20 08:41:43 brianwalenz Exp $";

// Module: AS_CGB_fga.c
// 
// A fragment graph analyser. This functional unit computes whether
// fragment overlaps are valid.
//
// Author: Clark M. Mobarry

#include "AS_CGB_all.h"
#include "AS_CGB_histo.h"

#define OUTPUT_DERIVED_QUANTITIES
#undef HIDE_OVERHANGS
#define DOVETAIL_OR_CONTAINED_FLAG


void view_fgb_chkpnt(char * Store_File_Prefix,
                     Tfragment frags[], 
                     Tedge edges[]) {

  /* We are creating two line-by-line sortable files that can be UNIX
     diff^ed to find the differences in two assemblies.  One of the
     files has the fragments and the other is the edges. */
  
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge     = GetNumEdges(edges);
  char thePath[FILENAME_MAX]={0};

  FILE *foutv = NULL;
  FILE *foute = NULL;

  if(NULL == Store_File_Prefix) { return;}
     
  sprintf(thePath,"%s%s",Store_File_Prefix,".fgv");
  foutv = fopen(thePath,"w");

  sprintf(thePath,"%s%s",Store_File_Prefix,".fge");
  foute = fopen(thePath,"w");

  { IntFragment_ID iv0;
    for(iv0=0; iv0<nfrag; iv0++) {
      fprintf(foutv,
              //"%6d: "
              "iid %9" F_IIDP " "
              "%c%c%c "
              "raw: %2d %2d %2d %2d "
              //"src %10ld "
              "prefix:"
              //"%8ld"
              "%5" F_S32P " %5" F_S32P " "
              "suffix:"
              //"%8ld"
              "%5" F_S32P " %5" F_S32P " lab:%d\n",
              //iv0,
              get_iid_fragment(frags,iv0),
              (char)get_typ_fragment(frags,iv0),
              (get_con_fragment(frags,iv0) ? 'C' : ' '),
              (get_del_fragment(frags,iv0) ? 'D' : ' '),

              get_raw_dvt_count_vertex(frags,iv0,FALSE),
              get_raw_dvt_count_vertex(frags,iv0,TRUE),
              get_raw_frc_count_fragment(frags,iv0),
              get_raw_toc_count_fragment(frags,iv0),

              //get_segstart_vertex(frags,iv0,FALSE),
              get_seglen_dvt_vertex(frags,iv0,FALSE),
              get_seglen_frc_vertex(frags,iv0,FALSE),
              //get_segstart_vertex(frags,iv0,TRUE),
              get_seglen_dvt_vertex(frags,iv0,TRUE),
              get_seglen_frc_vertex(frags,iv0,TRUE),
              get_lab_fragment(frags,iv0)
              );
    }}
  {
    IntEdge_ID ie0;
    for(ie0=0; ie0 < nedge; ie0++) {
      const IntFragment_ID avx = get_avx_edge(edges,ie0);
      const IntFragment_ID bvx = get_bvx_edge(edges,ie0);
      const IntFragment_ID afr = get_iid_fragment(frags,avx);
      const IntFragment_ID bfr = get_iid_fragment(frags,bvx);
      const Tnes nes = get_nes_edge(edges,ie0);
      if(
	 TRUE
	 ) {
	const int asx = get_asx_edge(edges,ie0);
	const int bsx = get_bsx_edge(edges,ie0);
	const int ahg = get_ahg_edge(edges,ie0);
	const int bhg = get_bhg_edge(edges,ie0);
	const int invalid = get_inv_edge(edges,ie0);
	fprintf(foute,
		"%9" F_IIDP " "
                "%2d "
#ifndef HIDE_OVERHANGS
                "%6d "
#endif                
		"%9" F_IIDP " "
                "%2d "
#ifndef HIDE_OVERHANGS
                "%6d "
#endif                
#ifdef DOVETAIL_OR_CONTAINED_FLAG
                "%c %c %c %c"
#endif                
#ifdef OUTPUT_DERIVED_QUANTITIES
		" : %5d "
                "%c %c "
#endif
                "\n"
		,afr
                ,asx
#ifndef HIDE_OVERHANGS
                ,ahg
#endif                
		,bfr
                ,bsx
#ifndef HIDE_OVERHANGS
                ,bhg
#endif                
#ifdef DOVETAIL_OR_CONTAINED_FLAG
                ,( is_a_dvt_edge(edges, ie0) ? ( get_thickest_dvt_edge(edges, ie0) ? 'D' : 'd' ) : 'C')
                ,( is_a_dvt_edge(edges, ie0) ? ( get_interchunk_dvt_edge(edges, ie0) ? 'E' : '-' ) : '*')
                ,( is_a_dvt_edge(edges, ie0) ? ( get_intrachunk_dvt_edge(edges, ie0) ? 'I' : '-' ) : '*')
                ,( is_a_dvt_edge(edges, ie0) ? (get_reflected_edge(edges,ie0) ? 'R' : '-') : '*')
                
#endif                
#ifdef OUTPUT_DERIVED_QUANTITIES
		,nes
		,(invalid ? 'T' : 'F')
		,( get_blessed_edge(edges, ie0) ? 'b' : 'u' )
#endif                
                );
      }
    }}
  fclose(foutv);
  fclose(foute);
}



static void analyze_the_fragment_overlap_graph(FILE *fout,
                                               Tfragment frags[],
                                               Tedge edges[]) {

  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

  IntFragment_ID ifrag;
  IntFragment_ID
    n_raw_noncontained_nonspur_thru_frag = 0,
    n_raw_noncontained_nonspur_prefix_hanging_frag = 0,
    n_raw_noncontained_nonspur_suffix_hanging_frag = 0,
    n_raw_noncontained_nonspur_solo_frag = 0,

    n_raw_noncontained_spur_thru_frag = 0,
    n_raw_noncontained_spur_prefix_hanging_frag = 0,
    n_raw_noncontained_spur_suffix_hanging_frag = 0,
    n_raw_noncontained_spur_solo_frag = 0,

    n_raw_contained_nonspur_thru_frag = 0,
    n_raw_contained_nonspur_prefix_hanging_frag = 0,
    n_raw_contained_nonspur_suffix_hanging_frag = 0,
    n_raw_contained_nonspur_solo_frag = 0,

    n_raw_contained_spur_thru_frag = 0,
    n_raw_contained_spur_prefix_hanging_frag = 0,
    n_raw_contained_spur_suffix_hanging_frag = 0,
    n_raw_contained_spur_solo_frag = 0,

    n_cur_noncontained_nonspur_thru_frag = 0,
    n_cur_noncontained_nonspur_prefix_hanging_frag = 0,
    n_cur_noncontained_nonspur_suffix_hanging_frag = 0,
    n_cur_noncontained_nonspur_solo_frag = 0,

    n_cur_noncontained_spur_thru_frag = 0,
    n_cur_noncontained_spur_prefix_hanging_frag = 0,
    n_cur_noncontained_spur_suffix_hanging_frag = 0,
    n_cur_noncontained_spur_solo_frag = 0,

    n_cur_contained_nonspur_thru_frag = 0,
    n_cur_contained_nonspur_prefix_hanging_frag = 0,
    n_cur_contained_nonspur_suffix_hanging_frag = 0,
    n_cur_contained_nonspur_solo_frag = 0,

    n_cur_contained_spur_thru_frag = 0,
    n_cur_contained_spur_prefix_hanging_frag = 0,
    n_cur_contained_spur_suffix_hanging_frag = 0,
    n_cur_contained_spur_solo_frag = 0,

    n_as_fgb_deleted_frag = 0;

  IntFragment_ID
    raw_dvt_count = 0,
    raw_frc_count = 0,
    raw_toc_count = 0,
    cur_dvt_count = 0,
    cur_frc_count = 0,
    cur_toc_count = 0;

  for(ifrag=0; ifrag<nfrag; ifrag++){
    const int deleted   = (TRUE == get_del_fragment(frags,ifrag));
    const int contained = (TRUE == get_con_fragment(frags,ifrag));
    const int spur      = (TRUE == get_spur_fragment(frags,ifrag));

    const int raw_prefix_dvt_degree = get_raw_dvt_count_vertex(frags,ifrag,FALSE);
    const int raw_suffix_dvt_degree = get_raw_dvt_count_vertex(frags,ifrag,TRUE);
    const int raw_frc_degree = get_raw_frc_count_fragment(frags,ifrag);
    const int raw_toc_degree = get_raw_toc_count_fragment(frags,ifrag);
    const int cur_prefix_dvt_degree = get_seglen_dvt_vertex(frags,ifrag,FALSE);
    const int cur_suffix_dvt_degree = get_seglen_dvt_vertex(frags,ifrag,TRUE);
    const int cur_prefix_frc_degree = get_seglen_frc_vertex(frags,ifrag,FALSE);
    const int cur_suffix_frc_degree = get_seglen_frc_vertex(frags,ifrag,TRUE);

    const int raw_thru 
      = ((raw_prefix_dvt_degree != 0) && 
         (raw_suffix_dvt_degree != 0) );
    const int raw_prefix_hanging
      = ((raw_prefix_dvt_degree == 0) && 
         (raw_suffix_dvt_degree != 0) );
    const int raw_suffix_hanging 
      = ((raw_prefix_dvt_degree != 0) && 
         (raw_suffix_dvt_degree == 0) );
    const int raw_solo
      = ((raw_prefix_dvt_degree == 0) && 
         (raw_suffix_dvt_degree == 0) );
    // The raw counters affected by the deleted and removed by breaker flags on fragments.
    // The raw counters not affected by the spur flag and contained flags on fragments.
      
    const int cur_thru 
      = ((cur_prefix_dvt_degree != 0) && 
         (cur_suffix_dvt_degree != 0) );
    const int cur_prefix_hanging
      = ((cur_prefix_dvt_degree == 0) && 
         (cur_suffix_dvt_degree != 0) );
    const int cur_suffix_hanging 
      = ((cur_prefix_dvt_degree != 0) && 
         (cur_suffix_dvt_degree == 0) );
    const int cur_solo
      = ((cur_prefix_dvt_degree == 0) && 
         (cur_suffix_dvt_degree == 0) );

    if(!deleted) {
      raw_dvt_count += raw_prefix_dvt_degree + raw_suffix_dvt_degree;
      raw_frc_count += raw_frc_degree;
      raw_toc_count += raw_toc_degree;
      cur_dvt_count += cur_prefix_dvt_degree + cur_suffix_dvt_degree;
      cur_frc_count += cur_prefix_frc_degree + cur_suffix_frc_degree;
      cur_toc_count += 0;

      assert( !(raw_prefix_dvt_degree == 0) || (cur_prefix_dvt_degree == 0));
      assert( !(raw_suffix_dvt_degree == 0) || (cur_suffix_dvt_degree == 0));
      // If the fragment-end has a zero dvt degree in the raw graph,
      // then it has a zero dvt degree in the current graph.

      if(
         ((cur_prefix_dvt_degree == 0)&&(raw_prefix_dvt_degree > 0)) ||
         ((cur_suffix_dvt_degree == 0)&&(raw_suffix_dvt_degree > 0))
         ) {
        // If the fragment-end has a zero dvt degree in the current graph,
        // then it should have a zero dvt degree in the raw graph.
        const IntFragment_ID iid = get_iid_fragment(frags,ifrag);
        fprintf(stdout,"REAPER DVT DISCONNECT: " F_IID " %d %d %d %d %d %d %d\n",
                iid,
                deleted, spur, contained,
                raw_prefix_dvt_degree,
                raw_suffix_dvt_degree,
                cur_prefix_dvt_degree,
                cur_suffix_dvt_degree
                );
      }
    }
      
              
    assert(deleted || (raw_thru + raw_prefix_hanging + raw_suffix_hanging + raw_solo == 1));
    assert(deleted || (cur_thru + cur_prefix_hanging + cur_suffix_hanging + cur_solo == 1));
    assert(deleted || (cur_prefix_dvt_degree <= raw_prefix_dvt_degree));
    assert(deleted || (cur_suffix_dvt_degree <= raw_suffix_dvt_degree));

    {
      const IntFragment_ID iid = get_iid_fragment(frags,ifrag);
	
      if(!(deleted || contained || (!raw_thru) || cur_thru)) {
        fprintf(stdout,"REAPER THRU DISCONNECT: "
                "iid=" F_IID " deleted=%d contained=%d raw_thru=%d cur_thru=%d\n",
                iid, deleted, contained, raw_thru, cur_thru);
      }
    }

    assert(deleted || (!cur_thru) || raw_thru);
    // Assert if cur_thru then raw_thru
      
    if(!(deleted || !spur || !raw_thru)) {
      const IntFragment_ID iid = get_iid_fragment(frags,ifrag);
      fprintf(stdout,"REAPER DVT DISCONNECT: " F_IID " %d %d %d %d %d %d %d\n",
              iid,
              deleted, spur, contained,
              raw_prefix_dvt_degree,
              raw_suffix_dvt_degree,
              cur_prefix_dvt_degree,
              cur_suffix_dvt_degree
              );
    }

    assert(deleted || !spur || !cur_thru); // if spur, then not cur_thru.
      
    if(deleted) {
      n_as_fgb_deleted_frag++;
    } else {

      if(contained) {
        if(spur) {
          if(raw_thru)           { n_raw_contained_spur_thru_frag++;}
          if(raw_prefix_hanging) { n_raw_contained_spur_prefix_hanging_frag++;}
          if(raw_suffix_hanging) { n_raw_contained_spur_suffix_hanging_frag++;}
          if(raw_solo)           { n_raw_contained_spur_solo_frag++;}
        } else {
          if(raw_thru)           { n_raw_contained_nonspur_thru_frag++;}
          if(raw_prefix_hanging) { n_raw_contained_nonspur_prefix_hanging_frag++;}
          if(raw_suffix_hanging) { n_raw_contained_nonspur_suffix_hanging_frag++;}
          if(raw_solo)           { n_raw_contained_nonspur_solo_frag++;}
        }
      } else {
        if(spur) {
          if(raw_thru)           { n_raw_noncontained_spur_thru_frag++;}
          if(raw_prefix_hanging) { n_raw_noncontained_spur_prefix_hanging_frag++;}
          if(raw_suffix_hanging) { n_raw_noncontained_spur_suffix_hanging_frag++;}
          if(raw_solo)           { n_raw_noncontained_spur_solo_frag++;}
        } else {
          if(raw_thru)           { n_raw_noncontained_nonspur_thru_frag++;}
          if(raw_prefix_hanging) { n_raw_noncontained_nonspur_prefix_hanging_frag++;}
          if(raw_suffix_hanging) { n_raw_noncontained_nonspur_suffix_hanging_frag++;}
          if(raw_solo)           { n_raw_noncontained_nonspur_solo_frag++;}
        }
      }

      if(contained) {
        if(spur) {
          if(cur_thru)           { n_cur_contained_spur_thru_frag++;}
          if(cur_prefix_hanging) { n_cur_contained_spur_prefix_hanging_frag++;}
          if(cur_suffix_hanging) { n_cur_contained_spur_suffix_hanging_frag++;}
          if(cur_solo)           { n_cur_contained_spur_solo_frag++;}
        } else {
          if(cur_thru)           { n_cur_contained_nonspur_thru_frag++;}
          if(cur_prefix_hanging) { n_cur_contained_nonspur_prefix_hanging_frag++;}
          if(cur_suffix_hanging) { n_cur_contained_nonspur_suffix_hanging_frag++;}
          if(cur_solo)           { n_cur_contained_nonspur_solo_frag++;}
        }
      } else {
        if(spur) {
          if(cur_thru)           { n_cur_noncontained_spur_thru_frag++;}
          if(cur_prefix_hanging) { n_cur_noncontained_spur_prefix_hanging_frag++;}
          if(cur_suffix_hanging) { n_cur_noncontained_spur_suffix_hanging_frag++;}
          if(cur_solo)           { n_cur_noncontained_spur_solo_frag++;}
        } else {
          if(cur_thru)           { n_cur_noncontained_nonspur_thru_frag++;}
          if(cur_prefix_hanging) { n_cur_noncontained_nonspur_prefix_hanging_frag++;}
          if(cur_suffix_hanging) { n_cur_noncontained_nonspur_suffix_hanging_frag++;}
          if(cur_solo)           { n_cur_noncontained_nonspur_solo_frag++;}
        }
      }
    }
  }

  fprintf(fout,"Deleted Fragment types\n");
  fprintf(fout,
          "%15" F_IIDP " : total number of fragments\n"
          "---Deleted Classification:\n"
          "%15" F_IIDP " : number of deleted fragments\n"
          "%15" F_IIDP " : number of non-deleted fragments\n"
          ,
          nfrag, 
          n_as_fgb_deleted_frag,
          nfrag - n_as_fgb_deleted_frag
          );

  fprintf(fout,"Non-Deleted Fragment types\n");

  //  The Raw overlap classifications:
  fprintf(fout,
          "---Raw Dovetail Degree Classification:\n"
          "%15" F_IIDP " : solo           fragments\n"
          "%15" F_IIDP " : prefix hanging fragments\n"
          "%15" F_IIDP " : suffix hanging fragments\n"
          "%15" F_IIDP " : thru           fragments\n"
          ,

          n_raw_noncontained_nonspur_solo_frag+
          n_raw_contained_nonspur_solo_frag+
          n_raw_noncontained_spur_solo_frag+
          n_raw_contained_spur_solo_frag,

          n_raw_noncontained_nonspur_prefix_hanging_frag+
          n_raw_contained_nonspur_prefix_hanging_frag+
          n_raw_noncontained_spur_prefix_hanging_frag+
          n_raw_contained_spur_prefix_hanging_frag,

          n_raw_noncontained_nonspur_suffix_hanging_frag+
          n_raw_contained_nonspur_suffix_hanging_frag+
          n_raw_noncontained_spur_suffix_hanging_frag+
          n_raw_contained_spur_suffix_hanging_frag,

          n_raw_noncontained_nonspur_thru_frag+
          n_raw_contained_nonspur_thru_frag+
          n_raw_noncontained_spur_thru_frag+
          n_raw_contained_spur_thru_frag

          );
  fprintf(fout,
          "---Raw Dovetail Degree and Spur Classification:\n"
          "%15" F_IIDP " : non-spur solo           fragments\n"
          "%15" F_IIDP " : non-spur prefix hanging fragments\n"
          "%15" F_IIDP " : non-spur suffix hanging fragments\n"
          "%15" F_IIDP " : non-spur thru           fragments\n"
          "%15" F_IIDP " :     spur solo           fragments\n"
          "%15" F_IIDP " :     spur prefix hanging fragments\n"
          "%15" F_IIDP " :     spur suffix hanging fragments\n"
          "%15" F_IIDP " :     spur thru           fragments\n"
          ,

          n_raw_noncontained_nonspur_solo_frag+
          n_raw_contained_nonspur_solo_frag,
          n_raw_noncontained_nonspur_prefix_hanging_frag+
          n_raw_contained_nonspur_prefix_hanging_frag,
          n_raw_noncontained_nonspur_suffix_hanging_frag+
          n_raw_contained_nonspur_suffix_hanging_frag,
          n_raw_noncontained_nonspur_thru_frag+
          n_raw_contained_nonspur_thru_frag,

          n_raw_noncontained_spur_solo_frag+
          n_raw_contained_spur_solo_frag,
          n_raw_noncontained_spur_prefix_hanging_frag+
          n_raw_contained_spur_prefix_hanging_frag,
          n_raw_noncontained_spur_suffix_hanging_frag+
          n_raw_contained_spur_suffix_hanging_frag,
          n_raw_noncontained_spur_thru_frag+
          n_raw_contained_spur_thru_frag
          );
  fprintf(fout,
          "---Raw Dovetail Degree, Spur, and Containment Degree Classification:\n"
          "%15" F_IIDP " : non-contained non-spur solo           fragments\n"
          "%15" F_IIDP " : non-contained non-spur prefix hanging fragments\n"
          "%15" F_IIDP " : non-contained non-spur suffix hanging fragments\n"
          "%15" F_IIDP " : non-contained non-spur thru           fragments\n"
          "%15" F_IIDP " : non-contained     spur solo           fragments\n"
          "%15" F_IIDP " : non-contained     spur prefix hanging fragments\n"
          "%15" F_IIDP " : non-contained     spur suffix hanging fragments\n"
          "%15" F_IIDP " : non-contained     spur thru           fragments\n"
          "%15" F_IIDP " :     contained non-spur solo           fragments\n"
          "%15" F_IIDP " :     contained non-spur prefix hanging fragments\n"
          "%15" F_IIDP " :     contained non-spur suffix hanging fragments\n"
          "%15" F_IIDP " :     contained non-spur thru           fragments\n"
          "%15" F_IIDP " :     contained     spur solo           fragments\n"
          "%15" F_IIDP " :     contained     spur prefix hanging fragments\n"
          "%15" F_IIDP " :     contained     spur suffix hanging fragments\n"
          "%15" F_IIDP " :     contained     spur thru           fragments\n"
          ,

          n_raw_noncontained_nonspur_solo_frag,
          n_raw_noncontained_nonspur_prefix_hanging_frag,
          n_raw_noncontained_nonspur_suffix_hanging_frag,
          n_raw_noncontained_nonspur_thru_frag,

          n_raw_noncontained_spur_solo_frag,
          n_raw_noncontained_spur_prefix_hanging_frag,
          n_raw_noncontained_spur_suffix_hanging_frag,
          n_raw_noncontained_spur_thru_frag,

          n_raw_contained_nonspur_solo_frag,
          n_raw_contained_nonspur_prefix_hanging_frag,
          n_raw_contained_nonspur_suffix_hanging_frag,
          n_raw_contained_nonspur_thru_frag,

          n_raw_contained_spur_solo_frag,
          n_raw_contained_spur_prefix_hanging_frag,
          n_raw_contained_spur_suffix_hanging_frag,
          n_raw_contained_spur_thru_frag
          );


  //  The Current overlap classifications:
  fprintf(fout,
          "---Current Dovetail Degree Classification:\n"
          "%15" F_IIDP " : solo           fragments\n"
          "%15" F_IIDP " : prefix hanging fragments\n"
          "%15" F_IIDP " : suffix hanging fragments\n"
          "%15" F_IIDP " : thru           fragments\n"
          ,

          n_cur_noncontained_nonspur_solo_frag+
          n_cur_contained_nonspur_solo_frag+
          n_cur_noncontained_spur_solo_frag+
          n_cur_contained_spur_solo_frag,

          n_cur_noncontained_nonspur_prefix_hanging_frag+
          n_cur_contained_nonspur_prefix_hanging_frag+
          n_cur_noncontained_spur_prefix_hanging_frag+
          n_cur_contained_spur_prefix_hanging_frag,

          n_cur_noncontained_nonspur_suffix_hanging_frag+
          n_cur_contained_nonspur_suffix_hanging_frag+
          n_cur_noncontained_spur_suffix_hanging_frag+
          n_cur_contained_spur_suffix_hanging_frag,

          n_cur_noncontained_nonspur_thru_frag+
          n_cur_contained_nonspur_thru_frag+
          n_cur_noncontained_spur_thru_frag+
          n_cur_contained_spur_thru_frag

          );
  fprintf(fout,
          "---Current Dovetail Degree and Spur Classification:\n"
          "%15" F_IIDP " : non-spur solo           fragments\n"
          "%15" F_IIDP " : non-spur prefix hanging fragments\n"
          "%15" F_IIDP " : non-spur suffix hanging fragments\n"
          "%15" F_IIDP " : non-spur thru           fragments\n"
          "%15" F_IIDP " :     spur solo           fragments\n"
          "%15" F_IIDP " :     spur prefix hanging fragments\n"
          "%15" F_IIDP " :     spur suffix hanging fragments\n"
          "%15" F_IIDP " :     spur thru           fragments\n"
          ,

          n_cur_noncontained_nonspur_solo_frag+
          n_cur_contained_nonspur_solo_frag,
          n_cur_noncontained_nonspur_prefix_hanging_frag+
          n_cur_contained_nonspur_prefix_hanging_frag,
          n_cur_noncontained_nonspur_suffix_hanging_frag+
          n_cur_contained_nonspur_suffix_hanging_frag,
          n_cur_noncontained_nonspur_thru_frag+
          n_cur_contained_nonspur_thru_frag,

          n_cur_noncontained_spur_solo_frag+
          n_cur_contained_spur_solo_frag,
          n_cur_noncontained_spur_prefix_hanging_frag+
          n_cur_contained_spur_prefix_hanging_frag,
          n_cur_noncontained_spur_suffix_hanging_frag+
          n_cur_contained_spur_suffix_hanging_frag,
          n_cur_noncontained_spur_thru_frag+
          n_cur_contained_spur_thru_frag
          );
  fprintf(fout,
          "---Current Dovetail Degree, Spur, and Containment Degree Classification:\n"
          "%15" F_IIDP " : non-contained non-spur solo           fragments\n"
          "%15" F_IIDP " : non-contained non-spur prefix hanging fragments\n"
          "%15" F_IIDP " : non-contained non-spur suffix hanging fragments\n"
          "%15" F_IIDP " : non-contained non-spur thru           fragments\n"
          "%15" F_IIDP " : non-contained     spur solo           fragments\n"
          "%15" F_IIDP " : non-contained     spur prefix hanging fragments\n"
          "%15" F_IIDP " : non-contained     spur suffix hanging fragments\n"
          "%15" F_IIDP " : non-contained     spur thru           fragments\n"
          "%15" F_IIDP " :     contained non-spur solo           fragments\n"
          "%15" F_IIDP " :     contained non-spur prefix hanging fragments\n"
          "%15" F_IIDP " :     contained non-spur suffix hanging fragments\n"
          "%15" F_IIDP " :     contained non-spur thru           fragments\n"
          "%15" F_IIDP " :     contained     spur solo           fragments\n"
          "%15" F_IIDP " :     contained     spur prefix hanging fragments\n"
          "%15" F_IIDP " :     contained     spur suffix hanging fragments\n"
          "%15" F_IIDP " :     contained     spur thru           fragments\n"
          ,

          n_cur_noncontained_nonspur_solo_frag,
          n_cur_noncontained_nonspur_prefix_hanging_frag,
          n_cur_noncontained_nonspur_suffix_hanging_frag,
          n_cur_noncontained_nonspur_thru_frag,

          n_cur_noncontained_spur_solo_frag,
          n_cur_noncontained_spur_prefix_hanging_frag,
          n_cur_noncontained_spur_suffix_hanging_frag,
          n_cur_noncontained_spur_thru_frag,

          n_cur_contained_nonspur_solo_frag,
          n_cur_contained_nonspur_prefix_hanging_frag,
          n_cur_contained_nonspur_suffix_hanging_frag,
          n_cur_contained_nonspur_thru_frag,

          n_cur_contained_spur_solo_frag,
          n_cur_contained_spur_prefix_hanging_frag,
          n_cur_contained_spur_suffix_hanging_frag,
          n_cur_contained_spur_thru_frag
          );

  //  Totals:
  fprintf(fout,
          "Fragment-end overlap edge degree counts\n"
          "%15" F_IIDP " : raw dovetail edges\n"
          "%15" F_IIDP " : raw from-contained edges\n"
          "%15" F_IIDP " : raw to-contained edges\n"
          "%15" F_IIDP " : current dovetail edges\n"
          "%15" F_IIDP " : current from-contained edges\n"
          "%15" F_IIDP " : current to-contained edges\n"
          ,
          raw_dvt_count,
          raw_frc_count,
          raw_toc_count,
          cur_dvt_count,
          cur_frc_count,
          cur_toc_count );

  {
    IntEdge_ID iedge;
    size_t 
      n_as_cgb_dovetail=0,
      n_as_cgb_thickest=0,
      n_as_cgb_interchunk=0,
      n_as_cgb_intrachunk=0,
      n_as_cgb_touches_contained=0,
      n_as_cgb_between_contained=0,
      n_as_cgb_containment_frc=0,
      n_as_cgb_containment_toc=0,
      n_as_cgb_touches_crappy_dvt=0,
      n_as_cgb_touches_crappy_frc=0,
      n_as_cgb_touches_crappy_toc=0,

      n_as_cgb_between_crappy_dvt=0,
      n_as_cgb_between_crappy_toc=0,
      n_as_cgb_between_crappy_frc=0,
        
      n_as_cgb_marked_by_branch_dvt=0,
      n_as_cgb_removed_by_transitivity_dvt=0,
      n_as_cgb_removed_by_transitivity_frc=0,
      n_as_cgb_removed_by_transitivity_toc=0,

      n_as_cgb_removed_by_threshold_dvt=0,
      n_as_cgb_removed_by_threshold_frc=0,
      n_as_cgb_removed_by_threshold_toc=0,
      n_as_cgb_removed_by_duplicate_dvt=0,
      n_as_cgb_removed_by_duplicate_frc=0,
      n_as_cgb_removed_by_duplicate_toc=0;
    
    for(iedge=0; iedge<nedge; iedge++){
      const Tnes ines = get_nes_edge(edges,iedge);
      
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
          n_as_cgb_touches_contained++; break;
        case AS_CGB_BETWEEN_CONTAINED_EDGE:
          n_as_cgb_between_contained++; break;
        case AS_CGB_TOUCHES_CRAPPY_DVT:
          n_as_cgb_touches_crappy_dvt++; break;
        case AS_CGB_BETWEEN_CRAPPY_DVT:
          n_as_cgb_between_crappy_dvt++; break;

        case AS_CGB_TOUCHES_CRAPPY_CON:
          if(is_a_frc_edge(edges,iedge)) {
            n_as_cgb_touches_crappy_frc++;
          } else {
            n_as_cgb_touches_crappy_toc++;
          }
          break;
        case AS_CGB_BETWEEN_CRAPPY_CON:
          if(is_a_frc_edge(edges,iedge)) {
            n_as_cgb_between_crappy_frc++;
          } else {
            n_as_cgb_between_crappy_toc++;
          }
          break;
        
        case AS_CGB_CONTAINED_EDGE:
          if(is_a_frc_edge(edges,iedge)) {
            n_as_cgb_containment_frc++; break;
          } else {
            n_as_cgb_containment_toc++; break;
          }
          break;

        case AS_CGB_MARKED_BY_BRANCH_DVT:
          n_as_cgb_marked_by_branch_dvt++; break;

        case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
          n_as_cgb_removed_by_transitivity_dvt++; break;
        case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
          if(is_a_frc_edge(edges,iedge)) {
            n_as_cgb_removed_by_transitivity_frc++;
          } else {
            n_as_cgb_removed_by_transitivity_toc++;
          }
          break;

        case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
          n_as_cgb_removed_by_duplicate_dvt++; break;
        case AS_CGB_REMOVED_BY_DUPLICATE_CON:
          if(is_a_frc_edge(edges,iedge)) {
            n_as_cgb_removed_by_duplicate_frc++;
          } else {
            n_as_cgb_removed_by_duplicate_toc++;
          }
          break;

        default:
          fprintf(stderr,"FGA: unexpected overlap edge type %d\n",ines);
          assert(FALSE);
      }
    }

    fprintf(fout,"Overlap types\n");
    fprintf(fout,
	    "%15" F_IIDP " : total number of fragment overlaps\n"
	    "%15" F_SIZE_TP " : dovetail \n"
	    "%15" F_SIZE_TP " : dovetail interchunk\n"
	    "%15" F_SIZE_TP " : dovetail intrachunk\n"
	    "%15" F_SIZE_TP " : dovetail touches contained\n"
	    "%15" F_SIZE_TP " : dovetail between contained\n"
	    "%15" F_SIZE_TP " : dovetail touches crappy\n"
	    "%15" F_SIZE_TP " : containment touches crappy\n"
	    "%15" F_SIZE_TP " : dovetail between crappy\n"
	    "%15" F_SIZE_TP " : containment between crappy\n"
	    "%15" F_SIZE_TP " : dovetail marked by branch\n"
	    "%15" F_SIZE_TP " : dovetail removed by transitivity\n"
	    "%15" F_SIZE_TP " : dovetail removed by threshold\n"
	    "%15" F_SIZE_TP " : pseudo-dovetail (containment)\n"
	    "%15" F_SIZE_TP " : pseudo-dovetail removed by transitivity\n"
	    "%15" F_SIZE_TP " : pseudo-dovetail removed by threshold\n"
	    "%15" F_SIZE_TP " : dovetail removed as duplicate\n"
	    "%15" F_SIZE_TP " : pseudo-dovetail removed as duplicate\n",
	    nedge/2,
	    n_as_cgb_dovetail/2,
	    n_as_cgb_interchunk/2,
	    n_as_cgb_intrachunk/2,
	    n_as_cgb_touches_contained/2,
	    n_as_cgb_between_contained/2,
            n_as_cgb_touches_crappy_dvt/2,
            (n_as_cgb_touches_crappy_frc+
             n_as_cgb_touches_crappy_toc)/2,
            n_as_cgb_between_crappy_dvt/2,
            (n_as_cgb_between_crappy_toc+
             n_as_cgb_between_crappy_frc)/2,
	    n_as_cgb_marked_by_branch_dvt/2,
	    n_as_cgb_removed_by_transitivity_dvt/2,
	    n_as_cgb_removed_by_threshold_dvt/2,
	    (n_as_cgb_containment_frc+
	     n_as_cgb_containment_toc)/2,
	    (n_as_cgb_removed_by_transitivity_frc+
	     n_as_cgb_removed_by_transitivity_toc)/2,
	    (n_as_cgb_removed_by_threshold_frc+
	     n_as_cgb_removed_by_threshold_toc)/2,
	    n_as_cgb_removed_by_duplicate_dvt/2,
	    (n_as_cgb_removed_by_duplicate_frc+
             n_as_cgb_removed_by_duplicate_toc)/2
	    );
  }

  {
    IntFragment_ID ifrag;
    const int nsample=500;
    const int nbucket=500;
    Histogram_t 
      * raw_dvt_count_histogram
      = create_histogram(nsample,nbucket,0,TRUE),
      * raw_toc_count_histogram
      = create_histogram(nsample,nbucket,0,TRUE),
      * raw_frc_count_histogram
      = create_histogram(nsample,nbucket,0,TRUE);
    for(ifrag=0;ifrag<nfrag;ifrag++) {
      { 
        add_to_histogram(raw_dvt_count_histogram,
                         get_raw_dvt_count_vertex(frags,ifrag,FALSE), NULL);
        add_to_histogram(raw_dvt_count_histogram,
                         get_raw_dvt_count_vertex(frags,ifrag,TRUE), NULL);
        add_to_histogram(raw_toc_count_histogram,
                         get_raw_toc_count_fragment(frags,ifrag), NULL);
        add_to_histogram(raw_frc_count_histogram,
                         get_raw_frc_count_fragment(frags,ifrag), NULL);
      }
    }
    fprintf(fout,"\n\nHistogram of the raw dvt degree of the fragment-ends\n");
    print_histogram(fout,raw_dvt_count_histogram, 0, 1);
    free_histogram(raw_dvt_count_histogram);

    fprintf(fout,"\n\nHistogram of the raw toc degree of the fragment-ends\n");
    print_histogram(fout,raw_toc_count_histogram, 0, 1);
    free_histogram(raw_toc_count_histogram);

    fprintf(fout,"\n\nHistogram of the raw frc degree of the fragment-ends\n");
    print_histogram(fout,raw_frc_count_histogram, 0, 1);
    free_histogram(raw_frc_count_histogram);
  }

  {
    IntFragment_ID ifrag;
    const int nsample=500;
    const int nbucket=500;
    Histogram_t 
      *edges_per_vertex_histogram
      = create_histogram(nsample,nbucket,0,TRUE);
    for(ifrag=0;ifrag<nfrag;ifrag++) {
      { 
        const int deleted = get_del_fragment(frags,ifrag);
        const int contained = get_con_fragment(frags,ifrag);
        const int spur = get_spur_fragment(frags,ifrag);
        int isuffix;
        if((!deleted)&&(!contained)&&(!spur)) for(isuffix=0;isuffix<2;isuffix++) {
	    int nnode;
	    nnode = get_seglen_dvt_vertex(frags,ifrag,isuffix);
	    add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	  }
      }
    }
    fprintf(fout,"\n\nHistogram of the dovetail degree of the non-contained non-spur fragment-ends\n");
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
      const int deleted = get_del_fragment(frags,ifrag);
      const int contained = get_con_fragment(frags,ifrag);
      const int spur = get_spur_fragment(frags,ifrag);
      int isuffix;
      if((!deleted)&&(!contained)&&(spur)) for(isuffix=0;isuffix<2;isuffix++) {
	  const IntEdge_ID nnode = get_seglen_dvt_vertex(frags,ifrag,isuffix);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
    }
    fprintf(fout,"\n\nHistogram of the dovetail degree of the non-contained spur fragment-ends\n");
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
      const int deleted = get_del_fragment(frags,ifrag);
      const int contained = get_con_fragment(frags,ifrag);
      const int spur = get_spur_fragment(frags,ifrag);
      int isuffix;
      if((!deleted)&&(contained)&&(!spur)) for(isuffix=0;isuffix<2;isuffix++) {
	  const IntEdge_ID nnode = get_seglen_dvt_vertex(frags,ifrag,isuffix);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
    }
    fprintf(fout,"\n\nHistogram of the dovetail degree of the contained non-spur fragment-ends\n");
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
      const int deleted = get_del_fragment(frags,ifrag);
      const int contained = get_con_fragment(frags,ifrag);
      const int spur = get_spur_fragment(frags,ifrag);
      int isuffix;
      if((!deleted)&&(contained)&&(spur)) for(isuffix=0;isuffix<2;isuffix++) {
	  const IntEdge_ID nnode = get_seglen_dvt_vertex(frags,ifrag,isuffix);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
    }
    fprintf(fout,"\n\nHistogram of the dovetail degree of the contained spur fragment-ends\n");
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
      const int deleted = get_del_fragment(frags,ifrag);
      const int contained = get_con_fragment(frags,ifrag);
      const int spur = get_spur_fragment(frags,ifrag);
      int isuffix;
      if((!deleted)&&(!contained)&&(!spur)) for(isuffix=0;isuffix<2;isuffix++) {
	  const IntEdge_ID nnode = get_seglen_frc_vertex(frags,ifrag,isuffix);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
    }
    fprintf(fout,"\n\nHistogram of the from-contained degree of the non-contained non-spur fragment-ends\n");
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
      const int deleted = get_del_fragment(frags,ifrag);
      const int contained = get_con_fragment(frags,ifrag);
      const int spur = get_spur_fragment(frags,ifrag);
      int isuffix;
      if((!deleted)&&(!contained)&&(spur)) for(isuffix=0;isuffix<2;isuffix++) {
	  const IntEdge_ID nnode = get_seglen_frc_vertex(frags,ifrag,isuffix);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
    }
    fprintf(fout,"\n\nHistogram of the from-contained degree of the non-contained spur fragment-ends\n");
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
      const int deleted = get_del_fragment(frags,ifrag);
      const int contained = get_con_fragment(frags,ifrag);
      const int spur = get_spur_fragment(frags,ifrag);
      int isuffix;
      if((!deleted)&&(contained)&&(!spur)) for(isuffix=0;isuffix<2;isuffix++) {
	  const IntEdge_ID nnode = get_seglen_frc_vertex(frags,ifrag,isuffix);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
    }
    fprintf(fout,"\n\nHistogram of the from-contained degree of the contained non-spur fragment-ends\n");
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
      const int deleted = get_del_fragment(frags,ifrag);
      const int contained = get_con_fragment(frags,ifrag);
      const int spur = get_spur_fragment(frags,ifrag);
      int isuffix;
      if((!deleted)&&(contained)&&(spur)) for(isuffix=0;isuffix<2;isuffix++) {
	  const IntEdge_ID nnode = get_seglen_frc_vertex(frags,ifrag,isuffix);
	  add_to_histogram(edges_per_vertex_histogram, nnode, NULL);
	}
    }
    fprintf(fout,"\n\nHistogram of the from-contained degree of the contained spur fragment-ends\n");
    print_histogram(fout,edges_per_vertex_histogram, 0, 1);
    free_histogram(edges_per_vertex_histogram);
  }
}




void fragment_graph_analysis(Tfragment frags[],
                             Tedge     edges[],
                             FILE      *ffga) {

  const IntFragment_ID nfrag = GetNumFragments(frags);

  assert(frags      != NULL);
  assert(edges      != NULL);

  analyze_the_fragment_overlap_graph ( ffga, frags, edges);

  IntEdge_ID ie1;
  const IntEdge_ID nedge = GetNumEdges(edges);
  const int nsample=500;
  const int nbucket=500;
  Histogram_t 
    *edges_locality_histogram
    = create_histogram(nsample,nbucket,0,TRUE);
  for(ie1=0; ie1 < nedge; ie1++) {
    const Tnes nes1 = get_nes_edge(edges,ie1);
    if((nes1 != AS_CGB_REMOVED_BY_DUPLICATE_DVT) &&
       (nes1 != AS_CGB_REMOVED_BY_DUPLICATE_CON)
       ) {
      const IntFragment_ID iv0 = get_avx_edge(edges,ie1);
      const IntFragment_ID iv1 = get_bvx_edge(edges,ie1);
      const int iv1_iv0 = iv1 - iv0;
      const int vdiff = (iv1_iv0 > 0 ? iv1_iv0 : -iv1_iv0);
      add_to_histogram(edges_locality_histogram, vdiff, NULL);
    }
  }
  fprintf(ffga,"\n\nHistogram of the locality of graph edges\n");
  print_histogram(ffga,edges_locality_histogram, 0, 1);
  free_histogram(edges_locality_histogram);
}
