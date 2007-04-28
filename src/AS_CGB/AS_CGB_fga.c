
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
= "$Id: AS_CGB_fga.c,v 1.9 2007-04-28 08:46:21 brianwalenz Exp $";
/*********************************************************************
 *
 * Module: AS_CGB_fga.c
 * 
 * Description: A fragment graph analyser. This functional unit computes
 * whether fragment overlaps are valid.
 *
 * Author: Clark M. Mobarry
 ********************************************************************/

#include "AS_CGB_all.h"

#define DEBUGGING
#undef DEBUGGING

#define OUTPUT_DERIVED_QUANTITIES
#undef HIDE_OVERHANGS
#define DOVETAIL_OR_CONTAINED_FLAG

#ifdef DEBUGGING
#undef DEBUG14
#undef DEBUG15
#define DEBUG27
#endif // DEBUGGING

#undef DEBUG54


static void processOneFragSourceString
(const IntFragment_ID frgiid,
 const char * const so,
 Afraginfo  * const fraginfo)
{

  /* Extract interval info from source string.  */
  char s[2048]={0};
  int HasAnnotation=FALSE;
  char *interval = NULL;

  /* the output of this routine: */
  BPTYPE genbgn=0,genend=0;
  CDS_COORD_t frglen=0; /* The length in bps of the fragment in the genome. */

  assert(strlen(so)<1023);
  strcpy(s,so); /* strtok() modifies the input string. */
  // fprintf(stderr,"Process strlen=%ld,<%s>\n",strlen(s),s);
  if(strlen(s)>0) {
    char *fragment_name = NULL, *genome_interval_source = NULL;

    /* The first line is the fragment number [fr ] */
    fragment_name = strtok(s,"\n\r");
    if(fragment_name != NULL) {
      // fprintf(stderr,"strlen=%ld,<%s>\n",
      // strlen(fragment_name),fragment_name);
      /* The second line is the genome interval info for the fragment. */
      genome_interval_source = strtok(NULL,"\n\r ");
      if(genome_interval_source != NULL) {
	// fprintf(stderr,"strlen=%ld,<%s>\n",
	// strlen(genome_interval_source),genome_interval_source);
	// genome_interval_source = strstr(s,"\n");
	// fprintf(stderr,"%ld\n",(long)(genome_interval_source-s));
	// assert(genome_interval_source);
	// Find the genome_interval_source info in the source string
	// interval = strstr(genome_interval_source,"[");
	interval = strstr(genome_interval_source,"[");
	if(genome_interval_source != NULL) {
	  //printf("interval=%p\n",interval);
	  HasAnnotation = TRUE;
	}
      }
    }
  }
  if(HasAnnotation) {
    FILE * f_simbpts = NULL; // The file of simulator induced branch-points.
    char *repeat_source = NULL;
    
    if(interval==NULL){
      fprintf(stderr,"Invalid source field for Fragment -- %s\n", s);
      assert(FALSE);
    }else {
      int iret2;
      iret2 = sscanf(interval,
		     "[" BPFORMAT "," BPFORMAT "]",
		     &(genbgn), &(genend));
      if(iret2 != 2 ){
	fprintf(stderr,
		"Warning: Couldn't parse Interval for Fragment -- %s\n",
		interval);
	assert(FALSE);
      }
      // fprintf(stderr,"processFragSource: ifrag,genbgn,genend=%d,%ld,%ld\n",
      // ifrag,genbgn,genend);

      frglen = genend-genbgn;
      frglen = (frglen > 0 ? frglen : -frglen);
    }


    if(f_simbpts != NULL){
      const int MINIMUM_OFFSET1=0;

      char pre_let='@',suf_let='@';
      int 
	pre_ins = 0,
	pre_brp = 0, 
	pre_end = 0,
	suf_ins = 0,
	suf_brp = 0,
	suf_end = 0;
      
      pre_brp = frglen;
      suf_end = frglen;
      
      // repeats = strstr(s,"]"); repeats +=2;
      if(frglen > 0) {
	char letter= (char)64;
	int  instance=0, copynum=0, idontknow=0;
	int repbgn=0,repend=0,frgbgn=0,frgend=0;
	int iret;

	/*
	  if (pre_brp >= pre_end ) then there is not a prefix branch point.
	  if (suf_brp <= suf_end ) then there is not a suffix branch point.
	*/
	
	while((repeat_source = strtok(NULL,"\n\r")) != NULL) {
	  /* The remaining lines are annotations about repeat regions
	     that overlap this fragment before sequencing errors and 
	     polymorphisms. (Created in AS_SIM_labelfrag.pl).
	     
	     The following fields describe the intra-fragment overlaps between
	     a repeat region of the genome and this fragment (before 
	     sampling errors).
	     
	     letter   = [A..Z] a labelled repeat region type of the genome
	     instance = the can be more than on instance of each type
	     (repbgn,repend) = the (5,3)-prime side of the repeat region
	     overlapping this fragment, in repeat coordinates.
	     (frgbgn,frgend) = the (5,3)-prime side of this fragment 
	     overlapping the repeat region, in fragment coordinates.
	     
	  */
	  iret = sscanf
	    (repeat_source,
	     "%c.%d.%d (%d) [%d,%d] [%d,%d]",
	     &(letter),&(instance),&(copynum),&(idontknow),
	     &(repbgn),&(repend),&(frgbgn),&(frgend));
	  assert(iret == 8);

#ifdef DEBUG14
	  fprintf(stdout,"ZIP: [" 
		  F_S64 "," F_S64
		  "] %c.%d.%d (%d) [%d,%d] [%d,%d]\n",
		  genbgn, genend,
		  letter,instance,copynum,idontknow,
		  repbgn,repend,frgbgn,frgend);
#endif
	  
	  if(frgbgn > frgend) {
	    int temp;
	    temp = frgbgn; frgbgn = frgend; frgend = temp;
	    temp = repbgn; repbgn = repend; repend = temp;
	  }

	  // (frgbgn,frgend) is the interval of the fragment covered
	  // by a paricular repeat copy.
	  
	  /* Check for branch points where the
	     repeat region extends to the end of the fragment. */
	  if(frgend >= frglen - MINIMUM_OFFSET1) { 
	    // The repeat region is overlaps the fragment-suffix.
	    pre_brp = MIN(pre_brp,frgbgn);
	    pre_end = MAX(pre_end,frgend);
	    pre_let = letter;
	    pre_ins = instance;
	    fprintf(f_simbpts,F_IID " %d %d %d\n",
                    frgiid, TRUE, (int)'S', frglen - pre_brp);
	  }
	  if(frgbgn <= MINIMUM_OFFSET1) {
	    // The repeat region is overlaps the fragment-suffix.
	    suf_brp = MAX(suf_brp,frgend);
	    suf_end = MIN(suf_end,frgbgn);
	    suf_let = letter;
	    suf_ins = instance;
	    fprintf(f_simbpts,F_IID " %d %d %d\n",
                    frgiid, FALSE, (int)'S', suf_brp);
	  }
	 
	}
	/*
	  To be a valid branch point we need:
	  (pre_brp < pre_end) or (suf_brp > suf_end)
	  A special case is when 
	  (pre_brp == 0) && (suf_brp == frglen)
	  which corresponds to a fragment contained in a repeat.
	*/
      }
    }
  }
 
  fraginfo->genbgn = genbgn;
  fraginfo->genend = genend;
}

#ifdef GENINFO
int check_overlap_with_simulator
(/* Input only */
 const IntFragment_ID nfrag,
 const Tfraginfo * const fraginfo,
 const IntFragment_ID iavx,
 const IntFragment_ID ibvx)
{
  /* Check whether two fragments overlap locally in the simulated
     genome.  The return value is ZERO if the overlap is a local
     overlap in the genome sequence.  The return value is 1 if the
     overlap is not a local overlap in the genome sequence.  */
  BPTYPE a_bgn,a_end,b_bgn,b_end;
  BPTYPE a_min,a_max,b_min,b_max,themin,themax;
  int iactual;

  // assert(iavx >= 0);
  // assert(ibvx >= 0);
  assert(iavx < nfrag);
  assert(ibvx < nfrag);
  
  a_bgn = get_genbgn_fraginfo(fraginfo,iavx);
  a_end = get_genend_fraginfo(fraginfo,iavx);
  b_bgn = get_genbgn_fraginfo(fraginfo,ibvx);
  b_end = get_genend_fraginfo(fraginfo,ibvx);
  
  /* For a linear genome, the fragments truly overlap if the length
     of a line segment that covers both of the fragments is less than
     the sum of the lengths of both fragments. */
  a_min = MIN(a_bgn,a_end);
  a_max = MAX(a_bgn,a_end);
  b_min = MIN(b_bgn,b_end);
  b_max = MAX(b_bgn,b_end);
  themin = MIN(a_min,b_min);
  themax = MAX(a_max,b_max);
  iactual = ((themax-themin) <= (a_max-a_min)+(b_max-b_min)); 
  /* The two fragments must overlap by at least one bp. But, I allow
     zero overlap as well, so that we can deal with no simulator
     information. In that case, a_min == a_max, etc.  */

  return (! iactual);
}

static int check_overlap_with_simulator2
(/* Input only */
 const IntFragment_ID nfrag,
 const Tfraginfo * const fraginfo,
 const IntFragment_ID iavx,
 const int iasx,
 const IntFragment_ID ibvx,
 const int ibsx)
{
  /* Check whether two fragments overlap locally in the simulated
     genome.  The return value is ZERO if the overlap is a local
     overlap in the genome sequence.  The return value is 1 if the
     overlap is not a local overlap in the genome sequence.  */
  BPTYPE a_bgn,a_end,b_bgn,b_end;
  BPTYPE a_min,a_max,b_min,b_max,themin,themax;
  int iactual;

  // assert(iavx >= 0);
  // assert(ibvx >= 0);
  assert(iavx < nfrag);
  assert(ibvx < nfrag);
  
  a_bgn = get_genbgn_fraginfo(fraginfo,iavx);
  a_end = get_genend_fraginfo(fraginfo,iavx);
  b_bgn = get_genbgn_fraginfo(fraginfo,ibvx);
  b_end = get_genend_fraginfo(fraginfo,ibvx);
  
  /* For a linear genome, the fragments truly overlap if the length
     of a line segment that covers both of the fragments is less than
     the sum of the lengths of both fragments. */
  a_min = MIN(a_bgn,a_end);
  a_max = MAX(a_bgn,a_end);
  b_min = MIN(b_bgn,b_end);
  b_max = MAX(b_bgn,b_end);
  themin = MIN(a_min,b_min);
  themax = MAX(a_max,b_max);
  iactual = ((themax-themin) < (a_max-a_min)+(b_max-b_min)); 
  /* The two fragments must overlap by at least one bp. */

  return (! iactual);
}
#endif /*GENINFO*/



static void processFragSourceString
(const VA_TYPE(char) * const fragsrc,
 const size_t isrc,
 Tfraginfo * const fraginfo, 
 const IntFragment_ID ifrag,
 const IntFragment_ID frgiid
 )
{
  /* Extract interval info from source string.  */
  char *so=NULL;
  Afraginfo the_fraginfo;

  so = GetVA_char(fragsrc,isrc); // Consult AS_UTL_Var.h
  assert(strlen(so)<1023);
  
  processOneFragSourceString(frgiid, so, &the_fraginfo);

#ifdef GENINFO
  set_genbgn_fraginfo(fraginfo,ifrag,the_fraginfo.genbgn);
  set_genend_fraginfo(fraginfo,ifrag,the_fraginfo.genend);
#endif /*GENINFO*/
}

static void mark_invalid_edges
(
 const IntFragment_ID nfrag,
 Tfraginfo *fraginfo,
 Tedge *edges
 )
{
  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID iedge;
  for(iedge=0;iedge<nedge;iedge++){
    const IntFragment_ID iavx = get_avx_edge(edges,iedge);
    const IntFragment_ID ibvx = get_bvx_edge(edges,iedge);
    // iasx = get_asx_edge(edges,iedge);
    // ibsx = get_bsx_edge(edges,iedge);
    // ines = get_nes_edge(edges,iedge);
    const int is_invalid = 
      check_overlap_with_simulator(nfrag,fraginfo,iavx,ibvx);
    set_inv_edge(edges,iedge,is_invalid);
  }
}


void view_fgb_chkpnt
(
 const char * const Store_File_Prefix,
 /* Input Only */
 Tfragment frags[], 
 Tedge edges[])
{
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
#ifndef OUTPUT_UID
              "iid %9" F_IIDP " "
#else
              "uid %16" F_UIDP " "
#endif
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
#ifndef OUTPUT_UID
              get_iid_fragment(frags,iv0),
#else
              get_uid_fragment(frags,iv0),
#endif	   
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
#ifndef OUTPUT_UID
      const IntFragment_ID afr = get_iid_fragment(frags,avx);
      const IntFragment_ID bfr = get_iid_fragment(frags,bvx);
#else
      const Fragment_ID afr = get_uid_fragment(frags,avx);
      const Fragment_ID bfr = get_uid_fragment(frags,bvx);
#endif
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
#ifndef OUTPUT_UID
		"%9" F_IIDP " "
#else // OUTPUT_UID
                "%16" F_UIDP " "
#endif // OUTPUT_UID
                "%2d "
#ifndef HIDE_OVERHANGS
                "%6d "
#endif                
#ifndef OUTPUT_UID
		"%9" F_IIDP " "
#else // OUTPUT_UID
                "%16 " F_UIDP " "
#endif // OUTPUT_UID
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
                //,( is_a_dvt_edge(edges, ie0) ? ( get_buddy_dvt_edge(edges, ie0) ? 'B' : '-' ) : '*')
                //,( is_a_dvt_edge(edges, ie0) ? ( get_tied_dvt_edge(edges,ie0) ? 'T' : '-' ) : '*')
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
#ifdef GENINFO
  { // 
    IntEdge_ID iedge;
    size_t
      count_invalid_dovetail_edges=0,
      count_invalid_thickest_edges=0,
      count_invalid_interchunk_edges=0,
      count_invalid_intrachunk_edges=0,
      count_invalid_containment_edges=0,
      count_invalid_touches_contained_edges=0,
      count_invalid_between_contained_edges=0,
      count_invalid_touches_crappy_dvt=0,
      count_invalid_between_crappy_dvt=0;
      
    for(iedge=0;iedge<nedge;iedge++){
      IntFragment_ID iavx,ibvx;
      int iasx,ibsx,iahg;
      Tnes ines;
      int invalid;
		
      iavx = get_avx_edge(edges,iedge);
      iasx = get_asx_edge(edges,iedge);
      ibvx = get_bvx_edge(edges,iedge);
      ibsx = get_bsx_edge(edges,iedge);

      ines = get_nes_edge(edges,iedge);
      iahg = get_ahg_edge(edges,iedge);

      invalid = get_inv_edge(edges,iedge);

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
	  case AS_CGB_CONTAINED_EDGE:
#ifdef DEBUGGING
	    printf("invalid_containment_edges: %d,%d %d,%d\n",
		   get_iid_fragment(frags,iavx),iasx,
		   get_iid_fragment(frags,ibvx),ibsx);
#endif // DEBUGGING
	    count_invalid_containment_edges++; break;

	  case AS_CGB_MARKED_BY_BRANCH_DVT:
	  case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
	  case AS_CGB_REMOVED_BY_THRESHOLD_DVT:
	  case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
	  case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
	  case AS_CGB_REMOVED_BY_THRESHOLD_CON:
	  case AS_CGB_REMOVED_BY_DUPLICATE_CON:
	    break;
	  default:
            //	    assert(FALSE);
            break;
        }
      }
    }
    fprintf(fout,"Fragment information from the simulator\n");

    fprintf(fout,
	    "%15" F_SIZE_TP " : number of invalid essential dovetail overlaps.\n"
	    "%15" F_SIZE_TP " : number of invalid essential interchunk dovetail overlaps.\n"
	    "%15" F_SIZE_TP " : number of invalid essential intrachunk dovetail overlaps.\n"
	    "%15" F_SIZE_TP " : number of invalid essential touches contained dovetail overlaps.\n"
	    "%15" F_SIZE_TP " : number of invalid essential between contained dovetail overlaps.\n"
	    "%15" F_SIZE_TP " : number of invalid essential touches crappy dovetail overlaps.\n"
	    "%15" F_SIZE_TP " : number of invalid essential between crappy dovetail overlaps.\n"
	    "%15" F_SIZE_TP " : number of invalid essential containment non-dovetail overlaps.\n"
            ,
	    count_invalid_dovetail_edges/2,
	    count_invalid_interchunk_edges/2,
	    count_invalid_intrachunk_edges/2,
	    count_invalid_touches_contained_edges/2,
	    count_invalid_between_contained_edges/2,
	    count_invalid_touches_crappy_dvt/2,
	    count_invalid_between_crappy_dvt/2,
	    count_invalid_containment_edges/2
            );
  }
#endif // GENINFO

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

void setup_fraginfo
(/* Input only */
 IntFragment_ID max_frag_iid,
 // The maximum fragment IID assigned by the Celera assembler
 // gatekeeper.
 Tfragment *frags,
 // The internal representation of fragments
 VA_TYPE(char) frag_annotations[],
 // The simulator fragment annotations.
 /* Output only */
 Tfraginfo *fraginfo
 )
{
  const IntFragment_ID nfrag = GetNumFragments(frags);
  // const IntEdge_ID nedge = GetNumEdges(edges);
  FragmentHashObject *afr_to_avx = create_FragmentHash((max_frag_iid+1));

#ifdef GENINFO
  /* Process the fragment source field. */ {
    IntFragment_ID ifrag;
    
#if 1
    fprintf(stderr,"FGA: max_frag_iid=" F_IID "\n", max_frag_iid);
    for(ifrag=0;ifrag<nfrag;ifrag++) { 
      const IntFragment_ID iid = get_iid_fragment(frags,ifrag);
      assert(iid >= 1);
      assert(iid <= max_frag_iid);
      set_vid_FragmentHash(afr_to_avx,iid,ifrag);
    }
#endif
    
    for(ifrag=0;ifrag<nfrag;ifrag++) { 
      const IntFragment_ID iid = get_iid_fragment(frags,ifrag);
      const size_t isrc = get_src_fragment(frags,ifrag);
      assert(iid >= 1);
      assert(iid <= max_frag_iid);
      set_vid_FragmentHash(afr_to_avx,iid,ifrag);
      processFragSourceString(frag_annotations,isrc,fraginfo,ifrag,iid);
    }
  }
  
  fprintf(stderr,"after processFragSourceString\n");

#endif // GENINFO
  destroy_FragmentHash(afr_to_avx);
}



void fragment_graph_analysis
(/* Input Only */
 const IntFragment_ID max_frag_iid,
 Tfragment frags[],
 // The internal representation of the fragment reads.
 Tedge     edges[],
 // The internal representation of the overlaps. 
 VA_TYPE(char) frag_annotations[], 
 // The simulator fragment annotations.
 const int ProcessFragmentAnnotationsForSimulatorCoordinates,
 /* Output only */
 FILE         *ffga
 ) 
{
  const IntFragment_ID nfrag = GetNumFragments(frags);
  Tfraginfo     *fraginfo   = CreateVA_Afraginfo(nfrag);

  assert(frags      != NULL);
  assert(edges      != NULL);
  assert(frag_annotations != NULL);
  assert(fraginfo != NULL);

  if( ProcessFragmentAnnotationsForSimulatorCoordinates ) {
    EnableRangeVA_Afraginfo(fraginfo,nfrag);
    setup_fraginfo( max_frag_iid, frags, frag_annotations, fraginfo);
    mark_invalid_edges( nfrag, fraginfo, edges);
  }
  
  analyze_the_fragment_overlap_graph
    ( ffga, max_frag_iid, frags, edges);

  { 
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
    fprintf(ffga,"\n\nHistogram of "
	    "the locality of graph edges\n");
    print_histogram(ffga,edges_locality_histogram, 0, 1);
    free_histogram(edges_locality_histogram);
    // fclose(fdiag);
  }
  DeleteVA_Afraginfo(fraginfo);

#ifdef DEBUG27
  view_fgb_chkpnt(frags, edges);
#endif // DEBUG27
}
