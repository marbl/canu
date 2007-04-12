
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
= "$Id: AS_CGB_iid2cid.c,v 1.10 2007-04-12 18:54:44 brianwalenz Exp $";
/* *******************************************************************
 *
 * Module: AS_CGB_fom2uom.c
 * 
 * Description: Chunk Graph Builder post-processor This functional
 * unit reads a *.cgc or *.cgi Celera Assembler i/o file an massages
 * into a *.cgb file.  The difference is that essential fragment
 * overlaps messages are converted into unitig overlap messages.
 *
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

/*********************************************************************/
/* System include files */
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*********************************************************************/
/* Local include files */
#include "AS_global.h"
#include "AS_UTL_version.h"  
#include "AS_MSG_pmesg.h"
#include "AS_CGB_histo.h"
#include "AS_CGB_all.h"

VA_DEF(int32)
VA_DEF(uint32)
VA_DEF(IntChunk_ID)

/****************************************************************************/
#define DEBUGGING

#define CMD_BUFFER_SIZE 1024

/*************************************************************************/

static void input_mesgs_pp
(int argc, char * argv [],
 FILE       *fcgi, 
 //FILE       *fcgb,
 int        *Pnadt,
 int        *Pnidt, 
 int        *Pnium,
 IntFragment_ID *Pnimp,
 IntEdge_ID *Pnuom,
 IntEdge_ID *Pnuom_dovetail, 
 IntEdge_ID *Pnuom_containment,
 IntChunk_ID *Pmin_unitig_cid,
 IntChunk_ID *Pmax_unitig_cid,
 IntFragment_ID *Pmin_frag_iid,
 IntFragment_ID *Pmax_frag_iid,
 int        analysis_flag,
 VA_TYPE(int32) * chunk_length,
 VA_TYPE(uint32) * chunk_num_frags,
 VA_TYPE(IntMultiPos) * imp_from_iid,
 VA_TYPE(IntChunk_ID) * cid_from_iid
 )
{ /* It is assumed that in the overlap records that new fragments
     point to old fragments.  */
  
  int nadt=0,nidt=0,nium=0,nimp=0;
  int nuom=0,nuom_dovetail=0,nuom_containment=0;
  GenericMesg *pmesg;

  IntFragment_ID lfrag = 0; // The global LID counter.
  
  const int nsample=500;
  const int nbucket=500;
  Histogram_t * uom_types_histogram
    = create_histogram(nsample,nbucket,TRUE,FALSE);

  while( EOF != ReadProtoMesg_AS(fcgi, &pmesg)) {
    const MessageType imesgtype = pmesg->t;
    
    switch(imesgtype) {
    case MESG_ADT: 
      {
	AuditMesg  *adt_mesg = (AuditMesg *)pmesg->m;
	nadt++;
#if 0
	AuditLine  auditLine;
	AppendAuditLine_AS(adt_mesg, &auditLine, time(NULL), "CGB_pp", 
			   CM_ID, "(empty)");
#else
	VersionStampADT(adt_mesg, argc, argv);
#endif
	//WriteProtoMesg_AS(fcgb,pmesg);
      }
      break;
    case MESG_IDT: 
      {
	InternalDistMesg  *idt_mesg;
	/*  Distance record--skip for now */
	idt_mesg = (InternalDistMesg *)(pmesg->m);
	nidt++;
	//WriteProtoMesg_AS(fcgb,pmesg);
      }
      break;
    case MESG_IUM: 
      {
	const IntUnitigMesg * const ium_mesg = (const IntUnitigMesg * const)(pmesg->m);
	// Do not change the value or address of the ium_mesg.
	const IntChunk_ID cid       = ium_mesg->iaccession;
	const IntFragment_ID  num_frags = ium_mesg->num_frags;
	const IntMultiPos * const f_list = ium_mesg->f_list;
	IntFragment_ID ifrag;

	const int unitig_len = ium_mesg->length;
	const char *chunk_sequence = ium_mesg->consensus;
	const int gapped_len = strlen(chunk_sequence);
	int sum_of_gaps_left_of[gapped_len+1]; 
	int ungapped_len = unitig_len;
	
	if ( (nium) % 1000 == 0) {
	  fprintf(stderr, "Unitig " F_IID "\r",cid);
	  //sleep(1);
	  fflush(stderr);
	}
	
	// If there is no seqence string (gapped_len==0), then do not
	// bother looking for gaps in it!
	assert( (unitig_len == gapped_len) || 
		(gapped_len == 0) );

	if(unitig_len == gapped_len) { 
	  // Compute the number of gaps in the consensus sequence to
	  // the left of an interface. Interface i is to the right of
	  // character i.

	  int ii;
	  // The SeqInt coordinates are between charaters.
	  sum_of_gaps_left_of[0]=0;
	  // coordinate 0 is to the left of the first character!
	  for(ii=0;ii<gapped_len;ii++){
	    sum_of_gaps_left_of[ii+1] = sum_of_gaps_left_of[ii]
	      + (chunk_sequence[ii] == '-');
	  }
	  ungapped_len = gapped_len - sum_of_gaps_left_of[gapped_len];
	}

	(*Pmax_unitig_cid) = MAX((*Pmax_unitig_cid),cid);
	SetVA_int32(chunk_length,cid,&ungapped_len);
	SetVA_uint32(chunk_num_frags,cid,&num_frags);

	//EnableRangeVA_IntChunk_ID(cid_array,nfrag+nimp);
	for( ifrag=0; ifrag<num_frags; ifrag++) {
	  const IntFragment_ID iid = f_list[ifrag].ident;
	  IntMultiPos imp = f_list[ifrag];
	  
	  if(unitig_len == gapped_len) { 
	    imp.position.bgn -= sum_of_gaps_left_of[imp.position.bgn];
	    imp.position.end -= sum_of_gaps_left_of[imp.position.end];
	  }

	  (*Pmax_frag_iid) = MAX((*Pmax_frag_iid),iid);
	  //fprintf(stderr,"iid=" F_IID " max_frag_iid=" F_IID "\n",iid, (*Pmax_frag_iid));
	  SetVA_IntMultiPos(imp_from_iid,iid,&imp);
	  SetVA_IntChunk_ID(cid_from_iid,iid,&cid);

	  lfrag ++;
	}
	
	//WriteProtoMesg_AS(fcgb,pmesg);
	// pass through the Unitig message
	
	nimp += num_frags;
	
	//fprintf(stderr,"nium,num_frags,nimp = %d," F_IID ",%d," F_IID "\n",
	// nium, num_frags, nimp, lfrag);
      }
      nium ++;
      break;
    case MESG_UOM: 
#if 0
      {
	/* Convert a FOM message into a UOM message. */
	UnitigOverlapMesg cea;

	UnitigOverlapMesg *cp = (UnitigOverlapMesg *)(pmesg->m);
	const IntFragment_ID iavx = cp->chunk1;
	const IntFragment_ID ibvx = cp->chunk2;
	const ChunkOrientationType  iorient = cp->orient;
	const UnitigOverlapType overlap_type = cp->overlap_type;
	const int32 best_overlap_length = cp->best_overlap_length;
	const int32 min_overlap_length = cp->min_overlap_length;
	const int32 max_overlap_length = cp->max_overlap_length;
	const float32 quality = cp->quality;

	if ( (nuom) % 1000 == 0) {
	  //fflush(fcgb);
	  fprintf(stderr, "Unitig overlap message %d\r",nuom);
	  //sleep(1);
	  fflush(stderr);
	}
	
	if(analysis_flag) {
	  add_to_histogram(uom_types_histogram, (int)overlap_type, NULL);
	}
  
	// assert(cp->source == NULL);

	{ // The following are chunk based edge vertices.

	  const int iasx = ( (iorient == AB_AB) || 
			     (iorient == AB_BA) );
	  const int ibsx = ( (iorient == AB_BA) || 
			     (iorient == BA_BA) );

	  const IntChunk_ID cavx = *(GetVA_IntChunk_ID(cid_from_iid,iavx));
	  const IntChunk_ID cbvx = *(GetVA_IntChunk_ID(cid_from_iid,ibvx));
	  const int32 ao5p = GetVA_IntMultiPos(imp_from_iid,iavx)->position.bgn;
	  const int32 ao3p = GetVA_IntMultiPos(imp_from_iid,iavx)->position.end;
	  const int32 bo5p = GetVA_IntMultiPos(imp_from_iid,ibvx)->position.bgn;
	  const int32 bo3p = GetVA_IntMultiPos(imp_from_iid,ibvx)->position.end;

	  const int casx = (iasx ^ (ao3p < ao5p) );
	  const int cbsx = (ibsx ^ (bo3p < bo5p) );

	  const ChunkOrientationType  orient
	    = ( ( casx) && ( cbsx) ? AB_BA : 0)
	    + ( ( casx) && (!cbsx) ? AB_AB : 0)
	    + ( (!casx) && ( cbsx) ? BA_BA : 0)
	    + ( (!casx) && (!cbsx) ? BA_AB : 0);
	  
	  const BPTYPE pa0 = ( iasx ? ao3p : ao5p );
	  const BPTYPE pa1 = 
	    ( casx 
	      ? *(GetVA_int32(chunk_length,cavx)) - pa0
	      : pa0 );
	  // pa1 is the gapped distance from the chunk-end of the
	  // A-chunk in the overlap to the fragment-end of the
	  // A-fragment in the overlap.

	  const BPTYPE pb0 = ( ibsx ? bo3p : bo5p );
	  const BPTYPE pb1
	    = ( cbsx 
		? *(GetVA_int32(chunk_length,cbvx)) - pb0
		: pb0 );
	  // pb1 is the gapped distance from the chunk-end of the
	  // B-chunk in the overlap to the fragment-end of the
	  // B-fragment in the overlap.
	  
	  assert( 0 != orient );

#if 0
	  if( ((cavx == 1) && (cbvx == 41)) ||
	      ((cavx == 41) && (cbvx == 1)) ) {
	    fprintf(stderr,
		    "cavx=" F_IID ",casx=%d,iavx=" F_IID ",iasx=%d\n"
		    "cbvx=" F_IID ",cbsx=%d,ibvx=" F_IID ",ibsx=%d\n"
		    "orient=%d,overlap_type=%d\n",
		    cavx,casx,iavx,iasx,
		    cbvx,cbsx,ibvx,ibsx,
		    orient,overlap_type);
	  }
#endif	  
	  cea.chunk1 = cavx;
	  cea.chunk2 = cbvx;
	  cea.orient = orient;
	  cea.overlap_type = overlap_type;
	  cea.best_overlap_length = pa1 + pb1 + best_overlap_length;
	  cea.min_overlap_length  = pa1 + pb1 + min_overlap_length;
	  cea.max_overlap_length  = pa1 + pb1 + max_overlap_length;
	  cea.quality = quality;
	  cea.source = NULL;
	  

	} // Convert to chunk overlaps from fragment overlaps

	pmesg->t = MESG_UOM;
	pmesg->m = &cea;

	//WriteProtoMesg_AS(fcgb,pmesg);
      }
#endif
      nuom++;
      break;
    default:
      {
	fprintf(stderr,"Unexpected message type %d\n",imesgtype);
	assert(FALSE);
      }
      break;
    }
  }
  fprintf(stderr, "Input Done\n");
  *Pnadt = nadt;
  *Pnidt = nidt;
  *Pnium = nium;
  *Pnimp = nimp;
  *Pnuom = nuom;
  *Pnuom_dovetail = nuom_dovetail;
  *Pnuom_containment = nuom_containment;

  if(analysis_flag) {
    fprintf(stderr,"\n\nHistogram of the UOM types\n");
    print_histogram(stderr, uom_types_histogram, 0, 1);
  }
  if(NULL != uom_types_histogram) {
    free_histogram(uom_types_histogram);
  }
}

int main(int argc, char * argv [])
{
  FILE       *fcgi = NULL;
  //FILE       *fcgb = NULL;
  int        nadt = 0;
  int        nidt = 0; 
  int        nium = 0;
  IntFragment_ID nimp = 0;
  IntEdge_ID nuom = 0;
  IntEdge_ID nuom_dovetail = 0;
  IntEdge_ID nuom_containment = 0;
  IntChunk_ID min_unitig_cid = 0;
  IntChunk_ID max_unitig_cid = 0;
  IntFragment_ID min_frag_iid = 0;
  IntFragment_ID max_frag_iid = 0;
  int analysis_flag=FALSE;
  int illegal=FALSE;

  VA_TYPE(int32) * chunk_length = NULL;
  VA_TYPE(uint32) * chunk_num_frags = NULL;
  VA_TYPE(IntMultiPos) * imp_from_iid = NULL;
  VA_TYPE(IntChunk_ID) * cid_from_iid = NULL;
  
  /**************** Process Command Line Arguments *********************/
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && 
	   ((ch = getopt(argc, argv, "APu:v:")) != EOF)) {
      switch(ch) {
      case 'A':
	analysis_flag = TRUE;
	break;
      case 'u':
	max_unitig_cid = atoi(optarg);
	break;
      case 'v':
	max_frag_iid = atoi(optarg);
	break;
      case '?':
      default :
	fprintf(stderr,"Unrecognized option -%c\n",optopt);
	errflg++;
      }
    }
    
    if((illegal == 1) // || (argc - optind < 2 )
       ){
      fprintf (stderr, "USAGE: %s \n"
	       "[-u <maximum unitig CID> ]\n"
	       "[-v <maximum fragment IID> ]\n"
	       "[-e <int>] Specify the maximum number of soft errors.\n"
	       "[-A] perform a data analysis\n"
	       "[-P] Specify ASCII output.\n"
	       " cgb-file < fragment_iid_list > unitig_id_list \n",
	       argv[0]);
      exit (EXIT_FAILURE);
    }
  }

  fprintf(stderr,"cgi file is %s\n", argv[optind]);
  fcgi = fopen(argv[optind++],"r");

  //fprintf(stderr,"cgb file is %s\n", argv[optind]);
  //fcgb = fopen(argv[optind++],"w");


  chunk_length = CreateVA_int32((max_unitig_cid)+1);
  chunk_num_frags = CreateVA_uint32((max_unitig_cid)+1);
  imp_from_iid = CreateVA_IntMultiPos((max_frag_iid)+1);
  cid_from_iid = CreateVA_IntChunk_ID((max_frag_iid)+1);


  input_mesgs_pp
    (argc, argv,
     fcgi, 
     //fcgb,
     &nadt,
     &nidt, 
     &nium,
     &nimp,
     &nuom,
     &nuom_dovetail, 
     &nuom_containment,
     &min_unitig_cid,
     &max_unitig_cid,
     &min_frag_iid,
     &max_frag_iid,
     analysis_flag,
     chunk_length,
     chunk_num_frags,
     imp_from_iid,
     cid_from_iid
     );

  fprintf(stderr,"nadt=%d\n",nadt);
  fprintf(stderr,"nidt=%d\n",nidt);
  fprintf(stderr,"nium=%d\n",nium);
  fprintf(stderr,"nuom=" F_IID "\n",nuom);
  fprintf(stderr,"nimp=" F_IID "\n",nimp);

  { 
    IntFragment_ID frag_iid;
    IntChunk_ID ichunk;
    while(EOF != fscanf(stdin,F_IID "\n",&frag_iid)) {
      ichunk = *(GetVA_IntChunk_ID(cid_from_iid,frag_iid));
      fprintf(stdout,F_IID "\n",ichunk);
    }
  }
  
  if(NULL != chunk_length) DeleteVA_int32(chunk_length);
  if(NULL != chunk_num_frags) DeleteVA_uint32(chunk_num_frags);
  if(NULL != imp_from_iid) DeleteVA_IntMultiPos(imp_from_iid);
  if(NULL != cid_from_iid) DeleteVA_IntChunk_ID(cid_from_iid);

  return 0;
}
