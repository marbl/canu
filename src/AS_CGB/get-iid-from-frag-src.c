
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
= "$Id: get-iid-from-frag-src.c,v 1.6 2006-02-13 22:16:31 eliv Exp $";
/* *******************************************************************
 *
 * Usage:
 * get-iid-from-frag-src target-string-in-fragment-comment < assembler-file > list-of-iids
 * 
 * Description: Reads a Celera Assembler i/o file with IUM messages
 * and outputs fragment IIDs to stdout of fragments that have a
 * substring in the "source" comment field that match
 * target-string-in-fragment-comment.
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

//#include "AS_CGB_all.h"
typedef CDS_IID_t IntEdge_ID;
#include "AS_CGB_histo.h"

VA_DEF(IntMultiPos)
VA_DEF(IntChunk_ID)

/****************************************************************************/
#define DEBUGGING

#define INPUT_EXTENSION        ".ovl"
#define INTERMEDIATE_EXTENSION ".ilk"
#define GRAPH_STORE_EXTENSION  ".fgb"
#define ANALYSIS_EXTENSION     ".cga"

#define CMD_BUFFER_SIZE 1024

/****************************************************************************/
/* Globals */

MesgWriter WriteMesg_AS, ErrorWriter_AS;

/*************************************************************************/

static void input_mesgs_pp
(int argc, char * argv [],
 FILE       *fcgi, 
 FILE       *fcgb,
 int        *Pnadt,
 int        *Pnidt,
 IntFragment_ID *Pnfrg,
 int        *Pnrpt,
 int        *Pnilk,
 int        *Pnium,
 IntFragment_ID *Pnimp,
 IntEdge_ID *Pnuom,
 IntEdge_ID *Pnuom_dovetail, 
 IntEdge_ID *Pnuom_containment,
 IntChunk_ID *Pmin_unitig_cid,
 IntChunk_ID *Pmax_unitig_cid,
 IntFragment_ID *Pmin_frag_iid,
 IntFragment_ID *Pmax_frag_iid,
 const int        analysis_flag,
 const char       target[]
 )
{ /* It is assumed that in the overlap records that new fragments
     point to old fragments.  */
  
  int nadt=0,nidt=0,nrpt=0,nilk=0,nium=0;
  IntFragment_ID nfrg=0,nimp=0;
  IntEdge_ID nuom=0,nuom_dovetail=0,nuom_containment=0;
  GenericMesg *pmesg;
  MesgReader ReadMesg_AS = (MesgReader)InputFileType_AS(fcgi);

  IntFragment_ID lfrag = 0; // The global LID counter.
  
  VA_TYPE(int32) * chunk_length =
    CreateVA_int32((*Pmax_unitig_cid)+1);
  VA_TYPE(int32) * chunk_num_frags =
    CreateVA_int32((*Pmax_unitig_cid)+1);
  VA_TYPE(IntMultiPos) * imp_from_iid =
    CreateVA_IntMultiPos((*Pmax_frag_iid)+1);
  VA_TYPE(IntChunk_ID) * cid_from_iid =
    CreateVA_IntChunk_ID((*Pmax_frag_iid)+1);


  const int nsample=500;
  const int nbucket=500;
  Histogram_t 
    *uom_types_histogram
    = create_histogram(nsample,nbucket,TRUE,FALSE);

  while( EOF != ReadMesg_AS(fcgi, &pmesg)) {
    const MessageType imesgtype = pmesg->t;
    
    switch(imesgtype) {
    case MESG_ADT: 
      {
	AuditMesg  *adt_mesg = (AuditMesg *)pmesg->m;
#if 0
	AuditLine  auditLine;
	AppendAuditLine_AS(adt_mesg, &auditLine, time(NULL), "CGB_pp", 
			   CM_ID, "(empty)");
#else
	VersionStampADT(adt_mesg, argc, argv);
#endif
	//WriteMesg_AS(fcgb,pmesg);
      }
      nadt++;
      break;
    case MESG_IDT: 
      {
	/*  Distance record--skip for now */
	// InternalDistMesg  *idt_mesg = (InternalDistMesg *) pmesg->m;

	//WriteMesg_AS(fcgb,pmesg);
      }
      nidt++;
      break;
    case MESG_ILK: 
      {
	//WriteMesg_AS(fcgb,pmesg);
      }
      nilk++;
      break;
    case MESG_IUM: 
      {
	const IntUnitigMesg * const ium_mesg = (const IntUnitigMesg * const) pmesg->m;
	// Do not change the value or address of the ium_mesg.
	const IntFragment_ID  num_frags = ium_mesg->num_frags;
	const IntMultiPos * const f_list = ium_mesg->f_list;
	IntFragment_ID ifrag;

	for( ifrag=0; ifrag<num_frags; ifrag++) {
	  const IntFragment_ID iid = f_list[ifrag].ident;
	  int present = f_list[ifrag].sourceInt;
	  lfrag++;
	  if(present != -1) {
	    fprintf(stdout,F_IID "\n",iid);
	  }
	}
	nimp += num_frags;
	
	//WriteMesg_AS(fcgb,pmesg);
	// pass through the Unitig message
	
	/*
	fprintf(stderr,"nium,num_frags,nimp,lfrag = %d," F_IID "," F_IID "," F_IID "\n",
                nium, num_frags, nimp, lfrag);
        */
      }
      nium ++;
      break;
    case MESG_IFG:
    case MESG_SFG:
    case MESG_OFG:
      //case MESG_FRG:
      {
	const ScreenedFragMesg * const frg_msg = (const ScreenedFragMesg * const) pmesg->m;
	const char * source = frg_msg->source;
	const char * present = strstr(source,target);
	//const Fragment_ID  eaccession = frg_msg->eaccession;
	const IntFragment_ID     iaccession = frg_msg->iaccession;
	if(present != NULL) {
	  fprintf(stdout,F_IID "\n",iaccession);
	}
      }
      nfrg ++;
      break;
    case MESG_RPT: 
      {
	//WriteMesg_AS(fcgb,pmesg);
      }
      nrpt++;
      break;
    case MESG_UOM: 
      {
      }
      nuom++;
      break;
    default:
      {
	//fprintf(stderr,"Unexpected message type %d\n",imesgtype);
	//assert(FALSE);
      }
      break;
    }
  }
  fprintf(stderr, "Input Done\n");
  *Pnadt = nadt;
  *Pnidt = nidt;
  *Pnfrg = nfrg;
  *Pnilk = nilk;
  *Pnium = nium;
  *Pnimp = nimp;
  *Pnuom = nuom;
  *Pnuom_dovetail = nuom_dovetail;
  *Pnuom_containment = nuom_containment;
  *Pnrpt = nrpt;

  if(NULL != chunk_length) DeleteVA_int32(chunk_length);
  if(NULL != chunk_num_frags) DeleteVA_int32(chunk_num_frags);
  if(NULL != imp_from_iid) DeleteVA_IntMultiPos(imp_from_iid);
  if(NULL != cid_from_iid) DeleteVA_IntChunk_ID(cid_from_iid);

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

  FILE       *fcgi = stdin;
  FILE       *fcgb = NULL;
  int        nadt = 0;
  int        nidt = 0; 
  IntFragment_ID nfrg = 0; 
  int        nrpt = 0;
  int        nilk = 0;
  int        nium = 0;
  IntFragment_ID nimp = 0;
  IntEdge_ID nuom = 0;
  IntEdge_ID nuom_dovetail = 0;
  IntEdge_ID nuom_containment = 0;
  IntChunk_ID min_unitig_cid = 0;
  IntChunk_ID max_unitig_cid = 0;
  IntFragment_ID min_frag_iid = 0;
  IntFragment_ID max_frag_iid = 0;
  int analysis_flag = FALSE;
  int illegal = FALSE;
  char *target = NULL;

  //VersionStamp(argc,argv);
  WriteMesg_AS = (MesgWriter)OutputFileType_AS(AS_BINARY_OUTPUT);
  ErrorWriter_AS = (MesgWriter)OutputFileType_AS(AS_PROTO_OUTPUT);
  
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
      case 'P':
	WriteMesg_AS = (MesgWriter)OutputFileType_AS(AS_PROTO_OUTPUT);
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
	       " < cgi-file > cgb-file \n",
	       argv[0]);
      exit (EXIT_FAILURE);
    }
  }

  if(optind == argc) {
    fprintf(stderr,"The target string is absent. It needs to be set. Try again.\n");
    exit(1);
  }
  target = argv[optind++];
  
  if(target == NULL) {
    fprintf(stderr,"The target string is NULL. It needs to be set. Try again.\n");
    exit(1);
  }

  input_mesgs_pp
    (argc, argv,
     fcgi, 
     fcgb,
     &nadt,
     &nidt, 
     &nfrg,
     &nrpt,
     &nilk,
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
     target
     );

  fprintf(stderr,"nadt=%d\n",nadt);
  fprintf(stderr,"nidt=%d\n",nidt);
  fprintf(stderr,"nfrg=" F_IID "\n",nfrg);
  fprintf(stderr,"nrpt=%d\n",nrpt);
  fprintf(stderr,"nilk=%d\n",nilk);
  fprintf(stderr,"nium=%d\n",nium);
  fprintf(stderr,"nuom=" F_IID "\n",nuom);
  fprintf(stderr,"nimp=" F_IID "\n",nimp);

  return 0;
}
