
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
static char CM_ID[] = "$Id: localUnitigging_CGW.c,v 1.7 2005-08-24 07:47:15 brianwalenz Exp $";


/*********************************************************************
 * Module:  localUnitigging_CGW
 * Description:
 *    New repeat resolution strategy which does unitigging on all potential
 *    fill material in a region of a scaffold with hopes of getting substantial
 *    pieces placed in that region.
 *
 * 
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "FbacREZ.h"
#include "PublicAPI_CNS.h"

#define CNS_DP_THRESH 1e-6
#define CNS_DP_MINLEN 30

VA_DEF(OFGMesg)
VA_DEF(OverlapMesg)
VA_DEF(IntUnitigMesg)

#include "AS_CGB_miniunitigger.h"
#include "localUnitigging_CGW.h"
#include "fragmentPlacement.h"

void OutputMergedMetaUnitig(CDS_CID_t sid,MultiAlignT *ma);

LengthT FindGapLength( ChunkInstanceT * lchunk, ChunkInstanceT * rchunk,
                       int verbose);
// externable variables for controlling use of Local_Overlap_AS_forCNS

// [initialized value is 12 -- no more than this many segments in the chain]
extern int MaxGaps;

// [ init value is 200; this could be set to the amount you extend the clear
// range of seq b, plus 10 for good measure]
extern int MaxBegGap;

// [ init value is 200; this could be set to the amount you extend the 
// clear range of seq a, plus 10 for good measure]
extern int MaxEndGap;

// [ initial value is 1000 (should have almost no effect)
// and defines the largest gap between segments in the chain]
// Also: allowed size of gap within the alignment
// -- forcing relatively good alignments, compared to those
// allowed in bubble-smoothing where indel polymorphisms are expected
extern int MaxInteriorGap;

// boolean to cause the size of an "end gap" to
// be evaluated with regard to the clear range extension
extern int asymmetricEnds;


static ReadStructp fsread = NULL;
int totalContigsBaseChange = 0;
MiniUnitiggerObject *muo = NULL;
static VA_TYPE(IntElementPos) *ContigPositions = NULL;

VA_DEF(ChunkPlacement)
static int OverlapPieceList(VA_TYPE(ChunkPlacement) *piece_list,
                            VA_TYPE(OFGMesg) *pieces,
                            VA_TYPE(OverlapMesg) *olaps);
static int CheckMetaUnitigOrderConsistency(VA_TYPE(ChunkPlacement) *piece_list,
                                           IntUnitigMesg *ium);

   
static void dumpFragsInChunk(FILE* stream, ChunkInstanceT* chunk);

static void DumpUnitigMultiAlignInfo ( VA_TYPE(ChunkPlacement) *piece_list,
                                       OFGMesg *ofgs, IntUnitigMesg *ium)
{
  int i;
  IntMultiPos *pos=ium->f_list;
  ChunkPlacement *piece;
  ContigTIterator cis;
  ChunkInstanceT *ci;
  
  
  fprintf( stderr, "in DumpUnitigMultiAlignInfo, unitig %8" F_IIDP "\n",
           ium->iaccession);

  for ( i = 0; i < ium->num_frags; i++,pos++)
  {
        piece = GetChunkPlacement(piece_list, pos->ident-1);
	fprintf( stderr, "	in DUMAI, %c fragment %8" F_UIDP ",%8" F_IIDP ", bgn: %10" F_COORDP ", end: %10" F_COORDP ", length: %10" F_COORDP ", source: %d\n", 
                 piece->status,
                 ofgs[pos->ident-1].eaccession,
                 pos->ident, pos->position.bgn, pos->position.end,
                 abs(pos->position.bgn - pos->position.end),
                 (int) pos->source);
	InitContigTIterator(ScaffoldGraph, ofgs[pos->ident-1].eaccession,
                            TRUE, FALSE, &cis);
	while(NULL != (ci = NextContigTIterator(&cis))){
	  CDS_CID_t cid = ci->id;
          fprintf(stderr,"		CI " F_CID ", cov: %d\n",
                  cid,ci->info.CI.coverageStat);
        }
  }
}

static void DumpUnitigVAInfo ( VA_TYPE(ChunkPlacement) *piece_list,
                               VA_TYPE(OFGMesg) *pieces,
                               VA_TYPE(IntUnitigMesg) *iumVA ) {
  int n_pieces = GetNumIntUnitigMesgs(iumVA);
  IntUnitigMesg *ium;
  OFGMesg *ofgs=GetOFGMesg(pieces,0);
  int i;
  for (i=0;i<n_pieces;i++) {
     ium=GetIntUnitigMesg(iumVA,i);
     DumpUnitigMultiAlignInfo(piece_list,ofgs,ium++);
  }
} 

static int CheckMetaUnitigOrderConsistency(VA_TYPE(ChunkPlacement) *piece_list,
                                           IntUnitigMesg *ium) {
   // If order is consistent with the scaffold order, return 1, else return 0
  int i;
  int rc=1;
  int num_frags=ium->num_frags;
  int monotone_incr=-1;
  int preo=-1;
  IntMultiPos *imp=ium->f_list;
  ChunkPlacement *piece; /*was:  =GetChunkPlacement(piece_list,imp->ident-1);*/
  for(i=0;i<num_frags;i++,imp++) {
     piece = GetChunkPlacement(piece_list,imp->ident-1);
     if ( piece->status == CGW_IN_SCAFFOLD) {
        if ( preo==-1) {
          preo=imp->ident; 
        } else if ( monotone_incr<0 ) {
          if ( imp->ident>preo ) {
            monotone_incr=1;
          } else {
            monotone_incr=0;
          }
          preo=imp->ident;
        } else {
          if ( monotone_incr && imp->ident < preo ) {
            rc =0;
            break;
          } else if ( ! monotone_incr && imp->ident > preo ) {
            rc =0;
            break;
          }
          preo=imp->ident;
        }
     }
   } 
   if ( !rc ) fprintf(stderr,"Input ordering of scaffolded contigs is changed upon meta-unitigging\n");
   return rc; // no conflicts
}

int CheckMetaUnitigOrientConsistency(VA_TYPE(ChunkPlacement) *piece_list,
                                     IntUnitigMesg *ium) {

   // true if orientation is consistent with scaffolded orientation and orientation implied by edges

  int i;
  int rc=1;
  int num_frags=ium->num_frags;
  NodeOrient metaOrient = A_B;
  IntMultiPos *imp=ium->f_list;
  ChunkPlacement *piece = GetChunkPlacement(piece_list,imp->ident-1);
  NodeOrient orientInMeta = ( imp->position.bgn<imp->position.end) ? A_B :  B_A;
  NodeOrient orientOfMeta = ( orientInMeta == piece->orient ) ? A_B : B_A;
  imp++;

  for(i=1;i<num_frags;i++,imp++) {
    orientInMeta =  ( imp->position.bgn<imp->position.end) ? A_B :  B_A;
    if(orientOfMeta == A_B){
      if( orientInMeta != piece->orient){
	rc = 0;
	break;
      }
    }
  }

  if ( !rc ) fprintf(stderr,"Input orientation of pieces is changed by meta-unitigging\n");

  return rc;

}

int CheckMetaUnitigPotential(VA_TYPE(ChunkPlacement) *piece_list,
                             IntUnitigMesg *ium) {
  if ( ium->num_frags < 2 ) return 0;
   // If orientation is consistent with the contigs' orientation in the scaffold, return 1, else return 0
  return 1;
}

// check whether the meta unitig contains any relevant contigs from original scaffold
// between anchor and max(anchor, rcontig-1) inclusive, reject
int MetaUnitigIncludesPiecesFromScaffold(VA_TYPE(ChunkPlacement) *piece_list,
					      IntUnitigMesg *ium,
					      int32 endCtgOrder){
  int i;
  int rc=0;
  int num_frags=ium->num_frags;
  IntMultiPos *imp=ium->f_list;
  ChunkPlacement *piece;
  for(i=0;i<num_frags;i++,imp++) {
     piece = GetChunkPlacement(piece_list,imp->ident-1);
     if(piece->order>-1 && piece->order <= endCtgOrder){
       rc=1;
       break;
     }
  }	
  return rc;
}


int MergeMetaUnitigIntoContig(VA_TYPE(ChunkPlacement) *piece_list,
                              IntUnitigMesg *ium) {
  // First of all, for each unitig within the meta-unitig, decide whether
  // it needs to be cloned as a surrogate, of it can go in as-is...
  // Single fragment unitigs can go in as is. Otherwise, previously unplaced
  // unitigs should be cloned and only those fragments which linked it into present
  // context should be instantiated.
  IntMultiPos *imp=ium->f_list;
  MultiAlignT  *newMultiAlign;
  IntMultiPos contigPos;
  ChunkPlacement *piece;
  ChunkPlacement *contained;
  // ContigT *scontig; /* not used? */
  int i;
  CDS_CID_t sid=NULLINDEX;

  // setup for contig merge
  contigPos.delta_length=0;
  contigPos.delta=NULL;
  if ( ContigPositions == NULL ) ContigPositions = CreateVA_IntMultiPos(20);
  ResetVA_IntMultiPos( ContigPositions );
  // start a list of the contigs which will be merged...
  //  for some of the pieces, we may want to instantiate new ci's
 
   for (i=0;i<ium->num_frags;i++,imp++ ) {
      piece = GetChunkPlacement(piece_list,imp->ident-1);

      contigPos.ident = piece->ident;
      contigPos.type = imp->type;
      if (imp->contained) {
          contained = GetChunkPlacement(piece_list,imp->contained-1);
          contigPos.contained = contained->ident; 
      } else {
          contigPos.contained = 0; 
      }
      contigPos.position.bgn = imp->position.bgn;
      contigPos.position.end = imp->position.end;
      AppendIntMultiPos(ContigPositions, &contigPos);
      if ( piece->status == CGW_PLACED ) {
        // consider here whether to instantiate one or more fragments
        fprintf(stderr,"Considering fragment from placed piece " F_CID "\n",piece->ident);
      } else if (piece->status == CGW_IN_SCAFFOLD) {
	fprintf(stderr,"Considering in-scaffold piece " F_CID " -- keep all fragments\n",
		piece->ident);
	sid = piece->scaff_id;
      } else {
	fprintf(stderr,"Considering fragment from unplaced piece " F_CID "\n",
		piece->ident);
      }
    }
    // here, we'll call MergeMultiAligns on the meta-unitig.
   newMultiAlign = MergeMultiAligns(ScaffoldGraph->sequenceDB,
                                    ScaffoldGraph->fragStore,
                                    ContigPositions, FALSE, TRUE,
                                    GlobalData->aligner,
                                    NULL);
   fprintf(stderr," Returned from call to MergeMultiAlign\n");
   PrintMultiAlignT(stderr,newMultiAlign,ScaffoldGraph->fragStore,0,0,0,0,0);


   OutputMergedMetaUnitig(sid,newMultiAlign);
   // then 
   return 1;
}

int  UnitigPieceList(VA_TYPE(OFGMesg) *pieces,VA_TYPE(OverlapMesg) *olaps,
                     VA_TYPE(IntUnitigMesg) *iums, VA_TYPE(IntMultiPos) *imps){


// these aren't used at present, but they need to be around for the interface
  static VA_TYPE(char)        * the_ofg_source=NULL;
  static VA_TYPE(char)        * the_ovl_source=NULL;
  static VA_TYPE(char)        * the_imp_source=NULL;
  static VA_TYPE(char)        * the_ium_source=NULL;

  muo = createMiniUnitigger(5000,5000,0);

  if ( the_ofg_source == NULL ) {
    // create the VA's for (nonexistent) char data
    the_ofg_source = CreateVA_char(0);
    the_ovl_source = CreateVA_char(0);
    the_imp_source = CreateVA_char(0);
    the_ium_source = CreateVA_char(0);
  }

  set_cgb_unique_cutoff_MiniUnitigger( muo, 5.0);
  set_overlap_error_threshold_MiniUnitigger( muo, 0.10);
  set_as_cgb_max_frag_iid_MiniUnitigger( muo, GetNumOFGMesgs(pieces) + 1);

  {
    RunMiniUnitiggerParams params =
    { pieces,
      olaps,
      the_ofg_source,
      the_ovl_source };
    RunMiniUnitiggerResults results =
    { imps,
      iums,
      the_imp_source,
      the_ium_source };

    run_MiniUnitigger( muo, &params, &results);
  }
  destroyMiniUnitigger(muo);
  return TRUE;
}

FILE *overlapfp;

int main( int argc, char *argv[])
{
  Global_CGW *data;
  char *outputPath = NULL;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int setSingleSid = FALSE;
  CDS_CID_t singleSid=0;
  int ckptNum = NULLINDEX;
  CDS_CID_t sid;
  int startingGap = 0, setStartingGap = FALSE;
  time_t t1;
  int *originalGaps;
  int *closedGap, *closedGapDelta, *lcontigIdGap, *rcontigIdGap, *lcontigLength, *rcontigLength;
  int *contigValid, *allContigLengths;
  int *closedGapAhang, *closedGapOlapLength, *closedGapBhang;
  int *closedGapLcontigBasesIntact, *closedGapRcontigBasesIntact;
  int overlapFound = 0;
  VA_TYPE(ChunkPlacement) *piece_list;
  VA_TYPE(OFGMesg) *pieces;
  VA_TYPE(OverlapMesg) *olaps;
  VA_TYPE(IntMultiPos) *imps;
  VA_TYPE(IntUnitigPos) *iums;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->stderro = stderr;
  data->stderrfp = fopen("localUnitigging.stderr","w");
  data->timefp = stderr;
  data->logfp = stderr;
  data->writer =  OutputFileType_AS(AS_PROTO_OUTPUT);
  sprintf(data->Output_File_Name,"localUnitigging.ICMs",outputPath);
  data->outfp = File_Open (data->Output_File_Name, "w", TRUE);     // cgw file

  overlapfp = fopen("localUnitigging.overlaps","w");
  
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "c:C:f:g:m:n:s:")) != EOF)){
      switch(ch) {
		case 'c':
		{
		  strcpy( data->File_Name_Prefix, argv[optind - 1]);
		  setPrefixName = TRUE;		  
		}
		break;
		case 'C':
		  startingGap = atoi(argv[optind - 1]);
		  setStartingGap = TRUE;
		  break;
		case 'f':
		{
		  strcpy( data->Frag_Store_Name, argv[optind - 1]);
		  setFragStore = TRUE;
		}
		break;
		case 'g':
		{
		  strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
		  setGatekeeperStore = TRUE;
		}
		break;	  
		case 'm':
		  MaxInteriorGap = atoi(argv[optind - 1]);
		  fprintf( stderr, "setting MaxInteriorGap to %d\n", MaxInteriorGap);
		  break;
		case 'n':
		  ckptNum = atoi(argv[optind - 1]);
		  break;
		case 's':
		  singleSid = atoi(argv[optind - 1]);
		  setSingleSid = TRUE;
		  fprintf( stderr, "setting singleSid to %d\n", singleSid);
		  break;
		case '?':
		  fprintf(stderr,"Unrecognized option -%c",optopt);
		default :
		  errflg++;
      }
    }
    if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0))
      {
	fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d outputPath = %s\n",
		argc, optind, setFragStore,setGatekeeperStore, outputPath);
	fprintf (stderr, "USAGE:  loadcgw -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum>\n");
	exit (EXIT_FAILURE);
      }
  }
  piece_list = CreateVA_ChunkPlacement(5000);
  pieces = CreateVA_OFGMesg(5000);
  olaps = CreateVA_OverlapMesg(5000);
  imps = CreateVA_IntMultiPos(5000);
  iums = CreateVA_IntUnitigMesg(5000);

  if (setStartingGap == TRUE)
	fprintf( stderr, "set starting gap to %d\n", startingGap);
  
  // hack
  {
	char temp_buf[1024];
	sprintf( temp_buf, "cp %s/db.frg.orig %s/db.frg", GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name);
	system( temp_buf );
	
	// system( "cp chrom21_dup.frgStore/db.frg.orig chrom21_dup.frgStore/db.frg" );
  }

  t1 = time(0);
  fprintf( stderr, "====> Starting at %s\n", ctime(&t1));

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, TRUE);
  GlobalData->aligner=Local_Overlap_AS_forCNS;

  // localeCam();

  originalGaps = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  closedGap = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  closedGapDelta = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  lcontigIdGap = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  rcontigIdGap = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  lcontigLength = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  rcontigLength = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  contigValid = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  allContigLengths = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  closedGapAhang = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  closedGapOlapLength = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  closedGapBhang = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  closedGapLcontigBasesIntact = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  closedGapRcontigBasesIntact = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));


  //
  // scan all the scaffolds
  //
  if (setSingleSid) { 
      // test the input sid and ensure that it is within range
      if ( singleSid >= GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) ) {
          fprintf(stderr," Specified single SID %d is out of range of the scaffold graph.\n",singleSid);
          assert(0);
      }
  }


  for (sid = singleSid; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++)
  {
        int scaff_order=0;
        CDS_CID_t anchor=0;
        ChunkPlacement cip;
#if 0
	NodeOrient lcontigOrientation;
#endif
        NodeOrient rcontigOrientation;
        int potential_fill_count=0;
        CDS_COORD_t potential_fill_length=0;
	CIScaffoldTIterator CIsTemp;
	CIScaffoldT * scaff;
	int icnt, lextension, rextension;
	CDS_CID_t rcontigID;
	ContigT *lcontig, *rcontig;
	
	lextension = rextension = 0;
        if ( setSingleSid && sid != singleSid ) break;
	
	scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
	// make sure the scaffold is there
	assert(scaff != NULL);
    
	// not interested in dead scaffold, not real scaffolds, or singleton scaffolds
    
	if ((isDeadCIScaffoldT(scaff)) ||
		(scaff->type != REAL_SCAFFOLD))
	{
	  continue;
	}
	fprintf(stderr,"\n=====================================================================\n");
	fprintf(stderr,  "=== examining scaffold " F_CID ", size %f\n", sid, scaff->bpLength.mean);
	  
	fprintf( stderr, "before unitigging in scaffold " F_CID "\n",
                 scaff->id);
	InitCIScaffoldTIterator( ScaffoldGraph, scaff, TRUE, FALSE, &CIsTemp);
	while ( NextCIScaffoldTIterator( &CIsTemp ) && 
			CIsTemp.next != GetNumGraphNodes( ScaffoldGraph->ContigGraph ) - 1)
	{
	  fprintf( stderr, "prev, curr, next: " F_CID "," F_CID "," F_CID "\n",
                   CIsTemp.prev, CIsTemp.curr, CIsTemp.next);
	}

	icnt = 0;
        ResetVA_ChunkPlacement(piece_list);
        ResetVA_OFGMesg(pieces);
        ResetVA_OverlapMesg(olaps);
        ResetVA_IntMultiPos(imps);
        ResetVA_IntUnitigMesg(iums);

	lcontig = GetGraphNode( ScaffoldGraph->ContigGraph, scaff->info.Scaffold.AEndCI);
        anchor = lcontig->id;
#if 0
        cip.ident = anchor;
        cip.order = scaff_order++;
        cip.status = CGW_IN_SCAFFOLD;
        cip.scaff_id = sid;
        cip.scaff_bpLength = (CDS_COORD_t) scaff->bpLength.mean;
        cip.scaff_numCI = scaff->info.Contig.numCI;
        cip.orient = lcontigOrientation;
        cip.AEndCI = lcontig->info.Contig.AEndCI;
        cip.BEndCI = lcontig->info.Contig.BEndCI;
        cip.edge = NULL;
        AppendVA_ChunkPlacement(piece_list,&cip);
#endif


        potential_fill_count=0;
        potential_fill_length=0;
	//rcontigID = lcontig->BEndNext;
        rcontigID = anchor;
	while ( rcontigID != -1 )
	{
          GraphEdgeIterator  ci_edges;
          CIEdgeT  * edge;

	  rcontig = GetGraphNode( ScaffoldGraph->ContigGraph, rcontigID);
	  if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
		rcontigOrientation = A_B;
          else
		rcontigOrientation = B_A;


            cip.ident = rcontig->id;
            cip.order = scaff_order++;
            cip.status = CGW_IN_SCAFFOLD;
            cip.scaff_id = sid;
            cip.scaff_numCI = scaff->info.Contig.numCI;
            cip.scaff_bpLength = (CDS_COORD_t) scaff->bpLength.mean;
            cip.orient = rcontigOrientation;
            cip.AEndCI = rcontig->info.Contig.AEndCI;
            cip.BEndCI = rcontig->info.Contig.BEndCI;
            cip.edge = NULL;
	    cip.type = AS_CONTIG;
            AppendVA_ChunkPlacement(piece_list,&cip);

	  
          
           // collect edges from chunk to other chunks

           InitGraphEdgeIterator (ScaffoldGraph->RezGraph, rcontig -> id, ALL_END,
                                  ALL_EDGES, GRAPH_EDGE_DEFAULT,
                                  & ci_edges);
           while  ((edge = NextGraphEdgeIterator (& ci_edges)) != NULL)
             {
              ChunkInstanceT  * other_chunk;

              if  (edge -> idA == rcontig->id)
                  other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                   edge -> idB);
                else
                  other_chunk = GetGraphNode(ScaffoldGraph->RezGraph,
                                                   edge -> idA);

            {
	       CIScaffoldT * oscaff;
               EdgeCGW_T *rawEdge;
               cip.ident = other_chunk->id;
               cip.order = -1;
               if ( other_chunk->scaffoldID != -1 ) {
                   oscaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, other_chunk->scaffoldID);
                   if ( other_chunk->scaffoldID == sid ) {
                     cip.status = CGW_IN_SCAFFOLD;
                   } else {
                     cip.status = CGW_PLACED;

		     fprintf(stderr,"Chunk " F_CID " in region is CGW_PLACED status; type is %d; fragments are as follows:\n",cip.ident,other_chunk->type);
		     dumpFragsInChunk(stderr,other_chunk);
		     
                   }
                   cip.scaff_id=other_chunk->scaffoldID;
                   cip.scaff_bpLength= (CDS_COORD_t) oscaff->bpLength.mean;
                   cip.scaff_numCI=oscaff->info.Contig.numCI;
               } else {
                   cip.status = CGW_UNPLACED;
                   cip.scaff_id=-1;
                   cip.scaff_bpLength=0;
                   cip.scaff_numCI=0;

		   fprintf(stderr,"Chunk " F_CID " in region is CGW_UNPLACED status; type is %d; fragments are as follows:\n",cip.ident,other_chunk->type);
		   dumpFragsInChunk(stderr,other_chunk);

               }
               cip.AEndCI = cip.ident;
               cip.BEndCI = cip.ident;
               cip.edge = edge;
              
            if ( other_chunk->scaffoldID == sid ) {
               fprintf(stderr,"\t\t(saw link to chunk " F_CID " already in this scaffold " F_CID " from chunk " F_CID " in scaffold " F_CID ")\n", 
                       other_chunk -> id,
                       other_chunk->scaffoldID,
                       rcontig->id,
                       sid );
            } else {

#undef UNITIG_CHUNKS	       
#ifdef UNITIG_CHUNKS	       
               fprintf(stderr,"  Putting unplaced or other chunk " F_CID " onto stack (from chunk " F_CID ")\n", 
                               other_chunk -> id, rcontig->id );
               AppendVA_ChunkPlacement(piece_list,&cip);
               potential_fill_count+=1;
               potential_fill_length+= (CDS_COORD_t) other_chunk->bpLength.mean;

    #if 0
	       {
		 CIFragT *otherfrg;
		 CGWFragIterator frags;
		 fprintf(stderr,"Testing chunk frag iterator ...\n");
		 InitCIFragTInChunkIterator(&frags,other_chunk,TRUE);
		 while( NextCIFragTInChunkIterator(&frags,&otherfrg) ){
		   fprintf(stderr,
			   "\tChecking frag " F_CID " ... mate is%s in scaffold " F_CID "\n",
			   otherfrg->iid,
			   matePlacedIn(otherfrg,sid) ? "" : " not",
			   sid);
		 }
	       }
    #endif
#else
	       {
		 CIFragT *otherfrg;
		 ChunkPlacement frgcip;
		 CGWFragIterator frags;

		 frgcip.order = -1;
		 frgcip.status = cip.status;
		 frgcip.scaff_id = cip.scaff_id;

		 InitCIFragTInChunkIterator(&frags,other_chunk,TRUE);
		 while( NextCIFragTInChunkIterator(&frags,&otherfrg) ){
		   int beg,end,ori;

		   // if we get here, we should be looking at fragments
		   // in chunks outside the scaffold of interest
		   assert(scaffoldOf(otherfrg->iid) != sid);

		   // we want to restrict our attention to fragments
		   // that are mated into the scaffold of interest
		   if(! matePlacedIn(otherfrg,sid))continue;


		   beg = (int) otherfrg->contigOffset5p.mean;
		   end = (int) otherfrg->contigOffset3p.mean;
		   if(beg<end){
		     ori=1;
		   } else {
		     ori=-1;
		   }
		   frgcip.ident = otherfrg->iid;
		   frgcip.AEndCI = otherfrg->iid;
		   frgcip.BEndCI = otherfrg->iid;
		   frgcip.type = otherfrg->type;

		   frgcip.orient = cip.orient;
		   if(ori==-1){
		     switch(frgcip.orient){
		     case 'F':
		       frgcip.orient = 'R';
		       break;
		     case 'R':
		       frgcip.orient = 'F';
		       break;
		     default:
		       assert(0);
		     }
		   }

    #define NO_EDGE_ASSOCIATED_WITH_FRG
    #ifdef  NO_EDGE_ASSOCIATED_WITH_FRG
		   frgcip.edge=NULL;
    #else
		   frgcip.edge->orient = cip.edge->orient;
		   if(ori==-1){
		     switch(cip.edge->orient){
		     case AB_BA:
		     case BA_AB:
		     case XX_XX:
		       break;
		     case AB_AB:
		       frgcip.edge->orient = BA_BA; 
		       break;
		     case BA_BA:
		       frgcip.edge->orient = AB_AB;
		       break;
		     default:
		       assert(0);
		     }
		   }
    #endif


		   frgcip.scaff_numCI=0;
		   frgcip.scaff_bpLength = (ori==1 ? end-beg : beg-end);


		   fprintf(stderr,"  Putting fragment " F_CID " from unplaced or other chunk " F_CID " onto stack (based on edge from chunk " F_CID ")\n", 
			   frgcip.ident,other_chunk -> id, rcontig->id );

		   AppendVA_ChunkPlacement(piece_list,&frgcip);

		   potential_fill_count+=1;
		   potential_fill_length+= (CDS_COORD_t) frgcip.scaff_bpLength;

		 }
		 CleanupCIFragTInChunkIterator(&frags);
	       }
#endif


            }


               if ( edge->nextRawEdge > 0 ) {
                 for(rawEdge = GetGraphEdge(ScaffoldGraph->RezGraph,edge->nextRawEdge);
                   rawEdge != NULL;
                   rawEdge = GetGraphEdge(ScaffoldGraph->RezGraph,rawEdge->nextRawEdge)){
                   assert(!rawEdge->flags.bits.isDeleted);
                   PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,"*raw* ", rawEdge, rawEdge->idA );
                   if (rawEdge->distIndex > 0 ) {
                     DistT *dist=GetDistT(ScaffoldGraph->Dists,rawEdge->distIndex);
                     fprintf(stderr," based on dist information from library " F_CID ", mean %d\n",rawEdge->distIndex,
                                    (int)dist->mean);
                   }
                 }
               } else {
                   PrintGraphEdge(GlobalData->stderrc,ScaffoldGraph->RezGraph,"*raw* ", edge, edge->idA );
                   if (edge->distIndex > 0 ) {
                     DistT *dist=GetDistT(ScaffoldGraph->Dists,edge->distIndex);
                     fprintf(stderr," based on dist information from library " F_CID ", mean %d\n",edge->distIndex,
                                    (int)dist->mean);
                   }
               }
            }
               
             }

          if ( ( anchor != rcontigID && rcontig->bpLength.mean > 20000 ) || rcontig->BEndNext == -1 ) {
             fprintf(stderr, "\n\tFound a region from anchor contig " F_CID " to " F_CID "\n",anchor,rcontig->id);
             fprintf(stderr,"\t Should try local unitigging here\n");
	     if(GetNumChunkPlacements(piece_list)>1){
	       overlapFound= OverlapPieceList(piece_list,pieces,olaps);
	       fprintf(stderr,"\t OverlapPieceList returned %d overlaps\n",overlapFound);
	       if ( overlapFound > 0 ) {
		 UnitigPieceList(pieces,olaps,iums,imps);
		 DumpUnitigVAInfo ( piece_list, pieces, iums );
		 { size_t i;
                 IntUnitigMesg *ium; 
                 for (i=0;i<GetNumIntUnitigMesgs(iums);i++) {
                   ium=GetIntUnitigMesg(iums,i);

		   // if only one element, no merging done
                   if ( !CheckMetaUnitigPotential( piece_list,ium)) continue;

		   // if the unitig changes the order of contigs in scaffold, reject
                   if ( !CheckMetaUnitigOrderConsistency( piece_list,ium)){
#ifdef NO_CONTIG_REORDERING		   
		     continue; 
#else
		     fprintf(stderr,"DIRE WARNING: Contigs reordered -- did you want to allow that?\n");
#endif
		   }

		   // the overlap checker should have restricted orientations so that only allowed ones
		   // occur, but just in case ...
                   assert( CheckMetaUnitigOrientConsistency( piece_list,ium));

		   // if the unitig doesn't contain any contigs from original scaffold
		   // between anchor and and rcontig, reject
		   // Also, if it only contains rcontig and rcontig != anchor, then reject.
		   // This latter case should be rejected because we are going to reprocess
		   // this contig next time through the loop ...
		   if ( ! MetaUnitigIncludesPiecesFromScaffold(piece_list,ium,
							       scaff_order - (anchor==rcontigID ? 1 : 2))) {

		     if ( anchor!=rcontigID &&
			  MetaUnitigIncludesPiecesFromScaffold(piece_list,ium,
							       scaff_order - 1)){
		       fprintf(stderr,"DIRE WARNING: meta-unitig affects only last contig in region, which will be reprocessed next time through looping on anchor!\n");
		     }


		     continue;

		   }

		   // apply the unitig
                   MergeMetaUnitigIntoContig(piece_list,ium);
                 }
		 }
	       }
	     }
             anchor=rcontig->id;
             scaff_order=0;
             ResetVA_ChunkPlacement(piece_list);
             ResetVA_OFGMesg(pieces);
             ResetVA_OverlapMesg(olaps);
             ResetVA_IntMultiPos(imps);
             ResetVA_IntUnitigMesg(iums);
          }
	  rcontigID = rcontig->BEndNext;
	}
  }

  t1 = time(0);
  fprintf( stderr, "====> Ending at %s\n", ctime(&t1));

  exit(0);
}

int compCP( const void *i1, const void *i2)
{
  ChunkPlacement *c1=(ChunkPlacement *) i1;
  ChunkPlacement *c2=(ChunkPlacement *) i2;
  if ( c1->ident >  c2->ident )
        return 1;
  else if ( c1->ident < c2->ident )
        return -1;
  else if ( c1->type != c2->type )
        return c1->type - c2->type;
  else 
        return 0;
}


static void dumpFastaRecord(FILE *stream, char *header, char *sequence);

static int OverlapPieceList(VA_TYPE(ChunkPlacement) *piece_list, VA_TYPE(OFGMesg) *pieces,VA_TYPE(OverlapMesg) *olaps ) {
   static VA_TYPE(char) *consensusA=NULL, *qualityA=NULL;
   static VA_TYPE(char) *consensusB=NULL, *qualityB=NULL;
   MultiAlignT *ma;
   MultiAlignT **multialigns;
   char **seqs;
   size_t *seqlen;
   OverlapMesg *O;
   ChunkPlacement *cip;
   ChunkPlacement *left_ci,*right_ci;
   // first of all, sort the list and remove duplicates
   ChunkPlacement *unsorted_pieces;
   size_t n_pieces = GetNumints(piece_list);
   size_t i;
   CDS_CID_t this_chunk,prev_chunk;
   int32 this_ori,prev_ori;
   FragType this_type, prev_type;

   if(n_pieces==1)return 0;
   if (consensusA== NULL ) 
   {
    consensusA = CreateVA_char(0);
    qualityA = CreateVA_char(0);
    consensusB = CreateVA_char(0);
    qualityB = CreateVA_char(0);
   } 
   else 
   {
    ResetVA_char(consensusA);
    ResetVA_char(qualityA);
    ResetVA_char(consensusB);
    ResetVA_char(qualityB);
   }
   if ( n_pieces > 0 ) {
      ContigTIterator cis;
      ChunkInstanceT *ci;
      InternalFragMesg  A, B;
      MultiAlignT *afrag,*bfrag;
      CIStatus cistat;
      OFGMesg ofg;
      char *aseq,*bseq;
      size_t length=0;
      time_t currentTime = time(0);
      CDS_CID_t iuniq=0;
      int is_placed=0;
      ofg.sequence = NULL;
      ofg.quality = NULL;
      ofg.source = NULL;
      ofg.action = AS_ADD;
      ofg.elocale = 0;
      ofg.ebactig_id = 0;
      ofg.ibactig_id = 0;
      ofg.entry_time = currentTime;
      A.quality=NULL;
      B.quality=NULL;
      unsorted_pieces = GetChunkPlacement(piece_list,0);
      qsort( unsorted_pieces, n_pieces, sizeof( ChunkPlacement ), &compCP );
      multialigns = (MultiAlignT **) safe_malloc(n_pieces*sizeof(MultiAlignT *));
      seqs = (char **) safe_malloc(n_pieces*sizeof(char *));
      seqlen = (size_t *) safe_malloc(n_pieces*sizeof(size_t));
      prev_chunk=0;
      this_chunk=0;
      for (i=0;i<n_pieces;i++) {
        cip = GetChunkPlacement(piece_list,i);
        this_chunk=cip->ident;
        cistat=cip->status;
	this_type = cip->type;
	this_ori = cip->orient;
	//process, avoiding duplicates
        if ( this_chunk != prev_chunk || 
	     this_type != prev_type
#ifdef ALLOW_BOTH_ORI	  
	  || this_ori != prev_ori
#endif
	     ){ 
          char stat;


	  // get the sequence -- different ways for contigs and reads, of course
	  if(this_type == AS_CONTIG){
	    GetConsensus(ScaffoldGraph->RezGraph, this_chunk, consensusA, qualityA);
	    seqlen[iuniq] = GetNumchars(consensusA);
	    seqs[iuniq] = (char *) safe_malloc((seqlen[iuniq]+1)*sizeof(char));
	    strcpy(seqs[iuniq],Getchar(consensusA,0));
	    ma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, this_chunk, FALSE);
	  } else {
	    static ReadStructp fsread = NULL;
	    static char frgSeqBuffer[AS_BACTIG_MAX_LEN+1], frgQltbuffer[AS_BACTIG_MAX_LEN+1];
	    uint clr_bgn, clr_end;

	    assert(this_type = AS_READ);

	    if ( fsread == NULL ) fsread = new_ReadStruct();

	    // get the read and its sequence
	    getFragStore( ScaffoldGraph->fragStore, this_chunk, FRAG_S_ALL, fsread);
	    getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	    getSequence_ReadStruct( fsread, frgSeqBuffer, frgQltbuffer, AS_BACTIG_MAX_LEN);
	    frgSeqBuffer[clr_end]='\0';

	    // copy its sequence
	    seqlen[iuniq] = clr_end - clr_bgn;
	    seqs[iuniq] = (char *) safe_malloc((seqlen[iuniq]+1)*sizeof(char));
	    strcpy(seqs[iuniq],frgSeqBuffer+clr_bgn);
	    ma = NULL;
	  }
	  length+=seqlen[iuniq];
	  ofg.clear_rng.bgn = 0;
	  ofg.clear_rng.end = seqlen[iuniq];



          ofg.iaccession = iuniq+1;
          ofg.eaccession = this_chunk;
	  ofg.type = cip->type;
          AppendVA_OFGMesg(pieces,&ofg);
          fprintf(stderr,"Chunk " F_UID " ",ofg.eaccession);
          if ( cistat == CGW_IN_SCAFFOLD ) {
            fprintf(stderr," * "); 
            is_placed++;
          }
          if ( cistat == CGW_IN_SCAFFOLD ) {
            stat = 'S';
          } else if ( cistat == CGW_PLACED ) {
            stat = 'P';
          } else {
            stat = 'U';
          }
          fprintf(stderr," (local index " F_IID ")\n",ofg.iaccession); 
          fprintf(stderr,"GenGraphNODE: " F_UID " %c " F_CID " " F_CID " " F_COORD " " F_COORD " " F_SIZE_T " " F_SIZE_T "\n",
                   ofg.eaccession, stat, cip->scaff_id, cip->scaff_numCI, cip->scaff_bpLength,
                   ofg.clear_rng.end,ma!=NULL ? GetNumIntMultiPoss(ma->f_list) : 1,
                   ma!=NULL ? GetNumIntUnitigPoss(ma->u_list) : 0);
	  InitContigTIterator(ScaffoldGraph, ofg.eaccession, TRUE, FALSE, &cis);
	  while(NULL != (ci = NextContigTIterator(&cis))){
	    CDS_CID_t cid = ci->id;
            fprintf(stderr,"	CI " F_CID ", cov: %d\n",cid,ci->info.CI.coverageStat);
          }
          SetVA_ChunkPlacement(piece_list,iuniq,cip);
          multialigns[iuniq++] = ma; 
        }
        prev_chunk = this_chunk;
	prev_type  = this_type;
	prev_ori = this_ori;
      }
      ResetToRange_ChunkPlacement(piece_list,iuniq);
      fprintf(stderr,"\n" F_CID " pieces to unitig, combined length of " F_SIZE_T "\n",iuniq,length);
      for (i=0;i<iuniq-1;i++) {
         // do overlaps here
        int j;
        int i_placed=0,j_placed;
        int i_orient,where;
        left_ci=NULL;
        right_ci=NULL;
        afrag=multialigns[i];
        A.sequence = seqs[i];
        A.iaccession = i+1;
        A.eaccession = i+1;
        aseq=seqs[i];
        cip = GetChunkPlacement(piece_list,i);
	i_orient = cip->orient;
        if ( cip->status == CGW_IN_SCAFFOLD ) {
  	    assert(cip->order>=0);
	    i_placed=1;
            if ( cip->order == 0 ) {
               left_ci = cip;
            } else if (cip->order == is_placed-1) {
               right_ci = cip;
            }
        }
        for (j=i+1;j<iuniq;j++) {
            j_placed=0;
            bfrag = multialigns[j];
            bseq = seqs[j];
            B.sequence = seqs[j];
            B.iaccession = j+1;
            B.eaccession = j+1;
            cip = GetChunkPlacement(piece_list,j);
            if ( cip->status == CGW_IN_SCAFFOLD ) {
              j_placed=1;
              if ( cip->order == 0 ) {
                 left_ci = cip;
              } else if (cip->order == is_placed-1) {
                 right_ci = cip;
              }
            }
            if ( ! (j_placed && i_placed ) ) {
	      //           for (ori=0;ori<2;ori++){ 
	      {
              int ori = (cip->orient == i_orient ? 0 : 1);
              int band_bgn=-seqlen[j];
              int band_end=seqlen[i];
              if ( left_ci != NULL ) {
                 fprintf(stderr,"Examining potential overlap at left end of region (ci = " F_CID ")\n",left_ci->ident);
                 if (  seqlen[i] > 10000) {
                    band_bgn=band_end-10000;
                 } 
              }
              if ( right_ci != NULL ) {
                 fprintf(stderr,"Examining potential overlap at right end of region (ci = " F_CID ")\n",right_ci->ident);
                 if ( seqlen[i] > 10000 ) {
                    band_end=10000;
                 }
              }
	      O = Local_Overlap_AS(&A,&B,band_bgn,band_end, ori, .10, 1e-6, 30, AS_FIND_LOCAL_ALIGN, &where);
	      //O = Local_Overlap_AS_forCNS(aseq,bseq,-seqlen[j], seqlen[i],0, .10, 1e-6, 30, AS_FIND_LOCAL_ALIGN);
              if ( O != NULL ) {
               fprintf(overlapfp,"Found overlap between chunks " F_CID " and " F_CID " (orient %d)\n",
		       GetChunkPlacement(piece_list,i)->ident,
		       GetChunkPlacement(piece_list,j)->ident,
		       ori);
               O->min_offset=O->max_offset=O->ahg;
               Print_Overlap_AS(overlapfp,&A,&B,O);
               //Print_Overlap(stderr,aseq,bseq,O);
               AppendVA_OverlapMesg(olaps,O);
              }
            }
          }   
        }
        free(seqs[i]);
      }  
      free(multialigns);
      free(seqs);
      free(seqlen);
   }
   return GetNumOverlapMesgs(olaps);
}
   
/***********************************************/
#define SEGLEN 100
static void dumpFastaRecord(FILE *stream, char *header, char *sequence){
  size_t i;
  size_t FragLen = strlen(sequence);

  fprintf(stream,"> %s\n", header);
  for (i = 0; i < FragLen; i += SEGLEN)
    if (FragLen-i < SEGLEN)
      fprintf(stream,"%.*s\n",FragLen-i,sequence+i);
    else
      fprintf(stream,"%.*s\n",SEGLEN,sequence+i);
}
/***********************************************/


static void PrintMateInfo(FILE* stream,CDS_CID_t mateIID){
  CIFragT *mfrag = GetCIFragT(ScaffoldGraph->CIFrags, 
			      GetInfoByIID(ScaffoldGraph->iidToFragIndex,
					   mateIID)->fragIndex);
  DistT *dist;
  AssertPtr(mfrag);
  dist = GetDistT(ScaffoldGraph->Dists, mfrag->dist);

  fprintf(stream,"\t\tMate " F_CID " dist (%f,%f) Unitig " F_CID " ChunkInstance " F_CID " Contig " F_CID " Scaffold " F_CID "\n",
	  mfrag->iid,
	  dist->mean,
	  dist->stddev,
	  mfrag->cid,
	  mfrag->CIid,
	  mfrag->contigID,
	  (mfrag->contigID != NULLINDEX ? GetGraphNode(ScaffoldGraph->ContigGraph,mfrag->contigID)->scaffoldID : -1)
	  );
}

static void dumpFragsInCI(FILE* stream, ChunkInstanceT* chunk){

  static MultiAlignT *unitig=NULL;
  IntMultiPos *f_list;
  int num_frags;
  int i,rv;
  char *prefix, pfx1[] = "Surrogate ", pfx2[] = "";

  if(chunk->flags.bits.isScaffold||chunk->flags.bits.isContig){
    fprintf(stream,"Tried to dumpFragsInCI() on a scaffold or contig -- null op\n");
    return;
  }

  if(chunk->flags.bits.isSurrogate){
    fprintf(stderr, "\t\t *** Chunk is a surrogate! *** \n");
    prefix = pfx1;
  } else {
    prefix = pfx2;
  }

  if(unitig==NULL){
    unitig= CreateEmptyMultiAlignT();
  }

  rv = ReLoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, unitig, chunk->id, TRUE);
  assert(rv==0);
  num_frags = GetNumIntMultiPoss(unitig->f_list);
  f_list = GetIntMultiPos(unitig->f_list,0);
  for(i=0;i<num_frags;i++){
    //    int frgAEnd = f_list[i].position.bgn;
    //    int frgBEnd = f_list[i].position.end;
    int frgIdent = f_list[i].ident;

    fprintf(stream,"\t%sFrag " F_CID " is externally mated as follows:\n", prefix,frgIdent);
    {
      CGWMateIterator mates;
      CDS_CID_t linkIID;
      InitCGWMateIterator(&mates,frgIdent,EXTERNAL_MATES_ONLY,chunk);
      while(NextCGWMateIterator(&mates,&linkIID)){
        PrintMateInfo(stream,linkIID);
      }
    }
  }
}

static void dumpFragsInChunk(FILE* stream, ChunkInstanceT* chunk){
  if(chunk->flags.bits.isScaffold){
    fprintf(stream,"Tried to dumpFragsInChunk() on a scaffold -- null op\n");
    return;
  }
  if(chunk->flags.bits.isContig){
    static MultiAlignT *contig=NULL;
    IntUnitigPos *u_list;
    int num_unitigs, rv, i;

    if(contig==NULL){
      contig= CreateEmptyMultiAlignT();
    }
    
    rv = ReLoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, contig, chunk->id,FALSE);
    assert(rv==0);
    num_unitigs = GetNumIntUnitigPoss(contig->u_list);
    for(i=0;i<num_unitigs;i++){
      IntUnitigPos *u_list = GetIntUnitigPos(contig->u_list,0);
      cds_int32 utgID = u_list[i].ident;
      fprintf(stream,
	      "  Chunk (contig) contains chunk (unitig) " F_CID " -- dump it\n", utgID);
      dumpFragsInCI(stream,GetGraphNode(ScaffoldGraph->CIGraph,utgID));
    }
    return;

  } else {
    dumpFragsInCI(stream,chunk);
  }
}




#if 0
void MassageAndApplyNewMultAlign(MultiAlignT *newMultiAlign, ContigT *contig){

  // Sort the Unitigs from left to right
  MakeCanonicalMultiAlignT(newMultiAlign);

  assert(newMultiAlign->refCnt == 0);

  // NOTE: we keep the new multi-align in the cache for a bit, but free it at the end of this routine
  InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, contig->id,FALSE, newMultiAlign, TRUE);

  assert(newMultiAlign->refCnt == 1);
  AddReferenceMultiAlignT(newMultiAlign); // this makes sure we survive any cache flushes

  // Propagate Overlaps, Tandem Marks and Branch Points to the new Contig
  PropagateOverlapsToNewContig(contig, ContigPositions, aEndID, aEndEnd, bEndID, bEndEnd, scaffold->id, FALSE);

  // Insert the new contig into the scaffold, in lieu of the old contigs
  ReplaceContigsInScaffolds(scaffold, contig,  ContigPositions, newOffsetAEnd, newOffsetBEnd, deltaOffsetBEnd);

  // Mark all frags as being members of this Contig, and set their offsets
  UpdateNodeFragments(ScaffoldGraph->RezGraph,contig->id, TRUE, FALSE);
  // Mark all of the Unitigs of this CI and set their offsets
  UpdateNodeUnitigs(newMultiAlign,contig);

  // Update simulator coordinates
  UpdateContigSimCoordinates(contig);
  UpdateScaffoldSimCoordinates(scaffold);

  // Create the raw link-based edges
  { 
    GraphEdgeStatT stats;
    BuildGraphEdgesFromMultiAlign(ScaffoldGraph->ContigGraph, contig, newMultiAlign, &stats, TRUE);
  }

  if(GlobalData->debugLevel > 0){
    fprintf(GlobalData->stderrc,"* Here is the new contig before merging: [%g,%g]\n",
	    newOffsetAEnd.mean, newOffsetBEnd.mean);
    DumpContig(GlobalData->stderrc,ScaffoldGraph, contig,FALSE);
  }
  // Merge the edges incident on this contig
  MergeNodeGraphEdges(ScaffoldGraph->ContigGraph, contig, FALSE, TRUE, FALSE);

  RemoveReferenceMultiAlignT(newMultiAlign); // we added a reference previously
  if(newMultiAlign->refCnt == 0){// we are the owner
    DeleteMultiAlignT(newMultiAlign);
  }else{ // the cache owns the memory
    // Free up the cache space from the new multiAlignT
    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
  }
  return TRUE;

}
#endif

void OutputMergedMetaUnitig(CDS_CID_t sid,MultiAlignT *ma){

  // this function largely copied from OutputContigsFromMultiAligns() -- should decompose the reused part!!

  GenericMesg		pmesg;
  IntConConMesg		icm_mesg;
  IntUnitigPos		*uptr;
  CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
  static CDS_CID_t  ctgID = NULLINDEX;
  int32 ubufSize = 100;
  CDS_IID_t numFrag;
  CDS_IID_t numUnitig;
  CDS_IID_t * tmpSource;
  IntMultiPos *mp;
  IntUnitigPos *up;
  int i;

  if(ctgID==NULLINDEX){
    ctgID = GetNumGraphNodes(ScaffoldGraph->ContigGraph);
  } else {
    ctgID++;
  }

  pmesg.m = &icm_mesg;
  pmesg.t = MESG_ICM;
  
  icm_mesg.unitigs = (IntUnitigPos *) safe_malloc(ubufSize*sizeof(IntUnitigPos));

  numFrag = GetNumIntMultiPoss(ma->f_list);
  mp = GetIntMultiPos(ma->f_list,0);
  numUnitig = GetNumIntUnitigPoss(ma->u_list);
  up = GetIntUnitigPos(ma->u_list,0);
      
  tmpSource = safe_malloc((GetNumIntMultiPoss(ma->f_list) + 1) * sizeof(CDS_IID_t));

  if(numUnitig >= ubufSize){
    ubufSize = numUnitig * 2;
    icm_mesg.unitigs = (IntUnitigPos *) safe_realloc(icm_mesg.unitigs, ubufSize*sizeof(IntUnitigPos));
  }
  uptr = icm_mesg.unitigs;
  for(i = 0; i < numUnitig; i++){
    IntUnitigPos *iup = up + i;
    NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, iup->ident);
    if(unitig->type == DISCRIMINATORUNIQUECHUNK_CGW){
      uptr[i].type = AS_UNIQUE_UNITIG;
    }else{
      if(unitig->scaffoldID != NULLINDEX){
	if(!unitig->flags.bits.isSurrogate){
	  uptr[i].type = AS_ROCK_UNITIG;
	}else  if(unitig->flags.bits.isStoneSurrogate){
	  uptr[i].type = AS_STONE_UNITIG;
	}else{
	  uptr[i].type = AS_PEBBLE_UNITIG;
	}
      }else{
	uptr[i].type = AS_SINGLE_UNITIG;
      }
    }
    uptr[i].position = iup->position;
    uptr[i].delta_length = iup->delta_length;
    uptr[i].delta = iup->delta;
    if(unitig->type == RESOLVEDREPEATCHUNK_CGW){
      iup->ident = unitig->info.CI.baseID; // map back to the parent of this instance
    }
    uptr[i].ident = iup->ident;
  }
  // Null out the source field
  for(i = 0; i < numFrag; i++){
    IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,
			       (CDS_CID_t)mp_i->source);
    tmpSource[i] = (CDS_IID_t) mp_i->source;
    mp_i->source = NULL;
  }

  icm_mesg.placed = (scaffold && (scaffold->type == REAL_SCAFFOLD)?AS_PLACED:AS_UNPLACED);
  icm_mesg.iaccession = ctgID;
  icm_mesg.forced = 0;
  icm_mesg.num_pieces = numFrag;
  icm_mesg.pieces = mp;
  icm_mesg.num_unitigs = numUnitig;
  icm_mesg.length = GetMultiAlignLength(ma);
  icm_mesg.num_vars = 0;
  icm_mesg.v_list   = NULL;
  if(icm_mesg.num_unitigs > 1){
    icm_mesg.consensus = ""; // Getchar(ma->consensus,0);
    icm_mesg.quality = ""; // Getchar(ma->quality,0);
  }else{
    icm_mesg.consensus = Getchar(ma->consensus,0);
    icm_mesg.quality = Getchar(ma->quality,0);
  }
      
  if(icm_mesg.num_unitigs > 1){
    assert(sid != NULLINDEX);
    (GlobalData->writer)(GlobalData->outfp1,&pmesg);
  }else{
    if(sid == NULLINDEX) {// contig is not placed
      GenericMesg		mesg;
      IntDegenerateScaffoldMesg dsc_mesg;
      NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph,
				       GetIntUnitigPos(ma->u_list,0)->ident);
          
      assert(unitig != NULL);
      if(unitig->info.CI.numInstances == 0){ // If this unitig has been placed as a surrogate, don't output contig
	dsc_mesg.icontig = ctgID;
	mesg.m = &dsc_mesg;
	mesg.t = MESG_IDS;
            
	(GlobalData->writer)(GlobalData->outfp,&pmesg); // write the contig
	(GlobalData->writer)(GlobalData->outfp,&mesg);  // write the associated degenerate scaffold
      }else{
	// do nothing. The unitig in this contig appears as a surrogate elsewhere in the assembly
      }
    }else{ // Contig is placed
      (GlobalData->writer)(GlobalData->outfp,&pmesg); // write the contig
    }     
  }
      
  // Restore the source values
  for(i = 0; i < numFrag; i++){
    IntMultiPos *mp_i = GetIntMultiPos(ma->f_list,i);
    mp_i->source = (char *)tmpSource[i];
  }
  free(tmpSource);

free(icm_mesg.unitigs);
fflush(NULL);

}
