
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
static char CM_ID[] = "$Id: localUnitigging_CGW.c,v 1.2 2004-09-23 20:25:19 mcschatz Exp $";


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
  ChunkPlacement *piece=GetChunkPlacement(piece_list,imp->ident-1);
  for(i=0;i<num_frags;i++,imp++) {
     piece = GetChunkPlacement(piece_list,imp->ident-1);
     if ( piece->status == CGW_PLACED) {
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
   // If orientation is consistent with the contigs' orientation in the scaffold, return 1, else return 0
  return 1;
}

int CheckMetaUnitigPotential(VA_TYPE(ChunkPlacement) *piece_list,
                             IntUnitigMesg *ium) {
  if ( ium->num_frags < 2 ) return 0;
   // If orientation is consistent with the contigs' orientation in the scaffold, return 1, else return 0
  return 1;
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
  ContigT *scontig;
  int i;

  // setup for contig merge
  contigPos.delta_length=0;
  contigPos.delta=NULL;
  if ( ContigPositions == NULL ) ContigPositions = CreateVA_IntMultiPos(20);
  ResetVA_IntMultiPos( ContigPositions );
  // start a list of the contigs which will be merged...
  //  for some of the pieces, we may want to instantiate new ci's
 
   for (i=0;i<ium->num_frags;i++,imp++ ) {
      piece = GetChunkPlacement(piece_list,imp->ident-1);
      scontig = GetGraphNode( ScaffoldGraph->ContigGraph,piece->ident);
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
      if ( piece->status != CGW_PLACED ) {
        // consider here whether to instantiate one or more fragments
        fprintf(stderr,"Considering new piece " F_CID "\n",piece->ident);
      } else {
        fprintf(stderr,"Considering scaffolded piece " F_CID "\n",
                piece->ident);
      }
    }
    // here, we'll call MergeMultiAligns on the meta-unitig.
   newMultiAlign = MergeMultiAligns(ScaffoldGraph->sequenceDB,
                                    ScaffoldGraph->fragStore,
                                    ContigPositions, FALSE, TRUE,
                                    GlobalData->aligner);
   fprintf(stderr," Returned from call to MergeMultiAlign\n");
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
  if (originalGaps == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for originalGaps\n");
	assert(0);
  }
  closedGap = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (closedGap == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for closedGap\n");
	assert(0);
  }
  closedGapDelta = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (closedGapDelta == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for closedGapDelta\n");
	assert(0);
  }  
  lcontigIdGap = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (lcontigIdGap == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for lcontigIdGap\n");
	assert(0);
  }
  rcontigIdGap = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (rcontigIdGap == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for rcontigIdGap\n");
	assert(0);
  }
  lcontigLength = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (lcontigLength == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for lcontigLength\n");
	assert(0);
  }
  rcontigLength = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (rcontigLength == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for rcontigLength\n");
	assert(0);
  }
  contigValid = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (contigValid == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for contigValid\n");
	assert(0);
  }
  allContigLengths = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (allContigLengths == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for allContigLengths\n");
	assert(0);
  }
  closedGapAhang = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (closedGapAhang == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for closedGapAhang\n");
	assert(0);
  }
  closedGapOlapLength = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (closedGapOlapLength == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for closedGapOlapLength\n");
	assert(0);
  }
  closedGapBhang = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (closedGapBhang == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for closedGapBhang\n");
	assert(0);
  }
  closedGapLcontigBasesIntact = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (closedGapLcontigBasesIntact == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for closedGapLcontigBasesIntact\n");
	assert(0);
  }
  closedGapRcontigBasesIntact = (int *) safe_malloc( GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  if (closedGapRcontigBasesIntact == NULL)
  {
	fprintf( stderr, "Could not safe_malloc space for closedGapRcontigBasesIntact\n");
	assert(0);
  }


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
	fprintf(stderr,  "=== examing scaffold " F_CID ", size %f\n", sid, scaff->bpLength.mean);
	  
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
                   }
                   cip.scaff_id=other_chunk->scaffoldID;
                   cip.scaff_bpLength= (CDS_COORD_t) oscaff->bpLength.mean;
                   cip.scaff_numCI=oscaff->info.Contig.numCI;
               } else {
                   cip.status = CGW_UNPLACED;
                   cip.scaff_id=-1;
                   cip.scaff_bpLength=0;
                   cip.scaff_numCI=0;
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
               fprintf(stderr,"  Putting unplaced or other chunk " F_CID " onto stack (from chunk " F_CID ")\n", 
                               other_chunk -> id, rcontig->id );
               AppendVA_ChunkPlacement(piece_list,&cip);
               potential_fill_count+=1;
               potential_fill_length+= (CDS_COORD_t) other_chunk->bpLength.mean;
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
             overlapFound= OverlapPieceList(piece_list,pieces,olaps);
             fprintf(stderr,"\t OverlapPieceList returned %d overlaps\n",overlapFound);
             if ( overlapFound > 0 ) {
               UnitigPieceList(pieces,olaps,iums,imps);
               DumpUnitigVAInfo ( piece_list, pieces, iums );
               { size_t i;
                 IntUnitigMesg *ium; 
                 for (i=0;i<GetNumIntUnitigMesgs(iums);i++) {
                   ium=GetIntUnitigMesg(iums,i);
                   if ( !CheckMetaUnitigOrderConsistency( piece_list,ium)) continue; 
                   if ( !CheckMetaUnitigOrientConsistency( piece_list,ium)) continue; 
                   if ( !CheckMetaUnitigPotential( piece_list,ium)) continue;
                   MergeMetaUnitigIntoContig(piece_list,ium);
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
      ofg.type = AS_CONTIG;
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
        if ( this_chunk != prev_chunk ) { //process, avoiding duplicates
          char stat;
          GetConsensus(ScaffoldGraph->RezGraph, this_chunk, consensusA, qualityA);
          seqs[iuniq] = (char *) safe_malloc((GetNumchars(consensusA)+1)*sizeof(char));
          seqlen[iuniq] = GetNumchars(consensusA);
          length+=seqlen[iuniq];
          strcpy(seqs[iuniq],Getchar(consensusA,0));
          ma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, this_chunk, FALSE);

          ofg.iaccession = iuniq+1;
          ofg.eaccession = this_chunk;
          ofg.clear_rng.bgn = 0;
          ofg.clear_rng.end = GetNumchars(consensusA);
          AppendVA_OFGMesg(pieces,&ofg);
          fprintf(stderr,"Chunk " F_UID " ",ofg.eaccession);
          if ( cistat == CGW_PLACED ) {
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
                   ofg.clear_rng.end,GetNumIntMultiPoss(ma->f_list),
                   GetNumIntUnitigPoss(ma->u_list));
	  InitContigTIterator(ScaffoldGraph, ofg.eaccession, TRUE, FALSE, &cis);
	  while(NULL != (ci = NextContigTIterator(&cis))){
	    CDS_CID_t cid = ci->id;
            fprintf(stderr,"	CI " F_CID ", cov: %d\n",cid,ci->info.CI.coverageStat);
          }
          SetVA_ChunkPlacement(piece_list,iuniq,cip);
          multialigns[iuniq++] = ma; 
        }
        prev_chunk = this_chunk;
      }
      ResetToRange_ChunkPlacement(piece_list,iuniq);
      fprintf(stderr,"\n" F_CID " pieces to unitig, combined length of " F_SIZE_T "\n",iuniq,length);
      for (i=0;i<iuniq-1;i++) {
         // do overlaps here
        int j;
        int i_placed=0,j_placed;
        int ori,where;
        left_ci=NULL;
        right_ci=NULL;
        afrag=multialigns[i];
        A.sequence = seqs[i];
        A.iaccession = i+1;
        A.eaccession = i+1;
        aseq=seqs[i];
        cip = GetChunkPlacement(piece_list,i);
        if ( cip->status == CGW_PLACED ) {
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
            cip = GetChunkPlacement(piece_list,i);
            if ( cip->status == CGW_PLACED ) {
              j_placed=1;
              if ( cip->order == 0 ) {
                 left_ci = cip;
              } else if (cip->order == is_placed-1) {
                 right_ci = cip;
              }
            }
            if ( ! (j_placed && i_placed ) ) {
           for (ori=0;ori<2;ori++){ 
            //for (ori=0;ori<1;ori++){ 
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
               fprintf(overlapfp,"Found overlap between chunks " F_CID " and " F_CID " (orient %d)\n",afrag->id,bfrag->id,ori);
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
