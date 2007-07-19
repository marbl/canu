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

static const char CM_ID[] = "$Id: eCR.c,v 1.23 2007-07-19 09:50:33 brianwalenz Exp $";

#include "eCR.h"
#include "ScaffoldGraph_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "PublicAPI_CNS.h"
#include "GapWalkerREZ.h"  //  FindGapLength

#define MAX_EXTENDABLE_FRAGS   100
#define NUM_STDDEV_CUTOFF        5.0

#define  USE_UNGAPPED_CONSENSUS_FOR_UNITIG


typedef struct extendableFragT {
  int fragIid;
  int extension;
  int addedBases;
  int basesToNextFrag;
  int fragOnEnd;
  int unitigID;
} extendableFrag;


typedef struct fragPositionsT {
  int bgn;
  int end;
} fragPositions;


typedef struct eCRstats {
} eCRstats;



int findFirstExtendableFrags(ContigT *contig, extendableFrag *extFragsArray);
int findLastExtendableFrags(ContigT *contig, extendableFrag *extFragsArray);

int findFirstUnitig(ContigT *contig, int *unitigID);
int findLastUnitig(ContigT *contig, int *unitigID);

void extendContig(ContigT *contig, int extendAEnd);

int GetNewUnitigMultiAlign(NodeCGW_T *unitig,
                           fragPositions *fragPoss,
                           int extendedFragIid);

void extendClearRange(int fragIid, int frag3pDelta);
void revertClearRange(int fragIid);

void getAlteredFragPositions(NodeCGW_T *unitig,
                             fragPositions **fragPoss,
                             int alteredFragIid,
                             int extension);

void leftShiftIUM(IntMultiPos *f_list, int numFrags, int extendedFragIid);
void rightShiftIUM(IntMultiPos *f_list, int numFrags, int extendedFragIid);

void saveFragAndUnitigData(int lFragIid, int rFragIid);
void restoreFragAndUnitigData(int lFragIid, int rFragIid);

//  In eCR-examineGap.c
//
void saveDefaultLocalAlignerVariables(void);
void restoreDefaultLocalAlignerVariables(void);
int  examineGap(ContigT *lcontig, int lFragIid,
                ContigT *rcontig, int rFragIid, 
                int  gapNumber,
                int *ahang,
                int *olapLengthOut,
                int *bhang,
                int *currDiffs,
                int *lcontigBasesIntact,
                int *rcontigBasesIntact,
                int *closedGapDelta,
                int  lBasesToNextFrag,
                int  rBasesToNextFrag,
                int *leftFragFlapLength,
                int *rightFragFlapLength);

//  In eCR-diagnostic.c
//
void DumpContigMultiAlignInfo (char *label, MultiAlignT *cma, int contigID);
void DumpUnitigInfo(char *label, NodeCGW_T *unitig);
void DumpContigUngappedOffsets(char *label, int contigID);


int                      totalContigsBaseChange = 0;
fragRecord              *fsread = NULL;
VA_TYPE(char)           *reformed_consensus = NULL;
VA_TYPE(char)           *reformed_quality = NULL;
VA_TYPE(int32)          *reformed_deltas = NULL;
VA_TYPE(IntElementPos)  *ContigPositions = NULL;
VA_TYPE(IntElementPos)  *UnitigPositions = NULL;

static MultiAlignT      *savedLeftUnitigMA = NULL;
static MultiAlignT      *savedRightUnitigMA = NULL;
static MultiAlignT      *savedLeftContigMA = NULL;
static MultiAlignT      *savedRightContigMA = NULL;

int                      iterNumber = 1;

debugflags_t             debug = {0, 0L, 0, 0L, 0, 0L};



//  Returns TRUE if the new unitig will cause problems in CreateAContigInScaffold().
//
//  This function is based on the failure in the following stack trace:
//    AppendFragToLocalStore()
//    MergeMultiAligns()
//    MergeMultiAlignsFast_new()
//    CreateAContigInScaffold()
//    main()
//
int
CheckNewUnitigMultiAlign_AppendFragToLocalStore(FragType type, int32 iid) {
  MultiAlignT *uma;

  int num_columns = 0;
  int num_frags   = 0;
  int num_unitigs = 0;

  VA_TYPE(int32) *gapped_positions = CreateVA_int32(num_columns+1);
  int32 ifrag;

  IntMultiPos  *frag;
  IntUnitigPos *unitig;

  uma     =  LoadMultiAlignTFromSequenceDB(sequenceDB, iid, type == AS_UNITIG);

  num_columns = GetMultiAlignLength(uma);
  num_frags   = GetNumIntMultiPoss(uma->f_list);
  num_unitigs = GetNumIntUnitigPoss(uma->u_list);

  frag = GetIntMultiPos(uma->f_list,0);

  for (ifrag=0;ifrag<num_frags;ifrag++,frag++){
    SetVA_int32(gapped_positions,frag->position.bgn,&frag->position.bgn);
    SetVA_int32(gapped_positions,frag->position.end,&frag->position.end);
  }

  unitig = GetIntUnitigPos(uma->u_list,0);

  for (ifrag=0;ifrag<num_unitigs;ifrag++,unitig++){
    SetVA_int32(gapped_positions,unitig->position.bgn,&unitig->position.bgn);
    SetVA_int32(gapped_positions,unitig->position.end,&unitig->position.end);
  }

  if ( Getint32(gapped_positions,num_columns) == NULL ) {
    fprintf(stderr,"Misformed Multialign... fragment positions only extend to bp %d out of %d\n",
            (int) GetNumint32s(gapped_positions),num_columns+1);
    DeleteVA_int32(gapped_positions);
    return -1;
  }

  DeleteVA_int32(gapped_positions);
  return(0);
}



int
CheckNewUnitigMultiAlign(CIScaffoldT *scaffold,
                         VA_TYPE(IntElementPos) *ContigPositions) {

  //  Our call in eCR that we need to check is:
  //
  //  success = CreateAContigInScaffold(scaff, ContigPositions, newOffsetAEnd, newOffsetBEnd);
  //

  //  CreateAContigInScaffold does nothing interesting at the start,
  //  it says it normalizes positions, but the actual normalization is
  //  commented out.

  //  MergeMultiAlignsFast_new


  static VA_TYPE(IntMultiPos)  *positions = NULL;
  static IntMultiPos            mpos;

  IntElementPos                *epos = GetIntElementPos(ContigPositions,0);
  int                           npos = GetNumIntElementPoss(ContigPositions);
  int                           i;

  int                           num_contigs;
  IntMultiPos                  *cpositions; 
  int                           num_columns = 0;

  mpos.contained    = 0;
  mpos.delta_length = 0;
  mpos.delta        = NULL;

  if (positions == NULL ) {
    positions = CreateVA_IntMultiPos(npos);
  } else {
    ResetVA_IntMultiPos(positions);
  }

  for (i=0; i<npos; i++,epos++) {
    mpos.type=epos->type;
    mpos.ident=epos->ident;
    mpos.position=epos->position;
    AppendVA_IntMultiPos(positions,&mpos);
  }

  num_contigs = GetNumIntMultiPoss(positions);
  cpositions  = GetIntMultiPos(positions,0);

  USE_SDB           = 1;
  ALIGNMENT_CONTEXT = AS_MERGE;
  sequenceDB        = ScaffoldGraph->sequenceDB;

  for (i=0;i<num_contigs;i++) {
    num_columns = (cpositions[i].position.bgn>num_columns) ? cpositions[i].position.bgn : num_columns;
    num_columns = (cpositions[i].position.end>num_columns) ? cpositions[i].position.end : num_columns;
  }

  if (num_contigs == 1) {
    //  Shouldn't happen!
    assert(0);
  } else {
    for (i=0;i<num_contigs;i++) {
      if (CheckNewUnitigMultiAlign_AppendFragToLocalStore(cpositions[i].type, 
                                                          cpositions[i].ident))
        return(TRUE);
    }
  }
     
  return(FALSE);
}


void
SynchUnitigTWithMultiAlignT(NodeCGW_T *unitig) {
  unitig->bpLength.mean = GetMultiAlignUngappedLength(LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
                                                                                    unitig->id,
                                                                                    TRUE));
}



int
main(int argc, char **argv) {

  //  Command line args and processing
  int   scaffoldBegin    = -1;
  int   startingGap      = -1;
  int   scaffoldEnd      = -1;
  int   ckptNum          = -1;
  int   arg              = 1;
  int   err              = 0;

  //  Loop counters
  int   sid         = 0;
  int   gapNumber   = 0;

  //  Statistics (in no apparent order)
  //
  int    numExtendableGaps = 0;

  int    leftContigExtendable = 0;
  int    rightContigExtendable = 0;
  int    bothContigExtendable = 0;

  int    surrogatesOnEnd = 0;

  int    numGaps = 0;
  int    numGapsClosed = 0;
  int    numGapsVarTooSmall = 0;

  int    totalBasesInClosedGaps = 0;

  int    totalOlapLength = 0;
  int    totalOlapDiffs = 0;
  double totalOlapVariance = 0.0;

  int    numClosingsTried = 0;
  int    unitigMultiAlignFailures = 0;
  int    replaceEndUnitigFailures = 0;
  int    createAContigFailures = 0;
  int    unitigToContigFailures = 0;
  int    noOverlapFound = 0;

  int    numSmallGaps = 0;
  int    numSmallGapsClosed = 0;

  int    numLargeGaps = 0;
  int    numLargeGapsClosed = 0;

  double maxGapSizeClosed = 0.0;
  int    maxGapSizeClosedNumber = -1;


  debug.eCRmainFP    = stderr;
  debug.examineGapFP = stderr;
  debug.diagnosticFP = stderr;



  // save off whatever the rest of the world has for default values
  // for Local_Overlap_AS_forCNS
  //
  saveDefaultLocalAlignerVariables();

  GlobalData           = CreateGlobal_CGW();
  GlobalData->stderrc  = stderr;
  GlobalData->timefp   = stderr;


  while (arg < argc) {
    if        (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->File_Name_Prefix, argv[++arg]);
    } else if (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->Gatekeeper_Store_Name, argv[++arg]);
    } else if (strcmp(argv[arg], "-C") == 0) {
      startingGap = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-n") == 0) {
      ckptNum = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-b") == 0) {
      scaffoldBegin = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      scaffoldEnd   = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-i") == 0) {
      iterNumber = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-v") == 0) {
      debug.eCRmainLV    = 1;
      debug.examineGapLV = 1;
    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((iterNumber < 1) || (iterNumber > 2)) {
    fprintf(stderr, "%s: ERROR!  Invalid iteration %d.\n", argv[0], iterNumber+1);
    err++;
  }
  if ((GlobalData->File_Name_Prefix[0] == 0) ||
      (GlobalData->Gatekeeper_Store_Name[0] == 0) ||
      (err)) {
    fprintf(stderr, "usage: %s [opts] -c ckpName -n ckpNumber -g gkpStore\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c ckpName   Use ckpName as the checkpoint name\n");
    fprintf(stderr, "  -n ckpNumber The checkpoint to use\n");
    fprintf(stderr, "  -g gkpStore  The gatekeeper store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -C gap#      Start at a specific gap number\n");
    fprintf(stderr, "  -b scafBeg   Begin at a specific scaffold\n");
    fprintf(stderr, "  -e scafEnd   End at a specific scaffold\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -i iterNum   The iteration of ECR; either 1 or 2\n");
    exit(1);
  }
  

  //  This is used all over the place, do not remove it!
  //
  fsread = new_fragRecord();


  //  LoadScaffoldGraphFromCheckpoint wants to CheckCIScaffoldT()
  //  which can RecomputeOffsetsInScaffold(), which can eventually,
  //  try to get an overlap.  Unless we set 'alligner', that will
  //  crash.
  //
  //  After the graph is loaded, we reopen the gatekeeper store for
  //  read/write.
  //
  GlobalData->aligner=Local_Overlap_AS_forCNS;
  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(GlobalData->File_Name_Prefix, ckptNum, TRUE);

  closeGateKeeperStore(ScaffoldGraph->gkpStore);

  ScaffoldGraph->gkpStore = openGateKeeperStore(GlobalData->Gatekeeper_Store_Name, TRUE);
  if(ScaffoldGraph->gkpStore == NULL){
    fprintf(stderr, "%s: Failed to open the gatekeeper store '%s' for read/write.  Bye.\n", argv[0], GlobalData->Gatekeeper_Store_Name);
    exit(1);
  }

  //  Update the begin/end scaffold ids.
  //
  if (scaffoldBegin == -1)
    scaffoldBegin = 0;
  if (scaffoldEnd == -1)
    scaffoldEnd = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
  if (scaffoldEnd > GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph))
    scaffoldEnd = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);


  //  Intiialize the variable arrays
  //
  reformed_consensus = CreateVA_char(256 * 1024);
  reformed_quality   = CreateVA_char(256 * 1024);
  reformed_deltas    = CreateVA_int32(8192);

  //
  //  Scan all the scaffolds, closing gaps.  Go!
  //

  for (sid = scaffoldBegin; sid < scaffoldEnd; sid++) {
    CIScaffoldT    *scaff            = NULL;
    ContigT        *lcontig          = NULL;
    int             lcontigID        = 0;
    ContigT        *rcontig          = NULL;
    int             rcontigID        = 0;

    int             numSmallGapsThisScaff       = 0;
    int             numSmallGapsClosedThisScaff = 0;
    int             numLargeGapsThisScaff       = 0;
    int             numLargeGapsClosedThisScaff = 0;

    extendableFrag  leftExtFragsArray[MAX_EXTENDABLE_FRAGS];
    extendableFrag  rightExtFragsArray[MAX_EXTENDABLE_FRAGS];


    scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
    assert(scaff != NULL);

    // not interested in dead scaffold, not real scaffolds, or singleton scaffolds

    if ((isDeadCIScaffoldT(scaff)) ||
        (scaff->type != REAL_SCAFFOLD) ||
        (scaff->info.Scaffold.numElements < 2))
      continue;

    fprintf(stderr,"\n=====================================================================\n");
    fprintf(stderr,"examing scaffold %d, size %f\n", sid, scaff->bpLength.mean);

    lcontig   = GetGraphNode(ScaffoldGraph->ContigGraph, scaff->info.Scaffold.AEndCI);
    lcontigID = lcontig->id;
    rcontigID = lcontig->BEndNext;

    while (rcontigID != -1) {
      NodeOrient lcontigOrientation  = B_A;
      NodeOrient rcontigOrientation  = B_A;

      LengthT    gapSize;
      int        numLeftFrags        = 0;
      int        numRightFrags       = 0;

      int        lFragIid            = -1;
      int        rFragIid            = -1;

      int        lunitigID           = 0;
      int        runitigID           = 0;

      ContigT    lcontigBackup;
      ContigT    rcontigBackup;

      int        leftFragIndex       = 0;
      int        rightFragIndex      = 0;

      //  Mostly args to examineGap()
      int        ahang               = 0;
      int        currLength          = 0;
      int        bhang               = 0;
      int        currDiffs           = 0;
      int        lcontigBasesIntact  = 0;
      int        rcontigBasesIntact  = 0;
      int        closedGap           = FALSE;
      int        closedGapDelta      = 0;
      int        leftFragFlapLength  = 0;
      int        rightFragFlapLength = 0;



      numGaps++;

      rcontig = GetGraphNode(ScaffoldGraph->ContigGraph, rcontigID);

      gapSize = FindGapLength(lcontig, rcontig, FALSE);

      if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
        lcontigOrientation = A_B;

      if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
        rcontigOrientation = A_B;

      fprintf(stderr, "---------------------------------------------------------------\n");
      fprintf(stderr, "examining gap %d from lcontig %d (orient: %c, pos: %f, %f) to rcontig %d (orient: %c, pos: %f, %f), size: %lf \n", 
              gapNumber, 
              lcontig->id, lcontigOrientation, lcontig->offsetAEnd.mean, lcontig->offsetBEnd.mean,
              rcontig->id, rcontigOrientation, rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean,
              gapSize.mean);

      if (gapSize.mean < 100.0) {
        numSmallGaps++;
        numSmallGapsThisScaff++;
      } else {
        numLargeGaps++;
        numLargeGapsThisScaff++;
      }


      if (debug.eCRmainLV > 0)
        fprintf(debug.eCRmainFP, "\nexamining lcontig %d (orientation %c) \n", lcontig->id, lcontigOrientation);

      // find the extreme read on the correct end of the lcontig
      if (lcontigOrientation == A_B) {
        numLeftFrags = findLastExtendableFrags(lcontig, leftExtFragsArray);
        if (findLastUnitig(lcontig, &lunitigID)) {
          surrogatesOnEnd++;
          numLeftFrags = 0;
        }
      } else {
        numLeftFrags = findFirstExtendableFrags(lcontig, leftExtFragsArray);
        if (findFirstUnitig(lcontig, &lunitigID)) {
          surrogatesOnEnd++;
          numLeftFrags = 0;
        }
      }

      if (debug.eCRmainLV > 0)
        fprintf(debug.eCRmainFP, "finished examining lcontig %d (orientation %c) \n", lcontig->id, lcontigOrientation);
  


      if (debug.eCRmainLV > 0)
        fprintf(debug.eCRmainFP, "\nexamining rcontig %d (orientation %c) \n", rcontig->id, rcontigOrientation);

      // find the extreme read on the correct end of the rchunk
      if (rcontigOrientation == A_B) {
        numRightFrags = findFirstExtendableFrags(rcontig, rightExtFragsArray);
        if (findFirstUnitig(rcontig, &runitigID)) {
          surrogatesOnEnd++;
          numRightFrags = 0;
        }
      } else {
        numRightFrags = findLastExtendableFrags(rcontig, rightExtFragsArray);
        if (findLastUnitig(rcontig, &runitigID)) {
          surrogatesOnEnd++;
          numRightFrags = 0;
        }
      }

      if (debug.eCRmainLV > 0)
        fprintf(debug.eCRmainFP, "finished examining rcontig %d (orientation %c) \n", rcontig->id, rcontigOrientation);



      //  Skip this gap if it is before the one we want to start at.
      //
      if (gapNumber < startingGap)
        numLeftFrags = numRightFrags = 0;


      //  This is a good place to check rcontig->id, lcontig->id,
      //  gapNumber, etc and abort if that gap is causing problems.
      //
      //if (rcontig->id == 405527)
      //  numLeftFrags = numRightFrags = 0;


      numExtendableGaps++;

      if (numLeftFrags > 0)
        leftContigExtendable++;
      if (numRightFrags > 0)
        rightContigExtendable++;
      if (numLeftFrags > 0 && numRightFrags > 0)
        bothContigExtendable++;


      if (debug.eCRmainLV > 0)
        fprintf(debug.eCRmainFP, "gap (%d, %d) has extendable frag pair on left: %d, on right: %d\n", 
                lcontig->id, rcontig->id, numLeftFrags, numRightFrags);

		
      // set an extra member of the arrays for when the contig is not being extended
      leftExtFragsArray[numLeftFrags].fragIid   = -1;
      leftExtFragsArray[numLeftFrags].extension = 0;
      numLeftFrags++;

      rightExtFragsArray[numRightFrags].fragIid   = -1;		  
      rightExtFragsArray[numRightFrags].extension = 0;
      numRightFrags++;


      //  In case we need to back out changes late in the extension, we save copies of the two contigs.
      //
      memcpy(&lcontigBackup, lcontig, sizeof(ContigT));
      memcpy(&rcontigBackup, rcontig, sizeof(ContigT));

			
      for (leftFragIndex = 0; leftFragIndex < numLeftFrags && closedGap == FALSE; leftFragIndex++) {
        for (rightFragIndex = 0; rightFragIndex < numRightFrags && closedGap == FALSE; rightFragIndex++) {

          if (debug.eCRmainLV > 0)
            fprintf(debug.eCRmainFP, "examining frags %d and %d\n", leftExtFragsArray[leftFragIndex].fragIid, 
                    rightExtFragsArray[rightFragIndex].fragIid);

          if ((gapSize.mean - leftExtFragsArray[leftFragIndex].extension - rightExtFragsArray[rightFragIndex].extension) > (NUM_STDDEV_CUTOFF * sqrt(gapSize.variance)) &&
              (gapSize.mean > 100.0)) {

            if (debug.eCRmainLV > 0)
              fprintf(debug.eCRmainFP, "leftExtFragsArray[%d].extension: %10d, rightExtFragsArray[%d].extension: %10d\n",
                      leftFragIndex, leftExtFragsArray[leftFragIndex].extension, 
                      rightFragIndex, rightExtFragsArray[rightFragIndex].extension);

            if (debug.eCRmainLV > 0)
              fprintf(debug.eCRmainFP, "gap variance too large (gapSize - extensions: %.2f, %.1f * sqrt(gapSize.variance): %.2f\n",
                      gapSize.mean - leftExtFragsArray[leftFragIndex].extension - 
                      rightExtFragsArray[rightFragIndex].extension, 
                      NUM_STDDEV_CUTOFF, NUM_STDDEV_CUTOFF * sqrt(gapSize.variance));

            numGapsVarTooSmall++;
            closedGap = FALSE;

            continue;
          }

          lFragIid = leftExtFragsArray[leftFragIndex].fragIid;
          rFragIid = rightExtFragsArray[rightFragIndex].fragIid;

          // have to check and make sure that the frags belong to the correct unitig
          if (lFragIid != -1) {
            InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);
            assert(info->set);
            if (GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex)->cid != lunitigID)
              continue;
          }
			
          if (rFragIid != -1) {
            InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
            assert(info->set);
            if (GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex)->cid != runitigID)
              continue;
          }

          numClosingsTried++;
			  
          //dumpContigInfo(lcontig);
          //dumpContigInfo(rcontig);

          if (examineGap(lcontig, lFragIid,
                         rcontig, rFragIid, 
                         gapNumber,
                         &ahang,
                         &currLength,
                         &bhang,
                         &currDiffs,
                         &lcontigBasesIntact,
                         &rcontigBasesIntact,
                         &closedGapDelta,
                         leftExtFragsArray[leftFragIndex].basesToNextFrag, 
                         rightExtFragsArray[rightFragIndex].basesToNextFrag,
                         &leftFragFlapLength,
                         &rightFragFlapLength)) {

            int        keepGap             = TRUE;

            int        lextension          = 0;
            int        rextension          = 0;

            LengthT    newOffsetAEnd;
            LengthT    newOffsetBEnd;
            double     maxRContigOffset;
            int32      nextContigIndex;



            // save off copies of everything we might alter so we can
            // restore if gap closing fails
            //
            saveFragAndUnitigData(lFragIid, rFragIid);

            totalBasesInClosedGaps += (int) gapSize.mean;
				
            if (CONTIG_BASES < 2000) {
              if ((ahang + currLength + bhang - 1000) > (MIN(CONTIG_BASES, (int) lcontig->bpLength.mean) +
                                                         MIN(CONTIG_BASES, (int) rcontig->bpLength.mean))) {

                if (debug.eCRmainLV > 0)
                  fprintf(debug.eCRmainFP, "at gapNumber %d, ahang + currLength + bhang - 1000 = %d, back(A) + back(B) = %d\n",
                          gapNumber, ahang + currLength + bhang - 1000,
                          MIN(CONTIG_BASES, (int) lcontig->bpLength.mean) + MIN(CONTIG_BASES, (int) rcontig->bpLength.mean));

                keepGap = FALSE;
              }
				
              if ((ahang + currLength + bhang + 700) < (MIN(CONTIG_BASES, (int) lcontig->bpLength.mean) +
                                                        MIN(CONTIG_BASES, (int) rcontig->bpLength.mean))) {

                if (debug.eCRmainLV > 0)
                  fprintf(debug.eCRmainFP, "at gapNumber %d, ahang + currLength + bhang + 500 = %d, back(A) + back(B) = %d\n",
                          gapNumber, ahang + currLength + bhang + 500,
                          MIN(CONTIG_BASES, (int) lcontig->bpLength.mean) + MIN(CONTIG_BASES, (int) rcontig->bpLength.mean));

                keepGap = FALSE;
              }
            }

            // extend the clear ranges of the frags
            //
            if (keepGap) {

              //  The previous version merged the ExtFragsArray[]
              //  test with the extendClearRange() -- so if the left
              //  frag updated, but the right frag failed, we might
              //  leave things inconsistent.  Hopefully, we restore
              //  the original clear ranges at the end of all this.
              //  But we still fix it to not modify unless both
              //  frags are OK.

              if (lFragIid != -1)
                if (leftExtFragsArray[leftFragIndex].addedBases - leftFragFlapLength < 0)
                  keepGap = FALSE;

              if (rFragIid != -1)
                if (rightExtFragsArray[rightFragIndex].addedBases - rightFragFlapLength < 0)
                  keepGap = FALSE;

              if (keepGap) {
                if (lFragIid != -1) {
                  if (debug.eCRmainLV > 0)
                    fprintf(debug.eCRmainFP,"adjusting left frg clear range by %d - %d = %d bases\n",
                            leftExtFragsArray[leftFragIndex].addedBases, leftFragFlapLength,
                            leftExtFragsArray[leftFragIndex].addedBases - leftFragFlapLength);

                  extendClearRange(lFragIid, leftExtFragsArray[leftFragIndex].addedBases - leftFragFlapLength);
                }

                if (rFragIid != -1) {
                  if (debug.eCRmainLV > 0)
                    fprintf(debug.eCRmainFP,"adjusting right frg clear range by %d - %d = %d bases\n",
                            rightExtFragsArray[rightFragIndex].addedBases, rightFragFlapLength,
                            rightExtFragsArray[rightFragIndex].addedBases - rightFragFlapLength);

                  extendClearRange(rFragIid, rightExtFragsArray[rightFragIndex].addedBases - rightFragFlapLength); 
                }
              }
            }


            if (keepGap) { // the fragment extensions have succeeded
              int             gotNewLeftMA        = TRUE;
              int             gotNewRightMA       = TRUE;

              CIFragT        *frag                = NULL;
              NodeCGW_T      *unitig              = NULL;
              MultiAlignT    *new_cma             = NULL;
              fragPositions  *fragPoss            = NULL;
				  
              if (debug.diagnosticLV > 0) {
                DumpContigMultiAlignInfo ("before anything (left)", NULL, lcontig->id);
                DumpContigUngappedOffsets("before anything (left)", lcontig->id);
                DumpContigMultiAlignInfo ("before anything (right)", NULL, rcontig->id);
                DumpContigUngappedOffsets("before anything (right)", rcontig->id);
              }

              if (debug.eCRmainLV > 0) {
                fprintf(debug.eCRmainFP, "before altering, lctg: %12.0f, %12.0f\n",
                        (lcontigOrientation == A_B) ? lcontig->offsetAEnd.mean : lcontig->offsetBEnd.mean,
                        (lcontigOrientation == A_B) ? lcontig->offsetBEnd.mean : lcontig->offsetAEnd.mean);
                fprintf(debug.eCRmainFP, "                 rctg: %12.0f, %12.0f\n",
                        (rcontigOrientation == A_B) ? rcontig->offsetAEnd.mean : rcontig->offsetBEnd.mean,
                        (rcontigOrientation == A_B) ? rcontig->offsetBEnd.mean : rcontig->offsetAEnd.mean);
              }

              // save the max offset of the right contig so we know
              // how to adjust the offsets of the contigs further
              // along the scaffold later

              maxRContigOffset = MAX(rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean);
              nextContigIndex = rcontig->BEndNext;
				  
              // alter the unitigs in cgw memory struct land

              // left unitig
              if (lFragIid != -1) {

                //  This is right before the big pile of overlap
                //  output.  Put in some more prints to isolate
                //  where those are, then figure out where "for
                //  unitig " gets printed then find why the
                //  fragment isn't being extended to the end of
                //  the contig.

                InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);
                assert(info->set);
                frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

                unitig = GetGraphNode(ScaffoldGraph->CIGraph, frag->CIid);
                //extendUnitigs(unitig, lFragIid, leftExtFragsArray[leftFragIndex], TRUE);
                getAlteredFragPositions(unitig, &fragPoss, lFragIid, 
                                        leftExtFragsArray[leftFragIndex].extension - leftFragFlapLength);

                gotNewLeftMA = GetNewUnitigMultiAlign(unitig, fragPoss, lFragIid);

                safe_free(fragPoss);  // inefficient, let's just do a big array once that get's reused

                if (!gotNewLeftMA)
                  unitigMultiAlignFailures++;

                if (gotNewLeftMA) {
                  SynchUnitigTWithMultiAlignT(unitig);

                  if (debug.diagnosticLV > 0)
                    DumpContigMultiAlignInfo ("before ReplaceEndUnitigInContig (left)", NULL, lcontig->id);

                  new_cma = ReplaceEndUnitigInContig(ScaffoldGraph->sequenceDB,
                                                     ScaffoldGraph->gkpStore,
                                                     lcontig->id, unitig->id,
                                                     lcontig->offsetAEnd.mean >= lcontig->offsetBEnd.mean,
                                                     GlobalData->aligner,
                                                     NULL);

                  if (new_cma) {
                    if (debug.diagnosticLV > 0) {
                      DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (left)", NULL, lcontig->id);
                      DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (left)", new_cma, lcontig->id);
                    }

                    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, lcontig->id, FALSE);
                    InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, lcontig->id, 
                                                  FALSE, new_cma, TRUE);


                    if (debug.eCRmainLV > 0)
                      fprintf(debug.eCRmainFP, "strlen(Getchar (new_cma->consensus)): " F_SIZE_T "\n",
                              strlen(Getchar (new_cma->consensus, 0)));

                    if (debug.diagnosticLV > 0)
                      DumpContigMultiAlignInfo ("after updating store (left)", NULL, lcontig->id);

                  } else {
                    fprintf(stderr, "WARNING:  Failed to align the new and old contigs.  Will not use this extension.\n");
                    replaceEndUnitigFailures++;
                    gotNewLeftMA = FALSE;
                  }

                  if (lcontigOrientation == A_B)
                    extendContig(lcontig, FALSE);
                  else
                    extendContig(lcontig, TRUE);
                }
              }
				  
              // right unitig
              if (rFragIid != -1) {
                InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
                assert(info->set);
                frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
					
                unitig = GetGraphNode(ScaffoldGraph->CIGraph, frag->CIid);
                //extendUnitigs(unitig, rFragIid, rightExtFragsArray[rightFragIndex], FALSE);
					
                getAlteredFragPositions(unitig, &fragPoss, rFragIid, 
                                        rightExtFragsArray[rightFragIndex].extension- rightFragFlapLength);

                gotNewRightMA = GetNewUnitigMultiAlign(unitig, fragPoss, rFragIid);

                safe_free(fragPoss);  // inefficient, let's just do a big array once that get's reused

                if (!gotNewRightMA)
                  unitigMultiAlignFailures++;

                if (gotNewRightMA) {
                  SynchUnitigTWithMultiAlignT(unitig);

                  if (debug.diagnosticLV > 0)
                    DumpContigMultiAlignInfo ("before ReplaceEndUnitigInContig (right)", NULL, rcontig->id);

                  new_cma = ReplaceEndUnitigInContig(ScaffoldGraph->sequenceDB,
                                                     ScaffoldGraph->gkpStore,
                                                     rcontig->id, unitig->id,
                                                     rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean,
                                                     GlobalData->aligner,
                                                     NULL);

                  if (new_cma) {
                    if (debug.diagnosticLV > 0) {
                      DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (right) (original contig)", NULL, rcontig->id);
                      DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (right) (new contig)", new_cma, rcontig->id);
                    }

                    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, rcontig->id, FALSE);
                    InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, rcontig->id, 
                                                  FALSE, new_cma, TRUE);
                    if (debug.eCRmainLV > 0)
                      fprintf(debug.eCRmainFP, "strlen(Getchar (new_cma->consensus)): " F_SIZE_T "\n",
                              strlen(Getchar (new_cma->consensus, 0)));

                    if (debug.diagnosticLV > 0)
                      DumpContigMultiAlignInfo ("after updating store (right)", NULL, rcontig->id);

                  } else {
                    fprintf(stderr, "WARNING:  Failed to align the new and old contigs.  Will not use this extension.\n");
                    replaceEndUnitigFailures++;
                    gotNewRightMA = FALSE;
                  }

                  // updateIntUnitigPoss(rcontig);
                  if (rcontigOrientation == A_B)
                    extendContig(rcontig, TRUE);
                  else
                    extendContig(rcontig, FALSE);
                }
              }  //  end of right unitig


              //  unwind changes if one or the other fails
              //
              if (gotNewLeftMA == FALSE || gotNewRightMA == FALSE)  
                keepGap = FALSE;
            }

            if (keepGap) {
              double delta;
              int success;

              // now shift the right contig into place
              rcontig->offsetAEnd.mean += closedGapDelta;
              rcontig->offsetBEnd.mean += closedGapDelta;				  

              if (debug.eCRmainLV > 0) {
                fprintf(debug.eCRmainFP, "after altering, lctg: %12.0f, %12.0f\n",
                        (lcontigOrientation == A_B) ? lcontig->offsetAEnd.mean : lcontig->offsetBEnd.mean,
                        (lcontigOrientation == A_B) ? lcontig->offsetBEnd.mean : lcontig->offsetAEnd.mean);
                fprintf(debug.eCRmainFP, "                rctg: %12.0f, %12.0f\n",
                        (rcontigOrientation == A_B) ? rcontig->offsetAEnd.mean : rcontig->offsetBEnd.mean,
                        (rcontigOrientation == A_B) ? rcontig->offsetBEnd.mean : rcontig->offsetAEnd.mean);
              }

              // setup for contig merge
              if (ContigPositions == NULL)
                ContigPositions = CreateVA_IntElementPos(2);
              ResetVA_IntElementPos(ContigPositions);
				  
              if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
                delta = lcontig->offsetAEnd.mean;
              else
                delta = lcontig->offsetBEnd.mean;

              {
                IntElementPos   contigPos;

                contigPos.ident = lcontig->id;
                contigPos.type = AS_CONTIG;
                contigPos.position.bgn = lcontig->offsetAEnd.mean - delta;
                contigPos.position.end = lcontig->offsetBEnd.mean - delta;
                AppendIntElementPos(ContigPositions, &contigPos);
				  
                if (debug.eCRmainLV > 0)
                  fprintf(debug.eCRmainFP, "lcontig %8d positioned at %8d, %8d\n", 
                          lcontig->id,contigPos.position.bgn, contigPos.position.end);
              }

              {
                IntElementPos   contigPos;
    
                contigPos.ident = rcontig->id;
                contigPos.type = AS_CONTIG;
                contigPos.position.bgn = rcontig->offsetAEnd.mean - delta;
                contigPos.position.end = rcontig->offsetBEnd.mean - delta;
                AppendIntElementPos(ContigPositions, &contigPos);

                if (debug.eCRmainLV > 0)
                  fprintf(debug.eCRmainFP, "rcontig %8d positioned at %8d, %8d\n", 
                          rcontig->id,contigPos.position.bgn, contigPos.position.end);
              }

              if (lcontigOrientation == A_B)
                newOffsetAEnd.mean = lcontig->offsetAEnd.mean;
              else
                newOffsetAEnd.mean = lcontig->offsetBEnd.mean;
				  
              if (rcontigOrientation == A_B)
                newOffsetBEnd.mean = rcontig->offsetBEnd.mean;
              else
                newOffsetBEnd.mean = rcontig->offsetAEnd.mean;

              if (debug.diagnosticLV > 0) {
                DumpContigMultiAlignInfo ("before CreateAContigInScaffold", NULL, lcontig->id);
                DumpContigMultiAlignInfo ("before CreateAContigInScaffold", NULL, rcontig->id);
                //DumpContigUngappedOffsets(rcontig->id);
              }

              if (debug.eCRmainLV > 0)
                fprintf(debug.eCRmainFP, "CreateAContigInScaffold()-- newOffsetAEnd=%d newOffsetBEnd=%d\n",
                        (int)newOffsetAEnd.mean, (int)newOffsetBEnd.mean);

              // have to call this routine with normalized positions

              if (CheckNewUnitigMultiAlign(scaff, ContigPositions)) {
                fprintf(stderr, "CheckNewUnitigMultiAlign()-- The new unitig multialignment is messed up, will not close this gap.\n");
                unitigToContigFailures++;
                success = FALSE;
              } else {
                success = CreateAContigInScaffold(scaff, ContigPositions, newOffsetAEnd, newOffsetBEnd);
              }

              if (!success) {
                // why does this happen?
                // now we have to unwind everything we've done up above???

                keepGap = FALSE;
                fprintf(stderr, "overlap found in extendClearRanges but not by CNS!\n");

                createAContigFailures++;
              }

              //  XXX: BPW - why are these here?  is this an
              //  attempt to undo stuff?  just making sure they
              //  are up to date?
              //
              rcontig = GetGraphNode(ScaffoldGraph->ContigGraph, rcontigID);
              lcontig = GetGraphNode(ScaffoldGraph->ContigGraph, lcontigID);

              if (keepGap) {
                fprintf(stderr, "closed gap %8d, contigs %8d and %8d, fragIids %9d and %9d\n",
                        gapNumber, lcontig->id, rcontig->id, lFragIid, rFragIid);
					
                closedGap = TRUE;
					
                numGapsClosed++;
                totalOlapVariance += currLength * currLength;
                totalOlapLength += currLength;
                totalOlapDiffs += currDiffs;
					
                if (gapSize.mean < 100.0) {
                  numSmallGapsClosed++;
                  numSmallGapsClosedThisScaff++;
                } else {
                  numLargeGapsClosed++;
                  numLargeGapsClosedThisScaff++;
                  if (gapSize.mean > maxGapSizeClosed) {
                    maxGapSizeClosed = gapSize.mean;
                    maxGapSizeClosedNumber = gapNumber;
                  }
                }
					
                // The new contig takes place of what was the right contig.
                //
                rcontig = GetGraphNode(ScaffoldGraph->ContigGraph, GetNumGraphNodes(ScaffoldGraph->ContigGraph) - 1);
                //adjustUnitigCoords(rcontig);

                // now we need to adjust contigs past the new contig, if not on end
                if (nextContigIndex != NULLINDEX) {
                  LengthT scaffoldDelta;
					  
                  scaffoldDelta.mean = MAX(rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean) - maxRContigOffset;
                  scaffoldDelta.variance = ComputeFudgeVariance(scaffoldDelta.mean);
                  AddDeltaToScaffoldOffsets(ScaffoldGraph, scaff->id, nextContigIndex,
                                            TRUE, FALSE, scaffoldDelta);
                }
              }  //  keep gap
            }  //  keep gap

            // if we didn't close gap for whatever reason undo all the frag and unitig changes
            if (keepGap == FALSE) {
              //fprintf(stderr, "did not close gap %8d, contigs %8d and %8d\n", gapNumber, lcontig->id, rcontig->id);

              closedGap = FALSE;

              revertClearRange(lFragIid);
              revertClearRange(rFragIid);				  

              restoreFragAndUnitigData(lFragIid, rFragIid);

              memcpy(lcontig, &lcontigBackup, sizeof(ContigT));
              memcpy(rcontig, &rcontigBackup, sizeof(ContigT));
            }
            //  end of if (examineGap())
          } else {
            //  examine gap failed
            //fprintf(stderr, "did not close gap %8d, contigs %8d and %8d\n", gapNumber, lcontig->id, rcontig->id);
            closedGap = FALSE;
            noOverlapFound++;
          }
        }  //  over all right frags
      }  //  over all left frags


      gapNumber++;
      lcontig = rcontig;
      lcontigID = lcontig->id;
      rcontigID = lcontig->BEndNext;
    }  //  over all contigs in the scaffold

    fprintf(stderr, "scaffold stats, scaff %10d, smallGaps %8d closed %8d, largeGaps %8d closed %8d\n",
            scaff->id, numSmallGapsThisScaff, numSmallGapsClosedThisScaff,
            numLargeGapsThisScaff, numLargeGapsClosedThisScaff);

    fprintf(stderr, "after gapNumber %d:\n", gapNumber);
    fprintf(stderr, "             numSmallGaps: %d (closed %d, %.2f%%)\n", 
            numSmallGaps, numSmallGapsClosed,
            (numSmallGaps == 0) ? 0 : 100.0 * numSmallGapsClosed / numSmallGaps);
    fprintf(stderr, "             numLargeGaps: %d (closed %d, %.2f%%)\n", 
            numLargeGaps, numLargeGapsClosed,
            (numLargeGaps == 0) ? 0 : 100.0 * numLargeGapsClosed / numLargeGaps);
    fprintf(stderr, "                  allGaps: %d (closed %d, %.2f%%)\n", 
            numSmallGaps + numLargeGaps, numSmallGapsClosed + numLargeGapsClosed, 
            (numSmallGaps + numLargeGaps == 0) ? 0 : 100.0 * (numSmallGapsClosed + numLargeGapsClosed) / (numSmallGaps + numLargeGaps));



    if (numSmallGapsClosedThisScaff + numLargeGapsClosedThisScaff > 0) {
      int status = RECOMPUTE_SINGULAR;
      int recomputeIteration = 0;

      while ((recomputeIteration++ < 3) &&
             (status == RECOMPUTE_SINGULAR ||
              status == RECOMPUTE_CONTIGGED_CONTAINMENTS)) {

        // need to make sure scaffold is connected with trusted raw edges

        MarkInternalEdgeStatus(ScaffoldGraph,
                               GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid),
                               PAIRWISECHI2THRESHOLD_CGW,
                               SLOPPY_EDGE_VARIANCE_THRESHHOLD,
                               TRUE, TRUE, 0, TRUE);

        assert(IsScaffoldInternallyConnected(ScaffoldGraph,
                                             GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid),
                                             ALL_EDGES));
	      
        status = RecomputeOffsetsInScaffold(ScaffoldGraph,
                                            GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid),
                                            TRUE, TRUE, FALSE);
      }
    }


    //  Good place to checkpoint the scaffold graph -- we'd also need
    //  to checkpoint the gkpstore Eventually, we'll just write a log
    //  of the gkpstore changes, and apply that when we're all done.

    //  Clear out any cached multialigns.  We're all done with them.
    //
    ClearCacheSequenceDB(ScaffoldGraph->sequenceDB, TRUE);
    ClearCacheSequenceDB(ScaffoldGraph->sequenceDB, FALSE);

  }  //  over all scaffolds




  //  Loading a checkpoint implicitly calls these -- and the
  //  downstream consumer of our checkpoint shouldn't modify the
  //  checkpoint if it is querying it (dumpDistanceUpdates, for
  //  example) -- so we just call them before the checkpoint is
  //  written.
  //
  SetCIScaffoldTLengths(ScaffoldGraph, TRUE);
  CheckCIScaffoldTs(ScaffoldGraph);

  fprintf(GlobalData->timefp, "checkpoint %d written at end of extendClearRanges\n", ScaffoldGraph->checkPointIteration);
  CheckpointScaffoldGraph(ScaffoldGraph, 1);

  DestroyScaffoldGraph(ScaffoldGraph);
  DeleteGlobal_CGW(GlobalData);




  // Variance = mean(x^2) - (mean(x))^2
  totalOlapVariance = -0.0;
  if (numGapsClosed > 0)
    totalOlapVariance = ((totalOlapVariance / numGapsClosed) - 
                         ((double)totalOlapLength / numGapsClosed) *
                         ((double)totalOlapLength / numGapsClosed));

  fprintf(stderr, "\n");
  fprintf(stderr, "                  numGaps: %d\n", numGaps);
  fprintf(stderr, "        numExtendableGaps: %d (left: %d, right: %d, both: %d)\n", 
          numExtendableGaps, leftContigExtendable, rightContigExtendable, bothContigExtendable);
  fprintf(stderr, "            numGapsClosed: %d\n", numGapsClosed);
  fprintf(stderr, "             numSmallGaps: %d (closed %d, %.2f%%)\n", 
          numSmallGaps, numSmallGapsClosed, 100.0 * numSmallGapsClosed / numSmallGaps);
  fprintf(stderr, "             numLargeGaps: %d (closed %d, %.2f%%)\n", 
          numLargeGaps, numLargeGapsClosed, 100.0 * numLargeGapsClosed / numLargeGaps);
  fprintf(stderr, "                  allGaps: %d (closed %d, %.2f%%)\n", 
          numSmallGaps + numLargeGaps, numSmallGapsClosed + numLargeGapsClosed, 
          100.0 * (numSmallGapsClosed + numLargeGapsClosed) / (numSmallGaps + numLargeGaps));
  fprintf(stderr, "         maxGapSizeClosed: %.2f (gap number %d) \n", maxGapSizeClosed, maxGapSizeClosedNumber);
  fprintf(stderr, "       numGapsVarTooSmall: %d\n", numGapsVarTooSmall);
  fprintf(stderr, "        numClosingsTried : %d\n", numClosingsTried);
  fprintf(stderr, "unitigMultiAlignFailures : %d\n", unitigMultiAlignFailures);
  fprintf(stderr, "replaceEndUnitigFailures : %d\n", replaceEndUnitigFailures);
  fprintf(stderr, "   unitigToContigFailures: %d\n", unitigToContigFailures);
  fprintf(stderr, "    createAContigFailures: %d\n", createAContigFailures);
  fprintf(stderr, "           noOverlapFound: %d\n", noOverlapFound);

  if (numGapsClosed > 0)
    fprintf(stderr, "            avgOlapLength: %.2f\n", (float) totalOlapLength / numGapsClosed);
  fprintf(stderr, "         stddevOlapLength: %.2f\n", (float) sqrt(totalOlapVariance));
  if (numGapsClosed > 0)
    fprintf(stderr, "             avgOlapDiffs: %.2f\n", (float) totalOlapDiffs / numGapsClosed);
  fprintf(stderr, "          surrogatesOnEnd: %d\n", surrogatesOnEnd);
  fprintf(stderr, "\n");
#if 0
  //  oops, now private to eCR-examineGap.c
  fprintf(stderr, "                  MaxGaps: %d\n", MaxGaps);
  fprintf(stderr, "           MaxInteriorGap: %d\n", MaxInteriorGap);
#endif
  fprintf(stderr, "             CONTIG_BASES: %d\n", CONTIG_BASES);
  fprintf(stderr, "   totalContigsBaseChange: %d\n", totalContigsBaseChange);
  fprintf(stderr, "   totalBasesInClosedGaps: %d\n", totalBasesInClosedGaps);
  // scaffold base change should be the diff between originalsScaffolds.cgm and alteredScaffolds.cgm
  fprintf(stderr, "     scaffold base change: %d\n", totalContigsBaseChange - totalBasesInClosedGaps);
  fprintf(stderr, "\n");

  exit(0);
}




static
int
compExtendableFrags(const void *s1, const void *s2) {
  const extendableFrag * t1 = s1;
  const extendableFrag * t2 = s2;
  assert(t1 == s1);
  assert(t2 == s2);
  
  if (t1->extension > t2->extension)
    return -1;
  else if (t1->extension < t2->extension)
    return 1;
  else 
    return 0;
}



// findFirstFrag looks for a 3p->5p frag at the low end of a contig
// basesToNextFrag has meaning only when the first frag is the end frag, since basesToNextFrag
// marks where we start to have 2x coverage and is then used to determine MaxBegGap or MaxEndGap
int
findFirstExtendableFrags(ContigT *contig, extendableFrag *extFragsArray) {
  MultiAlignT *ma;
  IntMultiPos *mp;
  CIFragT *frag;
  int i, numFrags, firstUnitigID;
  int extendableFragCount = 0;
  int secondFragStart = 10000;

  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numFrags = GetNumIntMultiPoss(ma->f_list);
  firstUnitigID = GetIntUnitigPos(ma->u_list, 0)->ident;

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "in findFirstExtendableFrags, firstUnitigID: %d\n", firstUnitigID);
  
  for (i = 0; i < numFrags; i++) {
    mp = GetIntMultiPos(ma->f_list, i);
    frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) mp->sourceInt);

    if (frag->contigOffset3p.mean < 100.0 &&      // frag is within a cutoff of the low end of the contig
        frag->locale == -1 &&                     // and is a read
        frag->contigOffset3p.mean < frag->contigOffset5p.mean &&  // and points in the right direction
        frag->cid == firstUnitigID) { // and is in the first unitig
      unsigned int clr_bgn, clr_end, seq_len;
      int frag3pExtra, extension;
	  
      getFrag(ScaffoldGraph->gkpStore, frag->iid, fsread, FRAG_S_SEQ);

      switch (iterNumber) {
        case 1:
          clr_bgn = getFragRecordClearRegionBegin(fsread, AS_READ_CLEAR_ECR1 - 1);
          clr_end = getFragRecordClearRegionEnd  (fsread, AS_READ_CLEAR_ECR1 - 1);
          break;
        case 2:
          clr_bgn = getFragRecordClearRegionBegin(fsread, AS_READ_CLEAR_ECR2 - 1);
          clr_end = getFragRecordClearRegionEnd  (fsread, AS_READ_CLEAR_ECR2 - 1);
          break;
        default:
          assert(0);
          break;
      }

      seq_len = getFragRecordSequenceLength  (fsread);

      //                 <--------------------------------------------------------------------------- contig
      // 3p <------------------|---------------------------------------------|----- 5p frag
      //                    clr_end                                       clr_bgn
      //    |-----------|
      //      extension

      //             <--------------------------------------------------------------------------- contig
      //                     <-|---------------------------------------------|----- 5p frag
      //                    clr_end                                       clr_bgn
      //             |------|
      //               extension (negative)

      if (debug.eCRmainLV > 0) {
        fprintf(debug.eCRmainFP, "contig->bpLength.mean: %f\n", contig->bpLength.mean);
        fprintf(debug.eCRmainFP, "frag iid: %d, frag->contigOffset5p.mean: %f, frag->contigOffset3p.mean: %f\n",
                frag->iid, frag->contigOffset5p.mean, frag->contigOffset3p.mean);
        fprintf(debug.eCRmainFP, "frag length: " F_SIZE_T ", 3p past clr_end length: " F_SIZE_T "\n", seq_len, 
                seq_len - clr_end);
        fprintf(debug.eCRmainFP, "extension: " F_SIZE_T "\n", seq_len - clr_end - (int) frag->contigOffset3p.mean);
      }
	  
      frag3pExtra = seq_len - clr_end;
      extension = frag3pExtra - frag->contigOffset3p.mean;

      // ask Granger what min extension we should accept
      if (extension > 30) {
        extFragsArray[extendableFragCount].fragIid = frag->iid;
        extFragsArray[extendableFragCount].extension = extension;
        extFragsArray[extendableFragCount].addedBases = frag3pExtra;

        if (debug.eCRmainLV > 0)
          fprintf(debug.eCRmainFP, "for frag %d, extension: %8d, frag3pExtra: %8d\n",
                  frag->iid, extension, frag3pExtra);

        if (frag->contigOffset3p.mean == 0)
          extFragsArray[extendableFragCount].fragOnEnd = TRUE;
        else
          extFragsArray[extendableFragCount].fragOnEnd = FALSE;

        if (debug.eCRmainLV > 0) {
          fprintf(debug.eCRmainFP, "in contig %d, frag %d is at %f -> %f (5p->3p) \n", 
                  contig->id, frag->iid,
                  frag->contigOffset5p.mean, frag->contigOffset3p.mean);
          fprintf(debug.eCRmainFP, "extension ratio: %.2f\n", extension / (float) (1.0 + frag3pExtra - extension));
        }

        extendableFragCount++;
        if (extendableFragCount > MAX_EXTENDABLE_FRAGS) {
          fprintf(stderr, "extendableFragCount (%d) is greater than MAX_EXTENDABLE_FRAGS, aborting...\n",
                  extendableFragCount);
          assert(0);
        }
      }
    }

    // secondFragStart is where the next to end frag starts, and thus where we start 2x coverage
    // we don't care if it's the 3p or 5p end
    if (frag->contigOffset3p.mean > 0 && frag->contigOffset3p.mean < secondFragStart) {
      secondFragStart = frag->contigOffset3p.mean;
      // fprintf(stderr, "secondFragStart %d set by frag %d (3p)\n", secondFragStart, frag->iid);
    }
    if (frag->contigOffset5p.mean < secondFragStart) {
      secondFragStart = frag->contigOffset5p.mean;
      // fprintf(stderr, "secondFragStart %d set by frag %d (5p)\n", secondFragStart, frag->iid);
    }
  }

  // now sort the extendable frags by their extendability
  qsort(extFragsArray, extendableFragCount, sizeof(extendableFrag), compExtendableFrags);

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "extendableFragCount: %d\n", extendableFragCount);

  for (i = 0; i < extendableFragCount; i++) {
    if (extFragsArray[i].fragOnEnd == TRUE)
      extFragsArray[i].basesToNextFrag = secondFragStart;
    else
      extFragsArray[i].basesToNextFrag = 0;

    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "contig %8d, frag %8d can extend %8d bases into the gap\n",
              contig->id, extFragsArray[i].fragIid, extFragsArray[i].extension);
  }

  return extendableFragCount;
}



// findLastFrag looks for a 5p->3p frag at the high end of a contig
// basesToNextFrag has meaning only when the first frag is the end frag, since basesToNextFrag
// marks where we start to have 2x coverage and is then used to determine MaxBegGap or MaxEndGap
int
findLastExtendableFrags(ContigT *contig, extendableFrag *extFragsArray) {
  MultiAlignT *ma;
  IntMultiPos *mp;
  CIFragT *frag;
  int i, numFrags, lastUnitigID;
  double maxContigPos;
  int secondFragEnd = 0;
  int extendableFragCount = 0;
  
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numFrags = GetNumIntMultiPoss(ma->f_list);
  lastUnitigID = GetIntUnitigPos(ma->u_list, GetNumIntUnitigPoss(ma->u_list) - 1)->ident;  
  maxContigPos = contig->bpLength.mean - 1.0;

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "in FindLastExtendableFrags, lastUnitigID: %d\n", lastUnitigID);
  
  for (i = 0; i < numFrags; i++) {
    mp = GetIntMultiPos(ma->f_list, i);
    frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) mp->sourceInt);

    if (frag->contigOffset3p.mean > maxContigPos - 100.0 &&      // frag is within a cutoff of the high end of the contig
        frag->locale == -1 &&                                    // and is a read
        frag->contigOffset5p.mean < frag->contigOffset3p.mean && // and points in the right direction
        frag->cid == lastUnitigID) {                             // and is in the last unitig
      unsigned int clr_bgn, clr_end, seq_len;
      int frag3pExtra, extension;

      getFrag(ScaffoldGraph->gkpStore, frag->iid, fsread, FRAG_S_SEQ);

      switch (iterNumber) {
        case 1:
          clr_bgn = getFragRecordClearRegionBegin(fsread, AS_READ_CLEAR_ECR1 - 1);
          clr_end = getFragRecordClearRegionEnd  (fsread, AS_READ_CLEAR_ECR1 - 1);
          break;
        case 2:
          clr_bgn = getFragRecordClearRegionBegin(fsread, AS_READ_CLEAR_ECR2 - 1);
          clr_end = getFragRecordClearRegionEnd  (fsread, AS_READ_CLEAR_ECR2 - 1);
          break;
        default:
          assert(0);
          break;
      }

      seq_len = getFragRecordSequenceLength  (fsread);

      //    contig ----------------------------------------------------------------------------------->
      //                                  5p -------|---------------------------------------------|------------> 3p 
      //                                         clr_bgn                                       clr_end
      //                                                                                               |-------|
      //                                                                                             extension
      frag3pExtra = seq_len - clr_end;
      extension = frag3pExtra - (int) (contig->bpLength.mean - frag->contigOffset3p.mean);

      if (debug.eCRmainLV > 0) {
        fprintf(debug.eCRmainFP, "contig->bpLength.mean: %f\n", contig->bpLength.mean);
        fprintf(debug.eCRmainFP, "frag iid: %d, frag->contigOffset5p.mean: %f, frag->contigOffset3p.mean: %f\n",
                frag->iid, frag->contigOffset5p.mean, frag->contigOffset3p.mean);
        fprintf(debug.eCRmainFP, "frag length: " F_SIZE_T ", 3p past clr_end length: %d\n", seq_len, frag3pExtra);
        fprintf(debug.eCRmainFP, "extension: %d\n", extension);
      }
	  
      if (extension > 30) {
        extFragsArray[extendableFragCount].fragIid = frag->iid;
        extFragsArray[extendableFragCount].extension = extension;
        extFragsArray[extendableFragCount].addedBases = frag3pExtra;

        if (debug.eCRmainLV > 0)
          fprintf(debug.eCRmainFP, "for frag %d, extension: %8d, frag3pExtra: %8d\n",
                  frag->iid, extension, frag3pExtra);

        if (frag->contigOffset3p.mean == contig->bpLength.mean)
          extFragsArray[extendableFragCount].fragOnEnd = TRUE;
        else
          extFragsArray[extendableFragCount].fragOnEnd = FALSE;

        if (debug.eCRmainLV > 0) {
          fprintf(debug.eCRmainFP, "in contig %d, frag %d is at %f -> %f (5p->3p) maxContigPos: %f\n", 
                  contig->id, frag->iid,
                  frag->contigOffset5p.mean, frag->contigOffset3p.mean, maxContigPos);
          fprintf(debug.eCRmainFP, "extension ratio: %.2f\n", extension / (float) (1.0 + frag3pExtra - extension));
        }

        extendableFragCount++;
        if (extendableFragCount > MAX_EXTENDABLE_FRAGS) {
          fprintf(stderr, "extendableFragCount (%d) is greater than MAX_EXTENDABLE_FRAGS, aborting...\n",
                  extendableFragCount);
          assert(0);
        }
      }
    }

    // secondFragEnd is where the next to end frag ends, and thus where we have 2x coverage
    // we don't care if it's the 3p or 5p end
    if (frag->contigOffset3p.mean < contig->bpLength.mean && (int) frag->contigOffset3p.mean > secondFragEnd) {
      secondFragEnd = (int) frag->contigOffset3p.mean;
      // fprintf(stderr, "secondFragEnd %d set by frag %d (3p)\n", secondFragEnd, frag->iid);
    }
    if ((int) frag->contigOffset5p.mean > secondFragEnd) {
      secondFragEnd = frag->contigOffset5p.mean;
      // fprintf(stderr, "secondFragEnd %d set by frag %d (5p)\n", secondFragEnd, frag->iid);
    }
  }

  // now sort the extendable frags by their extendability
  qsort(extFragsArray, extendableFragCount, sizeof(extendableFrag), compExtendableFrags);

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "extendableFragCount: %d\n", extendableFragCount);

  for (i = 0; i < extendableFragCount; i++) {
    if (extFragsArray[i].fragOnEnd == TRUE)
      extFragsArray[i].basesToNextFrag = maxContigPos - secondFragEnd;
    else
      extFragsArray[i].basesToNextFrag = 0;	  
	
    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "contig %8d, frag %8d can extend %8d bases into the gap\n",
              contig->id, extFragsArray[i].fragIid, extFragsArray[i].extension);
  }

  return extendableFragCount;
}

// findLastUnitig looks for a surrogate at the high end of a contig
int
findLastUnitig(ContigT *contig, int *unitigID) {
  MultiAlignT *ma;
  int i, numUnitigs, isSurrogate = FALSE;
  double maxContigPos = 0.0;
  NodeCGW_T *unitig = NULL;
  
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numUnitigs = GetNumIntUnitigPoss(ma->u_list);
  
  // can't just jump to last unitig since unitigs are arranged by starting position, not ending
  for (i = 0; i < numUnitigs; i++) {
    IntUnitigPos *upos = GetIntUnitigPos(ma->u_list, i);
    unitig = GetGraphNode(ScaffoldGraph->CIGraph, upos->ident);
	
    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "findLastUnitig()-- in contig %d, unitig %d is at %f -> %f maxContigPos: %f, isSurrogate: %d, baseID: %d\n", 
              contig->id, unitig->id,
              unitig->offsetAEnd.mean, unitig->offsetBEnd.mean, maxContigPos,
              isSurrogate, unitig->info.CI.baseID);

    if (unitig->offsetAEnd.mean >= maxContigPos || unitig->offsetBEnd.mean >= maxContigPos) {
      maxContigPos = MAX(unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);
      *unitigID = unitig->id;
      isSurrogate = unitig->flags.bits.isSurrogate;
    }
  }
  
  if (isSurrogate) {
    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "findLastUnitig()-- unitig %d on high end of contig %d is surrogate!\n",
              unitig->id, contig->id);
    return 1;
  }
  return 0;
}

// findFirstUnitig looks for a surrogate at the low end of a contig
int findFirstUnitig(ContigT *contig, int *unitigID) {
  MultiAlignT *ma;
  int numUnitigs;
  IntUnitigPos *upos;
  NodeCGW_T *unitig;
  int isSurrogate;

  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numUnitigs = GetNumIntUnitigPoss(ma->u_list);

  upos = GetIntUnitigPos(ma->u_list, 0);
  unitig = GetGraphNode(ScaffoldGraph->CIGraph, upos->ident);
  isSurrogate = unitig->flags.bits.isSurrogate;

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "findFirstUnitig()-- in contig %d, unitig %d is at %f -> %f, isSurrogate: %d\n", 
            contig->id, unitig->id,
            unitig->offsetAEnd.mean, unitig->offsetBEnd.mean,
            isSurrogate);

  *unitigID = unitig->id;

  if (isSurrogate) {
    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "findFirstUnitig()-- unitig %d on low end of contig %d is surrogate!\n",
              unitig->id, contig->id);
    return 1;
  }
  return 0;
}











void extendContig(ContigT *contig, int extendAEnd) {
  // have to alter the following fields in a NodeCGW_T: bpLength, offsetAEnd, offsetBEnd
  int contigOrientation = contig->offsetAEnd.mean < contig->offsetBEnd.mean ? A_B : B_A;
  MultiAlignT *cma = NULL;
  double lengthDelta;
  
  cma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);

  lengthDelta = strlen(Getchar(cma->consensus, 0)) - contig->bpLength.mean;

  if (debug.eCRmainLV > 0) {
    if (contigOrientation == A_B)
      fprintf(debug.eCRmainFP, "extendContig()-- contig %8d original pos: %.0f, %.0f A_B\n",
              contig->id, contig->offsetAEnd.mean, contig->offsetBEnd.mean);
    else
      fprintf(debug.eCRmainFP, "extendContig()-- contig %8d original pos: %.0f, %.0f B_A\n",
              contig->id, contig->offsetBEnd.mean, contig->offsetAEnd.mean);
  }
  
  contig->bpLength.mean += lengthDelta;
  
  if (contigOrientation == A_B) {
    if (extendAEnd == TRUE)
      contig->offsetAEnd.mean -= lengthDelta;
    else  // extend the B end
      contig->offsetBEnd.mean += lengthDelta;
  } else {
    // contigOrientation is B_A
    if (extendAEnd == TRUE)
      contig->offsetAEnd.mean += lengthDelta;
    else  // extend the B end
      contig->offsetBEnd.mean -= lengthDelta;
  }

  if (debug.eCRmainLV > 0) {
    if (contigOrientation == A_B)
      fprintf(debug.eCRmainFP, "extendContig()-- contig %8d altered pos: %.0f, %.0f A_B\n",
              contig->id, contig->offsetAEnd.mean, contig->offsetBEnd.mean);
    else
      fprintf(debug.eCRmainFP, "extendContig()-- contig %8d altered pos: %.0f, %.0f B_A\n",
              contig->id, contig->offsetBEnd.mean, contig->offsetAEnd.mean);  
  }
}



// stole most of this from OutputUnitigsFromMultiAligns
/********************************************************************************/
int GetNewUnitigMultiAlign(NodeCGW_T *unitig, fragPositions *fragPoss, int extendedFragIid) {
  GenericMesg			pmesg;
  IntUnitigMesg			ium_mesg;
  int i;
  int numCIs = GetNumGraphNodes(ScaffoldGraph->CIGraph);
  MultiAlignT *macopy = NULL;
  UnitigStatus   status;

#ifdef USE_UNGAPPED_CONSENSUS_FOR_UNITIG
  static VA_TYPE(char) *ungappedSequence=NULL,*ungappedQuality=NULL;
  if (ungappedSequence== NULL) {
    ungappedSequence = CreateVA_char(0);
    ungappedQuality = CreateVA_char(0);
  } else {
    ResetVA_char(ungappedSequence);
    ResetVA_char(ungappedQuality);
  }
#endif

  pmesg.m = &ium_mesg;
  pmesg.t = MESG_IUM;

  assert(unitig->id>=0 && unitig->id< numCIs);

  switch(unitig->type) {
    case DISCRIMINATORUNIQUECHUNK_CGW:
      status = AS_UNIQUE;
      //fprintf(GlobalData->stderrc,"* Unitig %d is UNIQUE: DISCRIMINATOR output %d \n",unitig->id, cid);
      break;
    case UNIQUECHUNK_CGW:
      status = AS_UNIQUE;
      //fprintf(GlobalData->stderrc,"* Unitig %d is UNIQUE: output %d \n",unitig->id, cid);
      break;
    case UNRESOLVEDCHUNK_CGW:
      if (unitig->info.CI.numInstances > 0) {
        //fprintf(GlobalData->stderrc, "* Unitig %d has %d instances--- output %d SEP\n",unitig->id, unitig->info.CI.numInstances,cid);
        status = AS_SEP;
        assert(!unitig->flags.bits.isUnique);
      } else if (unitig->scaffoldID != NULLINDEX) {
        //fprintf(GlobalData->stderrc, "* Unitig %d has %d instances--- output %d UNIQUE\n",unitig->id, unitig->info.CI.numInstances,cid);
        status = AS_UNIQUE;
      } else {
        //fprintf(GlobalData->stderrc, "* Unitig %d has %d instances--- output %d NOTREZ\n",unitig->id, unitig->info.CI.numInstances,cid);
        status = AS_NOTREZ;
      }
      break;
    case RESOLVEDREPEATCHUNK_CGW:
      /* SKIP THESE */
      //fprintf(GlobalData->stderrc,"* Skipping unitig %d --- RESOLVEDREPEAT\n",unitig->id);
      //continue;
    default:
      assert(0);
  }

  //  ReLoad...() makes a copy of the multialign in the store.  This is needed
  //  because at the end of this function, we delete the multialign in the store,
  //  and insert a new one.  The new one is based on the copy.
  //  
  macopy = CreateEmptyMultiAlignT();
  ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, macopy, unitig->id, TRUE);

  if (debug.eCRmainLV > 0) {
    fprintf(debug.eCRmainFP, "for unitig %d, before reforming, strlen(macopy->consensus) = " F_SIZE_T "\n", unitig->id, strlen(Getchar(macopy->consensus, 0)));
    fprintf(debug.eCRmainFP, "for unitig %d, before reforming, consensus:\n%s\n", unitig->id, Getchar(macopy->consensus, 0));
  }

  {
    int extendedFragLeftward, aligned;
    assert (unitig->type != CONTIG_CGW);

    ium_mesg.iaccession     = unitig->id;
#ifdef AS_ENABLE_SOURCE
    ium_mesg.source         = NULL;
#endif
    ium_mesg.coverage_stat  = unitig->info.CI.coverageStat;
    ium_mesg.status         = status;

#ifdef USE_UNGAPPED_CONSENSUS_FOR_UNITIG
    GetMultiAlignUngappedConsensus(macopy, ungappedSequence, ungappedQuality);
    ium_mesg.consensus = Getchar(ungappedSequence,0);
    ium_mesg.quality = Getchar(ungappedQuality,0);
    ium_mesg.length = GetMultiAlignUngappedLength(macopy);
#else
    ium_mesg.length = GetMultiAlignLength(macopy);
    ium_mesg.consensus = Getchar(macopy->consensus,0);
    ium_mesg.quality = Getchar(macopy->quality,0);
#endif
    ium_mesg.forced = 0;
    ium_mesg.num_frags = GetNumIntMultiPoss(macopy->f_list);
    ium_mesg.num_vars = 0;

    // replace the positions in the f_list with the adjusted positions
    extendedFragLeftward = FALSE;

    for (i = 0; i < GetNumIntMultiPoss(macopy->f_list); i++) {
      IntMultiPos *tempPos = GetIntMultiPos(macopy->f_list, i);

#if 0
      fprintf(stderr, "in GetNewUnitigMultiAlign, changing frag %10d from (%8d, %8d) to (%8d, %8d)\n",
              tempPos->ident, tempPos->position.bgn, tempPos->position.end, fragPoss[i].bgn, fragPoss[i].end);
#endif

      if (tempPos->ident == extendedFragIid) {
        tempPos->contained = 0;    // the extended frag is no longer contained by any other frag
        if (fragPoss[i].bgn == 0 || fragPoss[i].end == 0)
          extendedFragLeftward = TRUE;  // if at the beginning of the unitig have to reorder frags
      }
      tempPos->position.bgn = fragPoss[i].bgn;
      tempPos->position.end = fragPoss[i].end;
    }
    ium_mesg.f_list = GetIntMultiPos(macopy->f_list, 0);

    if (extendedFragLeftward)
      // definitely reorder frags in f_list if we extended a frag leftward
      leftShiftIUM(ium_mesg.f_list, GetNumIntMultiPoss(macopy->f_list), extendedFragIid);
    else
      // might need to reorder frags in f_list if we extended a frag rightward
      rightShiftIUM(ium_mesg.f_list, GetNumIntMultiPoss(macopy->f_list), extendedFragIid);

    //  The Local_Overlap_AS_forCNS below also had a commented out
    //  DP_Compare.
    //
    {
      CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                              CNS_OPTIONS_SMOOTH_WIN_DEFAULT,
                              CNS_OPTIONS_MAX_NUM_ALLELES };
      ALIGNMENT_CONTEXT=AS_CONSENSUS;


      cnslog = stderr;
      aligned = MultiAlignUnitig(&ium_mesg, 
                                 ScaffoldGraph->gkpStore,
                                 reformed_consensus,
                                 reformed_quality,
                                 reformed_deltas,
                                 CNS_QUIET,  //  CNS_VERBOSE, CNS_STATS_ONLY
                                 1,
                                 Local_Overlap_AS_forCNS,
                                 &options);
    }

    if (aligned == EXIT_FAILURE) {
      fprintf(stderr, "GetNewUnitigMultiAlign()-- MultiAlignUnitig failure on unitig %d\n", unitig->id);
      return FALSE;
      //assert(0);
    }

    if (debug.eCRmainLV > 0) {
      fprintf(debug.eCRmainFP, "for unitig %d, after reforming, strlen(reformed_consensus) = " F_SIZE_T "\n", unitig->id, strlen(Getchar(reformed_consensus, 0)));
      fprintf(debug.eCRmainFP, "for unitig %d, after reforming, consensus:\n%s\n", unitig->id, Getchar(reformed_consensus, 0));
    }

    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg.iaccession, TRUE);

    //  Set the keepInCache flag, since CreateMultiAlignTFromIUM() is
    //  returning an allocated multialign, and if we don't cache it,
    //  we have a giant memory leak.
    //
    InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg.iaccession, TRUE,
                                  CreateMultiAlignTFromIUM(&ium_mesg, -2, FALSE),
                                  TRUE);
  }

  DeleteMultiAlignT(macopy);

  return TRUE;
}




void
extendClearRange(int fragIid, int frag3pDelta) {
  unsigned int clr_bgn, clr_end;

  if (frag3pDelta < 0)
    fprintf(stderr, "extendClearRange()-- WARNING: frag3pDelta less than zero: %d\n", frag3pDelta);

  if (fragIid != -1) {
    getFrag(ScaffoldGraph->gkpStore, fragIid, fsread, FRAG_S_INF);

    switch (iterNumber) {
      case 1:
        clr_bgn = getFragRecordClearRegionBegin(fsread, AS_READ_CLEAR_ECR1 - 1);
        clr_end = getFragRecordClearRegionEnd  (fsread, AS_READ_CLEAR_ECR1 - 1) + frag3pDelta;
        setFragRecordClearRegion(fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR1);
        break;
      case 2:
        clr_bgn = getFragRecordClearRegionBegin(fsread, AS_READ_CLEAR_ECR2 - 1);
        clr_end = getFragRecordClearRegionEnd  (fsread, AS_READ_CLEAR_ECR2 - 1) + frag3pDelta;
        setFragRecordClearRegion(fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR2);
        break;
      default:
        assert(0);
        break;
    }

    setFrag(ScaffoldGraph->gkpStore, fragIid, fsread);

    fprintf(stderr, "extendClearRange()-- changed frag %d clr_end from %d to %d\n",
            fragIid, clr_end, clr_end + frag3pDelta);
  }
}


void
revertClearRange(int fragIid) {
  unsigned int clr_bgn, clr_end;
  
  if (fragIid != -1) {
    getFrag(ScaffoldGraph->gkpStore, fragIid, fsread, FRAG_S_INF);

    switch (iterNumber) {
      case 1:
        clr_bgn = getFragRecordClearRegionBegin(fsread, AS_READ_CLEAR_ECR1 - 1);
        clr_end = getFragRecordClearRegionEnd  (fsread, AS_READ_CLEAR_ECR1 - 1);
        setFragRecordClearRegion(fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR1);
        break;
      case 2:
        clr_bgn = getFragRecordClearRegionBegin(fsread, AS_READ_CLEAR_ECR2 - 1);
        clr_end = getFragRecordClearRegionEnd  (fsread, AS_READ_CLEAR_ECR2 - 1);
        setFragRecordClearRegion(fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR2);
        break;
      default:
        assert(0);
        break;
    }

    setFrag(ScaffoldGraph->gkpStore, fragIid, fsread);
  }
}





void
getAlteredFragPositions(NodeCGW_T *unitig, fragPositions **fragPoss, int alteredFragIid, int extension) {
  fragPositions *localFragPoss;
  int i, alteredFragIndex = NULLINDEX, orientation = 0, delta;
  MultiAlignT *uma = NULL;

  // currently code does not handle negative extensions, ie, trimming
  if (extension <= 0)
    fprintf(stderr, "getAlteredFragPositions()-- negative extension: %d\n", extension);

  uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);

  localFragPoss = (fragPositions *) safe_malloc( GetNumIntMultiPoss(uma->f_list) * sizeof(fragPositions));

  // get the current positions
  for (i = 0; i < GetNumIntMultiPoss(uma->f_list); i++) {
    IntMultiPos *pos = GetIntMultiPos(uma->f_list, i);
  
    localFragPoss[i].bgn = pos->position.bgn;
    localFragPoss[i].end = pos->position.end;

    if (pos->ident == alteredFragIid) {
      alteredFragIndex = i;
 
      if (pos->position.bgn < pos->position.end)
        orientation = 0;
      else
        orientation = 1;
    }

    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "getAlteredFragPositions()-- %2d] (%d,%d)\n", i, localFragPoss[i].bgn, localFragPoss[i].end);
  }

  // alteredFragDelta is the amount by which the altered frag gets extended

  if (orientation == 0) {
    // nothing much to do, just extend the altered fragment

    // first find out how much must be added to the frag to get to the end of the untig
    // int alteredFragDelta = unitig->bpLength.mean - localFragPoss[alteredFragIndex].end;

    int alteredFragDelta = strlen(Getchar(uma->consensus, 0)) - localFragPoss[alteredFragIndex].end;

    // then add how far it extends out into the gap
    alteredFragDelta += extension;

    localFragPoss[alteredFragIndex].end += alteredFragDelta;

    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "getAlteredFragPositions()-- adjust end fragment %d by %d\n", alteredFragIndex, alteredFragDelta);
  } else {
    // all frag positions get bumped by extension, that's how far the altered frag extends into the gap

    // this is the new minimum position in the unitig

    localFragPoss[alteredFragIndex].end = -extension;

    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "getAlteredFragPositions()-- adjust first fragment %d to (%d,%d)\n",
              alteredFragIndex, localFragPoss[alteredFragIndex].bgn, localFragPoss[alteredFragIndex].end);

    // if he extends off the front of the unitig, adjust everybody upward

    delta = localFragPoss[alteredFragIndex].end;

    if (delta < 0) {
      for (i = 0; i < GetNumIntMultiPoss(uma->f_list); i++) {
        localFragPoss[i].bgn -= delta;
        localFragPoss[i].end -= delta;
        if (debug.eCRmainLV > 0)
          fprintf(debug.eCRmainFP, "getAlteredFragPositions()-- %2d] (%d,%d) (adjusted by -delta = %d)\n", i, localFragPoss[i].bgn, localFragPoss[i].end, -delta);
      }
    }
  }
  *fragPoss = localFragPoss;
}





// this routine shifts the frag with iid extendedFragIid to the front of the f_list
void
leftShiftIUM(IntMultiPos *f_list, int numFrags, int extendedFragIid) {
  int i, currPos = 0, numShiftedInUnitig = 0;
  IntMultiPos tempIMP;
  
  // find out the extended frag's current position
  for (i = 0; i < numFrags; i++)
    if (f_list[i].ident == extendedFragIid)
      currPos = i;
  
  // save the extended frag off to the side
  memcpy(&tempIMP, &f_list[currPos], sizeof(IntMultiPos));

  // now shift everybody below currPos up a position
  for (i = currPos - 1; i >= 0; i--) {
    memcpy(&f_list[i+1], &f_list[i], sizeof(IntMultiPos));
    numShiftedInUnitig++;
  }
  
  // and place tempPos in position 0 of the unitig
  memcpy(&f_list[0], &tempIMP, sizeof(IntMultiPos));
  
#if 0
  // now do the same shifting for the contig
  {
    NodeCGW_T *contig = GetGraphNode(ScaffoldGraph->ContigGraph, contigID);
    MultiAlignT *cma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contigID, FALSE);
    IntMultiPos *cf_list = cma->f_list;
    int icnt, numContigFrags = GetNumIntMultiPoss(cma->f_list);
	
    // find out the extended frag's current position
    for (i = 0; i < numContigFrags; i++)
      if (cf_list[i].ident == extendedFragIid)
        currPos = i;
	
    // save the extended frag off to the side
    memcpy(&tempIMP, &cf_list[currPos], sizeof(IntMultiPos));
	
    // now shift the appropriate number of frags up a position starting at currPos
    for (i = currPos - 1, icnt = 0; icnt < numShiftedInUnitig; i--, icnt++)
      {
        memcpy(&cf_list[i+1], &cf_list[i], sizeof(IntMultiPos));
        numShiftedInUnitig++;
      }
	
    // and place tempPos in position 0
    memcpy(&cf_list[currPos - numShiftedInUnitig], &tempIMP, sizeof(IntMultiPos));
  }
#endif
}

// this routine shifts the frag with iid extendedFragIid to the appropriate place in the f_list
void
rightShiftIUM(IntMultiPos *f_list, int numFrags, int extendedFragIid) {
  int i, currPos = 0;
  IntMultiPos tempIMP;
  int numShifted, shiftedFrag, fragToMovePos = NULLINDEX, numToShift;

  // find out the extended frag's current position
  for (i = 0; i < numFrags; i++)
    if (f_list[i].ident == extendedFragIid)
      currPos = i;  

  // we need to move all frags that 
  // i) start to the left of the extended frag and 
  // ii) are not contained by another frag 
  // before the extended frag in the array.  We could move the contained ones, but don't have to

  numShifted = 0;
  shiftedFrag = TRUE;
  while (shiftedFrag) {
    shiftedFrag = FALSE;
	
    // look for a candidate frag
    i = currPos + numShifted + 1;
    while (i < numFrags) {
      if ((MIN(f_list[i].position.bgn, f_list[i].position.end) < 
           MIN(f_list[currPos + numShifted].position.bgn, f_list[currPos + numShifted].position.end)) &&
          f_list[i].contained == FALSE) {
        fragToMovePos = i;
        shiftedFrag = TRUE;
        break;
      }
      i++;
    }

    if (shiftedFrag) {
      // save the extended frag to be shifted off to the side
      memcpy(&tempIMP, &f_list[fragToMovePos], sizeof(IntMultiPos));
	  
      // now shift all the frags between currPos + numShifted up a position in the array
      numToShift = fragToMovePos - (currPos + numShifted);
      i = 0;
      while (i < numToShift) {
        memcpy(&f_list[fragToMovePos - i], &f_list[fragToMovePos - i - 1], sizeof(IntMultiPos));
        i++;
      }

      // and tempPos into open position
      memcpy(&f_list[currPos + numShifted], &tempIMP, sizeof(IntMultiPos));

      numShifted++;
    }
  }
}


void
saveFragAndUnitigData(int lFragIid, int rFragIid) {
  CIFragT *leftFrag, *rightFrag;
  InfoByIID *info;
  
  if (savedLeftUnitigMA == NULL) {
    savedLeftUnitigMA  = CreateEmptyMultiAlignT();
    savedRightUnitigMA = CreateEmptyMultiAlignT();
    savedLeftContigMA  = CreateEmptyMultiAlignT();
    savedRightContigMA = CreateEmptyMultiAlignT();
  }

  // grab the multialignments
  if (lFragIid != -1) {
    info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);
    assert(info->set);
    leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

    ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, savedLeftUnitigMA, leftFrag->cid, TRUE);
    ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, savedLeftContigMA, leftFrag->contigID, FALSE);
  }
  
  if (rFragIid != -1) {
    info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
    assert(info->set);
    rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

    ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, savedRightUnitigMA, rightFrag->cid, TRUE);
    ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, savedRightContigMA, rightFrag->contigID, FALSE);
  }
}

void
restoreFragAndUnitigData(int lFragIid, int rFragIid) {
  NodeCGW_T *unitig;
  CIFragT *leftFrag, *rightFrag;

  //  We do not own the saved*MA's that are inserted into the
  //  SequenceDB; do not set the keepInCache flag to
  //  InsertMultiAlignTInSequenceDB().

  // first the left frag
  if (lFragIid != -1) {
    InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);

    assert(info->set);
    leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
    unitig = GetGraphNode(ScaffoldGraph->CIGraph, leftFrag->cid);
    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);
    InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE, savedLeftUnitigMA, FALSE);
    SynchUnitigTWithMultiAlignT(unitig);

    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, leftFrag->contigID, FALSE);
    InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, leftFrag->contigID, FALSE, savedLeftContigMA, FALSE);

    if (debug.diagnosticLV > 1) {
      DumpContigMultiAlignInfo ("restoreFragAndUnitigData()--left", NULL, leftFrag->contigID);
      DumpContigUngappedOffsets("restoreFragAndUnitigData()--left", leftFrag->contigID);
    }
  }
  
  // now the right frag
  if (rFragIid != -1) {
    InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);

    assert(info->set);
    rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
    unitig = GetGraphNode(ScaffoldGraph->CIGraph, rightFrag->cid);
    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);
    InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE, savedRightUnitigMA, FALSE);
    SynchUnitigTWithMultiAlignT(unitig);
	
    UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, rightFrag->contigID, FALSE);
    InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, rightFrag->contigID, FALSE, savedRightContigMA, FALSE);

    if (debug.diagnosticLV > 1) {
      DumpContigMultiAlignInfo ("restoreFragAndUnitigData()--right", NULL, rightFrag->contigID);
      DumpContigUngappedOffsets("restoreFragAndUnitigData()--right", rightFrag->contigID);
    }
  }
}
