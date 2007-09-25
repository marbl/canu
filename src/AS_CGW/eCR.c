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

static const char CM_ID[] = "$Id: eCR.c,v 1.31 2007-09-25 05:46:28 brianwalenz Exp $";

#include "eCR.h"
#include "ScaffoldGraph_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "MultiAlignment_CNS.h"
#include "GapWalkerREZ.h"  //  FindGapLength

#define MAX_EXTENDABLE_FRAGS   100
#define NUM_STDDEV_CUTOFF        5.0

typedef struct {
  int fragIid;
  int ctgMaxExt;
  int frgMaxExt;
  int basesToNextFrag;
  int fragOnEnd;
} extendableFrag;

typedef struct {
  int bgn;
  int end;
} fragPositions;


int findFirstExtendableFrags(ContigT *contig, extendableFrag *extFragsArray);
int findLastExtendableFrags(ContigT *contig, extendableFrag *extFragsArray);

int findFirstUnitig(ContigT *contig, int *unitigID);
int findLastUnitig(ContigT *contig, int *unitigID);

void extendContig(ContigT *contig, int extendAEnd);

int GetNewUnitigMultiAlign(NodeCGW_T *unitig,
                           fragPositions *fragPoss,
                           int extFragIID,
                           int extLength);

void extendClearRange(int fragIid, int frag3pDelta);
void revertClearRange(int fragIid);

fragPositions *getAlteredFragPositions(NodeCGW_T *unitig,
                                       int alteredFragIid,
                                       int extension);

void leftShiftIUM(IntMultiPos *f_list, int numFrags, int extFragIID);
void rightShiftIUM(IntMultiPos *f_list, int numFrags, int extFragIID);

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




void
SynchUnitigTWithMultiAlignT(NodeCGW_T *unitig) {
  unitig->bpLength.mean = GetMultiAlignUngappedLength(loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
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

  argc = AS_configure(argc, argv);

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
  if (scaffoldBegin == scaffoldEnd) {
    fprintf(stderr, "%s: Nothing to do, no scaffolds selected: -b == -e\n", argv[0]);
    exit(0);
  }

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

#if 1
    //  Reset all the fragments in the scaffold to the correct clear
    //  range.  This only does something if eCR was restarted, if eCR
    //  had already processed a scaffold in out set.  If it's a fresh
    //  run, the clear ranges will all be 'correct' and this does
    //  nothing.
    //
    //  If we don't do this, the stuff in AS_CNS gets terribly confused.
    //
    {
      int        contigID = scaff->info.Scaffold.AEndCI;
      fragRecord fsread;

      while (contigID != -1) {
        ContigT     *contig = GetGraphNode(ScaffoldGraph->ContigGraph, contigID);
        MultiAlignT *ma     = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
        int          i;

        for (i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
          CDS_IID_t       iid = GetCIFragT(ScaffoldGraph->CIFrags,
                                           GetIntMultiPos(ma->f_list, i)->sourceInt)->iid;
          unsigned int    clr_bgn;
          unsigned int    clr_end;

          getFrag(ScaffoldGraph->gkpStore, iid, &fsread, FRAG_S_INF);

          switch (iterNumber) {
            case 1:
              clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR1 - 1);
              clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR1 - 1);
              if ((clr_bgn != getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR1)) ||
                  (clr_end != getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR1))) {
                fprintf(stderr, "WARNING:  Reset clear for frag "F_IID","F_UID" (ctg %d scf %d) from (%d,%d) to (%d,%d)\n",
                        getFragRecordIID(&fsread),
                        getFragRecordUID(&fsread),
                        contigID, sid,
                        getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR1),
                        getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR1),
                        clr_bgn,
                        clr_end);
                setFragRecordClearRegion(&fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR1);
                setFrag(ScaffoldGraph->gkpStore, iid, &fsread);
              }
              break;
            case 2:
              clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR2 - 1);
              clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR2 - 1);
              if ((clr_bgn != getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR2)) ||
                  (clr_end != getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR2))) {
                fprintf(stderr, "WARNING:  Reset clear for frag "F_IID","F_UID" (ctg %f scf %d) from (%d,%d) to (%d,%d)\n",
                        getFragRecordIID(&fsread),
                        getFragRecordUID(&fsread),
                        contigID, sid,
                        getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR2),
                        getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR2),
                        clr_bgn,
                        clr_end);
                setFragRecordClearRegion(&fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR2);
                setFrag(ScaffoldGraph->gkpStore, iid, &fsread);
              }
              break;
            default:
              assert(0);
              break;
          }
        }

        contigID = contig->BEndNext;
      }
    }
#endif

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

      ContigT    lcontigBackup = {0};
      ContigT    rcontigBackup = {0};

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

      if (debug.eCRmainLV > 0) {
        fprintf(debug.eCRmainFP, "finished examining lcontig %d (orientation %c) \n", lcontig->id, lcontigOrientation);
        fprintf(debug.eCRmainFP, "\nexamining rcontig %d (orientation %c) \n", rcontig->id, rcontigOrientation);
      }

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

		
      //  Add an extra "fragment" to tell us to also try the contig sequence itself.
      //
      leftExtFragsArray[numLeftFrags].fragIid           = -1;
      leftExtFragsArray[numLeftFrags].ctgMaxExt         = 0;
      leftExtFragsArray[numLeftFrags].frgMaxExt         = 0;
      leftExtFragsArray[numLeftFrags].basesToNextFrag   = 0;
      leftExtFragsArray[numLeftFrags].fragOnEnd         = FALSE;

      numLeftFrags++;

      rightExtFragsArray[numRightFrags].fragIid         = -1;		  
      rightExtFragsArray[numRightFrags].ctgMaxExt       = 0;
      rightExtFragsArray[numRightFrags].frgMaxExt       = 0;
      rightExtFragsArray[numRightFrags].basesToNextFrag = 0;
      rightExtFragsArray[numRightFrags].fragOnEnd       = FALSE;

      numRightFrags++;


      //  In case we need to back out changes late in the extension, we save copies of the two contigs.
      //
      memcpy(&lcontigBackup, lcontig, sizeof(ContigT));
      memcpy(&rcontigBackup, rcontig, sizeof(ContigT));

			
      for (leftFragIndex = 0; leftFragIndex < numLeftFrags && closedGap == FALSE; leftFragIndex++) {
        for (rightFragIndex = 0; rightFragIndex < numRightFrags && closedGap == FALSE; rightFragIndex++) {

          lFragIid = leftExtFragsArray[leftFragIndex].fragIid;
          rFragIid = rightExtFragsArray[rightFragIndex].fragIid;

          closedGap = FALSE;

          //fprintf(stderr, "attempting to close gap %d, contigs %d and %d, fragIids %9d and %9d\n", gapNumber, lcontig->id, rcontig->id, lFragIid, rFragIid);

          if (debug.eCRmainLV > 0)
            fprintf(debug.eCRmainFP, "examining frags %d and %d\n",
                    leftExtFragsArray[leftFragIndex].fragIid, 
                    rightExtFragsArray[rightFragIndex].fragIid);

          if ((gapSize.mean - leftExtFragsArray[leftFragIndex].ctgMaxExt - rightExtFragsArray[rightFragIndex].ctgMaxExt) > (NUM_STDDEV_CUTOFF * sqrt(gapSize.variance)) &&
              (gapSize.mean > 100.0)) {

            if (debug.eCRmainLV > 0) {
              fprintf(debug.eCRmainFP, "leftExtFragsArray[%d].ctgMaxExt: %10d, rightExtFragsArray[%d].ctgMaxExt: %10d\n",
                      leftFragIndex,  leftExtFragsArray[leftFragIndex].ctgMaxExt, 
                      rightFragIndex, rightExtFragsArray[rightFragIndex].ctgMaxExt);
              fprintf(debug.eCRmainFP, "gap variance too large (gapSize - extensions = %.2f > %.1f * sqrt(gapSize.variance) = %.2f\n",
                      gapSize.mean - leftExtFragsArray[leftFragIndex].ctgMaxExt - rightExtFragsArray[rightFragIndex].ctgMaxExt, 
                      NUM_STDDEV_CUTOFF, NUM_STDDEV_CUTOFF * sqrt(gapSize.variance));
            }

            numGapsVarTooSmall++;
            continue;
          }

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
                         &rightFragFlapLength) == FALSE) {
            noOverlapFound++;
            continue;
          }

          //
          //  Abort the extension?
          //

          if (CONTIG_BASES < 2000) {
            int  ctglen = (MIN(CONTIG_BASES, (int) lcontig->bpLength.mean) +
                           MIN(CONTIG_BASES, (int) rcontig->bpLength.mean));

            if ((ahang + currLength + bhang - 1000) > ctglen)
              continue;
            if ((ahang + currLength + bhang + 700) < ctglen)
              continue;
          }


          //  our alignment didn't include any of the added bases,
          //  don't proceed.
          //
          if ((lFragIid != -1) &&
              (leftExtFragsArray[leftFragIndex].frgMaxExt < leftFragFlapLength))
            continue;

          if ((rFragIid != -1) &&
              (rightExtFragsArray[rightFragIndex].frgMaxExt < rightFragFlapLength))
            continue;

          //
          // extend the clear ranges of the frags
          //

          saveFragAndUnitigData(lFragIid, rFragIid);

          int lFragExt = 0;
          int rFragExt = 0;

          if (lFragIid != -1) {
            if (debug.eCRmainLV > 0)
              fprintf(debug.eCRmainFP,"adjusting left frg clear range by %d - %d = %d bases\n",
                      leftExtFragsArray[leftFragIndex].frgMaxExt,
                      leftFragFlapLength,
                      leftExtFragsArray[leftFragIndex].frgMaxExt - leftFragFlapLength);

            lFragExt = leftExtFragsArray[leftFragIndex].frgMaxExt - leftFragFlapLength;

            extendClearRange(lFragIid, lFragExt);
          }

          if (rFragIid != -1) {
            if (debug.eCRmainLV > 0)
              fprintf(debug.eCRmainFP,"adjusting right frg clear range by %d - %d = %d bases\n",
                      rightExtFragsArray[rightFragIndex].frgMaxExt,
                      rightFragFlapLength,
                      rightExtFragsArray[rightFragIndex].frgMaxExt - rightFragFlapLength);

            rFragExt = rightExtFragsArray[rightFragIndex].frgMaxExt - rightFragFlapLength;

            extendClearRange(rFragIid, rFragExt);
          }


          // the fragment extensions have succeeded

          int             gotNewLeftMA        = TRUE;
          int             gotNewRightMA       = TRUE;

          CIFragT        *frag                = NULL;
          NodeCGW_T      *unitig              = NULL;
				  
          if (debug.diagnosticLV > 0) {
            DumpContigMultiAlignInfo ("before anything (left)", NULL, lcontig->id);
            DumpContigUngappedOffsets("before anything (left)", lcontig->id);
            DumpContigMultiAlignInfo ("before anything (right)", NULL, rcontig->id);
            DumpContigUngappedOffsets("before anything (right)", rcontig->id);
          }

          // save the max offset of the right contig so we know
          // how to adjust the offsets of the contigs further
          // along the scaffold later

          double maxRContigOffset = MAX(rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean);
          int    nextContigIndex  = rcontig->BEndNext;
				  
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

            frag  = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
            unitig = GetGraphNode(ScaffoldGraph->CIGraph, frag->CIid);

            gotNewLeftMA = GetNewUnitigMultiAlign(unitig,
                                                  getAlteredFragPositions(unitig, lFragIid, leftExtFragsArray[leftFragIndex].ctgMaxExt - leftFragFlapLength),
                                                  lFragIid,
                                                  lFragExt);

            if (!gotNewLeftMA) {
              unitigMultiAlignFailures++;
            } else {
              SynchUnitigTWithMultiAlignT(unitig);

              if (debug.diagnosticLV > 0)
                DumpContigMultiAlignInfo ("before ReplaceEndUnitigInContig (left)", NULL, lcontig->id);

              MultiAlignT *newLeftMA = ReplaceEndUnitigInContig(ScaffoldGraph->sequenceDB,
                                                                ScaffoldGraph->gkpStore,
                                                                lcontig->id, unitig->id,
                                                                lcontig->offsetAEnd.mean >= lcontig->offsetBEnd.mean,
                                                                GlobalData->aligner,
                                                                NULL);
              if (newLeftMA) {
                if (debug.diagnosticLV > 0) {
                  DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (left)", NULL, lcontig->id);
                  DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (left)", newLeftMA, lcontig->id);
                }

                newLeftMA->maID = lcontig->id;
                updateMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, lcontig->id, FALSE, newLeftMA, TRUE);

                if (debug.eCRmainLV > 0)
                  fprintf(debug.eCRmainFP, "strlen(Getchar (newLeftMA->consensus)): " F_SIZE_T "\n",
                          strlen(Getchar (newLeftMA->consensus, 0)));

                if (debug.diagnosticLV > 0)
                  DumpContigMultiAlignInfo ("after updating store (left)", NULL, lcontig->id);

                if (lcontigOrientation == A_B)
                  extendContig(lcontig, FALSE);
                else
                  extendContig(lcontig, TRUE);
              } else {
                fprintf(stderr, "WARNING:  Failed to align the new and old (left) contigs.  Will not use this extension.\n");
                replaceEndUnitigFailures++;
                gotNewLeftMA = FALSE;
              }
            }
          }
				  
          // right unitig
          if (rFragIid != -1) {
            InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
            assert(info->set);

            frag  = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
            unitig = GetGraphNode(ScaffoldGraph->CIGraph, frag->CIid);

            gotNewRightMA = GetNewUnitigMultiAlign(unitig,
                                                   getAlteredFragPositions(unitig, rFragIid, rightExtFragsArray[rightFragIndex].ctgMaxExt - rightFragFlapLength),
                                                   rFragIid,
                                                   rFragExt);

            if (!gotNewRightMA) {
              unitigMultiAlignFailures++;
            } else {
              SynchUnitigTWithMultiAlignT(unitig);

              if (debug.diagnosticLV > 0)
                DumpContigMultiAlignInfo ("before ReplaceEndUnitigInContig (right)", NULL, rcontig->id);

              MultiAlignT *newRightMA = ReplaceEndUnitigInContig(ScaffoldGraph->sequenceDB,
                                                                 ScaffoldGraph->gkpStore,
                                                                 rcontig->id, unitig->id,
                                                                 rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean,
                                                                 GlobalData->aligner,
                                                                 NULL);
              if (newRightMA) {
                if (debug.diagnosticLV > 0) {
                  DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (right) (original contig)", NULL, rcontig->id);
                  DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (right) (new contig)", newRightMA, rcontig->id);
                }

                newRightMA->maID = rcontig->id;
                updateMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, rcontig->id, FALSE, newRightMA, TRUE);

                if (debug.eCRmainLV > 0)
                  fprintf(debug.eCRmainFP, "strlen(Getchar (newRightMA->consensus)): " F_SIZE_T "\n",
                          strlen(Getchar (newRightMA->consensus, 0)));

                if (debug.diagnosticLV > 0)
                  DumpContigMultiAlignInfo ("after updating store (right)", NULL, rcontig->id);

                if (rcontigOrientation == A_B)
                  extendContig(rcontig, TRUE);
                else
                  extendContig(rcontig, FALSE);
              } else {
                fprintf(stderr, "WARNING:  Failed to align the new and old (right) contigs.  Will not use this extension.\n");
                replaceEndUnitigFailures++;
                gotNewRightMA = FALSE;
              }
            }
          }  //  end of right unitig


          //  unwind changes if one or the other fails
          //
          if (gotNewLeftMA == FALSE || gotNewRightMA == FALSE)  
            goto dontkeepgap;


          double delta;

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

          LengthT    newOffsetAEnd = {0};
          LengthT    newOffsetBEnd = {0};

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

          if (FALSE == CreateAContigInScaffold(scaff, ContigPositions, newOffsetAEnd, newOffsetBEnd)) {
            fprintf(stderr, "overlap found in extendClearRanges but not by CNS!\n");
            createAContigFailures++;
            goto dontkeepgap;
          }

          //  We might have reallocated the contig graph; lookup
          //  the contig pointers again.
          //
          rcontig = GetGraphNode(ScaffoldGraph->ContigGraph, rcontigID);
          lcontig = GetGraphNode(ScaffoldGraph->ContigGraph, lcontigID);

          closedGap = TRUE;
					
          numGapsClosed++;

          totalBasesInClosedGaps += (int) gapSize.mean;

          totalOlapVariance += currLength * currLength;
          totalOlapLength += currLength;
          totalOlapDiffs += currDiffs;

          if (gapSize.mean < 100.0) {
            numSmallGapsClosed++;
            numSmallGapsClosedThisScaff++;
          } else {
            numLargeGapsClosed++;
            numLargeGapsClosedThisScaff++;
          }

          if (gapSize.mean > maxGapSizeClosed) {
            maxGapSizeClosed = gapSize.mean;
            maxGapSizeClosedNumber = gapNumber;
          }

          //  Print before we reassign to get the new contig
          fprintf(stderr, "closed gap %d, contigs %d and %d, fragIids %9d and %9d\n", gapNumber, lcontig->id, rcontig->id, lFragIid, rFragIid);

          // The new contig takes place of what was the right contig.
          //
          rcontig   = GetGraphNode(ScaffoldGraph->ContigGraph, GetNumGraphNodes(ScaffoldGraph->ContigGraph) - 1);
          rcontigID = GetNumGraphNodes(ScaffoldGraph->ContigGraph) - 1;

          // now we need to adjust contigs past the new contig, if not on end
          //
          if (nextContigIndex != NULLINDEX) {
            LengthT scaffoldDelta;

            scaffoldDelta.mean = MAX(rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean) - maxRContigOffset;
            scaffoldDelta.variance = ComputeFudgeVariance(scaffoldDelta.mean);

            AddDeltaToScaffoldOffsets(ScaffoldGraph, scaff->id, nextContigIndex, TRUE, FALSE, scaffoldDelta);
          }

          continue;


        dontkeepgap:
          // if we didn't close gap for whatever reason undo all
          // the frag and unitig changes

          fprintf(stderr, "DONTKEEPGAP!  did not close gap %d, contigs %d and %d, fragIids %9d and %9d\n", gapNumber, lcontig->id, rcontig->id, lFragIid, rFragIid);

          closedGap = FALSE;

          restoreFragAndUnitigData(lFragIid, rFragIid);

          memcpy(lcontig, &lcontigBackup, sizeof(ContigT));
          memcpy(rcontig, &rcontigBackup, sizeof(ContigT));

        }  //  over all right frags
      }  //  over all left frags

      gapNumber++;

      lcontig   = rcontig;
      lcontigID = lcontig->id;
      rcontigID = lcontig->BEndNext;
    }  //  over all contigs in the scaffold


    fprintf(stderr, "scaffold stats, scaff %10d, smallGaps %8d closed %8d, largeGaps %8d closed %8d\n",
            scaff->id,
            numSmallGapsThisScaff, numSmallGapsClosedThisScaff,
            numLargeGapsThisScaff, numLargeGapsClosedThisScaff);


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

    //  Clear out any cached multialigns.  We're all done with them.
    //
    clearCacheSequenceDB(ScaffoldGraph->sequenceDB);
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

  if (t1->ctgMaxExt > t2->ctgMaxExt)
    return -1;
  if (t1->ctgMaxExt < t2->ctgMaxExt)
    return 1;
  return 0;
}



void
getExtendableClearRange(CDS_IID_t     iid,
                        unsigned int *clr_bgn,
                        unsigned int *clr_end,
                        unsigned int *len) {
  static fragRecord fsread;  //  static for performance only

  getFrag(ScaffoldGraph->gkpStore, iid, &fsread, FRAG_S_INF);

  *len = getFragRecordSequenceLength(&fsread);

  switch (iterNumber) {
    case 1:
      *clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR1 - 1);
      *clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR1 - 1);
      break;
    case 2:
      *clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR2 - 1);
      *clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR2 - 1);
      break;
    default:
      assert(0);
      break;
  }
}


// findFirstFrag looks for a 3p->5p frag at the low end of a contig
// basesToNextFrag has meaning only when the first frag is the end frag, since basesToNextFrag
// marks where we start to have 2x coverage and is then used to determine MaxBegGap or MaxEndGap
int
findFirstExtendableFrags(ContigT *contig, extendableFrag *extFragsArray) {

  MultiAlignT *ma                  = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  int          numFrags            = GetNumIntMultiPoss(ma->f_list);
  int          firstUnitigID       = GetIntUnitigPos(ma->u_list, 0)->ident;

  int          extendableFragCount = 0;
  int          secondFragStart     = contig->bpLength.mean;

  int          i;

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "in findFirstExtendableFrags, contig->bpLength.mean: %f, firstUnitigID: %d\n",
            contig->bpLength.mean, firstUnitigID);
  
  //                 <--------------------------------------------------------------------------- contig
  // 3p <---frgMaxExt------|---------------------------------------------|----- 5p frag
  //                    clr_end                                       clr_bgn
  //    |-----------|
  //      ctgMaxExt

  //             <--------------------------------------------------------------------------- contig
  //                     <-|---------------------------------------------|----- 5p frag
  //                    clr_end                                       clr_bgn
  //             |------|
  //               ctgMaxExt (negative)

  for (i=0; i<numFrags; i++) {
    IntMultiPos *mp   = GetIntMultiPos(ma->f_list, i);
    CIFragT     *frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) mp->sourceInt);

    if ((frag->contigOffset3p.mean < 100.0) &&                      // frag is within a cutoff of the low end of the contig
        (frag->contigOffset3p.mean < frag->contigOffset5p.mean) &&  // and points in the right direction
        (frag->cid == firstUnitigID)) {                             // and is in the first unitig
      unsigned int clr_bgn;
      unsigned int clr_end;
      unsigned int seq_len;

      getExtendableClearRange(frag->iid, &clr_bgn, &clr_end, &seq_len);

      int ext = seq_len - clr_end - frag->contigOffset3p.mean;

      if (ext > 30) {
        extFragsArray[extendableFragCount].fragIid           = frag->iid;
        extFragsArray[extendableFragCount].ctgMaxExt   = ext;
        extFragsArray[extendableFragCount].frgMaxExt  = seq_len - clr_end;
        extFragsArray[extendableFragCount].basesToNextFrag   = 0;
        extFragsArray[extendableFragCount].fragOnEnd         = FALSE;

        if ((int)frag->contigOffset3p.mean == 0)
          extFragsArray[extendableFragCount].fragOnEnd = TRUE;

        if (debug.eCRmainLV > 0)
          fprintf(debug.eCRmainFP, "frstExt: in contig %d, frag %d is at %f -> %f (5p->3p) -- ctgMaxExt %d, frgMaxExt %d\n", 
                  contig->id,
                  frag->iid,
                  frag->contigOffset5p.mean,
                  frag->contigOffset3p.mean,
                  extFragsArray[extendableFragCount].ctgMaxExt,
                  extFragsArray[extendableFragCount].frgMaxExt);

        extendableFragCount++;

        if (extendableFragCount >= MAX_EXTENDABLE_FRAGS)
          fprintf(stderr, "extendableFragCount (%d) is greater than MAX_EXTENDABLE_FRAGS, aborting...\n",
                  extendableFragCount);
        assert(extendableFragCount < MAX_EXTENDABLE_FRAGS);
      }
    }

    // secondFragStart is where the next to end frag starts, and thus
    // where we start 2x coverage

    if (((int)frag->contigOffset3p.mean > 0) && 
        ((int)frag->contigOffset3p.mean < secondFragStart))
      secondFragStart = (int)frag->contigOffset3p.mean;

    if ((int)frag->contigOffset5p.mean < secondFragStart)
      secondFragStart = (int)frag->contigOffset5p.mean;
  }


  // now sort the extendable frags by their extendability
  qsort(extFragsArray, extendableFragCount, sizeof(extendableFrag), compExtendableFrags);


  for (i=0; i<extendableFragCount; i++) {
    if (extFragsArray[i].fragOnEnd == TRUE)
      extFragsArray[i].basesToNextFrag = secondFragStart;

    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "contig %8d, frag %8d can extend %8d bases into the gap\n",
              contig->id, extFragsArray[i].fragIid, extFragsArray[i].ctgMaxExt);
  }


  return extendableFragCount;
}



// findLastFrag looks for a 5p->3p frag at the high end of a contig
// basesToNextFrag has meaning only when the first frag is the end frag, since basesToNextFrag
// marks where we start to have 2x coverage and is then used to determine MaxBegGap or MaxEndGap
int
findLastExtendableFrags(ContigT *contig, extendableFrag *extFragsArray) {

  MultiAlignT  *ma           = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  int           numFrags     = GetNumIntMultiPoss(ma->f_list);
  int           lastUnitigID = GetIntUnitigPos(ma->u_list, GetNumIntUnitigPoss(ma->u_list) - 1)->ident;  
  double        maxContigPos = contig->bpLength.mean - 1.0;

  int           extendableFragCount = 0;
  int           secondFragEnd       = 0;

  int           i;

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "in findLastExtendableFrags, contig->bpLength.mean: %f, firstUnitigID: %d\n",
            contig->bpLength.mean, lastUnitigID);
  
  //    contig ----------------------------------------------------------------------------------->
  //                                  5p -------|---------------------------------------------|------------> 3p 
  //                                         clr_bgn                                       clr_end
  //                                                                                               |-------|
  //                                                                                                ctgMaxExt

  for (i = 0; i < numFrags; i++) {
    IntMultiPos *mp   = GetIntMultiPos(ma->f_list, i);
    CIFragT     *frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) mp->sourceInt);

    if ((frag->contigOffset3p.mean > maxContigPos - 100.0) &&      // frag is within a cutoff of the high end of the contig
        (frag->contigOffset5p.mean < frag->contigOffset3p.mean) && // and points in the right direction
        (frag->cid == lastUnitigID)) {                             // and is in the last unitig

      unsigned int clr_bgn;
      unsigned int clr_end;
      unsigned int seq_len;

      getExtendableClearRange(frag->iid, &clr_bgn, &clr_end, &seq_len);

      int ext = seq_len - clr_end - (int) (contig->bpLength.mean - frag->contigOffset3p.mean);

      if (ext > 30) {
        extFragsArray[extendableFragCount].fragIid           = frag->iid;
        extFragsArray[extendableFragCount].ctgMaxExt   = ext;
        extFragsArray[extendableFragCount].frgMaxExt  = seq_len - clr_end;
        extFragsArray[extendableFragCount].basesToNextFrag   = 0;
        extFragsArray[extendableFragCount].fragOnEnd         = FALSE;

        if ((int)frag->contigOffset3p.mean == (int)contig->bpLength.mean)
          extFragsArray[extendableFragCount].fragOnEnd = TRUE;

        if (debug.eCRmainLV > 0)
          fprintf(debug.eCRmainFP, "lastExt: in contig %d, frag %d is at %f -> %f (5p->3p) -- ctgExt %d, frgExt %d\n", 
                  contig->id,
                  frag->iid,
                  frag->contigOffset5p.mean,
                  frag->contigOffset3p.mean,
                  extFragsArray[extendableFragCount].ctgMaxExt,
                  extFragsArray[extendableFragCount].frgMaxExt);

        extendableFragCount++;

        if (extendableFragCount >= MAX_EXTENDABLE_FRAGS)
          fprintf(stderr, "extendableFragCount (%d) is greater than MAX_EXTENDABLE_FRAGS, aborting...\n",
                  extendableFragCount);
        assert(extendableFragCount < MAX_EXTENDABLE_FRAGS);
      }
    }

    // secondFragEnd is where the next to end frag ends, and thus where we have 2x coverage
    // we don't care if it's the 3p or 5p end

    if (((int)frag->contigOffset3p.mean < (int)contig->bpLength.mean) &&
        ((int)frag->contigOffset3p.mean > secondFragEnd))
      secondFragEnd = (int) frag->contigOffset3p.mean;

    if ((int)frag->contigOffset5p.mean > secondFragEnd)
      secondFragEnd = frag->contigOffset5p.mean;
  }


  // now sort the extendable frags by their extendability
  qsort(extFragsArray, extendableFragCount, sizeof(extendableFrag), compExtendableFrags);


  for (i=0; i<extendableFragCount; i++) {
    if (extFragsArray[i].fragOnEnd == TRUE)
      extFragsArray[i].basesToNextFrag = maxContigPos - secondFragEnd;
	
    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "contig %8d, frag %8d can extend %8d bases into the gap\n",
              contig->id, extFragsArray[i].fragIid, extFragsArray[i].ctgMaxExt);
  }

  return extendableFragCount;
}



// findLastUnitig looks for a surrogate at the high end of a contig
int
findLastUnitig(ContigT *contig, int *unitigID) {
  int    i;
  int    isSurrogate = FALSE;
  double maxContigPos = 0.0;
  
  MultiAlignT *ma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  
  // can't just jump to last unitig since unitigs are arranged by starting position, not ending

  for (i=0; i<GetNumIntUnitigPoss(ma->u_list); i++) {
    IntUnitigPos *upos   = GetIntUnitigPos(ma->u_list, i);
    NodeCGW_T    *unitig = GetGraphNode(ScaffoldGraph->CIGraph, upos->ident);
	
    if ((unitig->offsetAEnd.mean >= maxContigPos) ||
        (unitig->offsetBEnd.mean >= maxContigPos)) {
      maxContigPos = MAX(unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);
      *unitigID    = unitig->id;
      isSurrogate  = unitig->flags.bits.isSurrogate;

      if (debug.eCRmainLV > 0)
        fprintf(debug.eCRmainFP, "findLastUnitig()-- in contig %d, unitig %d is at %f -> %f maxContigPos: %f, isSurrogate: %d, baseID: %d\n", 
                contig->id,
                unitig->id,
                unitig->offsetAEnd.mean,
                unitig->offsetBEnd.mean,
                maxContigPos,
                isSurrogate,
                unitig->info.CI.baseID);
    }
  }

  return(isSurrogate);
}



// findFirstUnitig looks for a surrogate at the low end of a contig
int
findFirstUnitig(ContigT *contig, int *unitigID) {

  MultiAlignT  *ma          = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  IntUnitigPos *upos        = GetIntUnitigPos(ma->u_list, 0);
  NodeCGW_T    *unitig      = GetGraphNode(ScaffoldGraph->CIGraph, upos->ident);
  int           isSurrogate = unitig->flags.bits.isSurrogate;

  *unitigID = unitig->id;

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "findFirstUnitig()-- in contig %d, unitig %d is at %f -> %f, isSurrogate: %d\n", 
            contig->id,
            unitig->id,
            unitig->offsetAEnd.mean,
            unitig->offsetBEnd.mean,
            isSurrogate);

  return(isSurrogate);
}






void
extendContig(ContigT *contig, int extendAEnd) {
  int           contigOrientation = (contig->offsetAEnd.mean < contig->offsetBEnd.mean) ? A_B : B_A;
  MultiAlignT  *cma               = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
  double        lengthDelta       = strlen(Getchar(cma->consensus, 0)) - contig->bpLength.mean;

  // have to alter the following fields in a NodeCGW_T: bpLength, offsetAEnd, offsetBEnd

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "extendContig()-- contig %8d original pos: %.0f, %.0f %s\n",
            contig->id,
            (contigOrientation == A_B) ? contig->offsetAEnd.mean : contig->offsetBEnd.mean,
            (contigOrientation == A_B) ? contig->offsetBEnd.mean : contig->offsetAEnd.mean,
            (contigOrientation == A_B) ? "A_B" : "B_A");
  
  contig->bpLength.mean += lengthDelta;
  
  if (contigOrientation == A_B) {
    if (extendAEnd == TRUE)
      contig->offsetAEnd.mean -= lengthDelta;
    else
      contig->offsetBEnd.mean += lengthDelta;
  } else {
    if (extendAEnd == TRUE)
      contig->offsetAEnd.mean += lengthDelta;
    else
      contig->offsetBEnd.mean -= lengthDelta;
  }

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "extendContig()-- contig %8d original pos: %.0f, %.0f %s\n",
            contig->id,
            (contigOrientation == A_B) ? contig->offsetAEnd.mean : contig->offsetBEnd.mean,
            (contigOrientation == A_B) ? contig->offsetBEnd.mean : contig->offsetAEnd.mean,
            (contigOrientation == A_B) ? "A_B" : "B_A");
}






fragPositions *
getAlteredFragPositions(NodeCGW_T *unitig, int alteredFragIid, int ctgExt) {
  fragPositions *fp    = NULL;
  int            i     = 0;
  int            afidx = NULLINDEX;

  if (ctgExt < 0)
    fprintf(stderr, "getAlteredFragPositions()-- negative contig extension (%d); hope you're in the first frag!\n", ctgExt);

  MultiAlignT *uma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);

  fp = (fragPositions *)safe_malloc(GetNumIntMultiPoss(uma->f_list) * sizeof(fragPositions));

  // get the current positions
  for (i=0; i<GetNumIntMultiPoss(uma->f_list); i++) {
    IntMultiPos *pos = GetIntMultiPos(uma->f_list, i);
  
    fp[i].bgn = pos->position.bgn;
    fp[i].end = pos->position.end;

    if (pos->ident == alteredFragIid)
      afidx = i;

    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "getAlteredFragPositions()-- %2d] (%d,%d)\n", i, fp[i].bgn, fp[i].end);
  }

  assert(afidx != NULLINDEX);

  if (fp[afidx].bgn < fp[afidx].end) {

    //  Nothing much to do, just extend the altered fragment -- it's
    //  on the end.  First find out how much must be added to the frag
    //  to get to the end of the untig, then add in how far it extends
    //  into the gap, then add that to the end.
    //
    //  -------------
    //     <---------+++++
    //
    fp[afidx].end = strlen(Getchar(uma->consensus, 0)) + ctgExt;

    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "getAlteredFragPositions()-- adjusted extended fragment %d to (%d,%d)\n",
              afidx, fp[afidx].bgn, fp[afidx].end);
  } else {
    //  All frag positions get bumped by 'ctgExt', that's how far
    //  the altered frag extends into the gap.  This is the new
    //  minimum position in the unitig.
    //
    //        -------------
    //   +++++------->
    //         <--------
    //

    for (i = 0; i < GetNumIntMultiPoss(uma->f_list); i++) {
      fp[i].bgn += ctgExt;
      fp[i].end += ctgExt;
      if (debug.eCRmainLV > 0)
        fprintf(debug.eCRmainFP, "getAlteredFragPositions()-- %2d] (%d,%d) (adjusted by %d)\n", i, fp[i].bgn, fp[i].end, ctgExt);
    }

    //  We're now the start of the unitig.
    fp[afidx].end = 0;

    if (debug.eCRmainLV > 0)
      fprintf(debug.eCRmainFP, "getAlteredFragPositions()-- adjusted extended fragment %d to (%d,%d)\n",
              afidx, fp[afidx].bgn, fp[afidx].end);
  }

  return(fp);
}




int
GetNewUnitigMultiAlign(NodeCGW_T *unitig,
                       fragPositions *fragPoss,
                       int extFragIID,
                       int extLength) {
  GenericMesg	    pmesg;
  IntUnitigMesg	    ium_mesg = {0};
  int               i;
  MultiAlignT      *macopy = NULL;

  //  If fragPoss is empty, then getAlteredFragPositions() failed and
  //  we should abort.  It currently never does, but we should check.
  //
  if (fragPoss == NULL)
    return(FALSE);

  static VA_TYPE(char) *ungappedSequence = NULL;
  static VA_TYPE(char) *ungappedQuality  = NULL;

  if (ungappedSequence== NULL) {
    ungappedSequence = CreateVA_char(0);
    ungappedQuality = CreateVA_char(0);
  } else {
    ResetVA_char(ungappedSequence);
    ResetVA_char(ungappedQuality);
  }

  pmesg.m = &ium_mesg;
  pmesg.t = MESG_IUM;

  assert((unitig->id >= 0) && (unitig->id < GetNumGraphNodes(ScaffoldGraph->CIGraph)));

  //  Make a copy of the multialign in the store.  This is needed
  //  because at the end of this function, we delete the multialign in
  //  the store, and insert a new one.  The new one is based on the
  //  copy....the construction of the new one uses an IUM message, and
  //  the IUM message uses pieces of macopy.
  //  
  macopy = CreateEmptyMultiAlignT();
  copyMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, macopy, unitig->id, TRUE);

  if (debug.eCRmainLV > 0) {
    //PrintMultiAlignT(stderr, loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE), ScaffoldGraph->gkpStore, 0, 0, AS_READ_CLEAR_ECR1);
    fprintf(debug.eCRmainFP, "for unitig %d, before reforming, strlen(macopy->consensus) = " F_SIZE_T "\n", unitig->id, strlen(Getchar(macopy->consensus, 0)));
    fprintf(debug.eCRmainFP, "for unitig %d, before reforming, consensus:\n%s\n", unitig->id, Getchar(macopy->consensus, 0));
  }

  assert (unitig->type != CONTIG_CGW);

  ium_mesg.iaccession     = unitig->id;
#ifdef AS_ENABLE_SOURCE
  ium_mesg.source         = NULL;
#endif
  ium_mesg.coverage_stat  = unitig->info.CI.coverageStat;
  ium_mesg.unique_rept    = unitig->unique_rept;

  switch(unitig->type) {
    case DISCRIMINATORUNIQUECHUNK_CGW:
      ium_mesg.status = AS_UNIQUE;
      break;
    case UNIQUECHUNK_CGW:
      ium_mesg.status = AS_UNIQUE;
      break;
    case UNRESOLVEDCHUNK_CGW:
      if (unitig->info.CI.numInstances > 0) {
        ium_mesg.status = AS_SEP;
        assert(!unitig->flags.bits.isUnique);
      } else if (unitig->scaffoldID != NULLINDEX) {
        ium_mesg.status = AS_UNIQUE;
      } else {
        ium_mesg.status = AS_NOTREZ;
      }
      break;
    default:
      assert(0);
  }

  GetMultiAlignUngappedConsensus(macopy, ungappedSequence, ungappedQuality);
  ium_mesg.consensus  = Getchar(ungappedSequence,0);
  ium_mesg.quality    = Getchar(ungappedQuality,0);
  ium_mesg.length     = GetMultiAlignUngappedLength(macopy);
  ium_mesg.forced     = 0;
  ium_mesg.num_frags  = GetNumIntMultiPoss(macopy->f_list);

  // replace the positions in the f_list with the adjusted positions
  int extendedFragLeftward = FALSE;

  for (i = 0; i < GetNumIntMultiPoss(macopy->f_list); i++) {
    IntMultiPos *tempPos = GetIntMultiPos(macopy->f_list, i);

    if (tempPos->ident == extFragIID) {
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
    leftShiftIUM(ium_mesg.f_list, GetNumIntMultiPoss(macopy->f_list), extFragIID);
  else
    // might need to reorder frags in f_list if we extended a frag rightward
    rightShiftIUM(ium_mesg.f_list, GetNumIntMultiPoss(macopy->f_list), extFragIID);

  //  The Local_Overlap_AS_forCNS below also had a commented out
  //  DP_Compare.
  //
  CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                          CNS_OPTIONS_SMOOTH_WIN_DEFAULT,
                          CNS_OPTIONS_MAX_NUM_ALLELES };
  ALIGNMENT_CONTEXT=AS_CONSENSUS;

  if (EXIT_FAILURE == MultiAlignUnitig(&ium_mesg, 
                                       ScaffoldGraph->gkpStore,
                                       reformed_consensus,
                                       reformed_quality,
                                       reformed_deltas,
                                       CNS_QUIET,  //  CNS_VERBOSE, CNS_STATS_ONLY
                                       1,
                                       Local_Overlap_AS_forCNS,
                                       &options)) {
    fprintf(stderr, "GetNewUnitigMultiAlign()-- MultiAlignUnitig failure on unitig %d\n", unitig->id);
    DeleteMultiAlignT(macopy);  //  We own this.
    return FALSE;
  }

  {
    MultiAlignT  *nma = CreateMultiAlignTFromIUM(&ium_mesg, -2, FALSE);

    //  If our new multialign differs wildly from the expected length,
    //  we assume the fragment extension was bogus --
    //  MultiAlignUnitig, instead of aligning the extension to the
    //  unitig, just inserted gaps (where x's are gaps).
    //
    //  Before:    ----------------   After:  -------xxxx-------
    //             <------                       <---xxxx---
    //           +++++++<-------              +++xxxx++++<-------
    //
    if (GetMultiAlignLength(nma) > GetMultiAlignLength(macopy) + 1.1 * extLength) {
      DeleteMultiAlignT(macopy);  //  We own this.
      DeleteMultiAlignT(nma);     //  We own this.
      fprintf(stderr, "WARNING:  new multialign is too long -- suspect alignment, don't use it.\n");
      return(FALSE);
    }

    fprintf(stderr, "WARNING:  GetNewUnitigMultiAlign()-- new length %d old length %d + %f\n",
            GetMultiAlignLength(nma),
            GetMultiAlignLength(macopy),
            1.1 * extLength);
        
    //  Set the keepInCache flag, since CreateMultiAlignTFromIUM() is
    //  returning an allocated multialign, and if we don't cache it,
    //  we have a giant memory leak.
    //
    updateMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg.iaccession, TRUE, nma, TRUE);
  }

  DeleteMultiAlignT(macopy);  //  We own this.

  if (debug.eCRmainLV > 0) {
    //PrintMultiAlignT(stderr, loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE), ScaffoldGraph->gkpStore, 0, 0, AS_READ_CLEAR_ECR1);
    fprintf(debug.eCRmainFP, "for unitig %d, after reforming, strlen(reformed_consensus) = " F_SIZE_T "\n", unitig->id, strlen(Getchar(reformed_consensus, 0)));
    fprintf(debug.eCRmainFP, "for unitig %d, after reforming, consensus:\n%s\n", unitig->id, Getchar(reformed_consensus, 0));
  }

  return(TRUE);
}




void
extendClearRange(int fragIid, int frag3pDelta) {
  unsigned int clr_bgn, clr_end;
  fragRecord fsread;

  if (fragIid == -1)
    return;

  if (frag3pDelta < 0)
    fprintf(stderr, "extendClearRange()-- WARNING: frag3pDelta less than zero: %d\n", frag3pDelta);
  assert(frag3pDelta >= 0);

  getFrag(ScaffoldGraph->gkpStore, fragIid, &fsread, FRAG_S_INF);

  switch (iterNumber) {
    case 1:
      clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR1 - 1);
      clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR1 - 1) + frag3pDelta;
      setFragRecordClearRegion(&fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR1);
      break;
    case 2:
      clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR2 - 1);
      clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR2 - 1) + frag3pDelta;
      setFragRecordClearRegion(&fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR2);
      break;
    default:
      assert(0);
      break;
  }

  setFrag(ScaffoldGraph->gkpStore, fragIid, &fsread);

  fprintf(stderr, "WARNING:  Set clear for frag "F_IID","F_UID" from (%d,%d) to (%d,%d)\n",
          getFragRecordIID(&fsread),
          getFragRecordUID(&fsread),
          clr_bgn,
          clr_end - frag3pDelta,
          clr_bgn,
          clr_end);
}


void
revertClearRange(int fragIid) {
  unsigned int clr_bgn, clr_end;
  fragRecord fsread;

  if (fragIid == -1)
    return;

  getFrag(ScaffoldGraph->gkpStore, fragIid, &fsread, FRAG_S_INF);

  switch (iterNumber) {
    case 1:
      clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR1 - 1);
      clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR1 - 1);
      setFragRecordClearRegion(&fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR1);
      break;
    case 2:
      clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_ECR2 - 1);
      clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_ECR2 - 1);
      setFragRecordClearRegion(&fsread, clr_bgn, clr_end, AS_READ_CLEAR_ECR2);
      break;
    default:
      assert(0);
      break;
  }

  setFrag(ScaffoldGraph->gkpStore, fragIid, &fsread);

  fprintf(stderr, "WARNING:  Revert clear for frag "F_IID","F_UID" to (%d,%d)\n",
          getFragRecordIID(&fsread),
          getFragRecordUID(&fsread),
          clr_bgn,
          clr_end);
}









// this routine shifts the frag with iid extFragIID to the front of the f_list
void
leftShiftIUM(IntMultiPos *f_list, int numFrags, int extFragIID) {
  int i, currPos = 0, numShiftedInUnitig = 0;
  IntMultiPos tempIMP;
  
  // find out the extended frag's current position
  for (i = 0; i < numFrags; i++)
    if (f_list[i].ident == extFragIID)
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
}


// this routine shifts the frag with iid extFragIID to the appropriate place in the f_list
void
rightShiftIUM(IntMultiPos *f_list, int numFrags, int extFragIID) {
  int i, currPos = 0;
  IntMultiPos tempIMP;
  int numShifted, shiftedFrag, fragToMovePos = NULLINDEX, numToShift;

  // find out the extended frag's current position
  for (i = 0; i < numFrags; i++)
    if (f_list[i].ident == extFragIID)
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
  
  if (savedLeftUnitigMA == NULL) {
    savedLeftUnitigMA  = CreateEmptyMultiAlignT();
    savedRightUnitigMA = CreateEmptyMultiAlignT();
    savedLeftContigMA  = CreateEmptyMultiAlignT();
    savedRightContigMA = CreateEmptyMultiAlignT();
  }

  if (lFragIid != -1) {
    InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);
    assert(info->set);

    CIFragT *leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

    copyMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, savedLeftUnitigMA, leftFrag->cid, TRUE);
    copyMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, savedLeftContigMA, leftFrag->contigID, FALSE);
  }
  
  if (rFragIid != -1) {
    InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
    assert(info->set);

    CIFragT *rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

    copyMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, savedRightUnitigMA, rightFrag->cid, TRUE);
    copyMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, savedRightContigMA, rightFrag->contigID, FALSE);
  }
}

void
restoreFragAndUnitigData(int lFragIid, int rFragIid) {

  revertClearRange(lFragIid);
  revertClearRange(rFragIid);				  

  if (lFragIid != -1) {
    InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);
    assert(info->set);

    CIFragT    *leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
    NodeCGW_T  *unitig   = GetGraphNode(ScaffoldGraph->CIGraph, leftFrag->cid);

    updateMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE, savedLeftUnitigMA, FALSE);
    updateMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, leftFrag->contigID, FALSE, savedLeftContigMA, FALSE);

    SynchUnitigTWithMultiAlignT(unitig);
  }
  
  if (rFragIid != -1) {
    InfoByIID *info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
    assert(info->set);

    CIFragT    *rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
    NodeCGW_T  *unitig    = GetGraphNode(ScaffoldGraph->CIGraph, rightFrag->cid);

    updateMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE, savedRightUnitigMA, FALSE);
    updateMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, rightFrag->contigID, FALSE, savedRightContigMA, FALSE);

    SynchUnitigTWithMultiAlignT(unitig);
  }
}
