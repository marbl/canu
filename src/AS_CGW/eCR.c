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

const char *mainid = "$Id: eCR.c,v 1.59 2011-11-29 11:50:00 brianwalenz Exp $";

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

static MultiAlignT      *savedLeftUnitigMA = NULL;
static MultiAlignT      *savedRightUnitigMA = NULL;
static MultiAlignT      *savedLeftContigMA = NULL;
static MultiAlignT      *savedRightContigMA = NULL;

int                      iterNumber = 1;

debugflags_t             debug = {0, 0L, 0, 0L, 0, 0L};




void
SynchUnitigTWithMultiAlignT(NodeCGW_T *unitig) {
  unitig->bpLength.mean = GetMultiAlignUngappedLength(ScaffoldGraph->tigStore->loadMultiAlign(unitig->id,
                                                                                              TRUE));
}



int
main(int argc, char **argv) {

  //  Command line args and processing
  int   scaffoldBegin    = -1;
  int   startingGap      = -1;
  int   scaffoldEnd      = -1;
  int   ckptNum          = -1;
  int   gkpPart          = 0;
  int   arg              = 1;
  int   err              = 0;

  //  Scaffold lists
  int  *scfSkip = NULL, scfSkipLen = 0;
  int  *scfOnly = NULL, scfOnlyLen = 0;
  int  *gapSkip = NULL, gapSkipLen = 0;
  int  *gapOnly = NULL, gapOnlyLen = 0;

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

  GlobalData           = new Globals_CGW();

  //  Could be cleaner (allocate only if options are present) but not
  //  simpler.
  scfSkip = (int *)safe_malloc(sizeof(int) * argc);
  scfOnly = (int *)safe_malloc(sizeof(int) * argc);
  gapSkip = (int *)safe_malloc(sizeof(int) * argc);
  gapOnly = (int *)safe_malloc(sizeof(int) * argc);

  while (arg < argc) {
    if        (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->outputPrefix, argv[++arg]);

    } else if (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->gkpStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      strcpy(GlobalData->tigStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-C") == 0) {
      startingGap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-n") == 0) {
      ckptNum = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-p") == 0) {
      gkpPart = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {
      scaffoldBegin = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      scaffoldEnd   = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-i") == 0) {
      iterNumber = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-v") == 0) {
      debug.eCRmainLV    = 1;
      debug.examineGapLV = 1;
      debug.diagnosticLV = 1;

    } else if (strcmp(argv[arg], "-o") == 0) {
      scfOnly[scfOnlyLen++] = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      scfSkip[scfSkipLen++] = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-O") == 0) {
      gapOnly[gapOnlyLen++] = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-S") == 0) {
      gapSkip[gapSkipLen++] = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      arg++;
      if (strcmp(argv[arg], "verbosemultialign") == 0) {
        VERBOSE_MULTIALIGN_OUTPUT = 1;
      }

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
  if ((GlobalData->outputPrefix[0] == 0) ||
      (GlobalData->gkpStoreName[0] == 0) ||
      (err)) {
    fprintf(stderr, "usage: %s [opts] -c ckpName -n ckpNumber -g gkpStore\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c ckpName     Use ckpName as the checkpoint name\n");
    fprintf(stderr, "  -n ckpNumber   The checkpoint to use\n");
    fprintf(stderr, "  -g gkpStore    The gatekeeper store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -C gap#        Start at a specific gap number\n");
    fprintf(stderr, "  -b scafBeg     Begin at a specific scaffold\n");
    fprintf(stderr, "  -e scafEnd     End after a specific scaffold (INCLUSIVE)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o scafIID     Process only this scaffold\n");
    fprintf(stderr, "  -s scafIID     Skip this scaffold\n");
    fprintf(stderr, "  -O gap#        Process only this gap\n");
    fprintf(stderr, "  -S gap#        Skip this gap\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -i iterNum     The iteration of ECR; either 1 or 2\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -p partition   Load a gkpStore partition into memory\n");
    exit(1);
  }

  LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix, ckptNum, TRUE);

  //  Create the starting clear range backup, if we're the first
  //  iteration.  We copy this from the LATEST clear, which will
  //  either by CLR or OBTCHIMERA.  This range does not exist before
  //  the first run of ECR, and we must create it.
  //
  if (iterNumber == 1)
    ScaffoldGraph->gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_ECR_0);

  //  Now, enable the range that we're going to be creating with this
  //  run.  Similar to ECR_0, this range doesn't exist before the
  //  first run of ECR.
  //
  ScaffoldGraph->gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_ECR_0 + iterNumber);

  //  If we're partitioned, load the partition.
  //
  //  THIS IS BROKEN.  GKPSTORE DOES NOT ALLOW UPDATES TO PARTITIONS.
  //
  //if (gkpPart)
  //  ScaffoldGraph->gkpStore->gkStore_loadPartition(gkpPart);

  //  Update the begin/end scaffold ids.
  //
  if (scaffoldBegin == -1)
    scaffoldBegin = 0;

  if ((scaffoldEnd == -1) ||
      (scaffoldEnd >= GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph)))
    scaffoldEnd = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) - 1;

  if (scaffoldBegin > scaffoldEnd) {
    fprintf(stderr, "%s: Nothing to do, no scaffolds selected: -b = %d > -e = %d\n", argv[0], scaffoldBegin, scaffoldEnd);
    exit(0);
  }

  //
  //  Scan all the scaffolds, closing gaps.  Go!
  //

  for (sid = scaffoldBegin; sid <= scaffoldEnd; sid++) {
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

    //  Are we supposed to do this scaffold?
    //
    {
      int  skipScaffold = 0;
      int  s;

      if (scfOnlyLen > 0) {
        int  doScaffold = 0;

        for (s=0; s<scfOnlyLen; s++)
          if (scfOnly[s] == sid)
            doScaffold = 1;

        if (doScaffold == 0)
          skipScaffold = 1;
      }

      if (scfSkipLen > 0) {
        for (s=0; s<scfSkipLen; s++)
          if (scfSkip[s] == sid)
            skipScaffold = 1;
      }

      if (skipScaffold) {
        fprintf(stderr,"\n=====================================================================\n");
        fprintf(stderr,"skipping scaffold %d, size %f (command line told me to)\n", sid, scaff->bpLength.mean);
        continue;
      }
    }


    fprintf(stderr,"\n=====================================================================\n");
    fprintf(stderr,"examining scaffold %d, size %f\n", sid, scaff->bpLength.mean);

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
      gkFragment fr;

      while (contigID != -1) {
        ContigT     *contig = GetGraphNode(ScaffoldGraph->ContigGraph, contigID);
        MultiAlignT *ma     = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);

        for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
          AS_IID          iid = GetIntMultiPos(ma->f_list, i)->ident;
          uint32  bgnOld, bgnCur;
          uint32  endOld, endCur;

          ScaffoldGraph->gkpStore->gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);

          fr.gkFragment_getClearRegion(bgnOld, endOld, AS_READ_CLEAR_ECR_0 + iterNumber - 1);
          fr.gkFragment_getClearRegion(bgnCur, endCur, AS_READ_CLEAR_ECR_0 + iterNumber);

          //fprintf(stderr, "CLR %d %d,%d %d,%d\n", iid, bgnOld, endOld, bgnCur, endCur);

          //  These are violated if the clear range is not configured.
          assert(bgnOld < endOld);
          assert(bgnCur < endCur);

          if ((bgnOld != bgnCur) || (endOld != endCur)) {
            fprintf(stderr, "WARNING:  Reset clear for frag "F_IID",%s (ctg %d scf %d) from (%d,%d) to (%d,%d)\n",
                    fr.gkFragment_getReadIID(),
                    AS_UID_toString(fr.gkFragment_getReadUID()),
                    contigID, sid,
                    bgnCur, endCur,
                    bgnOld, endOld);
            fr.gkFragment_setClearRegion(bgnOld, endOld, AS_READ_CLEAR_ECR_0 + iterNumber);

            ScaffoldGraph->gkpStore->gkStore_setFragment(&fr);
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
      SequenceOrient lcontigOrientation;
      SequenceOrient rcontigOrientation;

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
        lcontigOrientation.setIsForward();
      else
        lcontigOrientation.setIsReverse();

      if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
        rcontigOrientation.setIsForward();
      else
        rcontigOrientation.setIsReverse();


      fprintf(stderr, "---------------------------------------------------------------\n");
      fprintf(stderr, "examining gap %d from lcontig %d (orient: %c, pos: %f, %f) to rcontig %d (orient: %c, pos: %f, %f), size: %lf \n",
              gapNumber,
              lcontig->id, lcontigOrientation.toLetter(), lcontig->offsetAEnd.mean, lcontig->offsetBEnd.mean,
              rcontig->id, rcontigOrientation.toLetter(), rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean,
              gapSize.mean);

      if (gapSize.mean < 100.0) {
        numSmallGaps++;
        numSmallGapsThisScaff++;
      } else {
        numLargeGaps++;
        numLargeGapsThisScaff++;
      }


      if (debug.eCRmainLV > 0)
        fprintf(debug.eCRmainFP, "\nexamining lcontig %d (orientation %c) \n", lcontig->id, lcontigOrientation.toLetter());

      // find the extreme read on the correct end of the lcontig
      if (lcontigOrientation.isForward()) {
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
        fprintf(debug.eCRmainFP, "finished examining lcontig %d (orientation %c) \n", lcontig->id, lcontigOrientation.toLetter());
        fprintf(debug.eCRmainFP, "\nexamining rcontig %d (orientation %c) \n", rcontig->id, rcontigOrientation.toLetter());
      }

      // find the extreme read on the correct end of the rchunk
      if (rcontigOrientation.isForward()) {
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
        fprintf(debug.eCRmainFP, "finished examining rcontig %d (orientation %c) \n", rcontig->id, rcontigOrientation.toLetter());


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


      //  Are we supposed to do this gap?
      //
      {
        int  skipGap = 0;
        int  s;

        if (gapOnlyLen > 0) {
          int  doGap = 0;

          for (s=0; s<gapOnlyLen; s++)
            if (gapOnly[s] == gapNumber)
              doGap = 1;

          if (doGap == 0)
            skipGap = 1;
        }

        if (gapSkipLen > 0) {
          for (s=0; s<gapSkipLen; s++)
            if (gapSkip[s] == gapNumber)
              skipGap = 1;
        }

        if (skipGap) {
          fprintf(stderr, "skipping gap %d (command line told me to)\n", gapNumber);
          numLeftFrags = numRightFrags = 0;
        }
      }


      //  Skip this gap if it is before the one we want to start at.
      //
      if (gapNumber < startingGap)
        numLeftFrags = numRightFrags = 0;


      //  This is a good place to check rcontig->id, lcontig->id,
      //  gapNumber, etc and abort if that gap is causing problems.
      //
      //if (rcontig->id == 405527)
      //  numLeftFrags = numRightFrags = 0;


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
            if (GetCIFragT(ScaffoldGraph->CIFrags, lFragIid)->cid != lunitigID)
              continue;
          }

          if (rFragIid != -1) {
            if (GetCIFragT(ScaffoldGraph->CIFrags, rFragIid)->cid != runitigID)
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
            frag  = GetCIFragT(ScaffoldGraph->CIFrags, lFragIid);
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

              MultiAlignT *newLeftMA = ReplaceEndUnitigInContig(lcontig->id, unitig->id,
                                                                lcontig->offsetAEnd.mean >= lcontig->offsetBEnd.mean,
                                                                NULL);
              if (newLeftMA) {
                if (debug.diagnosticLV > 0) {
                  DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (left)", NULL, lcontig->id);
                  DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (left)", newLeftMA, lcontig->id);
                }

                newLeftMA->maID = lcontig->id;
                ScaffoldGraph->tigStore->insertMultiAlign(newLeftMA, FALSE, TRUE);

                if (debug.eCRmainLV > 0)
                  fprintf(debug.eCRmainFP, "strlen(Getchar (newLeftMA->consensus)): " F_SIZE_T "\n",
                          strlen(Getchar (newLeftMA->consensus, 0)));

                if (debug.diagnosticLV > 0)
                  DumpContigMultiAlignInfo ("after updating store (left)", NULL, lcontig->id);

                if (lcontigOrientation.isForward())
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
            frag  = GetCIFragT(ScaffoldGraph->CIFrags, rFragIid);
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

              MultiAlignT *newRightMA = ReplaceEndUnitigInContig(rcontig->id, unitig->id,
                                                                 rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean,
                                                                 NULL);
              if (newRightMA) {
                if (debug.diagnosticLV > 0) {
                  DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (right) (original contig)", NULL, rcontig->id);
                  DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (right) (new contig)", newRightMA, rcontig->id);
                }

                newRightMA->maID = rcontig->id;
                ScaffoldGraph->tigStore->insertMultiAlign(newRightMA, FALSE, TRUE);

                if (debug.eCRmainLV > 0)
                  fprintf(debug.eCRmainFP, "strlen(Getchar (newRightMA->consensus)): " F_SIZE_T "\n",
                          strlen(Getchar (newRightMA->consensus, 0)));

                if (debug.diagnosticLV > 0)
                  DumpContigMultiAlignInfo ("after updating store (right)", NULL, rcontig->id);

                if (rcontigOrientation.isForward())
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
                    (lcontigOrientation.isForward()) ? lcontig->offsetAEnd.mean : lcontig->offsetBEnd.mean,
                    (lcontigOrientation.isForward()) ? lcontig->offsetBEnd.mean : lcontig->offsetAEnd.mean);
            fprintf(debug.eCRmainFP, "                rctg: %12.0f, %12.0f\n",
                    (rcontigOrientation.isForward()) ? rcontig->offsetAEnd.mean : rcontig->offsetBEnd.mean,
                    (rcontigOrientation.isForward()) ? rcontig->offsetBEnd.mean : rcontig->offsetAEnd.mean);
          }

          // setup for contig merge

          if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
            delta = lcontig->offsetAEnd.mean;
          else
            delta = lcontig->offsetBEnd.mean;

          {
            VA_TYPE(IntElementPos) *ContigPositions = CreateVA_IntElementPos(2);

            IntElementPos   contigPos;

            contigPos.ident = lcontig->id;
            contigPos.type = AS_CONTIG;
            contigPos.position.bgn = lcontig->offsetAEnd.mean - delta;
            contigPos.position.end = lcontig->offsetBEnd.mean - delta;
            AppendIntElementPos(ContigPositions, &contigPos);

            if (debug.eCRmainLV > 0)
              fprintf(debug.eCRmainFP, "lcontig %8d positioned at %8d, %8d\n",
                      lcontig->id,contigPos.position.bgn, contigPos.position.end);

            contigPos.ident = rcontig->id;
            contigPos.type = AS_CONTIG;
            contigPos.position.bgn = rcontig->offsetAEnd.mean - delta;
            contigPos.position.end = rcontig->offsetBEnd.mean - delta;
            AppendIntElementPos(ContigPositions, &contigPos);

            if (debug.eCRmainLV > 0)
              fprintf(debug.eCRmainFP, "rcontig %8d positioned at %8d, %8d\n",
                      rcontig->id,contigPos.position.bgn, contigPos.position.end);

            LengthT    newOffsetAEnd = {0, 0};
            LengthT    newOffsetBEnd = {0, 0};

            if (lcontigOrientation.isForward())
              newOffsetAEnd.mean = lcontig->offsetAEnd.mean;
            else
              newOffsetAEnd.mean = lcontig->offsetBEnd.mean;

            if (rcontigOrientation.isForward())
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

            closedGap = CreateAContigInScaffold(scaff, ContigPositions, newOffsetAEnd, newOffsetBEnd);

            Delete_VA(ContigPositions);
          }

          //  We might have reallocated the contig graph; lookup
          //  the contig pointers again (yes, even if we fail).
          //
          rcontig = GetGraphNode(ScaffoldGraph->ContigGraph, rcontigID);
          lcontig = GetGraphNode(ScaffoldGraph->ContigGraph, lcontigID);

          if (closedGap == FALSE) {
            fprintf(stderr, "overlap found in extendClearRanges but not by CNS!\n");
            createAContigFailures++;
            goto dontkeepgap;
          }

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
          // the frag and unitig changes.

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
                                            sid,
                                            TRUE, TRUE, FALSE);
      }
    }

    //  Clear out any cached multialigns.  We're all done with them.
    //
    if ((sid % 1000) == 0)
      ScaffoldGraph->tigStore->flushCache();
  }  //  over all scaffolds




  //  Loading a checkpoint implicitly calls these -- and the
  //  downstream consumer of our checkpoint shouldn't modify the
  //  checkpoint if it is querying it (dumpDistanceUpdates, for
  //  example) -- so we just call them before the checkpoint is
  //  written.
  //
  SetCIScaffoldTLengths(ScaffoldGraph, TRUE);
  CheckCIScaffoldTs(ScaffoldGraph);

  CheckpointScaffoldGraph("extendClearRanges", "after extendClearRanges");

  DestroyScaffoldGraph(ScaffoldGraph);
  delete GlobalData;



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
  const extendableFrag * t1 = (extendableFrag *)s1;
  const extendableFrag * t2 = (extendableFrag *)s2;

  if (t1->ctgMaxExt > t2->ctgMaxExt)
    return -1;
  if (t1->ctgMaxExt < t2->ctgMaxExt)
    return 1;
  return 0;
}



void
getExtendableClearRange(AS_IID  iid,
                        uint32& clr_bgn,
                        uint32& clr_end,
                        uint32& len) {
  gkFragment fr;

  ScaffoldGraph->gkpStore->gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);

  len = fr.gkFragment_getSequenceLength();

  fr.gkFragment_getClearRegion(clr_bgn, clr_end, AS_READ_CLEAR_ECR_0 + iterNumber - 1);
  assert(clr_bgn < clr_end);
}


// findFirstFrag looks for a 3p->5p frag at the low end of a contig
// basesToNextFrag has meaning only when the first frag is the end frag, since basesToNextFrag
// marks where we start to have 2x coverage and is then used to determine MaxBegGap or MaxEndGap
int
findFirstExtendableFrags(ContigT *contig, extendableFrag *extFragsArray) {

  MultiAlignT *ma                  = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);
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
    CIFragT     *frag = GetCIFragT(ScaffoldGraph->CIFrags, mp->ident);

    if ((frag->contigOffset3p.mean < 100.0) &&                      // frag is within a cutoff of the low end of the contig
        (frag->contigOffset3p.mean < frag->contigOffset5p.mean) &&  // and points in the right direction
        (frag->cid == firstUnitigID)) {                             // and is in the first unitig
      uint32 clr_bgn;
      uint32 clr_end;
      uint32 seq_len;

      getExtendableClearRange(frag->read_iid, clr_bgn, clr_end, seq_len);

      int ext = seq_len - clr_end - frag->contigOffset3p.mean;

      if (ext > 30) {
        extFragsArray[extendableFragCount].fragIid           = frag->read_iid;
        extFragsArray[extendableFragCount].ctgMaxExt   = ext;
        extFragsArray[extendableFragCount].frgMaxExt  = seq_len - clr_end;
        extFragsArray[extendableFragCount].basesToNextFrag   = 0;
        extFragsArray[extendableFragCount].fragOnEnd         = FALSE;

        if ((int)frag->contigOffset3p.mean == 0)
          extFragsArray[extendableFragCount].fragOnEnd = TRUE;

        if (debug.eCRmainLV > 0)
          fprintf(debug.eCRmainFP, "frstExt: in contig %d, frag %d is at %f -> %f (5p->3p) -- ctgMaxExt %d, frgMaxExt %d\n",
                  contig->id,
                  frag->read_iid,
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

  MultiAlignT  *ma           = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);
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
    CIFragT     *frag = GetCIFragT(ScaffoldGraph->CIFrags, mp->ident);

    if ((frag->contigOffset3p.mean > maxContigPos - 100.0) &&      // frag is within a cutoff of the high end of the contig
        (frag->contigOffset5p.mean < frag->contigOffset3p.mean) && // and points in the right direction
        (frag->cid == lastUnitigID)) {                             // and is in the last unitig

      uint32 clr_bgn;
      uint32 clr_end;
      uint32 seq_len;

      getExtendableClearRange(frag->read_iid, clr_bgn, clr_end, seq_len);

      int ext = seq_len - clr_end - (int) (contig->bpLength.mean - frag->contigOffset3p.mean);

      if (ext > 30) {
        extFragsArray[extendableFragCount].fragIid           = frag->read_iid;
        extFragsArray[extendableFragCount].ctgMaxExt   = ext;
        extFragsArray[extendableFragCount].frgMaxExt  = seq_len - clr_end;
        extFragsArray[extendableFragCount].basesToNextFrag   = 0;
        extFragsArray[extendableFragCount].fragOnEnd         = FALSE;

        if ((int)frag->contigOffset3p.mean == (int)contig->bpLength.mean)
          extFragsArray[extendableFragCount].fragOnEnd = TRUE;

        if (debug.eCRmainLV > 0)
          fprintf(debug.eCRmainFP, "lastExt: in contig %d, frag %d is at %f -> %f (5p->3p) -- ctgExt %d, frgExt %d\n",
                  contig->id,
                  frag->read_iid,
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

  MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);

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

  MultiAlignT  *ma          = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);
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
  int           contigIsForward   = (contig->offsetAEnd.mean < contig->offsetBEnd.mean) ? true : false;
  MultiAlignT  *cma               = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);
  double        lengthDelta       = strlen(Getchar(cma->consensus, 0)) - contig->bpLength.mean;

  // have to alter the following fields in a NodeCGW_T: bpLength, offsetAEnd, offsetBEnd

  if (debug.eCRmainLV > 0)
    fprintf(debug.eCRmainFP, "extendContig()-- contig %8d original pos: %.0f, %.0f %s\n",
            contig->id,
            (contigIsForward) ? contig->offsetAEnd.mean : contig->offsetBEnd.mean,
            (contigIsForward) ? contig->offsetBEnd.mean : contig->offsetAEnd.mean,
            (contigIsForward) ? "A_B" : "B_A");

  contig->bpLength.mean += lengthDelta;

  if (contigIsForward) {
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
            (contigIsForward) ? contig->offsetAEnd.mean : contig->offsetBEnd.mean,
            (contigIsForward) ? contig->offsetBEnd.mean : contig->offsetAEnd.mean,
            (contigIsForward) ? "A_B" : "B_A");
}






fragPositions *
getAlteredFragPositions(NodeCGW_T *unitig, int alteredFragIid, int ctgExt) {
  fragPositions *fp    = NULL;
  int            i     = 0;
  int            afidx = NULLINDEX;

  if (ctgExt < 0)
    fprintf(stderr, "getAlteredFragPositions()-- negative contig extension (%d); hope you're in the first frag!\n", ctgExt);

  MultiAlignT *uma = ScaffoldGraph->tigStore->loadMultiAlign(unitig->id, TRUE);

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



//  Update containment relationships.  Fragment 'oldparent' is now
//  contained in 'newparent'.  For any other fragments that have a
//  parent of 'oldparent', make their parent be 'newparent'.  If the
//  fragment is already contained, we can skip the change though.
//
void
FixContainmentRelationships(int          oldparent,
                            int          newparent,
                            IntMultiPos *flst,
                            int          flen) {

  IntMultiPos  *oldpfrg = flst + oldparent;
  IntMultiPos  *newpfrg = flst + newparent;
  int           t;

  if (oldpfrg->contained == 0)
    return;

  for (t=oldparent+1; t<flen; t++) {

    //  This fragments parent is the 'oldparent' (now contained).  If
    //  we are not already contained, switch the parent to
    //  'newparent'...and recurse.
    //
    if ((flst[t].parent == oldpfrg->ident) &&
        (flst[t].contained == 0)) {
      int oldp = flst[t].parent;
      int olda = flst[t].ahang;
      int oldb = flst[t].bhang;

      flst[t].parent    = newpfrg->ident;
      flst[t].ahang     = MIN(flst[t].position.bgn, flst[t].position.end) - MIN(flst[newparent].position.bgn, flst[newparent].position.end);
      flst[t].bhang     = MAX(flst[t].position.bgn, flst[t].position.end) - MAX(flst[newparent].position.bgn, flst[newparent].position.end);
      flst[t].contained = ((flst[t].ahang > 0) && (flst[t].bhang < 0)) ? flst[t].parent : 0;

#if 0
      fprintf(stderr, "RESET fix-contain for id %d from %d,%d,%d to %d,%d,%d container=%d\n",
              flst[t].ident,
              oldp, olda, oldb, 
              flst[t].parent, flst[t].ahang, flst[t].bhang,
              flst[t].contained);
#endif

      if (flst[t].contained)
        FixContainmentRelationships(t, newparent, flst, flen);
    }
  }
}



int
GetNewUnitigMultiAlign(NodeCGW_T *unitig,
                       fragPositions *fragPoss,
                       int extFragIID,
                       int extLength) {

  //  If fragPoss is empty, then getAlteredFragPositions() failed and we should abort.  It currently
  //  never does, but we should check.
  //
  if (fragPoss == NULL)
    return(FALSE);

  assert(unitig->type != CONTIG_CGW);
  assert(unitig->id >= 0);
  assert(unitig->id < GetNumGraphNodes(ScaffoldGraph->CIGraph));

  MultiAlignT  *maorig = ScaffoldGraph->tigStore->loadMultiAlign(unitig->id, TRUE);
  MultiAlignT  *macopy = CopyMultiAlignT(NULL, maorig);

  if (debug.eCRmainLV > 0) {
    fprintf(debug.eCRmainFP, "for unitig %d, before reforming, strlen(macopy->consensus) = " F_SIZE_T "\n", unitig->id, strlen(Getchar(macopy->consensus, 0)));
    fprintf(debug.eCRmainFP, "for unitig %d, before reforming, consensus:\n%s\n", unitig->id, Getchar(macopy->consensus, 0));
  }

  //  Update the f_list in the maccopy.  We do this in several passes, for clarity mostly.

  IntMultiPos *flst = GetIntMultiPos(macopy->f_list, 0);
  int          flen = GetNumIntMultiPoss(macopy->f_list);
  int          efrg = 0;

  //
  //  Update the position of all fragments.
  //
  for (int32 i=0; i<flen; i++) {
#if 0
    fprintf(stderr, "RESET position for frag i=%d ident=%d from %d,%d to %d,%d\n",
            i,
            flst[i].ident,
            flst[i].position.bgn, flst[i].position.end,
            fragPoss[i].bgn, fragPoss[i].end);
#endif
    flst[i].position.bgn = fragPoss[i].bgn;
    flst[i].position.end = fragPoss[i].end;

    if (flst[i].ident == extFragIID)
      efrg = i;
  }

  safe_free(fragPoss);
  fragPoss = NULL;

  //
  //  If the extended fragment now starts the unitig, make it be the first fragment in the f_list.
  //
  for (int32 i=0; i<flen; i++) {
    if (flst[i].ident == extFragIID) {
      if ((i != 0) &&
          ((flst[i].position.bgn == 0) ||
           (flst[i].position.end == 0))) {
#if 0
        fprintf(stderr, "RESET swap-first frag i=%d (ident=%d pos=%d,%d) is now first frag\n",
                i, flst[i].ident, flst[i].position.bgn, flst[i].position.end);
#endif
        IntMultiPos  newfirst = flst[i];
        memmove(flst + 1, flst, sizeof(IntMultiPos) * i);
        flst[0] = newfirst;
        efrg = 0;
      }
      break;
    }
  }


  //
  //  Reset the hangs of the new first fragment.  The second fragment MUST have a parent of the
  //  first fragment.  We blindly set it, leaving fixing the hangs for the final block.
  //
  if (efrg == 0) {
    flst[0].parent    = 0;
    flst[0].ahang     = 0;
    flst[0].bhang     = 0;
    flst[0].contained = 0;

    flst[1].parent    = extFragIID;
  }


  //
  //  Check all the fragments to see if anyone needs fixing.  The first two are already done, and we
  //  skip them.
  //
  //  If this fragment is contained, we need to propagate the fact that fragment i is now contained
  //  in fragment 0 to future fragments, so they can update any dovetail relationships.
  //
  for (int32 i=1; i<flen; i++) {
    if (flst[i].parent == extFragIID) {
      int oldp = flst[i].parent;
      int olda = flst[i].ahang;
      int oldb = flst[i].bhang;

      flst[i].ahang     = MIN(flst[i].position.bgn, flst[i].position.end) - MIN(flst[efrg].position.bgn, flst[efrg].position.end);;
      flst[i].bhang     = MAX(flst[i].position.bgn, flst[i].position.end) - MAX(flst[efrg].position.bgn, flst[efrg].position.end);;
      flst[i].contained = ((flst[i].ahang > 0) && (flst[i].bhang < 0)) ? flst[i].parent : 0;;
    
#if 0
      fprintf(stderr, "RESET extended-parent for id %d from %d,%d,%d to %d,%d,%d container=%d\n",
              flst[i].ident,
              oldp, olda, oldb, 
              flst[i].parent, flst[i].ahang, flst[i].bhang,
              flst[i].contained);
#endif

      if (flst[i].contained)
        FixContainmentRelationships(i, efrg, flst, flen);
    }
  }



  CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                          CNS_OPTIONS_MIN_ANCHOR_DEFAULT };
  int         success = FALSE;

  int32       firstFailed = 0;

  if (!success)
    success = MultiAlignUnitig(macopy, ScaffoldGraph->gkpStore, CNS_QUIET, &options, firstFailed);

  if (!success) {
    fprintf(stderr, "WARNING: MultiAlignUnitig failure on unitig %d at fragment index %d\n",
            unitig->id, firstFailed);
    DeleteMultiAlignT(macopy);  //  We own this.
    //exit(1);
    return FALSE;
  }

  //  If our new multialign differs wildly from the expected length, we assume the fragment
  //  extension was bogus -- MultiAlignUnitig, instead of aligning the extension to the unitig, just
  //  inserted gaps (where x's are gaps).
  //
  //  Before:    ----------------   After:  -------xxxx-------
  //             <------                       <---xxxx---
  //           +++++++<-------              +++xxxx++++<-------
  //
  if (GetMultiAlignLength(macopy) > GetMultiAlignLength(maorig) + 1.1 * extLength) {
    PrintMultiAlignT(stderr, macopy, ScaffoldGraph->gkpStore, 0, 1, AS_READ_CLEAR_LATEST);
    fprintf(stderr, "WARNING:  new multialign is too long -- expected %d got %d -- suspect alignment, don't use it.\n",
            GetMultiAlignLength(maorig) + (int)1.1 * extLength,
            GetMultiAlignLength(macopy));
    DeleteMultiAlignT(macopy);     //  We own this.
    return(FALSE);
  }

  assert(macopy->maID == maorig->maID);

  //  We own macopy, but release ownership to the store.
  //
  ScaffoldGraph->tigStore->insertMultiAlign(macopy, TRUE, TRUE);

  if (debug.eCRmainLV > 0) {
    fprintf(debug.eCRmainFP, "for unitig %d, after reforming, strlen(reformed_consensus) = " F_SIZE_T "\n", unitig->id, strlen(Getchar(macopy->consensus, 0)));
    fprintf(debug.eCRmainFP, "for unitig %d, after reforming, consensus:\n%s\n", unitig->id, Getchar(macopy->consensus, 0));
  }

  return(TRUE);
}




void
extendClearRange(int fragIid, int frag3pDelta) {
  uint32 clr_bgn, clr_end;
  gkFragment fr;

  if (fragIid == -1)
    return;

  if (frag3pDelta < 0)
    fprintf(stderr, "extendClearRange()-- WARNING: frag3pDelta less than zero: %d\n", frag3pDelta);
  assert(frag3pDelta >= 0);

  ScaffoldGraph->gkpStore->gkStore_getFragment(fragIid, &fr, GKFRAGMENT_INF);

  fr.gkFragment_getClearRegion(clr_bgn, clr_end, AS_READ_CLEAR_ECR_0 + iterNumber - 1);
  assert(clr_bgn < clr_end);
  clr_end += frag3pDelta;
  fr.gkFragment_setClearRegion(clr_bgn, clr_end, AS_READ_CLEAR_ECR_0 + iterNumber);

  ScaffoldGraph->gkpStore->gkStore_setFragment(&fr);

#if 0
  fprintf(stderr, "WARNING:  Set clear for frag "F_IID",%s from (%d,%d) to (%d,%d)\n",
          fr.gkFragment_getReadIID(),
          AS_UID_toString(fr.gkFragment_getReadUID()),
          clr_bgn,
          clr_end - frag3pDelta,
          clr_bgn,
          clr_end);
#endif
}


void
revertClearRange(int fragIid) {
  uint32     clr_bgn;
  uint32     clr_end;
  gkFragment fr;

  if (fragIid == -1)
    return;

  ScaffoldGraph->gkpStore->gkStore_getFragment(fragIid, &fr, GKFRAGMENT_INF);

  fr.gkFragment_getClearRegion(clr_bgn, clr_end, AS_READ_CLEAR_ECR_0 + iterNumber - 1);
  assert(clr_bgn < clr_end);
  fr.gkFragment_setClearRegion(clr_bgn, clr_end, AS_READ_CLEAR_ECR_0 + iterNumber);

  ScaffoldGraph->gkpStore->gkStore_setFragment(&fr);

#if 0
  fprintf(stderr, "WARNING:  Revert clear for frag "F_IID",%s to (%d,%d)\n",
          fr.gkFragment_getReadIID(),
          AS_UID_toString(fr.gkFragment_getReadUID()),
          clr_bgn,
          clr_end);
#endif
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
    CIFragT *leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, lFragIid);

    ScaffoldGraph->tigStore->copyMultiAlign(leftFrag->cid, TRUE, savedLeftUnitigMA);
    ScaffoldGraph->tigStore->copyMultiAlign(leftFrag->contigID, FALSE, savedLeftContigMA);
  }

  if (rFragIid != -1) {
    CIFragT *rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, rFragIid);

    ScaffoldGraph->tigStore->copyMultiAlign(rightFrag->cid, TRUE, savedRightUnitigMA);
    ScaffoldGraph->tigStore->copyMultiAlign(rightFrag->contigID, FALSE, savedRightContigMA);
  }
}

void
restoreFragAndUnitigData(int lFragIid, int rFragIid) {

  revertClearRange(lFragIid);
  revertClearRange(rFragIid);

  if (lFragIid != -1) {
    CIFragT    *leftFrag = GetCIFragT(ScaffoldGraph->CIFrags, lFragIid);
    NodeCGW_T  *unitig   = GetGraphNode(ScaffoldGraph->CIGraph, leftFrag->cid);

    ScaffoldGraph->tigStore->insertMultiAlign(savedLeftUnitigMA, TRUE, FALSE);
    ScaffoldGraph->tigStore->insertMultiAlign(savedLeftContigMA, FALSE, FALSE);

    SynchUnitigTWithMultiAlignT(unitig);
  }

  if (rFragIid != -1) {
    CIFragT    *rightFrag = GetCIFragT(ScaffoldGraph->CIFrags, rFragIid);
    NodeCGW_T  *unitig    = GetGraphNode(ScaffoldGraph->CIGraph, rightFrag->cid);

    ScaffoldGraph->tigStore->insertMultiAlign(savedRightUnitigMA, TRUE, FALSE);
    ScaffoldGraph->tigStore->insertMultiAlign(savedRightContigMA, FALSE, FALSE);

    SynchUnitigTWithMultiAlignT(unitig);
  }
}
