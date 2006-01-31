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

static const char CM_ID[] = "$Id: eCR.c,v 1.1 2006-01-31 22:03:20 brianwalenz Exp $";

#include "eCR.h"

int                     totalContigsBaseChange = 0;
ReadStructp             fsread = NULL;
VA_TYPE(char)          *lContigConsensus = NULL;
VA_TYPE(char)          *rContigConsensus = NULL;
VA_TYPE(char)          *lContigQuality = NULL;
VA_TYPE(char)          *rContigQuality = NULL;
VA_TYPE(char)          *reformed_consensus = NULL;
VA_TYPE(char)          *reformed_quality = NULL;
VA_TYPE(int32)         *reformed_deltas = NULL;
VA_TYPE(IntElementPos) *ContigPositions = NULL;
VA_TYPE(IntElementPos) *UnitigPositions = NULL;

static MultiAlignT *savedLeftUnitigMA = NULL;
static MultiAlignT *savedRightUnitigMA = NULL;
static MultiAlignT *savedLeftContigMA = NULL;
static MultiAlignT *savedRightContigMA = NULL;



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
   int32 ipos;

   IntMultiPos  *frag;
   IntUnitigPos *unitig;


   //  type should be AS_UNITIG

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





int
main(int argc, char **argv) {

  Global_CGW *data;
  char *outputPath = NULL;
  int setFragStoreName = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;

  int ckptNum = NULLINDEX;
  int backupFrgStore = TRUE;
  int i; // , index;
  int sid, startingGap = 0, setStartingGap = FALSE;
  int numExtendableGaps = 0, numGaps = 0, numGapsClosed = 0, totalBasesInClosedGaps = 0;
  int leftContigExtendable = 0, rightContigExtendable = 0, bothContigsExtendable = 0;
  int gapNumber = 0, numSmallGaps = 0, numSmallGapsClosed = 0, numLargeGaps = 0, numLargeGapsClosed = 0;
  int numSmallGapsThisScaff, numSmallGapsClosedThisScaff, numLargeGapsThisScaff, numLargeGapsClosedThisScaff;
  float maxGapSizeClosed = 0.0;
  int maxGapSizeClosedNumber = -1, numGapsVarTooSmall = 0;
  int totalOlapLength = 0, totalOlapDiffs = 0;
  float totalOlapVariance = 0.0;
  int surrogateOnLeftEnd, surrogateOnRightEnd, surrogatesOnEnd = 0;

  int scaffoldBegin = -1;
  int scaffoldEnd   = -1;

  int *closedGap, *closedGapDelta;

  int numClosingsTried = 0;
  int unitigMultiAlignFailures = 0;
  int replaceEndUnitigFailures = 0;
  int createAContigFailures = 0;
  int unitigToContigFailures = 0;
  int noOverlapFound = 0;

  int doRevertFirst = FALSE;

  double sumScaffoldLengths = 0, sumScaffoldLengthsLastCkp = 0;
  int *numGapsInScaffold, *numGapsClosedInScaffold;
  int *numSmallGapsInScaffold, *numSmallGapsClosedInScaffold;
  int *numLargeGapsInScaffold, *numLargeGapsClosedInScaffold;

  // save off whatever the rest of the world has for default values
  // for Local_Overlap_AS_forCNS
  //
  saveDefaultLocalAlignerVariables();
  
#if 0
  //  moved into examinegap
  // set some variables to control Local_Overlap_AS_forCNS
  MaxGaps        = 5;
  MaxInteriorGap = 30;
  asymmetricEnds = TRUE;
#endif

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc   = stderr;
  data->stderro   = stderr;
  data->stderrfp  = stderr;
  data->timefp    = stderr;
  data->logfp     = stderr;

  //  This is used all over the place, do not remove it!
  //
  fsread = new_ReadStruct();


#if 0
#ifdef X86_GCC_LINUX
  /*
  ** Set the x86 FPU control word to force double
  ** precision rounding rather than `extended'
  ** precision rounding. This causes base
  ** calls and quality values on x86 GCC-Linux
  ** (tested on RedHat Linux) machines to be
  ** identical to those on IEEE conforming UNIX
  ** machines.
  */
  fpu_control_t fpu_cw;

  fpu_cw = ( _FPU_DEFAULT & ~_FPU_EXTENDED ) | _FPU_DOUBLE;

  _FPU_SETCW( fpu_cw );
#endif
#endif

  {
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "Bb:c:C:e:f:g:m:n:Rs:")) != EOF)){
      switch(ch) {
        case 'B':
          backupFrgStore = FALSE;
          break;
        case 'c':
          strcpy(data->File_Name_Prefix, argv[optind - 1]);
          setPrefixName = TRUE;		  
          break;
        case 'C':
          startingGap = atoi(argv[optind - 1]);
          setStartingGap = TRUE;
          break;
        case 'f':
          strcpy(data->Frag_Store_Name, argv[optind - 1]);
          setFragStoreName = TRUE;
          break;
        case 'g':
          strcpy(data->Gatekeeper_Store_Name, argv[optind - 1]);
          setGatekeeperStore = TRUE;
          break;	  
        case 'm':
          fprintf(stderr, "option m removed.\n");
          //MaxInteriorGap = atoi(argv[optind - 1]);
          //fprintf(stderr, "setting MaxInteriorGap to %d\n", MaxInteriorGap);
          exit(1);
          break;
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case 'R':
          doRevertFirst = TRUE;
          break;
        case 'b':
          scaffoldBegin = atoi(argv[optind - 1]);
          break;
        case 'e':
          scaffoldEnd   = atoi(argv[optind - 1]);
          break;
        case 's':
          fprintf(stderr, "singleSid is no longer a valid option; use -b and -e.\n");
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if ((setPrefixName == FALSE) ||
        (setFragStoreName == 0) ||
        (setGatekeeperStore == 0)) {
      fprintf(stderr, "usage: %s [opts] -c ckpName -n ckpNumber -f frgStore -g gkpStore\n", argv[0]);
      fprintf(stderr, "  -B           DISABLE fragStore backups\n");
      fprintf(stderr, "  -c ckpName   Use ckpName as the checkpoint name\n");
      fprintf(stderr, "  -C gap#      Start at a specific gap number\n");
      fprintf(stderr, "  -f frgStore  The fragment store\n");
      fprintf(stderr, "  -g gkpStore  The gatekeeper store\n");
      fprintf(stderr, "  -n ckpNumber The checkpoint to use\n");
      fprintf(stderr, "  -R           Revert the fragment store to the consensus clear range (SLOW!)\n");
      fprintf(stderr, "  -b scafBeg   Begin at a specific scaffold\n");
      fprintf(stderr, "  -e scafEnd   End at a specific scaffold\n");
      exit(1);
    }
  }

  if (setStartingGap == TRUE)
    fprintf(stderr, "set starting gap to %d\n", startingGap);


  // check to see if db.frg.orig exists, and if not, copy db.frg to it
  if (backupFrgStore) {
    char temp_buf[1024];
    int sysReturn;
    FILE *fileExistenceCheck;
	  
    sprintf(temp_buf, "%s/db.frg.orig", GlobalData->Frag_Store_Name);
    fileExistenceCheck = fopen(temp_buf, "r");
    if (fileExistenceCheck == NULL) {
      fprintf(stderr, "file %s/db.frg.orig does not exist, creating upon start from ckpt %d.\n",
              GlobalData->Frag_Store_Name, ckptNum);
		  
      sprintf(temp_buf, "cp %s/db.frg %s/db.frg.orig", 
              GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name);
      sysReturn = system(temp_buf);
      if (sysReturn != -1)
        fprintf(stderr, "copied %s/db.frg to %s/db.frg.orig\n",
                GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name);
      else {
        fprintf(stderr, "error copying %s/db.frg to %s/db.frg.orig\n",
                GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name);
        assert(0);
      }
		  
      sprintf(temp_buf, "chmod 444 %s/db.frg.orig", GlobalData->Frag_Store_Name);
      sysReturn = system(temp_buf);
      if (sysReturn == -1) {
        fprintf(stderr, "error doing chmod on %s/db.frg.orig\n", GlobalData->Frag_Store_Name);
        assert(0);
      }
    }
  }


  //  LoadScaffoldGraphFromCheckpoint wants to CheckCIScaffoldT()
  //  which can RecomputeOffsetsInScaffold(), which can eventually,
  //  try to get an overlap.  Unless this is set, it bombs.
  //
  GlobalData->aligner=Local_Overlap_AS_forCNS;

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix, ckptNum, TRUE);

  //  The ScaffoldGraph, by default now, opens the fragstore
  //  read-only.  We reopen it read-write here.
  //
  closeFragStore(ScaffoldGraph->fragStore);
  ScaffoldGraph->fragStore = openFragStore(GlobalData->Frag_Store_Name, "r+");
  if(ScaffoldGraph->fragStore == NULLSTOREHANDLE){
    fprintf(stderr,"**** Failure to re-open frag store %s ...exiting\n", GlobalData->Frag_Store_Name);
    exit(1);
  }


  //  Intiialize the variable arrays
  //
  lContigConsensus   = CreateVA_char(1024);
  rContigConsensus   = CreateVA_char(1024);
  lContigQuality     = CreateVA_char(1024);
  rContigQuality     = CreateVA_char(1024);
  reformed_consensus = CreateVA_char(200000);
  reformed_quality   = CreateVA_char(200000);
  reformed_deltas    = CreateVA_int32(1);



  //  Revert back to the CNS clear range for all frags.  Useful if you
  //  happen to screw stuff up, but too slow to be enabled by default.
  //
  if (doRevertFirst) {
    int  firstFrag=0, lastFrag=0;
    int  fragIid=0;
    int  cnsBeg=0, cnsEnd=0, cgwBeg=0, cgwEnd=0;
    int  modified=0;

    fprintf(stderr, "Reverting back to the CNS clear range!\n");

    firstFrag = getFirstElemFragStore(ScaffoldGraph->fragStore);
    lastFrag  = getLastElemFragStore(ScaffoldGraph->fragStore) + 1;

    for (fragIid=firstFrag; fragIid<lastFrag; fragIid++) {
      getFragStore(ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	

      getClearRegion_ReadStruct(fsread, &cnsBeg, &cnsEnd, READSTRUCT_CNS);
      getClearRegion_ReadStruct(fsread, &cgwBeg, &cgwEnd, READSTRUCT_CGW);

      if ((cnsBeg != cgwBeg) || (cnsEnd != cgwEnd)) {
        fprintf(stderr, "Fragment %d modified CNS (%d,%d)  CGW (%d,%d).\n",
                fragIid, cnsBeg, cnsEnd, cgwBeg, cgwEnd);
        setClearRegion_ReadStruct(fsread, cnsBeg, cnsEnd, READSTRUCT_CGW);
        setFragStore(ScaffoldGraph->fragStore, fragIid, fsread);
        modified++;
      }
    }

    fprintf(stderr, "Reverted %d frags back to the CNS clear range.\n", modified);
    fprintf(stderr, "Bye.\n");
    exit(0);
  }



  closedGap                    = (int *) safe_malloc(GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));
  closedGapDelta               = (int *) safe_malloc(GetNumGraphNodes(ScaffoldGraph->ContigGraph) * sizeof(int));

  // following arrays are used in checkpointing
  numGapsInScaffold            = (int *) safe_malloc(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  numGapsClosedInScaffold      = (int *) safe_malloc(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  numSmallGapsInScaffold       = (int *) safe_malloc(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  numSmallGapsClosedInScaffold = (int *) safe_malloc(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  numLargeGapsInScaffold       = (int *) safe_malloc(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));
  numLargeGapsClosedInScaffold = (int *) safe_malloc(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) * sizeof(int));

  loadEcrCheckpoint(ckptNum, numGapsInScaffold, numGapsClosedInScaffold,
                    numSmallGapsInScaffold, numSmallGapsClosedInScaffold,
                    numLargeGapsInScaffold, numLargeGapsClosedInScaffold);

  for (i = 0; i < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); i++) {
    numGapsInScaffold[i]            = 0;
    numGapsClosedInScaffold[i]      = 0;
    numSmallGapsInScaffold[i]       = 0;
    numSmallGapsClosedInScaffold[i] = 0;
    numLargeGapsInScaffold[i]       = 0;
    numLargeGapsClosedInScaffold[i] = 0;
  }
  

  //
  // scan all the scaffolds
  //

  if (scaffoldBegin == -1)
    scaffoldBegin = 0;
  if (scaffoldEnd == -1)
    scaffoldEnd = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);
  if (scaffoldEnd > GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph))
    scaffoldEnd = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);

  for (sid = scaffoldBegin; sid < scaffoldEnd; sid++) {
    CIScaffoldT * scaff;
    int icnt, lFragIid, rFragIid, lextension, rextension;
    extendableFrag leftExtFragsArray[ MAX_EXTENDABLE_FRAGS ], rightExtFragsArray[ MAX_EXTENDABLE_FRAGS ];
    int numLeftFrags, numRightFrags, rcontigID,lcontigID;
    IntElementPos contigPos;
    ContigT *lcontig, *rcontig, *newContig;
    ContigT  lcontigBackup;
    ContigT  rcontigBackup;
    int lunitigID, runitigID;

    lextension = rextension = 0;
	
    scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
    // make sure the scaffold is there
    assert(scaff != NULL);
    


    // not interested in dead scaffold, not real scaffolds, or singleton scaffolds
    

    if ((isDeadCIScaffoldT(scaff)) ||
        (scaff->type != REAL_SCAFFOLD) ||
        (scaff->info.Scaffold.numElements < 2))
      continue;


    fprintf(stderr,"\n=====================================================================\n");
    fprintf(stderr,"=== examing scaffold %d, size %f\n", sid, scaff->bpLength.mean);

    numSmallGapsThisScaff = 0;
    numSmallGapsClosedThisScaff = 0;
    numLargeGapsThisScaff = 0;
    numLargeGapsClosedThisScaff = 0;

    icnt = 0;

    lcontig = GetGraphNode(ScaffoldGraph->ContigGraph, scaff->info.Scaffold.AEndCI);
    lcontigID = lcontig->id;
    rcontigID = lcontig->BEndNext;
    while (rcontigID != -1) {
      NodeOrient lcontigOrientation, rcontigOrientation;
      LengthT gapSize, newOffsetAEnd, newOffsetBEnd;
      double maxRContigOffset;
      int32 nextContigIndex;

#ifdef DEBUG_ECR
      fprintf(stderr, "at top of loop: lcontig->BEndNext: %d\n", lcontig->BEndNext);
#endif
      rcontig = GetGraphNode(ScaffoldGraph->ContigGraph, rcontigID);

      // lcontigIdGap[ gapNumber ] = lcontig->id;
      // rcontigIdGap[ gapNumber ] = rcontig->id;	  
      // lcontigLength[ gapNumber ] = (int) lcontig->bpLength.mean;
      // rcontigLength[ gapNumber ] = (int) rcontig->bpLength.mean;

      assert(lcontig != NULL);
      assert(rcontig != NULL);		

      surrogatesOnEnd = 0;
      gapSize = FindGapLength(lcontig, rcontig, FALSE);
      // originalGaps[numGaps] = (int) gapSize.mean;
      numGaps++;

      if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
        lcontigOrientation = A_B;
      else
        lcontigOrientation = B_A;

      if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
        rcontigOrientation = A_B;
      else
        rcontigOrientation = B_A;

      fprintf(stderr, "\n\n\n---------------------------------------------------------------\n");
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


      fprintf(stderr, "\nexamining lcontig %d (orientation %c) \n", lcontig->id, lcontigOrientation);
      lFragIid = - 1;
      numLeftFrags = numRightFrags = 0;

      // find the extreme read on the correct end of the lcontig
      if (lcontigOrientation == A_B) {
        numLeftFrags = findLastExtendableFrags(lcontig, leftExtFragsArray);
        surrogateOnLeftEnd = findLastUnitig(lcontig, &lunitigID);
      } else {
        numLeftFrags = findFirstExtendableFrags(lcontig, leftExtFragsArray);
        surrogateOnLeftEnd = findFirstUnitig(lcontig, &lunitigID);
      }
      surrogatesOnEnd += surrogateOnLeftEnd;
      if (surrogateOnLeftEnd)
        numLeftFrags = 0;
      fprintf(stderr, "finished examining lcontig %d (orientation %c) \n", lcontig->id, lcontigOrientation);

	  
      fprintf(stderr, "\nexamining rcontig %d (orientation %c) \n", rcontig->id, rcontigOrientation);
      rFragIid = - 1;
      // find the extreme read on the correct end of the rchunk
      if (rcontigOrientation == A_B) {
        numRightFrags = findFirstExtendableFrags(rcontig, rightExtFragsArray);
        surrogateOnRightEnd = findFirstUnitig(rcontig, &runitigID);
      } else {
        numRightFrags = findLastExtendableFrags(rcontig, rightExtFragsArray);
        surrogateOnRightEnd = findLastUnitig(rcontig, &runitigID);
      }
      surrogatesOnEnd += surrogateOnRightEnd;
      if (surrogateOnRightEnd)
        numRightFrags = 0;
      fprintf(stderr, "finished examining rcontig %d (orientation %c) \n", rcontig->id, rcontigOrientation);


      if (setStartingGap == TRUE && gapNumber < startingGap)
        numLeftFrags = numRightFrags = 0;

      //  This is a good place to check rcontig->id and abort if that gap is causing problems.
      //
      //if (rcontig->id == 405527)
      //  numLeftFrags = numRightFrags = 0;


      if (numLeftFrags > 0)
        leftContigExtendable++;
      if (numRightFrags > 0)
        rightContigExtendable++;
      if (numLeftFrags > 0 && numRightFrags > 0)
        bothContigsExtendable++;

      {
        // this test should be applied after the extensions have been factored into gap size
        // original:
        // if (gapSize.mean - lextension - rextension > NUM_STDDEV_CUTOFF * sqrt(gapSize.variance) && gapSize.mean > 100.0)
        // hacking
        // compute here using {left. right}ExtFragsArray[ index ].extension, not lextension and rextension

        int ahang, currLength, bhang, lcontigBasesIntact, rcontigBasesIntact;
        int currDiffs;
        int leftFragIndex, rightFragIndex;
        int leftFragFlapLength, rightFragFlapLength;
        int gotNewLeftMA, gotNewRightMA;
        InfoByIID *info;
		
        fprintf(stderr, "gap (%d, %d) has extendable frag pair on left: %d, on right: %d\n", 
                lcontig->id, rcontig->id, numLeftFrags, numRightFrags);
        numExtendableGaps++;
		
        // set an extra member of the arrays for when the contig is not being extended
        leftExtFragsArray[ numLeftFrags ].fragIid = -1;
        leftExtFragsArray[ numLeftFrags++ ].extension = 0;
        rightExtFragsArray[ numRightFrags ].fragIid = -1;		  
        rightExtFragsArray[ numRightFrags++ ].extension = 0;
        closedGap[ gapNumber ] = FALSE;

        //  Incase we need to back out changes late in the extension, we save copies of the two contigs.
        //
        memcpy(&lcontigBackup, lcontig, sizeof(ContigT));
        memcpy(&rcontigBackup, rcontig, sizeof(ContigT));

			
        for (leftFragIndex = 0; leftFragIndex < numLeftFrags && closedGap[ gapNumber ] == FALSE; leftFragIndex++) {
          for (rightFragIndex = 0; rightFragIndex < numRightFrags && closedGap[ gapNumber ] == FALSE; rightFragIndex++) {
            fprintf(stderr, "examining frags %d and %d\n", leftExtFragsArray[ leftFragIndex ].fragIid, 
                    rightExtFragsArray[ rightFragIndex ].fragIid);


            if (gapSize.mean - leftExtFragsArray[ leftFragIndex ].extension - 
                rightExtFragsArray[ rightFragIndex ].extension > 
                NUM_STDDEV_CUTOFF * sqrt(gapSize.variance) && gapSize.mean > 100.0) {
              fprintf(stderr, "leftExtFragsArray[ %d ].extension: %10d, rightExtFragsArray[ %d ].extension: %10d\n",
                      leftFragIndex, leftExtFragsArray[ leftFragIndex ].extension, 
                      rightFragIndex, rightExtFragsArray[ rightFragIndex ].extension);
              // numGapsVarTooSmall++;
              closedGap[ gapNumber ] = FALSE;
              fprintf(stderr, "gap variance too large (gapSize - extensions: %.2f, %.1f * sqrt(gapSize.variance): %.2f\n",
                      gapSize.mean - leftExtFragsArray[ leftFragIndex ].extension - 
                      rightExtFragsArray[ rightFragIndex ].extension, 
                      NUM_STDDEV_CUTOFF, NUM_STDDEV_CUTOFF * sqrt(gapSize.variance));
              continue;
            }

            lFragIid = leftExtFragsArray[ leftFragIndex ].fragIid;
            rFragIid = rightExtFragsArray[ rightFragIndex ].fragIid;

            // have to check and make sure that the frags belong to the correct unitig
            if (lFragIid != -1) {
              info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);
              assert(info->set);
              if (GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex)->cid != lunitigID)
                continue;
            }
			
            if (rFragIid != -1) {
              info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
              assert(info->set);
              if (GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex)->cid != runitigID)
                continue;
            }

            numClosingsTried++;
			  
            // dumpContigInfo(lcontig);
            // dumpContigInfo(rcontig);

            if (examineGap(lcontig, lFragIid, rcontig, rFragIid, 
                           gapNumber, &ahang, &currLength, &bhang, &currDiffs,
                           &lcontigBasesIntact, &rcontigBasesIntact, &closedGapDelta[gapNumber],
                           leftExtFragsArray[ leftFragIndex ].basesToNextFrag, 
                           rightExtFragsArray[ rightFragIndex ].basesToNextFrag,
                           &leftFragFlapLength, &rightFragFlapLength)) {
              int keepGap = TRUE;

              // save off copies of everything we might alter so we can restore if gap closing fails
              saveFragAndUnitigData(lFragIid, rFragIid);

              totalBasesInClosedGaps += (int) gapSize.mean;
				
              if (CONTIG_BASES < 2000) {

                // these checks Granger suggested
                if (ahang + currLength + bhang - 1000 > min(CONTIG_BASES, (int) lcontig->bpLength.mean) +
                    min(CONTIG_BASES, (int) rcontig->bpLength.mean)) {
                  fprintf(stderr, "at gapNumber %d, ahang + currLength + bhang - 1000 = %d\n",
                          gapNumber, ahang + currLength + bhang - 1000);
                  fprintf(stderr, "at gapNumber %d, back(A) + back(B) = %d\n",
                          gapNumber, min(CONTIG_BASES, (int) lcontig->bpLength.mean) +
                          min(CONTIG_BASES, (int) rcontig->bpLength.mean));
                  keepGap = FALSE;
                }
				
                if (ahang + currLength + bhang + 700 < min(CONTIG_BASES, (int) lcontig->bpLength.mean) +
                    min(CONTIG_BASES, (int) rcontig->bpLength.mean)) {
                  fprintf(stderr, "at gapNumber %d, ahang + currLength + bhang + 500 = %d\n",
                          gapNumber, ahang + currLength + bhang + 500);
                  fprintf(stderr, "at gapNumber %d, back(A) + back(B) = %d\n",
                          gapNumber, min(CONTIG_BASES, (int) lcontig->bpLength.mean) +
                          min(CONTIG_BASES, (int) rcontig->bpLength.mean));
                  keepGap = FALSE;
                }
              }

              gotNewLeftMA = gotNewRightMA = TRUE;

              // extend the clear ranges of the frags
              //
              if (keepGap) {

                //  The previous version merged the ExtFragsArray[]
                //  test with the extendCgwClearRange() -- so if the
                //  left frag updated, but the right frag failed, we
                //  might leave things inconsistent.  Hopefully, we
                //  restore the original clear ranges at the end of
                //  all this.  But we still fix it to not modify
                //  unless both frags are OK.

                if (lFragIid != -1)
                  if (leftExtFragsArray[ leftFragIndex ].addedBases - leftFragFlapLength < 0)
                    keepGap = FALSE;

                if (rFragIid != -1)
                  if (rightExtFragsArray[ rightFragIndex ].addedBases - rightFragFlapLength < 0)
                    keepGap = FALSE;

                if (keepGap) {
                  if (lFragIid != -1) {
                    fprintf(stderr,"adjusting left frg clear range by %d - %d = %d bases\n",
                            leftExtFragsArray[ leftFragIndex ].addedBases,leftFragFlapLength,
                            leftExtFragsArray[ leftFragIndex ].addedBases - leftFragFlapLength);
                    extendCgwClearRange(lFragIid,
                                        leftExtFragsArray[leftFragIndex].addedBases - leftFragFlapLength);
                  }

                  if (rFragIid != -1) {
                    fprintf(stderr,"adjusting right frg clear range by %d - %d = %d bases\n",
                            rightExtFragsArray[ rightFragIndex ].addedBases,rightFragFlapLength,
                            rightExtFragsArray[ rightFragIndex ].addedBases - rightFragFlapLength);
                    extendCgwClearRange(rFragIid,
                                        rightExtFragsArray[rightFragIndex].addedBases - rightFragFlapLength); 
                  }
                }
              }


              if (keepGap) { // the fragment extensions have succeeded
                // InfoByIID *info;
                CIFragT *frag;
                NodeCGW_T *unitig;
                MultiAlignT *new_cma;
                fragPositions *fragPoss;
				  
                DumpContigMultiAlignInfo ("before anything (left)", NULL, lcontig->id);
                DumpContigUngappedOffsets("before anything (left)", lcontig->id);
                DumpContigMultiAlignInfo ("before anything (right)", NULL, rcontig->id);
                DumpContigUngappedOffsets("before anything (right)", rcontig->id);

                fprintf(stderr, "before altering, lctg: %12.0f, %12.0f\n",
                        (lcontigOrientation == A_B) ? lcontig->offsetAEnd.mean : lcontig->offsetBEnd.mean,
                        (lcontigOrientation == A_B) ? lcontig->offsetBEnd.mean : lcontig->offsetAEnd.mean);
                fprintf(stderr, "                 rctg: %12.0f, %12.0f\n\n",
                        (rcontigOrientation == A_B) ? rcontig->offsetAEnd.mean : rcontig->offsetBEnd.mean,
                        (rcontigOrientation == A_B) ? rcontig->offsetBEnd.mean : rcontig->offsetAEnd.mean);

                // save the max offset of the right contig so we know
                // how to adjust the offsets of the contigs further
                // along the scaffold later

                maxRContigOffset = max(rcontig->offsetAEnd.mean, rcontig->offsetBEnd.mean);
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

                  info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, lFragIid);
                  assert(info->set);
                  frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);

                  unitig = GetGraphNode(ScaffoldGraph->CIGraph, frag->CIid);
                  //extendUnitigs(unitig, lFragIid, leftExtFragsArray[ leftFragIndex ], TRUE);
                  getAlteredFragPositions(unitig, &fragPoss, lFragIid, 
                                          leftExtFragsArray[ leftFragIndex ].extension - leftFragFlapLength);

                  gotNewLeftMA = GetNewUnitigMultiAlign(unitig, fragPoss, lFragIid);

                  free(fragPoss);  // inefficient, let's just do a big array once that get's reused

                  if (!gotNewLeftMA)
                    unitigMultiAlignFailures++;

                  if (gotNewLeftMA) {
                    int extendToLeft;
                    SynchUnitigTWithMultiAlignT(unitig);

                    DumpContigMultiAlignInfo ("before ReplaceEndUnitigInContig (left)", NULL, lcontig->id);

                    if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
                      extendToLeft = FALSE;
                    else
                      extendToLeft = TRUE;

                    new_cma = ReplaceEndUnitigInContig(ScaffoldGraph->sequenceDB,
                                                       ScaffoldGraph->fragStore,
                                                       lcontig->id, unitig->id, extendToLeft,
                                                       GlobalData->aligner,
                                                       NULL);

                    if (new_cma) {
                      DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (left)", NULL, lcontig->id);
                      DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (left)", new_cma, lcontig->id);

                      UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, lcontig->id, FALSE);
                      InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, lcontig->id, 
                                                    FALSE, new_cma, TRUE);

                      fprintf(stderr, "strlen(Getchar (new_cma->consensus)): " F_SIZE_T "\n",
                              strlen(Getchar (new_cma->consensus, 0)));

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
                  info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, rFragIid);
                  assert(info->set);
                  frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
					
                  unitig = GetGraphNode(ScaffoldGraph->CIGraph, frag->CIid);
                  //extendUnitigs(unitig, rFragIid, rightExtFragsArray[ rightFragIndex ], FALSE);
					
                  getAlteredFragPositions(unitig, &fragPoss, rFragIid, 
                                          rightExtFragsArray[ rightFragIndex ].extension- rightFragFlapLength);

                  gotNewRightMA = GetNewUnitigMultiAlign(unitig, fragPoss, rFragIid);

                  free(fragPoss);  // inefficient, let's just do a big array once that get's reused

                  if (!gotNewRightMA)
                    unitigMultiAlignFailures++;

                  if (gotNewRightMA) {
                    int extendToLeft;
                    SynchUnitigTWithMultiAlignT(unitig);

                    DumpContigMultiAlignInfo ("before ReplaceEndUnitigInContig (right)", NULL, rcontig->id);

                    if (rcontig->offsetAEnd.mean < rcontig->offsetBEnd.mean)
                      extendToLeft = TRUE;
                    else
                      extendToLeft = FALSE;

                    new_cma = ReplaceEndUnitigInContig(ScaffoldGraph->sequenceDB,
                                                       ScaffoldGraph->fragStore,
                                                       rcontig->id, unitig->id, extendToLeft,
                                                       GlobalData->aligner,
                                                       NULL);

                    if (new_cma) {
                      DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (right) (original contig)", NULL, rcontig->id);
                      DumpContigMultiAlignInfo("after ReplaceEndUnitigInContig (right) (new contig)", new_cma, rcontig->id);

                      UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, rcontig->id, FALSE);
                      InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, rcontig->id, 
                                                    FALSE, new_cma, TRUE);
                      fprintf(stderr, "strlen(Getchar (new_cma->consensus)): " F_SIZE_T "\n",
                              strlen(Getchar (new_cma->consensus, 0)));

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
                rcontig->offsetAEnd.mean += closedGapDelta[gapNumber];
                rcontig->offsetBEnd.mean += closedGapDelta[gapNumber];				  

                fprintf(stderr, "after altering, lctg: %12.0f, %12.0f\n",
                        (lcontigOrientation == A_B) ? lcontig->offsetAEnd.mean : lcontig->offsetBEnd.mean,
                        (lcontigOrientation == A_B) ? lcontig->offsetBEnd.mean : lcontig->offsetAEnd.mean);
                fprintf(stderr, "                rctg: %12.0f, %12.0f\n\n",
                        (rcontigOrientation == A_B) ? rcontig->offsetAEnd.mean : rcontig->offsetBEnd.mean,
                        (rcontigOrientation == A_B) ? rcontig->offsetBEnd.mean : rcontig->offsetAEnd.mean);

                /* temp hack!!*/  // if (lFragIid == -1 && rFragIid == -1)
                // this effectively kills the contigs and rebuilds from the unitigs
					
                // setup for contig merge
                if (ContigPositions == NULL)
                  ContigPositions = CreateVA_IntElementPos(2);
                ResetVA_IntElementPos(ContigPositions);
				  
                if (lcontig->offsetAEnd.mean < lcontig->offsetBEnd.mean)
                  delta = lcontig->offsetAEnd.mean;
                else
                  delta = lcontig->offsetBEnd.mean;

                contigPos.ident = lcontig->id;
                contigPos.type = AS_CONTIG;
                contigPos.position.bgn = lcontig->offsetAEnd.mean - delta;
                contigPos.position.end = lcontig->offsetBEnd.mean - delta;
                AppendIntElementPos(ContigPositions, &contigPos);
				  
                fprintf(stderr, "lcontig %8d positioned at %8d, %8d\n", 
                        lcontig->id,contigPos.position.bgn, contigPos.position.end);
				  
                contigPos.ident = rcontig->id;
                contigPos.type = AS_CONTIG;
                contigPos.position.bgn = rcontig->offsetAEnd.mean - delta;
                contigPos.position.end = rcontig->offsetBEnd.mean - delta;
                AppendIntElementPos(ContigPositions, &contigPos);
				  
                fprintf(stderr, "rcontig %8d positioned at %8d, %8d\n", 
                        rcontig->id,contigPos.position.bgn, contigPos.position.end);

                if (lcontigOrientation == A_B)
                  newOffsetAEnd.mean = lcontig->offsetAEnd.mean;
                else
                  newOffsetAEnd.mean = lcontig->offsetBEnd.mean;
				  
                if (rcontigOrientation == A_B)
                  newOffsetBEnd.mean = rcontig->offsetBEnd.mean;
                else
                  newOffsetBEnd.mean = rcontig->offsetAEnd.mean;

                DumpContigMultiAlignInfo ("before CreateAContigInScaffold", NULL, lcontig->id);
                DumpContigMultiAlignInfo ("before CreateAContigInScaffold", NULL, rcontig->id);
                //DumpContigUngappedOffsets(rcontig->id);

                fprintf(stderr, "CreateAContigInScaffold()-- newOffsetAEnd=%d newOffsetBEnd=%d\n",
                        newOffsetAEnd, newOffsetBEnd);

                // have to call this routine with normalized positions

                if (CheckNewUnitigMultiAlign(scaff, ContigPositions)) {
                  fprintf(stderr, "CheckNewUnitigMultiAlign()-- The new unitig multialignment is messed up, will not close this gap.\n");
                  unitigToContigFailures++;
                  success = NULL;
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
					
                  closedGap[ gapNumber ] = TRUE;
					
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
					
                  // newContig takes the place of what was the right contig
                  rcontig = newContig = GetGraphNode(ScaffoldGraph->ContigGraph, 
                                                     GetNumGraphNodes(ScaffoldGraph->ContigGraph) - 1);
                  //adjustUnitigCoords(newContig);

                  // now we need to adjust contigs past the new contig, if not on end
                  if (nextContigIndex != NULLINDEX) {
                    LengthT scaffoldDelta;
					  
                    scaffoldDelta.mean = max(newContig->offsetAEnd.mean, newContig->offsetBEnd.mean) 
                      - maxRContigOffset;
                    scaffoldDelta.variance = ComputeFudgeVariance(scaffoldDelta.mean);
                    AddDeltaToScaffoldOffsets(ScaffoldGraph, scaff->id, nextContigIndex,
                                              TRUE, FALSE, scaffoldDelta);
                  }
                }  //  keep gap
              }  //  keep gap

              // if we didn't close gap for whatever reason undo all the frag and unitig changes
              if (keepGap == FALSE) {
                fprintf(stderr, "did not close gap %8d, contigs %8d and %8d\n",
                        gapNumber, lcontig->id, rcontig->id);
                closedGap[ gapNumber ] = FALSE;
                revertToCnsClearRange(lFragIid);
                revertToCnsClearRange(rFragIid);				  
                restoreFragAndUnitigData(lFragIid, rFragIid);

                memcpy(lcontig, &lcontigBackup, sizeof(ContigT));
                memcpy(rcontig, &rcontigBackup, sizeof(ContigT));
              }
              //  end of if (examineGap())
            } else {
              fprintf(stderr, "did not close gap %8d, contigs %8d and %8d\n",
                      gapNumber, lcontig->id, rcontig->id);
              closedGap[ gapNumber ] = FALSE;
              noOverlapFound++;
            }
          }
        }
      }

      fprintf(stderr, "after gapNumber %d:\n", gapNumber);
      fprintf(stderr, "             numSmallGaps: %d (closed %d, %.2f%%)\n", 
              numSmallGaps, numSmallGapsClosed, 100.0 * numSmallGapsClosed / numSmallGaps);
      fprintf(stderr, "             numLargeGaps: %d (closed %d, %.2f%%)\n", 
              numLargeGaps, numLargeGapsClosed, 100.0 * numLargeGapsClosed / numLargeGaps);
      fprintf(stderr, "                  allGaps: %d (closed %d, %.2f%%)\n", 
              numSmallGaps + numLargeGaps, numSmallGapsClosed + numLargeGapsClosed, 
              100.0 * (numSmallGapsClosed + numLargeGapsClosed) / (numSmallGaps + numLargeGaps));

      gapNumber++;
      lcontig = rcontig;
      lcontigID = lcontig->id;
      rcontigID = lcontig->BEndNext;

      fprintf(stderr, "at bottom of loop: lcontig->BEndNext: %d\n", lcontig->BEndNext);
    }
	
    // checkpointing
    sumScaffoldLengths += scaff->bpLength.mean;

    fprintf(stderr, "sumScaffoldLengths: %f, sumScaffoldLengthsLastCkp: %f\n",
            sumScaffoldLengths, sumScaffoldLengthsLastCkp);
    fprintf(stderr, "scaffold stats, scaff %10d, smallGaps %8d closed %8d, largeGaps %8d closed %8d\n",
            scaff->id, numSmallGapsThisScaff, numSmallGapsClosedThisScaff,
            numLargeGapsThisScaff, numLargeGapsClosedThisScaff);

    numSmallGapsInScaffold[sid] = numSmallGapsThisScaff;
    numSmallGapsClosedInScaffold[sid] = numSmallGapsClosedThisScaff;
    numLargeGapsInScaffold[sid] = numLargeGapsThisScaff;
    numLargeGapsClosedInScaffold[sid] = numLargeGapsClosedThisScaff;
    numGapsInScaffold[sid] = numSmallGapsInScaffold[sid] + numLargeGapsInScaffold[sid];
    numGapsClosedInScaffold[sid] = numSmallGapsClosedInScaffold[sid] + numLargeGapsClosedInScaffold[sid];	


    if(numSmallGapsClosedInScaffold[sid]+numLargeGapsClosedInScaffold[sid]>0){
      int status = RECOMPUTE_SINGULAR;
      int recomputeIteration = 0;
      while(recomputeIteration++ < 3 &&
            (status == RECOMPUTE_SINGULAR ||
             status == RECOMPUTE_CONTIGGED_CONTAINMENTS)) {
        // need to make sure scaffold is connected with trusted raw edges
        MarkInternalEdgeStatus(ScaffoldGraph,
                               GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                            sid),
                               PAIRWISECHI2THRESHOLD_CGW,
                               1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
                               TRUE, TRUE, 0, TRUE);

        assert(IsScaffoldInternallyConnected(ScaffoldGraph,
                                             GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                          sid),
                                             ALL_EDGES));
	      
        status =
          RecomputeOffsetsInScaffold(ScaffoldGraph,
                                     GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                  sid),
                                     TRUE, TRUE, FALSE);
      }
    }

    if (sumScaffoldLengths - sumScaffoldLengthsLastCkp > 90000000) {
      sumScaffoldLengthsLastCkp = sumScaffoldLengths;
      fprintf(stderr, "checkpoint %d written during extendClearRanges, sumScaffoldLengths: %f\n",
              ScaffoldGraph->checkPointIteration, sumScaffoldLengths);
      fprintf(GlobalData->timefp, "checkpoint %d written during extendClearRanges, sumScaffoldLengths: %f\n",
              ScaffoldGraph->checkPointIteration, sumScaffoldLengths);
      CheckpointScaffoldGraph(ScaffoldGraph, 1);
      if (backupFrgStore) {
        char temp_buf[1024];
        int sysReturn;
        sprintf(temp_buf, "cp %s/db.frg %s/db.frg.%d", 
                GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration - 1);
        sysReturn = system(temp_buf);
        if (sysReturn != -1)
          fprintf(stderr, "copied frgStore to %s/db.frg.%d\n",
                  GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration - 1);
        else
          fprintf(stderr, "error encountered copying frgStore to %s/db.frg.%d\n",
                  GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration);		  
      }
      writeEcrCheckpoint(numGapsInScaffold, numGapsClosedInScaffold,
                         numSmallGapsInScaffold, numSmallGapsClosedInScaffold,
                         numLargeGapsInScaffold, numLargeGapsClosedInScaffold);
    }
  }

  // Variance = mean(x^2) - (mean(x))^2
  if (numGapsClosed > 0)
    totalOlapVariance = (totalOlapVariance / numGapsClosed) - 
      ((float) totalOlapLength / numGapsClosed) * ((float) totalOlapLength / numGapsClosed);
  else
    totalOlapVariance = -0.0;

  fprintf(stderr, "\n");
  fprintf(stderr, "                  numGaps: %d\n", numGaps);
  fprintf(stderr, "        numExtendableGaps: %d (left: %d, right: %d, both: %d)\n", 
          numExtendableGaps, leftContigExtendable, rightContigExtendable, bothContigsExtendable);
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

  fprintf(GlobalData->timefp, "checkpoint %d written at end of extendClearRanges, sumScaffoldLengths: %f",
          ScaffoldGraph->checkPointIteration, sumScaffoldLengths);
  CheckpointScaffoldGraph(ScaffoldGraph, 2);

  if (backupFrgStore) {
    char temp_buf[1024];
    int sysReturn;
    sprintf(temp_buf, "cp %s/db.frg %s/db.frg.ecr.%d", 
            GlobalData->Frag_Store_Name, GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration - 1);
    sysReturn = system(temp_buf);
    if (sysReturn != -1)
      fprintf(stderr, "copied frgStore to %s/db.frg.%d\n",
              GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration - 1);
    else
      fprintf(stderr, "error encountered copying frgStore to %s/db.frg.%d\n",
              GlobalData->Frag_Store_Name, ScaffoldGraph->checkPointIteration);		  
  }
  writeEcrCheckpoint(numGapsInScaffold, numGapsClosedInScaffold,
                     numSmallGapsInScaffold, numSmallGapsClosedInScaffold,
                     numLargeGapsInScaffold, numLargeGapsClosedInScaffold);
  
  exit(0);
}





// since we mucked with the unitigs multialignment, reset the offsets
// and bpLengths of all the unitigs in the contig the ones in the
// ScaffoldGraph are not valid anymore
//
//  Technically dead, but might be useful later
//
#if 0
void
adjustUnitigCoords(NodeCGW_T *contig) {
  MultiAlignT *ma;
  int i;
  
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
  
  for (i = 0; i < GetNumIntUnitigPoss(ma->u_list); i++) {
    IntUnitigPos *pos = GetIntUnitigPos(ma->u_list, i);
    NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, pos->ident);
    MultiAlignT *uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);

    fprintf(stderr, "before unitig %8d, bgn: %10d, end: %10d, length: %10d\n", 
            unitig->id, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end));
	
    fprintf(stderr, "in adjustUnitigCoords, for unitig %d strlen(ma->consensus) = " F_SIZE_T "\n",
            unitig->id, strlen(Getchar(uma->consensus, 0)));

    unitig->bpLength.mean = pos->position.end - pos->position.bgn;      // set length
    if (unitig->offsetAEnd.mean < unitig->offsetBEnd.mean) {
      // ordering info is okay
        unitig->offsetAEnd.mean = pos->position.bgn;
        unitig->offsetBEnd.mean = pos->position.end;
      } else {
        unitig->offsetAEnd.mean = pos->position.end;
        unitig->offsetBEnd.mean = pos->position.bgn;
      }

    fprintf(stderr, " after unitig %8d, bgn: %10d, end: %10d, length: %10d\n", 
            unitig->id, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end));
  }
}
#endif




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

  fprintf(stderr, "in findFirstExtendableFrags, firstUnitigID: %d\n", firstUnitigID);
  
  // foundFrag = FALSE;
  // *extensionOut = 0;
  for (i = 0; i < numFrags; i++) {
    mp = GetIntMultiPos(ma->f_list, i);
    frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) mp->source);

    if (frag->contigOffset3p.mean < 100.0 &&      // frag is within a cutoff of the low end of the contig
        frag->locale == -1 &&                     // and is a read
        frag->contigOffset3p.mean < frag->contigOffset5p.mean &&  // and points in the right direction
        frag->cid == firstUnitigID) { // and is in the first unitig
        char seqbuffer[AS_BACTIG_MAX_LEN+1], qltbuffer[AS_BACTIG_MAX_LEN+1];
        unsigned int clr_bgn, clr_end;
        int frag3pExtra, extension;
	  
        getFragStore(ScaffoldGraph->fragStore, frag->iid, FRAG_S_ALL, fsread);
        getClearRegion_ReadStruct(fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
        getSequence_ReadStruct(fsread, seqbuffer, qltbuffer, AS_BACTIG_MAX_LEN);


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

        fprintf(stderr, "contig->bpLength.mean: %f\n", contig->bpLength.mean);
        fprintf(stderr, "frag iid: %d, frag->contigOffset5p.mean: %f, frag->contigOffset3p.mean: %f\n",
                frag->iid, frag->contigOffset5p.mean, frag->contigOffset3p.mean);
        fprintf(stderr, "frag length: " F_SIZE_T ", 3p past clr_end length: " F_SIZE_T "\n", strlen(seqbuffer), 
                strlen(seqbuffer) - clr_end);
        fprintf(stderr, "extension: " F_SIZE_T "\n", strlen(seqbuffer) - clr_end - (int) frag->contigOffset3p.mean);
	  
        frag3pExtra = strlen(seqbuffer) - clr_end;
        extension = frag3pExtra - frag->contigOffset3p.mean;

        // ask Granger what min extension we should accept
        if (extension > 30) {
          // foundFrag = TRUE;

          extFragsArray[ extendableFragCount ].fragIid = frag->iid;
          extFragsArray[ extendableFragCount ].extension = extension;
          extFragsArray[ extendableFragCount ].addedBases = frag3pExtra;

          fprintf(stderr, "for frag %d, extension: %8d, frag3pExtra: %8d\n",
                  frag->iid, extension, frag3pExtra);

          if (frag->contigOffset3p.mean == 0)
            extFragsArray[ extendableFragCount ].fragOnEnd = TRUE;
          else
            extFragsArray[ extendableFragCount ].fragOnEnd = FALSE;

#ifdef DEBUG_ECR
              fprintf(stderr, "in contig %d, frag %d is at %f -> %f (5p->3p) \n", 
                      contig->id, frag->iid,
                      frag->contigOffset5p.mean, frag->contigOffset3p.mean);
              fprintf(stderr, "extension ratio: %.2f\n", extension / (float) (1.0 + frag3pExtra - extension));
#endif

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
  qsort(extFragsArray, extendableFragCount, sizeof(extendableFrag), &compExtendableFrags);

  fprintf(stderr, "extendableFragCount: %d\n", extendableFragCount);
  for (i = 0; i < extendableFragCount; i++) {
    if (extFragsArray[i].fragOnEnd == TRUE)
      extFragsArray[i].basesToNextFrag = secondFragStart;
    else
      extFragsArray[i].basesToNextFrag = 0;

    fprintf(stderr, "contig %8d, frag %8d can extend %8d bases into the gap\n",
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
  float maxContigPos;
  int secondFragEnd = 0;
  int extendableFragCount = 0;
  // int basesToNextFrag;
  
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numFrags = GetNumIntMultiPoss(ma->f_list);
  lastUnitigID = GetIntUnitigPos(ma->u_list, GetNumIntUnitigPoss(ma->u_list) - 1)->ident;  
  maxContigPos = contig->bpLength.mean - 1.0;

  fprintf(stderr, "in FindLastExtendableFrags, lastUnitigID: %d\n", lastUnitigID);
  
  // *extensionOut = 0;
  for (i = 0; i < numFrags; i++) {
    mp = GetIntMultiPos(ma->f_list, i);
    frag = GetCIFragT(ScaffoldGraph->CIFrags, (int32) mp->source);

    // if (frag->contigOffset3p.mean == maxContigPos && frag->locale == -1)
    if (frag->contigOffset3p.mean > maxContigPos - 100.0 &&      // frag is within a cutoff of the high end of the contig
        frag->locale == -1 &&                                    // and is a read
        frag->contigOffset5p.mean < frag->contigOffset3p.mean && // and points in the right direction
        frag->cid == lastUnitigID) {                             // and is in the last unitig
        char seqbuffer[AS_BACTIG_MAX_LEN+1], qltbuffer[AS_BACTIG_MAX_LEN+1];
        unsigned int clr_bgn, clr_end;
        int frag3pExtra, extension;

        getFragStore(ScaffoldGraph->fragStore, frag->iid, FRAG_S_ALL, fsread);
        getClearRegion_ReadStruct(fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
        getSequence_ReadStruct(fsread, seqbuffer, qltbuffer, AS_BACTIG_MAX_LEN);

        //    contig ----------------------------------------------------------------------------------->
        //                                  5p -------|---------------------------------------------|------------> 3p 
        //                                         clr_bgn                                       clr_end
        //                                                                                               |-------|
        //                                                                                             extension
        frag3pExtra = strlen(seqbuffer) - clr_end;
        extension = frag3pExtra - (int) (contig->bpLength.mean - frag->contigOffset3p.mean);

        fprintf(stderr, "contig->bpLength.mean: %f\n", contig->bpLength.mean);
        fprintf(stderr, "frag iid: %d, frag->contigOffset5p.mean: %f, frag->contigOffset3p.mean: %f\n",
                frag->iid, frag->contigOffset5p.mean, frag->contigOffset3p.mean);
        fprintf(stderr, "frag length: " F_SIZE_T ", 3p past clr_end length: %d\n", strlen(seqbuffer), frag3pExtra);
        fprintf(stderr, "extension: %d\n", extension);
	  
        if (extension > 30) {
          // foundFrag = TRUE;

          extFragsArray[ extendableFragCount ].fragIid = frag->iid;
          extFragsArray[ extendableFragCount ].extension = extension;
          extFragsArray[ extendableFragCount ].addedBases = frag3pExtra;

          fprintf(stderr, "for frag %d, extension: %8d, frag3pExtra: %8d\n",
                  frag->iid, extension, frag3pExtra);

          if (frag->contigOffset3p.mean == contig->bpLength.mean)
            extFragsArray[ extendableFragCount ].fragOnEnd = TRUE;
          else
            extFragsArray[ extendableFragCount ].fragOnEnd = FALSE;

#ifdef DEBUG_ECR
              fprintf(stderr, "in contig %d, frag %d is at %f -> %f (5p->3p) maxContigPos: %f\n", 
                      contig->id, frag->iid,
                      frag->contigOffset5p.mean, frag->contigOffset3p.mean, maxContigPos);
              fprintf(stderr, "extension ratio: %.2f\n", extension / (float) (1.0 + frag3pExtra - extension));
#endif

          extendableFragCount++;
          if (extendableFragCount > MAX_EXTENDABLE_FRAGS)
            {
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
  qsort(extFragsArray, extendableFragCount, sizeof(extendableFrag), &compExtendableFrags);

  fprintf(stderr, "extendableFragCount: %d\n", extendableFragCount);
  for (i = 0; i < extendableFragCount; i++) {
    if (extFragsArray[i].fragOnEnd == TRUE)
      extFragsArray[i].basesToNextFrag = maxContigPos - secondFragEnd;
    else
      extFragsArray[i].basesToNextFrag = 0;	  
	
    fprintf(stderr, "contig %8d, frag %8d can extend %8d bases into the gap\n",
            contig->id, extFragsArray[i].fragIid, extFragsArray[i].extension);
  }

  return extendableFragCount;
}

// findLastUnitig looks for a surrogate at the high end of a contig
int
 findLastUnitig(ContigT *contig, int *unitigID) {
  MultiAlignT *ma;
  int i, numUnitigs, isSurrogate = FALSE;
  float maxContigPos = 0.0;
  NodeCGW_T *unitig = NULL;
  
  fprintf(stderr, "in FindLastUnitig\n");
  
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numUnitigs = GetNumIntUnitigPoss(ma->u_list);
  
  // maxContigPos = contig->bpLength.mean - 1.0;
  
  // can't just jump to last unitig since unitigs are arranged by starting position, not ending
  for (i = 0; i < numUnitigs; i++) {
    IntUnitigPos *upos = GetIntUnitigPos(ma->u_list, i);
    unitig = GetGraphNode(ScaffoldGraph->CIGraph, upos->ident);
    // int isSurrogate = unitig->flags.bits.isSurrogate; // (unitig->info.CI.baseID > 0);
	
#ifdef DEBUG_ECR
      fprintf(stderr, "in contig %d, unitig %d is at %f -> %f maxContigPos: %f, isSurrogate: %d, baseID: %d\n", 
              contig->id, unitig->id,
              unitig->offsetAEnd.mean, unitig->offsetBEnd.mean, maxContigPos,
              isSurrogate, unitig->info.CI.baseID);
#endif

    if (unitig->offsetAEnd.mean >= maxContigPos || unitig->offsetBEnd.mean >= maxContigPos) {
        maxContigPos = max(unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);
        *unitigID = unitig->id;
        isSurrogate = unitig->flags.bits.isSurrogate;
      }
  }
  
  if (isSurrogate) {
    fprintf(stderr, "unitig %d on high end of contig %d is surrogate!\n",
            unitig->id, contig->id);
    return 1;
  } else {
    return 0;
  }
}

// findFirstUnitig looks for a surrogate at the low end of a contig
int findFirstUnitig(ContigT *contig, int *unitigID) {
  MultiAlignT *ma;
  int numUnitigs;
  IntUnitigPos *upos;
  NodeCGW_T *unitig;
  int isSurrogate;
  
  fprintf(stderr, "in FindFirstUnitig\n");
  
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE); 
  numUnitigs = GetNumIntUnitigPoss(ma->u_list);
  
  upos = GetIntUnitigPos(ma->u_list, 0);
  unitig = GetGraphNode(ScaffoldGraph->CIGraph, upos->ident);
  isSurrogate = unitig->flags.bits.isSurrogate; // (unitig->info.CI.baseID > 0);
  
#ifdef DEBUG_ECR
    fprintf(stderr, "in contig %d, unitig %d is at %f -> %f, isSurrogate: %d\n", 
            contig->id, unitig->id,
            unitig->offsetAEnd.mean, unitig->offsetBEnd.mean,
            isSurrogate);
#endif
  
  *unitigID = unitig->id;
  if (isSurrogate) {
    fprintf(stderr, "unitig %d on low end of contig %d is surrogate!\n",
            unitig->id, contig->id);
    return 1;
  }
  else
    return 0;
}











void extendContig(ContigT *contig, int extendAEnd) {
  // have to alter the following fields in a NodeCGW_T: bpLength, offsetAEnd, offsetBEnd
  int contigOrientation = contig->offsetAEnd.mean < contig->offsetBEnd.mean ? A_B : B_A;
  MultiAlignT *cma = NULL;
  float lengthDelta;
  
  cma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);

  lengthDelta = strlen(Getchar(cma->consensus, 0)) - contig->bpLength.mean;

  if (contigOrientation == A_B)
    fprintf(stderr, "in extendContig, contig %8d original pos: %.0f, %.0f A_B\n",
            contig->id, contig->offsetAEnd.mean, contig->offsetBEnd.mean);
  else
    fprintf(stderr, "in extendContig, contig %8d original pos: %.0f, %.0f B_A\n",
            contig->id, contig->offsetBEnd.mean, contig->offsetAEnd.mean);
  
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

  if (contigOrientation == A_B)
    fprintf(stderr, "in extendContig, contig %8d altered pos: %.0f, %.0f A_B\n",
            contig->id, contig->offsetAEnd.mean, contig->offsetBEnd.mean);
  else
    fprintf(stderr, "in extendContig, contig %8d altered pos: %.0f, %.0f B_A\n",
            contig->id, contig->offsetBEnd.mean, contig->offsetAEnd.mean);  
}



// stole most of this from OutputUnitigsFromMultiAligns
/********************************************************************************/
int GetNewUnitigMultiAlign(NodeCGW_T *unitig, fragPositions *fragPoss, int extendedFragIid)
{
  GenericMesg			pmesg;
  IntUnitigMesg			ium_mesg;
  int i;
  int numCIs = GetNumGraphNodes(ScaffoldGraph->CIGraph);
  MultiAlignT *ma = CreateEmptyMultiAlignT();
  UnitigStatus   status;

  //fprintf(stderr, "GetNewUnitigMultiAlign()--\n");

#define  USE_UNGAPPED_CONSENSUS_FOR_UNITIG
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

  //    MultiAlignT *ma = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, unitig->id);

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
      //      fprintf(GlobalData->stderrc,"* Skipping unitig %d --- RESOLVEDREPEAT\n",unitig->id);
      // continue;
    default:
      assert(0);
  }

  ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, unitig->id, TRUE);

#if 1
  fprintf(stderr, "for unitig %d, before reforming, strlen(ma->consensus) = " F_SIZE_T "\n", unitig->id, strlen(Getchar(ma->consensus, 0)));
  fprintf(stderr, "for unitig %d, before reforming, consensus:\n%s\n", unitig->id, Getchar(ma->consensus, 0));
#endif

  {
    int extendedFragLeftward, aligned;
    assert (unitig->type != CONTIG_CGW);

    ium_mesg.iaccession = unitig->id;

    ium_mesg.source = NULL;

    ium_mesg.coverage_stat  = unitig->info.CI.coverageStat;
    ium_mesg.status         = status;
    ium_mesg.a_branch_point = unitig->info.CI.branchPointA;
    ium_mesg.b_branch_point = unitig->info.CI.branchPointB;

#ifdef USE_UNGAPPED_CONSENSUS_FOR_UNITIG
    GetMultiAlignUngappedConsensus(ma, ungappedSequence, ungappedQuality);
    ium_mesg.consensus = Getchar(ungappedSequence,0);
    ium_mesg.quality = Getchar(ungappedQuality,0);
    ium_mesg.length = GetMultiAlignUngappedLength(ma);
#else
    ium_mesg.length = GetMultiAlignLength(ma);
    ium_mesg.consensus = Getchar(ma->consensus,0);
    ium_mesg.quality = Getchar(ma->quality,0);
#endif
    ium_mesg.forced = 0;
    ium_mesg.num_frags = GetNumIntMultiPoss(ma->f_list);


    // replace the positions in the f_list with the adjusted positions
    extendedFragLeftward = FALSE;

    for (i = 0; i < GetNumIntMultiPoss(ma->f_list); i++) {
      IntMultiPos *tempPos = GetIntMultiPos(ma->f_list, i);

#if 1
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
    ium_mesg.f_list = GetIntMultiPos(ma->f_list, 0);

    if (extendedFragLeftward)
      // definitely reorder frags in f_list if we extended a frag leftward
      leftShiftIUM(ium_mesg.f_list, GetNumIntMultiPoss(ma->f_list), extendedFragIid);
    else
      // might need to reorder frags in f_list if we extended a frag rightward
      rightShiftIUM(ium_mesg.f_list, GetNumIntMultiPoss(ma->f_list), extendedFragIid);

    //  The Local_Overlap_AS_forCNS below also had a commented out
    //  DP_Compare.
    //
    {
    CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                            CNS_OPTIONS_SMOOTH_WIN_DEFAULT,
                            CNS_OPTIONS_MAX_NUM_ALLELES };
    ALIGNMENT_CONTEXT=AS_CONSENSUS;


    //  Added options, as in consensus.

    cnslog = stderr;
    aligned = MultiAlignUnitig(&ium_mesg, 
                               ScaffoldGraph->fragStore,
                               reformed_consensus,
                               reformed_quality,
                               reformed_deltas,
                               CNS_STATS_ONLY,  //  CNS_VERBOSE
                               1,
                               Local_Overlap_AS_forCNS,
                               &options);
    }

    if (aligned == -1) {
      fprintf(stderr, "MultiAlignUnitig failure on unitig %d\n", unitig->id);
      return FALSE;   // assert(0);
    }

#ifdef DEBUG_ECR
    fprintf(stderr, "for unitig %d, after reforming, strlen(reformed_consensus) = " F_SIZE_T "\n", unitig->id, strlen(Getchar(reformed_consensus, 0)));
    fprintf(stderr, "for unitig %d, after reforming, consensus:\n%s\n", unitig->id, Getchar(reformed_consensus, 0));
#endif

    {

      // "-2" tells CreateMultiAlignTFromIUM to copy over the source
      // field in an IntMultiPos from the ium_mesg
      //
      MultiAlignT *new_ma = CreateMultiAlignTFromIUM(&ium_mesg, -2, FALSE);

      // if (new_ma == NULL) ...   handle failure, but must always unload multialign since we changed it

      // This adds a reference only if keepInCache is true...
      UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg.iaccession, TRUE);
      InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg.iaccession, TRUE, new_ma, FALSE);

      //PrintMultiAlignT(GlobalData->stderrc, new_ma, ScaffoldGraph->fragStore, 0x0, 0, TRUE, 0, READSTRUCT_LATEST);
    }
  }

  // DeleteMultiAlignT(ma);
  fflush(NULL);
  return TRUE;
}




int
extendCgwClearRange(int fragIid, int frag3pDelta) {
  unsigned int clr_bgn, clr_end;
  int setStatus = 0;

  if (frag3pDelta < 0)
    fprintf(stderr, "Warning: frag3pDelta less than zero: %d\n", frag3pDelta);

  if (fragIid != -1) {
    getFragStore(ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	
    getClearRegion_ReadStruct(fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
    setClearRegion_ReadStruct(fsread, clr_bgn, clr_end + frag3pDelta, READSTRUCT_CGW);
    setStatus = setFragStore(ScaffoldGraph->fragStore, fragIid, fsread);

    fprintf(stderr, "extendCgwClearRange, changed frag %d clr_end from %d to %d\n",
            fragIid, clr_end, clr_end + frag3pDelta);
  }

  return (setStatus); 
}


int
revertToCnsClearRange(int fragIid) {
  unsigned int clr_bgn, clr_end;
  int setStatus = 0;
  
  if (fragIid != -1) {
    getFragStore(ScaffoldGraph->fragStore, fragIid, FRAG_S_ALL, fsread);	
    getClearRegion_ReadStruct(fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
    setClearRegion_ReadStruct(fsread, clr_bgn, clr_end, READSTRUCT_CGW);
    setStatus = setFragStore(ScaffoldGraph->fragStore, fragIid, fsread);
  }
  return (setStatus); 
}





void
getAlteredFragPositions(NodeCGW_T *unitig, fragPositions **fragPoss, int alteredFragIid, int extension) {
  fragPositions *localFragPoss;
  int i, alteredFragIndex = NULLINDEX, orientation = 0, delta;
  MultiAlignT *uma = NULL;

  // currently code does not handle negative extensions, ie, trimming
  if (extension <= 0)
    fprintf(stderr, "negative extension: %d\n", extension);

  // assert (extension > 0);

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

      fprintf(stderr, "getAlteredFragPositions()-- %2d] (%d,%d)\n", i, localFragPoss[i].bgn, localFragPoss[i].end);
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

      fprintf(stderr, "getAlteredFragPositions()-- adjust end fragment %d by %d\n",
              alteredFragIndex, alteredFragDelta);
    } else {
      // all frag positions get bumped by extension, that's how far the altered frag extends into the gap

      // this is the new minimum position in the unitig

      localFragPoss[alteredFragIndex].end = -extension;

      fprintf(stderr, "getAlteredFragPositions()-- adjust first fragment %d to (%d,%d)\n",
              alteredFragIndex, localFragPoss[alteredFragIndex].bgn, localFragPoss[alteredFragIndex].end);

      // if he extends off the front of the unitig, adjust everybody upward

      delta = localFragPoss[alteredFragIndex].end;

      if (delta < 0) {
        for (i = 0; i < GetNumIntMultiPoss(uma->f_list); i++) {
          localFragPoss[i].bgn -= delta;
          localFragPoss[i].end -= delta;
          fprintf(stderr, "getAlteredFragPositions()-- %2d] (%d,%d) (adjusted by -delta = %d)\n", i, localFragPoss[i].bgn, localFragPoss[i].end, -delta);
        }
      }
    }
  *fragPoss = localFragPoss;
}



void
SynchUnitigTWithMultiAlignT(NodeCGW_T *unitig) {
  MultiAlignT *uma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);

  // unitig->bpLength.mean = strlen(Getchar(uma->consensus, 0));      // set length
  unitig->bpLength.mean = GetMultiAlignUngappedLength(uma);

  fprintf(stderr, "in SynchUnitigTWithMultiAlignT, for unitig %d strlen(ma->consensus)=%d  ungapped=%d\n",
          unitig->id,
          strlen(Getchar(uma->consensus, 0)),
          (int)unitig->bpLength.mean);

#if 0
  if (unitig->offsetAEnd.mean < unitig->offsetBEnd.mean) {         // ordering info is okay
    unitig->offsetAEnd.mean = pos->position.bgn;
    unitig->offsetBEnd.mean = pos->position.end;
  } else {
    unitig->offsetAEnd.mean = pos->position.end;
    unitig->offsetBEnd.mean = pos->position.bgn;
  }
#endif  
}



// this routine shifts the frag with iid extendedFragIid to the front of the f_list
void
leftShiftIUM(IntMultiPos *f_list, int numFrags, int extendedFragIid) {
  int i, currPos = 0, numShiftedInUnitig = 0;
  IntMultiPos tempIMP;
  
  fprintf(stderr, "leftShiftIUM()-- \n");

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

  fprintf(stderr, "rightShiftIUM()-- \n");
  
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
        if ((min(f_list[i].position.bgn, f_list[i].position.end) < 
             min(f_list[currPos + numShifted].position.bgn, f_list[currPos + numShifted].position.end)) &&
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

#if 0
      fprintf(stderr, "in restoreFragAndUnitigData, left contig %d:\n", leftFrag->contigID);
      DumpContigMultiAlignInfo (NULL, NULL, leftFrag->contigID);
      DumpContigUngappedOffsets(NULL, NULL, leftFrag->contigID);
#endif
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

#if 0
      fprintf(stderr, "in restoreFragAndUnitigData, right contig %d:\n", rightFrag->contigID);
      DumpContigMultiAlignInfo (NULL, NULL, rightFrag->contigID);
      DumpContigUngappedOffsets(NULL, NULL, rightFrag->contigID);
#endif
    }
}







int
writeEcrCheckpoint(int *numGapsInScaffold, int *numGapsClosedInScaffold,
                       int *numSmallGapsInScaffold, int *numSmallGapsClosedInScaffold,
                       int *numLargeGapsInScaffold, int *numLargeGapsClosedInScaffold) {
  char ckpFileName[1024];
  FILE *ckpFile;
  int i;
  
  sprintf(ckpFileName, "%s.ecr.ckp.%d", GlobalData->File_Name_Prefix, ScaffoldGraph->checkPointIteration - 1);

  ckpFile = fopen(ckpFileName, "w+");
  if (ckpFile == NULL) {
    fprintf(stderr, "Could not open checkpoint file %s, aborting.\n", ckpFileName);
    assert(0);
  }
  
  fprintf(ckpFile, "%d\n", GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));
  for (i = 0; i < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); i++) {
    CIScaffoldT *scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, i);
    assert(scaff != NULL);
	
    //if ((isDeadCIScaffoldT(scaff)) || (scaff->type != REAL_SCAFFOLD) ||	(scaff->info.Scaffold.numElements < 2))
    //continue;
	
    fprintf(ckpFile, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", scaff->id, 
            numGapsInScaffold[ scaff->id ], numGapsClosedInScaffold[ scaff->id ],
            numSmallGapsInScaffold[ scaff->id ], numSmallGapsClosedInScaffold[ scaff->id ],
            numLargeGapsInScaffold[ scaff->id ], numLargeGapsClosedInScaffold[ scaff->id ]);
  }
  fclose(ckpFile);
  return 0;
}

int
loadEcrCheckpoint(int ckptNum, int *numGapsInScaffold, int *numGapsClosedInScaffold,
                      int *numSmallGapsInScaffold, int *numSmallGapsClosedInScaffold,
                      int *numLargeGapsInScaffold, int *numLargeGapsClosedInScaffold) {
  char ckpFileName[1024];
  FILE *ckpFile;
  int i, scaffid;
  int totalGaps = 0, totalGapsClosed = 0, 
    totalSmallGaps = 0, totalSmallGapsClosed = 0, totalLargeGaps = 0, totalLargeGapsClosed = 0;
  int numGraphNodes;
  
  sprintf(ckpFileName, "%s.ecr.ckp.%d", GlobalData->File_Name_Prefix, ckptNum);

  ckpFile = fopen(ckpFileName, "r");
  if (ckpFile == NULL) {
    fprintf(stderr, "Warning: could not open checkpoint file %s for reading, continuing.\n", ckpFileName);
    return -1;
  }
  
  fscanf(ckpFile, "%d\n", &numGraphNodes);
  for (i = 0; i < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); i++) {
    // CIScaffoldT* scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, i);
    // assert(scaff);
    //if ((isDeadCIScaffoldT(scaff)) || (scaff->type != REAL_SCAFFOLD) ||	(scaff->info.Scaffold.numElements < 2))
    //continue;
	
    fscanf(ckpFile, "%d %d %d %d %d %d %d\n", &scaffid, 
           &numGapsInScaffold[ i ], &numGapsClosedInScaffold[ i ],
           &numSmallGapsInScaffold[ i ], &numSmallGapsClosedInScaffold[ i ],
           &numLargeGapsInScaffold[ i ], &numLargeGapsClosedInScaffold[ i ]);

    //if (scaffid != i)
    // fprintf(stderr, "Warning: scaffid (%d) != i (%d) in loadEcrCheckpoint!\n",
    //		   scaffid, i);

    totalGaps += numGapsInScaffold[ scaffid ];
    totalGapsClosed += numGapsClosedInScaffold[ scaffid ];
    totalSmallGaps += numSmallGapsInScaffold[ scaffid ];
    totalSmallGapsClosed += numSmallGapsClosedInScaffold[ scaffid ];
    totalLargeGaps += numLargeGapsInScaffold[ scaffid ];
    totalLargeGapsClosed += numLargeGapsClosedInScaffold[ scaffid ];
  }
  fclose(ckpFile);

  fprintf(stderr, "Stats after loadEcrCheckpoint for ckptNum %d:\n", ckptNum);
  fprintf(stderr, "           totalGaps: %d\n", totalGaps);
  fprintf(stderr, "     totalGapsClosed: %d\n", totalGapsClosed);
  fprintf(stderr, "      totalSmallGaps: %d\n", totalSmallGaps);
  fprintf(stderr, "totalSmallGapsClosed: %d\n", totalSmallGapsClosed);
  fprintf(stderr, "      totalLargeGaps: %d\n", totalLargeGaps);
  fprintf(stderr, "totalLargeGapsClosed: %d\n", totalLargeGapsClosed);
  
  return 0;
}
