
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

static char *rcsid = "$Id: Input_CGW.c,v 1.66 2009-09-14 16:09:04 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_CGW_dataTypes.h"
#include "AS_PER_gkpStore.h"
#include "ScaffoldGraph_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "Input_CGW.h"


static int32 DiscriminatorUniques = 0;
static int32 ShortDiscriminatorUniques = 0;

static int32 totalReadFrags = 0;
static int32 inUniqueReadFrags = 0;
static int32 inRepeatReadFrags = 0;
static int32 inTeenyUnitigReadFrags = 0;
static int32 inSingletonUnitigReadFrags = 0;
static int32 onEndReadFrags = 0;

static int32 TouchesContained = 0;
static int32 TransChunk = 0;
static int32 Containment = 0;
static int32 DoveTail = 0;
static int32 Tandem = 0;
static int32 BetweenContained = 0;
static int32 ContainStack = 0;
static int32 BadQuality = 0;


int ProcessInput(int optind, int argc, char *argv[]){
  GenericMesg   *pmesg;
  FILE *infp;
  int i,j = 0;
  int32 numIUM = 0;

  EnableRange_VA(ScaffoldGraph->CIFrags, ScaffoldGraph->gkpStore->gkStore_getNumFragments() + 1);

  gkFragment  fr;
  gkStream   *fs = new gkStream(ScaffoldGraph->gkpStore, 0, 0, GKFRAGMENT_INF);

  //  There is no zero fragment.  Make the data junk (this also sets isDeleted, so if we do ever
  //  happen to use the fragment, it'll generally be skipped.

  memset(GetCIFragT(ScaffoldGraph->CIFrags, 0), 0xff, sizeof(CIFragT));

  while (fs->next(&fr)) {
    CIFragT   *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, fr.gkFragment_getReadIID());

    cifrag->read_iid                              = fr.gkFragment_getReadIID();
    cifrag->mate_iid                              = fr.gkFragment_getMateIID();
    cifrag->dist                                  = fr.gkFragment_getLibraryIID();

    cifrag->cid                                   = NULLINDEX;
    cifrag->CIid                                  = NULLINDEX;
    cifrag->offset5p.mean                         = 0.0;
    cifrag->offset5p.variance                     = 0.0;
    cifrag->offset3p.mean                         = 0.0;
    cifrag->offset3p.variance                     = 0.0;

    cifrag->contigID                              = NULLINDEX;
    cifrag->contigOffset5p.mean                   = 0.0;
    cifrag->contigOffset5p.variance               = 0.0;
    cifrag->contigOffset3p.mean                   = 0.0;
    cifrag->contigOffset3p.variance               = 0.0;

    cifrag->flags.bits.hasInternalOnlyCILinks     = FALSE;
    cifrag->flags.bits.hasInternalOnlyContigLinks = FALSE;
    cifrag->flags.bits.isPlaced                   = FALSE;
    cifrag->flags.bits.isSingleton                = FALSE;
    cifrag->flags.bits.isChaff                    = FALSE;
    cifrag->flags.bits.isDeleted                  = fr.gkFragment_getIsDeleted();
    cifrag->flags.bits.innieMate                  = (fr.gkFragment_getOrientation() == AS_READ_ORIENT_INNIE);
    cifrag->flags.bits.hasMate                    = (fr.gkFragment_getMateIID()      > 0);

    cifrag->flags.bits.edgeStatus                 = (fr.gkFragment_getMateIID()      > 0) ? UNKNOWN_EDGE_STATUS : INVALID_EDGE_STATUS;
    cifrag->flags.bits.chunkLabel                 = AS_SINGLETON;

    cifrag->flags.bits.mateDetail                 = UNASSIGNED_MATE;
  }


  int              extraUnitigsMax = 0;
  int              extraUnitigsLen = 0;
  IntUnitigMesg   *extraUnitigs    = NULL;

  for(i = optind; i < argc; i++){
    infp = fopen(argv[i],"r");

    while ((EOF != ReadProtoMesg_AS(infp, &pmesg))) {
      if (pmesg->t == MESG_IUM) {
        IntUnitigMesg *ium_mesg = (IntUnitigMesg *)pmesg->m;

        //  Insert both a unitig and a contig -- these MUST be saved
        //  in the cache, else a huge memory leak.
        //
        //  If one already exists for that iacc, save it for later and
        //  give it a new iacc.  fixUnitigs generates these guys.

        if (ium_mesg->iaccession < 1000000000) {
          MultiAlignT *uma = CreateMultiAlignTFromIUM(ium_mesg, 0, FALSE);
          MultiAlignT *cma = CopyMultiAlignT(NULL, uma);

          insertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, TRUE,  uma, TRUE);
          insertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, ium_mesg->iaccession, FALSE, cma, TRUE);

          ProcessIUM_ScaffoldGraph(ium_mesg, GetMultiAlignUngappedLength(uma), FALSE);

          numIUM++;
          if ((numIUM % 10000) == 0) {
            fprintf(stderr, "processed "F_S32" IUM messages.\n", numIUM);
            clearCacheSequenceDB(ScaffoldGraph->sequenceDB);
          }
        } else {
          int l = strlen(ium_mesg->consensus) + 1;

          if (extraUnitigsLen >= extraUnitigsMax) {
            extraUnitigsMax = (extraUnitigsMax == 0) ? 1024 : (extraUnitigsMax * 2);
            extraUnitigs    = (IntUnitigMesg *)safe_realloc(extraUnitigs, sizeof(IntUnitigMesg) * extraUnitigsMax);
          }

          extraUnitigs[extraUnitigsLen] = *ium_mesg;

          extraUnitigs[extraUnitigsLen].consensus = (char *)safe_malloc(sizeof(char) * l);
          extraUnitigs[extraUnitigsLen].quality   = (char *)safe_malloc(sizeof(char) * l);
          extraUnitigs[extraUnitigsLen].f_list    = (IntMultiPos *)safe_malloc(sizeof(IntMultiPos) * ium_mesg->num_frags);

          memcpy(extraUnitigs[extraUnitigsLen].consensus, ium_mesg->consensus, sizeof(char) * l);
          memcpy(extraUnitigs[extraUnitigsLen].quality,   ium_mesg->quality,   sizeof(char) * l);
          memcpy(extraUnitigs[extraUnitigsLen].f_list,    ium_mesg->f_list,    sizeof(IntMultiPos) * ium_mesg->num_frags);

          extraUnitigsLen++;
        }
      }
    }
    fclose(infp);
  }


  for (i=0; i<extraUnitigsLen; i++) {
    fprintf(stderr, "WARNING: assign new accession to duplicate IUM "F_U32"; now "F_U32"\n", extraUnitigs[i].iaccession, numIUM);

    extraUnitigs[i].iaccession = numIUM++;

    MultiAlignT *uma = CreateMultiAlignTFromIUM(extraUnitigs + i, 0, FALSE);
    MultiAlignT *cma = CopyMultiAlignT(NULL, uma);

    insertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, extraUnitigs[i].iaccession, TRUE,  uma, TRUE);
    insertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, extraUnitigs[i].iaccession, FALSE, cma, TRUE);

    ProcessIUM_ScaffoldGraph(extraUnitigs + i, GetMultiAlignUngappedLength(uma), FALSE);

    safe_free(extraUnitigs[i].consensus);
    safe_free(extraUnitigs[i].quality);
    safe_free(extraUnitigs[i].f_list);
  }

  safe_free(extraUnitigs);

  fprintf(stderr,"Processed %d IUM messages (max IUM acc = %d) with %d fragments\n",
          numIUM,
          (int)GetNumGraphNodes(ScaffoldGraph->CIGraph),
	  (int)GetNumCIFragTs(ScaffoldGraph->CIFrags));

  ScaffoldGraph->numLiveCIs     = GetNumGraphNodes(ScaffoldGraph->CIGraph);
  ScaffoldGraph->numOriginalCIs = GetNumGraphNodes(ScaffoldGraph->CIGraph);


  fprintf(stderr,"* Total Long Discriminator Uniques : %d   Short Uniques: %d\n",
          DiscriminatorUniques,
          ShortDiscriminatorUniques);
  fprintf(stderr,"* Total Reads:%d in discriminator unique:%d in other:%d ; in teeny: %d in singles:%d on ends:%d\n",
	  totalReadFrags, inUniqueReadFrags, inRepeatReadFrags, inTeenyUnitigReadFrags, inSingletonUnitigReadFrags,
	  onEndReadFrags);

  return(0);
}





void ProcessIUM_ScaffoldGraph(IntUnitigMesg *ium_mesg, int32 length, int sequenceOnly){
  CDS_CID_t cfr;
  ChunkInstanceT CI;

  memset(&CI, 0, sizeof(ChunkInstanceT));

  CI.id = ium_mesg->iaccession;
  CI.bpLength.mean = length;
  CI.bpLength.variance = MAX(1.0,ComputeFudgeVariance(CI.bpLength.mean));
  CI.edgeHead = NULLINDEX;
  CI.setID = NULLINDEX;
  CI.scaffoldID = NULLINDEX;
  CI.indexInScaffold = NULLINDEX;
  CI.prevScaffoldID = NULLINDEX;
  CI.numEssentialA = 0;
  CI.numEssentialB = 0;
  CI.essentialEdgeA = NULLINDEX;
  CI.essentialEdgeB = NULLINDEX;
  CI.smoothExpectedCID = NULLINDEX;
  CI.BEndNext = CI.AEndNext = NULLINDEX;
  CI.info.CI.numFragments = ium_mesg->num_frags;
  CI.info.CI.coverageStat = (ium_mesg->coverage_stat < -1000.0? -1000:ium_mesg->coverage_stat);
  CI.info.CI.microhetProb = ium_mesg->microhet_prob;
  CI.info.CI.forceUniqueRepeat = ium_mesg->unique_rept;
  CI.info.CI.contigID = NULLINDEX;
  CI.info.CI.numInstances = 0;
  CI.info.CI.instances.in_line.instance1 = 0;
  CI.info.CI.instances.in_line.instance2 = 0;
  CI.info.CI.instances.va = NULL;
  CI.info.CI.source = NULLINDEX;
  CI.flags.all = 0;
  CI.offsetAEnd.mean = 0.0;
  CI.offsetAEnd.variance = 0.0;
  CI.offsetBEnd = CI.bpLength;

  if(ium_mesg->coverage_stat >= GlobalData->cgbUniqueCutoff){
    if(length < CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH ||
       ium_mesg->num_frags < CGW_MIN_READS_IN_UNIQUE){
      ShortDiscriminatorUniques++;
    }else{
      DiscriminatorUniques++;
    }
  }


  {
    int isUnique = FALSE;
    if(ium_mesg->coverage_stat >= GlobalData->cgbUniqueCutoff &&
       length >= CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH &&
       ium_mesg->num_frags >= CGW_MIN_READS_IN_UNIQUE){
      // microhet probability is actually the probability of the
      // sequence being UNIQUE, based on microhet considerations.
      // Falling below threshhold makes something a repeat.

      if( CI.info.CI.microhetProb < GlobalData->cgbMicrohetProb){
	if(ium_mesg->coverage_stat < GlobalData->cgbApplyMicrohetCutoff){
	  //fprintf(stderr,"* CI " F_CID " with astat: %g classified as repeat based on microhet unique prob of %g < %g\n",
          //        CI.id, ium_mesg->coverage_stat, CI.info.CI.microhetProb, GlobalData->cgbMicrohetProb);
	  isUnique = FALSE;
	  CI.type = UNRESOLVEDCHUNK_CGW;
	}else{
	  isUnique = TRUE;
	  //fprintf(stderr,"* WARNING: CI " F_CID " with coverage %g WOULD HAVE BEEN classified as repeat based on microhet unique prob of %g < %g\n",
          //        CI.id, ium_mesg->coverage_stat, CI.info.CI.microhetProb, GlobalData->cgbMicrohetProb);
	}
      }else{
	isUnique = TRUE;
      }
    }else{
      isUnique = FALSE;
    }

    // allow flag to overwrite what the default behavior for a chunk and force it to be unique or repeat

    if (CI.info.CI.forceUniqueRepeat == AS_FORCED_UNIQUE) {
       isUnique = TRUE;
    }
    else if (CI.info.CI.forceUniqueRepeat == AS_FORCED_REPEAT) {
       isUnique = FALSE;
    }

    if(isUnique){
      ScaffoldGraph->numDiscriminatorUniqueCIs++;
      CI.flags.bits.isUnique = 1;
      CI.type = DISCRIMINATORUNIQUECHUNK_CGW;
    }else{
      CI.flags.bits.isUnique = 0;
      CI.type = UNRESOLVEDCHUNK_CGW;
    }
  }

  CI.flags.bits.smoothSeenAlready = FALSE;
  CI.flags.bits.isCI              = TRUE;
  CI.flags.bits.isChaff           = FALSE;
  CI.flags.bits.isClosure         = FALSE;

  if( ! sequenceOnly ) {
      CDS_CID_t extremalA = NULLINDEX;
      CDS_CID_t extremalB = NULLINDEX;
      int32 minOffset = INT32_MAX;
      int32 maxOffset = INT32_MIN;

      for(cfr = 0; cfr < ium_mesg->num_frags; cfr++){
	IntMultiPos *cfr_mesg = ium_mesg->f_list + cfr;
	int32 end = MAX( cfr_mesg->position.end, cfr_mesg->position.bgn);
	int32 beg = MIN( cfr_mesg->position.end, cfr_mesg->position.bgn);

	if(minOffset > beg){
	  minOffset = beg;
	  extremalA = cfr;
	}
	if(maxOffset < end){
	  maxOffset = end;
	  extremalB = cfr;
	}
      }


      for(cfr = 0; cfr < ium_mesg->num_frags; cfr++){
	IntMultiPos *cfr_mesg = ium_mesg->f_list + cfr;
        CIFragT     *cifrag   = GetCIFragT(ScaffoldGraph->CIFrags, cfr_mesg->ident);

	cifrag->cid      = ium_mesg->iaccession;
	cifrag->CIid     = ium_mesg->iaccession;

        //  Singleton chunks are chaff; singleton frags are chaff unless proven otherwise
        //
        if (ium_mesg->num_frags < 2) {
          CI.flags.bits.isChaff          = TRUE;
          cifrag->flags.bits.isSingleton = TRUE;
          cifrag->flags.bits.isChaff     = TRUE;
	}

	// Collect read stats
        if (AS_FA_READ(cfr_mesg->type)) {
	  totalReadFrags++;

	  if(CI.flags.bits.isUnique)
	    inUniqueReadFrags++;
	  else
	    inRepeatReadFrags++;

	  if(ium_mesg->num_frags <= 2)
            inTeenyUnitigReadFrags++;

          if(ium_mesg->num_frags < 2)
	    inSingletonUnitigReadFrags++;

	  if(cfr == extremalA || cfr == extremalB)
	    onEndReadFrags++;
	}
      }
    }

  // Insert the Chunk Instance
  SetChunkInstanceT(ScaffoldGraph->CIGraph->nodes, CI.id, &CI);

  // Mark all frags as being members of this CI, and set their offsets within
  // the CI
  if( ! sequenceOnly )
    UpdateNodeFragments(ScaffoldGraph->CIGraph,CI.id, CI.type == DISCRIMINATORUNIQUECHUNK_CGW, TRUE ); // mark unitigs and contigs
}




void
LoadDistData(void) {
  int32 numDists = ScaffoldGraph->gkpStore->gkStore_getNumLibraries();
  CDS_CID_t i;

  for(i = 1; i <= numDists; i++){
    DistT dist;
    gkLibrary  *gkpl = ScaffoldGraph->gkpStore->gkStore_getLibrary(i);

    dist.mu             = gkpl->mean;
    dist.sigma          = gkpl->stddev;
    dist.numSamples     = 0;
    dist.min            = INT32_MAX;
    dist.max            = INT32_MIN;
    dist.bnum           = 0;
    dist.bsize          = 0;
    dist.histogram      = NULL;
    dist.lower          = dist.mu - CGW_CUTOFF * dist.sigma;
    dist.upper          = dist.mu + CGW_CUTOFF * dist.sigma;
    dist.numReferences  = 0;
    dist.numBad         = 0;

    fprintf(stderr,"* Loaded dist %s,"F_CID" (%g +/- %g)\n",
            AS_UID_toString(gkpl->libraryUID), i, dist.mu, dist.sigma);

    SetDistT(ScaffoldGraph->Dists, i, &dist);
  }
}
