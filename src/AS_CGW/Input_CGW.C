
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

static char *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.H"
#include "AS_UTL_Var.H"
#include "AS_CGW_dataTypes.H"
#include "AS_PER_gkpStore.H"
#include "ScaffoldGraph_CGW.H"
#include "Globals_CGW.H"
#include "ScaffoldGraph_CGW.H"
#include "Output_CGW.H"
#include "Input_CGW.H"


#define CGW_MIN_READS_IN_UNIQUE  2

// Due to the FBAC fragments, we get some pathologically short U-Unitigs
// Set the following threshhold to eliminate the really short ones
//
#define CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH 1000


//  For each unitig, write logging on the repeat/unique decision.
#undef LOG_REPEAT_UNIQUE_LABELING


//  Stats on repeat labeling of input unitigs.
//
class ruLabelStat {
public:
  ruLabelStat() {
    num = 0;
    len = 0;
  };

  void     operator+=(uint64 len_) {
    num++;
    len += len_;
  };

  uint32   num;
  uint64   len;
};

ruLabelStat  repeat_LowReads;
ruLabelStat  repeat_LowCovStat;
ruLabelStat  repeat_Short;
ruLabelStat  repeat_MicroHet;

ruLabelStat  repeat_IsUnique;
ruLabelStat  repeat_IsRepeat;



static
void
ProcessInputUnitig(MultiAlignT *uma) {
  CDS_CID_t       cfr;
  ChunkInstanceT  CI;

  int32           length = GetMultiAlignUngappedLength(uma);

  memset(&CI, 0, sizeof(ChunkInstanceT));

  CI.id                                  = uma->maID;
  CI.bpLength.mean                       = length;
  CI.bpLength.variance                   = MAX(1.0,ComputeFudgeVariance(CI.bpLength.mean));
  CI.setID                               = NULLINDEX;
  CI.scaffoldID                          = NULLINDEX;
  CI.indexInScaffold                     = NULLINDEX;
  CI.prevScaffoldID                      = NULLINDEX;
  CI.numEssentialA                       = 0;
  CI.numEssentialB                       = 0;
  CI.essentialEdgeA                      = NULLINDEX;
  CI.essentialEdgeB                      = NULLINDEX;
  CI.smoothExpectedCID                   = NULLINDEX;
  CI.BEndNext                            = NULLINDEX;
  CI.AEndNext                            = NULLINDEX;

  if (ScaffoldGraph->tigStore->getUnitigCoverageStat(uma->maID) < -1000)
    ScaffoldGraph->tigStore->setUnitigCoverageStat(uma->maID, -1000);

  CI.info.CI.contigID                    = NULLINDEX;
  CI.info.CI.numInstances                = 0;
  CI.info.CI.instances.in_line.instance1 = 0;
  CI.info.CI.instances.in_line.instance2 = 0;
  CI.info.CI.instances.va                = NULL;
  CI.info.CI.source                      = NULLINDEX;

  CI.flags.all                           = 0;

  CI.offsetAEnd.mean                     = 0.0;
  CI.offsetAEnd.variance                 = 0.0;
  CI.offsetBEnd                          = CI.bpLength;

  int  isUnique     = TRUE;


  //fprintf(stderr, "unitig %d "F_SIZE_T" reads %d length\n",
  //        uma->maID, GetNumIntMultiPoss(uma->f_list), length);


  if (GetNumIntMultiPoss(uma->f_list) < CGW_MIN_READS_IN_UNIQUE) {
#ifdef LOG_REPEAT_UNIQUE_LABELING
    fprintf(stderr, "unitig %d not unique -- "F_SIZE_T" reads, need at least %d\n",
             uma->maID, GetNumIntMultiPoss(uma->f_list), CGW_MIN_READS_IN_UNIQUE);
#endif
    repeat_LowReads += length;
    isUnique = FALSE;
  }

  if (ScaffoldGraph->tigStore->getUnitigCoverageStat(uma->maID) < GlobalData->cgbUniqueCutoff) {
#ifdef LOG_REPEAT_UNIQUE_LABELING
    fprintf(stderr, "unitig %d not unique -- coverage stat %d, needs to be at least %d\n",
            uma->maID, ScaffoldGraph->tigStore->getUnitigCoverageStat(uma->maID), GlobalData->cgbUniqueCutoff);
#endif
    repeat_LowCovStat += length;
    isUnique = FALSE;
  }

#ifdef SHORT_HIGH_ASTAT_ARE_UNIQUE
  //  This is an attempt to not blindly call all short unitigs as non-unique.  It didn't work so
  //  well in initial limited testing.  The threshold is arbitrary; older versions used
  //  cgbDefinitelyUniqueCutoff.
   if ((ScaffoldGraph->tigStore->getUnitigCoverageStat(uma->maID) < GlobalData->cgbUniqueCutoff * 10) &&
      (length < CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH)) {
#ifdef LOG_REPEAT_UNIQUE_LABELING
    fprintf(stderr, "unitig %d not unique -- length %d too short, need to be at least %d AND coverage stat %d must be larger than %d\n",
            uma->maID, length, CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH,
            ScaffoldGraph->tigStore->getUnitigCoverageStat(uma->maID), GlobalData->cgbUniqueCutoff * 10);
#endif
    repeat_Short += length;
    isUnique = FALSE;
  }
#else
  if (length < CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH) {
#ifdef LOG_REPEAT_UNIQUE_LABELING
    fprintf(stderr, "unitig %d not unique -- length %d too short, need to be at least %d\n",
            uma->maID, length, CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH);
#endif
    repeat_Short += length;
    isUnique = FALSE;
  }
#endif

  //  MicroHet probability is actually the probability of the sequence being UNIQUE, based on
  //  microhet considerations.  Falling below threshhold makes something a repeat.
  //  Note that this is off by default (see options -e, -i)
  //  Defaults  cgbMicrohetProb        = 1.0 e-5
  //            cgbApplyMicrohetCutoff = -1
  if ((ScaffoldGraph->tigStore->getUnitigMicroHetProb(uma->maID) < GlobalData->cgbMicrohetProb) &&
      (ScaffoldGraph->tigStore->getUnitigCoverageStat(uma->maID) < GlobalData->cgbApplyMicrohetCutoff)) {
#ifdef LOG_REPEAT_UNIQUE_LABELING
    fprintf(stderr, "unitig %d not unique -- low microhetprob %f (< %f) and low coverage stat %d (< %d)\n",
            uma->maID,
            ScaffoldGraph->tigStore->getUnitigMicroHetProb(uma->maID), GlobalData->cgbMicrohetProb,
            ScaffoldGraph->tigStore->getUnitigCoverageStat(uma->maID), GlobalData->cgbApplyMicrohetCutoff);
#endif
    repeat_MicroHet += length;
    isUnique = FALSE;
  }

  // allow flag to overwrite what the default behavior for a chunk and force it to be unique or repeat

  if (ScaffoldGraph->tigStore->getUnitigFUR(CI.id) == AS_FORCED_UNIQUE) {
#ifdef LOG_REPEAT_UNIQUE_LABELING
    fprintf(stderr, "unitig %d forced unique\n",
            uma->maID);
#endif
    isUnique = TRUE;
  }

  if (ScaffoldGraph->tigStore->getUnitigFUR(CI.id) == AS_FORCED_REPEAT) {
#ifdef LOG_REPEAT_UNIQUE_LABELING
    fprintf(stderr, "unitig %d forced repeat\n",
            uma->maID);
#endif
    isUnique = FALSE;
  }

  if (isUnique) {
    CI.flags.bits.isUnique = 1;
    CI.type                = DISCRIMINATORUNIQUECHUNK_CGW;

    repeat_IsUnique += length;
  } else {
    CI.flags.bits.isUnique = 0;
    CI.type                = UNRESOLVEDCHUNK_CGW;

    repeat_IsRepeat += length;
  }

  CI.flags.bits.smoothSeenAlready = FALSE;
  CI.flags.bits.isCI              = TRUE;
  CI.flags.bits.isChaff           = FALSE;
  CI.flags.bits.isClosure         = FALSE;

  //  Assign each read to this unitig.

  for(cfr = 0; cfr < GetNumIntMultiPoss(uma->f_list); cfr++){
    IntMultiPos *imp    = GetIntMultiPos(uma->f_list, cfr);
    CIFragT     *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, imp->ident);

    cifrag->cid  = uma->maID;
    cifrag->CIid = uma->maID;
  }

  //  Decide if this unitig is a potential rock or stone (note, uses cid mark just set)

  uint32  numExternMates = 0;

  for(cfr = 0; cfr < GetNumIntMultiPoss(uma->f_list); cfr++){
    IntMultiPos *imp    = GetIntMultiPos(uma->f_list, cfr);
    CIFragT     *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, imp->ident);

    if (cifrag->flags.bits.hasMate == false)
      continue;

    CIFragT     *mifrag = GetCIFragT(ScaffoldGraph->CIFrags, cifrag->mate_iid);

    if (cifrag->cid != mifrag->cid)
      numExternMates++;
  }

  if        (numExternMates == 0) {
    CI.flags.bits.isPotentialRock  = FALSE;
    CI.flags.bits.isPotentialStone = TRUE;

  } else if (numExternMates == 1) {
    CI.flags.bits.isPotentialRock  = FALSE;
    CI.flags.bits.isPotentialStone = TRUE;

  } else {
    CI.flags.bits.isPotentialRock  = TRUE;
    CI.flags.bits.isPotentialStone = TRUE;
  }

  //  Singleton chunks are chaff; singleton frags are chaff unless proven otherwise

  if (GetNumIntMultiPoss(uma->f_list) < 2) {
    IntMultiPos *imp    = GetIntMultiPos(uma->f_list, 0);
    CIFragT     *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, imp->ident);

    CI.flags.bits.isChaff          = TRUE;
    cifrag->flags.bits.isSingleton = TRUE;
    cifrag->flags.bits.isChaff     = TRUE;
  }

  // Insert the Chunk Instance

  SetGraphNode(ScaffoldGraph->CIGraph, &CI);

  //  Ensure that there are no edges, and that the edgeList is allocated.
  assert(ScaffoldGraph->CIGraph->edgeLists[CI.id].empty() == true);

  // Mark all frags as being members of this CI, and set their offsets within the CI
  // mark unitigs and contigs
  UpdateNodeFragments(ScaffoldGraph->CIGraph,CI.id, CI.type == DISCRIMINATORUNIQUECHUNK_CGW, TRUE);
}


void
ReloadMatesFromGatekeeper(void) {

  fprintf(stderr, "Reading fragments - just for mate pairs.\n");

  uint32      nAlreadySet = 0;
  uint32      nNotMated   = 0;
  uint32      nAdded      = 0;

  gkFragment  fr;
  gkStream   *fs = new gkStream(ScaffoldGraph->gkpStore, 0, 0, GKFRAGMENT_INF);

  while (fs->next(&fr)) {
    CIFragT   *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, fr.gkFragment_getReadIID());

    if (cifrag->mate_iid > 0) {
      //  Mate already set, do not modify the CIFragT.
      nAlreadySet++;
      continue;
    }

    if (fr.gkFragment_getMateIID() == 0) {
      //  Not mated, do not modify the CIFragT.
      nNotMated++;
      continue;
    }

    nAdded++;

    assert(cifrag->read_iid == fr.gkFragment_getReadIID());

    cifrag->mate_iid                              = fr.gkFragment_getMateIID();
    cifrag->dist                                  = fr.gkFragment_getLibraryIID();

    cifrag->flags.bits.hasInternalOnlyCILinks     = FALSE;
    cifrag->flags.bits.hasInternalOnlyContigLinks = FALSE;

    cifrag->flags.bits.innieMate                  = (fr.gkFragment_getOrientation() == AS_READ_ORIENT_INNIE);
    cifrag->flags.bits.hasMate                    = 1;

    cifrag->flags.bits.edgeStatus                 = UNKNOWN_EDGE_STATUS;

    //  Set in ComputeMatePairDetailedStatus() appears unused though
    cifrag->flags.bits.mateDetail                 = UNASSIGNED_MATE;
  }

  delete fs;

  fprintf(stderr, "Reading fragments - just for mate pairs - "F_U32" reads already mated; "F_U32" reads unmated; "F_U32" reads with new mate added.\n",
          nAlreadySet, nNotMated, nAdded);

  //  Decide if this unitig is a potential rock or stone (note, uses cid mark just set)

  for (int32 i=0; i<ScaffoldGraph->tigStore->numUnitigs(); i++) {
    MultiAlignT   *uma = ScaffoldGraph->tigStore->loadMultiAlign(i, TRUE);

    if (uma == NULL)
      continue;

    uint32  numExternMates = 0;

    for(int32 cfr=0; cfr<GetNumIntMultiPoss(uma->f_list); cfr++){
      IntMultiPos *imp    = GetIntMultiPos(uma->f_list, cfr);
      CIFragT     *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, imp->ident);

      if (cifrag->flags.bits.hasMate == false)
        continue;

      CIFragT     *mifrag = GetCIFragT(ScaffoldGraph->CIFrags, cifrag->mate_iid);

      if (cifrag->cid != mifrag->cid)
        numExternMates++;
    }

    NodeCGW_T *CI = GetGraphNode(ScaffoldGraph->CIGraph, uma->maID);

    if        (numExternMates == 0) {
      CI->flags.bits.isPotentialRock  = FALSE;
      CI->flags.bits.isPotentialStone = TRUE;
      
    } else if (numExternMates == 1) {
      CI->flags.bits.isPotentialRock  = FALSE;
      CI->flags.bits.isPotentialStone = TRUE;

    } else {
      CI->flags.bits.isPotentialRock  = TRUE;
      CI->flags.bits.isPotentialStone = TRUE;
    }
  }
}


void
ProcessInput(int optind, int argc, char *argv[]){
  int32  numFRG = 0;
  int32  numUTG = 0;

  fprintf(stderr, "Reading fragments.\n");

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
    cifrag->contigID                              = NULLINDEX;

    cifrag->offset5p.mean                         = 0.0;
    cifrag->offset5p.variance                     = 0.0;
    cifrag->offset3p.mean                         = 0.0;
    cifrag->offset3p.variance                     = 0.0;

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
    cifrag->flags.bits.hasMate                    = (fr.gkFragment_getMateIID() > 0);

    cifrag->flags.bits.edgeStatus                 = (fr.gkFragment_getMateIID() > 0) ? UNKNOWN_EDGE_STATUS : INVALID_EDGE_STATUS;
    cifrag->flags.bits.chunkLabel                 = AS_SINGLETON;

    cifrag->flags.bits.mateDetail                 = UNASSIGNED_MATE;

    if ((++numFRG % 1000000) == 0) {
      fprintf(stderr, "...processed "F_S32" fragments.\n", numFRG);
    }
  }

  delete fs;

  fprintf(stderr, "Reading unitigs.\n");

  //  We flush the cache while loading, as each unitig is immediately copied to a contig.  Even
  //  thought the immediate next step (of checking for chimeric unitigs) will reload all unitigs
  //  again, the flush is useful (it gets rid of that second copy in a unitig).  During assembly, we
  //  usually never need ALL this stuff loaded at the same time.

  for (int32 i=0; i<ScaffoldGraph->tigStore->numUnitigs(); i++) {
    MultiAlignT   *uma = ScaffoldGraph->tigStore->loadMultiAlign(i, TRUE);

    if (uma == NULL) {
      ChunkInstanceT  CI;

      //fprintf(stderr, "WARNING:  Unitig %d does not exist.\n", i);

      memset(&CI, 0, sizeof(ChunkInstanceT));

      //  BEWARE!  Magic.  We try to fake the process of ProcessInputUnitig() then DeleteGraphNode().

      CI.id                                  = i;                    //  uma->maID, if only it existed

      CI.bpLength.mean                       = 0;          //  Boilerplate from below, to NULLINDEX
      CI.bpLength.variance                   = 1.0;
      CI.setID                               = NULLINDEX;
      CI.scaffoldID                          = NULLINDEX;
      CI.indexInScaffold                     = NULLINDEX;
      CI.prevScaffoldID                      = NULLINDEX;
      CI.essentialEdgeA                      = NULLINDEX;
      CI.essentialEdgeB                      = NULLINDEX;
      CI.smoothExpectedCID                   = NULLINDEX;
      CI.BEndNext                            = NULLINDEX;
      CI.AEndNext                            = NULLINDEX;

      CI.info.CI.contigID                    = NULLINDEX;
      CI.info.CI.source                      = NULLINDEX;

      CI.flags.bits.isCI                     = TRUE;  //  Crashes NextGraphNodeIterator() if not set.
      CI.flags.bits.isDead                   = TRUE;  //  DOA.

      //  Mark this as a unique CI, to prevent it from being used in stones.  Stones crashes on
      //  these empty unitigs.

      CI.flags.bits.isUnique                 = 1;
      CI.type                                = DISCRIMINATORUNIQUECHUNK_CGW;

      //  And while we're at it, just mark it chaff too.

      CI.flags.bits.isChaff                  = TRUE;

      //  Add it to the graph.

      SetGraphNode(ScaffoldGraph->CIGraph, &CI);

      //  Ensure that there are no edges, and that the edgeList is allocated.
      assert(ScaffoldGraph->CIGraph->edgeLists[CI.id].empty() == true);

      continue;
    }

    assert(i == uma->maID);

    if (1 != GetNumIntUnitigPoss(uma->u_list))
      fprintf(stderr, "ERROR:  Unitig %d has no placement; probably not run through consensus.\n", i);
    assert(1 == GetNumIntUnitigPoss(uma->u_list));

    if (i != GetIntUnitigPos(uma->u_list, 0)->ident)
      fprintf(stderr, "ERROR:  Unitig %d has incorrect unitig ident; error in utgcns or utgcnsfix or manual fixes to unitig.\n", i);
    assert(i == GetIntUnitigPos(uma->u_list, 0)->ident);

    MultiAlignT   *cma = CopyMultiAlignT(NULL, uma);

    ScaffoldGraph->tigStore->insertMultiAlign(cma, FALSE, TRUE);

    ProcessInputUnitig(uma);

    if ((++numUTG % 100000) == 0) {
      fprintf(stderr, "...processed "F_S32" unitigs.\n", numUTG);
      ScaffoldGraph->tigStore->flushCache();
    }
  }

  ScaffoldGraph->tigStore->flushCache();

  fprintf(stderr, "Checking sanity of loaded fragments.\n");
 
  uint32  numErrors = 0;

  for (int32 i=1, s=GetNumCIFragTs(ScaffoldGraph->CIFrags); i<s; i++) {
    CIFragT     *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, i);
 
    if (cifrag->flags.bits.isDeleted)
      continue;

    //  We could instead delete these fragments from the assembly.  That is somewhat difficult to do
    //  here, since we (a) don't have the gkpStore opened for writing and (b) don't have permission to
    //  do that anyway.

    if ((cifrag->cid == NULLINDEX) || (cifrag->CIid == NULLINDEX)) {
      fprintf(stderr, "ERROR:  Frag %d has null cid or CIid.  Fragment is not in an input unitig!\n", i);
      numErrors++;
    }
  }

  fprintf(stderr, "Processed %d unitigs with %d fragments.\n", numUTG, numFRG);
  fprintf(stderr, "  unique:        "F_U32" unitigs with total length "F_U64"\n", repeat_IsUnique.num,   repeat_IsUnique.len);
  fprintf(stderr, "  repeat:        "F_U32" unitigs with total length "F_U64"\n", repeat_IsRepeat.num,   repeat_IsRepeat.len);
  fprintf(stderr, "  few reads:     "F_U32" unitigs with total length "F_U64"\n", repeat_LowReads.num,   repeat_LowReads.len);
  fprintf(stderr, "  low cov stat:  "F_U32" unitigs with total length "F_U64"\n", repeat_LowCovStat.num, repeat_LowCovStat.len);
  fprintf(stderr, "  too short:     "F_U32" unitigs with total length "F_U64"\n", repeat_Short.num,      repeat_Short.len);
  fprintf(stderr, "  microhet:      "F_U32" unitigs with total length "F_U64"\n", repeat_MicroHet.num,   repeat_MicroHet.len);

  if (numErrors > 0)
    fprintf(stderr, "ERROR:  %u fragments are not in unitigs.\n", numErrors);
  assert(numErrors == 0);

  ScaffoldGraph->numLiveCIs     = GetNumGraphNodes(ScaffoldGraph->CIGraph);
  ScaffoldGraph->numOriginalCIs = GetNumGraphNodes(ScaffoldGraph->CIGraph);

  //  Load the distances.

  int32 numDists = ScaffoldGraph->gkpStore->gkStore_getNumLibraries();

  for (int32 i=1; i<=numDists; i++) {
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
    dist.numBad         = 0;
    dist.allowUpdate    = (gkpl->constantInsertSize == false);

    fprintf(stderr,"* Loaded dist %s,"F_CID" (%g +/- %g)\n",
            AS_UID_toString(gkpl->libraryUID), i, dist.mu, dist.sigma);

    SetDistT(ScaffoldGraph->Dists, i, &dist);
  }
}
