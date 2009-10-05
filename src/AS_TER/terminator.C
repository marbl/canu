
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

const char *mainid = "$Id: terminator.C,v 1.1 2009-10-05 22:49:42 brianwalenz Exp $";

//  Assembly terminator module. It is the backend of the assembly pipeline and replaces internal
//  accession numbers by external accession numbers.

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <unistd.h>
#include  <assert.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "SYS_UIDclient.h"

#include "AS_PER_gkpStore.h"
#include "MultiAlignStore.h"

#include "Globals_CGW.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "ChiSquareTest_CGW.h"



class IIDtoUIDmap {
public:
  IIDtoUIDmap() {
    len = 0;
    max = 1048576;
    map = new AS_UID [max];
    memset(map, 0, sizeof(AS_UID) * max);
  };

  ~IIDtoUIDmap() {
    delete [] map;
  };

  void   add(AS_IID iid, AS_UID uid) {
    if (iid >= max) {
      max *= 2;
      AS_UID *M = new AS_UID [max];
      memcpy(M, map, sizeof(AS_UID) * len);
      memset(M + len, 0, sizeof(AS_UID) * (max - len));
      delete [] map;
      map = M;
    }
    if (iid >= len) {
      len = iid+1;
    }
    map[iid] = uid;
  };

  AS_UID lookup(AS_IID iid) {
    if (iid < len)
      return(map[iid]);
    fprintf(stderr, "Unknown IID %d\n", iid);  //  What a terrible error message
    exit(1);
  };

  void   dump(char *label, FILE *F) {
    for(int32 i=0; i<=len; i++) {
      if (AS_UID_isDefined(map[i]))
        fprintf(F,"%s\t%d\t%s\n", label, i, AS_UID_toString(map[i]));
    }
  };

private:
  int32    len;
  int32    max;
  AS_UID  *map;
};



UIDserver     *uidServer = NULL;
uint64         uidMin    = 0;

IIDtoUIDmap    MDImap;
IIDtoUIDmap    FRGmap;
IIDtoUIDmap    UTGmap;
IIDtoUIDmap    CCOmap;
IIDtoUIDmap    SCFmap;



void
writeMDI(FILE *asmFile, bool doWrite) {
  GenericMesg           pmesg;
  SnapMateDistMesg      mdi;

  assert(ScaffoldGraph->doRezOnContigs);

  pmesg.m = &mdi;
  pmesg.t = MESG_MDI;

  for (int32 i=1; i<GetNumDistTs(ScaffoldGraph->Dists); i++){
    DistT *dptr = GetDistT(ScaffoldGraph->Dists, i);

    //  Believe whatever estimate is here.  We used to reset to zero and the input (except we had
    //  already munged the input stddev) if there were 30 or fewer samples.

    mdi.erefines    = ScaffoldGraph->gkpStore->gkStore_getLibrary(i)->libraryUID;
    mdi.irefines    = i;
    mdi.mean        = dptr->mu;
    mdi.stddev      = dptr->sigma;
    mdi.min         = INT32_MIN;
    mdi.max         = INT32_MAX;
    mdi.num_buckets = 0;
    mdi.histogram   = NULL;

    //  The histogram does not get stored in a checkpoint.  If the current run of CGW did not have
    //  enough samples to recompute the histogram, we have to live without it

    if (dptr->bnum > 0) {
      mdi.min         = dptr->min;
      mdi.max         = dptr->max;
      mdi.num_buckets = dptr->bnum;
      mdi.histogram   = dptr->histogram;
    }

    if (doWrite)
      WriteProtoMesg_AS(asmFile, &pmesg);

    MDImap.add(mdi.irefines, mdi.erefines);

    if ((AS_UID_isString(mdi.erefines) == FALSE) &&
        (uidMin <= AS_UID_toInteger(mdi.erefines)))
      uidMin = AS_UID_toInteger(mdi.erefines) + 1;

    safe_free(dptr->histogram);

    dptr->histogram  = NULL;
    dptr->numSamples = 0;
    dptr->bnum       = 0;
  }
}



void
writeAFG(FILE *asmFile, bool doWrite) {
  GenericMesg       pmesg;
  AugFragMesg       afg;

  pmesg.m = &afg;
  pmesg.t = MESG_AFG;

  gkFragment     fr;
  gkStream      *fs = new gkStream(ScaffoldGraph->gkpStore, 0, 0, GKFRAGMENT_INF);

  for (int32 i=1; i<GetNumCIFragTs(ScaffoldGraph->CIFrags); i++) {
    CIFragT  *cifrag = GetCIFragT(ScaffoldGraph->CIFrags, i);

    fs->next(&fr);

    if (cifrag->flags.bits.isDeleted)
      continue;

    assert(cifrag->read_iid == i);
    assert(cifrag->read_iid == fr.gkFragment_getReadIID());

    afg.eaccession     = fr.gkFragment_getReadUID();
    afg.iaccession     = i;
    afg.mate_status    = cifrag->flags.bits.mateDetail;
    afg.chaff          = cifrag->flags.bits.isChaff;
    afg.clear_rng.bgn  = fr.gkFragment_getClearRegionBegin();
    afg.clear_rng.end  = fr.gkFragment_getClearRegionEnd  ();

    if (doWrite)
      WriteProtoMesg_AS(asmFile, &pmesg);

    FRGmap.add(afg.iaccession, afg.eaccession);

    if ((AS_UID_isString(afg.eaccession) == FALSE) &&
        (uidMin <= AS_UID_toInteger(afg.eaccession)))
      uidMin = AS_UID_toInteger(afg.eaccession) + 1;
  }

  delete fs;
}



void
writeAMP(FILE *asmFile, bool doWrite) {
  GenericMesg           pmesg;
  AugMatePairMesg       amp;

  pmesg.m = &amp;
  pmesg.t = MESG_AMP;

  for (int32 i=1; i<GetNumCIFragTs(ScaffoldGraph->CIFrags); i++) {
    CIFragT            *cif1 = GetCIFragT(ScaffoldGraph->CIFrags, i);
    CIFragT            *cif2 = NULL;

    if (cif1->flags.bits.isDeleted)
      continue;

    if (cif1->mate_iid == 0)
      continue;

    cif2 = GetCIFragT(ScaffoldGraph->CIFrags, cif1->mate_iid);

    if (cif2->flags.bits.isDeleted)
      continue;

    if (cif1->read_iid > cif2->read_iid)
      continue;

    assert(cif1->flags.bits.edgeStatus == cif2->flags.bits.edgeStatus);
    assert(cif1->flags.bits.mateDetail == cif2->flags.bits.mateDetail);

    amp.fragment1   = FRGmap.lookup(cif1->read_iid);
    amp.fragment2   = FRGmap.lookup(cif2->read_iid);
    amp.mate_status = cif1->flags.bits.mateDetail;

    if (doWrite)
      WriteProtoMesg_AS(asmFile, &pmesg);
  }
}



void
writeUTG(FILE *asmFile, bool doWrite) {
  GenericMesg         pmesg;
  SnapUnitigMesg      utg;

  pmesg.m = &utg;
  pmesg.t = MESG_UTG;

  GraphNodeIterator   unitigs;
  ChunkInstanceT     *ci;

  InitGraphNodeIterator(&unitigs, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((ci = NextGraphNodeIterator(&unitigs)) != NULL) {
    assert(ci->id >= 0);
    assert(ci->id < GetNumGraphNodes(ScaffoldGraph->CIGraph));

    if (ci->flags.bits.isChaff)
      //  Don't write chaff
      continue;

    if (ci->type == RESOLVEDREPEATCHUNK_CGW)
      //  Don't write surrogate instances
      continue;

    MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(ci->id, TRUE);

    utg.eaccession    = AS_UID_fromInteger(getUID(uidServer));
    utg.iaccession    = ci->id;
    utg.coverage_stat = ScaffoldGraph->tigStore->getUnitigCoverageStat(ci->id);
    utg.microhet_prob = ScaffoldGraph->tigStore->getUnitigMicroHetProb(ci->id);
    utg.status        = ScaffoldGraph->tigStore->getUnitigStatus(ci->id);
    utg.length        = GetMultiAlignLength(ma);
    utg.consensus     = Getchar(ma->consensus, 0);
    utg.quality       = Getchar(ma->quality, 0);
    utg.forced        = 0;
    utg.num_frags     = GetNumIntMultiPoss(ma->f_list);
    utg.num_vars      = 0;
    utg.f_list        = (SnapMultiPos*)safe_malloc(utg.num_frags * sizeof(SnapMultiPos));
    utg.v_list        = NULL;

    for (int32 i=0; i<utg.num_frags; i++) {
      IntMultiPos  *imp = GetIntMultiPos(ma->f_list, i);

      utg.f_list[i].type          = imp->type;
      utg.f_list[i].eident        = FRGmap.lookup(imp->ident);
      utg.f_list[i].position      = imp->position;
      utg.f_list[i].delta_length  = imp->delta_length;
      utg.f_list[i].delta         = imp->delta;
    }

    if (doWrite)
      WriteProtoMesg_AS(asmFile, &pmesg);

    safe_free(utg.f_list);

    UTGmap.add(utg.iaccession, utg.eaccession);
  }
}



void
writeULK(FILE *asmFile, bool doWrite) {
  GenericMesg          pmesg;
  SnapUnitigLinkMesg   ulk;

  pmesg.m = &ulk;
  pmesg.t = MESG_ULK;

  GraphNodeIterator  nodes;
  ChunkInstanceT    *ci;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((ci = NextGraphNodeIterator(&nodes)) != NULL) {
    assert(ci->type != CONTIG_CGW);

    if (ci->type == RESOLVEDREPEATCHUNK_CGW)
      continue;

    if (ci->flags.bits.isChaff)
      continue;

    GraphEdgeIterator	edges;
    CIEdgeT		*edge, *redge;

    InitGraphEdgeIterator(ScaffoldGraph->CIGraph, ci->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
    while ((edge = NextGraphEdgeIterator(&edges)) != NULL) {

      if (edge->idA != ci->id ||
          edge->flags.bits.isInferred ||
          edge->flags.bits.isInferredRemoved ||
          edge->flags.bits.isMarkedForDeletion)
        continue;

      ChunkInstanceT *mi = GetGraphNode(ScaffoldGraph->CIGraph, edge->idB);

      if (mi->flags.bits.isChaff)
        continue;

      ulk.eunitig1 = UTGmap.lookup(edge->idA);  //  == ci->id
      ulk.eunitig2 = UTGmap.lookup(edge->idB);

      ulk.orientation = edge->orient;  //  Don't need to map orientation, always using canonical orientation

      ulk.overlap_type = (isOverlapEdge(edge)) ? AS_OVERLAP : AS_NO_OVERLAP;

      ulk.is_possible_chimera = edge->flags.bits.isPossibleChimera;
      ulk.includes_guide      = FALSE;
      ulk.mean_distance       = edge->distance.mean;
      ulk.std_deviation       = sqrt(edge->distance.variance);
      ulk.num_contributing    = edge->edgesContributing;
      ulk.status              = AS_UNKNOWN_IN_ASSEMBLY;

      uint32  edgeCount = 0;
      uint32  edgeTotal = ulk.num_contributing;

      if ((edgeTotal == 1) && (ulk.overlap_type == AS_OVERLAP))
	// don't output pure overlap edges
        continue;

      // Look through the fragment pairs in this edge to decide the status of the link.

      redge = (edge->flags.bits.isRaw) ? edge : GetGraphEdge(ScaffoldGraph->CIGraph, edge->nextRawEdge);

      int numBad     = 0;
      int numGood    = 0;
      int numUnknown = 0;

      for (; redge != NULL; redge = GetGraphEdge(ScaffoldGraph->CIGraph, redge->nextRawEdge)) {
        if(isOverlapEdge(redge))
          continue;

        CIFragT *fragA = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragA);
        CIFragT *fragB = GetCIFragT(ScaffoldGraph->CIFrags, redge->fragB);

        assert(fragA->flags.bits.edgeStatus == fragB->flags.bits.edgeStatus);

        if ((fragA->flags.bits.edgeStatus == UNTRUSTED_EDGE_STATUS) ||
            (fragA->flags.bits.edgeStatus == TENTATIVE_UNTRUSTED_EDGE_STATUS))
          numBad++;

        else if ((fragA->flags.bits.edgeStatus == TRUSTED_EDGE_STATUS) ||
                 (fragA->flags.bits.edgeStatus == TENTATIVE_TRUSTED_EDGE_STATUS))
          numGood++;

        else
          numUnknown++;
        }

      if (numBad > 0)
        ulk.status = AS_BAD;

      else if (numGood > 0)
        ulk.status = AS_IN_ASSEMBLY;

      else
        ulk.status = AS_UNKNOWN_IN_ASSEMBLY;

      ulk.jump_list = (SnapMate_Pairs *)safe_malloc(sizeof(SnapMate_Pairs) * edgeTotal);

      if (edge->flags.bits.isRaw) {
        assert(edgeTotal == 1);

        ulk.jump_list[edgeCount].in1  = FRGmap.lookup(edge->fragA);
        ulk.jump_list[edgeCount].in2  = FRGmap.lookup(edge->fragB);
        ulk.jump_list[edgeCount].type = AS_MATE;

        edgeCount++;
      } else {
        assert(edgeTotal > 0);

        redge = edge;

        assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge

        while (redge->nextRawEdge != NULLINDEX) {
          redge = GetGraphEdge(ScaffoldGraph->CIGraph, redge->nextRawEdge);

          if (isOverlapEdge(redge)) {
            // overlap edges don't count
            edgeTotal--;
            continue;
          }

          ulk.jump_list[edgeCount].in1  = FRGmap.lookup(redge->fragA);
          ulk.jump_list[edgeCount].in2  = FRGmap.lookup(redge->fragB);
          ulk.jump_list[edgeCount].type = AS_MATE;

          edgeCount++;
        }
      }

      assert(edgeCount == edgeTotal);

      if (doWrite)
        WriteProtoMesg_AS(asmFile, &pmesg);

      safe_free(ulk.jump_list);
    }
  }
}



void
writeCCO(FILE *asmFile, bool doWrite) {
  GenericMesg         pmesg;
  SnapConConMesg      cco;

  pmesg.m = &cco;
  pmesg.t = MESG_CCO;

  GraphNodeIterator   contigs;
  ContigT	     *contig;

  InitGraphNodeIterator(&contigs, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while ((contig = NextGraphNodeIterator(&contigs)) != NULL) {
    assert(contig->id >= 0);
    assert(contig->id < GetNumGraphNodes(ScaffoldGraph->ContigGraph));

    if (contig->flags.bits.isChaff)
      continue;

    NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

    if ((ScaffoldGraph->tigStore->getNumUnitigs(contig->id, FALSE) == 1) &&
        (contig->scaffoldID == NULLINDEX) &&
        (unitig->info.CI.numInstances > 0))
      //  Contig is a surrogate instance
      continue;

    MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);

    cco.eaccession  = AS_UID_fromInteger(getUID(uidServer));
    cco.iaccession  = contig->id;
    cco.placed      = ScaffoldGraph->tigStore->getContigStatus(contig->id);
    cco.length      = GetMultiAlignLength(ma);
    cco.consensus   = Getchar(ma->consensus, 0);
    cco.quality     = Getchar(ma->quality, 0);
    cco.forced      = 0;
    cco.num_pieces  = GetNumIntMultiPoss(ma->f_list);
    cco.num_unitigs = GetNumIntMultiPoss(ma->u_list);
    cco.num_vars    = GetNumIntMultiPoss(ma->v_list);
    cco.pieces      = NULL;
    cco.unitigs     = NULL;
    cco.vars        = NULL;

    if (cco.num_pieces > 0) {
      cco.pieces = (SnapMultiPos *)safe_malloc(cco.num_pieces * sizeof(SnapMultiPos));

      for(int32 i=0; i<cco.num_pieces; i++) {
        IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

        cco.pieces[i].type         = imp->type;
        cco.pieces[i].eident       = FRGmap.lookup(imp->ident);
        cco.pieces[i].delta_length = imp->delta_length;
        cco.pieces[i].position     = imp->position;
        cco.pieces[i].delta        = imp->delta;
      }
    }

    if (cco.num_unitigs > 0) {
      cco.unitigs = (UnitigPos *)safe_malloc(cco.num_unitigs * sizeof(UnitigPos));

      for(int32 i=0; i<cco.num_unitigs; i++) {
        IntUnitigPos *imp = GetIntUnitigPos(ma->u_list, i);

        cco.unitigs[i].type         = imp->type;
        cco.unitigs[i].eident       = UTGmap.lookup(imp->ident);
        cco.unitigs[i].position     = imp->position;
        cco.unitigs[i].delta        = imp->delta;
        cco.unitigs[i].delta_length = imp->delta_length;
      }
    }

    if (cco.num_vars > 0) {
      cco.vars = (IntMultiVar *)safe_malloc(cco.num_vars * sizeof(IntMultiVar));

      for(int32 i=0; i<cco.num_vars; i++) {
        IntMultiVar *imv = GetIntMultiVar(ma->v_list, i);

        cco.vars[i].var_id                = imv->var_id;
        cco.vars[i].phased_id             = imv->phased_id;

        cco.vars[i].position              = imv->position;
        cco.vars[i].num_reads             = imv->num_reads;
        cco.vars[i].num_alleles           = imv->num_alleles;
        cco.vars[i].num_alleles_confirmed = imv->num_alleles_confirmed;
        cco.vars[i].min_anchor_size       = imv->min_anchor_size;
        cco.vars[i].var_length            = imv->var_length;

        cco.vars[i].alleles               = imv->alleles;
        cco.vars[i].var_seq_memory        = imv->var_seq_memory;
        cco.vars[i].read_id_memory        = imv->read_id_memory;

        cco.vars[i].enc_num_reads         = NULL;
        cco.vars[i].enc_weights           = NULL;
        cco.vars[i].enc_var_seq           = NULL;
        cco.vars[i].enc_read_ids          = NULL;
      }
    }

    if (doWrite)
      WriteProtoMesg_AS(asmFile, &pmesg);

    safe_free(cco.pieces);
    safe_free(cco.unitigs);
    safe_free(cco.vars);

    CCOmap.add(cco.iaccession, cco.eaccession);
  }
}



int
SurrogatedSingleUnitigContig(NodeCGW_T* contig) {

  if (contig->info.Contig.numCI > 1)
    //  Contig has multiple unitigs
    return(FALSE);

  if (contig->scaffoldID != NULLINDEX)
    //  Contig is placed
    return(FALSE);

  NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

  if (unitig->info.CI.numInstances == 0)
    //  Unitig has not been placed as a surrogate
    return(FALSE);

  //  Else, the unitig in this contig appears as a surrogate elsewhere in the assembly
  return(TRUE);
}



void
writeCLK(FILE *asmFile, bool doWrite) {
  GenericMesg			pmesg;
  SnapContigLinkMesg		clk;

  pmesg.m = &clk;
  pmesg.t = MESG_CLK;

  GraphNodeIterator      nodes;
  ContigT               *ctg;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while ((ctg = NextGraphNodeIterator(&nodes)) != NULL) {

    if (ctg->flags.bits.isChaff)
      continue;

    if (SurrogatedSingleUnitigContig(ctg))
      continue;

    GraphEdgeIterator	edges;
    CIEdgeT		*edge, *redge;

    InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, ctg->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
    while((edge = NextGraphEdgeIterator(&edges)) != NULL){

      if (edge->idA != ctg->id)
	continue;

      ContigT *mate = GetGraphNode(ScaffoldGraph->ContigGraph, edge->idB);

      if(mate->flags.bits.isChaff)
	continue;

      if (SurrogatedSingleUnitigContig(mate))
        continue;

      clk.econtig1 = CCOmap.lookup(edge->idA);
      clk.econtig2 = CCOmap.lookup(edge->idB);

      clk.orientation = edge->orient;  //  Don't need to map orientation, always using canonical orientation

      clk.overlap_type = (isOverlapEdge(edge)) ? AS_OVERLAP : AS_NO_OVERLAP;

      switch (GetEdgeStatus(edge)) {
        case LARGE_VARIANCE_EDGE_STATUS:
        case UNKNOWN_EDGE_STATUS:
        case INTER_SCAFFOLD_EDGE_STATUS:
          clk.status = AS_UNKNOWN_IN_ASSEMBLY;
          break;

        case TENTATIVE_TRUSTED_EDGE_STATUS:
        case TRUSTED_EDGE_STATUS:
          clk.status = AS_IN_ASSEMBLY;
          break;

        case TENTATIVE_UNTRUSTED_EDGE_STATUS:
        case UNTRUSTED_EDGE_STATUS:
          clk.status = AS_BAD;
          break;

        default:
          assert(0 /* Invalid edge status */);
      }

      clk.is_possible_chimera = edge->flags.bits.isPossibleChimera;
      clk.includes_guide      = FALSE;
      clk.mean_distance       = edge->distance.mean;
      clk.std_deviation       = sqrt(edge->distance.variance);
      clk.num_contributing    = edge->edgesContributing;

      uint32 edgeCount = 0;
      uint32 edgeTotal = clk.num_contributing;

      if ((edgeTotal == 1) &&
          (clk.overlap_type == AS_OVERLAP) &&
          (GlobalData->outputOverlapOnlyContigEdges == FALSE))
          // don't output pure overlap edges
	continue;

      clk.jump_list = (SnapMate_Pairs *)safe_malloc(sizeof(SnapMate_Pairs) * edgeTotal);

      if (edge->flags.bits.isRaw) {
	assert(edgeTotal == 1);

	if (clk.overlap_type == AS_NO_OVERLAP) {
	  clk.jump_list[edgeCount].in1  = FRGmap.lookup(edge->fragA);
	  clk.jump_list[edgeCount].in2  = FRGmap.lookup(edge->fragB);
          clk.jump_list[edgeCount].type = AS_MATE;
	} else {
	  assert(GlobalData->outputOverlapOnlyContigEdges);
	  clk.jump_list[edgeCount].in1  = AS_UID_undefined();
          clk.jump_list[edgeCount].in2  = AS_UID_undefined();
	  clk.jump_list[edgeCount].type = AS_UNMATED;
	}

        edgeCount++;

      } else {
	redge = edge;

	assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge

	while (redge->nextRawEdge != NULLINDEX) {
	  redge = GetGraphEdge(ScaffoldGraph->ContigGraph, redge->nextRawEdge);

	  if (isOverlapEdge(redge)) {
            // overlap edges don't count
            edgeTotal--;
	    continue;
          }

	  clk.jump_list[edgeCount].in1  = FRGmap.lookup(redge->fragA);
	  clk.jump_list[edgeCount].in2  = FRGmap.lookup(redge->fragB);
          clk.jump_list[edgeCount].type = AS_MATE;

          edgeCount++;
	}
      }

      assert(edgeCount == edgeTotal);

      if (doWrite)
        WriteProtoMesg_AS(asmFile, &pmesg);

      safe_free(clk.jump_list);
    }
  }
}



void
writeSCF(FILE *asmFile, bool doWrite) {
  GenericMesg			pmesg;
  SnapScaffoldMesg		scf;

  pmesg.m = &scf;
  pmesg.t = MESG_SCF;

  GraphNodeIterator   scaffolds;
  CIScaffoldT        *scaffold;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    if(scaffold->type != REAL_SCAFFOLD)
      continue;

    assert(scaffold->info.Scaffold.numElements > 0);

    scf.eaccession       = AS_UID_fromInteger(getUID(uidServer));
    scf.iaccession       = scaffold->id;
    scf.num_contig_pairs = scaffold->info.Scaffold.numElements - 1;
    scf.contig_pairs     = (SnapContigPairs *)safe_malloc(sizeof(SnapContigPairs) * scaffold->info.Scaffold.numElements);

    CIScaffoldTIterator	    contigs;
    ChunkInstanceT         *contigCurr;
    ChunkInstanceT         *contigLast;

    InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &contigs);
    contigLast = NextCIScaffoldTIterator(&contigs);

    ChunkOrient  orientLast = (contigLast->offsetAEnd.mean < contigLast->offsetBEnd.mean) ? A_B : B_A;
    ChunkOrient  orientCurr;

    assert(contigLast->scaffoldID == scaffold->id);

    if (scf.num_contig_pairs == 0) {
      scf.contig_pairs[0].econtig1 = CCOmap.lookup(contigLast->id);
      scf.contig_pairs[0].econtig2 = CCOmap.lookup(contigLast->id);
      scf.contig_pairs[0].mean     = 0.0;
      scf.contig_pairs[0].stddev   = 0.0;
      scf.contig_pairs[0].orient   = AB_AB; // got to put something

    } else {
      int32 pairCount = 0;

      while ((contigCurr = NextCIScaffoldTIterator(&contigs)) != NULL) {

        assert(pairCount < scf.num_contig_pairs);
        assert(contigCurr->scaffoldID == scaffold->id);

        scf.contig_pairs[pairCount].econtig1 = CCOmap.lookup(contigLast->id);
        scf.contig_pairs[pairCount].econtig2 = CCOmap.lookup(contigCurr->id);

        ChunkOrient orientCurr = (contigCurr->offsetAEnd.mean < contigCurr->offsetBEnd.mean) ? A_B : B_A;

        if (orientLast == A_B) {
          if (orientCurr == A_B) {
            scf.contig_pairs[pairCount].mean   = contigCurr->offsetAEnd.mean - contigLast->offsetBEnd.mean;
            scf.contig_pairs[pairCount].stddev = sqrt(contigCurr->offsetAEnd.variance -
                                                      contigLast->offsetBEnd.variance);
            scf.contig_pairs[pairCount].orient = AB_AB;
          } else {  //orientCurr == B_A
            scf.contig_pairs[pairCount].mean   = contigCurr->offsetBEnd.mean - contigLast->offsetBEnd.mean;
            scf.contig_pairs[pairCount].stddev = sqrt(contigCurr->offsetBEnd.variance -
                                                      contigLast->offsetBEnd.variance);
            scf.contig_pairs[pairCount].orient = AB_BA;
          }
        } else {  //orientLast == B_A
          if (orientCurr == A_B) {
            scf.contig_pairs[pairCount].mean   = contigCurr->offsetAEnd.mean - contigLast->offsetAEnd.mean;
            scf.contig_pairs[pairCount].stddev = sqrt(contigCurr->offsetAEnd.variance -
                                                      contigLast->offsetAEnd.variance);
            scf.contig_pairs[pairCount].orient = BA_AB;
          } else {  //orientCurr == B_A
            scf.contig_pairs[pairCount].mean   = contigCurr->offsetBEnd.mean - contigLast->offsetAEnd.mean;
            scf.contig_pairs[pairCount].stddev = sqrt(contigCurr->offsetBEnd.variance -
                                                      contigLast->offsetAEnd.variance);
            scf.contig_pairs[pairCount].orient = BA_BA;
          }
        }

        contigLast = contigCurr;
        orientLast = orientCurr;

        ++pairCount;
      }
    }

    if (doWrite)
      WriteProtoMesg_AS(asmFile, &pmesg);

    SCFmap.add(scf.iaccession, scf.eaccession);

    safe_free(scf.contig_pairs);
  }
}



void
writeSLK(FILE *asmFile, bool doWrite) {
  SnapScaffoldLinkMesg slk;
  GenericMesg          pmesg;

  pmesg.m = &slk;
  pmesg.t = MESG_SLK;

  GraphNodeIterator scaffolds;
  CIScaffoldT      *scaffold;
  CIScaffoldT      *scafmate;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    GraphEdgeIterator	 edges;
    CIEdgeT		*edge;
    CIEdgeT             *redge;

    InitGraphEdgeIterator(ScaffoldGraph->ScaffoldGraph, scaffold->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
    while((edge = NextGraphEdgeIterator(&edges)) != NULL) {
      if (edge->idA != scaffold->id)
        continue;

      scafmate = GetGraphNode(ScaffoldGraph->ScaffoldGraph, edge->idB);

      assert(!isOverlapEdge(edge));

      slk.escaffold1       = SCFmap.lookup(scaffold->id);
      slk.escaffold2       = SCFmap.lookup(scafmate->id);

      slk.orientation      = edge->orient;

      slk.includes_guide   = FALSE;
      slk.mean_distance    = edge->distance.mean;
      slk.std_deviation    = sqrt(edge->distance.variance);
      slk.num_contributing = edge->edgesContributing;

      int edgeTotal = slk.num_contributing;
      int edgeCount = 0;

      if(edgeTotal < 2)
        continue;

      slk.jump_list = (SnapMate_Pairs *)safe_malloc(sizeof(SnapMate_Pairs) * slk.num_contributing);

      if (edge->flags.bits.isRaw) {
        assert(edgeTotal <= 1);		// sanity check

        if (edgeTotal == 1) {
          slk.jump_list[edgeCount].in1 = FRGmap.lookup(edge->fragA);
          slk.jump_list[edgeCount].in2 = FRGmap.lookup(edge->fragB);
        }else{
          slk.jump_list[edgeCount].in1 = AS_UID_undefined();
          slk.jump_list[edgeCount].in2 = AS_UID_undefined();
        }

        slk.jump_list[edgeCount].type = AS_MATE;

        edgeCount++;

      } else {
        redge = edge;

        assert(redge->flags.bits.isRaw == FALSE);

        assert(redge->nextRawEdge != NULLINDEX); // must have >= 1 raw edge

        while (redge->nextRawEdge != NULLINDEX) {
          redge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph,redge->nextRawEdge);

          assert(!isOverlapEdge(redge));

          slk.jump_list[edgeCount].in1  = FRGmap.lookup(redge->fragA);
          slk.jump_list[edgeCount].in2  = FRGmap.lookup(redge->fragB);
          slk.jump_list[edgeCount].type = AS_MATE;

          edgeCount++;
        }
      }

      assert(edgeCount == edgeTotal);

      if (doWrite)
        WriteProtoMesg_AS(asmFile, &pmesg);

      safe_free(slk.jump_list);
    }
  }
}



int main (int argc, char *argv[]) {
  FILE       *asmFile                  = NULL;
  char       *outputPrefix             = NULL;
  char        outputName[FILENAME_MAX] = {0};
  int32       checkpointVers           = 0;
  int32       tigStoreVers             = 0;
  uint64      uidStart                 = 0;

  GlobalData = new Globals_CGW();

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->gkpStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      strcpy(GlobalData->tigStoreName, argv[++arg]);
      tigStoreVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->outputPrefix, argv[++arg]);
      checkpointVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-s") == 0) {
      uidStart = strtoul(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-n") == 0) {
      SYS_UIDset_euid_namespace(argv[++arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      SYS_UIDset_euid_server(argv[++arg]);

    } else if (strcmp(argv[arg], "-h") == 0) {
      err++;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((GlobalData->gkpStoreName[0] == 0) ||
      (GlobalData->tigStoreName[0] == 0) ||
      (GlobalData->outputPrefix[0] == 0) ||
      (err)) {
    fprintf(stderr, "usage: %s -g gkpStore [-o prefix] [-s firstUID] [-n namespace] [-E server] [-h]\n", argv[0]);
    fprintf(stderr, "  -g gkpStore             mandatory path to the gkpStore\n");
    fprintf(stderr, "  -t tigStore version     mandatory path to the tigStore and version\n");
    fprintf(stderr, "  -c checkpoint version   mandatory path to a checkpoint and version\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o prefix               write the output here\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -s firstUID      don't use real UIDs, but start counting from here\n");
    fprintf(stderr, "  -n namespace     use this UID namespace\n");
    fprintf(stderr, "  -E server        use this UID server\n");
    exit(1);
  }

  sprintf(outputName, "%s.asm", outputPrefix);
  errno = 0;
  asmFile = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "%s: Couldn't open '%s' for write: %s\n", outputName, strerror(errno)), exit(1);

  LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix, checkpointVers, FALSE);

  //  Reopen the tigStore used for consensus.
  delete ScaffoldGraph->tigStore;
  ScaffoldGraph->tigStore = new MultiAlignStore(GlobalData->tigStoreName, tigStoreVers, 0, 0, FALSE);

  fprintf(stderr, "Writing assembly file\n");

  writeMDI(asmFile, true);
  writeAFG(asmFile, true);

  uidServer = UIDserverInitialize(256, MAX(uidMin, uidStart));

  writeAMP(asmFile, true);
  writeUTG(asmFile, true);
  writeULK(asmFile, true);
  writeCCO(asmFile, true);
  writeCLK(asmFile, true);
  writeSCF(asmFile, true);
  writeSLK(asmFile, true);

  fclose(asmFile);

  fprintf(stderr, "Assembly file complete.\n");
  fprintf(stderr, "Writing IID to UID mapping files.\n");

  sprintf(outputName, "%s.iidtouid", outputPrefix);
  errno = 0;
  asmFile = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "%s: Couldn't open '%s' for write: %s\n", argv[0], outputName, strerror(errno)), exit(1);

  FRGmap.dump("FRG", asmFile);
  UTGmap.dump("UTG", asmFile);
  CCOmap.dump("CTG", asmFile);
  SCFmap.dump("SCF", asmFile);

  fclose(asmFile);

  fprintf(stderr, "IID to UID mapping files complete.\n");

  DestroyScaffoldGraph(ScaffoldGraph);

  return(0);
}
