
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
static char *rcsid = "$Id: Output_CGW.c,v 1.54 2012-09-20 19:18:40 brianwalenz Exp $";

#include <assert.h>
#include <math.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ChiSquareTest_CGW.h"



//  Mark the trustedness of the intra-scaffold, inter-contig edges
//
void
MarkContigEdges(void) {

  CIScaffoldT *scaffold;
  GraphNodeIterator scaffolds;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    if(scaffold->type != REAL_SCAFFOLD)
      continue;
    //  OPERATES ON RAW EDGES
    MarkInternalEdgeStatus(ScaffoldGraph, scaffold, 0, FALSE, PAIRWISECHI2THRESHOLD_010, 100000000000.0);
  }

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    ContigT *contig;
    CIScaffoldTIterator Contigs;

    InitCIScaffoldTIterator(ScaffoldGraph, scaffold,TRUE, FALSE, &Contigs);
    while((contig = NextCIScaffoldTIterator(&Contigs)) != NULL){
      GraphEdgeIterator edges(ScaffoldGraph->ContigGraph, contig->id, ALL_END, ALL_EDGES);
      EdgeCGW_T        *edge;

      while((edge = edges.nextRaw()) != NULL){
        ContigT *mcontig;

        assert(edge->flags.bits.isRaw);
        if((edge->idA != contig->id) || isSingletonOverlapEdge(edge))
          continue;

        mcontig = GetGraphNode(ScaffoldGraph->ContigGraph, edge->idB);

        if(contig->scaffoldID != mcontig->scaffoldID)
          SetEdgeStatus(ScaffoldGraph->ContigGraph, edge, INTER_SCAFFOLD_EDGE_STATUS);

        PropagateEdgeStatusToFrag(ScaffoldGraph->ContigGraph, edge);
      }
    }
  }
}





void
OutputUnitigsFromMultiAligns(void) {
  GraphNodeIterator   nodes;
  ContigT	     *ci;

  InitGraphNodeIterator(&nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((ci = NextGraphNodeIterator(&nodes)) != NULL) {
    UnitigStatus   status = AS_UNASSIGNED;

    assert(ci->id >= 0);
    assert(ci->id < GetNumGraphNodes(ScaffoldGraph->CIGraph));

    if (ci->flags.bits.isChaff)
      continue;

    if (ci->type == RESOLVEDREPEATCHUNK_CGW)
      continue;

    if      (ci->type == DISCRIMINATORUNIQUECHUNK_CGW)
      status = AS_UNIQUE;

    else if (ci->type == UNIQUECHUNK_CGW)
      status = AS_UNIQUE;

    else if (ci->type == UNRESOLVEDCHUNK_CGW)
      if (ci->info.CI.numInstances > 0)
        status = AS_SEP;
      else if (ci->scaffoldID != NULLINDEX)
        status = AS_UNIQUE;
      else
        status = AS_NOTREZ;

    else {
      fprintf(stderr, "Unknown ci->type %d\n", ci->type);
      assert(0);
    }

    ScaffoldGraph->tigStore->setUnitigStatus(ci->id, status);
  }
}









void
OutputContigsFromMultiAligns(int32 outputFragsPerPartition) {
  GraphNodeIterator      nodes;
  ContigT		*ctg;

  uint32                 partitionNum     = 1;
  uint32                 partitionLimit   = outputFragsPerPartition;
  uint32                 partitionSize    = 0;
  uint32                 partitionContigs = 0;

  uint32                *partmap  = new uint32 [ScaffoldGraph->tigStore->numContigs() + 1];

  char                   partName[FILENAME_MAX];
  FILE                  *part;
  FILE                  *pari;

  //  Figure out how many fragments to place in each partition.  We can generally do this better
  //  than runCA, since we can ignore singletons.....too bad we don't know the correct parameters
  //  (being the number of partitions desired (128) and the minimum partition size (75000)).

#if 0
  if (partitionLimit == 0) {
    InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
    while ((ctg = NextGraphNodeIterator(&nodes)) != NULL) {
      if (ctg->flags.bits.isChaff)
        continue;

      partitionLimit += ScaffoldGraph->tigStore->getNumFrags(ctg->id, FALSE);
    }

    partitionLimit = partitionLimit / 128 + 1;

    if (partitionLimit < 75000)
      partitionLimit = 75000;
  }
#endif

  //  Open the partition mapping output file -- we do NOT fail here, as we're so close to finishing
  //  scaffolding.  Better to let runCA fail, and have a complete scaffolder output.

  sprintf(partName, "%s.partitioning", GlobalData->outputPrefix);
  errno = 0;
  part = fopen(partName, "w");
  if (errno) {
    fprintf(stderr, "Couldn't open '%s' for writing; NO CONTIGS WRITTEN: %s\n", partName, strerror(errno));
    return;
  }

  sprintf(partName, "%s.partitionInfo", GlobalData->outputPrefix);
  errno = 0;
  pari = fopen(partName, "w");
  if (errno) {
    fprintf(stderr, "Couldn't open '%s' for writing; NO CONTIGS WRITTEN: %s\n", partName, strerror(errno));
    return;
  }

  //  Build the partition mapping

  memset(partmap, 0xff, sizeof(uint32) * (ScaffoldGraph->tigStore->numContigs() + 1));

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while ((ctg = NextGraphNodeIterator(&nodes)) != NULL) {
    if (ctg->flags.bits.isChaff)
      continue;

    MultiAlignT   *ma        = ScaffoldGraph->tigStore->loadMultiAlign(ctg->id, FALSE);

    uint32         numFrag   = GetNumIntMultiPoss(ma->f_list);
    uint32         numUnitig = GetNumIntUnitigPoss(ma->u_list);

    IntMultiPos   *imp       = GetIntMultiPos(ma->f_list, 0);
    IntUnitigPos  *iup       = GetIntUnitigPos(ma->u_list, 0);


    if ((partitionSize + numFrag >= partitionLimit) &&
        (partitionSize           >  0)) {
      fprintf(pari, "Partition %d has %d contigs and %d fragments.\n",
              partitionNum, partitionContigs, partitionSize);

      partitionNum++;
      partitionContigs = 0;
      partitionSize    = 0;
    }

    for (uint32 i=0; i<numFrag; i++)
      fprintf(part, "%d\t%d\n", partitionNum, imp[i].ident);

    //fprintf(stderr, "contig %d into partition %d\n", ma->maID, partitionNum);

    assert(ma->maID < ScaffoldGraph->tigStore->numContigs() + 1);
    partmap[ma->maID]  = partitionNum;

    partitionContigs += 1;
    partitionSize    += numFrag;
  }

  fprintf(pari, "Partition %d has %d contigs and %d fragments.\n",
          partitionNum, partitionContigs, partitionSize);

  fclose(part);
  fclose(pari);

  //  Reset the tigStore for partitioning

  ScaffoldGraph->tigStore->writeToPartitioned(NULL, 0, partmap, ScaffoldGraph->tigStore->numContigs() + 1);

  //  Finalize contigs, remove existing alignments and write to the partitioned store.

  InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while ((ctg = NextGraphNodeIterator(&nodes)) != NULL) {
    if (ctg->flags.bits.isChaff)
      continue;

    CIScaffoldT   *scaffold  = GetGraphNode(ScaffoldGraph->ScaffoldGraph, ctg->scaffoldID);

    MultiAlignT   *ma        = ScaffoldGraph->tigStore->loadMultiAlign(ctg->id, FALSE);

    uint32         numFrag   = GetNumIntMultiPoss(ma->f_list);
    uint32         numUnitig = GetNumIntUnitigPoss(ma->u_list);

    IntMultiPos   *imp       = GetIntMultiPos(ma->f_list, 0);
    IntUnitigPos  *iup       = GetIntUnitigPos(ma->u_list, 0);

    NodeCGW_T     *unitig    = GetGraphNode(ScaffoldGraph->CIGraph, iup[0].ident);

    if ((numUnitig == 1) &&
        (ctg->scaffoldID == NULLINDEX) &&
        (unitig->info.CI.numInstances > 0))
      //  Skip surrogate instances.
      continue;

    for (uint32 i=0; i<numUnitig; i++) {
      NodeCGW_T    *unitig        = GetGraphNode(ScaffoldGraph->CIGraph, iup[i].ident);
      int32         num_instances = unitig->info.CI.numInstances;
 
      if (unitig->type == DISCRIMINATORUNIQUECHUNK_CGW) {
        iup[i].type = AS_UNIQUE_UNITIG;

      } else if (unitig->scaffoldID == NULLINDEX) {
        iup[i].type = AS_SINGLE_UNITIG;

      } else if (unitig->flags.bits.isSurrogate == FALSE) {
        iup[i].type = AS_ROCK_UNITIG;

      } else  if (unitig->flags.bits.isStoneSurrogate) {
        iup[i].type = AS_STONE_UNITIG;

      } else {
        iup[i].type = AS_PEBBLE_UNITIG;
      }

      iup[i].delta_length = 0;
      iup[i].delta        = NULL;

      if (unitig->type == RESOLVEDREPEATCHUNK_CGW) {
        iup[i].ident         = unitig->info.CI.baseID; // map back to the parent of this instance
        iup[i].num_instances = GetGraphNode(ScaffoldGraph->CIGraph, unitig->info.CI.baseID)->info.CI.numInstances;
      }

      //fprintf(stderr, "CTG %d UTG %d %d-%d\n",
      //        ma->maID, iup[i].ident, iup[i].position.bgn, iup[i].position.end);
    }

    ResetVA_char(ma->consensus);
    ResetVA_char(ma->quality);

    //  Important: we want to keep this ma in the cache, only so that we don't have to explicitly
    //  delete it right here.  When the store is deleted (below) the cache will get flushed.

    ScaffoldGraph->tigStore->setContigStatus(ctg->id, ((scaffold != NULL) && (scaffold->type == REAL_SCAFFOLD)) ? AS_PLACED : AS_UNPLACED);
    ScaffoldGraph->tigStore->insertMultiAlign(ma, FALSE, TRUE);
  }

  //  Flush and close the tigStore.  We are no longer able to do anything with the scaffold graph.

  delete ScaffoldGraph->tigStore;
  ScaffoldGraph->tigStore = NULL;

  delete [] partmap;
}
