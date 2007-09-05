
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
static char CM_ID[] = "$Id: Stats_CGW.c,v 1.15 2007-09-05 11:22:12 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "ScaffoldGraphIterator_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "GraphCGW_T.h"
#include "AS_UTL_interval.h"
#include "AS_UTL_fileIO.h"
#include "ChunkOverlap_CGW.h"
#include "MultiAlignment_CNS.h"

VA_DEF(int);

/* Generate statistics on the U-Unitig induced subgraph of
   the CIGraph */
void GenerateCIGraph_U_Stats(void){
  char buffer[2048];
  FILE *fout;
  GraphCGW_T *graph = ScaffoldGraph->CIGraph;
  VA_TYPE(int) *aEndDegree = CreateVA_int(GetNumGraphNodes(graph));
  VA_TYPE(int) *bEndDegree = CreateVA_int(GetNumGraphNodes(graph));
  GraphNodeIterator Nodes;
  NodeCGW_T *node;
  int unodes = 0;
  int uedges = 0;

  AS_UTL_mkdir("stat");

  fprintf(GlobalData->stderrc,"**** GenerateCIGraph_U_Stats ****\n");

  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    GraphEdgeIterator Edges;
    EdgeCGW_T *edge;
    int cnt = 0;

    if(node->type != DISCRIMINATORUNIQUECHUNK_CGW)
      continue;
    unodes++;

    InitGraphEdgeIterator(graph, node->id, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &Edges);
    while(NULL != (edge = NextGraphEdgeIterator(&Edges))){
      NodeCGW_T *other = GetGraphNode(graph,edge->idA == node->id? edge->idB: edge->idA);
      if(other->type != DISCRIMINATORUNIQUECHUNK_CGW)
	continue;
      if(edge->idA == node->id)
	uedges++;
      cnt++;
    }

    Setint(aEndDegree, node->id, &cnt);

    InitGraphEdgeIterator(graph, node->id, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &Edges);
    while(NULL != (edge = NextGraphEdgeIterator(&Edges))){
      NodeCGW_T *other = GetGraphNode(graph,edge->idA == node->id? edge->idB: edge->idA);
      if(other->type != DISCRIMINATORUNIQUECHUNK_CGW)
	continue;
      if(edge->idA == node->id)
	uedges++;
      cnt++;
    }

    Setint(bEndDegree, node->id, &cnt);
  }

  sprintf(buffer,"stat/CIGraph_U.nodeoutdegree.cgm");


  fout = fopen(buffer,"w");
  AssertPtr(fout);

  fprintf(fout,"CIGraph_U (V:%d  E:%d) OutDegrees\n",
	  unodes, uedges);

  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_UNIQUE_ONLY);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    int aDegree = *Getint(aEndDegree, node->id);
    int bDegree = *Getint(bEndDegree, node->id);
    fprintf(fout,"%d\n",aDegree + bDegree);
  }
  fclose(fout);

  sprintf(buffer,"stat/CIGraph_U.nodeendoutdegree.cgm");
  fout = fopen(buffer,"w");
  AssertPtr(fout);

  fprintf(fout,"CIGraph_U (V:%d  E:%d) End OutDegrees\n",
	  2*unodes, uedges);

  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_UNIQUE_ONLY);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    int aDegree = *Getint(aEndDegree, node->id);
    int bDegree = *Getint(bEndDegree, node->id);
    fprintf(fout,"%d\n",aDegree);
    fprintf(fout,"%d\n",bDegree);
  }
  fclose(fout);

  DeleteVA_int(aEndDegree);
  DeleteVA_int(bEndDegree);

}

/* Generate statistics on CIGraph */
void GenerateCIGraphStats(void){
  GraphCGW_T *graph = ScaffoldGraph->CIGraph;
  GraphNodeIterator Nodes;
  NodeCGW_T *node;
  int nu_unitigs = 0;
  int tfrags_nobf = 0;
  int tfrags_nolinks = 0;
  int tfrags_nobf_nolinks= 0;
  int nu_unitigs_no_bac_fragments = 0;
  int nu_unitigs_no_links = 0;
  int nu_unitigs_no_links_no_bac_fragments = 0;
  int n_unitigs = 0;

  fprintf(GlobalData->stderrc,"**** GenerateCIGraphStats ****\n");

  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    GraphEdgeIterator Edges;
    EdgeCGW_T *edge;
    int cnt = 0;

    // Filter surrogates
    if(node->info.CI.numFragments == 0)
      continue;

    n_unitigs++;
    if(node->flags.bits.isUnique)
      continue;
    nu_unitigs++;

    InitGraphEdgeIterator(graph, node->id, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &Edges);
    while(NULL != (edge = NextGraphEdgeIterator(&Edges))){
      cnt++;
    }

    //if(node->flags.bits.includesFinishedBacFragments){
    //  nu_unitigs_no_bac_fragments++;
    //  tfrags_nobf += node->info.CI.numFragments;
    //}
    if(cnt == 0){
      nu_unitigs_no_links++;
      tfrags_nolinks += node->info.CI.numFragments;
    }

    //if(cnt == 0 && node->flags.bits.includesFinishedBacFragments){
    //  nu_unitigs_no_links_no_bac_fragments++;
    //  tfrags_nobf_nolinks += node->info.CI.numFragments;
    //}
  }
  fprintf(GlobalData->stderrc,"*@ Graph has %d unitigs of which %d are non-unique\n",
	  n_unitigs, nu_unitigs);
  fprintf(GlobalData->stderrc,"*@ %d unitigs have no external data comprising %d total fragments\n",
	  nu_unitigs_no_bac_fragments, tfrags_nobf);
  fprintf(GlobalData->stderrc,"*@ %d unitigs have no links comprising %d total fragments\n",
	  nu_unitigs_no_links, tfrags_nolinks);
  fprintf(GlobalData->stderrc,"*@ %d unitigs have no external data AND no links comprising %d total fragments\n",
	  nu_unitigs_no_links_no_bac_fragments, tfrags_nobf_nolinks);
}


/* Generate stats for the subgraph of the Contig graph induced by placed contigs */
void GeneratePlacedContigGraphStats(char *label,int iteration){

  char buffer[2048];
  FILE *fout;
  FILE *fEndOut;
  FILE *fLengthOut;
  FILE *fUnitigOut;
  GraphCGW_T *graph = ScaffoldGraph->ContigGraph;
  VA_TYPE(int) *aEndDegree = CreateVA_int(GetNumGraphNodes(graph));
  VA_TYPE(int) *bEndDegree = CreateVA_int(GetNumGraphNodes(graph));
  GraphNodeIterator Nodes;
  NodeCGW_T *node;
  int unodes = 0;
  int uedges = 0;

  if(GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph) < 1){
    return;
  }

  AS_UTL_mkdir("stat");

  fprintf(GlobalData->stderrc,"**** GeneratePlacedContigStats %s %d****\n", label, iteration);

  
  sprintf(buffer,"stat/%s%d.PlacedContig.nodelength.cgm", label, iteration);
  fLengthOut = fopen(buffer,"w");
  AssertPtr(fLengthOut);

  sprintf(buffer,"stat/%s%d.PlacedContig.unitigs.cgm", label,iteration);
  fUnitigOut = fopen(buffer,"w");
  AssertPtr(fUnitigOut);

  sprintf(buffer,"stat/%s%d.PlacedContig.nodeoutdegree.cgm", label,iteration);
  fout = fopen(buffer,"w");
  AssertPtr(fout);

  sprintf(buffer,"stat/%s%d.PlacedContig.nodeendoutdegree.cgm",label,iteration);
  fEndOut = fopen(buffer,"w");
  AssertPtr(fEndOut);

  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    GraphEdgeIterator Edges;
    EdgeCGW_T *edge;
    int cnt = 0;

    if(node->scaffoldID == NULLINDEX){
      Setint(aEndDegree, node->id, &cnt);
      Setint(bEndDegree, node->id, &cnt);
      continue;
    }

    unodes++;

    InitGraphEdgeIterator(graph, node->id, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &Edges);
    while(NULL != (edge = NextGraphEdgeIterator(&Edges))){
      NodeCGW_T *other = GetGraphNode(graph,edge->idA == node->id? edge->idB: edge->idA);
      if(other->scaffoldID == NULLINDEX)
	continue;
      if(edge->idA == node->id)
	uedges++;
      cnt++;
    }

    Setint(aEndDegree, node->id, &cnt);

    cnt = 0;

    InitGraphEdgeIterator(graph, node->id, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &Edges);
    while(NULL != (edge = NextGraphEdgeIterator(&Edges))){
      NodeCGW_T *other = GetGraphNode(graph,edge->idA == node->id? edge->idB: edge->idA);
      if(other->scaffoldID == NULLINDEX)
	continue;
      if(edge->idA == node->id)
	uedges++;
      cnt++;
    }

    Setint(bEndDegree, node->id, &cnt);
  }


  fprintf(fout,"%s (V:%d E:%d) OutDegrees\n",
	  label, unodes, uedges);
  fprintf(fEndOut,"PContig (V:%d  E:%d) End OutDegrees\n",
	  2*unodes, uedges);
  fprintf(fLengthOut,"%s (V:%d E:%d) Length\n",
	  label, unodes, uedges);
  fprintf(fUnitigOut,"%s (V:%d E:%d) Unitigs\n",
	  label, unodes, uedges);

  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    int aDegree;
    int bDegree;

    if(node->scaffoldID == NULLINDEX)
      continue;

    aDegree = *Getint(aEndDegree, node->id);
    bDegree = *Getint(bEndDegree, node->id);

    fprintf(fLengthOut,"%d\n", (int)node->bpLength.mean);
    fprintf(fUnitigOut,"%d\n", node->info.Contig.numCI);
    fprintf(fout,"%d\n",aDegree + bDegree);
    fprintf(fEndOut,"%d\n",aDegree);
    fprintf(fEndOut,"%d\n",bDegree);
  }

  fclose(fout);
  fclose(fEndOut);
  fclose(fLengthOut);
  fclose(fUnitigOut);
  DeleteVA_int(aEndDegree);
  DeleteVA_int(bEndDegree);

}



/* Generate stats for the scaffold graph  */
/* # of scaffolds
   # scaffold edges (weight 2 or more?)
   # scaffold edges BAC only (weight 2 or more?)
   Scaffold lengths
   Contigs/Scaffold
   Intra-Scaffold Gaps
   Links/Scaffold edge w/o BACs   (weight 2 or more?)
   Links/Scaffold edge w/BACs     (weight 2 or more?)
   Nature of 2/10k scaffolding (output to text file)
*/
void GenerateScaffoldGraphStats(char *label, int iteration){

  char buffer[2048];
  FILE *fout;
  FILE *fEndOut;
  FILE *fLengthOut, *fLengthSinglesOut;
  FILE *fContigOut;
  FILE *fLinksPerEdge_WBacs;
  FILE *fLinksPerEdge_WOBacs;
  FILE *fScaffoldGapMeans, *fScaffoldGapStds;
  FILE *fScaffoldNature;
  GraphCGW_T *graph = ScaffoldGraph->ScaffoldGraph;
  VA_TYPE(int) *aEndDegree = CreateVA_int(GetNumGraphNodes(graph));
  VA_TYPE(int) *bEndDegree = CreateVA_int(GetNumGraphNodes(graph));
  GraphNodeIterator Nodes;
  NodeCGW_T *node;
  int nodes = 0;
  int edges = 0;
  int bacOnlyEdges = 0;
  int numberInferred = 0;
  int numberRemoved = 0;

  if(GetNumGraphNodes(graph) < 1){
    return;
  }

  AS_UTL_mkdir("stat");

  fprintf(GlobalData->stderrc,"**** GeneratePlacedContigStats %s %d ****\n", label,iteration);

  
  sprintf(buffer,"stat/%s%d.Scaffolds.nodelength.cgm", label,iteration);
  fLengthOut = fopen(buffer,"w");
  AssertPtr(fLengthOut);

  sprintf(buffer,"stat/%s%d.SingleScaffolds.nodelength.cgm", label,iteration);
  fLengthSinglesOut = fopen(buffer,"w");
  AssertPtr(fLengthSinglesOut);

  sprintf(buffer,"stat/%s%d.Scaffolds.Nature.txt", label,iteration);
  fScaffoldNature = fopen(buffer,"w");
  AssertPtr(fScaffoldNature);

  sprintf(buffer,"stat/%s%d.Scaffolds.contigs.cgm", label,iteration);
  fContigOut = fopen(buffer,"w");
  AssertPtr(fContigOut);

  sprintf(buffer,"stat/%s%d.Scaffolds.nodeoutdegree.cgm", label,iteration);
  fout = fopen(buffer,"w");
  AssertPtr(fout);

  sprintf(buffer,"stat/%s%d.Scaffolds.nodeendoutdegree.cgm", label,iteration);
  fEndOut = fopen(buffer,"w");
  AssertPtr(fEndOut);

  sprintf(buffer,"stat/%s%d.Scaffolds.links_per_edge_w_bac.cgm",label,iteration);
  fLinksPerEdge_WBacs = fopen(buffer,"w");
  AssertPtr(fLinksPerEdge_WBacs);

  sprintf(buffer,"stat/%s%d.Scaffolds.links_per_edge_wo_bac.cgm",label,iteration);
  fLinksPerEdge_WOBacs = fopen(buffer,"w");
  AssertPtr(fLinksPerEdge_WOBacs);

  sprintf(buffer,"stat/%s%d.Scaffolds.intra_scaffold_gap_means.cgm",label,iteration);
  fScaffoldGapMeans = fopen(buffer,"w");
  AssertPtr(fScaffoldGapMeans);

  sprintf(buffer,"stat/%s%d.Scaffolds.intra_scaffold_gap_stds.cgm",label,iteration);
  fScaffoldGapStds = fopen(buffer,"w");
  AssertPtr(fScaffoldGapStds);

  fprintf(fLinksPerEdge_WBacs, "Links Per Edge including BACs\n");
  fprintf(fLinksPerEdge_WOBacs, "Links Per Edge excluding BACs\n");

  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    GraphEdgeIterator Edges;
    EdgeCGW_T *edge;
    int cnt = 0;
    int noBacDegree = 0, totalDegree = 0;

    nodes++;

    InitGraphEdgeIterator(graph, node->id, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &Edges);
    while(NULL != (edge = NextGraphEdgeIterator(&Edges))){
      if(edge->idA == node->id){
	if(edge->flags.bits.isInferred){
	  numberInferred++;
	}
	if(edge->flags.bits.isTransitivelyRemoved){
	  numberRemoved++;
	}
	if(edge->edgesContributing > 1){
          edges++;
          totalDegree = EdgeDegree(graph,edge);
          assert(edge->edgesContributing == totalDegree);
          if(totalDegree && (totalDegree == noBacDegree))
            bacOnlyEdges++;
          fprintf(fLinksPerEdge_WBacs,  "%d\n",totalDegree);
          fprintf(fLinksPerEdge_WOBacs,  "%d\n",noBacDegree);
        }
        cnt++;
      }
    }

    Setint(aEndDegree, node->id, &cnt);

    cnt = 0;
    InitGraphEdgeIterator(graph, node->id, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &Edges);
    while(NULL != (edge = NextGraphEdgeIterator(&Edges))){
      if(edge->idA == node->id){
	if( edge->edgesContributing > 1){
          edges++;
          totalDegree = EdgeDegree(graph,edge);
          assert(edge->edgesContributing == totalDegree);
          if(totalDegree && (totalDegree == noBacDegree))
            bacOnlyEdges++;
          fprintf(fLinksPerEdge_WBacs,  "%d\n",totalDegree);
          fprintf(fLinksPerEdge_WOBacs,  "%d\n",noBacDegree);
	}
	cnt++;
      }
    }

    Setint(bEndDegree, node->id, &cnt);
  }

  fprintf(fScaffoldNature,"Number removed by transitive reduction: %d\n", numberRemoved);
  fprintf(fScaffoldNature,"Number Edges marked inferred by end-tuck: %d\n", numberInferred);

  fprintf(fout,"%s (V:%d E:%d (%d)) OutDegrees\n",
	  label, nodes, edges, bacOnlyEdges);
  fprintf(fEndOut,"Scaffolds (V:%d  E:%d (%d)) End OutDegrees\n",
	  2*nodes, edges, bacOnlyEdges);
  fprintf(fLengthOut,"Scaffold (V:%d E:%d (%d)) Length\n",
	  nodes, edges,bacOnlyEdges);
  fprintf(fLengthSinglesOut,"Scaffold (V:%d E:%d (%d)) Length\n",
	  nodes, edges,bacOnlyEdges);
  fprintf(fContigOut,"Scaffold (V:%d E:%d (%d)) Contigs\n",
	  nodes, edges, bacOnlyEdges);
  fprintf(fScaffoldGapMeans, "Mean IntraScaffold Gap Sizes\n");
  fprintf(fScaffoldGapStds, "Std IntraScaffold Gap Sizes\n");

  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){

    int aDegree = *Getint(aEndDegree, node->id);
    int bDegree = *Getint(bEndDegree, node->id);


    if(node->info.Scaffold.numElements > 1)
      fprintf(fLengthOut,"%d\n", (int)node->bpLength.mean);
    else
      fprintf(fLengthSinglesOut,"%d\n", (int)node->bpLength.mean);

    fprintf(fContigOut, "%d\n", node->info.Scaffold.numElements);
    fprintf(fout, "%d\n",aDegree + bDegree);
    fprintf(fEndOut, "%d\n",aDegree);
    fprintf(fEndOut, "%d\n",bDegree);

    {
      CIScaffoldTIterator Contigs;
      NodeCGW_T *prev = NULL;
      NodeCGW_T *next;
      double mean, std;

      InitCIScaffoldTIterator(ScaffoldGraph, node, TRUE, FALSE, &Contigs);
      for(next = NextCIScaffoldTIterator(&Contigs), prev = NULL; next != NULL; prev = next, next = NextCIScaffoldTIterator(&Contigs)){
	if(prev && next){
	  if(prev->offsetAEnd.mean < prev->offsetBEnd.mean){
	    if(next->offsetAEnd.mean < next->offsetBEnd.mean){
	      mean = next->offsetAEnd.mean - prev->offsetBEnd.mean;
	      std = next->offsetAEnd.variance - prev->offsetBEnd.variance;
	    }else{
	      mean = next->offsetBEnd.mean - prev->offsetBEnd.mean;
	      std = next->offsetBEnd.variance - prev->offsetBEnd.variance;
	    }
	  }else{
	    if(next->offsetAEnd.mean < next->offsetBEnd.mean){
	      mean = next->offsetAEnd.mean - prev->offsetAEnd.mean;
	      std = next->offsetAEnd.variance - prev->offsetAEnd.variance;
	    }else{
	      mean = next->offsetBEnd.mean - prev->offsetAEnd.mean;
	      std = next->offsetBEnd.variance - prev->offsetAEnd.variance;
	    }

	  }
	  if(std > 0)
	    std = sqrt(std);
	  else
	    std = 0.0;
	  fprintf(fScaffoldGapMeans, "%d\n", (int)mean);
	  fprintf(fScaffoldGapStds, "%d\n", (int)std);

	}

      }
    }
  }


  fclose(fout);
  fclose(fEndOut);
  fclose(fLengthOut);
  fclose(fLengthSinglesOut);
  fclose(fContigOut);
  fclose(fScaffoldNature);
  fclose(fScaffoldGapMeans);
  fclose(fScaffoldGapStds);
  DeleteVA_int(aEndDegree);
  DeleteVA_int(bEndDegree);

}



/* Compute statistics on links
   # of mates per link edge
   std of link edges by
   all
   overlap confirmed (non tandem)
   not confirmed (no overlap)
*/
void GenerateLinkStats(GraphCGW_T *graph, char *label, int iteration){
  int i;
  char buffer[2048];
  FILE *mates_per_link;
  FILE *linkstd_all;
  FILE *linkstd_w_overlap;
  FILE *linkstd_no_overlap;
  int mates;
  int cgbOverlap = 0;
  int nonCGBOverlap = 0;
  const char *graphName = (graph->type == CI_GRAPH?"CI":"Contig");

  AS_UTL_mkdir("stat");

  sprintf(buffer,"stat/%s%s%d.mates_per_link.cgm",graphName,label,iteration);
  mates_per_link = fopen(buffer,"w");
  AssertPtr(mates_per_link);
  fprintf(mates_per_link,"Mates per link %s\n", label);

  sprintf(buffer,"stat/%s%s%d.linkstd_all.cgm",graphName,label,iteration);
  linkstd_all = fopen(buffer,"w");
  AssertPtr(linkstd_all);
  fprintf(linkstd_all,"Link Standard Deviation -- All %s\n", label);

  sprintf(buffer,"stat/%s%s%d.linkstd_w_overlap.cgm",graphName,label,iteration);
  linkstd_w_overlap = fopen(buffer,"w");
  AssertPtr(linkstd_w_overlap);
  fprintf(linkstd_w_overlap,"Link Standard Deviation -- Confirmed %s\n", label);

  sprintf(buffer,"stat/%s%s%d.linkstd_no_overlap.cgm",graphName,label, iteration);
  linkstd_no_overlap = fopen(buffer,"w");
  AssertPtr(linkstd_no_overlap);
  fprintf(linkstd_no_overlap,"Link Standard Deviation -- Unconfirmed %s\n", label);



  for(i = 0; i < GetNumGraphEdges(graph); i++){
    EdgeCGW_T *edge = GetGraphEdge(graph,i);
    NodeCGW_T *nodeA, *nodeB;
    CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, edge);
    int std = (edge->distance.variance > 0.0 ? (int)sqrt(edge->distance.variance):-1);
    if(edge->flags.bits.isDeleted || edge->topLevelEdge != eid)
      continue;

    mates = edge->edgesContributing;
    nodeA = GetGraphNode(graph, edge->idA);
    nodeB = GetGraphNode(graph, edge->idB);

    if(graph->type == CONTIG_GRAPH &&
       (nodeA->scaffoldID == NULLINDEX ||
	nodeB->scaffoldID == NULLINDEX ))
      continue;

    if(isOverlapEdge(edge)){
      mates--;
      
      if(mates){  // we only want edges that are NOT overlap only
        fprintf(linkstd_w_overlap,"%d\n", std);
        if(graph->type == CONTIG_GRAPH){
          ChunkOverlapCheckT olap = {0};
          int overlapFound = LookupOverlap(graph, edge->idA, edge->idB, edge->orient, &olap);
          if(overlapFound && olap.fromCGB){
            cgbOverlap++;
          }else{
            nonCGBOverlap++;
          }
        }
      }
    }else{
      fprintf(linkstd_no_overlap,"%d\n", std);
    }

    if(mates){
      fprintf(linkstd_all,"%d\n", std);
      fprintf(mates_per_link,"%d\n", mates);
    }
  }
  if(graph->type == CONTIG_GRAPH)
    fprintf(GlobalData->stderrc,"*** Links confirmed by cgbOlaps %d  onCGBOlaps %d\n",
	    cgbOverlap, nonCGBOverlap);

  fclose(linkstd_all);
  fclose(linkstd_w_overlap);
  fclose(linkstd_no_overlap);
  fclose(mates_per_link);
}


int32 ApproximateUnitigCoverage(NodeCGW_T *unitig){
  CDS_COORD_t length=0;
  int i;
  MultiAlignT *ma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);

  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++){
    IntMultiPos *pos = GetIntMultiPos(ma->f_list, i);
    length += abs(pos->position.bgn - pos->position.end);
  }

  if(length == 0)
    return 0;

  return (int32) ( length / unitig->bpLength.mean);
}


void GenerateSurrogateStats(char *phase){
  FILE *surrogPer, *surrogSize, *surrogFrags, *surrogRatio;
  FILE *surrogCreated;
  char buffer[1024];
  GraphNodeIterator Nodes;
  NodeCGW_T *node;
  int stoneSurrogs=0, 
    walkSurrogs=0;

  AS_UTL_mkdir("stat");
  
  sprintf(buffer,"stat/%s.surrogates_per_repeatCI.cgm",phase);
  surrogPer = fopen(buffer,"w");
  AssertPtr(surrogPer);

  sprintf(buffer,"stat/%s.surrogates_Created.cgm", phase);
  surrogCreated = fopen(buffer,"w");
  AssertPtr(surrogCreated);

  fprintf(surrogPer,"Surrogates Per Repeat CI\n");

  sprintf(buffer,"stat/%s.surrogates_size.cgm", phase);
  surrogSize = fopen(buffer,"w");
  AssertPtr(surrogSize);

  fprintf(surrogSize,"Surrogate Sizes\n");

  sprintf(buffer,"stat/%s.surrogates_fragsPer.cgm", phase);
  surrogFrags = fopen(buffer,"w");
  AssertPtr(surrogFrags);

  fprintf(surrogFrags,"Fragments per Repeat CI\n");

  sprintf(buffer,"stat/%s.surrogates_ratio.cgm",phase);
  surrogRatio = fopen(buffer,"w");
  AssertPtr(surrogRatio);

  fprintf(surrogRatio,"Approx Cov / Num Copies\n");

  InitGraphNodeIterator(&Nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    char type;
    switch(node->type){
      case DISCRIMINATORUNIQUECHUNK_CGW:
        type = 'U';
        break;
      case UNIQUECHUNK_CGW:
        type = 'u';
        break;
      case UNRESOLVEDCHUNK_CGW:
        type = '?';
        break;
      case RESOLVEDREPEATCHUNK_CGW:
        type = 'R';
        if(node->flags.bits.isWalkSurrogate){
          walkSurrogs++;
        }else if(node->flags.bits.isStoneSurrogate){
          stoneSurrogs++;
        }
        break;
      default:
        type = 'X';
        break;
    }
#ifdef DEBUG_DETAILED
    fprintf(GlobalData->stderrc,"* Node " F_CID " %c contig:" F_CID "  numInstances %d\n", 
	    node->id, type, node->info.CI.contigID, node->info.CI.numInstances);
#endif
    if((node->type != UNRESOLVEDCHUNK_CGW) ||    // is not a surrogate parent
       (node->info.CI.numInstances == 0))        // has no surrogates
      continue;
    fprintf(surrogSize,"%d\n", (int)node->bpLength.mean);
    fprintf(surrogPer, "%d\n", (int)node->info.CI.numInstances);
    fprintf(surrogFrags, "%d\n", node->info.CI.numFragments);
    fprintf(surrogRatio, "%d\n", ApproximateUnitigCoverage(node)/node->info.CI.numInstances );
  }
  fprintf(GlobalData->stderrc,"* Stones: %d  Walks:%d\n",
	  stoneSurrogs, walkSurrogs);


  fclose(surrogCreated);
  fclose(surrogPer);
  fclose(surrogSize);
  fclose(surrogRatio);
  fclose(surrogFrags);
}

////////////////////////////////////////////////////////////////////////////////

typedef enum {
  CNS_INSERT = (int) 'I',
  CNS_DELETE = (int) 'D',    
  CNS_SUBSTITUTE = (int) 'S'
}ErrorType;

typedef struct{
  int position;
  ErrorType type;
}ErrorStruct;

VA_DEF(ErrorStruct)

static
void 
CollectStats(MultiAlignT *ma,
             GateKeeperStore *frag_store, 
             FILE *column_stats, 
             FILE *frag_stats,
             uint32 clrrng_flag)
{
  /*  
      Need to append to column_stats and frag_stats the following:
      To column_stats:
      Foreach column in multialignment, print contigID, column index, coverage, quality value

      To frag_stats:
      Foreach fragment in multialignment, print fragIID, fragUID,  clr_bgn, clr_end, errors (in apos,type pairs)
  */
  int32 readptr;
  int32 delptr;
  int32 left;
  int32 right;
  int32 flen;
  int32 ungapped=0;
  uint clrbgn;
  uint clrend;
  CDS_UID_t accession;
  IntMultiPos *reads=GetIntMultiPos(ma->f_list,0);
  int i; // tracks reads
  int j; // tracks columns
  int32 ma_len = GetMultiAlignLength(ma);
  int num_reads = GetNumIntMultiPoss(ma->f_list);
  int *column_cov;
  int *column_mm;
  char column_call;
  int num_errors;
  VA_TYPE(ErrorStruct) *errors; 
  ErrorStruct frag_error;
  char seqdata[AS_READ_MAX_LEN+2];
  char qvdata[AS_READ_MAX_LEN+2];
  fragRecord *rsp = new_fragRecord();
    
  column_cov = (int *) safe_malloc(ma_len*sizeof(int));
  column_mm = (int *) safe_malloc(ma_len*sizeof(int));
  errors = CreateVA_ErrorStruct(250);

  assert(column_cov && column_mm);

  // special case for singletons
  if (num_reads == 1) {
    getFrag(frag_store,reads[0].ident,rsp,FRAG_S_QLT);
    clrbgn = getFragRecordClearRegionBegin(rsp, clrrng_flag);
    clrend = getFragRecordClearRegionEnd  (rsp, clrrng_flag);
    fprintf(frag_stats,F_IID "  " F_UID " %c %d %d\n",
            reads[0].ident,accession,
            reads[0].type,(int) clrbgn,(int) clrend);
    flen = clrend - clrbgn;
    // capture only the clear range for analysis
    // reverse complement if necessary:
    memcpy(seqdata, getFragRecordSequence(rsp) + clrbgn, (flen+1)*sizeof(char));
    memcpy(qvdata,  getFragRecordQuality(rsp)  + clrbgn, (flen+1)*sizeof(char));
    seqdata[flen] = '\0';
    qvdata[flen] = '\0';
    for (j=0;j<ma_len;j++) {
      fprintf(column_stats,"%d %d %d %d %d %c %d %d %d\n",ma->maID,j,1,
              0,
              (reads[0].type != AS_READ &&
               reads[0].type != AS_EXTR &&
               reads[0].type != AS_TRNR)?1:0, 
              seqdata[j],qvdata[j] - '0',0,j);
    }
  } else {
        
    // initialize column coverage to zero
    for ( j=0;j<ma_len;j++) {
      column_cov[j]=0;
      column_mm[j]=0;
    }

    for(i=0;i<num_reads;i++) {
      left = (reads[i].position.bgn < reads[i].position.end)? reads[i].position.bgn : reads[i].position.end;
      right= (reads[i].position.bgn > reads[i].position.end)?reads[i].position.bgn:reads[i].position.end;
      getFrag(frag_store,reads[i].ident,rsp,FRAG_S_QLT);
      clrbgn = getFragRecordClearRegionBegin(rsp, clrrng_flag);
      clrend = getFragRecordClearRegionEnd  (rsp, clrrng_flag);
      flen = clrend - clrbgn;
      assert(flen < AS_READ_MAX_LEN);
      assert(flen > 0);
      // capture only the clear range for analysis
      // reverse complement if necessary:
      memcpy(seqdata, getFragRecordSequence(rsp) + clrbgn, (flen+1)*sizeof(char));
      memcpy(qvdata,  getFragRecordQuality(rsp)  + clrbgn, (flen+1)*sizeof(char));
      seqdata[flen] = '\0';
      qvdata[flen] = '\0';
      if (reads[i].position.bgn > reads[i].position.end) {
        SequenceComplement(seqdata,qvdata);
      }
      accession = getFragRecordUID(rsp);
      ResetErrorStruct(errors);
      
      readptr= 0;
      delptr = 0;
      for ( j=left;j<right;j++) {
        assert (j < ma_len );
        column_call=*Getchar(ma->consensus,j);
        if ( delptr < reads[i].delta_length ) {
          if ( readptr != *(reads[i].delta + delptr)) {
            // non gap coverage for this fragment in this column
            // compare base at readptr[i] in frag sequence to column_call
            if ( seqdata[readptr] != column_call ) {
              // record the error
              frag_error.position = readptr;
              if (column_call == '-') { // insertion
                frag_error.type = CNS_INSERT; 
              } else {
                frag_error.type = CNS_SUBSTITUTE; 
              } 
              AppendErrorStruct(errors,&frag_error);
              column_mm[j]+=1;
            }
            readptr++;  column_cov[j]+=1;
          } else {
            // gap for this fragment in this column
            if ( '-' != column_call ) {
              // record the error
              frag_error.position = readptr;
              frag_error.type = CNS_DELETE;
              AppendErrorStruct(errors,&frag_error);
              column_mm[j]+=1;
              column_cov[j]+=1; //adding this so that intra-fragment gaps count as coverage
            }
            delptr++;
          }
        } else {
          // non gap coverage for this fragment in this column
          // compare base at readptr[i] in frag sequence to column_call
          if ( seqdata[readptr] != column_call ) {
            // record the error
            frag_error.position = readptr;
            if (column_call == '-') { // insertion
              frag_error.type = CNS_INSERT; 
            } else {
              frag_error.type = CNS_SUBSTITUTE; 
            } 
            AppendErrorStruct(errors,&frag_error);
            column_mm[j]+=1;
          }
          readptr++;  column_cov[j]+=1;
        }
      }
      fprintf(frag_stats,F_IID " " F_UID " %c %d %d",
              reads[i].ident,accession,
              reads[i].type,(int) clrbgn,(int) clrend);
      num_errors = 0;
      if ( GetNumErrorStructs(errors) > 75 ) {
        fprintf(frag_stats," misaligned fragment with %d mismatches\n",
                (int) GetNumErrorStructs(errors));
      } else {
        while ( GetErrorStruct(errors,num_errors) ) {
          frag_error = *GetErrorStruct(errors,num_errors);
          fprintf(frag_stats," %d %c",frag_error.position, frag_error.type);
          num_errors++;
        }
        fprintf(frag_stats,"\n"); 
      }
    }
    for (j=0;j<ma_len;j++) {
      fprintf(column_stats,"%d %d %d %c %d %d %d\n",ma->maID,j,
              column_cov[j],
              *Getchar(ma->consensus,j),
              (int) *Getchar(ma->quality,j) - '0',
              column_mm[j],ungapped);
      if (*Getchar(ma->consensus,j) != '-') ungapped++;
    }
  }
  fflush(column_stats);
  fflush(frag_stats);
  safe_free(column_cov);
  safe_free(column_mm);
  DeleteVA_ErrorStruct(errors);
  del_fragRecord(rsp);
}

void GenerateContigAlignmentStats(char *phase){
  GraphCGW_T *graph = ScaffoldGraph->ContigGraph;
  FILE *pcs = NULL;
  FILE *pfs = NULL;
  FILE *dcs = NULL;
  FILE *dfs = NULL;
  char buffer[256];
  GraphNodeIterator     nodes;
  ContigT		*ctg;
  sprintf(buffer,"stat/%s.CNS.placed_column_stats",phase);
  pcs = fopen(buffer,"w");
  sprintf(buffer,"stat/%s.CNS.placed_frag_stats",phase);
  pfs = fopen(buffer,"w");
  sprintf(buffer,"stat/%s.CNS.dregs_column_stats",phase);
  dcs = fopen(buffer,"w");
  sprintf(buffer,"stat/%s.CNS.dregs_frag_stats",phase);
  dfs = fopen(buffer,"w");
  assert(pfs && pcs && dfs && dcs);
  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  /* 1st get min and max values */
  while(NULL != (ctg = NextGraphNodeIterator(&nodes))){
    CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, ctg->scaffoldID);
    MultiAlignT *ma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ctg->id, graph->type == CI_GRAPH);
    if (scaffold && (scaffold->type == REAL_SCAFFOLD)) { // contig is placed
      CollectStats(ma, ScaffoldGraph->gkpStore, pcs, pfs, AS_READ_CLEAR_LATEST);
    } else {
      CollectStats(ma, ScaffoldGraph->gkpStore, dcs, dfs, AS_READ_CLEAR_LATEST);
    }
  }
  fclose(pcs);
  fclose(pfs);
  fclose(dcs);
  fclose(dfs);
}

void ComputeFragmentMembershipStats(void){
  int numFrags = GetNumCIFragTs(ScaffoldGraph->CIFrags);
  int i;
  int reads = 0;
  int uniqueReads = 0;
  int repeatReads = 0;
  int smallReads = 0;
  int bigReads = 0;
  CIFragT *frag;

  for(i = 0,  frag = GetCIFragT(ScaffoldGraph->CIFrags, 0); 
      i < numFrags; i++, 
	frag++){
    NodeCGW_T *ci = NULL;

    if (! AS_FA_READ(frag->type)) 
      continue;
    reads++;
    ci = GetGraphNode(ScaffoldGraph->CIGraph, frag->cid);
    if(ci->flags.bits.isUnique){
      uniqueReads++;
    }else{
      repeatReads++;
    }
    if(ci->bpLength.mean >= CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH){
      bigReads++;
    }else{
      smallReads++;
    }
  }
  fprintf(GlobalData->stderrc," There are %d reads of which %d are in unique unitigs and %d are in repeat unitigs\n",
	  reads, uniqueReads, repeatReads);
  fprintf(GlobalData->stderrc,"    %d reads are in unitigs of length < %d\n",
	  smallReads,CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH );
  fflush(GlobalData->stderrc);
}
