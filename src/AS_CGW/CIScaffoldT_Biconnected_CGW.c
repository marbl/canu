
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
static char CM_ID[] = "$Id: CIScaffoldT_Biconnected_CGW.c,v 1.6 2006-10-08 08:47:39 brianwalenz Exp $";

//#define DEBUG 1
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
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "UnionFind_AS.h"
#include "UtilsREZ.h"
#include "ChiSquareTest_CGW.h"
#include "MultiAlignment_CNS.h"
#include "DataTypesREZ.h"
#include "CommonREZ.h"
#include "Stats_CGW.h"   // for collecting scaffold merging stats

/****************************************************************************/
/* IsScaffold2EdgeConnected is intended to identify the case where a single mate link
   is all that holds a scaffold together.  Basically we are looking for a 'bridge' edge,
   where the bridge is a single mate link.
   To identify this situation, we use a variant on the biconnected component
   implementation from LEDA.  This dfs-based algorithm groups the edges in a graph into disjoint
   components.  Any component containing a single edge is a bridge edge.
   Bridge edges are marked using a flag bit, and ignored by the code in CheckScaffoldConnectivityAndSplit.
   So, marking an edge as a bridge edge SHOULD disconnect the scaffold into two components.
*/


/* Create a parametrized stack of PtrT variables */

STACK_DEF(PtrT)

     /****** Stolen from LEDA _bicomponents.c
             For a detailed reference, see "Computer Algorithms" , Horowitz,Sahni, Rajesekaran pg 335.
     ********/
     static void bcc_dfs(ScaffoldGraphT *sgraph,
                         ContigT **contigs, ContigT *contig,
                         int32 *dfsnum, int32 *lowpt, int32 *father,
                         Stack_PtrT *stack, int32 *count1, int32 *count2,
                         int32 *numBridges);


     static Stack_PtrT *Stack = NULL;

     int IsScaffold2EdgeConnected(ScaffoldGraphT *graph, CIScaffoldT *scaffold){
       int numElements = scaffold->info.Scaffold.numElements;
       NodeCGW_T *contigs[numElements], *contig1;
       int32    dfsnum[numElements];
       int32    father[numElements];
       int32    lowpt[numElements];
       CIScaffoldTIterator Nodes;
       int count1 = 0;
       int count2 = 0;
       CDS_CID_t i = 0;
       int num_isolated = 0;
       int32 numBridges = 0;


       if(Stack == NULL){
         Stack = CreateStack_PtrT(numElements);
       }else{
         ResetStack_PtrT(Stack);
       }



       InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &Nodes);

       while((contig1 = NextCIScaffoldTIterator(&Nodes)) != NULL){
         assert(i < numElements);
         dfsnum[i] = NULLINDEX;
         father[i] = NULLINDEX;
         lowpt[i] = NULLINDEX;
         contig1->microhetScore = (float)(i) + 0.01; // this field is unused for contigs
#ifdef DEBUG
         fprintf(stderr,"* contig " F_CID " is " F_CID " (%g)\n",
                 i, contig1->id, contig1->microhetScore);
#endif
         contigs[i++] = contig1;
       }

       assert(i == numElements);

#ifdef DEBUG
       fprintf(stderr,
               "* IsScaffoldBiconnected on scaffold " F_CID " with %d elements\n",
               scaffold->id, numElements);
#endif

       for(i = 0; i < numElements; i++){
         NodeCGW_T *contig = contigs[i];

         if(dfsnum[i] == NULLINDEX){
           int isIsolated = 1;
           EdgeCGW_T *edge;
           GraphEdgeIterator Edges;

           // Iterate over raw edges
           InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, contig->id, ALL_END, ALL_TRUSTED_EDGES, GRAPH_EDGE_DEFAULT, &Edges);
      
           dfsnum[i] = ++count1;
#ifdef DEBUG
           fprintf(stderr,"* i = " F_CID "\n", i);
           fprintf(stderr,"* Contig " F_CID " dfsnum %d\n", contig->id, dfsnum[i]);
#endif
      
           // Does this node have any edges that are scaffold internal
           while (NULL != (edge = NextGraphEdgeIterator(&Edges))){
             NodeCGW_T *nodeA = GetGraphNode(ScaffoldGraph->ContigGraph, edge->idA);
             NodeCGW_T *nodeB = GetGraphNode(ScaffoldGraph->ContigGraph, edge->idB);

             // filter out edges that are not internal to this scaffold
             if(nodeA->scaffoldID != nodeB->scaffoldID ||
                isSingletonOverlapEdge(edge))
               continue;
	
             assert(nodeA->scaffoldID == scaffold->id);
             assert(nodeB->scaffoldID == scaffold->id);

             isIsolated = 0;
             break;
           }

           if(isIsolated){
#ifdef DEBUG
             fprintf(stderr,"* Node " F_CID " is isolated\n",contig->id);
#endif
             num_isolated++;
           }else{
#ifdef DEBUG
             fprintf(stderr,
                     "* Pushing on stack and calling dfs for node " F_CID "\n",
                     contig->id);
#endif
             PushStack_PtrT(Stack, contig);
             bcc_dfs(ScaffoldGraph, contigs, contig,  dfsnum, lowpt, father, Stack, &count1, &count2, &numBridges);
#ifdef DEBUG
             fprintf(stderr,"* Popping stack for contig " F_CID "\n", contig->id);
#endif
             PopStack_PtrT(Stack);
           }
         }
       }

       return (numBridges == 0);

     }


// bcc_dfs
//   This is the main dfs-based function for determining the bicomponents and labeling the bridges
//
static void bcc_dfs(ScaffoldGraphT *sgraph,
                    ContigT **contigs, ContigT *contig,
                    int32 *dfsnum, int32 *lowpt, int32 *father,
                    Stack_PtrT *stack, 
		    int32 *count1, int32 *count2, int32 *numBridges){
  CDS_CID_t contigNum = (CDS_CID_t)contig->microhetScore;
  GraphEdgeIterator Edges;
  EdgeCGW_T *edge;

  lowpt[contigNum] = dfsnum[contigNum];
#ifdef DEBUG
  fprintf(stderr, "contig " F_CID " contigNum %d lowpt[contigNum] = %d\n",
          contig->id,contigNum,lowpt[contigNum]);
#endif
  
  // Iterate over raw edges
  InitGraphEdgeIterator(sgraph->ContigGraph, contig->id, ALL_END, ALL_TRUSTED_EDGES, GRAPH_EDGE_DEFAULT, &Edges);

  while (NULL != (edge = NextGraphEdgeIterator(&Edges))){
    NodeCGW_T *otherContig = GetGraphNode(ScaffoldGraph->ContigGraph, (edge->idA == contig->id? edge->idB:edge->idA));
    CDS_CID_t otherNum = (CDS_CID_t)otherContig->microhetScore;

    // filter out edges that are not internal to this scaffold
    if(otherContig->scaffoldID != contig->scaffoldID ||
       isSingletonOverlapEdge(edge))
      continue;
    
#ifdef DEBUG
    fprintf(stderr," (" F_CID "," F_CID ",%c) wght %d\n",
	    edge->idA, edge->idB, edge->orient, edge->edgesContributing);
#endif
    
    if(dfsnum[otherNum] == NULLINDEX){ // haven't seen it yet
      dfsnum[otherNum] = ++(*count1);
      PushStack_PtrT(stack, otherContig);
      father[otherNum] = contigNum;

#ifdef DEBUG
      fprintf(stderr,"Pushing ctg " F_CID " contigNum %d dfsnum %d father %d and calling dfs\n",
              otherContig->id,otherNum,dfsnum[otherNum],father[otherNum]);
#endif
      
      bcc_dfs(sgraph, contigs, otherContig,  dfsnum, lowpt, father,stack, count1, count2, numBridges);
      lowpt[contigNum] = min(lowpt[contigNum], lowpt[otherNum]);
    }else{
      lowpt[contigNum] = min(lowpt[contigNum], dfsnum[otherNum]);
    }
#ifdef DEBUG
    fprintf(stderr,"lowpt[contigNum] = %d\n",lowpt[contigNum]);
#endif
  } 

  if(father[contigNum] != NULLINDEX &&
     (lowpt[contigNum] == dfsnum[father[contigNum]])){
    CDS_CID_t wNum;
    int cnt = 0;
    EdgeCGW_T *lastEdge = NULL;

    do{
      GraphEdgeIterator Edges2;
      NodeCGW_T *w = PopStack_PtrT(stack);
      EdgeCGW_T *edge2;

      wNum = (CDS_CID_t)w->microhetScore;
      InitGraphEdgeIterator(sgraph->ContigGraph, w->id, ALL_END, ALL_TRUSTED_EDGES, GRAPH_EDGE_DEFAULT, &Edges2);

      while (NULL != (edge2 = NextGraphEdgeIterator(&Edges2))){
        NodeCGW_T *otherContig = GetGraphNode(ScaffoldGraph->ContigGraph, (edge2->idA == w->id? edge2->idB:edge2->idA));
        CDS_CID_t otherNum = (CDS_CID_t)otherContig->microhetScore;

        if(otherContig->scaffoldID != w->scaffoldID ||
           isSingletonOverlapEdge(edge2))
          continue;
	  
#ifdef DEBUG
        fprintf(stderr,
                " w:" F_CID " cnt:%d (" F_CID "," F_CID ",%c) wght %d\n",
                w->id, cnt, edge2->idA, edge2->idB,
                edge2->orient, edge2->edgesContributing);
#endif
        if(dfsnum[wNum] > dfsnum[otherNum]){
          edge2->flags.bits.isBridge = FALSE;
          lastEdge = edge2;
          cnt++;
          w->setID = otherContig->setID = (*count2);
        }
      }
	
    }while(wNum != contigNum);

    if(cnt == 1){
      int numLinks = lastEdge->edgesContributing - (isOverlapEdge(lastEdge)?1:0);
      if(numLinks == 1){
        (*numBridges)++;
        lastEdge->flags.bits.isBridge = TRUE;
#ifdef DEBUG
        fprintf(stderr,"* Edge (" F_CID "," F_CID ",%c) is a bridge\n",
                lastEdge->idA, lastEdge->idB, lastEdge->orient);
#endif
        PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, "Bridge ", lastEdge, lastEdge->idA);

      }
    }
      
    (*count2)++;
  }

}
