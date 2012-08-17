
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2012, J. Craig Venter Institute.
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

#ifndef GRAPH_EDGE_ITERATOR_H
#define GRAPH_EDGE_ITERATOR_H

static const char *rcsid_GRAPH_EDGE_ITERATOR_H = "$Id: GraphEdgeIterator.H,v 1.4 2012-08-17 05:25:28 brianwalenz Exp $";

//  includes are broken here
//#include "GraphCGW_T.h"
//#include "AS_CGW_dataTypes.h"
//#include "InputDataTypes_CGW.h"

//  We can iterate over three types of edges:
//    graph->SEdges          SEdgeT
//    graph->CIEdges         CIEdgeT
//    graph->ContigEdges     CIEdgeT
//
//  The ScaffoldGraphT has three edge lists:
//    VA_TYPE(CIEdgeT)        *CIEdges;
//    VA_TYPE(CIEdgeT)        *ContigEdges;
//    VA_TYPE(SEdgeT)         *SEdges;
//
//  Which edge list is used depends on which GraphCGW_T is passed in.  The ScaffoldGraphT
//  has three graphs:
//    GraphCGW_T             *CIGraph;
//    GraphCGW_T             *ContigGraph;
//    GraphCGW_T             *ScaffoldGraph;



class GraphEdgeIterator {
public:
  GraphEdgeIterator(GraphCGW_T *graph_,
                    CDS_CID_t   cid_,
                    int         end_,
                    int         edgeStatusSet_);

  ~GraphEdgeIterator() {
  };

  void         excludeContainment(void) {
    noContains = true;
  };

  EdgeCGW_T   *nextRaw(void);
  EdgeCGW_T   *nextMerged(void);
  EdgeCGW_T   *prevMerged(void)   { return(GetGraphEdge(graph, prevM)); };

  CDS_CID_t    prevMergedID(void) { return(prevM); };  //  Used when adding edges to the graph.
  CDS_CID_t    currMergedID(void) { return(currM); };  //  Quite ugly.  Same for prevMerged().

  void         reset(void) {
    prevM              = NULLINDEX;
    currM              = NULLINDEX;
    nextM              = GetGraphNode(graph, cid)->edgeHead;

    nextR              = NULLINDEX;
  };

private:
  void         sortEdges(void);

  GraphCGW_T  *graph;
  CDS_CID_t    cid;
  int          end;
  uint32       edgeStatusSet;

  bool         noContains;

  CDS_CID_t    prevM;
  CDS_CID_t    currM;
  CDS_CID_t    nextM;

  CDS_CID_t    prevR;
  CDS_CID_t    currR;
  CDS_CID_t    nextR;
};


bool  edgeCompare(EdgeCGW_T const *edge1, EdgeCGW_T const *edge2);



#endif  //  GRAPH_EDGE_ITERATOR_H
