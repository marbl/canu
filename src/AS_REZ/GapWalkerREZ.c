
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
/**********************************************************************

        Module:  GapWalkerREZ.c

   Description:  GapWalker is a collection of routines to walk in
                 subgraphs of the chunk graph. It features
		 - depth first searching
		 - greedy path finding
		 - Dijkstra single source shortest path

		 log info is sent to <inputFile>.gwlog

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)
                 K. Reinert
		 M. Flanigan

       Written:  20 May 99

 **********************************************************************/


static char fileID[] = "$Id: GapWalkerREZ.c,v 1.4 2005-03-22 19:49:22 jason_miller Exp $";


#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
//
// AS_CGW
//
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ChiSquareTest_CGW.h"

//
// AS_REZ
//
#include "DataTypesREZ.h"
#include "UtilsREZ.h"
#include "CommonREZ.h"
#include "ConsistencyChecksREZ.h"
#include "GapWalkerREZ.h"
#include "BccREZ.h"
#include "SubgraphREZ.h"
#include "GWDriversREZ.h"
#include "FbacREZ.h"

//
// vars
//

int
  start_index[4] = {0, 2, 0, 0}, // remember that A_END == 1, B_END == 2, ALL_END = 3 (see AS_CTW_dataTypes.h)
  end_index[4] = {0, NUM_ORIENTATIONS, 2, NUM_ORIENTATIONS},
  new_end[NUM_ORIENTATIONS] = {B_END, A_END, B_END, A_END},
  orient[NUM_ORIENTATIONS] = {OR2NUM_AB_AB, OR2NUM_AB_BA, OR2NUM_BA_AB, OR2NUM_BA_BA};

//
// externs
//
#if CREATE_CAM_FILE > 0
extern FILE
  * walker_cam_file;
#endif

extern char
  * GW_Filename_Prefix;


//
// filters
//


int Stop_At_Unique(ChunkInstanceT * to, CDS_CID_t tid, chunk_subgraph * s) {
  //
  // termination condition used in Find_Greedy_Path
  //
  return Is_Unique(to);
}



int Stop_At_The_Other_End(ChunkInstanceT * to, CDS_CID_t tid, chunk_subgraph * s) { 
  //
  // termination condition used in Find_Greedy_Path
  //
  return (tid == to->id);
}



int Stop_At_Path_Bit_Set(ChunkInstanceT * to, CDS_CID_t tid, chunk_subgraph * s) {
  //
  // termination condition used in Find_Greedy_Path
  //
  assert(s->table[to->id] != NULL);
  return s->table[to->id]->path_bit;
}


//
// edge quality
//

float Edge_Quality_GW(CIEdgeT * edge, CDS_CID_t cid) {
  //
  // Return a real value that indicates the reliability of the edge
  // mate (* edge) going out from <cid>
  //
  // we need a function to be MINimised: that is better
  // edges have lower score
  //
  float
    val;
  
  val = 3.0 * edge->edgesContributing;
 
  if (isProbablyBogusEdge(edge)) 
    return 50.0;

  if (edge->flags.bits.isPossibleChimera)
    val -= 0.2;
  if (edge->flags.bits.hasContributingOverlap)
    val -= 1.0;
  if (edge -> flags.bits.hasTandemOverlap)
    val -= 3.0;
  if (edge->flags.bits.hasRepeatOverlap)
    val -= 2.0;
  if (edge->flags.bits.hasGuide)
    val -= 0.5;

# if 0
  {
    CDS_CID_t
      other_cid;
    ChunkInstanceT 
      * other_chunk;
    
    //
    // get the other end
    //
    if (cid == edge->idA)
      other_cid = edge->idB;
    else
      other_cid = edge->idA;
    
    other_chunk = GetGraphNode(ScaffoldGraph->RezGraph, other_cid);
    assert(other_chunk != NULL);
    
    if (Is_Unique(other_chunk)) {
     val += 25.0;
     //fprintf(stderr,"> " F_CID " to unique " F_CID "\n", cid, other_cid);
    }
  }
# endif

  return max(0.0, 50.0 - val);
}




float No_Quality(CIEdgeT * edge, CDS_CID_t cid) {
  //
  // self-explanatory
  //
  return 1.0;
}

float Bayesian_Quality(CIEdgeT *edge, CDS_CID_t cid)
{
  ChunkOverlapCheckT olap;
  float quality;
  GraphCGW_T *graph = ScaffoldGraph->RezGraph;

  // Is there an overlap ?? If yes, return its quality
  if( TRUE == LookupQualityOverlap(graph,edge,
				   edge->orient,&olap,BAYESIAN,
				   &quality,stderr) )
    {
      // now we adapt the edge->distance.mean to our computed overlap
#if DEBUG_GAP_WALKER > 3
      PrintGraphEdge(stderr,ScaffoldGraph->RezGraph,"before adapting ", edge, edge->idA);
#endif
      edge->flags.bits.MeanChangedByWalking = TRUE;
      edge->minDistance   =  edge->distance.mean;
      edge->distance.mean = -olap.overlap;
      /* save the value of the distance mean into the field minDistance for retrival
	 in the cgw output */

#if DEBUG_GAP_WALKER > 3
      PrintGraphEdge(stderr,ScaffoldGraph->RezGraph,"after adapting ", edge, edge->idA);
#endif
      return(quality);
    }
  else
    {
      return 1.0;
    }
}

float Bac_Walking_Quality(CIEdgeT *edge, CDS_CID_t cid)
{
  ChunkOverlapSpecT spec;
  ChunkOverlapperT *chunkOverlapper = ScaffoldGraph->RezGraph->overlapper;
  ChunkOverlapCheckT *lookup;
  int isCanonical;

  isCanonical = InitCanonicalOverlapSpec(edge->idA, edge->idB, edge->orient, &spec);
  lookup = LookupCanonicalOverlap(chunkOverlapper, &spec);
  if(!lookup)  // We didn't find anything
    return FALSE;

  // There was something in the table but the overlap length is 0
  if( lookup->overlap == 0 )
    return FALSE;

  // There is an overlap, return its quality
  return(lookup->quality);
}





static int Compare_Edge_Quality(const void *e1,
				const void *e2) {
  //
  // The function that compares two edge_quality (for the qsort())
  //
  edge_quality *eq1 = (edge_quality *)e1;
  edge_quality *eq2 = (edge_quality *)e2;

  assert((eq1 != NULL) && (eq2 != NULL));

  if (eq1->quality > eq2->quality)
    return 1;
  else if (eq1->quality < eq2->quality)
    return -1;
  else
    return 0;
}


//
// Print_Cam_Subgraph() prints all the chunks not used by any path or
// not visited and the uniques in a .cam file (the other are
// supposedly printed by Visit_Subgraph())
//

static void Print_Cam_Subgraph(chunk_subgraph * s,
			       FILE * wfile) {
  //
  // print all the chunks not used by any path or not visited
  // and the uniques in a .cam file (the others are printed
  // by Visit_Subgraph())
  //
  CDS_CID_t
    a_end,
    b_end,
    i,
    id;
  int
    color;
  ChunkInstanceT 
    * chunk;
  char
    orient[STR_LEN],
    unique[STR_LEN],
    visited[STR_LEN],
    path[STR_LEN];

  assert (wfile != NULL);

  for (i = 0; i < s->size; i++) {
    id = s->node[i].cid;
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, id);
    assert(chunk != NULL);
    if (chunk->aEndCoord <= chunk->bEndCoord) {
      a_end = chunk->aEndCoord;
      b_end = chunk->bEndCoord;
      strcpy(orient, "direct");
    } else {
      a_end = chunk->bEndCoord;
      b_end = chunk->aEndCoord;
      strcpy(orient, "reverse");
    }

    if (Is_Unique(chunk)) {
      strcpy(unique, "unique");
      color = MY_UNIQUE_COLOUR;
    } else {
      strcpy(unique, "not unique");
      color = 0;
    }

    if (s->node[i].visited) {
      strcpy(visited, "visited");
      color = 0;
    } else {
      strcpy(visited, "not visited");
      color = NOT_VISITED_COLOUR;
    }

    if (s->node[i].path_id) {
      strcpy(path, "path");
      color = 0;
    } else {
      strcpy(path, "not path");
      color = NOT_PATH_COLOUR;
    }

    if (Is_Unique(chunk) ||
	(! s->node[i].path_id) ||
	(! s->node[i].visited))
      fprintf(wfile,
	      "%d: %d A%d %d # (%s - %s - %s) chunk %d\n",
	      id,
	      a_end,
	      color,
	      b_end,
	      unique,
	      visited,
	      path,
	      id);
  }
}


//
// Print_Path() prints all the chunk id that are on some path
//
static void Print_Path (chunk_subgraph * s) {
  //
  // prints all the chunk id that are on some path
  //
  int
    i,
    cnt = 0;

  assert(s != NULL);
  for (i = 0; i < s->size; i++) {
    if (s->node[i].path_id) {
#     if DEBUG_GAP_WALKER > 0
      fprintf(GlobalData->gwlogfp,
	      "cid %5d, path_id %4d, path_parent %4d, (%5.2f stdDev %5.2f)\n",
	      s->node[i].cid,
	      s->node[i].path_id,
	      s->node[i].path_parent,
	      s->node[i].distance.mean,
	      sqrt(s->node[i].distance.variance));
#     endif
      cnt++;
    }
  }
# if DEBUG_GAP_WALKER > 0
  fprintf(GlobalData->gwlogfp,"\n Total nodes on some path: %d\n", cnt);
# endif
}



//
// walking code
//


void Visit_Subgraph(chunk_subgraph * s,
		    CDS_CID_t begin_cid,
		    CDS_CID_t end_cid,
		    int end,
		    LengthT * max,
		    float (* quality)(CIEdgeT *, CDS_CID_t),
		    int (* terminate)(ChunkInstanceT *, CDS_CID_t, chunk_subgraph *)) {
  //
  // it finds and prints all the paths from <begin_cid> to <end_cid>
  // into a "*.<begin_cid>.<end_cid>gw.cam" file
  // it uses Find_Greedy_Path() as a subroutine
  //
  // NOTE: this version is directional, that is if it enters a chunk
  // from the A end it goes out from the B end and viceversa
  //
# if CREATE_CAM_FILE > 0
  char
    * Colour[PATH_NUM_COLOURS] = {
      "CFFFF00 T4 S # start/end",
      "C888888 T2 S # not_visited",
      "CBBBBBB T2 S # not_path",
      "C0000FF T2 S # path_1",
      "C00FF00 T2 S # path_2",
      "CFF0000 T2 S # path_3",
      "C00FFFF T2 S # path_4", 
      "CFFFF00 T2 S # path_5",
      "CFF00FF T2 S # path_6",
      "CFF00AA T2 S # path_7",
      "CAA00FF T2 S # path_8",
      "C00AAFF T2 S # path_9"};
  char
    filename[256];
  FILE
    * walker_cam_file;
# endif
  LengthT
    dist;
  int
    i,
    c = FIRST_PATH_COLOUR;
  DIR *dbDir;
  mode_t mode = S_IRWXU | S_IRWXG | S_IROTH;
  long currentCalls=0;
  float tooShort = 0.0;
  float tooLong = 0.0;

#include <sys/types.h>
  assert(s != NULL);
  assert(quality != NULL);

  fprintf(stderr,"** open cam dir\n");
  dbDir = opendir("cam");
  
  if(dbDir == NULL)
    {
      fprintf(stderr,"Could not open cam dir\n");
      if(mkdir("cam",mode))
	{
	  exit(1);
	}
    }
  //
  // if the start/end chunk are not in <s>, add them
  //
  Add_Node(s, begin_cid, Is_Not_Bogus);
  Add_Node(s, end_cid, Is_Not_Bogus);

# if CREATE_CAM_FILE > 0
  //
  // open the file where the graphical path information
  // will be stored
  //
  sprintf(filename, "./cam/%s.%d.%d.gw.cam",
	  GW_Filename_Prefix,
	  begin_cid, end_cid);



  walker_cam_file = file_open (filename, "w");
  assert(walker_cam_file != NULL);

  //
  // spit out the colors
  //
  for (i = 0; i < PATH_NUM_COLOURS; i ++)
    fprintf(walker_cam_file, "%d: %s\n",
	    i,
	    Colour[i]);
# endif

  //
  // clear path/visited bits
  //
  Clear_All_Path_Bit(s);
  Clear_All_Visited_Bit(s);

  //
  // initialize the distance
  //
  dist.mean = 0.0;
  dist.variance = 0.0;

  //
  // initialize the distances
  //
  for (i = 0; i < s->size; i++)
    s->node[i].d = 0.0;

  //
  // search for the paths and print them out with different colours
  //
  while (Find_Greedy_Path(s, begin_cid, begin_cid, end_cid, max, dist,
			  0, end, quality, NO_QUALITY_THRESH,  terminate, c, FALSE, &currentCalls,MAXCALLS,&tooShort,&tooLong)) { 
#   if CREATE_CAM_FILE > 0
    ChunkInstanceT
      * to_chunk;
#endif

#   if DEBUG_GAP_WALKER > 2
    fprintf(GlobalData->gwlogfp,"\nNodes on the path (ordered by chunk id):\n");
    Print_Path(s);
#   endif
    Clear_All_Path_Bit(s);          // clear the path bit 
#   if CREATE_CAM_FILE > 0
    //
    // print the path in the camfile (backwards)
    //
    i = end_cid;
    while (i != begin_cid) {
      assert(s->table[i] != NULL);
      if (s->table[i]->best_edge == NULL)
	break;
      fprintf (walker_cam_file, "LNK: %d %d A%d # edge quality %5.2f contrib %d\n",
	       s->table[i]->path_parent,
	       i,
	       c,
	       quality(s->table[i]->best_edge, begin_cid),
	       s->table[i]->best_edge->edgesContributing);
      to_chunk = GetGraphNode(ScaffoldGraph->RezGraph, i);
      if (! Is_Unique(to_chunk)) 
	if (to_chunk->aEndCoord <= to_chunk->bEndCoord)
	  fprintf (walker_cam_file,
		   "%d: %d A%d %d # chunk %d, mean distance %5.2f, var %5.2f, cumul quality %5.2f\n",
		   i,
		   to_chunk->aEndCoord,
		   c,
		   to_chunk->bEndCoord,
		   i,
		   s->table[i]->distance.mean,
		   s->table[i]->distance.variance,
		   s->table[i]->d);
	else
	  fprintf (walker_cam_file,
		   "%d: %d A%d %d # reverse chunk %d, mean distance %5.2f, var %5.2f, cumul quality %5.2f\n",
		   i,
		   to_chunk->bEndCoord,
		   c,
		   to_chunk->aEndCoord,
		   i,
		   s->table[i]->distance.mean,
		   s->table[i]->distance.variance,
		   s->table[i]->d);

      //
      // follow the parent
      //
      i = s->table[i]->path_parent;
    }

    //
    // choose the next color. if used all them
    // start again from the beginning
    //
    c++;
    if (c == PATH_NUM_COLOURS)
      c = FIRST_PATH_COLOUR; 
#   endif
  }

  //
  // restore the status of the end chunk
  //
  Set_Visited_Bit(s, end_cid);
  Set_Visited_Bit(s, begin_cid);
# if DEBUG_GAP_WALKER > 2
  Print_Subgraph(s);
# endif

# if CREATE_CAM_FILE > 0
  //
  // print all the other chunks
  //
  Print_Cam_Subgraph(s, walker_cam_file);

  //
  // good boys close the files
  //
  fclose (walker_cam_file);
# endif
}




// FindGapLength provides the functionality of Compute_Gap_Length but 
// is much more straightforward
LengthT FindGapLength( ChunkInstanceT * lchunk,
					   ChunkInstanceT * rchunk,
					   int verbose)
{
  LengthT gapSize;
  float lchunkMaxOffset, lchunkMaxVariance;
  float rchunkMinOffset, rchunkMinVariance;
  
  if ( lchunk->offsetAEnd.mean < lchunk->offsetBEnd.mean)
  {
	lchunkMaxOffset = lchunk->offsetBEnd.mean;
	lchunkMaxVariance = lchunk->offsetBEnd.variance;
  }
  else
  {
	lchunkMaxOffset = lchunk->offsetAEnd.mean;
	lchunkMaxVariance = lchunk->offsetAEnd.variance;
  }

  if ( rchunk->offsetAEnd.mean < rchunk->offsetBEnd.mean)
  {
	rchunkMinOffset = rchunk->offsetAEnd.mean;
	rchunkMinVariance = rchunk->offsetAEnd.variance;
  }
  else
  {
	rchunkMinOffset = rchunk->offsetBEnd.mean;
	rchunkMinVariance = rchunk->offsetBEnd.variance;
  }

  gapSize.mean = rchunkMinOffset - lchunkMaxOffset;
  gapSize.variance = rchunkMinVariance - lchunkMaxVariance;

  return gapSize;
}	

LengthT Compute_Gap_Length(ChunkInstanceT * lchunk,
			   ChunkInstanceT * rchunk,
			   int verbose) {
  //
  // find the gap size between two chunks: handles the case of intra
  // scaffolds gaps
  //
  NodeOrient
    lorient = GetNodeOrient(lchunk),
    rorient = GetNodeOrient(rchunk);
  LengthT
    gap = {0.0, 0.0};
  
  assert(lchunk != NULL);
  assert(rchunk != NULL);
  assert(lchunk->scaffoldID == lchunk->scaffoldID);
  if (min(lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean) >
      min(rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean)) { 
    fprintf(stderr,
	    "                       Warning: chunks %d and %d are out of order (scaffold %d)\n",
	    lchunk->id,
	    rchunk->id,
	    lchunk->scaffoldID);
  }
  
# if DEBUG_GAP_WALKER > -1
  if( verbose )
    {
      MultiAlignT *lma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, lchunk->id, ScaffoldGraph->RezGraph->type == CI_GRAPH);
      MultiAlignT *rma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, rchunk->id, ScaffoldGraph->RezGraph->type == CI_GRAPH);
      //      MultiAlignT *lma = GetMultiAlignInStore(ScaffoldGraph->RezGraph->maStore,lchunk->id);
      //      MultiAlignT *rma = GetMultiAlignInStore(ScaffoldGraph->RezGraph->maStore,rchunk->id);
      
      VA_TYPE(char) *consensus = CreateVA_char(2048);
      VA_TYPE(char) *quality   = CreateVA_char(2048);
					    
      long ullength = GetMultiAlignUngappedLength(lma);
      long urlength = GetMultiAlignUngappedLength(rma);

	fprintf(GlobalData->gwlogfp,"\n*--------------------------------------------------*\n means: lchunk [%5.2f, %5.2f] o:%s - rchunk [%5.2f, %5.2f] o:%s\n", lchunk->offsetAEnd.mean, lchunk->offsetBEnd.mean,(lorient == A_B) ? "A_B" : "B_A",rchunk->offsetAEnd.mean, rchunk->offsetBEnd.mean,(rorient == A_B) ? "A_B" : "B_A");
      fprintf(GlobalData->gwlogfp,
	      "       vars: lchunk [%5.2f, %5.2f] - rchunk [%5.2f, %5.2f]\n",
	      lchunk->offsetAEnd.variance, lchunk->offsetBEnd.variance,
	      rchunk->offsetAEnd.variance, rchunk->offsetBEnd.variance);
      
      {
	GraphEdgeIterator edges;
	CIEdgeT *edge;
	SeqInterval interval;
	IntUnitigPos *pos;
	if( lorient == A_B )
	  {	    
	    pos = GetBendUnitigPos(lma);
	    fprintf(GlobalData->gwlogfp,"\nRightmost Unitig %d with position (%d,%d)\n",pos->ident,pos->position.bgn,pos->position.end);
	    fprintf(GlobalData->gwlogfp,"\nLast 1000 bps of consensus and quality of the B end of lchunk %d\n",lchunk->id);
	    interval.bgn = max(0,ullength - 1000);
	    interval.end = ullength;
	    GetMultiAlignUngappedConsensusFromInterval(lma,interval,consensus,quality);
	    fprintf(GlobalData->gwlogfp,"Consensus \n%s\n",Getchar(consensus,0));
	    fprintf(GlobalData->gwlogfp,"Quality   \n%s\n",Getchar(quality,0));


	    fprintf(GlobalData->gwlogfp,"\nEdges going out from the B end of lchunk %d\n",lchunk->id);

	    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, lchunk->id, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);

	    while(( edge = NextGraphEdgeIterator(&edges) ) != NULL)
	      PrintGraphEdge(GlobalData->gwlogfp,ScaffoldGraph->RezGraph,"relativ from idA ", edge, edge->idA);
	  }
	else
	  {
	    pos = GetAendUnitigPos(lma);
	    fprintf(GlobalData->gwlogfp,"\nLeftmost Unitig %d with position (%d,%d)\n",pos->ident,pos->position.bgn,pos->position.end);
	    fprintf(GlobalData->gwlogfp,"\nLast 1000 bps of consensus and quality of the A end of lchunk %d\n",lchunk->id);
	    interval.bgn = 0;
	    interval.end = min(ullength,1000);
	    GetMultiAlignUngappedConsensusFromInterval(lma,interval,consensus,quality);
	    fprintf(GlobalData->gwlogfp,"Consensus \n%s\n",Getchar(consensus,0));
	    fprintf(GlobalData->gwlogfp,"Quality   \n%s\n",Getchar(quality,0));

	    fprintf(GlobalData->gwlogfp,"\nEdges going out from the A end of lchunk %d\n",lchunk->id);

	    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, lchunk->id, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);

	    while(( edge = NextGraphEdgeIterator(&edges) ) != NULL)
	      PrintGraphEdge(GlobalData->gwlogfp,ScaffoldGraph->RezGraph,"relativ from idA ", edge, edge->idA);
	  }
	if( rorient == A_B )
	  {
	    pos = GetAendUnitigPos(rma);
	    fprintf(GlobalData->gwlogfp,"\nLeftmost Unitig %d with position (%d,%d)\n",pos->ident,pos->position.bgn,pos->position.end);
	    fprintf(GlobalData->gwlogfp,"\nLast 1000 bps of consensus and quality of the A end of lchunk %d\n",lchunk->id);   
	    fprintf(GlobalData->gwlogfp,"\nLast 1000 bps of consensus and quality of the A end of rchunk %d\n",rchunk->id);
	    interval.bgn = 0;
	    interval.end = min(urlength,1000);
	    GetMultiAlignUngappedConsensusFromInterval(rma,interval,consensus,quality);
	    fprintf(GlobalData->gwlogfp,"Consensus \n%s\n",Getchar(consensus,0));
	    fprintf(GlobalData->gwlogfp,"Quality   \n%s\n",Getchar(quality,0));

	    fprintf(GlobalData->gwlogfp,"\nEdges going out from the A end of rchunk %d\n",rchunk->id);

	    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, rchunk->id, A_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);

	    while(( edge = NextGraphEdgeIterator(&edges) ) != NULL)
	      PrintGraphEdge(GlobalData->gwlogfp,ScaffoldGraph->RezGraph,"relativ from idA ", edge, edge->idA);
	  }
	else
	  {
	    pos = GetBendUnitigPos(rma);
	    fprintf(GlobalData->gwlogfp,"\nRightmost Unitig %d with position (%d,%d)\n",pos->ident,pos->position.bgn,pos->position.end);
	    fprintf(GlobalData->gwlogfp,"\nLast 1000 bps of consensus and quality of the B end of rchunk %d\n",rchunk->id);
	    interval.bgn = max(0,urlength-1000);
	    interval.end = urlength;
	    GetMultiAlignUngappedConsensusFromInterval(rma,interval,consensus,quality);
	    fprintf(GlobalData->gwlogfp,"Consensus \n%s\n",Getchar(consensus,0));
	    fprintf(GlobalData->gwlogfp,"Quality   \n%s\n",Getchar(quality,0));

	    fprintf(GlobalData->gwlogfp,"\nEdges going out from the B end of rchunk %d\n",rchunk->id);

	    InitGraphEdgeIterator(ScaffoldGraph->RezGraph, rchunk->id, B_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
	    while(( edge = NextGraphEdgeIterator(&edges) ) != NULL)
	      PrintGraphEdge(GlobalData->gwlogfp,ScaffoldGraph->RezGraph,"relativ from idA ", edge, edge->idA);
	  }
      }
      DeleteVA_char(consensus);
      DeleteVA_char(quality);
    }
# endif
  
  switch (lorient) {
  case A_B:
    switch (rorient) {
    case A_B: // ---l--->  ----r--->
      gap.mean = rchunk->offsetAEnd.mean - lchunk->offsetBEnd.mean;
      gap.variance = rchunk->offsetAEnd.variance - lchunk->offsetBEnd.variance;
      break;
    case B_A: // ---l--->  <----r---
      gap.mean = rchunk->offsetBEnd.mean - lchunk->offsetBEnd.mean;
      gap.variance = rchunk->offsetBEnd.variance - lchunk->offsetBEnd.variance;
      break;
    default:
      assert(0);
      break;
    }
    break;
  case B_A:
    switch (rorient) {
    case A_B: //  <---l---  ----r--->
      gap.mean = rchunk->offsetAEnd.mean - lchunk->offsetAEnd.mean;
      gap.variance = rchunk->offsetAEnd.variance - lchunk->offsetAEnd.variance;
      break;
    case B_A: // <---l---  <----r---
      gap.mean = rchunk->offsetBEnd.mean - lchunk->offsetAEnd.mean;
      gap.variance = rchunk->offsetBEnd.variance - lchunk->offsetAEnd.variance;
      break;
    default:
      assert(0);
      break;
    }
    break;
  default:
    assert(0);
    break;
  }

# if DEBUG_GAP_WALKER > 1
  if( verbose)
    fprintf(GlobalData->gwlogfp,
	    "        gap:  mean %5.2f, variance %5.2f\n\n",
	    gap.mean, gap.variance);
# endif

# if DEBUG_GAP_WALKER > 0
  if( verbose )
    if (gap.variance < 0.0) {
      fprintf(stderr,
	      "* warning: the variance of the gap between chunk %d and %d is negative (%.3f).\n",
	      lchunk->id,
	      rchunk->id,
	      gap . variance);
      fflush(stderr);
	  DumpContig(stderr, ScaffoldGraph, lchunk, FALSE);
	  DumpContig(stderr, ScaffoldGraph, rchunk, FALSE);
      assert(gap.variance >= 0.0);
    }
# endif
  return gap;
}



chunk_subgraph * Intra_Scaffold_Gap_Walker(chunk_subgraph * full,
					   CDS_CID_t sid,
					   CIScaffoldTIterator * CIs,
					   float (* quality)(CIEdgeT *, CDS_CID_t),
					   float qualityThresh, int *hops, 
					   long* calls,float *tooShort,
					   float *tooLong) {
  //
  // Intra_Scaffold_Gap_walker() takes a chunk_subgraph of the graph we
  // want to analyze, a scaff_id, and two adjacents chunks in the
  // scaffold (CIs.curr and CIs.next), then:
  //
  // * it estimates the gap length between the two chunks
  //
  // * it tryes to isolates the chunks in the gap calling Build_Subgraph_Path()
  //
  // * it returns the subgraph
  //
  // * it checks if the chunks created by the previous step actually
  //   belongs to the gap (print debug info/a .dot file of that subgraph)
  //
  CDS_CID_t
    left_cid = CIs->curr,
    right_cid = CIs->next;
  ChunkInstanceT
    * lchunk = GetGraphNode(ScaffoldGraph->RezGraph, left_cid),
    * rchunk = GetGraphNode(ScaffoldGraph->RezGraph, right_cid);
  CIScaffoldT
    * scaff = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, sid);
  LengthT
    gap = {0.0, 0.0};
  chunk_subgraph
    * f0;

  assert(quality != NULL);
  assert(scaff != NULL);
  assert(! isDeadCIScaffoldT(scaff));
  assert(scaff->type == REAL_SCAFFOLD);
  assert(lchunk != NULL);
  assert(rchunk != NULL);
  assert(full != NULL);

  //
  // if the scaffold has one element (or less?)
  // return (this is already checked in the calling function
  // but does not any harm here)
  //
  if (scaff->info.Scaffold.numElements < 2)
    {
# if DEBUG_GAP_WALKER > 0
      fprintf(GlobalData->gwlogfp,"Scaffold has less than 2 elements\n");
#endif
      return NULL;
    }
  //
  // compute the gap size from the offset of the chunks in the
  // scaffold
  //
  // gap = Compute_Gap_Length(lchunk, rchunk, FALSE);
  gap = FindGapLength(lchunk, rchunk, FALSE);

  /*** mjf ***/ 
#if DEBUG_GAP_WALKER > 0
  fprintf( GlobalData->gwlogfp, 
			 "\nlooking for an intra-scaffold walk from %d ---> |%.2f| ---> %d\n", 
			 left_cid, gap.mean, right_cid);
#endif

  //
  // if one of the two chunks is not in the "full" graph, ie, the graph we passed to this function
  // return
  //
  if ((full->table[left_cid] == NULL) ||
      (full->table[right_cid] == NULL))
    {
      fprintf( GlobalData->gwlogfp,
               "\nNULL: full->table[%d] = %p, full->table[%d] = %p\n", 
               left_cid, full->table[left_cid],
               right_cid, full->table[right_cid]);
      return NULL;
    }

  //
  // if the gap is big enough to be reasonably "walked", then we
  // build the subgraph of all the chunks reachable from the two
  // flanking chunks
  //
  f0 = Build_Subgraph_Path(full, CIs, &gap, quality, qualityThresh,
						   Is_Not_Bogus_Not_Contained_Is_Overlap, TRUE, hops, calls, tooShort, tooLong);
  

  if (f0 == NULL) {

#   if DUMP_PATHS_NOT_CROSSED > 0
    //
    // let's see why we where not able to cross ...
    //
    Build_Subgraph_Path(full, CIs, &gap, quality,qualityThresh,
    			// Is_Not_Bogus_Not_Guide_Not_Pure_Overlap, TRUE);
			/*** mjf ***/ // OVERLAPWALKING
			Is_Not_Bogus_Not_Guide_Is_Overlap_Not_Tandem, TRUE);
    Clear_All_Path_Bit(full);
    Print_Subgraph_Cam(full, full, left_cid, right_cid, TRUE);
#   endif

    return NULL;
  }
  
# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
	  "    Subgraph of size %d on the gap of size %5.2f between %d and %d\n",
	  f0->size,
	  gap.mean,
	  left_cid,
	  right_cid);
  fflush(GlobalData->gwlogfp);
# endif
  
# if 1
  //
  // check if I am doing ok
  //
  {
    static int
      total = 0,
      bad = 0,
      bad_due_to_rr = 0,
      bad_due_to_ru = 0,
      bad_due_to_ur = 0,
      bad_due_to_other_combinations = 0;
    static double
      total_gap = 0.0;
    int
      i,
      uu,
      ru,
      ur,
      rr,
      cid,
      uniques,
      low,
      mistake = FALSE,
      high;
    ChunkInstanceT
      * chunk;
    
    assert(lchunk != NULL);
    assert(rchunk != NULL);
    low = min(min(lchunk->aEndCoord, lchunk->bEndCoord),
	      min(rchunk->aEndCoord, rchunk->bEndCoord));
    high = max(max(lchunk->aEndCoord, lchunk->bEndCoord),
	       max(rchunk->aEndCoord, rchunk->bEndCoord));
    uu = 0;
    ru = 0;
    rr = 0;
    ur = 0;
    uniques = 0;
    for (i = 0; i < f0->size; i++) {
      cid = f0->node[i].cid;
      chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
      assert(chunk != NULL);
      //
      // check only non-uniques
      //
      if (Is_Unique(chunk)) {
	uniques++;
	continue;
      }

      switch (chunk->flags.bits.cgbType) {
      case UU_CGBTYPE: uu++;
	break;
      case RR_CGBTYPE: rr++;
	break;
      case UR_CGBTYPE: ur++;
	break;
      case RU_CGBTYPE: ru++;
	break;
      default:
	assert(0);
	break;
      }

      if (((chunk->aEndCoord >= low) && (chunk->aEndCoord <= high)) ||
	  ((chunk->bEndCoord >= low) && (chunk->bEndCoord <= high)))   
	continue;
      else {
#       if DEBUG_GAP_WALKER > 1
	fprintf(GlobalData->gwlogfp,
		"    The chunk %d has been erroneusly selected\n    range [%d,%d] (flanking chunks %d, %d)\n",
		cid,
		low,
		high,
		left_cid,
		right_cid);
#       endif
	mistake = TRUE;
      }
    }
    if (mistake) {
      bad++;
      if ((rr > 0) && (ur == 0) && (ru == 0))
	bad_due_to_rr++;
      else if ((ur > 0) && (rr == 0) && (ru == 0))
	bad_due_to_ur++;
      else if ((ru > 0) && (ur == 0) && (rr == 0))
	bad_due_to_ru++;
      else
	bad_due_to_other_combinations++;
# if DEBUG_GAP_WALKER > 0
      fprintf(GlobalData->gwlogfp, "    This subgraph is NOT ok and contains %d uu, %d ur, %d ru, %d rr\n",
	      uu,
	      ur,
	      ru,
	      rr);
      fprintf(GlobalData->gwlogfp, "    bad due to rr %d, bad due to ur %d, bad due to ru %d, other %d\n",
	      bad_due_to_rr,
	      bad_due_to_ur,
	      bad_due_to_ru,
	      bad_due_to_other_combinations);
#endif
    } else if (uu + uniques != f0->size) {
# if DEBUG_GAP_WALKER > 0
      fprintf(GlobalData->gwlogfp, "    This subgraph is ok but contains %d uu, %d ur, %d ru, %d rr\n",
	      uu,
	      ur,
	      ru,
	      rr);
      fprintf(GlobalData->gwlogfp, "    bad due to rr %d, bad due to ur %d, bad due to ru %d, other %d\n",
	      bad_due_to_rr,
	      bad_due_to_ur,
	      bad_due_to_ru,
	      bad_due_to_other_combinations);
#endif
    }
    total++;
#   if DEBUG_GAP_WALKER > 1
    fprintf(GlobalData->gwlogfp,"    bad %d, total %d\n",
	    bad,
	    total);
#   endif
    total_gap += gap.mean;
#   if DEBUG_GAP_WALKER > 1
    fprintf(GlobalData->gwlogfp,
	    "    avg_gap_size %5.2f\n",
	    total_gap/(double)total);
#   endif
  }
# endif
    
# if 0
  //
  // check if I am doing ok (prints a .dot file of the "bad" gap)
  //
  {
    int
      i,
      cid,
      a_end,
      b_end,
      low,
      high;
    ChunkInstanceT
      * chunk;
    
    assert(lchunk != NULL);
    assert(rchunk != NULL);
    low = min(min(lchunk->aEndCoord, lchunk->bEndCoord),
	      min(rchunk->aEndCoord, rchunk->bEndCoord));
    high = max(max(lchunk->aEndCoord, lchunk->bEndCoord),
	       max(rchunk->aEndCoord, rchunk->bEndCoord));
    for (i = 0; i < f0->size; i++) {
      cid = f0->node[i].cid;
      chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
      assert(chunk != NULL);
      
      if (((chunk->aEndCoord >= low) && (chunk->aEndCoord <= high)) ||
	  ((chunk->bEndCoord >= low) && (chunk->bEndCoord <= high)))   
	continue;
      else {
	//
	// print f in .dot format
	//
#       if DEBUG_GAP_WALKER > 1
	fprintf(GlobalData->gwlogfp,"Oops: the chunk %d has been erroneusly selected. Range [%d,%d] (chunks %d, %d)\n",
		cid,
		low,
		high,
		left_cid,
		right_cid);
#       endif
	{
	  //
	  // print the subgraph in .dot format (it can produce the same graph
	  // multiple times, but it will overwrite the old file)
	  //
	  FILE *
	    dot_file;
	  char
	    filename[256];
	  
	  sprintf(filename, "./dot/%s.%d.%d.dot", GW_Filename_Prefix, left_cid, right_cid);
	  dot_file = file_open(filename, "w");
	  assert(dot_file != NULL);
	  Print_Dot_Subgraph(f, dot_file, filename, All_Edges);
	  fclose(dot_file);
	}
      }
    }
  }
# endif

  return f0;
}



int Inter_Scaffold_Gap_Walker(CDS_CID_t a_sid,
			      CDS_CID_t b_sid,
			      float (* quality)(CIEdgeT *, CDS_CID_t)) {
  //
  // Inter_Scaffold_Gap_walker() takes a scaffold id <A>, a scaffold
  // id <B>
  //
  // * it finds the gap length between the last chunk of <A> and
  //   the first chunk in <B> by looking for some scaffold edge
  //
  // * if none, then return FALSE (we don't have a way to know the gap size)
  //
  // * it builds the subgraph of all the chunks in the gap
  //
  // * it looks for a path from <A> to <B> in the subgraph
  //   such that the length of the path is consistent with
  //   the gap length
  //

  ChunkInstanceT
    * chunk;
  CIScaffoldT
    * a_s = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, a_sid),
    * b_s = GetGraphNode(ScaffoldGraph->ScaffoldGraph, b_sid);
  GraphEdgeIterator   SEdges;
  SEdgeT
    * edge;
  CDS_CID_t
    a_cid,
    b_cid,
    a_left_cid,
    a_right_cid,
    b_left_cid,
    b_right_cid,
    other_sid;
  int
    end;
  LengthT
    gap = {0.0, 0.0};

  assert(quality != NULL);
  assert(a_s != NULL);
  assert(b_s != NULL);
  assert(! isDeadCIScaffoldT(a_s));
  assert(! isDeadCIScaffoldT(b_s));
  assert(a_s->type == REAL_SCAFFOLD);
  assert(b_s->type == REAL_SCAFFOLD);

  //
  // try to estimate the distance between <A> and <B> by going through
  // the scaffold mate edges: what if there is no scaffold mate edge? 
  //
# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
	  "Looking to outgoing edges of scaffold %d\n",
	  a_sid);
# endif
  InitGraphEdgeIterator(ScaffoldGraph->ScaffoldGraph, a_sid, 
		     ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &SEdges);
  while((edge = NextGraphEdgeIterator(&SEdges))!= NULL) {
    assert(edge != NULL);

    //
    // if bogus ... discard!
    //
    if (isProbablyBogusEdge(edge))
      continue;

    //
    // get the other end
    //
    if (a_sid == edge->idA)
      other_sid = edge->idB;
    else
      other_sid = edge->idA;

#   if DEBUG_GAP_WALKER > 1
    fprintf(GlobalData->gwlogfp,
	    "* Edge %d to %d, length %5.2f, variance %5.2f, items %d\n",
	    a_sid, 
	    other_sid,
	    edge->distance.mean, 
	    edge->distance.variance,
	    edge->edgesContributing);
#   endif
    if ((other_sid == b_sid) &&
	(edge->distance.mean > MIN_SCAFF_EDGE)) {
      gap.mean = edge->distance.mean;
      assert(edge->distance.variance >= 0.0);
      gap.variance = edge->distance.variance;
      break;
    }
  }

  if (gap.mean == 0.0) {
#   if DEBUG_GAP_WALKER > 1
    fprintf(GlobalData->gwlogfp,
	    "no Sedge ... return\n");
#   endif
    return FALSE;
  }

  //
  // find the "beginning" of the gap (it could be the last or the
  // first chunk in scaffold <A>, it depends on the orientation)
  //
  a_left_cid = a_s->info.Scaffold.AEndCI;
  a_right_cid = a_s->info.Scaffold.BEndCI;
  
  //
  // find the "end" the gap (it could be the last or the
  // first chunk in scaffold <B>, it depends on the orientation)
  //
  b_left_cid = b_s->info.Scaffold.AEndCI;
  b_right_cid = b_s->info.Scaffold.AEndCI;
  
  //
  // now get the orientation of the two scaffolds
  // from the Sedge
  //
  assert(edge != NULL);
  switch (GetEdgeOrientationWRT(edge, a_sid)) {
  case AB_AB :
    a_cid = a_right_cid;
    b_cid = b_left_cid;
    break;
  case AB_BA :
    a_cid = a_right_cid;
    b_cid = b_right_cid;
    break;
  case BA_AB :
    a_cid = a_left_cid;
    b_cid = b_left_cid;
    break;
  case BA_BA :
    a_cid = a_left_cid;
    b_cid = b_right_cid;
    break;
  default :
    assert(0);
  }

# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
	  "a_left_cid %d, a_right_cid %d\n",
	  a_left_cid,
	  a_right_cid);
  fprintf(GlobalData->gwlogfp,
	  "b_left_cid %d, b_right_cid %d\n",
	  b_left_cid,
	  b_right_cid);
  fprintf(GlobalData->gwlogfp,
	  "orientation %s\n",
	  Orientation_As_String(GetEdgeOrientationWRT(edge, a_sid)));
  fprintf(GlobalData->gwlogfp,
	  "A is %d, B is %d \n",
	  a_cid,
	  b_cid);
# endif
 
  //
  // compute the orientation of chunk <a_cid>
  //
  chunk = GetGraphNode(ScaffoldGraph->RezGraph, a_cid);
  assert(chunk != NULL);
  if (chunk->aEndCoord <= chunk->bEndCoord)
    end = B_END;
  else
    end = A_END;
  
# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
	  "\nGap size mean %5.2f, variance %5.2f\n",
	  gap.mean,
	  gap.variance);
# endif

  //
  // TO DO 
  //

  return TRUE;
}

int Bac_Inter_Scaffold_Gap_Walker(CDS_CID_t a_sid, ChunkInstanceT *lchunk,
                                  CDS_CID_t b_sid, ChunkInstanceT *rchunk,
                                  float (* quality)(CIEdgeT *, CDS_CID_t),
                                  int walkLocale) {
  //
  // Bac_Inter_Scaffold_Gap_walker() takes a scaffold id <A>, a scaffold
  // id <B>
  //
  // * it finds the gap length between the last chunk of <A> and
  //   the first chunk in <B> by looking for some scaffold edge
  //
  // * it builds the subgraph of all the chunks that contain fragments from locale "walkLocale"
  //
  // * it looks for a path from <A> to <B> in the subgraph
  //   such that the length of the path is consistent with
  //   the estimated gap length (if any)
  //

  ChunkInstanceT
    * chunk;
  CIScaffoldT
    * a_scaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, a_sid),
    * b_scaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, b_sid);
  GraphEdgeIterator   SEdges;
  SEdgeT
    * edge;
  CDS_CID_t
    a_cid,
    b_cid,
    a_left_cid,
    a_right_cid,
    b_left_cid,
    b_right_cid,
    other_sid;
  int
    end;
  LengthT
    gap = {0.0, 0.0};

  assert(quality != NULL);
  assert(a_scaffold != NULL);
  assert(b_scaffold != NULL);
  assert(! isDeadCIScaffoldT(a_scaffold));
  assert(! isDeadCIScaffoldT(b_scaffold));
  assert(a_scaffold->type == REAL_SCAFFOLD);
  assert(b_scaffold->type == REAL_SCAFFOLD);

  //
  // try to estimate the distance between <A> and <B> by going through
  // the scaffold mate edges: what if there is no scaffold mate edge? 
  //
# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
	  "Looking to outgoing edges of scaffold %d\n",
	  a_sid);
# endif
  InitGraphEdgeIterator(ScaffoldGraph->ScaffoldGraph, a_sid, 
		     ALL_END, ALL_EDGES, ITERATOR_VERBOSE, &SEdges);
  while((edge = NextGraphEdgeIterator(&SEdges)) != NULL){
    assert(edge != NULL);

    //
    // if bogus ... discard!
    //
    if (isProbablyBogusEdge(edge))
      continue;

    //
    // get the other end
    //
    if (a_sid == edge->idA)
      other_sid = edge->idB;
    else
      other_sid = edge->idA;

#   if DEBUG_GAP_WALKER > 1
    fprintf(GlobalData->gwlogfp,
	    "* Edge %d to %d, length %5.2f, variance %5.2f, items %d\n",
	    a_sid, 
	    other_sid,
	    edge->distance.mean, 
	    edge->distance.variance,
	    edge->edgesContributing);
#   endif
    if ((other_sid == b_sid) &&
	(edge->distance.mean > MIN_SCAFF_EDGE)) {
      gap.mean = edge->distance.mean;
      assert(edge->distance.variance >= 0.0);
      gap.variance = edge->distance.variance;
      break;
    }
  }

  //
  // find the "beginning" of the gap (it could be the last or the
  // first chunk in scaffold <A>, it depends on the orientation)
  //
  a_left_cid = a_scaffold->info.Scaffold.AEndCI;
  a_right_cid = a_scaffold->info.Scaffold.BEndCI;
  
  //
  // find the "end" of the gap (it could be the last or the
  // first chunk in scaffold <B>, it depends on the orientation)
  //
  b_left_cid = b_scaffold->info.Scaffold.AEndCI;
  b_right_cid = b_scaffold->info.Scaffold.AEndCI;
  
  //
  // now get the orientation of the two scaffolds
  // from the Sedge
  //
  assert(edge != NULL);
  switch (GetEdgeOrientationWRT(edge, a_sid)) {
  case AB_AB :
    a_cid = a_right_cid;
    b_cid = b_left_cid;
    break;
  case AB_BA :
    a_cid = a_right_cid;
    b_cid = b_right_cid;
    break;
  case BA_AB :
    a_cid = a_left_cid;
    b_cid = b_left_cid;
    break;
  case BA_BA :
    a_cid = a_left_cid;
    b_cid = b_right_cid;
    break;
  default :
    assert(0);
  }

# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
	  "a_left_cid %d, a_right_cid %d\n",
	  a_left_cid,
	  a_right_cid);
  fprintf(GlobalData->gwlogfp,
	  "b_left_cid %d, b_right_cid %d\n",
	  b_left_cid,
	  b_right_cid);
  fprintf(GlobalData->gwlogfp,
	  "orientation %s\n",
	  Orientation_As_String(GetEdgeOrientationWRT(edge, a_sid)));
  fprintf(GlobalData->gwlogfp,
	  "A is %d, B is %d \n",
	  a_cid,
	  b_cid);
# endif

  // check to make sure these are the chunks we passed in as flanking the gap
  if ((a_cid != lchunk->id) || (b_cid != rchunk->id))
  {
	fprintf( stderr, "Sedge does not match locale information!\n");
	fprintf( stderr, "Sedge: a_cid: %d, b_cid: %d\n", a_cid, b_cid);
	fprintf( stderr, "Passed: lchunk->id: %d, rchunk->id: %d\n", lchunk->id, rchunk->id);
	assert(0);
  }
 
  //
  // compute the orientation of chunk <a_cid>
  //
  chunk = GetGraphNode(ScaffoldGraph->RezGraph, a_cid);
  assert(chunk != NULL);
  if (chunk->aEndCoord <= chunk->bEndCoord)
    end = B_END;
  else
    end = A_END;
  
# if DEBUG_GAP_WALKER > 1
  fprintf(GlobalData->gwlogfp,
	  "\nGap size mean %5.2f, variance %5.2f\n",
	  gap.mean,
	  gap.variance);
# endif

  //
  // TO DO 
  //

  return TRUE;
}



void Compute_Tentative_Position(chunk_subgraph * s,
                                CDS_CID_t start_cid, CDS_CID_t end_cid) {
  /* a simple minded way to compute a tentative position */
  /* on the chunks */
  /* */
  /* write the position in the A_end/B_end fields of the */
  /* chunk_subgraph */
  /* */
  ChunkInstanceT
    * u_chunk,
    * v_chunk,
    * start_chunk, 
    * end_chunk;
  int
    u,
    v,
    k,
    id,
    next_u;
  CIEdgeT
    * edge;
  LengthT gap_length,  
    gap_start;
#if 0
  int i;
#endif
  
  assert(s != NULL);
  
  start_chunk = GetGraphNode(ScaffoldGraph->RezGraph, start_cid);
  assert(start_chunk != NULL);
  end_chunk = GetGraphNode(ScaffoldGraph->RezGraph, end_cid);
  assert(end_chunk != NULL);
  // gap_length = Compute_Gap_Length(start_chunk, end_chunk, FALSE);
  gap_length = FindGapLength(start_chunk, end_chunk, FALSE);

  edge = s->table[start_cid]->best_edge;
  assert(edge != NULL);

  // don't worry about orientation, max is end going into gap
  gap_start.mean = max(start_chunk->offsetAEnd.mean, start_chunk->offsetBEnd.mean);
  gap_start.variance = max(start_chunk->offsetAEnd.variance, start_chunk->offsetBEnd.variance);

# if DEBUG_GAP_WALKER > 1
  fprintf(stderr, "in Compute_Tentative_Position, start_cid = %d, end_cid = %d\n", start_cid, end_cid);
  fprintf(stderr, "gap_start.mean = %f\n", gap_start.mean);
  fprintf(stderr, "gap_start.variance = %f\n", gap_start.variance);
  fprintf(stderr, "gap_length.mean = %f\n", gap_length.mean);
  fprintf(stderr, "gap_length.variance = %f\n", gap_length.variance);
#endif  
  Clear_All_Done_Bit(s);
  Set_Done_Bit(s, start_cid);  

  /* we are going to build from below, so start at start_cid  */

#if 0
  fprintf(stderr, "\ngraph size: %d\n", s->size);
  for (i = 0; i < s->size; i++) {
    u = s->node[i].cid;
    fprintf(stderr, "node[%d]: %d\n", i, u);
  }
#endif

  u = start_cid;
  while (u != end_cid) {
    float64 u_AMean, u_BMean, u_AVariance, u_BVariance;
    
    u_chunk = GetGraphNode(ScaffoldGraph->RezGraph, u);
    assert(u_chunk != NULL);

# if DEBUG_GAP_WALKER > 0
    fprintf(stderr, "\n******************Processing node: %d\n", u);
    for (k = 0; k < NUM_ORIENTATIONS; k++)
      {
	fprintf(stderr, "s->table[%d]->num_edges[orient[%d]] = %d\n", u, k, s->table[u]->num_edges[orient[k]]);
	for (id = 0; id < s->table[u]->num_edges[orient[k]]; id++) {
	  edge = s->table[u]->edge[orient[k]][id];
	  fprintf(stderr, "orientation: %d edge: %d to %d\n", k, edge->idA, edge->idB);
	}
      }
#endif
	
    /* get the edge */
    edge = s->table[u]->best_edge;
    assert(edge != NULL);
    
    /* get the other end */
    if (u == edge->idA)
      v = edge->idB;
    else
      v = edge->idA;
    next_u = v;
    
    v_chunk = GetGraphNode(ScaffoldGraph->RezGraph, v);
    assert(v_chunk != NULL);

# if DEBUG_GAP_WALKER > 1
    fprintf(stderr, "u = %d, v = %d\n", u, v);
    fprintf(stderr, "edge->idA = %d, edge->idB = %d, orient = %s, length = %.1f, variance = %.1f\n",
	    edge->idA, edge->idB, Orientation_As_String(edge->orient), 
	    edge->distance.mean, edge->distance.variance); 
    if (edge->flags.bits.hasContributingOverlap) 
      fprintf(stderr, "edge->flags.bits.hasContributingOverlap\n"); 
    if (edge->flags.bits.hasRepeatOverlap) 
      fprintf(stderr, "edge->flags.bits.hasRepeatOverlap\n"); 
    if (edge->flags.bits.hasTandemOverlap) 
      fprintf(stderr, "edge->flags.bits.hasTandemOverlap\n"); 
#endif	    
    
    if ( u_chunk == start_chunk ) { /* get the "real" means and variances */
      u_AMean = u_chunk->offsetAEnd.mean;
      u_AVariance = u_chunk->offsetAEnd.variance;
      u_BMean = u_chunk->offsetBEnd.mean;
      u_BVariance = u_chunk->offsetBEnd.variance;
      s->table[u]->A_end.variance = u_AVariance;
      s->table[u]->B_end.variance = u_BVariance;
      s->table[u]->A_end.mean = u_AMean;
      s->table[u]->B_end.mean = u_BMean;

    }
    else {	/* use the ones we have previously calculated for this chunk */
      u_AMean = s->table[u]->A_end.mean;
      u_AVariance = s->table[u]->A_end.variance;
      u_BMean = s->table[u]->B_end.mean;
      u_BVariance = s->table[u]->B_end.variance;
    }
    
    /* all the orientations are from idA's point of view */
    /* k = 0 => AB_AB, k = 1 => AB_BA, k = 2 => BA_AB, k = 3 => BA_BA */
# if DEBUG_GAP_WALKER > 0
    if (u == edge->idA) {
      fprintf(stderr, "uv %d	%d	%d	%f	%f", 
	      u, v, k, edge->distance.mean, edge->distance.variance); 
      fprintf(stderr, "	%d	%d\n", 
	      (s->table[u_chunk->id]->done) ? 1 : 0, (s->table[v_chunk->id]->done) ? 1 : 0);
    }
    else {
      fprintf(stderr, "vu %d	%d	%d	%f	%f", 
	      v, u, k, edge->distance.mean, edge->distance.variance); 
      fprintf(stderr, "	%d	%d\n", 
	      (s->table[v_chunk->id]->done) ? 1 : 0, (s->table[u_chunk->id]->done) ? 1 : 0);
    }
#endif	    
    /* k is computed based on node u's point of view */
    if (TRUE /*!Is_Unique(v_chunk)*/ ) {
      if (u == edge->idA) {
	switch (edge->orient) {
	case 'N':		/* AB_AB */
	  if (u_AMean < u_BMean) {
	    s->table[v]->A_end.mean = u_BMean + edge->distance.mean;  
	    s->table[v]->B_end.mean =  s->table[v]->A_end.mean + v_chunk->bpLength.mean;
	  }
	  else {
	    s->table[v]->A_end.mean = u_BMean - edge->distance.mean;  
	    s->table[v]->B_end.mean =  s->table[v]->A_end.mean - v_chunk->bpLength.mean;
	  }
	  break;
	case 'I':		/* AB_BA */
	  if (u_AMean < u_BMean) {
	    s->table[v]->B_end.mean = u_BMean + edge->distance.mean;   
	    s->table[v]->A_end.mean =  s->table[v]->B_end.mean + v_chunk->bpLength.mean;
	  }
	  else {
	    s->table[v]->B_end.mean = u_BMean - edge->distance.mean; 
	    s->table[v]->A_end.mean =  s->table[v]->B_end.mean - v_chunk->bpLength.mean;
	  }
	  break;
	case 'O':		/* BA_AB */
	  if (u_AMean < u_BMean) {
	    s->table[v]->A_end.mean = u_AMean - edge->distance.mean;  
	    s->table[v]->B_end.mean =  s->table[v]->A_end.mean - v_chunk->bpLength.mean;
	  }
	  else {
	    s->table[v]->A_end.mean = u_AMean + edge->distance.mean;  
	    s->table[v]->B_end.mean = s->table[v]->A_end.mean + v_chunk->bpLength.mean;
	  }
	  break;
	case 'A':		/* BA_BA */
	  if (u_AMean < u_BMean) {
	    s->table[v]->B_end.mean = u_AMean - edge->distance.mean; 
	    s->table[v]->A_end.mean = s->table[v]->B_end.mean - v_chunk->bpLength.mean;
	  }
	  else {
	    s->table[v]->B_end.mean = u_AMean + edge->distance.mean; 
	    s->table[v]->A_end.mean = s->table[v]->B_end.mean + v_chunk->bpLength.mean;
	  }
	  break;
        default:
          assert(0);
          break;
	}
      }
      else { /* v == edge->idA */
	switch (edge->orient) {
	case 'N':		/* AB_AB == v_u */
	  if (u_AMean < u_BMean) {
	    s->table[v]->B_end.mean = u_AMean - edge->distance.mean;
	    s->table[v]->A_end.mean =  s->table[v]->B_end.mean - v_chunk->bpLength.mean;
	  }
	  else {
	    s->table[v]->B_end.mean = u_AMean + edge->distance.mean; 
	    s->table[v]->A_end.mean =  s->table[v]->B_end.mean + v_chunk->bpLength.mean;
	  }
	  break;
	case 'I':		/* AB_BA */
	  if (u_AMean < u_BMean) {
	    s->table[v]->B_end.mean = u_BMean + edge->distance.mean;   
	    s->table[v]->A_end.mean =  s->table[v]->B_end.mean + v_chunk->bpLength.mean;
	  }
	  else {
	    s->table[v]->B_end.mean = u_BMean - edge->distance.mean;  
	    s->table[v]->A_end.mean =  s->table[v]->B_end.mean - v_chunk->bpLength.mean;
	  }
	  break;
	case 'O':		/* BA_AB */
	  if (u_AMean < u_BMean) {
	    s->table[v]->A_end.mean = u_AMean - edge->distance.mean; 
	    s->table[v]->B_end.mean =  s->table[v]->A_end.mean - v_chunk->bpLength.mean;
	  }
	  else {
	    s->table[v]->A_end.mean = u_AMean + edge->distance.mean; 
	    s->table[v]->B_end.mean = s->table[v]->A_end.mean + v_chunk->bpLength.mean;
	  }
	  break;
	case 'A':		/* BA_BA == v_u */
	  if (u_AMean < u_BMean) {
	    s->table[v]->A_end.mean = u_BMean + edge->distance.mean;  
	    s->table[v]->B_end.mean = s->table[v]->A_end.mean + v_chunk->bpLength.mean;
	  }
	  else {
	    s->table[v]->A_end.mean = u_BMean - edge->distance.mean;  
	    s->table[v]->B_end.mean = s->table[v]->A_end.mean - v_chunk->bpLength.mean;
	  }
	  break;
        default:
          assert(0);
          break;
	}
      }
     

      {
	double edgeVariance = ComputeFudgeVariance(edge->distance.mean); 

	if( s->table[v]->A_end.mean < s->table[v]->B_end.mean )
	  {
	    if( s->table[u]->A_end.mean < s->table[u]->B_end.mean )
	      {
		/*
		  u                              v
		  ------------------>          --------------------->

		*/
		s->table[v]->A_end.variance = s->table[u]->B_end.variance + edgeVariance;
		s->table[v]->B_end.variance = s->table[v]->A_end.variance + v_chunk->bpLength.variance;
	      }
	    else
	      {
		/*
		  u                              v
		  <------------------          --------------------->


		*/
		s->table[v]->A_end.variance = s->table[u]->A_end.variance + edgeVariance;
		s->table[v]->B_end.variance = s->table[v]->A_end.variance + v_chunk->bpLength.variance;
	      }
	  }
	else
	  {
	    if( s->table[u]->A_end.mean < s->table[u]->B_end.mean )
	      {
		/*
		  u                              v
		  ------------------>          <---------------------

		*/
		s->table[v]->B_end.variance = s->table[u]->B_end.variance + edgeVariance;
		s->table[v]->A_end.variance = s->table[v]->B_end.variance + v_chunk->bpLength.variance;  
	      }
	    else
	      {
		/*
		  u                              v
		  <-----------------           <---------------------

		*/
		s->table[v]->B_end.variance = s->table[u]->A_end.variance + edgeVariance;
		s->table[v]->A_end.variance = s->table[v]->B_end.variance + v_chunk->bpLength.variance;
	      }
	  }

      }
       /*
       s->table[v]->A_end.variance = gap_start.variance + ComputeFudgeVariance(s->table[v]->A_end.mean - gap_start.mean);

      assert (s->table[v]->A_end.variance > 0);

      s->table[v]->B_end.variance = gap_start.variance + ComputeFudgeVariance(s->table[v]->B_end.mean - gap_start.mean);
      */

    }

# if DEBUG_GAP_WALKER > 0	 
    fprintf(stderr, "placing u = %d, s->table[%d]->A_end.mean = %f, variance = %f\n", 
	    u, u, s->table[u]->A_end.mean, s->table[u]->A_end.variance);
    fprintf(stderr, "placing u = %d, s->table[%d]->B_end.mean = %f, variance = %f\n", 
	    u, u, s->table[u]->B_end.mean, s->table[u]->B_end.variance);

    fprintf(stderr, "placing v = %d, s->table[%d]->A_end.mean = %f, variance = %f\n", 
	    v, v, s->table[v]->A_end.mean, s->table[v]->A_end.variance);
    fprintf(stderr, "placing v = %d, s->table[%d]->B_end.mean = %f, variance = %f\n", 
	    v, v, s->table[v]->B_end.mean, s->table[v]->B_end.variance);
	    
#endif
    assert (s->table[v]->B_end.variance > 0);


    if (next_u == end_cid)
      {
	// figure out where to put end_cid based on walk to this point
	// it should start at the endpoint of the current u_chunk minus the overlap length of the current edge
	LengthT new_gap_length, delta;
# if DEBUG_GAP_WALKER > 0
	fprintf(stderr, "abcedf, edge->orient = %c, edge->distance.mean = %f, gap_start.mean = %f\n, gap_start.variance = %f\n", 
		edge->orient, edge->distance.mean, gap_start.mean,gap_start.variance); 
	
	
	fprintf(stderr, "NEWGAP, A_END = %f, B_END = %f,  edgedist = %f,\n",s->table[u]->A_end.mean,s->table[u]->B_end.mean,edge->distance.mean);
#endif

	new_gap_length.mean = max(s->table[u]->A_end.mean + edge->distance.mean - gap_start.mean,
				  s->table[u]->B_end.mean + edge->distance.mean - gap_start.mean);

	delta.mean = new_gap_length.mean - gap_length.mean;
	
# if DEBUG_GAP_WALKER > 0
	fprintf(stderr, "abcedf, new_gap_length.mean = %f, gap_length.mean = %f\n", 
		new_gap_length.mean, gap_length.mean);
#endif 
	// figure out what the variance of the destination chunk should be, subtract off current to get the delta


	delta.variance = max(s->table[u]->A_end.variance,s->table[u]->B_end.variance)+
	                 ComputeFudgeVariance(edge->distance.mean)-
	                 min(v_chunk->offsetAEnd.variance,v_chunk->offsetBEnd.variance);
	/* delta.variance = gap_start.variance + ComputeFudgeVariance(new_gap_length.mean) - min(v_chunk->offsetAEnd.variance,v_chunk->offsetBEnd.variance);
	 */
# if DEBUG_GAP_WALKER > 0
	{
	  fprintf(stderr, "abcedf, %d->offsetAEnd.mean = %f\n", v, v_chunk->offsetAEnd.mean);
	  fprintf(stderr, "abcedf, %d->offsetAEnd.variance = %f\n", v, v_chunk->offsetAEnd.variance);
	  fprintf(stderr, "abcedf, %d->offsetBEnd.mean = %f\n", v, v_chunk->offsetBEnd.mean);
	  fprintf(stderr, "abcedf, %d->offsetBEnd.variance = %f\n", v, v_chunk->offsetBEnd.variance);

	  fprintf(stderr, "abcedf, delta.mean = %f\n", delta.mean);
	  fprintf(stderr, "abcedf, delta.variance = %f\n", delta.variance);
	}
#endif

	assert(v_chunk->offsetBEnd.variance != 0);
	assert(v_chunk->offsetAEnd.variance != 0);

	
	{
# if DEBUG_GAP_WALKER > 0
	  // CIScaffoldT *scaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds,v_chunk->scaffoldID);
	  fprintf(stderr, "BEFORE AddDelta %d \n",v_chunk->scaffoldID);

	  /*	  DumpCIScaffold(stderr,ScaffoldGraph,scaffold, FALSE);*/
	  fprintf(stderr,"**** FORCE variances before AddDelta in from = %d, to = %d  ****\n",start_cid,end_cid );	 
#endif

	  Force_Increasing_Variances();
	  AddDeltaToScaffoldOffsets(ScaffoldGraph,
				    v_chunk->scaffoldID,
				    v,
				    1,
				    0,
				    delta);
# if DEBUG_GAP_WALKER > 0
#if 1
	  fprintf(stderr, "AFTER AddDelta\n");
	  fprintf(stderr,"**** FORCE variances after AddDelta in from = %d, to = %d  ****\n",start_cid,end_cid );
#endif
#endif

	  Force_Increasing_Variances();
	  /*  DumpCIScaffold(stderr,ScaffoldGraph,scaffold, FALSE);*/
	}

# if DEBUG_GAP_WALKER > 0
	{
	  fprintf(stderr, "abcedf, %d->offsetAEnd.mean = %f\n", v, v_chunk->offsetAEnd.mean);
	  fprintf(stderr, "abcedf, %d->offsetAEnd.variance = %f\n", v, v_chunk->offsetAEnd.variance);
	  fprintf(stderr, "abcedf, %d->offsetBEnd.mean = %f\n", v, v_chunk->offsetBEnd.mean);
	  fprintf(stderr, "abcedf, %d->offsetBEnd.variance = %f\n", v, v_chunk->offsetBEnd.variance);
	  fprintf(stderr, "abcedf\n");
	}
#endif
      }
    u = next_u;
  }
}



int Count_TransChunk_Edges(chunk_subgraph * s)
{
  ChunkInstanceT
    * u_chunk;
  int
    i,
    u,
    k,
    id,
    num_transchunk_edges = 0;
  CIEdgeT
    * edge;
  
  assert(s != NULL);
  
  for (i = 0; i < s->size; i++) {
    /* grab the next node */
    u = s->node[i].cid;

    u_chunk = GetGraphNode(ScaffoldGraph->RezGraph, u);
    assert(u_chunk != NULL);

    /* scan the adjacency of u */
    for (k = 0; k < NUM_ORIENTATIONS; k++) {
      for (id = 0; id < s->table[u]->num_edges[orient[k]]; id++) {

	/* get the edge */
	edge = s->table[u]->edge[orient[k]][id];
	assert(edge != NULL);

	if (isTransChunkEdge(edge))
	  num_transchunk_edges++;
      }
    }
  }
  return num_transchunk_edges;
}
	  


int Count_Edges(chunk_subgraph * s)
{
  ChunkInstanceT
    * u_chunk;
  int
    i,
    u,
    k,
    id,
    num_edges = 0;
  CIEdgeT
    * edge;
  
  assert(s != NULL);
  
  for (i = 0; i < s->size; i++) {
    /* grab the next node */
    u = s->node[i].cid;

    u_chunk = GetGraphNode(ScaffoldGraph->RezGraph, u);
    assert(u_chunk != NULL);

    /* scan the adjacency of u */
    for (k = 0; k < NUM_ORIENTATIONS; k++) {
      for (id = 0; id < s->table[u]->num_edges[orient[k]]; id++) {

	/* get the edge */
	edge = s->table[u]->edge[orient[k]][id];
	assert(edge != NULL);

	num_edges++;
      }
    }
  }
  return num_edges;
}
	  


void Shortest_Path(chunk_subgraph * s,
		   CDS_CID_t begin_cid,
		   CDS_CID_t end_cid,
		   float (* quality)(CIEdgeT *, CDS_CID_t)) {
  //
  // Shortest_Path() finds the shortest paths from <begin_cid> to any
  // node in <s> by calling Dijkstra() and it prints the informations
  // into a "*.sp.cam" file
  //
# if CREATE_CAM_FILE > 0
  char
    * Colour[SP_NUM_COLOURS] = {
      "CFFFF00 T4 S # start/end",         //  0
      "CFF0000 T2 S # links",
      "CFFFFFF T2 S #   0.0_<_d_<__30.0", //  2
      "CEEEEEE T2 S #  30.0_<_d_<__60.0",
      "CDDDDDD T2 S #  60.0_<_d_<__90.0", //  4
      "CCCCCCC T2 S #  90.0_<_d_<_120.0",
      "CBBBBBB T2 S # 120.0_<_d_<_150.0", //  6
      "CAAAAAA T2 S # 150.0_<_d_<_180.0",
      "C999999 T2 S # 180.0_<_d_<_210.0", //  8
      "C888888 T2 S # 210.0_<_d_<_240.0",
      "C777777 T2 S # 240.0_<_d_<_270.0", // 10
      "C666666 T2 S # 270.0_<_d_<_300.0",
      "C555555 T2 S # 300.0_<_d_<_330.0", // 12
      "C444444 T2 S # 330.0_<_d_<_360.0",
      "C333333 T2 S # 360.0_<_d_<_360.0", // 14
      "C222222 T2 S # 390.0_<_d_<_410.0",
      "C00AA00 T2 S # d_>_410.0",         // 16
      "CFFFFFF T4 S #   0.0_<_d_<__30.0__unique",
      "CEEEEEE T4 S #  30.0_<_d_<__60.0__unique",
      "CDDDDDD T4 S #  60.0_<_d_<__90.0__unique",
      "CCCCCCC T4 S #  90.0_<_d_<_120.0__unique",
      "CBBBBBB T4 S # 120.0_<_d_<_150.0__unique", 
      "CAAAAAA T4 S # 150.0_<_d_<_180.0__unique",
      "C999999 T4 S # 180.0_<_d_<_210.0__unique",
      "C888888 T4 S # 210.0_<_d_<_240.0__unique",
      "C777777 T4 S # 240.0_<_d_<_270.0__unique",
      "C666666 T4 S # 270.0_<_d_<_300.0__unique",
      "C555555 T4 S # 300.0_<_d_<_330.0__unique",
      "C444444 T4 S # 330.0_<_d_<_360.0__unique",
      "C333333 T4 S # 360.0_<_d_<_360.0__unique",
      "C222222 T4 S # 390.0_<_d_<_410.0__unique",
      "C00AA00 T4 S # d_>_410.0__unique"};
  CDS_CID_t
    i,
    id,
    a_end,
    b_end;
  int
    colour;
  FILE
    * walker_cam_file;
  char
    unique[STR_LEN],
    orient[STR_LEN],
    filename[STR_LEN];
# endif
  ChunkInstanceT  
    * begin_chunk = GetGraphNode(ScaffoldGraph->RezGraph, begin_cid),
    * end_chunk = GetGraphNode(ScaffoldGraph->RezGraph, end_cid);
# if CREATE_CAM_FILE > 0
  ChunkInstanceT  
    * chunk;
#endif

  assert(s != NULL);
  assert(quality != NULL);
  assert(begin_chunk != NULL);
  assert(end_chunk != NULL);

# if CREATE_CAM_FILE > 0
  //
  // open the cam file
  //
  sprintf(filename, "./cam/%s.%d.%d.sp.cam", GW_Filename_Prefix, begin_cid, end_cid);
  walker_cam_file = file_open (filename, "w");
  assert(walker_cam_file != NULL);

  //
  // output the colours
  //
  for (i = 0; i < SP_NUM_COLOURS; i ++)
    fprintf(walker_cam_file, "%d: %s\n", i, Colour[i]);
# endif

  //
  // if the start/end chunk are not in <s>, add them
  //
  Add_Node(s, begin_cid, Is_Not_Bogus);
  Add_Node(s, end_cid, Is_Not_Bogus);

  //
  // search for the shortest path
  //
  Dijkstra(s, begin_cid, end_cid, ALL_END, quality);

# if CREATE_CAM_FILE > 0
  //
  // now print all the chunks
  //
  for (i = 0; i < s->size; i++) {
    //
    // compute the color
    //
    if (s->node[i].d >= 410.0)
      colour = 16;
    else
      colour = (int)(s->node[i].d / 30.0) + 2;
    id = s->node[i].cid;
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, id);
    assert(chunk != NULL);
    if (Is_Unique(chunk)) {
      strcpy(unique, "unique");
      colour += 15;
    } else {
      strcpy(unique, "not unique");
    }
    if ((id == begin_cid) || (id == end_cid))
      colour = 0;
    if  (chunk->aEndCoord <= chunk->bEndCoord) {
      strcpy(orient, "direct");
      a_end = chunk->aEndCoord;
      b_end = chunk->bEndCoord;
    } else {
      a_end = chunk->bEndCoord;
      b_end = chunk->aEndCoord;
      strcpy(orient, "reverse");
    }

    fprintf(walker_cam_file,
	    "%d: %d A%d %d # %s %s chunk %d - bn distance from %d %5.2f, sow distance %5.2f, bp distance %5.2f\n",
	    id,
	    a_end,
	    colour,
	    b_end,
	    orient,
	    unique,
	    id,
	    begin_cid,
	    s->node[i].d_neck,
	    s->node[i].d,
	    s->node[i].distance.mean);
  }

  //
  // and print the best path (backward)
  //
  i = end_cid;
  while (i != begin_cid) {
    assert(s->table[i] != NULL);
    if (s->table[i]->best_edge == NULL)
      break;
    //if (! isOverlapEdge(s->table[i]->best_edge)) {
      fprintf (walker_cam_file, "LNK: %d %d A%d # edge quality %5.2f, weight %d, mean %f, var %f\n",
	       s->table[i]->path_parent,
	       i,
	       1,
	       quality(s->table[i]->best_edge, begin_cid),
               s->table[i]->best_edge->edgesContributing,
	       s->table[i]->best_edge->distance.mean,
	       s->table[i]->best_edge->distance.variance);
      //}
    //
    // follow the parent
    //
    i = s->table[i]->path_parent;
  }


  fclose (walker_cam_file);
# endif
}



void Dijkstra(chunk_subgraph * s,
	      CDS_CID_t source_id,
	      CDS_CID_t sink_id,
	      int end,
	      float (* quality)(CIEdgeT *, CDS_CID_t)) {
  //
  // Dijkstra's single source shortest path (see, e.g., CLR page 527)
  //
  // The substantial difference is the definition of path lenght that
  // we assume as the maximum edge cost on the edges of the path
  // (bottleneck distance); we break ties by using the standard
  // distance (sum of edges weights on the path)
  //
  // NOTE: this version is directional, that is if enters a chunk
  // from the A end it goes out from the B end and viceversa
  //
  // The sink is used only to avoid to pick the trivial path source->sink
  //
  int
    i,
    k,
    u = -1,
    v,
    id,
    new_u = -1,
    num_nodes;
  double
    q,
    min;
  CIEdgeT
    * edge;
  ChunkInstanceT 
    * v_chunk;
  assert(s != NULL);
  assert(quality != NULL);

  //
  // initialization
  //
  for (i = 0; i < s->size; i++) {
    s->node[i].d = FLT_MAX;
    s->node[i].d_neck = FLT_MAX;
    s->node[i].distance.mean = 0.0;
  }
  assert(source_id < s->max);
  assert(s->table[source_id] != NULL);
  s->table[source_id]->d = 0.0;
  s->table[source_id]->d_neck = 0.0;
  s->table[source_id]->end = end;
  Clear_All_Path_Bit(s);

# if DEBUG_GAP_WALKER > 2
  fprintf(GlobalData->gwlogfp,"** initialized ok\n");
# endif

  //
  // repeat for all the nodes
  //
  num_nodes = s->size;
  while (num_nodes > 0) {
#   if DEBUG_GAP_WALKER > 2
    fprintf(GlobalData->gwlogfp, "\nnum_nodes = %d \n",
	    num_nodes);
#   endif

    //
    // find the minimum in the vertices not marked
    //
    min = FLT_MAX;
    //new_u = u;
    for (i = 0; i < s->size; i++)
      if ((! s->node[i].path_bit) && (s->node[i].d < min)) {
	min = s->node[i].d;
	new_u = s->node[i].cid;
      }

    //
    // avoid looping
    //
    if (u == new_u)
      break;
    u = new_u;

#   if DEBUG_GAP_WALKER > 2
    fprintf(GlobalData->gwlogfp,"At node u = %d, mindist = %5.2f, end = %d\n",
	    u,
	    min,
	    s->table[u]->end);
#   endif

    //
    // mark the node <u>
    //
    Set_Path_Bit(s, u);
    num_nodes--;

    //
    // scan the adjancency of <u>
    //
    for (k = start_index[s->table[u]->end]; k < end_index[s->table[u]->end]; k++)
      for (id = 0; id < s->table[u]->num_edges[orient[k]]; id++) {
	//
	// get the edge
	//
	edge = s->table[u]->edge[orient[k]][id];
	assert(edge != NULL);
	//
	// get the other end
	//
	if (u == edge->idA)
	  v = edge->idB;
	else
	  v = edge->idA;
	v_chunk = GetGraphNode(ScaffoldGraph->RezGraph, v);
	assert(v_chunk != NULL);

	//
	// if trivial path ... get the next edge
	// 
	if ((u == source_id) && (v == sink_id))
	  continue;

	//
	// avoid self loops (if any)
	//
	if (v == u)
	  continue;
	
#       if DEBUG_GAP_WALKER > 3
	fprintf(GlobalData->gwlogfp,"edge from %d to %d (quality %5.2f) orientation %s ",
		u,
		v,
		quality(edge,u),
	        Orientation_As_String(GetEdgeOrientationWRT(edge, u)));
#       endif

	//
        // now .... relax the node <v> where the
	// path length is max edge cost on the edges of the path
	//
        assert(s->table[v] != NULL);
	q = quality(edge, u);

	//
	// if the new d_neck distance is better, I choose the edge
	// if they are equal I break the ties with the standard sum of weights
	//
	if ((s->table[v]->d_neck > max(s->table[u]->d_neck, q)) ||
	    ((s->table[v]->d_neck == max(s->table[u]->d_neck, q)) &&
	     (s->table[v]->d > s->table[u]->d + q))) {
	  s->table[v]->end = new_end[orient[k]];
	  s->table[v]->d_neck = max(s->table[u]->d_neck, q);
	  s->table[v]->d = s->table[u]->d + q;
	  s->table[v]->path_parent = u;
	  s->table[v]->best_edge = edge;
	  s->table[v]->distance.mean = s->table[u]->distance.mean +
	    edge->distance.mean + v_chunk->bpLength.mean;
	  s->table[v]->distance.variance = s->table[u]->distance.variance +
	    edge->distance.variance + v_chunk->bpLength.variance;
	}
    }
  }
}

int Find_Greedy_Path(chunk_subgraph * subgraph,
		     CDS_CID_t from_cid,
		     CDS_CID_t source,
		     CDS_CID_t destination,
		     LengthT * max_distance,
		     LengthT travelled_distance,
		     int level,
		     int end,
		     float (* quality)(CIEdgeT *, CDS_CID_t),
		     float qualityThresh,
		     int (* terminate)(ChunkInstanceT *, CDS_CID_t, chunk_subgraph *),
		     int colour,
		     int verbose,
		     long *currentCalls,
		     long maxCalls,
		     float *tooShort,
		     float *tooLong) {
  //
  // Find_Greedy_Path() searches for a path in a subgraph <subgraph> of
  // chunks from the chunk <source_cid> to the chunk <destination> (or
  // to a unique, depending on the function terminate() that you are
  // passing) which has maximal length <max_distance> (if don't want
  // to check for distance use a huge value).  The search goes in the
  // direction specified by <end>.  At each node, Find_Greedy_Path()
  // evaluates the quality of each outgoing edge, and pick greedly the
  // best one that has not been visited already
  //
  // Notes:
  //
  //   * In the initial call, use from_cid == source
  //
  //   * it will NOT consider a path the trivial path source->destination
  //
  //   * <level> should be set to 0 (correspond to the recursion level)
  //
  //   * the <quality()> function takes an edge and should evaluate to a float score
  //     (see Edge_Quality2() as example)
  //
  //   * the <terminate()> function tells the procedure when to stop
  //     (see Stop_At_Unique() or Stop_At_The_Other_End() as examples)
  //
  //   * the path is directional, that is if enters a chunk
  //     from the A end it goes out from the B end and viceversa
  //
  //
  // Returns:
  //   * 0 (==FALSE) if a path has been found
  //   * l, the number of hops
  //
  // Modifies:
  //   * the subgraph->table[]->path* fields to record the path(s) found (see
  //     the defn in GapWalkerREZ.h) and subgraph->table[]->distance where we store
  //     the travelled distance
  //
  //

  // Knut Reinert: (use of GlobalData->walkLevel)
  // here in short
  // level 1 : standard gap walking
  // level 2 : agressive gap walking
  // level 3 : force walking
  // level 4 : min bottleneck (not done yet)


  CIEdgeT
    * e;
  ChunkInstanceT
    * from_chunk = GetGraphNode(ScaffoldGraph->RezGraph, from_cid),
    * source_chunk = GetGraphNode(ScaffoldGraph->RezGraph, source),
    * to_chunk   = NULL; 
  edge_quality
    * eq2,
    * eq;
  CDS_CID_t
    id,
    tid;
  int
    i,
    j,
    k,
    o,
    hops,
    num_edges = 0,
    edge_start_index = 0;
  float64
    neg_num_std_deviations,
    pos_num_std_deviations;
#   if DEBUG_GAP_WALKER > 5
  int
    do_check_distance;
  float64
    add_mean,
    add_variance;
#endif

  assert(subgraph != NULL);
  assert(end != 0);
  assert(from_chunk != NULL);
  
# if DEBUG_GAP_WALKER > 2
  if (verbose) {
    fprintf(GlobalData->gwlogfp,
	    "\n");
    for (i = 0; i < level; i++) 
      fprintf(GlobalData->gwlogfp,
	      " ");
    if (end == B_END)
      fprintf(GlobalData->gwlogfp,
	      "* at the B");
    else if (end == A_END)
      fprintf(GlobalData->gwlogfp,
	      "* at the A");
    else 
      fprintf(GlobalData->gwlogfp,
	      "* at the A/B");
    fprintf(GlobalData->gwlogfp,
	    " end of node %d (mean distance %5.2f, var %5.2f)\n",
	    from_cid,
	    travelled_distance.mean,
	    travelled_distance.variance);
  }

  // FragInfo(rchunk);
# endif


  // We increment the counter of the overall edges explored
  (*currentCalls)++;

# if DEBUG_GAP_WALKER > 2
  fprintf(GlobalData->gwlogfp,
	  "\nCurrent Calls = %ld maxCalls = %ld\n",*currentCalls,maxCalls);
#endif

  // if we have explored to much we return;
  if( *currentCalls > maxCalls)
    {
      source_chunk->flags.bits.walkMaxedOut = TRUE;
      return FALSE;
    }

  //
  // compute the size of the edge quality array
  //

  for (k = start_index[end]; k < end_index[end]; k++)
    num_edges += subgraph->table[from_cid]->num_edges[orient[k]];
  
  //
  // no outgoing edge? => return
  //
  if (num_edges == 0)
    {
      // if( level == 0 )
#if DEBUG_GAP_WALKER > -1
	fprintf(GlobalData->gwlogfp,
	      "$$$ No outgoing edge\n");
#endif
      return FALSE;
    }

  //
  // allocate the array
  //
  eq2 = eq = (edge_quality *)safe_calloc(num_edges, sizeof(edge_quality));

  //
  // get the edges with the correct orientation and store them
  // in the <eq> array
  //
  j = 0;
  for (k = start_index[end]; k < end_index[end]; k++) {
    o = orient[k];
    for (id = 0; id < subgraph->table[from_cid]->num_edges[o]; id++) {
      e = subgraph->table[from_cid]->edge[o][id];
      assert(e != NULL);
      eq[j].orientation = o;
      eq[j].edge = e;
      eq[j].quality = quality(e, from_cid);
      j++;
    }
  }
  assert(j == num_edges);

  //
  // sort w.r.t. the quality
  //
  qsort(eq, num_edges, sizeof(edge_quality), Compare_Edge_Quality); 

# if DEBUG_GAP_WALKER > 2
  if (verbose) {
    for (j = 0; j < num_edges; j++) {
      for (i = 0; i < level; i++) 
	fprintf(GlobalData->gwlogfp," ");
      fprintf(GlobalData->gwlogfp,
	      "> edge (%d,%d), quality %5.2f, orient %d, distance %5.2f containment relation aCb %d bCa %d\n",
	      eq[j].edge->idA,
	      eq[j].edge->idB,
	      eq[j].quality,
	      eq[j].orientation,
	      eq[j].edge->distance.mean,
	      eq[j].edge->flags.bits.aContainsB,
	      eq[j].edge->flags.bits.bContainsA);
    }
  }

  if (verbose) {
    int temp_end;
    if (end == B_END) 
      temp_end = A_END;
    else
      temp_end = B_END;
      
    for (k = start_index[temp_end]; k < end_index[temp_end]; k++) {
      o = orient[k];
      for (id = 0; id < subgraph->table[from_cid]->num_edges[o]; id++) {
	e = subgraph->table[from_cid]->edge[o][id];
	assert(e != NULL);
	for (i = 0; i < level; i++) 
	  fprintf(GlobalData->gwlogfp, " ");
	fprintf(GlobalData->gwlogfp,
		"> unused edge (idA, idB) (%d,%d), quality %5.2f, orient %d, distance %5.2f\n",
		e->idA,
		e->idB,
		quality(e, from_cid),
		o,
		e->distance.mean);
      }
    }
  }
#endif

  edge_start_index = 0;
  //
  // iterate through all the edges looking for the destination
  //
#if 0
  for (j = 0; j < num_edges; j++) {
    //
    // get the edges 
    //
    e = eq[j].edge;
    o = eq[j].orientation;
    assert(e != NULL);


    // check to see if other end is destination
    if( (e->idA == destination &&  e->idB != source)  
	|| (e->idB == destination && e->idA != source) )
      {
	float64 gap_estimate;
	float64 lower_bound;

	if ( GlobalData->walkLevel == 1 )
	  neg_num_std_deviations = 3.0;
	else 
	  neg_num_std_deviations = 300000.0;  
	
	// can't be too short a walk if walkLevel != 0
	// only shortcut if we have travelled far enough
	
	gap_estimate = travelled_distance.mean+e->distance.mean;
	lower_bound  = max_distance->mean-neg_num_std_deviations*sqrt(e->distance.variance);

#if DEBUG_GAP_WALKER > 2
	    fprintf(GlobalData->gwlogfp,"gap estimate+SLOP = %f, lower bound = %f\n",gap_estimate+SLOP,lower_bound);
#endif

      if ( gap_estimate+SLOP > lower_bound )
	// if sufficient quality, set this edge as first index
	if( eq[j].quality < qualityThresh )
	  {
	    edge_start_index = j;
#if DEBUG_GAP_WALKER > 2
	    if (j != 0) fprintf(GlobalData->gwlogfp, "able to detect destination early for dest: %d\n", destination);
#endif
	  }
      }
 
  }
#endif
  //
  // iterate through all the edges in decreasing order of quality
  //
  for (j = edge_start_index; j < num_edges; j++) {
    //
    // get the edges 
    //

    // look only at the top MAXOUTDEGREE edges
    if( j >= edge_start_index + MAXOUTDEGREE )
      {
#if DEBUG_GAP_WALKER > 2
	fprintf(GlobalData->gwlogfp,"Skipping rest\n");
#endif
	break;
      }

    e = eq[j].edge;
    o = eq[j].orientation;
    assert(e != NULL);

    // continue if the quality is not good enough
    if(  eq[j].quality > qualityThresh )
      continue;

    //
    // filter the trivial path
    // HERE WE HAVE TO ADAPT THE CALL SUCH THAT
    // TRIVIAL WALKS ARE WALKED
    if (((e->idA == source) && (e->idB == destination)) ||
	((e->idB == source) && (e->idA == destination))) {

      source_chunk->flags.bits.walkedTrivial = TRUE;

#if DEBUG_GAP_WALKER > 2
      if (verbose) {
	for (i = 0; i < level; i++) 
	  fprintf(GlobalData->gwlogfp," ");
	fprintf(GlobalData->gwlogfp," (trivial path, edge skipped)\n");
      }
# endif
      continue;
    }

    //
    // get the other end
    //
    if (from_cid == e->idA)
      tid = e->idB;
    else
      tid = e->idA;
    to_chunk = GetGraphNode(ScaffoldGraph->RezGraph, tid);
    assert(to_chunk != NULL);

    // filter some of the containments
    // we ignore containment edges that have a shorter
    // overlap than the chunk length (we allow some slop)

    if (0)  // don't want to do this if walking external bacs
	{
	  if( to_chunk->bpLength.mean+CGW_DP_MINLEN <= -e->distance.mean)
      {
#if DEBUG_GAP_WALKER > 2
		fprintf(GlobalData->gwlogfp,
				"-> IGNORE non contributing containment EDGE (%d,%d), quality %5.2f, orientation %s\n",
				from_cid,
				tid,
				quality(e, from_cid),
				Orientation_As_String(GetEdgeOrientationWRT(e, from_cid)));
#endif
		continue;
      }
	}

#if DEBUG_GAP_WALKER > 2
    if (verbose) {
      for (i = 0; i < level; i++) 
	fprintf(GlobalData->gwlogfp," ");
      fprintf(GlobalData->gwlogfp,
	      "-> chosen EDGE (%d,%d), quality %5.2f, orientation %s\n",
	      from_cid,
	      tid,
	      quality(e, from_cid),
	      Orientation_As_String(GetEdgeOrientationWRT(e, from_cid)));
    }
#endif

    //
    // if we have already been there and the end vertex is not the one
    // we are looking for (except the case of a direct link) ... get
    // the next edge
    //
    assert(subgraph->table[tid] != NULL);
    if ( ( (subgraph->table[tid]->visited) || (tid == source ) )   &&
	 ((! Stop_At_The_Other_End(to_chunk, destination, subgraph)) ))  {
#   if DEBUG_GAP_WALKER > 2
      fprintf(GlobalData->gwlogfp,"Already been here, but not at the other end\n");
#endif
      continue;
    }



#   if DEBUG_GAP_WALKER > 2
    if (verbose) {
      for (i = 0; i < level; i++)
	fprintf(GlobalData->gwlogfp," ");
      fprintf(GlobalData->gwlogfp,"travelled distance mean %5.2f, variance %5.2f\n",
	      travelled_distance.mean,
	      travelled_distance.variance);
      for (i = 0; i < level; i++)
	fprintf(GlobalData->gwlogfp," ");
      fprintf(GlobalData->gwlogfp,"edge mean %5.2f, edge variance %5.2f\n",
	      e->distance.mean,
	      e->distance.variance);
      for (i = 0; i < level; i++) 
	fprintf(GlobalData->gwlogfp," ");
      fprintf(GlobalData->gwlogfp,"chunk mean %5.2f, chunk variance %5.2f\n",
	      to_chunk->bpLength.mean,
	      to_chunk->bpLength.variance);
    }
#   endif


    // add the overlap distance (which is negative for overlaps), 
    // this puts at the trailing edge of the current chunk
    travelled_distance.mean     += e->distance.mean;
    travelled_distance.variance += e->distance.variance;

    assert(travelled_distance.variance >= 0.0);
    
    if ( GlobalData->walkLevel == 1 )
    {
      neg_num_std_deviations = 3.0;
      pos_num_std_deviations = 3.0;
    }
    else 
    {
      neg_num_std_deviations = 5.0;
      pos_num_std_deviations = 5.0;
    }
 

#   if DEBUG_GAP_WALKER > 5
    if (verbose) {
      for (i = 0; i < level; i++) 
	fprintf(GlobalData->gwlogfp,
		" ");
      fprintf(GlobalData->gwlogfp,
	      "adding mean %5.2f, adding variance %5.2f, checked distance? %d\n",
	      add_mean,
	      add_variance,
	      do_check_distance);
    }
#   endif

    // adjust the travelled distance so far, mark the to_chunk visited
    // and save the distance
    // add the chunk length to get ot the end of the current chunk
    travelled_distance.mean     += to_chunk->bpLength.mean;
    travelled_distance.variance += to_chunk->bpLength.variance;

#   if DEBUG_GAP_WALKER > 2
    if (verbose) {
      for (i = 0; i < level; i++)
	fprintf(GlobalData->gwlogfp," ");
      fprintf(GlobalData->gwlogfp,"travelled after adding to_chunk distance mean %5.2f, variance %5.2f\n",
	      travelled_distance.mean,
	      travelled_distance.variance);
    }
#endif

    assert(subgraph->table[tid] != NULL);
    if (subgraph->table[from_cid])
      subgraph->table[tid]->d = subgraph->table[from_cid]->d + quality(e, from_cid);
    else
      subgraph->table[tid]->d = quality(e, from_cid);
    Set_Visited_Bit(subgraph, tid);
    subgraph->table[tid]->distance = travelled_distance;

    //
    // if we have reached the other end or this node will bring us to
    // the end and we have travelled far enough and not too far.
    //
    {
      float64 gap_estimate = travelled_distance.mean-to_chunk->bpLength.mean;
      float64 upper_bound  = max_distance->mean + pos_num_std_deviations * sqrt(max_distance->variance);
      float64 lower_bound  = max_distance->mean - neg_num_std_deviations * sqrt(max_distance->variance);
#   if DEBUG_GAP_WALKER > 2
      fprintf(GlobalData->gwlogfp,"gap estimate+SLOP = %f, upper bound = %f, lower bound = %f\n",gap_estimate+SLOP,upper_bound,lower_bound);
#endif
      
      fprintf(GlobalData->gwlogfp, 
			  "*** terminate(to_chunk, destination, subgraph): %d\n",
			  terminate(to_chunk, destination, subgraph));
      fprintf(GlobalData->gwlogfp, 
			  "*** gap_estimate+SLOP > lower_bound: %d\n",
			  gap_estimate+SLOP > lower_bound);
      fprintf(GlobalData->gwlogfp, 
			  "*** gap_estimate < upper_bound + SLOP: %d\n",
			  gap_estimate < upper_bound + SLOP);
      fprintf(GlobalData->gwlogfp, 
			  "*** SLOP: %d\n",
			  SLOP);


      if( terminate(to_chunk, destination, subgraph) &&( gap_estimate+SLOP <= lower_bound ) )
	{
	  // if we reached the destination but walked too short
	  // we undo the to_chunk
	  source_chunk->flags.bits.walkedTooShort = TRUE;
	  // we also record the closest miss 
	  if( *tooShort > lower_bound-gap_estimate )
	    *tooShort = lower_bound-gap_estimate;

	  travelled_distance.mean     -= to_chunk->bpLength.mean;
	  travelled_distance.variance -= to_chunk->bpLength.variance;
	  travelled_distance.mean     -= e->distance.mean;
	  travelled_distance.variance -= e->distance.variance;
	  continue;

	}

      if( terminate(to_chunk, destination, subgraph) &&( gap_estimate >= upper_bound + SLOP) )
	{
	  source_chunk->flags.bits.walkedTooLong = TRUE;
	  // we also record the closest miss 
	  if( *tooLong > gap_estimate-upper_bound )
	    *tooLong = gap_estimate-upper_bound;	  
	}

      if( terminate(to_chunk, destination, subgraph) &&
	  ( gap_estimate+SLOP > lower_bound ) && 
	  ( gap_estimate < upper_bound + SLOP)) 
	{
	  hops = level + 1;
	} else
	  if( gap_estimate >= upper_bound )
	    {
	      // undo the to_chunk if we walked too far
	      travelled_distance.mean     -= to_chunk->bpLength.mean;
	      travelled_distance.variance -= to_chunk->bpLength.variance;
	      travelled_distance.mean     -= e->distance.mean;
	      travelled_distance.variance -= e->distance.variance;
	      continue;
	    }
	  else
	    hops = Find_Greedy_Path(subgraph, tid, source, destination, max_distance,
				    travelled_distance, level + 1,
				    new_end[o], quality, qualityThresh, terminate, colour, verbose,currentCalls,maxCalls,tooShort,tooLong);
      if (hops) {
#if DEBUG_GAP_WALKER > -1
	if (verbose) {
	  int i;
	  for (i = 0; i < level; i++)
	    fprintf(GlobalData->gwlogfp,
		    " ");
	  fprintf(GlobalData->gwlogfp,
		  "Found path of length %5.2f (stddev %5.2f) with %d hops using edge (%d,%d)\n",
		  travelled_distance.mean,
		  sqrt(travelled_distance.variance),
		  hops,
		  from_cid,
		  tid);
	}
#endif
	//
	// ... mark the path
	//
	Set_Path_Bit(subgraph, tid);
	subgraph->table[tid]->path_id = colour - FIRST_PATH_COLOUR + 1;
	subgraph->table[tid]->path_parent = from_cid;
	subgraph->table[from_cid]->best_edge = e;
	
	assert(eq == eq2); 
	free(eq);
	
	return hops;
      }
      else
	{
#if DEBUG_GAP_WALKER > 0
	  fprintf(GlobalData->gwlogfp,"Got FALSE returned, undoing to_chunk\n");
#endif
	}
    }
    //
    // discard the changes to the distance 
    //
    travelled_distance.mean     -= to_chunk->bpLength.mean;
    travelled_distance.variance -= to_chunk->bpLength.variance;
    travelled_distance.mean     -= e->distance.mean;
    travelled_distance.variance -= e->distance.variance;
  }
  assert(eq == eq2);
  free(eq);
  return FALSE;
}


 void Check_Edge_Distance(chunk_subgraph * subgraph,
			 int (* filter)(CIEdgeT *)) {
  //
  // compute the error of the CI edges by comparing them with the
  // simulator coordinates. It will consider only edges such that
  // filter(e) == TRUE
  //
  CDS_CID_t
    cid,
    tid,
    eid;
  int
    i,
    o,
    /*
    n_a,
    n_b,
    */
    or,
    end = ALL_END,
    chunks = 0,
    bad_edges = 0,
    bad_edges_10 = 0,
    bad_edges_100 = 0,
    bad_edges_1000 = 0,
    tot_edges = 0;
  CIEdgeT
    * e;
  ChunkInstanceT  
    * from,
    * to;
  double
    d,
    delta,
    left,
    right,
    delta_max = 0.0,
    delta_avg = 0.0;

  for (i = 0; i < subgraph->size; i++) {
    cid = subgraph->node[i].cid;
    /*    
    //
    // Art's variation
    //
    for (n_a = 0, o = start_index[A_END]; o < end_index[A_END]; o++)
      n_a += subgraph->table[cid]->num_edges[orient[o]];
    for (n_b = 0, o = start_index[B_END]; o < end_index[B_END]; o++)
      n_b += subgraph->table[cid]->num_edges[orient[o]];

    //fprintf(GlobalData->gwlogfp,
    //	    "() %d %d\n", n_a, n_b);

    if ((n_a == 1) && (n_b == 1))
      end = ALL_END;
    else 
      if (n_a == 1)
	end = A_END;
      else
	if (n_b == 1)
	  end = B_END;
	else continue;
    */
    //
    // get the chunk
    //
    chunks++;
    from = GetGraphNode(ScaffoldGraph->RezGraph, cid);
    assert(from != NULL);

    //
    // scan the edges
    //
    for (or = start_index[end]; or < end_index[end]; or++) {
      o = orient[or];
      for (eid = 0; eid < subgraph->table[cid]->num_edges[o]; eid++) {

	//
	// get the edge
	//
	e = subgraph->table[cid]->edge[o][eid];
	assert(e != NULL);

	//
	// check the edge
	//
	if (! filter(e))
	  continue;

	//
	// get the other node
	//
	if (cid == e->idA)
	  tid = e->idB;
	else
	  tid = e->idA;
	to = GetGraphNode(ScaffoldGraph->RezGraph, tid);
	assert(to != NULL);

	//
	// disregard repeat chunks
	//
	if ((from->flags.bits.cgbType != (unsigned int)UU_CGBTYPE)||
	    (to->flags.bits.cgbType != (unsigned int)UU_CGBTYPE)) {
#         if DEBUG_GAP_WALKER > 2
	  fprintf(GlobalData->gwlogfp,
		  "> Found a repeat chunk %d(%d,%d) or %d(%d,%d)\n",
		  cid,
		  from->aEndCoord,
		  from->bEndCoord,
		  tid,
		  to->aEndCoord,
		  to->bEndCoord);
#         endif
	  continue;
	}

	//
	// compute the actual distance between the
	// chunk <from> and the chunk <to> from
	// the simulator coordinates
	//
	switch (GetEdgeOrientationWRT(e,cid)) {
	case AB_AB :
	  d = (double)(to->aEndCoord - from->bEndCoord);
	  break;
	case AB_BA :
	  d = (double)(to->bEndCoord - from->bEndCoord);
	  break;
	case BA_AB :
	  d = (double)(to->aEndCoord - from->aEndCoord);
	  break;
	case BA_BA :
	  d = (double)(to->bEndCoord - from->aEndCoord);
	  break;
	default :
	  assert(0);
	}

	//
	// adjust the sign
	//
	if (min(from->aEndCoord, from->bEndCoord) > min(to->aEndCoord, to->bEndCoord))
	  d = -d;
	
	//
	// compute the interval
	//
	assert(e->distance.variance >= 0.0);
	left = (double)e->distance.mean - 3.0 * sqrt(e->distance.variance);
	right = (double)e->distance.mean + 3.0 * sqrt(e->distance.variance);
	
	//
	// compute delta (error)
	//
	if (d < left)
	  delta = left - d;
	else if (d > right)
	  delta = d - right;
	else 
	  delta = 0.0;

	if (delta > 0)
	  bad_edges++;

	if (delta > 10)
	  bad_edges_10++;

	if (delta > 100)
	  bad_edges_100++;

	if (delta > 1000)
	  bad_edges_1000++;

	//
	// add up
	//
	delta_avg += delta;
	if (delta > delta_max)
	  delta_max = delta;
	tot_edges++;
      }	
    }
  }

# if DEBUG_GAP_WALKER > 0
  assert(tot_edges != 0);
  fprintf(GlobalData->gwlogfp,
	  "chunks %d, total %d\ncorrect %d, bad %d, total %d\n  %3d  |  %3d  |  %3d  |  %3d  |  %3d  |  %3d  |\nmax_err %5.2f\navg_err %5.2f\n",
	  chunks,
	  subgraph->size,
	  tot_edges - bad_edges,
	  bad_edges,
	  tot_edges,
	  tot_edges - bad_edges,
	  bad_edges,
	  bad_edges_10,
	  bad_edges_100,
	  bad_edges_1000,
	  tot_edges,
	  delta_max,
	  delta_avg/(double)tot_edges);
# endif
}



//
// Path_Degree() returns the number of unique chunks reachable from
// the chunk <cid>
//
static int Path_Degree(chunk_subgraph * s,
		       CDS_CID_t cid,
		       int end,
		       int (* filter)(CIEdgeT *)) {
  //
  // Path_Degree() returns the number of unique chunks reachable from
  // the chunk <cid>
  //
  int
    i,
    c = 0;
  LengthT
    dist = {0.0, 0.0},
    max = {1000000000000000.0, 0.0};
    
  long currentCalls = 0;
  // initialize the distances
  //
  float tooShort = 0.0;
  float tooLong = 0.0;

  for (i = 0; i < s->size; i++)
    s->node[i].d = 0.0;

  //
  // search for the paths
  //
  while (Find_Greedy_Path(s, cid, cid, 0, &max, dist,
			  0, end, No_Quality, NO_QUALITY_THRESH, Stop_At_Unique, c, FALSE,&currentCalls,MAXCALLS,&tooShort,&tooLong)) { 
#   if DEBUG_GAP_WALKER > 0
    fprintf(GlobalData->gwlogfp,"\nNodes on the path (ordered by chunk id):\n");
    Print_Path(s);
#   endif
    Clear_All_Path_Bit(s);          // clear the path bit 
    c++;
  }

  //
  // store the path degree
  //
  return c;
}



// --------------------
// bit and path methods
// --------------------



void Mark_Unique(Scaffold_Fill_t * fill_chunks) {
  //
  // set the isUnique flag for all chunks in the struct <fill_chunks>
  // (not used anymore)
  //
  Gap_Fill_t
    * fc;
  ChunkInstanceT  
    * chunk;
  int
    sid,
    gapid,
    k;

  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++) {
    CIScaffoldT
      * scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
    if ((isDeadCIScaffoldT(scaff)) || (scaff->type != REAL_SCAFFOLD))
      continue;
    for (gapid = 0; gapid < fill_chunks[sid].num_gaps; gapid++) {
      fc = &(fill_chunks[sid].gap[gapid]);
      for (k = 0; k < fc->num_chunks; k++) 
	if ((fc->chunk[k].chunk_id) && (fc->chunk[k].keep)) {
	  chunk = GetGraphNode(ScaffoldGraph->RezGraph, fc->chunk[k].chunk_id);
	  assert(chunk != NULL);
	  chunk->type = UNIQUECHUNK_CGW;
	}
    }
  }
}



void Compute_Outdegree (chunk_subgraph * s) {
  //
  // computes the outdegree of each side of each chunk
  // stores into a celagram file
  //
  int
    i,
    n_a,
    n_b,
    cid,
    k;
  ChunkInstanceT 
    * chunk;
  FILE
    * log1,
    * log2;
  char
    filename1[256],
    filename2[256];
  
  //
  // file handling
  //
  sprintf(filename1,
	  "%s.A.dat",
	  GW_Filename_Prefix);
  sprintf(filename2,
	   "%s.B.dat",
	   GW_Filename_Prefix);
  log1 = file_open(filename1,"w");
  log2 = file_open(filename2,"w");
  assert(log1 != NULL);
  assert(log2 != NULL);
  fprintf(log1,
	  "%s: outdegree A end\n",
	  GW_Filename_Prefix);
  fprintf(log2,
	  "%s: outdegree B end\n",
	  GW_Filename_Prefix);
  
  //
  // iterate all over the chunks
  //
  for (i = 0; i < s->size; i++) {
    cid = s->node[i].cid;
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
    assert(chunk != NULL);
    for (n_a = 0, k = start_index[A_END]; k < end_index[A_END]; k++)
      n_a += s->table[cid]->num_edges[orient[k]];
      fprintf(log1, "%d ", n_a);
      for (n_b = 0, k = start_index[B_END]; k < end_index[B_END]; k++)
	n_b += s->table[cid]->num_edges[orient[k]];
      fprintf(log2, "%d ", n_b);
  }
}



void Compute_Path_Outdegree(chunk_subgraph * s)  {
  //
  // compute the number of uniques reachable from each
  // non-unique in the subgraph <s>
  //
  int
    i,
    a_out,
    a_end,
    b_out,
    b_end,
    color,
    cid;
  FILE
    * log1,
    * log2;
  char
    filename1[STR_LEN],
    filename2[STR_LEN],
    orient[STR_LEN];
  ChunkInstanceT
    * chunk;
# if CREATE_CAM_FILE > 0
  FILE
    * cam;
  char
    cam_filename[STR_LEN],
    * Colour[12] = {
      "CFFFF00 T2 S # uniques",
      "CFF0000 T2 S # zero",
      "CEEEEEE T2 S # one",
      "CDDDDDD T2 S # two",
      "CCCCCCC T2 S # three",
      "CBBBBBB T2 S # four",
      "CAAAAAA T2 S # five",
      "C999999 T2 S # six",
      "C888888 T2 S # seven",
      "C777777 T2 S # eight",
      "C666666 T2 S # nine",
      "C555555 T2 S # ten_or_more"};
# endif
  
# if CREATE_CAM_FILE > 0
  //
  // open the cam file
  //
  sprintf(cam_filename, "%s.po.cam", GW_Filename_Prefix);
  cam = file_open (cam_filename, "w");
  assert(cam != NULL);

  //
  // output the colours
  //
  for (i = 0; i < 12; i ++)
    fprintf(cam, "%d: %s\n", i, Colour[i]);
# endif

  //
  // open the two celagram files
  //
  sprintf(filename1,
	  "%s.a.dat",
	  GW_Filename_Prefix);
  sprintf(filename2,
	  "%s.b.dat",
	  GW_Filename_Prefix);
  log1 = file_open(filename1,"w");
  log2 = file_open(filename2,"w");
  assert(log1 != NULL);
  assert(log2 != NULL);

  //
  // print titles
  //
  fprintf(log1,
	  "%s: path outdegree A end\n",
	  GW_Filename_Prefix);
  fprintf(log2,
	  "%s: path outdegree B end\n",
	  GW_Filename_Prefix);
  
  //
  // compute the path outdegree of each one
  //
  for (i = 0; i < s->size; i++) {  
    cid = s->node[i].cid;
    chunk = GetGraphNode(ScaffoldGraph->RezGraph, cid);
    if  (chunk->aEndCoord <= chunk->bEndCoord) {
      strcpy(orient, "direct");
      a_end = chunk->aEndCoord;
      b_end = chunk->bEndCoord;
    } else {
      a_end = chunk->bEndCoord;
      b_end = chunk->aEndCoord;
      strcpy(orient, "reverse");
    }
    if (! Is_Unique(chunk)) {
      Clear_All_Path_Bit(s);
      Clear_All_Visited_Bit(s);
      a_out = Path_Degree(s, cid, A_END, Is_Not_Guide);
      fprintf (log1, "%d ",
	       a_out);
      Clear_All_Path_Bit(s);
      Clear_All_Visited_Bit(s);
      b_out = Path_Degree(s, cid, B_END, Is_Not_Guide);
      fprintf (log2, "%d ",
	       b_out);
      Clear_All_Path_Bit(s);
      Clear_All_Visited_Bit(s);

      if (max(a_out, b_out) >= 10) {
	color = 11; 
      } else {
	color = max(a_out, b_out) + 1;
      }
#     if CREATE_CAM_FILE > 0
      fprintf(cam,
	      "%d: %d A%d %d # %s chunk, %d A end outdegree, %d B end outdegree\n",
	      cid,
	      a_end,
	      color,
	      b_end,
	      orient,
	      a_out,
	      b_out);
#     endif
    } else {
#    if CREATE_CAM_FILE > 0
     fprintf(cam,
	      "%d: %d A0 %d # %s unique chunk\n",
	      cid,
	      a_end,
	      b_end,
	      orient);
#   endif
    }
  }
# if CREATE_CAM_FILE > 0
  fclose(cam);
# endif
  fclose(log1);
  fclose(log2);
}



int Annotate_Edges_For_Consistent_Path(chunk_subgraph * s, chunk_subgraph * ref) {
  //
  // annotates the edges in s if there is a path that confirm them
  // f is the full subgraph to use a reference (the one with more edges)
  //
  int
    i,
    o,
    k,
    cid,
    hops,
    id;
  static int
    not_confirmed_but_good = 0,
    not_confirmed_and_bad = 0,
    confirmed_but_bad = 0,
    confirmed_and_good = 0,
    total = 0;
  CIEdgeT
    * e;
  LengthT
    dist = {0.0, 0.0};
  long currentCalls=0;
  float tooShort = 0.0;
  float tooLong  = 0.0;

  //
  // make sure that the distance are zero
  //
  for (i = 0; i < s->size; i++) {
    s->node[i].distance.mean = 0.0;
    s->node[i].distance.variance = 0.0;
  }
  
  //
  // now scan all the chunks
  //
  for (i = 0; i < s->size; i++) {
    cid = s->node[i].cid;
    for (k = 0; k < NUM_ORIENTATIONS; k++) {
      o = orient[k];
      for (id = 0; id < s->table[cid]->num_edges[o]; id++) {
	e = s->table[cid]->edge[o][id];
	assert(e != NULL);
	if ((e->idA < e->idB)) { //&&(e->edgesContributing == 1)	
#         if DEBUG_GAP_WALKER > 1
	  double
	    d,
	    left,
	    right,
	    delta;
	  ChunkInstanceT  
	    * from = GetGraphNode(ScaffoldGraph->RezGraph, e->idA),
	    * to = GetGraphNode(ScaffoldGraph->RezGraph, e->idB);
	  assert(from != NULL);
	  assert(to != NULL);
	  
	  switch (GetEdgeOrientationWRT(e, e->idA)) {
	  case AB_AB :
	    d = (double)(to->aEndCoord - from->bEndCoord);
	    break;
	  case AB_BA :
	    d = (double)(to->bEndCoord - from->bEndCoord);
	    break;
	  case BA_AB :
	    d = (double)(to->aEndCoord - from->aEndCoord);
	    break;
	  case BA_BA :
	    d = (double)(to->bEndCoord - from->aEndCoord);
	    break;
	  default :
	    assert(0);
	  }
	  
	  //
	  // adjust the sign
	  //
	  if (min(from->aEndCoord, from->bEndCoord) > min(to->aEndCoord, to->bEndCoord))
	    d = -d;
	  
	  //
	  // compute the interval
	  //
	  assert(e->distance.variance >= 0.0);
	  left = (double)e->distance.mean - 3.0 * sqrt(e->distance.variance);
	  right = (double)e->distance.mean + 3.0 * sqrt(e->distance.variance);
	  
	  //
	  // compute delta (error)
	  //
	  if (d < left)
	    delta = left - d;
	  else if (d > right)
	    delta = d - right;
	  else 
	    delta = 0.0;

	  if ((from->flags.bits.cgbType != (unsigned int)UU_CGBTYPE)||
	      (to->flags.bits.cgbType != (unsigned int)UU_CGBTYPE)) {
	    fprintf(GlobalData->gwlogfp,
		    "One of the two end is a chunk %d(%d,%d) or %d(%d,%d)\n",
		    e->idA,
		    from->aEndCoord,
		    from->bEndCoord,
		    e->idB,
		    to->aEndCoord,
		    to->bEndCoord);
	  }

	  fprintf(GlobalData->gwlogfp,">> Edge %d->%d, weight %d, distance %5.2f\n",
		  e->idA, e->idB,
		  e->edgesContributing,
		  e->distance.mean);
#         endif

	  Clear_All_Visited_Bit(ref);
	  //
	  // try to find a path from idA to idB
	  //
	  hops = Find_Greedy_Path(ref, e->idA, e->idA, e->idB, &(e->distance), dist,
				  0, ALL_END, No_Quality, NO_QUALITY_THRESH, Stop_At_The_Other_End, 0, FALSE,&currentCalls,MAXCALLS,&tooShort,&tooLong);
	  if (hops) {
#           if DEBUG_GAP_WALKER > 1
	    if (delta > 0.0) {
	      fprintf(GlobalData->gwlogfp,"Edge is confirmed by a path of lenght %5.2f, hops %d, but delta %5.2f\n",
		      ref->table[e->idB]->distance.mean,
		      hops,
		      delta);
	      fprintf(GlobalData->gwlogfp,
		      "Chunks %d(%d,%d) -> %d(%d,%d)\n",
		      e->idA,
		      from->aEndCoord,
		      from->bEndCoord,
		      e->idB,
		      to->aEndCoord,
		      to->bEndCoord);
	      confirmed_but_bad++;
	    } else {
	      confirmed_and_good++;
	    }
	    fprintf(GlobalData->gwlogfp,"confirmed: bad/good %d/%d, not_confirmed: good/bad %d/%d, total %d\n\n",
		    confirmed_but_bad,
		    confirmed_and_good,
		    not_confirmed_but_good,
		    not_confirmed_and_bad,
		    total);
#           endif
	    e->flags.bits.hasConfirmingPath = TRUE;
	  } else {
	    Clear_All_Visited_Bit(ref);
	    //
	    // try to find a path from idB to idA
	    //
	    tooShort = tooLong = 0.0;

	    hops = Find_Greedy_Path(ref, e->idB, e->idB, e->idA, &(e->distance), dist,
				    0, ALL_END, No_Quality, NO_QUALITY_THRESH, Stop_At_The_Other_End, 0, FALSE,&currentCalls,MAXCALLS,&tooShort,&tooLong);
	    if (hops) {
#             if DEBUG_GAP_WALKER > 1
	      if (delta > 0.0) {
		fprintf(GlobalData->gwlogfp,"Edge is confirmed by a path of lenght %5.2f, hops %d, but delta %5.2f\n",
			ref->table[e->idB]->distance.mean,
			hops,
			delta);
		confirmed_but_bad++;
	      } else {
		confirmed_and_good++;
	      }
	      fprintf(GlobalData->gwlogfp,
		      "confirmed: bad/good %d/%d, not_confirmed: good/bad %d/%d, total %d\n\n",
		      confirmed_but_bad,
		      confirmed_and_good,
		      not_confirmed_but_good,
		      not_confirmed_and_bad,
		      total);
#             endif
	      e->flags.bits.hasConfirmingPath = TRUE;
	    } else {
#             if DEBUG_GAP_WALKER > 1
	      if (delta == 0.0) {
		fprintf(GlobalData->gwlogfp,
			"Edge is NOT confirmed by a path, but the real distance %5.2f was ok\n",
			d);
		not_confirmed_but_good++;
	      } else {
		not_confirmed_and_bad++;
	      }
	      fprintf(GlobalData->gwlogfp,
		      "confirmed: bad/good %d/%d, not_confirmed: good/bad %d/%d, total %d\n\n",
		      confirmed_but_bad,
		      confirmed_and_good,
		      not_confirmed_but_good,
		      not_confirmed_and_bad,
		      total);
#             endif
	      e->flags.bits.hasConfirmingPath = FALSE;
	    }
	  }
	  total++;
	}
      }
    }
  }
  return 0;
}
