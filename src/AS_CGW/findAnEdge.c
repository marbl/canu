
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

#include "cds.h"
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "FbacREZ.h"
#include "PublicAPI_CNS.h"
#include "AS_ALN_aligners.h"
#include "AS_ALN_forcns.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_Hash.h"
#include "OlapStoreOVL.h"
#include "Instrument_CGW.h"

void usage(char *pgm){
  fprintf (stderr, "USAGE:  %s -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum>\n",
           pgm);
}

int main (int argc , char * argv[] ) {

  Global_CGW *data;
  char *prefix;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setOvlStore = FALSE;
  int setPrefixName = FALSE;
  int ckptNum = NULLINDEX;
  int i, index;
  char subset_map[1000];
  char full_map[1000];
  char setSubsetMap=0;
  char setFullMap=0;
  char ovlPath[1000];
  int setFullOvl=0;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,"c:f:g:n:")) != EOF)){
      switch(ch) {
        case 'c':
          strcpy( data->File_Name_Prefix, argv[optind - 1]);
          setPrefixName = TRUE;		  
          break;
        case 'f':
          strcpy( data->Frag_Store_Name, argv[optind - 1]);
          setFragStore = TRUE;
          break;
        case 'g':
          strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
          setGatekeeperStore = TRUE;
          break;	  
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0) ){
      fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d\n",
              argc, optind, setFragStore,setGatekeeperStore);

      usage(argv[0]);
      exit (-1);
    }
  }

  ScaffoldGraph = 
    LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);


  fprintf(stdout,"Examining RezGraph ...\n");
  {
    int32 cid;
    ChunkInstanceT *thisCI;
    GraphNodeIterator nodes;

    InitGraphNodeIterator(&nodes, ScaffoldGraph->RezGraph, GRAPH_NODE_DEFAULT);
    while(thisCI = NextGraphNodeIterator(&nodes)){
      GraphEdgeIterator edges;
      CIEdgeT *edge;
      
      InitGraphEdgeIterator(ScaffoldGraph->RezGraph, thisCI->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
      while((edge = NextGraphEdgeIterator(&edges)) != NULL){
	if((edge->idA == 707 && edge->idB == 1303 )||
	   (edge->idB == 707 && edge->idA == 1303 )){
	  PrintGraphEdge(stdout,ScaffoldGraph->RezGraph,"",edge,edge->idA);
	}
      }
    }
  }
  fprintf(stdout,"Examining ScaffoldGraph ...\n");
  {
    int32 cid;
    ChunkInstanceT *thisCI;
    GraphNodeIterator nodes;

    InitGraphNodeIterator(&nodes, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
    while(thisCI = NextGraphNodeIterator(&nodes)){
      GraphEdgeIterator edges;
      CIEdgeT *edge;
      
      InitGraphEdgeIterator(ScaffoldGraph->ScaffoldGraph, thisCI->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
      while((edge = NextGraphEdgeIterator(&edges)) != NULL){
	if((edge->idA == 707 && edge->idB == 1303 )||
	   (edge->idB == 707 && edge->idA == 1303 )){
	  PrintGraphEdge(stdout,ScaffoldGraph->ScaffoldGraph,"",edge,edge->idA);
	}
      }
    }
  }
  fprintf(stdout,"Examining ContigGraph ...\n");
  {
    int32 cid;
    ChunkInstanceT *thisCI;
    GraphNodeIterator nodes;

    InitGraphNodeIterator(&nodes, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
    while(thisCI = NextGraphNodeIterator(&nodes)){
      GraphEdgeIterator edges;
      CIEdgeT *edge;
      
      InitGraphEdgeIterator(ScaffoldGraph->ContigGraph, thisCI->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);// Use merged edges
      while((edge = NextGraphEdgeIterator(&edges)) != NULL){
	if((edge->idA == 707 && edge->idB == 1303 )||
	   (edge->idB == 707 && edge->idA == 1303 )){
	  PrintGraphEdge(stdout,ScaffoldGraph->ContigGraph,"",edge,edge->idA);
	}
      }
    }
  }
  return 0;
}

