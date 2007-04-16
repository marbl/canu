
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


static char CM_ID[] = "$Id: assemblyStructure.c,v 1.8 2007-04-16 17:36:31 brianwalenz Exp $";


/*********************************************************************/

//#define DEBUG 1
//#define DEBUG_BUCIS 1
//#define DEBUG_MERGE_SCAF 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>


#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "PublicAPI_CNS.h"
#include "MultiAlignStore_CNS.h"

static int numCIs;
static MultiAlignT *ma;
static int cid = 0;

void init_frg_processing(void ){
  numCIs = GetNumGraphNodes(ScaffoldGraph->CIGraph);
  ma = CreateEmptyMultiAlignT();
}

void process_frags_in_UTG(ChunkInstanceT *ci, int ciIsAtoB){
  IntUnitigMesg			ium_mesg;
  int numFrag;
  UnitigStatus   status;

  switch(ci->type){
    case DISCRIMINATORUNIQUECHUNK_CGW:
      status = AS_UNIQUE;
      fprintf(stdout,"* Unitig %d is UNIQUE: DISCRIMINATOR output %d \n",ci->id, cid);
      break;
    case UNIQUECHUNK_CGW:
      status = AS_UNIQUE;
      fprintf(stdout,"* Unitig %d is UNIQUE: output %d \n",ci->id, cid);
      break;
    case UNRESOLVEDCHUNK_CGW:
      if(ci->info.CI.numInstances > 0){
        assert(!ci->flags.bits.isUnique);
        status = AS_SEP;
        fprintf(stdout,"* Unitig %d has %d instances--- output %d SEP\n",ci->id, ci->info.CI.numInstances,cid);
      }else{
        if(ci->scaffoldID != NULLINDEX){
          status = AS_UNIQUE;
          fprintf(stdout,"* Unitig %d has %d instances--- output %d UNIQUE\n",ci->id, ci->info.CI.numInstances,cid);
        }else{
          status = AS_NOTREZ;
          fprintf(stdout,"* Unitig %d has %d instances--- output %d NOTREZ\n",ci->id, ci->info.CI.numInstances,cid);
        }
      }
      break;
    case RESOLVEDREPEATCHUNK_CGW:
      /* SKIP THESE */
      fprintf(stdout,"* Skipping unitig %d --- RESOLVEDREPEAT\n",ci->id);
      return;
    default:
      assert(0);
  }

  ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, ci->id, TRUE);

  numFrag = GetNumIntMultiPoss(ma->f_list);
  assert (ci->type != CONTIG_CGW);

  ium_mesg.coverage_stat = ci->info.CI.coverageStat;
  ium_mesg.status = status;
  ium_mesg.a_branch_point = ci->info.CI.branchPointA;
  ium_mesg.b_branch_point = ci->info.CI.branchPointB;
  ium_mesg.length = GetMultiAlignLength(ma);
  ium_mesg.consensus = Getchar(ma->consensus,0);
  ium_mesg.quality = Getchar(ma->quality,0);
  ium_mesg.forced = 0;
  ium_mesg.num_frags = GetNumIntMultiPoss(ma->f_list);
  ium_mesg.f_list = GetIntMultiPos(ma->f_list,0);

  fprintf(stdout,"     sta:%c cov:%.3f nfr:%d\n",ium_mesg.status,ium_mesg.coverage_stat,ium_mesg.num_frags);
  { int i;
  for (i=0; i < numFrag; ++i){
    fprintf(stdout,"      fragment %d [ %d , %d ] typ:%c\n",
	    ium_mesg.f_list[i].ident,
	    ium_mesg.f_list[i].position.bgn,
	    ium_mesg.f_list[i].position.end,
	    ium_mesg.f_list[i].type
	    );
  }
  }
}


int main( int argc, char *argv[])
{
  int32 restartFromCheckpoint = NULLINDEX;
  Global_CGW *data;
  char *inputPath;
  char *prefix;
  char *outputPath = NULL;
  int setPrefixName = FALSE;
  int setSingleSid = FALSE, singleSid;
  int ckptNum = NULLINDEX;
  int i, index;
  int sid, startingGap, setStartingGap = FALSE;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "c:n:s:")) != EOF)){
      switch(ch) {
        case 'c':
          strcpy( data->File_Name_Prefix, argv[optind - 1]);
          setPrefixName = TRUE;		  
          break;
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case 's':
          singleSid = atoi(argv[optind - 1]);
          setSingleSid = TRUE;
          fprintf( stderr, "setting singleSid to %d\n", singleSid);
          break;
        default :
          errflg++;
      }
    }
  }

  if(data->File_Name_Prefix == NULL || ckptNum == NULLINDEX){
    fprintf(stderr,"Usage: %s -n ckptNum -c prefix [-s singleSid]\n",argv[0]);
    exit(-1);
  }
  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);
  init_frg_processing();

  // over all scfs in graph
  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++){

    CIScaffoldTIterator CIs;
    CIScaffoldT * scaff;
    ContigT *contig;
	
    scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
    assert(scaff != NULL);
    if ((isDeadCIScaffoldT(scaff)) ||
	(scaff->type != REAL_SCAFFOLD) ||
	(scaff->info.Scaffold.numElements < 2))
      continue;

    if ( setSingleSid == TRUE && sid != singleSid) continue;

    fprintf(stdout,"Working on scaffold %d\n",sid);
    // over all ctgs in scf
    InitCIScaffoldTIterator( ScaffoldGraph, scaff, TRUE, FALSE, &CIs);
    while ( (contig = NextCIScaffoldTIterator( &CIs )) != NULL){
      {
	int ctgIsAtoB;
	ContigTIterator UTGs;
	ChunkInstanceT *utg;

	ctgIsAtoB = ( contig->offsetAEnd.mean < contig->offsetBEnd.mean ? 1 : 0 ) ;

	fprintf(stdout,"  Working on contig %d [ %.0f , %.0f ] ori in scf %s\n",contig->id,
		contig->offsetAEnd.mean,contig->offsetBEnd.mean,
		(ctgIsAtoB ? "fwd" : "rev"));

	// over all utgs in contig
	InitContigTIterator(ScaffoldGraph,contig->id,ctgIsAtoB,FALSE,&UTGs);
	while( ( utg = NextContigTIterator( &UTGs) ) != NULL){
	  int utgIsAtoB;
	  utgIsAtoB = (utg->offsetAEnd.mean < utg->offsetBEnd.mean ? 1 : 0 );
	  fprintf(stdout,"    Working on unitig %d [ %.0f , %.0f ] ori in contig %s\n",utg->id,
                  utg->offsetAEnd.mean,utg->offsetBEnd.mean,
		  (utgIsAtoB ? "fwd" : "rev"));
	  process_frags_in_UTG(utg,utgIsAtoB);
	}
      }
    }
  }

  exit(0);
}
