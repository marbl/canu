
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


int DistFromGap = 20;

HashTable_AS *CreateHashTable_AS(int numItemsToHash, HashFn_AS hash, HashCmpFn_AS cmp); /*);*/
int StringHashFn_AS(const void *pointerToString, int length);

static OVL_Store_t  * my_store = NULL;
static char verbose = 0;

int strhashcmp(const void *a , const void *b){
  char *A = (char *) a;
  char *B = (char *) b;
  return strcmp(A,B);
}

void setup_ovlStore(char *OVL_Store_Path){

  assert(my_store==NULL);
  my_store = New_OVL_Store ();
  Open_OVL_Store (my_store, OVL_Store_Path);

}

void print_olap(Long_Olap_Data_t olap){
  printf ("    %8d %8d %c %5d %5d %4.1f %4.1f\n",
          olap . a_iid,
          olap . b_iid,
          olap . flipped ? 'I' : 'N',
          olap . a_hang, olap . b_hang,
          olap . orig_erate / 10.0, olap . corr_erate / 10.0);
}

int count_overlaps_off_end(int id, int offAEnd){

  Long_Olap_Data_t  olap;
  OVL_Stream_t  * my_stream = New_OVL_Stream ();
  int retval=0;
  assert(my_store!=NULL);

  Init_OVL_Stream (my_stream, id, id, my_store);

  while  (Next_From_OVL_Stream (& olap, my_stream)){
    //    print_olap(olap);
    if(offAEnd){
      if  (olap . a_hang < 0)
	retval++;
    } else {
      if  (olap . b_hang > 0)
	retval++;
    }
  }
 
  Free_OVL_Stream (my_stream);
  return(retval);
}


void finished_with_ovlStore(void){

  Free_OVL_Store (my_store);

}


void usage(char *pgm){
  fprintf (stderr, "USAGE:  %s -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum> -d distance_from_end -o <full_ovlStore>\n",
           pgm);
}

void find_first_and_last_unitigs(
                                 int ctgID,
                                 int ctgIsAtoB,
                                 int *firstUTGid, 
                                 int *lastUTGid, 
                                 int *firstUTGisAtoB, 
                                 int *lastUTGisAtoB)
{
  ContigTIterator UTGs;
  ChunkInstanceT *utg;
  int beginpos = 999999999;
  int endingpos = -1;
	  
  // over all utgs in contig
  InitContigTIterator(ScaffoldGraph,ctgID,ctgIsAtoB,FALSE,&UTGs);
  while( ( utg = NextContigTIterator( &UTGs) ) != NULL){
    int begpos, endpos;
    int utgIsAtoB = (utg->offsetAEnd.mean < utg->offsetBEnd.mean ? 1 : 0);
    begpos = (utgIsAtoB ? utg->offsetAEnd.mean : utg->offsetBEnd.mean);
    endpos = (utgIsAtoB ? utg->offsetBEnd.mean : utg->offsetAEnd.mean);
    if(begpos<beginpos){
      beginpos=begpos;
      *firstUTGid = utg->id;
      *firstUTGisAtoB = utgIsAtoB;
    }
    if(endpos>endingpos){
      endingpos=endpos;
      *lastUTGid = utg->id;
      *lastUTGisAtoB = utgIsAtoB;
    }
    if (verbose) 
      fprintf(stdout,"    Processing unitig %d [ %.0f , %.0f ] ori in contig %s\n",utg->id,
              utg->offsetAEnd.mean,utg->offsetBEnd.mean,
              (utg->offsetAEnd.mean < utg->offsetBEnd.mean ? "fwd" : "rev"));
  }
  if (verbose) 
    fprintf(stdout,"    First utg %d (ori %s), last %d (ori %s)\n",
            *firstUTGid,
            *firstUTGisAtoB ? "fwd" : "rev",
            *lastUTGid,
            *lastUTGisAtoB ? "fwd" : "rev");
}


static int numCIs;
static MultiAlignT *ma;
static int cid = 0;

void init_frg_processing(void){
  numCIs = GetNumGraphNodes(ScaffoldGraph->CIGraph);
  ma = CreateEmptyMultiAlignT();
}

void getTipFrag(int utgID,int *theFrg, int *theFrgIsAtoB, int wantAEnd){
  int numFrag;
  IntMultiPos *imp;
  ChunkInstanceT *ci = GetGraphNode(ScaffoldGraph->CIGraph, utgID);

  if(ci->flags.bits.isSurrogate){
    ci = GetGraphNode(ScaffoldGraph->CIGraph, ci->info.CI.baseID);
    fprintf(stdout,"Replaced surrogate utg %d with baseCI %d\n",
	    utgID,ci->id);
  }


  ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, ci->id, TRUE);

  numFrag = GetNumIntMultiPoss(ma->f_list);

  assert (ci->type != CONTIG_CGW);

  imp = GetIntMultiPos(ma->f_list,0);
  if(wantAEnd){
    int best=999999999;
    int i;
    for (i=0; i < numFrag; ++i){
      int frgIsAtoB = ( imp[i].position.bgn < imp[i].position.end ?
			1 : 0 );
      int beg = ( frgIsAtoB ? 
		  imp[i].position.bgn : imp[i].position.end);

      if(beg<best){
	best=beg;
	*theFrg = imp[i].ident;
	*theFrgIsAtoB = frgIsAtoB;
      }
    }
  } else {
    int best=-1;
    int i;
    for (i=0; i < numFrag; ++i){
      int frgIsAtoB = ( imp[i].position.bgn < imp[i].position.end ?
			1 : 0 );
      int tail = ( frgIsAtoB ? 
                   imp[i].position.end : imp[i].position.bgn);

      if(best<tail){
	best=tail;
	*theFrg = imp[i].ident;
	*theFrgIsAtoB = frgIsAtoB;
      }
    }
  }
}

void explore_ending_of_contig(ContigT *contig, int *frontCnt, int *tailCnt){

  VA_TYPE(IntElementPos) *positions = CreateVA_IntElementPos(5000);
  MultiAlignT *ma =  LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
  int num_tigs = GetNumIntUnitigPoss(ma->u_list);
  int i;
  for (i=0;i<num_tigs;i++) {
    IntUnitigPos *upos = GetIntUnitigPos(ma->u_list,i);
    IntElementPos pos;
    pos.type = AS_UNITIG;
    pos.ident = upos->ident;
    pos.position.bgn = upos->position.bgn;
    pos.position.end = upos->position.end;
    SetVA_IntElementPos(positions,i,&pos);
  }
   
  {
    MultiAlignT *newma = MergeMultiAlignsFast_new(ScaffoldGraph->sequenceDB, NULL, positions, 0, 1, NULL, NULL);

    int nfr = GetNumIntMultiPoss(newma->f_list);
    int len = GetMultiAlignLength(newma);
    int i;
    *frontCnt=0;
    *tailCnt=0;

    //     fprintf(stdout,"Numfrgs in contig %d , len = %d\n",nfr,len);

    for(i=0;i<nfr;i++){
      IntMultiPos *mpos = GetIntMultiPos(newma->f_list,i);
      int beg = mpos->position.bgn;
      int end = mpos->position.end;
      if(end< beg){
        int tmp=end;
        end=beg;
        beg=tmp;
      }
      //       fprintf(stdout,"  frg pos %d %d\n",beg,end);
      if(beg<DistFromGap){
        (*frontCnt)++;
        //	 fprintf(stdout,"    front incremented to %d\n",*frontCnt);
      }
      if(len-end<DistFromGap){
        (*tailCnt)++;
        //	 fprintf(stdout,"    tail incremented to %d\n",*tailCnt);
      }
    }
  } 
  DeleteVA_IntElementPos(positions);
}

void explore_gap_structure_for_a_scaffold(int sid){
  CIScaffoldTIterator CIs;
  CIScaffoldT * scaff;
  ContigT *contig;
  float endPrev;
  float varPrev;
	
  scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);
  assert(scaff != NULL);
  if ((isDeadCIScaffoldT(scaff)) ||
      (scaff->type != REAL_SCAFFOLD))
    return;

  fprintf(stdout,"Working on scaffold %d\n",sid);
  // over all ctgs in scf
  InitCIScaffoldTIterator( ScaffoldGraph, scaff, TRUE, FALSE, &CIs);
  endPrev=-1;
  while ( (contig = NextCIScaffoldTIterator( &CIs )) != NULL){
    float begThis,varThis;
    if(contig->offsetAEnd.mean<contig->offsetBEnd.mean){
      begThis=contig->offsetAEnd.mean;
      varThis=contig->offsetAEnd.variance;
    } else {
      begThis=contig->offsetBEnd.mean;
      varThis=contig->offsetBEnd.variance;
    }
      
    if(endPrev!=-1){
      fprintf(stdout,"              gap size %f , %f\n",
              begThis-endPrev,
              sqrt(varThis-varPrev));
    }
		
    { 
      int frontCnt=0;
      int tailCnt=0;
      explore_ending_of_contig(contig,&frontCnt,&tailCnt);
      fprintf(stdout,"      %d  <-- frgs within %dbp of end --> %d\n",
              frontCnt,DistFromGap,tailCnt);
    }

    if(contig->offsetAEnd.mean<contig->offsetBEnd.mean){
      endPrev=contig->offsetBEnd.mean;
      varPrev=contig->offsetBEnd.variance;
    } else {
      endPrev=contig->offsetAEnd.mean;
      varPrev=contig->offsetAEnd.variance;
    }

  }
}

int main (int argc , char * argv[] ) {

  Global_CGW *data;
  char *prefix;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int ckptNum = NULLINDEX;
  int i, index;
  char full_ovlPath[1000];
  int setFullOvl=0;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,"c:d:f:g:n:o:v")) != EOF)){
      switch(ch) {
        case 'c':
          strcpy( data->File_Name_Prefix, argv[optind - 1]);
          setPrefixName = TRUE;		  
          break;
        case 'd':
          if (strlen(argv[optind - 1]) > 4) {
            fprintf(stderr,"-d option too large.\n");
            exit(1);
          }
          DistFromGap = atoi(argv[optind - 1]);
          break;
        case 'g':
          strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
          setGatekeeperStore = TRUE;
          break;	  
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case 'v':
          verbose = 1; break;
        case 'o':
          strcpy(full_ovlPath,argv[optind-1]);
          setFullOvl=1;
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if((setPrefixName == FALSE) || (setGatekeeperStore == 0) || !setFullOvl)
      {
	fprintf(stderr,"* argc = %d optind = %d setGatekeeperStore = %d\n",
		argc, optind, setGatekeeperStore);

	usage(argv[0]);
	exit (-1);
      }
  }

  ScaffoldGraph = 
    LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);

  setup_ovlStore(full_ovlPath);

  init_frg_processing();

  // over all scfs in graph
  { int sid;
  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++){
    explore_gap_structure_for_a_scaffold(sid);
  }}
  finished_with_ovlStore();  
}

