
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


//#define VERBOSE

HashTable_AS *CreateHashTable_AS(int numItemsToHash, HashFn_AS hash, HashCmpFn_AS cmp); /*);*/
int StringHashFn_AS(const void *pointerToString, int length);

static HashTable_AS *subset_name2iid;
static HashTable_AS *full_name2iid;
static char subset_names[2000000][20];
static char full_names[2000000][20];
static VA_TYPE(uint64) *subset_iid2name;
static VA_TYPE(uint64) *full_iid2name;
static OVL_Store_t  * my_store = NULL;

typedef enum { SCAFFOLD_BEGIN,SCAFFOLD_END } scfEnds;

int strhashcmp(const void *a , const void *b){
  char *A = (char *) a;
  char *B = (char *) b;
  return strcmp(A,B);
}

static int useFrgMaps = 1;

void setup_frg_maps(const char *subset_map, const char *full_map){
  FILE *subset_frg_map = fopen(subset_map,"r");
  FILE *full_frg_map = fopen(full_map,"r");
  int iid,uid,i=0,j=0;
  char name[100];

  assert(subset_frg_map !=NULL);
  assert(full_frg_map !=NULL);

  subset_iid2name = CreateVA_uint64(2000000);
  full_iid2name = CreateVA_uint64(2000000);

  subset_name2iid = CreateHashTable_AS(2000000,StringHashFn_AS,strhashcmp);
  full_name2iid = CreateHashTable_AS(2000000,StringHashFn_AS,strhashcmp);

  while(fscanf(subset_frg_map,"%d %d %s\n",&iid,&uid,name)==3){

    strcpy(&(subset_names[i][0]),name);

    assert(LookupInHashTable_AS(subset_name2iid,name,strlen(name))==NULL);

    // N.B.: the hash table actually stores the IID, instead of a pointer; but it will have to be cast

    InsertInHashTable_AS(subset_name2iid,&(subset_names[i][0]),strlen(name),(void*)iid);

    // demonstration of how to get value back from hash table:
    //    fprintf(stdout,"post-insert Hash table contains %d\n",
    //            (int) LookupInHashTable_AS(subset_name2iid,name,strlen(name)));

    { 
      uint64 dummy = (uint64) &(subset_names[i][0]);
      SetElement_VA(subset_iid2name,iid, &dummy);
    }
    // demonstration of how to get value back from variable array
    //    fprintf(stdout,"post-insert iid2name gives %s\n",
    //            (char*)(* (uint64*)GetElement_VA(subset_iid2name,iid)));

    i++;
  }

  while(fscanf(full_frg_map,"%d %d %s\n",&iid,&uid,name)==3){
    strcpy(&(full_names[j][0]),name);
    //    fprintf(stdout,"Working on %d %d %s\n",iid,uid,subset_names[j]);
    assert(LookupInHashTable_AS(full_name2iid,name,strlen(name))==NULL);
    // N.B.: the hash table actually stores the IID, instead of a pointer; but it will have to be cast
    InsertInHashTable_AS(full_name2iid,&(full_names[j][0]),strlen(name),(void*)iid);
    SetElement_VA(full_iid2name,iid, &(full_names[j][0]));
    j++;
  }
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
      if  (olap . a_hang < 0){
	//	print_olap(olap);
	retval++;
      }
    } else {
      if  (olap . b_hang > 0){
	//	print_olap(olap);
	retval++;
      }
     }
  }
 
  Free_OVL_Stream (my_stream);
  return(retval);
}


void finished_with_ovlStore(void){

  Free_OVL_Store (my_store);

}


void usage(char *pgm){
	fprintf (stderr, "USAGE:  %s -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum> -1 <subset_map> [ -2 <full_map> ] -o <full_ovlStore>\n"
		 "(set subset_map to NULL and omit full_map if using single asm)\n",
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
#ifdef VERBOSE
    fprintf(stdout,"    Processing unitig %d [ %.0f , %.0f ] ori in contig %s\n",utg->id,
	    utg->offsetAEnd.mean,utg->offsetBEnd.mean,
	    (utg->offsetAEnd.mean < utg->offsetBEnd.mean ? "fwd" : "rev"));
#endif
  }
#ifdef VERBOSE
  fprintf(stdout,"    First utg %d (ori %s), last %d (ori %s)\n",
	  *firstUTGid,
	  *firstUTGisAtoB ? "fwd" : "rev",
	  *lastUTGid,
	  *lastUTGisAtoB ? "fwd" : "rev");
#endif
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
     MultiAlignT *newma = MergeMultiAlignsFast_new(ScaffoldGraph->sequenceDB, NULLFRAGSTOREHANDLE, positions, 0, 1, NULL);

     int nfr = GetNumIntMultiPoss(newma->f_list);
     int len = GetMultiAlignLength(newma);
     *frontCnt=0;
     *tailCnt=0;

     //     fprintf(stdout,"Numfrgs in contig %d , len = %d\n",nfr,len);
     {int i;
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
       if(beg<40){
	 (*frontCnt)++;
	 //	 fprintf(stdout,"    front incremented to %d\n",*frontCnt);
       }
       if(len-end<40){
	 (*tailCnt)++;
	 //	 fprintf(stdout,"    tail incremented to %d\n",*tailCnt);
       }
     }}
   } 
   DeleteVA_IntElementPos(positions);
}

void explore_end_of_contig(ContigT *contig,int whichEnd){
  int ctgIsAtoB;
  int firstUTGid;
  int lastUTGid;
  int firstUTGisAtoB;
  int lastUTGisAtoB;
  int ctgID = contig->id;
  int firstFrag;
  int firstFragisAtoB;
  int firstFrag_full;
  int lastFrag;
  int lastFragisAtoB;
  int lastFrag_full;
  int contigIsAtoB = ( contig->offsetAEnd.mean < contig->offsetBEnd.mean ? 1 : 0 ) ;
  char *name;
  int ctgAendCnt,ctgBendCnt;

#ifdef VERBOSE
  fprintf(stdout,
	  "  Working on contig %d [ %.0f , %.0f ] ori in scf %s\n",contig->id,
	  contig->offsetAEnd.mean,contig->offsetBEnd.mean,
	  (contigIsAtoB ? "fwd" : "rev"));
#endif

  find_first_and_last_unitigs(ctgID,contigIsAtoB,
			      &firstUTGid,
			      &lastUTGid,
			      &firstUTGisAtoB,
			      &lastUTGisAtoB);

  if(! contigIsAtoB ){
    int tmp;
    tmp=firstUTGid;
    firstUTGid=lastUTGid;
    lastUTGid=tmp;
    tmp=firstUTGisAtoB;
    firstUTGisAtoB = ( 1 - lastUTGisAtoB );
    lastUTGisAtoB = ( 1 - tmp );
	  
#ifdef VERBOSE
    fprintf(stdout,
	    "   First utg %d (ori %s), last %d (ori %s)  (after adjusting for contig ori)\n",
	    firstUTGid,
	    firstUTGisAtoB ? "fwd" : "rev",
	    lastUTGid,
	    lastUTGisAtoB ? "fwd" : "rev");
#endif
  }

  
  explore_ending_of_contig(contig,&ctgAendCnt,&ctgBendCnt);  

  if(whichEnd == 0){

    // at this point, firstUTGisAtoB is relative to scaffold as a whole;
    // if yes, then we want the first frag of the unitig; else, we want
    // the last frag of the unitig

    getTipFrag(firstUTGid,&firstFrag,&firstFragisAtoB,firstUTGisAtoB);
    
    // fragment orientation: if !firstUTGisAtoB, then invert the fragment

    if(!firstUTGisAtoB){
      firstFragisAtoB = 1 - firstFragisAtoB;
    }

    if(!useFrgMaps){
      firstFrag_full = firstFrag;
    } else {
      name = (char *) (* (uint64 *)  GetElement_VA(subset_iid2name,firstFrag));
      firstFrag_full = (int) LookupInHashTable_AS(full_name2iid,name,strlen(name));
    }
#ifdef VERBOSE
    fprintf(stdout,"First frg of contig is %d ori %s\n",firstFrag,
	    firstFragisAtoB ? "fwd" : "rev");
#endif
    fprintf(stdout,"  overlaps off scaffold BEGIN (contig %d %c end, frg %d %c end): %d\n",
	    contig->id,contigIsAtoB ? 'A' : 'B',
	    firstFrag_full,firstFragisAtoB ? 'A' : 'B',
	    count_overlaps_off_end(firstFrag_full,firstFragisAtoB));
    fprintf(stdout,"    fragments ending near end: %d\n",
	    contigIsAtoB ? ctgAendCnt : ctgBendCnt);


  } else {

    // at this point, firstUTGisAtoB is relative to scaffold as a whole;
    // if yes, then we want the last frag of the unitig; else, we want
    // the first frag of the unitig

    getTipFrag(lastUTGid,&lastFrag,&lastFragisAtoB,1-lastUTGisAtoB);

    // fragment orientation: if !lastUTGisAtoB, then invert the fragment

    if(!lastUTGisAtoB){
      lastFragisAtoB = 1 - lastFragisAtoB;
    }

    if(!useFrgMaps){
      lastFrag_full = lastFrag;
    } else {
      name = (char *) (* (uint64 *)  GetElement_VA(subset_iid2name,lastFrag));
      lastFrag_full = (int) LookupInHashTable_AS(full_name2iid,name,strlen(name));
    }
    
#ifdef VERBOSE
    fprintf(stdout,"Last frg of contig is %d ori %s\n",lastFrag,
	  lastFragisAtoB ? "fwd" : "rev");
#endif
    fprintf(stdout,"  overlaps off scaffold END (contig %d %c end frag %d %c end): %d\n",
	    contig->id,contigIsAtoB ? 'B' : 'A',
	    lastFrag_full,lastFragisAtoB ? 'B' : 'A',
	    count_overlaps_off_end(lastFrag_full,1-lastFragisAtoB));
    fprintf(stdout,"    fragments ending near end: %d\n",
	    contigIsAtoB ? ctgBendCnt : ctgAendCnt);

  }

}

void explore_ends_of_a_scaffold(int sid){
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

  {
    int ori,ctgAendCnt, ctgBendCnt;
    float begThis,varThis;
    ContigT *lastCtg = NULL;
    contig = NextCIScaffoldTIterator(&CIs);
    
    if(contig->offsetAEnd.mean<contig->offsetBEnd.mean){
      begThis=contig->offsetAEnd.mean;
      varThis=contig->offsetAEnd.variance;
      ori=A_B;
    } else {
      begThis=contig->offsetBEnd.mean;
      varThis=contig->offsetBEnd.variance;
      ori=B_A;
    }
    explore_end_of_contig(contig,SCAFFOLD_BEGIN);

    lastCtg=contig;
    while ( (contig = NextCIScaffoldTIterator( &CIs )) != NULL){
      lastCtg=contig;
    }
    contig=lastCtg;
    if(contig->offsetAEnd.mean<contig->offsetBEnd.mean){
      begThis=contig->offsetAEnd.mean;
      varThis=contig->offsetAEnd.variance;
      ori=A_B;
    } else {
      begThis=contig->offsetBEnd.mean;
      varThis=contig->offsetBEnd.variance;
      ori=B_A;
    }
    explore_end_of_contig(contig,SCAFFOLD_END);  
  }
}

int main (int argc , char * argv[] ) {

  Global_CGW *data;
  char *prefix;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
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
  data->stderro = stderr;
  data->stderrfp = stderr;
  data->timefp = stderr;
  data->logfp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,"c:f:g:n:1:2:o:")) != EOF)){
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
      case 'o':
	strcpy(ovlPath,argv[optind-1]);
	setFullOvl=1;
	break;
      case '1':
	strcpy(subset_map,argv[optind-1]);
	setSubsetMap=1;
	break;
      case '2':
	strcpy(full_map,argv[optind-1]);
	setFullMap=1;
	break;
      case '?':
	fprintf(stderr,"Unrecognized option -%c",optopt);
      default :
	errflg++;
      }
    }

    if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0) || !setSubsetMap || (!setFullMap&& strcmp(subset_map,"NULL") != 0) || !setFullOvl)
      {
	fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d\n",
		argc, optind, setFragStore,setGatekeeperStore);

	usage(argv[0]);
	exit (-1);
      }
  }

  ScaffoldGraph = 
    LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);

  if(strcmp(subset_map,"NULL") != 0){
    setup_frg_maps(subset_map, full_map);
  } else {
    useFrgMaps=0;
  }

  setup_ovlStore(ovlPath);

  init_frg_processing();

  // over all scfs in graph
  { int sid;
  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++){
    explore_ends_of_a_scaffold(sid);
  }
  }
  finished_with_ovlStore();  
}

