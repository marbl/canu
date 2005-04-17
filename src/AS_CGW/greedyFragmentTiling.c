
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
#include "AS_PER_fragStore.h"
#include "AS_PER_ReadStruct.h"


HashTable_AS *CreateHashTable_AS(int numItemsToHash, HashFn_AS hash, HashCmpFn_AS cmp); /*);*/
int StringHashFn_AS(const void *pointerToString, int length);

static OVL_Store_t  * my_store = NULL;
static OVL_Store_t  * my_second_store = NULL;


void setup_ovlStores(char *OVL_Store_Path){

  assert(my_store==NULL);
  my_store = New_OVL_Store ();
  Open_OVL_Store (my_store, OVL_Store_Path);

  assert(my_second_store==NULL);
  my_second_store = New_OVL_Store ();
  Open_OVL_Store (my_second_store, OVL_Store_Path);

}

void print_olap(Long_Olap_Data_t olap){
        printf ("    %8d %8d %c %5d %5d %4.1f %4.1f\n",
		olap . a_iid,
		olap . b_iid,
		olap . flipped ? 'I' : 'N',
		olap . a_hang, olap . b_hang,
		olap . orig_erate / 10.0, olap . corr_erate / 10.0);
}


typedef enum {
  BEST_MEANS_THICKEST=1,
  BEST_MEANS_LOWEST_ERROR=2,
  BEST_MEANS_BEST_UPPER_BOUND=3
}BestMeasure;



uint64 iid2uid(uint32 iid){
  static uint64 *uid=NULL;
  if(uid==NULL){
    int i;
    int last = getLastElemFragStore (my_frg_store);
    uid = (uint64 *)safe_malloc(sizeof(uint64)*(last+1));
    for(i=0;i<=last;i++){
      uid[i]=0;
    }
  }

  if(uid[iid]!=0){
    return (uid[iid]);
  } else {
    GateKeeperFragmentRecord gkpFrag;
    if(getGateKeeperFragmentStore(my_gkp_store.frgStore,iid,&gkpFrag)!=0)
      assert(0);
    uid[iid] = gkpFrag.readUID;
  }
  return uid[iid];
}


int get_clr_len(uint32 iid){
  static int *len=NULL;
  if(len==NULL){
    int i;
    int last = getLastElemFragStore (my_frg_store);
    len = (int *)safe_malloc(sizeof(int)*(last+1));
    for(i=0;i<=last;i++){
      len[i]=-1;
    }
  }

  if(len[iid]!=-1){
    return (len[iid]);
  } else {
    uint clr_bgn,clr_end;
    getFragStore(my_frg_store,iid,FRAG_S_ALL,fsread);
    getClearRegion_ReadStruct(fsread, &clr_bgn,&clr_end, READSTRUCT_LATEST);
    len[iid]= clr_end-clr_bgn;
  }
  return len[iid];
}



void olap_ref_swap(Long_Olap_Data_t *o){

  /*
                ---->  x
                  ----> y

                  ----> y
                ---->  x


			 
           ---->  x
             <---- y
			

           ----> y
             <---- x  
  */            
  int t=o->a_iid;
  o->a_iid=o->b_iid;
  o->a_iid=t;
  if(o->flipped){
    t=o->a_hang;
    o->a_hang=o->b_hang;
    o->b_hang=t;
  } else {
    o->a_hang*=-1;
    o->b_hang*=-1;
  }
}

void olap_reverse(Long_Olap_Data_t *o){

  /* if we flip the orientation of a fragment then ...

          ---> x
           ---> y

           ---> x
         <--- y


		 or
 
          ---> x
           <--- y

           ---> x
         ---> y

       so "flipped" gets inverted and we also have to let ah' <- -bh and bh' <- -ah
  */
  int t=-o->a_hang;
  o->a_hang=-o->b_hang;
  o->b_hang=t;
  o->flipped=1-o->flipped;
}

double binom(int i,int n,double p){
  //  (n choose i)*p^i*(1-p)*(n-i)
  return exp ( lgamma(n+1.) - lgamma(n-i+1.) - lgamma(i+1.) 
	       + ((double)i)*log(p) + ((double)(n-i))*log(1.-p));
}

double get_95pct_upper(int L, double obsrate){

  // given an observed rate of event, in L observations, we want to find an upper bound on the underlying rate;
  // i.e. a rate that is so high that 95% of the time, we'd expect to observe more events than we really do
  
  //  absolute brute force:

  double lop = 0;
  double tol=.001;
  double hip = 1-tol;
  double delta=1.;
  double try;
  int n = (int) (obsrate * L + .1); // number of observed events

  /*  fprintf(stderr,"Trying to find 95%% confidence upper bound for observed rate of %f out of %d aligned bases ...\n",
      obsrate,L);*/

  do{
    int i; 
    double cum=0;
    try = (hip+lop)/2.;
    for(i=0;i<=n;i++){
      /*      fprintf(stderr,"binom(%d,%d,%f)=%e\n",
	      i,L,try,binom(i,L,try));*/
      cum+=binom(i,L,try);
    }
    if(cum>.05){
      lop=try;
      delta = cum-.05;
    } else {
      hip=try;
      delta = .05-cum;
    }
  } while(delta>tol);

  //  fprintf(stderr,"... decided on %f\n",try);

  return (try);
}  
    
				  

int isDeadEnd(int frg, int offAEnd, double erate, int useCorrected, int skipContaining, int minlen,int frglen){

  OVL_Stream_t  * my_stream = New_OVL_Stream ();
  int retval=0;
  Long_Olap_Data_t  olap;

  assert(my_second_store!=NULL);

  Init_OVL_Stream (my_stream, frg, frg, my_second_store);

  while  (Next_From_OVL_Stream (& olap, my_stream)){

    // exclude too-sloppy overlaps

    if(useCorrected&&olap.corr_erate>erate*1000)continue;

    if(!useCorrected&&olap.orig_erate>erate*1000)continue;

    assert(olap.a_iid == frg);

    // exclude contained overlaps
    if( olap.a_hang > 0 && olap.b_hang < 0)continue;

    // exclude containing overlaps
    if( skipContaining && olap.a_hang < 0 && olap.b_hang > 0)continue;

    if(frglen-max(0,olap.a_hang)+min(0,olap.b_hang) < minlen)continue;

    // if it's off the correct end ...
    if( (offAEnd ? (olap . a_hang < 0) : (olap.b_hang > 0))){
      Free_OVL_Stream (my_stream);
      return 0;
    }
  } 
  Free_OVL_Stream (my_stream);
  return 1;
}



Long_Olap_Data_t better_olap(Long_Olap_Data_t a, Long_Olap_Data_t b, int startingFrg, int offAEnd, BestMeasure bestType, int frglen, 
			     int useCorrected, int skipContaining, int minlen, int avoidDeadEnds, double erate ){

  int aThick, bThick;
  int aMis, bMis;
  int aLen, bLen;

  Long_Olap_Data_t copy_a =a;
  Long_Olap_Data_t copy_b =b;

  if(avoidDeadEnds){ // check for extension off far end of other fragment

    int aIsDead=0;
    int bIsDead=0;

    int aOtherAEnd;  // whether the far end of the other fragment of overlap "a" is its A-end or not
    int bOtherAEnd;  // whether the far end of the other fragment of overlap "b" is its A-end or not

    assert(a.a_iid==startingFrg);
    assert(b.a_iid==startingFrg);

    if(offAEnd){
      aOtherAEnd = !a.flipped;
      bOtherAEnd = !b.flipped;
    } else {
      aOtherAEnd = a.flipped;
      bOtherAEnd = b.flipped;
    }

    aIsDead = isDeadEnd(a.b_iid,aOtherAEnd,erate,useCorrected,skipContaining,minlen,frglen);
    bIsDead = isDeadEnd(b.b_iid,bOtherAEnd,erate,useCorrected,skipContaining,minlen,frglen);
    
    // if status is different, then prefer the live one ...
    if(aIsDead<bIsDead)return a;
    if(aIsDead>bIsDead)return b;
  }

  switch (bestType){

  case BEST_MEANS_LOWEST_ERROR:
    if(useCorrected){ 
      if(a.corr_erate>=b.corr_erate){
	return b;
      } else {
	return a;
      }
    } else {
      if(a.orig_erate>=b.orig_erate){
	return b;
      } else {
	return a;
      }
    }
    break;

  case BEST_MEANS_THICKEST:

    /*  Suppose we want the thickest; considering the fragment of interest to be the reference, and the end we want to extend
	off as the 3' end, then what we want is the minimal ahang, subject to a positive bhang */

    // so, start the canonicalization relative to overlap a;

    // first, if startingFrg not reference, then swap:
    if(copy_a.a_iid!=startingFrg){
      olap_ref_swap(&copy_a);
    }
    // next, if startingFrg orientation wrong, then reverse:
    if(offAEnd){
      olap_reverse(&copy_a);
    }


    // now, repeat for overlap b:

    // first, if startingFrg not reference, then swap:
    if(copy_b.a_iid!=startingFrg){
      olap_ref_swap(&copy_b);
    }
    // next, if startingFrg orientation wrong, then reverse:
    if(offAEnd){
      olap_reverse(&copy_b);
    }


    // now, verify that the bhangs
    assert(copy_a.b_hang >0 );
    assert(copy_b.b_hang >0 );
    if(copy_a.a_hang>=copy_b.a_hang){
      return b;
    } else {
      return a;
    }
    break;

  
  case BEST_MEANS_BEST_UPPER_BOUND:

    // we will take the length of the overlap plus the number of errors observed and compute the
    // 95% confidence interval on the lower bound of the true mismatch rate, favoring whichever
    // gives a lower result.

    {  double upperA,errorsA;
    double upperB,errorsB;
    // compute length of overlap
    // overlap = frglen - max(0,ahang) + min(0,bhang)
    int  olenA = frglen - max(0,a.a_hang) + min(0,a.b_hang);
    int  olenB = frglen - max(0,b.a_hang) + min(0,b.b_hang);

    if(useCorrected){
      upperA = get_95pct_upper(olenA, (a.corr_erate+.9)/1000.);
      upperB = get_95pct_upper(olenB, (b.corr_erate+.9)/1000.);
    } else {
      upperA = get_95pct_upper(olenA, (a.orig_erate+.9)/1000.);
      upperB = get_95pct_upper(olenB, (b.orig_erate+.9)/1000.);
    }
    

    }
    break;

  default:
    fprintf(stderr,"Unknown criterion for best overlap! exit \n");
    exit(-1);
      break;
  }
}


int best_overlap_off_end(int id, int offAEnd,BestMeasure bestType,Long_Olap_Data_t *bestovl,double erate,int useCorrected,int skipContaining,int minlen, int frglen,int avoidDeadEnds){
  Long_Olap_Data_t  olap, bestolap;
  int goodOlap=0;
  OVL_Stream_t  * my_stream = New_OVL_Stream ();
  int retval=0;

  assert(my_store!=NULL);

  Init_OVL_Stream (my_stream, id, id, my_store);

  while  (Next_From_OVL_Stream (& olap, my_stream)){

    // exclude too-sloppy overlaps

    if(useCorrected&&olap.corr_erate>erate*1000)continue;

    if(!useCorrected&&olap.orig_erate>erate*1000)continue;

    assert(olap.a_iid == id);

    // exclude contained overlaps
    if( olap.a_hang > 0 && olap.b_hang < 0)continue;

    // exclude containing overlaps
    if( skipContaining && olap.a_hang < 0 && olap.b_hang > 0)continue;

    // exclude too-short overlaps
    if(frglen-max(0,olap.a_hang)+min(0,olap.b_hang) < minlen)continue;

    // if it's off the correct end ...
    if( (offAEnd ? (olap . a_hang < 0) : (olap.b_hang > 0))){

      //    print_olap(olap);

      if(goodOlap==0){
	goodOlap=1;
	bestolap=olap; 
      } else {
	bestolap=better_olap(bestolap,olap, id, offAEnd,bestType,frglen,useCorrected,skipContaining,minlen,avoidDeadEnds,erate);
      }
    } else {
      // does not stick off the correct end
    }
  }
 
  Free_OVL_Stream (my_stream);

  *bestovl=bestolap;
  if(goodOlap==0)return 0;
  return 1;
}


struct dfs_node_tag {
  int orderReached;
  int finished;
  int parent;
  int height;
  int bestChild;
} dfsGreedyFrg;

void init_dfsGreedyFrg(dfsGreedyFrg *f){
  f->orderReached=-1;
  f->finished=0;
  f->parent=0;
  f->height=0;
  f->bestChild=0;
}

struct ovlfilt_tag {
  BestMeasure type;
  double erate;
  int useCorrected;
  int skipContaining;
  int minlen;
} overlapFilters;



      
       
  

void DFS_longest_chain(dfsGreedyFrg *frgNodes,overlapFilters filter){


/* THE PLAN:

set current = first fragment
while (current <= last fragment )
  if current not already marked finished
    mark current visit time
    init overlap iterator
    set finished = 1
    foreach overlap
      if target not yet visited, 
        set finished = 0 
        set target's parent as current
        set current to target
        break
      if target currently in stack (i.e. visit time valid but not marked finished), continue
      else (i.e. if target is already marked finished)
         if current's height through this one is an improvement
           set height
           set best child
    end
    if finished (i.e. we didn't go down farther)
      mark current finished
      if current has a parent
        set current = current's parent
      else 
        advance current to next unfinished node
end

*/

  dfsGreedyFrg *frgNodes=NULL;
  OVL_Stream_t  * my_stream = New_OVL_Stream ();
  Long_Olap_Data_t olap;

  int curr=0;
  int visited=0;
  int last =  Last_Frag_In_OVL_Store (my_ovl_store);

  frgNodes = (dfsGreedyFrg *) safe_malloc(sizeof(dfsGreedyFrg)*(last+1));
  for(int i=1;i<=last;i++){
    init_dfsGreedyFrg(frgNodes+i);
  }

  //set current = first fragment
  curr=1;

  // while (current <= last fragment )
  while(curr<=last){

    // set deadEnd = 1 -- signal that there were no (non-cycle) overlaps available
    int deadEnd = 1;

    // set finished = 1 -- signal that we did not hit an unvisited target and descend
    int finished = 1

    assert(!frgNodes[curr].finished);

    //  mark current visit time
    frgNodes[startfrg].orderReached = ++visited;

    // init overlap iterator
    Init_OVL_Stream (my_stream, curr, curr, my_ovl_store);

    
    // foreach overlap
    while  (Next_From_OVL_Stream (& olap, my_stream)){

      int target = olap.b_iid;
      assert(olap.a_iid=curr);

      // if target currently in stack (i.e. visit time valid but not marked finished), continue
      if (frgNodes[target].orderReached>0&&frgNodes[target].finished==0){
	continue;
      }

      // if we got here, there is at least one path to follow from current
      deadEnd=0;

      // else if target not yet visited, 
      if(frgNodes[target].visited==0){
	// set finished = 0 
	finished=0;
        // set target's parent as current
	frgNodes[target].parent=curr;
        // set current to target
	curr=target;
        // break out of loop over overlaps -- we are going down to a lower node!
	break;
      }

      // else (i.e. if target is already marked finished)
      {

	int wouldBeHeight = frgNodes[target].height + ???extra_height;

	// if current's height through this one is an improvement
	if(wouldBeHeight > frgNodes[curr].height){
	  // set height
	  frgNodes[curr].height = wouldBeHeight;
	  // set best child
	  frgNodes[curr].bestChild = target;
	}

      }

    }

    // if finished (i.e. we didn't go down farther)
    if(finished){

      if(deadEnd){
	frgNodes[curr].height = get_clr_len(curr);
      }

      // mark current finished
      frgNodes[curr].finished=1;

      // if current has a parent
      if(frgNodes[curr].parent!=0){

        // set current = current's parent
	curr=frgNodes[curr].parent;

      } else {
	
        // advance current to next unfinished node
	while(curr<=last&&frgNodes[curr].finished){
	  curr++;
	}

      }
    }
  }
}

void DFS_longest_chain(){
  /* For each fragment, record
     on the way down:
     - initial order
     - parent
     - depth
     on the way up:
     - favorite child
     - height
     - pop order (or at least pop status)

     At a node, consider all of its children;
     for each child,
       if the child has not been marked finished, continue  // break cycles
       if not already visited,
         order_visited=visisted_nodes++;
         process (recurse)
         mark finished
       if height > highest
          favorite = this child
          highest = height

     We will either be 
       - inefficient, having to loop through the N overlaps of a node in O(N*N) steps ... plus a *lot* of time setting up iterators
       - memory bloated (if we add all overlaps to the stack)
       - tricky (if we let each node know how many of its overlaps we have already processed and we can jump to the i'th overlap)
            but still we have a lot of overlap iterators to set up
  */
}


void finished_with_ovlStore(void){

  Free_OVL_Store (my_store);
  Free_OVL_Store (my_second_store);

}


void usage(char *pgm){
	fprintf (stderr, 
		 "USAGE:  %s -f <FragStoreName> -n <startingFrgNum> -o <full_ovlStore> [-C] [-e <erate cutoff>] [-E] [-m <minlen>] [-Q] [-D]\n"
		 "\t-n startingFrgNum = fragment to walk out from\n"
		 "\t-C specifies to not use containing fragments, i.e. only dovetails\n"
		 "\t-e specifies the maximum mismatch rate (as a fraction, i.e. .01 means one percent)\n"
		 "\t-E specifies that the corrected error rate rather than the original is to be used\n"
		 "\t-Q specifies that the lowest error rate overlap will be used (in place of thickest)\n"
		 "\t-B specifies that the overlap with the lowest upper bound on the mismatch rate will be used\n"
		 "\t-D specifies that best quality should win over avoiding dead ends; by default, a fragment that is not a dead end wins over one that is\n",
		 pgm);
}


int main (int argc , char * argv[] ) {

  Global_CGW *data;
  char *prefix;
  int setFragStore = FALSE;
  int ckptNum = NULLINDEX;
  int i, index;
  char subset_map[1000];
  char full_map[1000];
  char setSubsetMap=0;
  char setFullMap=0;
  char full_ovlPath[1000];
  int setFullOvl=0;
  int bestType = BEST_MEANS_THICKEST;
  int startingFrgNum=-1;
  int stillGoing=0;
  int Aend=-1; 
  int currFrg=-1;
  int firstExtend=-1;

  FragStoreHandle  frag_store;
  int last_stored_frag;
  char *seen;
  int ahang;
  ReadStructp fsread = new_ReadStruct();
  uint frglen,clr_bgn,clr_end;
  double maxError=.3;
  int useCorrectedErate=0;
  int skipContaining=1;
  int minlen=40;
  int avoidDeadEnds=1;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->stderro = stderr;
  data->stderrfp = fopen("findMissedOverlaps.stderr","w");
  data->timefp = stderr;
  data->logfp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,"BCDe:Ef:m:n:o:Q")) != EOF)){
      switch(ch) {
      case 'B':
	bestType=BEST_MEANS_BEST_UPPER_BOUND;
	break;
      case 'C':
	skipContaining = 0;
	break;
      case 'D':
	avoidDeadEnds=0;
	break;
      case 'e':
	maxError = atof(optarg);
	break;
      case 'E':
	useCorrectedErate=1;
	break;
      case 'f':
	strcpy( data->Frag_Store_Name, argv[optind - 1]);
	setFragStore = TRUE;
	break;
      case 'm':
	minlen=atoi(optarg);
	break;
      case 'n':
	startingFrgNum = atoi(argv[optind - 1]);
	break;
      case 'o':
	strcpy(full_ovlPath,argv[optind-1]);
	setFullOvl=1;
	break;
      case 'Q':
	bestType=BEST_MEANS_LOWEST_ERROR;
	break;
      case '?':
	fprintf(stderr,"Unrecognized option -%c",optopt);
      default :
	errflg++;
      }
    }

    if( (setFragStore == 0) || !setFullOvl || startingFrgNum < 0)
      {
	fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setFullOvl = %d startingFrgNum %d\n",
		argc, optind, setFragStore,setFullOvl,startingFrgNum);

	usage(argv[0]);
	exit (-1);
      }
  }
  
  frag_store = openFragStore (data->Frag_Store_Name, "r"); 

  last_stored_frag = getLastElemFragStore (frag_store);

  seen = (char*) malloc((last_stored_frag+1)*sizeof(char));
  assert(seen!=NULL);
  {int i;
   for(i=0;i<=last_stored_frag;i++)seen[i]='\0';
  }

  setup_ovlStores(full_ovlPath);

  stillGoing=1;
  Aend=1;
  currFrg=startingFrgNum;
  seen[currFrg]='\1';
  ahang=0;
  getFragStore(frag_store,currFrg,FRAG_S_ALL,fsread);

  while(stillGoing){
    Long_Olap_Data_t o;
    getClearRegion_ReadStruct(fsread, &clr_bgn,&clr_end, READSTRUCT_LATEST);
    frglen=clr_end-clr_bgn;
    printf("%d  %d %d %s\n",currFrg,ahang,ahang+frglen, Aend?"<--":"-->");
    stillGoing = best_overlap_off_end(currFrg, Aend, bestType,&o,maxError,useCorrectedErate,skipContaining,minlen,frglen,avoidDeadEnds);
    if(stillGoing){
      if(currFrg==startingFrgNum){
	firstExtend= ( o.a_iid==currFrg) ? o.b_iid : o.a_iid;
      }
      if(o.a_iid==currFrg){
	currFrg=o.b_iid;

	if(Aend){
	  ahang+=-o.b_hang;
	} else {
	  ahang+=o.a_hang;
	}

      } else {
	currFrg=o.a_iid;

	if(Aend){
	  if(o.flipped){
	    ahang+=-o.a_hang;
	  } else {
	    ahang+=o.b_hang;
	  }
	} else {
	  if(o.flipped){
	    ahang+=o.b_hang;
	  } else {
	    ahang+=-o.b_hang;
	  }
	}
      }

      if(Aend){
	Aend=1-o.flipped;
      } else {
	Aend=o.flipped;
      }

    }
    if(!seen[currFrg]){
      getFragStore(frag_store,currFrg,FRAG_S_ALL,fsread);
    } else {
      stillGoing=0;
    }
    seen[currFrg]++;
    assert(seen[currFrg]<(char)127);
  }

  if(firstExtend>=0){
    seen[firstExtend]--;
  }


    fprintf(stderr,"now for the other end...\n");

  stillGoing=1;
  Aend=0;
  currFrg=startingFrgNum;
  getFragStore(frag_store,currFrg,FRAG_S_ALL,fsread);
  getClearRegion_ReadStruct(fsread, &clr_bgn,&clr_end, READSTRUCT_LATEST);
  frglen=clr_end-clr_bgn;
  ahang=frglen;
  while(stillGoing){
    Long_Olap_Data_t o;
    if(currFrg!=startingFrgNum&&currFrg!=firstExtend)  printf("%d %d %d %s\n",currFrg,ahang-frglen,ahang,!Aend?"<--":"-->");

    stillGoing = best_overlap_off_end(currFrg, Aend, bestType,&o,maxError,useCorrectedErate,skipContaining,minlen,frglen,avoidDeadEnds);
    if(stillGoing){
      if(o.a_iid==currFrg){
	currFrg=o.b_iid;

	if(Aend){
	  ahang-=-o.b_hang;
	} else {
	  ahang-=o.a_hang;
	}

      } else {
	currFrg=o.a_iid;

	if(Aend){
	  if(o.flipped){

	    /*
                    ---> a_iid
	         <--- b_iid
                 |-| = ahang

	    */
	    ahang-=o.a_hang;
	  } else {
	    /*
                    ---> a_iid
	              ----> b_iid
                        |-| = bhang

	    */
	    ahang-=o.b_hang;
	  }
	} else {
	  if(o.flipped){
	    /*
                    ---> a_iid
	              <--- b_iid
                        |-| = bhang

	    */
	    ahang-=o.b_hang;
	  } else {
	    /*
                    ---> a_iid
	         ---> b_iid
                     |-| = -bhang

	    */
	    ahang-=-o.b_hang;
	  }
	}
      }
      if(Aend){
	Aend=1-o.flipped;
      } else {
	Aend=o.flipped;
      }
    }
    if(!seen[currFrg]){
      seen[currFrg]='\1';
      getFragStore(frag_store,currFrg,FRAG_S_ALL,fsread);
      getClearRegion_ReadStruct(fsread, &clr_bgn,&clr_end, READSTRUCT_LATEST);
      frglen=clr_end-clr_bgn;
    } else {
      stillGoing=0;
    }
  }

  finished_with_ovlStore();  
}

