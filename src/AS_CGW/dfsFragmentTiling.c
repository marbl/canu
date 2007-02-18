
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


#include "AS_global.h"
#include "OlapStoreOVL.h"
#include "AS_PER_gkpStore.h"
#include "UtilsREZ.h"


#undef DEBUG_DFS
#undef GRAPHIC_DEBUG_DFS
//#define GRAPHIC_DEBUG_DFS 0

static  GateKeeperStore *my_gkp_store ;
static  OVL_Store_t *my_ovl_store = NULL;
static ReadStructp fsread;



void setup_stores(char *OVL_Store_Path, char *Gkp_Store_Path){

  assert(my_ovl_store==NULL);

  assert(OVL_Store_Path!=NULL);
  assert(Gkp_Store_Path!=NULL);

  my_ovl_store = New_OVL_Store ();
  Open_OVL_Store (my_ovl_store, OVL_Store_Path);

  my_gkp_store = openGateKeeperStore( Gkp_Store_Path, "r");
}

void finished_with_stores(void){

  Free_OVL_Store (my_ovl_store);
  closeGateKeeperStore(&my_gkp_store);
  closeGateKeeperStore(my_gkp_store);

}

void print_olap(Long_Olap_Data_t olap,FILE*fs,char *space){
  fprintf (fs,"%s%8d %8d %c %5d %5d %4.1f %4.1f\n",
           space,
           olap . a_iid,
           olap . b_iid,
           olap . flipped ? 'I' : 'N',
           olap . a_hang, olap . b_hang,
           olap . orig_erate / 10.0, olap . corr_erate / 10.0);
}




uint64 iid2uid(uint32 iid){
  static uint64 *uid=NULL;
  if(uid==NULL){
    int i;
    int last = getLastElemFragStore (my_gkp_store);
    uid = (uint64 *)safe_malloc(sizeof(uint64)*(last+1));
    for(i=0;i<=last;i++){
      uid[i]=0;
    }
  }

  if(uid[iid]!=0){
    return (uid[iid]);
  } else {
    GateKeeperFragmentRecord gkpFrag;
    if(getGateKeeperFragmentStore(my_gkp_store.frg,iid,&gkpFrag)!=0)
      assert(0);
    uid[iid] = gkpFrag.UID;
  }
  return uid[iid];
}


int get_clr_len(uint32 iid){
  static int *len=NULL;
  if(len==NULL){
    int i;
    int last = getLastElemFragStore (my_gkp_store);
    len = (int *)safe_malloc(sizeof(int)*(last+1));
    for(i=0;i<=last;i++){
      len[i]=-1;
    }
  }

  if(len[iid]!=-1){
    return (len[iid]);
  } else {
    uint clr_bgn,clr_end;
    getFrag(my_gkp_store,iid,fsread,FRAG_S_INF);
    getClearRegion_ReadStruct(fsread, &clr_bgn,&clr_end, READSTRUCT_LATEST);
    len[iid]= clr_end-clr_bgn;
  }
  return len[iid];
}

				  
typedef struct dfs_node_tag {
  int orderReached;
  int finished;
  int parent;
  int height;
  int bestChild;
  int flipped;
  int numOvlsCompleted;
  int numOvls;
  Long_Olap_Data_t *ovls;
  double avgerr;
  double errvar;
} dfsGreedyFrg;


void init_dfsGreedyFrg(dfsGreedyFrg *f){
  f->orderReached=-1;
  f->finished=0;
  f->parent=0;
  f->height=0;
  f->bestChild=0;
  f->flipped=-1;
  f->numOvls=0;
  f->numOvlsCompleted=0;
  f->ovls=NULL;
}

typedef struct ovlfilt_tag {
  double erate;
  int useCorrected;
  int skipContaining;
  int minlen;
} overlapFilters;




int usefulOverlap(  Long_Olap_Data_t olap, int id, int offAEnd, overlapFilters filter){
  // exclude too-sloppy overlaps

  assert(olap.a_iid == id);

  if(filter.useCorrected&&olap.corr_erate>filter.erate*1000)return 0;

  if(!filter.useCorrected&&olap.orig_erate>filter.erate*1000)return 0;

  // exclude contained overlaps
  if( olap.a_hang > 0 && olap.b_hang < 0)return 0;

  // exclude containing overlaps
  if( filter.skipContaining && olap.a_hang < 0 && olap.b_hang > 0)return 0;

  // exclude too-short overlaps
  if(get_clr_len((uint)id)-MAX(0,olap.a_hang)+MIN(0,olap.b_hang) < filter.minlen)return 0;

  // if it's off the correct end ...
  if( (offAEnd ? (olap . a_hang < 0) : (olap.b_hang > 0))){
    return 1;
  }

  return 0;
}
       
int extra_height(int id, int offAEnd, Long_Olap_Data_t olap){
  if(offAEnd){
    // the extra is the part the B end of the top fragment sticks out -- which is
    //   > 0 only if the bhang is negative
    return MAX(0,-olap.b_hang);
  }

  // otherwise, the extra part is the A end stick out: the ahang
  return MAX(0,olap.a_hang);
}


int thickestSort(const void *A,const void *B){
  Long_Olap_Data_t *a = (  Long_Olap_Data_t *) A;
  Long_Olap_Data_t *b = (  Long_Olap_Data_t *) B;
  int thickA,thickB;
  
  assert(a->a_hang*a->b_hang >= 0);
  if(a->a_hang<0){
    thickA = get_clr_len(a->a_iid)+a->b_hang;
  } else {
    thickA = get_clr_len(a->a_iid)-a->a_hang;
  }

  assert(b->a_hang*b->b_hang >= 0);
  if(b->a_hang<0){
    thickB = get_clr_len(b->a_iid)+b->b_hang;
  } else {
    thickB = get_clr_len(b->a_iid)-b->a_hang;
  }
  return thickB-thickA;
} 

void DFS_longest_chain(dfsGreedyFrg *frgNodes,overlapFilters filter,int startingFrg){


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


  NOTES ON EFFICIENCY:

  We will either be 
  - inefficient, having to loop through the N overlaps of a node in O(N*N) steps ... plus a *lot* of time setting up iterators
  - memory bloated (if we add all overlaps to the stack)
  - tricky (if we let each node know how many of its overlaps we have already processed and we can jump to the i'th overlap)
  but still we have a lot of overlap iterators to set up

  */


  OVL_Stream_t  * my_stream = New_OVL_Stream ();
  Long_Olap_Data_t olap;

  //set current = first fragment
  int curr=startingFrg;
  int visited=0;
  int flipped=0;
  int last = Last_Frag_In_OVL_Store (my_ovl_store);
  int longest = 0;
  int maxlen = 0;

  int current_depth=1;


  // while (current <= last fragment )
  while(curr<=last){

    // set deadEnd = 1 -- signal that there were no (non-cycle) overlaps available
    int deadEnd = 1;

    // set finished = 1 -- signal that we did not hit an unvisited target and descend
    int finished = 1;

    assert(!frgNodes[curr].finished);

    //  mark current visit time
    if(    frgNodes[curr].orderReached == -1 ){
      frgNodes[curr].orderReached = ++visited;
      frgNodes[curr].flipped = flipped;

      { // set up overlaps
	int ovlsAlloced=0;
	Long_Olap_Data_t  *olaps=NULL;
	Long_Olap_Data_t  olap;
	int numOvls=0;
	double mean=0,var=0;

	// init overlap iterator
	Init_OVL_Stream_Intra_Frg (my_stream, curr, 0, my_ovl_store);

	// foreach overlap
	while  (Next_From_OVL_Stream (& olap, my_stream)){

	  if( !usefulOverlap(olap,curr,flipped,filter) ){
	    continue;
	  }
	  numOvls++;
	  if(numOvls>ovlsAlloced){
	    ovlsAlloced+=50;
	    if(frgNodes[curr].ovls==NULL){
	      frgNodes[curr].ovls=(Long_Olap_Data_t*)safe_malloc(sizeof(Long_Olap_Data_t)*ovlsAlloced);
	    } else {
	      frgNodes[curr].ovls=(Long_Olap_Data_t*)safe_realloc(frgNodes[curr].ovls,
							       sizeof(Long_Olap_Data_t)*ovlsAlloced);
	    }
	  }
	  frgNodes[curr].ovls[numOvls-1] = olap;
	  mean += olap.orig_erate;
	  var += olap.orig_erate * olap.orig_erate;
	}
	if(numOvls>0){
	  frgNodes[curr].ovls=(Long_Olap_Data_t*)safe_realloc(frgNodes[curr].ovls,
							   sizeof(Long_Olap_Data_t)*numOvls);
	  qsort(frgNodes[curr].ovls,numOvls,sizeof(Long_Olap_Data_t),thickestSort);
	  mean /= (double) numOvls;
	  frgNodes[curr].avgerr = mean;
	  frgNodes[curr].errvar = var/(double)numOvls - mean*mean;
	}else {
	  frgNodes[curr].avgerr = -1;
	  frgNodes[curr].errvar = -1;
	}

	frgNodes[curr].numOvls = numOvls;
      }


    } else {
      flipped = frgNodes[curr].flipped;
    }

    //#ifdef DEBUG_DFS
    fprintf(stderr,"dfs %d visit %d   %s   %f %f %d\n",curr,frgNodes[curr].orderReached,flipped ? "<--" : "-->",frgNodes[curr].avgerr,frgNodes[curr].errvar,frgNodes[curr].numOvls);
    //#endif



    while( frgNodes[curr].numOvls> frgNodes[curr].numOvlsCompleted){
      int target;
      int targetWouldBeFlippedGlobally;
      int wouldBeHeight=-1;
      
      Long_Olap_Data_t olap=frgNodes[curr].ovls[frgNodes[curr].numOvlsCompleted] ;

      //      memcpy((void*)&olap, (void*)(frgNodes[curr].ovls+frgNodes[curr].numOvlsCompleted),sizeof(Long_Olap_Data_t));

      assert(olap.a_iid==curr);

      target = olap.b_iid;
      targetWouldBeFlippedGlobally = (flipped != olap.flipped);

#ifdef DEBUG_DFS
      fprintf(stderr," overlaps %d -- visit %d finished %d height %d flipped %d\n",
              target,
              frgNodes[target].orderReached,
              frgNodes[target].finished,
              frgNodes[target].height,
              frgNodes[target].flipped ? "<--" : "-->");
      print_olap(olap,stderr,"   "); 
#endif
      
      // filter out overlaps we aren't interested in -- wrong end or not good enough
      if( (frgNodes[target].flipped !=-1 && frgNodes[target].flipped != targetWouldBeFlippedGlobally ) ){

	//#ifdef DEBUG_DFS
	fprintf(stderr," bad orientation! curr %d is %s olap.flipped=%d target %d is %s (visited %d)\n",
		curr,flipped?"<--":"-->",olap.flipped,target,frgNodes[target].flipped?"<--":"-->",frgNodes[target].orderReached); 
	//#endif
	frgNodes[curr].numOvlsCompleted++;
	continue;
      }

      // if target currently in stack (i.e. visit time valid but not marked finished), continue
      if (frgNodes[target].orderReached>0&&frgNodes[target].finished==0){
#ifdef DEBUG_DFS
	fprintf(stderr,"  cycle break!\n"); 
#endif
	frgNodes[curr].numOvlsCompleted++;
	continue;
      }

      // if we got here, there is at least one path to follow from current
      deadEnd=0;

      // else if target not yet visited, 
      if(frgNodes[target].orderReached==-1){

#ifdef DEBUG_DFS
	fprintf(stderr," drop down!\n"); 
#endif

	current_depth++;
#ifdef GRAPHIC_DEBUG_DFS
#if GRAPHIC_DEBUG_DFS > 1
	fprintf(stderr,".");
#else
	fprintf(stderr,"%d\n",current_depth);
#endif
#endif

	// set finished = 0 
	finished=0;

        // set target's parent as current
	frgNodes[target].parent=curr;

        // set current to target
	curr=target;

	// set the new orientation
	flipped = targetWouldBeFlippedGlobally;
 
        // break out of loop over overlaps -- we are going down to a lower node!
	break;

      }
    
      // else (i.e. if target is already marked finished)


#ifdef DEBUG_DFS
      fprintf(stderr," found a finished node -- evalute:\n"); 
#endif

      // if current's height through this one is an improvement
      wouldBeHeight = frgNodes[target].height + extra_height(curr,flipped,olap);
	

#ifdef DEBUG_DFS
      fprintf(stderr," height %d", wouldBeHeight);
#endif
	
#define MIN_ADVANTAGE 10
      if(wouldBeHeight > frgNodes[curr].height + MIN_ADVANTAGE){
	  
#ifdef DEBUG_DFS
	fprintf(stderr," new best!\n");
#endif
	// set height
	frgNodes[curr].height = wouldBeHeight;

	// set best child
	frgNodes[curr].bestChild = target;

      } else {
#ifdef DEBUG_DFS
        fprintf(stderr," no improvement\n");
#endif
      }
      frgNodes[curr].numOvlsCompleted++;
    }


    // if finished (i.e. we didn't go down farther)
    if(finished){

#ifdef DEBUG_DFS
      fprintf(stderr," finished working on %d\n",curr);
#endif
      if(frgNodes[curr].numOvls>0){
	free(frgNodes[curr].ovls);
	frgNodes[curr].ovls=NULL;
      }

      if(deadEnd){

	//#ifdef DEBUG_DFS
	fprintf(stderr," end of line!\n");
	//#endif

	frgNodes[curr].height = get_clr_len((uint)curr);
      }

      // mark current finished
      frgNodes[curr].finished=1;


      // if current has a parent
      if(frgNodes[curr].parent!=0){

	fprintf(stderr,"Node %d depth %d, height %d\n",curr, current_depth, frgNodes[curr].height);

        // set current = current's parent
	curr=frgNodes[curr].parent;

#ifdef DEBUG_DFS
	fprintf(stderr," pop up to parent %d!\n",curr);
#endif

	current_depth--;

#ifdef GRAPHIC_DEBUG_DFS
#if GRAPHIC_DEBUG_DFS > 1
	fprintf(stderr,"\n");
	{ 
	  int z;
	  for(z=1;z<current_depth;z++){
	    fprintf(stderr," ");
	  }
	}
#else
	fprintf(stderr,"%d\n",current_depth);
#endif
#endif

      } else {


	fprintf(stderr,"Node %d top level, height %d\n",curr, frgNodes[curr].height);
	if(frgNodes[curr].height>maxlen){
	  maxlen=frgNodes[curr].height;
	  longest=curr;
	}
	
        // advance current to next unfinished node
	while(curr<=last&&frgNodes[curr].finished){
	  curr++;
	}

	fprintf(stderr,"Advance to next unfinished node %d\n",curr);

      }
    }
  }


  curr=longest;
  assert(curr>0);
  while(curr!=0){
    printf("%d %d %d %s %f %f %d\n",
	   curr,
	   maxlen-frgNodes[curr].height,
	   maxlen-frgNodes[curr].height+get_clr_len(curr),
	   frgNodes[curr].flipped ? "<--" : "-->",
	   frgNodes[curr].avgerr,frgNodes[curr].errvar,
	   frgNodes[curr].numOvls);
    curr=frgNodes[curr].bestChild;
  }

}


void usage(char *pgm){
  fprintf (stderr, 
           "USAGE:  %s -g <GkpStoreName> -o <OvlStoreName> [-C] [-e <erate cutoff>] [-E] [-m <minlen>] [-N startingFrg]\n"
           "\t-C specifies to use containing fragments; default is only dovetails\n"
           "\t-e specifies the maximum mismatch rate (as a fraction, i.e. .01 means one percent)\n"
           "\t-E specifies that the corrected error rate rather than the original is to be used\n"
           "\t-N specifies the first fragment to examine; useful for partial ovlStore\n",
           pgm);
}


int main (int argc , char * argv[] ) {

  char full_gkpPath[1000];
  char full_ovlPath[1000];
  int setFullGkp=0;
  int setFullOvl=0;
  int i;
  int last;

  static dfsGreedyFrg *frgNodes=NULL;
  overlapFilters filter;
  int startingFrg=1; 

  setbuf(stdout,NULL);

  filter.erate = .3;
  filter.useCorrected = 0;
  filter.minlen = 40;
  filter.skipContaining = 1;


  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,"Ce:Ef:g:m:N:o:")) != EOF)){
      switch(ch) {
        case 'C':
          filter.skipContaining = 0;
          break;
        case 'e':
          filter.erate = atof(optarg);
          break;
        case 'E':
          filter.useCorrected=1;
          break;
        case 'g':
          strcpy(full_gkpPath, argv[optind - 1]);
          setFullGkp=1;
          break;
        case 'm':
          filter.minlen=atoi(optarg);
          break;
        case 'N':
          startingFrg = atoi(optarg);
          break;
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

    if(!setFullOvl || !setFullGkp )
      {
	usage(argv[0]);
	exit (-1);
      }
  }
  
  setup_stores(full_ovlPath,full_gkpPath);
  fsread = new_ReadStruct();

  // really ought to test store validity

  last =  Last_Frag_In_OVL_Store (my_ovl_store);
  frgNodes = (dfsGreedyFrg *) safe_malloc(sizeof(dfsGreedyFrg)*(last+1));
  for(i=1;i<=last;i++){
    init_dfsGreedyFrg(frgNodes+i);
  }

  DFS_longest_chain(frgNodes,filter,startingFrg);


  finished_with_stores();  
}

