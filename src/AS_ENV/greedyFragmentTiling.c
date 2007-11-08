
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
#include "AS_ALN_aligners.h"
#include "AS_ALN_forcns.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_Hash.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_gkpStore.h"
#include "UtilsREZ.h"
#include "AS_GKP_include.h"


// WAS: static OVL_Store_t  * my_store = NULL;
static OverlapStore  * my_store = NULL;
// WAS: static OVL_Store_t  * my_second_store = NULL;
static OverlapStore  * my_second_store = NULL;
static  GateKeeperStore *my_gkp_store;
static fragRecord fsread;
static int useCorrectedErate=0;
static int maxOvlsToConsider=-1;
static int *IIDusability;
static int *SeedUsability;
static int restrictIDs=0;
static int restrictSeeds=0;
static char iidListFile[250];
static char seedListFile[250];
static int *iid2sample=NULL;
static char *seen;
static int thickestOvlsCountAsSeen=0;
static int seedSample=-1;

#define DEFAULT_SAMPLE_ADVANTAGE .05

void setup_stores(char *OVL_Store_Path, char *Gkp_Store_Path){

  assert(OVL_Store_Path!=NULL);
  assert(Gkp_Store_Path!=NULL);

  assert(my_store==NULL);
  //WAS:   my_store = New_OVL_Store ();
  //WAS:  Open_OVL_Store (my_store, OVL_Store_Path);
  my_store = AS_OVS_openOverlapStore(OVL_Store_Path);

  assert(my_second_store==NULL);
  //WAS:  my_second_store = New_OVL_Store ();
  //WAS:  Open_OVL_Store (my_second_store, OVL_Store_Path);
  my_second_store = AS_OVS_openOverlapStore(OVL_Store_Path);

  my_gkp_store = openGateKeeperStore( Gkp_Store_Path, FALSE);
}

//WAS: void print_olap(Long_Olap_Data_t olap){
void print_olap(OVSoverlap olap){
  printf ("    %8d %8d %c %5d %5d %4.1f %4.1f\n",
	  olap . a_iid,
	  olap . b_iid,
	  olap . dat . ovl . flipped ? 'I' : 'N',
	  olap . dat . ovl . a_hang, olap . dat . ovl . b_hang,
	  AS_OVS_decodeQuality(olap . dat . ovl . orig_erate)*100.0,
	  AS_OVS_decodeQuality(olap . dat . ovl . corr_erate)*100.0);
}


typedef enum {
  BEST_MEANS_THICKEST=1,
  BEST_MEANS_LOWEST_ERROR=2,
  BEST_MEANS_BEST_UPPER_BOUND=3
}BestMeasure;



AS_UID iid2uid(uint32 iid){
  static AS_UID *uid=NULL;

  if(uid==NULL){
    int i;
    int last = getLastElemFragStore (my_gkp_store);
    uid = (AS_UID *)safe_malloc(sizeof(AS_UID)*(last+1));
    for(i=0;i<=last;i++){
      uid[i] = AS_UID_undefined();
    }
  }

  if(AS_UID_isDefined(uid[iid]) == FALSE){
    GateKeeperFragmentRecord gkpFrag;
    getGateKeeperFragment(my_gkp_store,iid,&gkpFrag);
    assert(gkpFrag.readIID == iid);
    uid[iid] = gkpFrag.readUID;
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
    getFrag(my_gkp_store,iid,&fsread,FRAG_S_INF);
    clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_LATEST);
    clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_LATEST);

    len[iid]= clr_end-clr_bgn;
  }
  return len[iid];
}



void olap_ref_swap(OVSoverlap *o){

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
  if(o-> dat . ovl . flipped){
    t=o-> dat . ovl .a_hang;
    o-> dat . ovl .a_hang=o-> dat . ovl .b_hang;
    o-> dat . ovl .b_hang=t;
  } else {
    o-> dat . ovl .a_hang*=-1;
    o-> dat . ovl .b_hang*=-1;
  }
}

void olap_reverse(OVSoverlap *o){

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
  int t=-o-> dat . ovl . a_hang;
  o-> dat . ovl .a_hang=-o-> dat . ovl .b_hang;
  o-> dat . ovl .b_hang=t;
  o-> dat . ovl .flipped = 1 - (o-> dat . ovl .flipped);
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

  //WAS:  OVL_Stream_t  * my_stream = New_OVL_Stream ();
  int retval=0;
  //WAS:  Long_Olap_Data_t  olap;
  OVSoverlap  olap;

  assert(my_second_store!=NULL);

  //WAS:  Init_OVL_Stream (my_stream, frg, frg, my_second_store);
  AS_OVS_setRangeOverlapStore(my_second_store,frg,frg);


  while  (AS_OVS_readOverlapFromStore (my_second_store, &olap, AS_OVS_TYPE_OVL )){

    // exclude too-sloppy overlaps

    if(useCorrected && olap . dat . ovl . corr_erate >erate)continue;

    if(!useCorrected&&olap . dat . ovl.orig_erate>erate)continue;

    assert(olap.a_iid == frg);

    // exclude contained overlaps
    if( olap . dat . ovl .a_hang > 0 && olap . dat . ovl .b_hang < 0)continue;

    // exclude containing overlaps
    if( skipContaining && olap. dat . ovl.a_hang < 0 && olap. dat . ovl.b_hang > 0)continue;

    if(frglen-MAX(0,olap. dat . ovl.a_hang)+MIN(0,olap. dat . ovl.b_hang) < minlen)continue;

    // if it's off the correct end ...
    if( (offAEnd ? (olap. dat . ovl . a_hang < 0) : (olap. dat . ovl.b_hang > 0))){
      return 0;
    }
  } 
  return 1;
}



OVSoverlap better_olap(OVSoverlap a, OVSoverlap b, int startingFrg, int offAEnd, BestMeasure bestType, int frglen, 
			     int useCorrected, int skipContaining, int minlen, int avoidDeadEnds, double erate ){

  int aThick, bThick;
  int aMis, bMis;
  int aLen, bLen;

  OVSoverlap copy_a =a;
  OVSoverlap copy_b =b;

  if(avoidDeadEnds){ // check for extension off far end of other fragment

    int aIsDead=0;
    int bIsDead=0;

    int aOtherAEnd;  // whether the far end of the other fragment of overlap "a" is its A-end or not
    int bOtherAEnd;  // whether the far end of the other fragment of overlap "b" is its A-end or not

    assert(a.a_iid==startingFrg);
    assert(b.a_iid==startingFrg);

    if(offAEnd){
      aOtherAEnd = !a . dat . ovl .flipped;
      bOtherAEnd = !b . dat . ovl.flipped;
    } else {
      aOtherAEnd = a . dat . ovl.flipped;
      bOtherAEnd = b . dat . ovl.flipped;
    }

    aIsDead = isDeadEnd(a.b_iid,aOtherAEnd,erate,useCorrected,skipContaining,minlen,get_clr_len(a.b_iid));
    bIsDead = isDeadEnd(b.b_iid,bOtherAEnd,erate,useCorrected,skipContaining,minlen,get_clr_len(b.b_iid));
    
    // if status is different, then prefer the live one ...
    if(aIsDead<bIsDead)return a;
    if(aIsDead>bIsDead)return b;
  }

  switch (bestType){

    case BEST_MEANS_LOWEST_ERROR:
      if(useCorrected){ 
        if(a.dat . ovl . corr_erate>=b.dat . ovl . corr_erate){
          return b;
        } else {
          return a;
        }
      } else {
        if(a.dat . ovl . orig_erate>=b.dat . ovl . orig_erate){
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
      assert(copy_a.dat . ovl . b_hang >0 );
      assert(copy_b.dat . ovl . b_hang >0 );
      if(copy_a. dat . ovl . a_hang>=copy_b. dat . ovl . a_hang){
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
      // overlap = frglen - MAX(0,ahang) + MIN(0,bhang)
      int  olenA = frglen - MAX(0,a. dat . ovl . a_hang) + MIN(0,a.dat . ovl . b_hang);
      int  olenB = frglen - MAX(0,b. dat . ovl . a_hang) + MIN(0,b.dat . ovl . b_hang);

      if(useCorrected){
        upperA = get_95pct_upper(olenA, AS_OVS_decodeQuality(a.dat . ovl . corr_erate)+.0009);
        upperB = get_95pct_upper(olenB, AS_OVS_decodeQuality(b.dat . ovl . corr_erate)+.0009);
      } else {
        upperA = get_95pct_upper(olenA, AS_OVS_decodeQuality(a.dat . ovl . orig_erate)+.0009);
        upperB = get_95pct_upper(olenB, AS_OVS_decodeQuality(b.dat . ovl . orig_erate)+.0009);
      }
    

      }
      break;

    default:
      fprintf(stderr,"Unknown criterion for best overlap! exit \n");
      exit(-1);
      break;
  }
}

int lowestErrorSort(const void *A,const void *B){
  OVSoverlap *a = (  OVSoverlap *) A;
  OVSoverlap *b = (  OVSoverlap *) B;
  int thickA,thickB;

  // sort by error rate
  if(useCorrectedErate){
    if(a->dat . ovl . corr_erate!=b->dat . ovl . corr_erate){
      return a->dat . ovl . corr_erate-b->dat . ovl . corr_erate;
    }
  } else {
    if(a->dat . ovl . orig_erate!=b->dat . ovl . orig_erate){
      return a->dat . ovl . orig_erate-b->dat . ovl . orig_erate;
    }
    /* in case of a tie in uncorrected ... */
    if(a->dat . ovl . corr_erate!=b->dat . ovl . corr_erate){
      return a->dat . ovl . corr_erate-b->dat . ovl . corr_erate;
    }
  }

  { // sort by thickness
    int laa = get_clr_len(a->a_iid);
    int  olenA = laa - MAX(0,a-> dat . ovl . a_hang) + MIN(0,a->dat . ovl . b_hang);
    int lba = get_clr_len(b->a_iid);
    int  olenB = lba - MAX(0,b-> dat . ovl . a_hang) + MIN(0,b->dat . ovl . b_hang);
    if(olenA!=olenB)return olenA-olenB;
  }

  return 0;
    
} 

// setupolaps() returns a count of usable olaps and sets *retolaps to point to
// a block of memory containing them ... but the memory is static (i.e. the 
// calling routine does not own the memory) and should not be freed!
int setupolaps(int id, int offAEnd,BestMeasure bestType,double erate,int useCorrected,int skipContaining,int minlen, int frglen,int avoidDeadEnds,OVSoverlap **retolaps, double favorSameSample,double favorSameSampleAsSeed){

  static OVSoverlap  *olaps=NULL;
  OVSoverlap  olap;
  static int ovlsAlloced=0;
  int numOvls=0;

  assert(my_store!=NULL);

  AS_OVS_setRangeOverlapStore(my_store,id,id);

  while  (AS_OVS_readOverlapFromStore(my_store, & olap, AS_OVS_TYPE_OVL)){

    assert(olap.a_iid == id);

    // exclude overlaps to fragments not in the use list
    if( restrictIDs && IIDusability[olap.b_iid] != 1)continue;

    // exclude too-sloppy overlaps
    if(useCorrected&&olap.dat . ovl . corr_erate>erate)continue;
    if(!useCorrected&&olap.dat . ovl . orig_erate>erate)continue;

    // exclude contained overlaps
    if( olap. dat . ovl . a_hang > 0 && olap.dat . ovl . b_hang < 0)continue;

    // exclude containing overlaps
    if( skipContaining && olap. dat . ovl . a_hang < 0 && olap.dat . ovl . b_hang > 0)continue;

    // exclude too-short overlaps
    if(frglen-MAX(0,olap. dat . ovl . a_hang)+MIN(0,olap.dat . ovl . b_hang) < minlen)continue;

    // if it's off the correct end ...
    if( (offAEnd ? (olap .  dat . ovl . a_hang < 0) : (olap.dat . ovl . b_hang > 0))){

      //    print_olap(olap);

      numOvls++;
      if(numOvls>ovlsAlloced){
	ovlsAlloced+=50;
	if(olaps==NULL){
	  olaps=(OVSoverlap*)safe_malloc(sizeof(OVSoverlap)*ovlsAlloced);
	} else {
	  olaps=(OVSoverlap*)safe_realloc(olaps,
						sizeof(OVSoverlap)*ovlsAlloced);
	}
      }
      if(bestType==BEST_MEANS_LOWEST_ERROR&&favorSameSample>0){
	if(iid2sample[olap.a_iid]==iid2sample[olap.b_iid]){
	  if(useCorrected){
	    if( ((int)olap.dat . ovl . corr_erate)-favorSameSample >= 0 ){
	      olap.dat . ovl . corr_erate-=favorSameSample;
	    } else {
	      olap.dat . ovl . corr_erate=0;
	    }
	  } else {
	    if( ((int)olap.dat . ovl . orig_erate)-favorSameSample >= 0 ){
	      olap.dat . ovl . orig_erate-=favorSameSample;
	    } else {
	      olap.dat . ovl . orig_erate=0;
	    }
	  }
	}
      }
      if(bestType==BEST_MEANS_LOWEST_ERROR&&favorSameSampleAsSeed>0){
	if(seedSample==iid2sample[olap.b_iid]){
	  if(useCorrected){
	    if( ((int)olap.dat . ovl . corr_erate)-favorSameSampleAsSeed >= 0 ){
	      olap.dat . ovl . corr_erate-=favorSameSampleAsSeed;
	    } else {
	      olap.dat . ovl . corr_erate=0;
	    }
	    if(seedSample!=iid2sample[olap.a_iid]){
	      if( ((int)olap.dat . ovl . corr_erate)-favorSameSampleAsSeed >= 0 ){
		olap.dat . ovl . corr_erate-=favorSameSampleAsSeed;
	      } else {
		olap.dat . ovl . corr_erate=0;
	      }
	    }
	  } else {
	    if( ((int)olap.dat . ovl . orig_erate)-favorSameSampleAsSeed >= 0 ){
	      olap.dat . ovl . orig_erate-=favorSameSampleAsSeed;
	    } else {
	      olap.dat . ovl . orig_erate=0;
	    }
	    if(seedSample!=iid2sample[olap.a_iid]){
	      if( ((int)olap.dat . ovl . orig_erate)-favorSameSampleAsSeed >= 0 ){
		olap.dat . ovl . orig_erate-=favorSameSampleAsSeed;
	      } else {
		olap.dat . ovl . orig_erate=0;
	      }
	    }
	  }
	}
      }
      olaps[numOvls-1] = olap;

    } else {
      // does not stick off the correct end
    }

  }

  qsort(olaps,numOvls,sizeof(OVSoverlap),lowestErrorSort);
  *retolaps=olaps;
  return numOvls;
}

int best_overlap_off_end(int id, int offAEnd,BestMeasure bestType,OVSoverlap *bestovl,double erate,int useCorrected,int skipContaining,int minlen, int frglen,int avoidDeadEnds,OVSoverlap *olaps, int numOvls){
  OVSoverlap  olap, bestolap;
  int i;
  int goodOlap=0;
  int retval=0;

  i=0;
  while(i<numOvls&&(maxOvlsToConsider == -1 
		    || (i < maxOvlsToConsider)
		    || (useCorrectedErate ? (olaps[i].dat . ovl . corr_erate == olaps[maxOvlsToConsider].dat . ovl . corr_erate) :  (olaps[i].dat . ovl . orig_erate == olaps[maxOvlsToConsider].dat . ovl . orig_erate)))){
		    
    /*
      while(i<numOvls&&(maxOvlsToConsider == -1 || i<maxOvlsToConsider)){
    */

    OVSoverlap olap = olaps[i];
    
    if(goodOlap==0){
      goodOlap=1;
      bestolap=olap; 
    } else {
      bestolap=better_olap(bestolap,olap, id, offAEnd,bestType,frglen,useCorrected,skipContaining,minlen,avoidDeadEnds,erate);
    }
    i++;
  }
  *bestovl=bestolap;

  if(goodOlap==0)return 0;
  return 1;
}



void finished_with_ovlStore(void){

  AS_OVS_closeOverlapStore (my_store);
  AS_OVS_closeOverlapStore (my_second_store);

}


void setUpRestrictions(int last_stored_frag,char *listfile,int *usabilityArray){
  FILE *frglist = fopen(listfile,"r");
  char range[256];
  int i;

  assert(frglist!=NULL);

  for(i=0;i<=last_stored_frag;i++){
    usabilityArray[i]=0;
  }

  // this should be made much more sophisticated, with token parsing etc
  while(fscanf(frglist,"%s",range)==1){
    if(strlen(range)>255){
      fprintf(stderr,"Sloppy coding: range string %s ran out of memory; exit!\n",range);
      exit(-4);
    }
    if(strstr(range,",")!=NULL){
      fprintf(stderr,"the IID restriction code doesn't really support comma-separated lists yet;exit!\n");
      exit(-4);
    }
    if(strstr(range,"-")!=NULL){
      int b,e,i;
      assert(strstr(range,"-")!=range);
      b = atoi(range);
      e = atoi(strstr(range,"-")+1);
      assert(b>0&&b<=e&&e<=last_stored_frag);
      for(i=b;i<=e;i++){
	usabilityArray[i]=1;
      }
    } else {
      int i = atoi(range); 
      assert(i>0&&i<=last_stored_frag);
      usabilityArray[i]=1;
    }
  }
  fclose(frglist);
}


void usage(char *pgm){
  fprintf (stderr, 
           "USAGE:  %s -g <GkpStoreName> -o <full_ovlStore> -n <startingFrgNum> [-C] [-e <erate cutoff>] [-E] [-m <minlen>] [-Q] [-D] [-N <maxovls>] [-i <file specifying IIDs to use> | -I <file specifying IIDs to use>] [-R] [-s <uid2sample file> [-S <same-sample bonus> | -T <same-sample-as-seed bonus>]] [-P] [-5 | -3]\n"
           "\t-n startingFrgNum = fragment to walk out from\n"
           "\t-e specifies the maximum mismatch rate (as a fraction, i.e. .01 means one percent)\n"
           "\t-B specifies that the overlap with the lowest upper bound on the mismatch rate will be used\n"
           "\t-C specifies to use containing fragments (default is only dovetails\n"
           "\t-D specifies that best quality should win over avoiding dead ends; by default, a fragment that is not a dead end wins over one that is\n"
           "\t-E specifies that the corrected error rate rather than the original is to be used\n"
           "\t-N specifies that only the maxovls lowest-error overlaps will be evaluated\n"
           "\t-Q specifies that the lowest error rate overlap will be used (in place of thickest)\n"
           "\t-i specifies a file specifying a comma or \\n separated list\n"
           "\t\tof either single IIDs or <start>-<end> ranges; IIDs thus\n"
           "\t\tspecified can be used in an assembly--others not!\n"
           "\t-I is like -i except that the restriction only applies to seeds,\n"
           "\t\tnot extensions\n"
           "\t-P causes error rate to be printed\n"
           "\t-R specifies that fragments involved in thicker overlaps than the one chosen at a given extension will count as used\n"
           "\t-s species the name of a file containing a uid to sample (integer) mapping; this implicitly has the effect\n"
           "\t-S specifies the amount (in percent error) that a same-sample overlap is preferred to a different-sample overlap\n"
           "\t-T specifies the amount (in percent error) that an overlap is preferred if it is to a fragment whose sample matches that of the seed\n"
           "\t-5 specifies that only the 5' end of the seed be extended\n"
           "\t-3 specifies that only the 3' end of the seed be extended\n"
           ,pgm);
}

AS_IID uid2iid(AS_UID uid){
  static int firstFailure=1;
  AS_IID     iid;
  iid = getGatekeeperUIDtoIID(my_gkp_store, uid, 0);
  if (iid == 0) {
    if (firstFailure) {
      fprintf(stderr,"Tried to look up iid of unknown uid: %s; this may reflect trying to use a deleted fragment; further instances will not be reported.\n", AS_UID_toString(uid));
      firstFailure=0;
    }
    return (-1);
  }
  return (iid);
}

int ovlThickness(OVSoverlap o, int frglen, int afrg, int Aend){

  /* thickness is afrg's length minus
     bhang (from afrg's perspective) if Aend
     ahang (from afrg's perspective) if Bend

     N.B. this makes thickness > frglen if relevant hang is negative
  */

  int change;

  if(o.a_iid == afrg){
    return frglen + (Aend ? o.dat . ovl . b_hang : -o. dat . ovl . a_hang);
  } 

  if(o.dat . ovl . flipped){
    if(Aend){
      change = o. dat . ovl . a_hang;
    } else {
      change = -o. dat . ovl . a_hang;
    }
  } else {
    if(Aend){
      change = -o.dat . ovl . b_hang;
    } else {
      change = o. dat . ovl . a_hang;
    }
  }
  return frglen + change;
}



int main (int argc , char * argv[] ) {

  Global_CGW *data;
  char *prefix;
  int ckptNum = NULLINDEX;
  int i, index;
  char subset_map[1000];
  char full_map[1000];
  char setSubsetMap=0;
  char setFullMap=0;
  char full_ovlPath[1000];  
  char full_gkpPath[1000];
  char sampleFileName[1000];
  double favorSameSample=0;
  double favorSameSampleAsSeed=0;
  FILE *sampleFile;  
  int setFullGkp=0;
  int setFullOvl=0;
  int bestType = BEST_MEANS_THICKEST;
  int startingFrgNum=0;
  int seediid;
  int stillGoing=0;
  int Aend=-1; 
  int currFrg=-1;
  int firstExtend=-1;

  int fivePonly=0;
  int threePonly=0;

  int printpid=0;
  double pid=-1.;

  int last_stored_frag;
  int ahang,bhang;
  int rightEnd,leftEnd;
  uint frglen,clr_bgn,clr_end;
  double maxError=.3;
  int skipContaining=1;
  int minlen=40;
  int avoidDeadEnds=1;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  sampleFileName[0]='\0';

  argc = AS_configure(argc, argv);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,"BCDe:Ef:g:m:n:N:o:PQs:S:T:i:I:R53")) != EOF)){
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
          maxError = AS_OVS_encodeQuality(atof(optarg));
          break;
        case 'E':
          useCorrectedErate=1;
          break;
        case 'g':
          strcpy( full_gkpPath, argv[optind - 1]);
          setFullGkp = TRUE;
          break;
        case 'i':
          restrictIDs=1;
          {
            int n ;
            char *c= strncpy(iidListFile,optarg,249);
            n = strlen(c);
            if(n==249&&optarg[249]!='\0'){
              fprintf(stderr,"Not enough memory for iid list path!\n");
              exit(-4);
            }
          }
          break;
        case 'I':
          restrictSeeds=1;
          {
            char *c = strncpy(seedListFile,optarg,249);
            int n;
            n = strlen(c);
            if(n==249&&optarg[249]!='\0'){
              fprintf(stderr,"Not enough memory for iid list path!\n");
              exit(-4);
            }
          }
          break;
        case 'm':
          minlen=atoi(optarg);
          break;
        case 'n':
          startingFrgNum = atoi(argv[optind - 1]);
          break;
        case 'N':
          maxOvlsToConsider = atoi(argv[optind-1]);
          break;
        case 'o':
          strcpy(full_ovlPath,argv[optind-1]);
          setFullOvl=1;
          break;
        case 'P':
          printpid=1;
          break;
        case 'Q':
          bestType=BEST_MEANS_LOWEST_ERROR;
          break;
        case 'R':
          thickestOvlsCountAsSeen=1;
          break;
        case 's':
          strcpy(sampleFileName,argv[optind-1]);
          assert(strlen(sampleFileName)<=999);
          break;
        case 'S':
          favorSameSample=AS_OVS_encodeQuality(atof(argv[optind-1]));
          assert(favorSameSample>0);
          break;
        case 'T':
          favorSameSampleAsSeed=AS_OVS_encodeQuality(atof(argv[optind-1]));
          assert(favorSameSampleAsSeed>0);
          break;
        case '3':
          threePonly=1;
          break;
        case '5':
          fivePonly=1;
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if( !setFullGkp || !setFullOvl )
      {
	usage(argv[0]);
	exit (-1);
      }
  }
  
  assert(! (fivePonly&&threePonly) );

  assert(AS_CNS_ERROR_RATE >= AS_OVS_decodeQuality(maxError));

  setup_stores(full_ovlPath, full_gkpPath);
  last_stored_frag = getLastElemFragStore (my_gkp_store);

  if(restrictIDs){
    IIDusability = (int *) safe_malloc(sizeof(int)*(last_stored_frag+1));
    setUpRestrictions(last_stored_frag,iidListFile,IIDusability);
  }

  if(restrictSeeds){
    SeedUsability = (int *) safe_malloc(sizeof(int)*(last_stored_frag+1));
    setUpRestrictions(last_stored_frag,seedListFile,SeedUsability);
  }

  if(favorSameSample>0||favorSameSampleAsSeed>0){
    if(favorSameSampleAsSeed>0){
      favorSameSample=0;
    }
    assert(sampleFileName[0]!='\0');
  }
#if 1
  fprintf(stderr, "sampleFileName is broken.  Blame BPW.\n");
  exit(1);
#else
  if(sampleFileName[0]!='\0'){
    AS_UID uid;
    uint32 smp;
    AS_IID iid;
    int i;
    iid2sample = (int *) safe_malloc(sizeof(int)*(last_stored_frag+1));
    for(i=0;i<=last_stored_frag;i++){
      iid2sample[i]=-1;
    }
    sampleFile = fopen(sampleFileName,"r");
    assert(sampleFile!=NULL);
#error broken
    while(fscanf(sampleFile, "%s "F_IID,&uid,&smp)==2){
      AS_IID iid=uid2iid(uid);
      if(iid>0) iid2sample[iid]=smp;
    }
    if(favorSameSample==0&&favorSameSampleAsSeed==0){
      //      favorSameSample=AS_OVS_encodeQuality(DEFAULT_SAMPLE_ADVANTAGE); 
      fprintf(stderr,"Unexpected settings in %s :  -s <read UID to sample file> specified without positive value for -S or -T -- abort!\n",argv[0]);
      exit(-1);
    }
  }
#endif

  seen = (char*) safe_malloc((last_stored_frag+1)*sizeof(char));
  {
    for(i=0;i<=last_stored_frag;i++){
      seen[i]='\0';
    }
  }



  if(startingFrgNum==0) {
    startingFrgNum=1;
  } else {
    assert(startingFrgNum>0&&startingFrgNum<last_stored_frag);
    last_stored_frag=startingFrgNum;
  }

  if(last_stored_frag>AS_OVS_lastFragInStore(my_store)){
    last_stored_frag=AS_OVS_lastFragInStore(my_store);
  }

  for(seediid=startingFrgNum;seediid<=last_stored_frag;seediid++){
    int numOvls=0,numFwdOvls=0;
    OVSoverlap *olaps;
    int built5End = 0;

    // skip deleted fragments
    GateKeeperFragmentRecord gkpFrag;
    getGateKeeperFragment(my_gkp_store,seediid,&gkpFrag);
    if(gkpFrag.deleted)continue;
    assert(gkpFrag.readIID);

    // if this fragment is not in the select list(s) ... skip it
    if(restrictSeeds&& SeedUsability[seediid]==0){continue;}
    if(restrictIDs&&IIDusability[seediid]==0){continue;}

    // if already seen, skip it
    if(seen[seediid]!='\0')continue;

    // otherwise, use it as a seed ...

    printf("Seed %d\n",seediid);

    // keep track of the sample the seed comes from
    if(favorSameSampleAsSeed){
      seedSample=iid2sample[seediid];
    }


    // do forward extensions
    stillGoing=1;
    Aend=1;
    currFrg=seediid;
    seen[currFrg]='\1';
    ahang=0;
    getFrag(my_gkp_store,currFrg,&fsread,FRAG_S_INF);

    clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_LATEST);
    clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_LATEST);
    frglen=clr_end-clr_bgn;
    numFwdOvls = numOvls = setupolaps(currFrg, Aend, bestType,maxError,useCorrectedErate,skipContaining,minlen,frglen,avoidDeadEnds,&olaps,favorSameSample,favorSameSampleAsSeed);
    if(numOvls==0||threePonly){
      stillGoing=0;
    }
    while(stillGoing){
      OVSoverlap o;
      if(printpid){
	if(iid2sample!=NULL){
	  printf("%d  %d %d %s %f %d\n",currFrg,ahang,ahang+frglen, Aend?"<--":"-->",AS_OVS_decodeQuality(pid),iid2sample[currFrg]);
	} else {
	  printf("%d  %d %d %s %f\n",currFrg,ahang,ahang+frglen, Aend?"<--":"-->",AS_OVS_decodeQuality(pid));
	}
      } else {
	printf("%d  %d %d %s\n",currFrg,ahang,ahang+frglen, Aend?"<--":"-->");
      }
      built5End = 1;
      
      if(numOvls>0){
	stillGoing = best_overlap_off_end(currFrg, Aend, bestType,&o,maxError,useCorrectedErate,skipContaining,minlen,frglen,avoidDeadEnds,olaps,numOvls);
	if(stillGoing){

	  pid=o.dat . ovl . orig_erate;

	  if(thickestOvlsCountAsSeen){
	    int idx;
	    int bestotheriid=( currFrg == o.a_iid ? o.b_iid : o.a_iid );
	    int bestthickness = ovlThickness(o,frglen,currFrg,Aend);
	    for(idx=0;idx<numOvls;idx++){
	      int otheriid = ( currFrg == olaps[idx].a_iid ? olaps[idx].b_iid : olaps[idx].a_iid );
	      if(otheriid==seediid||otheriid==bestotheriid)continue;
	      if(bestthickness<ovlThickness(olaps[idx],frglen,currFrg,Aend)){
		if(!seen[otheriid])seen[otheriid]++;
	      }
	    }
	  }

	  if(currFrg==seediid){
	    firstExtend= ( o.a_iid==currFrg) ? o.b_iid : o.a_iid;
	  }
	  if(o.a_iid==currFrg){
	    currFrg=o.b_iid;

	    if(Aend){
	      ahang+=-o.dat . ovl . b_hang;
	    } else {
	      ahang+=o. dat . ovl . a_hang;
	    }

	  } else {
	    currFrg=o.a_iid;

	    if(Aend){
	      if(o.dat . ovl . flipped){
		ahang+=-o. dat . ovl . a_hang;
	      } else {
		ahang+=o.dat . ovl . b_hang;
	      }
	    } else {
	      if(o.dat . ovl . flipped){
		ahang+=o.dat . ovl . b_hang;
	      } else {
		ahang+=-o.dat . ovl . b_hang;
	      }
	    }
	  }

	  if(Aend){
	    Aend=1-o.dat . ovl . flipped;
	  } else {
	    Aend=o.dat . ovl . flipped;
	  }

	}
	if(!seen[currFrg]){
	  getFrag(my_gkp_store,currFrg,&fsread,FRAG_S_INF);
	  clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_LATEST);
	  clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_LATEST);
	  frglen=clr_end-clr_bgn;
	  numOvls = setupolaps(currFrg, Aend, bestType,maxError,useCorrectedErate,skipContaining,minlen,frglen,avoidDeadEnds,&olaps,favorSameSample,favorSameSampleAsSeed);
	} else {
	  printf("SEEN: %d  %d %d %s\n",currFrg,ahang,ahang+frglen, Aend?"<--":"-->");
	  stillGoing=0;
	}   
	seen[currFrg]++;
   assert(seen[currFrg]<(char)127);
      } else {
	stillGoing=0;
      }
    }

    seen[seediid]--; // decrement the seen counter for the seed frg

    // if we hit the seed again, don't go back off the other end
    if(seen[seediid]!=0){
      continue;
    }

    //    fprintf(stderr,"now for the other end...\n");
    
    stillGoing=1;
    Aend=0;
    currFrg=seediid;
    getFrag(my_gkp_store,currFrg,&fsread,FRAG_S_INF);
    clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_LATEST);
    clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_LATEST);
    frglen=clr_end-clr_bgn;

    numOvls = setupolaps(currFrg, Aend, bestType,maxError,useCorrectedErate,skipContaining,minlen,frglen,avoidDeadEnds,&olaps,favorSameSample,favorSameSampleAsSeed);
    if(numOvls==0||fivePonly){stillGoing=0;}

    rightEnd=frglen; /* should be frglen of seed fragment */
    leftEnd=0;

    while(stillGoing){
      OVSoverlap o;
      // when we had no overlaps on the 5' end, we need to also output the seed info or else it is never mentioned
      if(currFrg!=seediid||firstExtend==-1||built5End==0){
	if(printpid){
	  if(iid2sample!=NULL){
	    printf("%d %d %d %s %f %d\n",currFrg,leftEnd,rightEnd,!Aend?"<--":"-->",AS_OVS_decodeQuality(pid),iid2sample[currFrg]);
	  } else {
	    printf("%d %d %d %s %f\n",currFrg,leftEnd,rightEnd,!Aend?"<--":"-->",AS_OVS_decodeQuality(pid));
	  }
	} else {
	  printf("%d %d %d %s\n",currFrg,leftEnd,rightEnd,!Aend?"<--":"-->");
	}
      }      

      if(numOvls>0){
	stillGoing = best_overlap_off_end(currFrg, Aend, bestType,&o,maxError,useCorrectedErate,skipContaining,minlen,frglen,avoidDeadEnds,olaps,numOvls);
	if(stillGoing){
	  pid=o.dat . ovl . orig_erate;

	  if(thickestOvlsCountAsSeen){
	    int idx;
	    int bestotheriid=( currFrg == o.a_iid ? o.b_iid : o.a_iid );
	    int bestthickness = ovlThickness(o,frglen,currFrg,Aend);
	    for(idx=0;idx<numOvls;idx++){
	      int otheriid = ( currFrg == olaps[idx].a_iid ? olaps[idx].b_iid : olaps[idx].a_iid );
	      if(otheriid==seediid||otheriid==bestotheriid)continue;
	      if(bestthickness<ovlThickness(olaps[idx],frglen,currFrg,Aend)){
		if(!seen[otheriid])seen[otheriid]++;
	      }
	    }
	  }
 

	  if(o.a_iid==currFrg){
	    currFrg=o.b_iid;

	    if(Aend){
	      rightEnd=leftEnd + frglen + o.dat . ovl . b_hang;
	      leftEnd+=o. dat . ovl . a_hang;
	    } else {
	      rightEnd=leftEnd + frglen - o.dat . ovl . a_hang;
	      leftEnd-=o.dat . ovl . b_hang;
	    }

	  } else {

	    fprintf(stderr,"Diagnostic note: I was not sure this would ever happen -- ALH\n");

	    currFrg=o.a_iid;

	    if(Aend==o.dat . ovl . flipped){
	      rightEnd+=o. dat . ovl . a_hang;
	      leftEnd+=o.dat . ovl . b_hang;
	    } else {
	      rightEnd-=o.dat . ovl . b_hang;
	      leftEnd-=o. dat . ovl . a_hang;
	    }
	  }
	  if(Aend){
	    Aend=1-o.dat . ovl . flipped;
	  } else {
	    Aend=o.dat . ovl . flipped;
	  }
	}
	if(!seen[currFrg]){
	  seen[currFrg]='\1';
	  getFrag(my_gkp_store,currFrg,&fsread,FRAG_S_INF);
	  clr_bgn = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_LATEST);
	  clr_end = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_LATEST);
	  frglen=clr_end-clr_bgn;
	  numOvls = setupolaps(currFrg, Aend, bestType,maxError,useCorrectedErate,skipContaining,minlen,frglen,avoidDeadEnds,&olaps,favorSameSample,favorSameSampleAsSeed);
	} else {
	  stillGoing=0;
	  printf("SEEN: %d  %d %d %s\n",currFrg,leftEnd,rightEnd, Aend?"<--":"-->");
	}
	seen[currFrg]++;
	assert(seen[currFrg]<(char)127);
      } else {
	stillGoing=0;
      }
    }
     
     // mark seed as seen so it cannot be used again
     seen[seediid]++;
  }

  finished_with_ovlStore();  
  exit(0);
}

