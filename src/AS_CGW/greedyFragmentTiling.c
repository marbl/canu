
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

static HashTable_AS *subset_name2iid;
static HashTable_AS *full_name2iid;
static char subset_names[2000000][20];
static char full_names[2000000][20];
static VA_TYPE(uint64) *subset_iid2name;
static VA_TYPE(uint64) *full_iid2name;
static OVL_Store_t  * my_store = NULL;


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


typedef enum {
  BEST_MEANS_THICKEST=1,
  BEST_MEANS_LOWEST_ERROR=2,
  BEST_MEANS_BEST_UPPER_BOUND=3
}BestMeasure;

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
    
				  





Long_Olap_Data_t better_olap(Long_Olap_Data_t a, Long_Olap_Data_t b, int startingFrg, int offAEnd, BestMeasure bestType, int frglen, int useCorrected){

  int aThick, bThick;
  int aMis, bMis;
  int aLen, bLen;

  Long_Olap_Data_t copy_a =a;
  Long_Olap_Data_t copy_b =b;

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


int best_overlap_off_end(int id, int offAEnd,BestMeasure bestType,Long_Olap_Data_t *bestovl,double erate,int useCorrected,int skipContaining,int minlen, int frglen){
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


    // if it's off the correct end ...
    if( (offAEnd ? (olap . a_hang < 0) : (olap.b_hang > 0))){

      //    print_olap(olap);

      if(goodOlap==0){
	goodOlap=1;
	bestolap=olap; 
      } else {
	bestolap=better_olap(bestolap,olap, id, offAEnd,bestType,frglen,useCorrected);
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


void finished_with_ovlStore(void){

  Free_OVL_Store (my_store);

}


void usage(char *pgm){
	fprintf (stderr, 
		 "USAGE:  %s -f <FragStoreName> -n <startingFrgNum> -o <full_ovlStore> [-C] [-e <erate cutoff>] [-E] [-m <minlen>] [-Q]\n"
		 "\t-n startingFrgNum = fragment to walk out from\n"
		 "\t-C specifies to not use containing fragments, i.e. only dovetails\n"
		 "\t-e specifies the maximum mismatch rate (as a fraction, i.e. .01 means one percent)\n"
		 "\t-E specifies that the corrected error rate rather than the original is to be used\n"
		 "\t-Q specifies that the lowest error rate overlap will be used (in place of thickest)\n"
		 "\t-B specifies that the overlap with the lowest upper bound on the mismatch rate will be used\n",
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
    while (!errflg && ((ch = getopt(argc, argv,"BCe:Ef:m:n:o:Q")) != EOF)){
      switch(ch) {
      case 'B':
	bestType=BEST_MEANS_BEST_UPPER_BOUND;
	break;
      case 'C':
	skipContaining = 0;
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

  setup_ovlStore(full_ovlPath);

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
    stillGoing = best_overlap_off_end(currFrg, Aend, bestType,&o,maxError,useCorrectedErate,skipContaining,minlen,frglen);
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

    stillGoing = best_overlap_off_end(currFrg, Aend, bestType,&o,maxError,useCorrectedErate,skipContaining,minlen,frglen);
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

