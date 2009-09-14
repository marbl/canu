
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

const char *mainid = "$Id: frgs2clones.c,v 1.37 2009-09-14 16:09:05 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>


#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "SYS_UIDclient.h"
#include "AS_UTL_reverseComplement.h"
#include "MultiAlignment_CNS.h"

#define MAXSEQLEN 20000


int
main( int argc, char *argv[]) {
  char *inputPath;
  char *prefix;

  int setGatekeeperStore = FALSE;
  int fragIID,mateIID;
  AS_UID fragUID,mateUID;
  char GKP_Store_Name[2000];
  gkStore *gkpStore;
  char *seq1,*seq2,*qul1,*qul2,*clear1,*clear2;
  uint clr_bgn1,clr_end1;
  uint clr_bgn2,clr_end2;
  int alloclen1=5000;
  int alloclen2=5000;
  int len1,len2,lastfrg;
  gkFragment fsread;
  gkFragment fsmate;
  uint64 UIDstart = 1230000;
  UIDserver   *uids              = NULL;

  Overlap *ovl;
  IntUnitigMesg ium;
  IntMultiPos the_imps[2];
  uint64  mergeUid;
  char seq[MAXSEQLEN], qlt[MAXSEQLEN];
  int clr_bgn,clr_end;
  VA_TYPE(int32) *deltas=CreateVA_int32(1);
  VA_TYPE(char) *sequence=CreateVA_char(200000);
  VA_TYPE(char) *quality=CreateVA_char(200000);
  int runConsensus=0;
  int Ngaps=0;

  //  setbuf(stdout,NULL);

  argc = AS_configure(argc, argv);

  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "g:NUC")) != EOF)){
      switch(ch) {
        case 'C':
          runConsensus=1;
          break;
        case 'g':
          strcpy( GKP_Store_Name, argv[optind - 1]);
          setGatekeeperStore = TRUE;
          break;
        case 'N':
          Ngaps=1;
          break;
        case 'U':
          UIDstart = 0;
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if((setGatekeeperStore == 0) || errflg>0)
      {
	fprintf(stderr,"* argc = %d optind = %d setGatekeeperStore = %d\n",
		argc, optind,setGatekeeperStore);
	fprintf (stderr, "USAGE:  %s -g <GatekeeperStoreName> [-U] [-C]\n",argv[0]);
	fprintf (stderr, "\t-U uses real UIDs\n");
	fprintf (stderr, "\t-C computes a consensus rather than splicing fragment seqs (slower, but better?)\n");
	exit (1);
      }

  }

  gkpStore = new gkStore(GKP_Store_Name, FALSE, FALSE);

  //  seq1=(char*)safe_malloc(sizeof(char)*alloclen1);
  //  qul1=(char*)safe_malloc(sizeof(char)*alloclen1);
  clear1=(char*)safe_malloc(sizeof(char)*alloclen1);
  //  assert(seq1!=NULL);
  //  assert(qul1!=NULL);
  assert(clear1!=NULL);
  //  seq2=(char*)safe_malloc(sizeof(char)*alloclen2);
  //  qul2=(char*)safe_malloc(sizeof(char)*alloclen2);
  clear2=(char*)safe_malloc(sizeof(char)*alloclen2);
  //  assert(seq2!=NULL);
  //  assert(qul2!=NULL);
  assert(clear2!=NULL);


  /*************************/
  // Set up UID server stuff
  /*************************/
  uids = UIDserverInitialize(256, UIDstart);
  assert(uids!=NULL);

  /*************************/
  // over all fragments, check for overlap with (previously unseen) mate
  /*************************/

  lastfrg = gkpStore->gkStore_getNumFragments () ;
  for (fragIID = 1; fragIID <= lastfrg; fragIID++){

    /*************************/
    // get the fragment
    /*************************/

    gkpStore->gkStore_getFragment(fragIID, &fsread, GKFRAGMENT_QLT);

    if (fsread.gkFragment_getIsDeleted())
      continue;

    fragUID = fsread.gkFragment_getReadUID();

    fsread.gkFragment_getClearRegion(clr_bgn1, clr_end1);

    while(alloclen1<=fsread.gkFragment_getSequenceLength()){
      alloclen1*=2;
      clear1=(char*)safe_realloc(clear1,alloclen1*sizeof(char));
    }
    seq1 = fsread.gkFragment_getSequence();
    qul1 = fsread.gkFragment_getQuality();
    strcpy(clear1,seq1+clr_bgn1);
    len1=clr_end1-clr_bgn1;
    clear1[len1]='\0';


    /*************************/
    // check for an appropriate mate
    /*************************/

    mateIID = fsread.gkFragment_getMateIID();

    if(mateIID == 0){

      // if no mate (or multiple mates), output fragment itself
      printf(">%s\n%s\n",AS_UID_toString(fragUID),clear1);

    } else { // there are links

      /*************************/
      // get (clear) sequence of mate
      /*************************/

      gkpStore->gkStore_getFragment(mateIID, &fsmate, GKFRAGMENT_QLT);

      if ((mateIID < fragIID) && (fsmate.gkFragment_getIsDeleted() == 0))
        continue;

      mateUID = fsmate.gkFragment_getReadUID();

      fsmate.gkFragment_getClearRegion(clr_bgn2, clr_end2);
      while(alloclen2<=fsmate.gkFragment_getSequenceLength()){
	alloclen2*=2;
	clear2=(char*)safe_realloc(clear2,alloclen2*sizeof(char));
      }
      seq2 = fsmate.gkFragment_getSequence();
      qul2 = fsmate.gkFragment_getQuality();
      strcpy(clear2,seq2+clr_bgn2);
      len2=clr_end2-clr_bgn2;
      clear2[len2]='\0';

      if(fsmate.gkFragment_getIsDeleted()){
	// if no mate (or multiple mates), output fragment itself
	printf(">%s\n%s\n",AS_UID_toString(mateUID),clear2);
	continue;
      }

      /*********************************************/
      // Create a UID for the clone
      /*********************************************/

      mergeUid = getUID(uids);

      /*************************/
      // check for an overlap
      /*************************/

      ovl = Local_Overlap_AS_forCNS(clear1, clear2,
                                    -len2, len1,
                                    len2,  len1,
                                    1,
                                    .06,
                                    1e-6,
                                    40,
                                    AS_FIND_LOCAL_ALIGN_NO_TRACE);

      if(ovl==NULL||
	 ( (ovl->begpos<0||ovl->endpos<0) &&
	   ((len1+len2)-abs(ovl->begpos)-abs(ovl->endpos))/2<100)
	 ){

	// if they don't overlap reasonably,

	if(Ngaps){
	  printf(">" F_S64 " from mated fragments %s and %s\n%sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN%s\n",
		 mergeUid,AS_UID_toString(fragUID),AS_UID_toString(mateUID),clear1,clear2);
	} else {
	  // output two sequences, but with a clone UID plus "a" or "b"

	  printf(">" F_S64 "a (fragment %s)\n%s\n>" F_S64 "b (fragment %s)\n%s\n",
		 mergeUid,AS_UID_toString(fragUID),clear1,mergeUid,AS_UID_toString(mateUID),clear2);
	}

      } else { // there is an overlap

	if(runConsensus){
	  /*************************/
	  // create a unitig
	  /*************************/

	  ium.consensus = "";
	  ium.quality = "";
	  ium.iaccession = 0;
	  ium.forced = FALSE;
	  ium.coverage_stat = 10;
          ium.microhet_prob = 1.01;
	  ium.status = (UnitigStatus)'U';
	  ium.num_frags = 2;
	  ium.f_list = &(the_imps[0]);
	  {
	    the_imps[0].type = (FragType)'R';
	    the_imps[0].ident = fragIID;
	    the_imps[0].contained = 0;
	    the_imps[0].position.bgn = (ovl->begpos >= 0 ) ? 0 : -ovl->begpos;
	    the_imps[0].position.end = (ovl->begpos >= 0 ) ? len1 : len1 - ovl->begpos;
	    the_imps[0].delta_length = 0;
	    the_imps[0].delta        = NULL;
	    the_imps[1].type = (FragType)'R';
	    the_imps[1].ident = mateIID;
	    the_imps[1].contained = 0;
	    // due to inversion of mate, note the following swap of end and beg
	    the_imps[1].position.end = (ovl->begpos >= 0) ? ovl->begpos : 0;
	    the_imps[1].position.bgn = (ovl->begpos >= 0) ? ovl->begpos + len2 : len2;
	    the_imps[1].delta_length = 0;
	    the_imps[1].delta        = NULL;
	  }
	  ium.length = ( ium.f_list[0].position.end >  ium.f_list[1].position.bgn ) ?
	    ium.f_list[0].position.end :  ium.f_list[1].position.bgn;


	  /*************************/
	  // run consensus on unitig
	  /*************************/
	  {
	    MultiAlignT *ma;
	    CNS_PrintKey printwhat=CNS_STATS_ONLY;
	    int i,j,len;
	    char *s,*q;

	    if(ovl->begpos<0){
	      allow_neg_hang=1;
	    }

	    if (MultiAlignUnitig(&ium,gkpStore,sequence,quality,deltas,printwhat, NULL) == 0 ) {
	      fprintf(stderr,"MultiAlignUnitig failed for overlap of fragments %d and %d\n",
                      fragIID,mateIID);
	      assert(FALSE);
	    }

	    if(ovl->begpos<0){
	      allow_neg_hang=1;
	    }


	    len = GetNumVA_char(sequence)-1;
	    assert(len<MAXSEQLEN);
	    j=0;
	    s = Getchar(sequence,0);
	    q = Getchar(quality,0);
	    for(i=0;i<len;i++){
	      if(s[i]!='-'){
		seq[j] = s[i];
		qlt[j] = q[i];
		j++;
	      }
	    }
	    seq[j]='\0';
	    qlt[j]='\0';
	    clr_bgn=0;
	    clr_end=j;
	  }

	} else { // do not run consensus

	  int into1 = len1;
	  int into2 = strlen(clear2)-ovl->endpos;

	  assert(len1+len2+50<MAXSEQLEN);
	  strcpy(seq,clear1);
	  reverseComplementSequence(clear2, strlen(clear2));
	  strcpy(seq+into1,clear2+into2);
	  assert(strlen(seq)<MAXSEQLEN);

	}

	printf(">" F_S64,mergeUid);
	printf(" merged sequence of mated fragments %s and %s\n",AS_UID_toString(fragUID),AS_UID_toString(mateUID));
	printf("%s\n",seq);

      }
    }
  }
  exit(0);
}
