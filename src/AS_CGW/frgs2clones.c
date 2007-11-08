
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
#include "AS_PER_gkpStore.h"
#include "SYS_UIDclient.h"
#include "MultiAlignment_CNS.h"

#define MAXSEQLEN 20000


void RevCompl(char *seq, char *qul)
{ static char WCinvert[256];
 static int Firstime = 1;

 if (Firstime)          /* Setup complementation array */
   { 
     int i;
     Firstime = 0;
     for(i = 0; i < 256;i++){
       WCinvert[i] = '?';
     }
     WCinvert['a'] = 't';
     WCinvert['c'] = 'g';
     WCinvert['g'] = 'c';
     WCinvert['t'] = 'a';
     WCinvert['n'] = 'n';
     WCinvert['A'] = 'T';
     WCinvert['C'] = 'G';
     WCinvert['G'] = 'C';
     WCinvert['T'] = 'A';
     WCinvert['N'] = 'N';
     WCinvert['-'] = '-'; // added this to enable alignment of gapped consensi
   }
      
 { int len;                    /* Complement and reverse sequence */
 len = strlen(seq);

 { register char *s, *t;
 int c;

 s = seq;
 t = seq + (len-1);
 while (s < t)
   { // Sanity Check!
     assert(WCinvert[(int) *t] != '?' &&
            WCinvert[(int) *s] != '?');

     c = *s;
     *s++ = WCinvert[(int) *t];
     *t-- = WCinvert[c];
   }
 if (s == t)
   *s = WCinvert[(int) *s];
 }

 if (qul != NULL)
   { register char *s, *t;   /* Reverse quality value array */
   int c;
    
   s = qul;
   t = qul + (len-1);
   while (s < t)
     { c = *s;
     *s++ = *t;
     *t-- = c;
     }
   }
 }
}



int main( int argc, char *argv[])
{
  char *inputPath;
  char *prefix;

  int setGatekeeperStore = FALSE;
  int fragIID,mateIID;
  AS_UID fragUID,mateUID;
  char GKP_Store_Name[2000];
  GateKeeperStore *gkpStore;
  GateKeeperFragmentRecord gkpFrag,gkpMate;
  uint64 uid, mateuid;
  char *seq1,*seq2,*qul1,*qul2,*clear1,*clear2;
  uint clr_bgn1,clr_end1;
  uint clr_bgn2,clr_end2;
  int alloclen1=5000;
  int alloclen2=5000;
  int len1,len2,lastfrg;
  fragRecord fsread;
  fragRecord fsmate;
  uint64 UIDstart = 1230000;
  UIDserver   *uids              = NULL;

  Overlap *ovl;
  IntUnitigMesg ium;
  IntMultiPos the_imps[2];
  AS_UID  mergeUid;
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
	exit (EXIT_FAILURE);
      }

  }

  gkpStore = openGateKeeperStore(GKP_Store_Name, FALSE);

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
  
  lastfrg = getLastElemFragStore (gkpStore) ;
  for (fragIID = 1; fragIID <= lastfrg; fragIID++){

    /*************************/
    // get the fragment
    /*************************/

    getGateKeeperFragment(gkpStore,fragIID,&gkpFrag);
    if(gkpFrag.deleted)continue;

    fragUID = gkpFrag.readUID;

      //    if(getFrag(gkpStore,fragIID,fsread,FRAG_S_ALL)!=0){
      //      fprintf(stderr,"Couldn't get fragment from gkpStore for iid %d\n",fragIID);
      //      assert(0);
      //    }
    getFrag(gkpStore,fragIID,&fsread,FRAG_S_ALL);
    //    getClearRegion_ReadStruct(fsread, &clr_bgn1,&clr_end1, READSTRUCT_LATEST);
    clr_bgn1 = getFragRecordClearRegionBegin(&fsread, AS_READ_CLEAR_LATEST);
    clr_end1 = getFragRecordClearRegionEnd  (&fsread, AS_READ_CLEAR_LATEST);

    while(alloclen1<=getFragRecordSequenceLength(&fsread)){
      alloclen1*=2;
      clear1=(char*)safe_realloc(clear1,alloclen1*sizeof(char));
    }
    seq1 = getFragRecordSequence(&fsread);
    qul1 = getFragRecordQuality(&fsread);
    strcpy(clear1,seq1+clr_bgn1);
    len1=clr_end1-clr_bgn1;
    clear1[len1]='\0';


    /*************************/
    // check for an appropriate mate
    /*************************/

    if(gkpFrag.mateIID == 0){

      // if no mate (or multiple mates), output fragment itself
      printf(">" F_S64 "\n%s\n",fragUID,clear1);

    } else { // there are links
      mateIID = gkpFrag.mateIID;

      /*************************/
      // get (clear) sequence of mate
      /*************************/

      getGateKeeperFragment(gkpStore,mateIID,&gkpFrag);
      mateUID = gkpFrag.readUID;
	
      if(mateIID<fragIID&&gkpFrag.deleted!=1)continue;

      //      if(getFrag(gkpStore,mateIID,fsmate,FRAG_S_ALL)!=0){
      //	fprintf(stderr,"Couldn't get fragment from gkpStore for iid %d\n",mateIID);
      //	assert(0);
      //      }
      getFrag(gkpStore,mateIID,&fsmate,FRAG_S_ALL);
      //      getClearRegion_ReadStruct(fsmate, &clr_bgn2,&clr_end2, READSTRUCT_LATEST);
      clr_bgn2 = getFragRecordClearRegionBegin(&fsmate, AS_READ_CLEAR_LATEST);
      clr_end2 = getFragRecordClearRegionEnd  (&fsmate, AS_READ_CLEAR_LATEST);
      while(alloclen2<=getFragRecordSequenceLength(&fsmate)){
	alloclen2*=2;
	clear2=(char*)safe_realloc(clear2,alloclen2*sizeof(char));
      }
      seq2 = getFragRecordSequence(&fsmate);
      qul2 = getFragRecordQuality(&fsmate);
      strcpy(clear2,seq2+clr_bgn2);
      len2=clr_end2-clr_bgn2;
      clear2[len2]='\0';

      if(gkpFrag.deleted){
	// if no mate (or multiple mates), output fragment itself
	printf(">" F_S64 "\n%s\n",mateUID,clear2);
	continue;
      }

      /*********************************************/
      // Create a UID for the clone
      /*********************************************/

      mergeUid = AS_UID_fromInteger(getUID(uids));

      /*************************/
      // check for an overlap
      /*************************/

      ovl = Local_Overlap_AS_forCNS(clear1, clear2, -len2,len1,1,.06,1e-6,40,AS_FIND_LOCAL_ALIGN_NO_TRACE);

      if(ovl==NULL||
	 ( (ovl->begpos<0||ovl->endpos<0) &&
	   ((len1+len2)-abs(ovl->begpos)-abs(ovl->endpos))/2<100)
	 ){
	
	// if they don't overlap reasonably,

	if(Ngaps){
	  printf(">" F_S64 " from mated fragments " F_S64 " and " F_S64 "\n%sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN%s\n",
		 mergeUid,fragUID,mateUID,clear1,clear2);
	} else {
	  // output two sequences, but with a clone UID plus "a" or "b"
	
	  printf(">" F_S64 "a (fragment " F_S64 ")\n%s\n>" F_S64 "b (fragment " F_S64 ")\n%s\n",
		 mergeUid,fragUID,clear1,mergeUid,mateUID,clear2);
	}

      } else { // there is an overlap

	if(runConsensus){
	  /*************************/
	  // create a unitig
	  /*************************/

	  ium.consensus = "";
	  ium.quality = "";
	  ium.iaccession = 0;
# ifdef AS_ENABLE_SOURCE
	  ium.source = "";
#endif
	  ium.forced = FALSE;
	  ium.coverage_stat = 10; 
	  ium.status = 'U';
	  ium.num_frags = 2;
	  ium.f_list = &(the_imps[0]);
	  {
	    the_imps[0].type = 'R';
	    the_imps[0].ident = fragIID;
	    the_imps[0].contained = 0;
	    the_imps[0].sourceInt = -1;
	    the_imps[0].position.bgn = (ovl->begpos >= 0 ) ? 0 : -ovl->begpos;
	    the_imps[0].position.end = (ovl->begpos >= 0 ) ? len1 : len1 - ovl->begpos;
	    the_imps[0].delta_length = 0;
	    the_imps[0].delta        = NULL;
	    the_imps[1].type = 'R';
	    the_imps[1].ident = mateIID;
	    the_imps[1].contained = 0;
	    the_imps[1].sourceInt = -1;
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
	    int printwhat=CNS_STATS_ONLY;
	    int do_rez=0;
	    int i,j,len;
	    char *s,*q;

	    Overlap *(*COMPARE_FUNC)(COMPARE_ARGS)=Local_Overlap_AS_forCNS;
	    //Overlap *(*COMPARE_FUNC)(COMPARE_ARGS)=DP_Compare;

	    
	    if(ovl->begpos<0){
	      allow_neg_hang=1;
	    }

	    //      fprintf(stderr,"Doing the multialignment\n");

	    if (MultiAlignUnitig(&ium,gkpStore,sequence,quality,deltas,printwhat,do_rez,
                                 COMPARE_FUNC, NULL)==-1 ) {
	      fprintf(stderr,"MultiAlignUnitig failed for overlap of fragments %d and %d\n",
                      fragIID,mateIID);
	      assert(FALSE);
	    }

	    //      fprintf(stderr,"Done with the multialignment\n");

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
	  RevCompl(clear2,NULL);
	  strcpy(seq+into1,clear2+into2);
	  assert(strlen(seq)<MAXSEQLEN);

	}

	printf(">" F_S64,mergeUid);
	printf(" merged sequence of mated fragments " F_S64 " and " F_S64 "\n",fragUID,mateUID);
	printf("%s\n",seq);

      }
    }
  }
  exit(0);
}
