
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

static char CM_ID[] = "$Id: combineMates.c,v 1.10 2006-10-08 08:47:39 brianwalenz Exp $";


/*********************************************************************/

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

#include "AS_PER_gkpStore.h"

#include "SYS_UIDcommon.h"
#include "SYS_UIDclient.h"

#include "MultiAlignment_CNS.h"

#define MAXSEQLEN 20000

/* Output text field item with 3-code field-name "tag". */

static void copyPutText(FILE *fout, const char * const tag, 
                        char * text, const int format)
{
  // Note that the data of "text" is modified!!!
  int i, len;
  fprintf(fout,"%s\n",tag);
  if(text != NULL){
    len = strlen(text);
    if (format)
      { for (i = 0; i < len; i += 70)
        { fprintf(fout,"%.*s\n",70,text);
	text += 70;
        }
      fprintf(fout,".\n");
      } else{ 
	if (text[len-1] == '\n')     /* Strip trailing new line if prez. */
	  text[len-1] = '\0';
	fprintf(fout,"%s\n.\n", text);
      }
  }else{
    fprintf(fout,"\n.\n");
  }
}


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

  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int fragIID,mateIID;
  Fragment_ID fragUID,mateUID;
  char Frag_Store_Name[2000];
  char GKP_Store_Name[2000];
  FragStoreHandle frgStore = 0;
  GateKeeperStore gkpStore;
  GateKeeperFragmentRecord gkpFrag,gkpMate;
  uint64 uid, mateuid;
  char *seq1,*seq2,*qul1,*qul2,*clear1,*clear2;
  uint clr_bgn1,clr_end1;
  uint clr_bgn2,clr_end2;
  int alloclen1=5000;
  int alloclen2=5000;
  int len1,len2,lastfrg;
  ReadStructp fsread=new_ReadStruct();
  ReadStructp fsmate=new_ReadStruct();
  int realUID=0;
  int UIDstart=1230000;
  int firstUID=1;
  CDS_UID_t       interval_UID[4];
  Overlap *ovl;
  IntUnitigMesg ium;
  IntMultiPos the_imps[2];
  CDS_UID_t mergeUid;
  char seq[MAXSEQLEN], qlt[MAXSEQLEN];
  int clr_bgn,clr_end;
  VA_TYPE(int32) *deltas=CreateVA_int32(1);
  VA_TYPE(char) *sequence=CreateVA_char(200000);
  VA_TYPE(char) *quality=CreateVA_char(200000);
  int firstIID =1;
  int lastIID = 0;
  int useNegHangs=0;
  int verbose=0;

  //  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "f:g:Us:e:nV")) != EOF)){
      switch(ch) {
        case 'e':
          lastIID = atoi(argv[optind-1]);
          assert(lastIID>=1);
          break;
        case 'f':
          strcpy( Frag_Store_Name, argv[optind - 1]);
          setFragStore = TRUE;
          break;
        case 'g':
          strcpy( GKP_Store_Name, argv[optind - 1]);
          setGatekeeperStore = TRUE;
          break;	  
        case 'n':
          useNegHangs=1;
          break;
        case 's':
          firstIID = atoi(argv[optind-1]);
          assert(firstIID>=1);
          break;
        case 'U':
          realUID=1;
          break;
        case 'V':
          verbose=1;
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if((setFragStore == 0) || (setGatekeeperStore == 0) || errflg>0)
      {
	fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d\n",
		argc, optind, setFragStore,setGatekeeperStore);
	fprintf (stderr, "USAGE:  %s -f <FragStoreName> -g <GatekeeperStoreName> [-U] [-V] [-s <firstIID>] [-e <lastIID>] [-n]\n",argv[0]);
	fprintf (stderr, "\t-U causes generation of real UIDs\n");
	fprintf (stderr, "\t-V causes info about used overlaps to be printed to stderr\n");
	fprintf (stderr, "\t-s specifies first IID to examine (default = 1)\n");
	fprintf (stderr, "\t-e specifies last IID to examine (default = last frag in stores)\n");
	fprintf (stderr, "\t-n specifies that even mates with negative hangs should be merged\n");
	exit (EXIT_FAILURE);
      }

  }

  assert(existsFragStore(Frag_Store_Name) == TRUE);
  frgStore = openFragStore(Frag_Store_Name,"r");

  InitGateKeeperStore(&gkpStore,GKP_Store_Name);
  assert(TestOpenReadOnlyGateKeeperStore(&gkpStore) == TRUE);
  OpenReadOnlyGateKeeperStore(&gkpStore);

  seq1=(char*)malloc(sizeof(char)*alloclen1);
  qul1=(char*)malloc(sizeof(char)*alloclen1);
  clear1=(char*)malloc(sizeof(char)*alloclen1);
  assert(seq1!=NULL);
  assert(qul1!=NULL);
  assert(clear1!=NULL);
  seq2=(char*)malloc(sizeof(char)*alloclen2);
  qul2=(char*)malloc(sizeof(char)*alloclen2);
  clear2=(char*)malloc(sizeof(char)*alloclen2);
  assert(seq2!=NULL);
  assert(qul2!=NULL);
  assert(clear2!=NULL);

  /*************************/
  // Construct a BAT message
  /*************************/
  {

    /*************************/
    // Get a UID to use
    /*************************/
    {
      int32 blockSize = 300;
      int32  uidStatus;
      CDS_UID_t interval_UID[4];
      if(firstUID){
	firstUID=0;
	set_start_uid(UIDstart); /* used if realUID == FALSE */
	get_uids(blockSize,interval_UID,realUID);
      }

      uidStatus = get_next_uid(&mergeUid,realUID);
      if( uidStatus != UID_CODE_OK )
	{
	  get_uids(blockSize,interval_UID,realUID);
	  uidStatus = get_next_uid(&mergeUid,realUID);
	}	  
      if( UID_CODE_OK != uidStatus )
	{ 
          fprintf(stderr, "Could not get UID \n");
          assert(0);
	}
    }
    /***********************/
    // Print a BAT message
    /***********************/
    printf("{BAT\n");
    printf("bna:(Batch name)\n");
    printf("crt:" F_TIME_T "\n",time(NULL));
    printf("acc:" F_UID "\n",mergeUid);
    printf("com:\nCreated by %s\n.\n",__FILE__);
    printf("}\n");
  }

  /*************************/
  // over all fragments, check for overlap with (previously unseen) mate
  /*************************/

  if(lastIID==0){
    lastfrg = getLastElemFragStore (frgStore) ;
  } else {
    lastfrg = lastIID;
    assert(lastfrg<=getLastElemFragStore (frgStore) );
  }
  assert(firstIID<=lastfrg);
  for (fragIID = firstIID; fragIID <= lastfrg; fragIID++){
    int rv1,rv2;

    /*************************/
    // get the fragment
    /*************************/


    //    fprintf(stderr,"Working on frgIID %d\n",fragIID);
    rv1 = getGateKeeperFragmentStore(gkpStore.frgStore,fragIID,&gkpFrag);

    assert(rv1==0);
    fragUID = gkpFrag.readUID;

    /*************************/
    // check for an appropriate mate
    /*************************/

    if(gkpFrag.numLinks!=1){
      continue;
    }
    {
      GateKeeperLinkRecordIterator iterator;
      GateKeeperLinkRecord link;
      CreateGateKeeperLinkRecordIterator(gkpStore.lnkStore, gkpFrag.linkHead,fragIID, &iterator);
      while(NextGateKeeperLinkRecordIterator(&iterator, &link))
	mateIID = (link.frag1 == fragIID) ? link.frag2 : link.frag1;
      if(mateIID<fragIID)continue;
    }

    /*************************/
    // get (clear) sequences of fragments
    /*************************/

    if(getFragStore(frgStore,fragIID,FRAG_S_ALL,fsread)!=0){
      fprintf(stderr,"Couldn't get fragment from frgStore for iid %d\n",fragIID);
      assert(0);
    }
    getClearRegion_ReadStruct(fsread, &clr_bgn1,&clr_end1, READSTRUCT_LATEST);
    while(getSequence_ReadStruct(fsread,seq1,qul1,alloclen1)!=0){
      alloclen1*=2;
      seq1=(char*)realloc(seq1,alloclen1*sizeof(char));
      qul1=(char*)realloc(qul1,alloclen1*sizeof(char));
      clear1=(char*)realloc(clear1,alloclen1*sizeof(char));
    }
    strcpy(clear1,seq1+clr_bgn1);
    len1=clr_end1-clr_bgn1;
    clear1[len1]='\0';
    //    fprintf(stderr, F_UID " %d : %s\n", 	    fragUID,len1,clear1);

    rv2 = getGateKeeperFragmentStore(gkpStore.frgStore,mateIID,&gkpFrag);
    assert(rv2==0);
    mateUID = gkpFrag.readUID;
    if(getFragStore(frgStore,mateIID,FRAG_S_ALL,fsmate)!=0){
      fprintf(stderr,"Couldn't get fragment from frgStore for iid %d\n",mateIID);
      assert(0);
    }
    getClearRegion_ReadStruct(fsmate, &clr_bgn2,&clr_end2, READSTRUCT_LATEST);
    while(getSequence_ReadStruct(fsmate,seq2,qul2,alloclen2)!=0){
      alloclen2*=2;
      seq2=(char*)realloc(seq2,alloclen2*sizeof(char));
      qul2=(char*)realloc(qul2,alloclen2*sizeof(char));
      clear2=(char*)realloc(clear2,alloclen2*sizeof(char));
    }
    strcpy(clear2,seq2+clr_bgn2);
    len2=clr_end2-clr_bgn2;
    clear2[len2]='\0';

    /*************************/
    // check for an overlap
    /*************************/

    //ovl = Local_Overlap_AS_forCNS(clear1, clear2, -len2,len1,1,.06,1e-6,40,AS_FIND_LOCAL_ALIGN_NO_TRACE);
    ovl = DP_Compare(clear1, clear2, -len2,len1,1,.06,1e-6,40,AS_FIND_LOCAL_ALIGN_NO_TRACE);

    if(ovl==NULL)continue;

    if(ovl->begpos<0||ovl->endpos<0){
      fprintf(stderr,"NEGATIVE HANG(S) ENCOUNTERED: [%d,%d] for fragments (" F_UID ",%d) and (" F_UID ",%d)\n",
	      ovl->begpos,ovl->endpos,
	      fragUID,fragIID,mateUID,mateIID);
      if(!useNegHangs){
	continue;
      }
    }

    if(verbose){
      fprintf(stderr,"Using overlap with length %d between %d and %d\n",ovl->length,fragIID,mateIID);
    }

#undef REALLY_USE_CONSENSUS
#ifdef REALLY_USE_CONSENSUS
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
    ium.a_branch_point = 0;
    ium.b_branch_point = 0;
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


      //      fprintf(stderr,"Doing the multialignment\n");

      if (MultiAlignUnitig(&ium,frgStore,sequence,quality,deltas,printwhat,do_rez,COMPARE_FUNC)==-1 ) {
	fprintf(stderr,"MultiAlignUnitig failed for overlap of fragments %d and %d\n",fragIID,mateIID);
	assert(FALSE);
      }

      //      fprintf(stderr,"Done with the multialignment\n");


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
#else
    {
      int into1 = clr_bgn1+len1;
      int into2 = strlen(seq2)-(ovl->endpos+clr_bgn2);
      strcpy(seq,seq1);
      strcpy(qlt,qul1);
      RevCompl(seq2,qul2);
      strcpy(seq+into1,seq2+into2);
      strcpy(qlt+into1,qul2+into2);
      clr_bgn=clr_bgn1;
      clr_end=clr_bgn1+len1+ovl->endpos;
    }
#endif

    /*********************************************/
    // Create a UID for the merged "fragment"
    /*********************************************/

    {
      int32 blockSize = 300;
      int32  uidStatus;
      CDS_UID_t interval_UID[4];
      if(firstUID){
	firstUID=0;
	set_start_uid(UIDstart); /* used if realUID == FALSE */
	get_uids(blockSize,interval_UID,realUID);
      }

      uidStatus = get_next_uid(&mergeUid,realUID);
      if( uidStatus != UID_CODE_OK )
	{
	  get_uids(blockSize,interval_UID,realUID);
	  uidStatus = get_next_uid(&mergeUid,realUID);
	}	  
      if( UID_CODE_OK != uidStatus )
	{ 
          fprintf(stderr, "Could not get UID \n");
          assert(0);
	}
    }


    /*************************/
    // generate messages to delete existing fragments and link
    /*************************/

    printf("{LKG\n");
    printf("act:D\n");
    printf("typ:M\n");
    printf("fg1:" F_UID "\n",fragUID);
    printf("fg2:" F_UID "\n",mateUID);
    printf("}\n");
    printf("{FRG\n");
    printf("act:D\n");
    printf("acc:" F_UID "\n",fragUID);
    printf("}\n");
    printf("{FRG\n");
    printf("act:D\n");
    printf("acc:" F_UID "\n",mateUID);
    printf("}\n");

#if 0 
    /*************************/
    // munge unitig into pseudo-fragment
    /*************************/

    /*************************/
    // output the fragment
    /*************************/
#else

    /*************************/
    // hack a new FRG message
    /*************************/

    
    printf("{FRG\n");
    printf("act:A\n");
    printf("acc:" F_S64 "\n",mergeUid);
    printf("typ:R\n");
    printf("src:\nmerged sequence of overlapping mated fragments " F_S64 " and " F_S64 "\n.\n",fragUID,mateUID);
    printf("etm:0\n");
    copyPutText(stdout,"seq:",seq,1);
    copyPutText(stdout,"qlt:",qlt,1);
    printf("clr:%d,%d\n",clr_bgn,clr_end);
    printf("}\n");

#endif


  }
  exit(0);
}
