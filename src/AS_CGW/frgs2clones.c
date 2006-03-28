static char CM_ID[] = "$Id:";


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

#include "cds.h"
#include "AS_global.h"

#include "AS_PER_gkpStore.h"

#include "SYS_UIDcommon.h"
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
  int runConsensus=0;
  int Ngaps=0;

  //  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "f:g:NUC")) != EOF)){
      switch(ch) {
      case 'C':
	runConsensus=1;
	break;
      case 'f':
	strcpy( Frag_Store_Name, argv[optind - 1]);
	setFragStore = TRUE;
	break;
      case 'g':
	strcpy( GKP_Store_Name, argv[optind - 1]);
	setGatekeeperStore = TRUE;
	break;	  
      case 'N':
	Ngaps=1;
	break;
      case 'U':
	realUID=1;
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
	fprintf (stderr, "USAGE:  %s -f <FragStoreName> -g <GatekeeperStoreName> [-U] [-C]\n",argv[0]);
	fprintf (stderr, "\t-U uses real UIDs\n");
	fprintf (stderr, "\t-C computes a consensus rather than splicing fragment seqs (slower, but better?)\n");
	exit (EXIT_FAILURE);
      }

  }

  assert(existsFragStore(Frag_Store_Name) == TRUE);
  frgStore = openFragStore(Frag_Store_Name,"r");

  InitGateKeeperStore(&gkpStore,GKP_Store_Name);
  assert(TestOpenGateKeeperStore(&gkpStore) == TRUE);
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
  // Set up UID server stuff
  /*************************/
  {

    {
      int32 blockSize = 300;
      int32  uidStatus;
      CDS_UID_t interval_UID[4];
      if(firstUID){
	firstUID=0;
	set_start_uid(UIDstart); /* used if readUID == FALSE */
	get_uids(blockSize,interval_UID,realUID);
      }

      uidStatus = get_next_uid(&mergeUid,realUID);
      if( uidStatus != UID_CODE_OK )
	{
	  uidStatus = get_uids(blockSize,interval_UID,realUID);
	  get_next_uid(&mergeUid,realUID);
	}	  
      if( UID_CODE_OK != uidStatus )
	{ 
          fprintf(stderr, "Could not get UID \n");
          assert(0);
	}
    }
  }

  /*************************/
  // over all fragments, check for overlap with (previously unseen) mate
  /*************************/
  
  lastfrg = getLastElemFragStore (frgStore) ;
  for (fragIID = 1; fragIID <= lastfrg; fragIID++){
    int rv1,rv2;

    /*************************/
    // get the fragment
    /*************************/

    rv1 = getGateKeeperFragmentStore(gkpStore.frgStore,fragIID,&gkpFrag);
    if(gkpFrag.deleted)continue;
    assert(rv1==0);
    fragUID = gkpFrag.readUID;

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


    /*************************/
    // check for an appropriate mate
    /*************************/

    if(gkpFrag.numLinks!=1){

      // if no mate (or multiple mates), output fragment itself
      printf(">" F_S64 "\n%s\n",fragUID,clear1);

    } else { // there are links
      GateKeeperLinkRecordIterator iterator;
      GateKeeperLinkRecord link;
      CreateGateKeeperLinkRecordIterator(gkpStore.lnkStore, gkpFrag.linkHead,fragIID, &iterator);
      while(NextGateKeeperLinkRecordIterator(&iterator, &link))
	mateIID = (link.frag1 == fragIID) ? link.frag2 : link.frag1;

      /*************************/
      // get (clear) sequence of mate
      /*************************/


      rv2 = getGateKeeperFragmentStore(gkpStore.frgStore,mateIID,&gkpFrag);
      assert(rv2==0);
      mateUID = gkpFrag.readUID;
	
      if(mateIID<fragIID&&gkpFrag.deleted!=1)continue;

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

      if(gkpFrag.deleted){
	// if no mate (or multiple mates), output fragment itself
	printf(">" F_S64 "\n%s\n",mateUID,clear2);
	continue;
      }

      /*********************************************/
      // Create a UID for the clone
      /*********************************************/

      {
	int32 blockSize = 300;
	int32  uidStatus;
	CDS_UID_t interval_UID[4];
	if(firstUID){
	  firstUID=0;
	  set_start_uid(UIDstart); /* used if readUID == FALSE */
	  get_uids(blockSize,interval_UID,realUID);
	}

	uidStatus = get_next_uid(&mergeUid,realUID);
	if( uidStatus != UID_CODE_OK )
	  {
	    uidStatus = get_uids(blockSize,interval_UID,realUID);
	    get_next_uid(&mergeUid,realUID);
	  }	  
	if( UID_CODE_OK != uidStatus )
	  { 
            fprintf(stderr, "Could not get UID \n");
            assert(0);
	  }
      }


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

	    if (MultiAlignUnitig(&ium,frgStore,sequence,quality,deltas,printwhat,do_rez,
                COMPARE_FUNC, NULL)==-1 ) {
	      fprintf(stderr,"MultiAlignUnitig failed for overlap of fragments %d and %d\n",
                  fragIID,mateIID);
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
