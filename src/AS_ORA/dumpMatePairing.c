/*********************************************************************
 *		  Confidential -- Do Not Distribute                  *
 *	    Copyright © 2005 The J. Craig Venter Institute           *
 *			 All rights Reserved.                        *
 *                                                                   *
 *********************************************************************/


static char CM_ID[] = "$Id";


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
#include "AS_PER_fragStore.h"



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
  int len1,len2,lastfrg;
  ReadStructp fsread=new_ReadStruct();
  ReadStructp fsmate=new_ReadStruct();
  int printUID=0;
  int firstIID =1;
  int lastIID = 0;

  //  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "f:g:U")) != EOF)){
      switch(ch) {
      case 'f':
	strcpy( Frag_Store_Name, argv[optind - 1]);
	setFragStore = TRUE;
	break;
      case 'g':
	strcpy( GKP_Store_Name, argv[optind - 1]);
	setGatekeeperStore = TRUE;
	break;	  
      case 'U':
	printUID=1;
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
	fprintf (stderr, "USAGE:  %s -f <FragStoreName> -g <GatekeeperStoreName> [-U]\n",argv[0]);
	fprintf (stderr, "\t-n specifies that even mates with negative hangs should be merged\n");
	fprintf (stderr, "\t-U prints UIDs instead of IIDs\n");
	exit (EXIT_FAILURE);
      }

  }

  assert(existsFragStore(Frag_Store_Name) == TRUE);
  frgStore = openFragStore(Frag_Store_Name,"r");

  InitGateKeeperStore(&gkpStore,GKP_Store_Name);
  assert(TestOpenGateKeeperStore(&gkpStore) == TRUE);
  OpenReadOnlyGateKeeperStore(&gkpStore);

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

    rv2 = getGateKeeperFragmentStore(gkpStore.frgStore,mateIID,&gkpFrag);
    assert(rv2==0);
    mateUID = gkpFrag.readUID;
    if(printUID){
      printf(F_UID "\t" F_UID "\n",fragUID,mateUID);
    } else {
      printf(F_IID "\t" F_IID "\n",fragIID,mateIID);
    }
  }
  exit(0);
}
