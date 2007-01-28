
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
static char CM_ID[] = "$Id: dumpGatekeeperRecords.c,v 1.6 2007-01-28 21:52:24 brianwalenz Exp $";

/* Dump the gatekeeper stores for debug */

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
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_version.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"

int  nerrs = 0;   // Number of errors in current run
int maxerrs = 10; // Number of errors allowed before we punt

static MesgReader reader;
static MesgWriter writer;


int  main(int argc, char * argv [])

{
  int  summary;
  char *gatekeeperStorePath;
  GateKeeperStore gkpStore;
  CDS_UID_t uid;
  int type = AS_IID_FRG;

  summary = 0;
  /**************** Process Command Line Arguments *********************/
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "bf")) != EOF))
      switch(ch) {
      case 'f':
	type = AS_IID_FRG;
	break;
      case 'b':
	type = AS_IID_BTG;
	break;
      case '?':
	fprintf(stderr,"Unrecognized option -%c",optopt);
      default :
	errflg++;
      }

     
    if(argc - optind != 1 )
      {
	fprintf (stderr, "USAGE:  dumpGatekeeperRecords [-fb] <gatekeeperStorePath> < <input-file>\n");
	fprintf(stderr," use -b option for bactig stores  or -f (default) for read stores\n");
	exit (EXIT_FAILURE);
      }

    gatekeeperStorePath = argv[optind++];

    /* End of command line parsing */
  }
   


  /**************** Open or Create Files *********************/
  fprintf(stderr,"* GatekeeperStorePath is %s\n",
	  gatekeeperStorePath);

  InitGateKeeperStore(&gkpStore, gatekeeperStorePath);
  OpenReadOnlyGateKeeperStore(&gkpStore);
  {
    PHashValue_AS value;
     
     
    while(EOF != scanf(F_UID, &uid)){
      if(HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore.hashTable, 
						   UID_NAMESPACE_AS,
						   uid, 
						   type, 
						   TRUE,
						   stderr,
						   &value)){
	fprintf(stderr,"* Couldn't find frag with uid " F_UID " ...\n", uid);
      }else{
	GateKeeperFragmentRecord gkf;
	GateKeeperLinkRecordIterator iterator;
	GateKeeperLinkRecord link;
	CDS_IID_t fragIID;
	fragIID = value.IID;
     
	if(type == AS_IID_FRG){

	getGateKeeperFragmentStore(gkpStore.frgStore, fragIID, &gkf);
       fprintf(stderr,"* uid:" F_UID " Fragment " F_IID ": UID:" F_UID " type%c refs: %d links:%d(" F_IID ") lID:" F_IID " sID:" F_IID " bID:" F_IID " batch(%u,%u)  \n",
	       uid,
	       fragIID, 
	       gkf.readUID, 
	       gkf.type,
	       value.refCount, gkf.numLinks, gkf.linkHead,
	       gkf.localeID, gkf.seqID, gkf.bactigID, gkf.birthBatch, gkf.deathBatch);
       fflush(stderr);
	if(gkf.numLinks > 0){
	  CreateGateKeeperLinkRecordIterator(gkpStore.lnkStore, gkf.linkHead,fragIID, &iterator);
	  while(NextGateKeeperLinkRecordIterator(&iterator, &link)){
	    CDS_UID_t mateUID;
	    int reversed = !(link.frag1 == fragIID);
	    getGateKeeperFragmentStore(gkpStore.frgStore, (reversed?link.frag1:link.frag2), &gkf);
	    mateUID =  gkf.readUID;
	    fprintf(stderr,"\tLink (" F_IID "," F_IID ") (" F_UID "," F_UID ") dist: " F_IID " type %d\n",
		    link.frag1, link.frag2, 
		    (reversed?mateUID:uid), 
		    (reversed?uid:mateUID), 
		    link.distance, link.type);
	  }
	}
	}else{
	  GateKeeperBactigRecord gkb;
	getGateKeeperBactigStore(gkpStore.btgStore, fragIID, &gkb);
       fprintf(stderr,"* uid:" F_UID " Bactig " F_IID ": UID:" F_UID " refs: %d lID:" F_IID " sID:" F_IID " bID:" F_IID "\n",
	       uid,
	       fragIID, 
	       gkb.UID, 
	       value.refCount,
	       gkb.bacID, gkb.seqID, fragIID);
       fflush(stderr);




	}
      }
    }
  }
  return 0;
}
