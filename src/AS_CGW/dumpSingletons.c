
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

static char CM_ID[] = "$Id: dumpSingletons.c,v 1.8 2006-09-21 21:34:01 brianwalenz Exp $";


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
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"

#include "SYS_UIDcommon.h"
#include "SYS_UIDclient.h"

int USE_SDB;
int USE_SDB_PART;

#ifndef LD
#ifdef linux
#define LD "%lld"
#else
#define LD "%ld"
#endif
#endif


static void Complement(char *seq, int len)
{ static char WCinvert[256];
 static int Firstime = 1;

 if (Firstime)          /* Setup complementation array */
   { int i;

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

 /* Complement and reverse sequence */

 { register char *s, *t;
 int c;

 s = seq;
 t = seq + (len-1);
 while (s < t)
   { c = *s;
   *s++ = WCinvert[(int) *t];
   *t-- = WCinvert[c];
   }
 if (s == t)
   *s = WCinvert[(int) *s];
 }
}

int main( int argc, char *argv[])
{
  int32 restartFromCheckpoint = NULLINDEX;
  Global_CGW *data;
  char *inputPath;
  char *prefix;
  MesgReader reader;
  MesgWriter writer;
  MesgWriter errorWriter;
  FILE *myerr = stderr; 
  FILE *myout = stdout; 
  char *outputPath = NULL;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int setSingleSid = FALSE, singleSid;
  int ckptNum = NULLINDEX;
  int mateIID;
  int ifrag;
  FragStoreHandle storeHandle = 0;
  GateKeeperStore gkpStore;
  GateKeeperFragmentRecord gkpFrag,gkpMate;
  CIFragT *frag,*mate;
  uint64 uid, mateuid;
  char *seq1,*seq2,*qul1,*qul2,*toprint1,*toprint2;
  uint clr_bgn1,clr_end1;
  uint clr_bgn2,clr_end2;
  int alloclen1=1000;
  int alloclen2=1000;
  ReadStructp fsread=new_ReadStruct();
  ReadStructp fsmate=new_ReadStruct();
  int realUID=0;
  int UIDstart=1230000;
  int firstUID=1;
  CDS_UID_t       interval_UID[4];

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "c:f:g:n:U")) != EOF)){
      switch(ch) {
        case 'c':
          strcpy( data->File_Name_Prefix, argv[optind - 1]);
          setPrefixName = TRUE;		  
          break;
        case 'f':
          strcpy( data->Frag_Store_Name, argv[optind - 1]);
          setFragStore = TRUE;
          break;
        case 'g':
          strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
          setGatekeeperStore = TRUE;
          break;	  
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
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

    if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0))
      {
	fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d outputPath = %s\n",
		argc, optind, setFragStore,setGatekeeperStore, outputPath);
	fprintf (stderr, "USAGE:  %s -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum> [-U]\n",argv[0]);
	exit (EXIT_FAILURE);
      }

  }
  seq1=(char*)malloc(sizeof(char)*alloclen1);
  qul1=(char*)malloc(sizeof(char)*alloclen1);
  toprint1=(char*)malloc(sizeof(char)*alloclen1);
  assert(seq1!=NULL);
  assert(qul1!=NULL);
  assert(toprint1!=NULL);
  seq2=(char*)malloc(sizeof(char)*alloclen2);
  qul2=(char*)malloc(sizeof(char)*alloclen2);
  toprint2=(char*)malloc(sizeof(char)*alloclen2);
  assert(seq2!=NULL);
  assert(qul2!=NULL);
  assert(toprint2!=NULL);

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);

  for (ifrag = 0; ifrag < GetNumVA_CIFragT( ScaffoldGraph->CIFrags ); ifrag++){

    frag = GetCIFragT( ScaffoldGraph->CIFrags, ifrag);
    assert(frag->cid!=NULLINDEX);

    if(GetGraphNode(ScaffoldGraph->CIGraph,frag->cid)->flags.bits.isChaff){
      InfoByIID * info;

      if(frag->numLinks>0){
	assert(frag->numLinks==1);
	mate = GetCIFragT(ScaffoldGraph->CIFrags,frag->mateOf);

        //  Hmmm, why don't we have a mate?!  Probably a bug somewhere
        //  (in the input, perhaps??)  Perhaps this is from OBT
        //  deleting fragments, but not deleting the link.
        //
        if (!mate) {
          frag->numLinks = 0;
        } else {
          if(mate->flags.bits.isChaff){
            if(frag->iid>mate->iid) {
              //	    printf("%d is chaff ",frag->iid);
              //	    printf(" would not print (should have been taken care of already)\n");
              continue;
            }
          }
	}

      }
	    
      if(getFragStore(ScaffoldGraph->fragStore,frag->iid,FRAG_S_ALL,fsread)!=0){
	fprintf(stderr,"Couldn't get fragment from frgStore for iid %d\n",frag->iid);
	assert(0);
      } else {
	int rv1;
	rv1 = getGateKeeperFragmentStore(ScaffoldGraph->gkpStore.frgStore,frag->iid,&gkpFrag);
	assert(rv1==0);
	getClearRegion_ReadStruct(fsread, &clr_bgn1,&clr_end1, READSTRUCT_LATEST);
	while(getSequence_ReadStruct(fsread,seq1,qul1,alloclen1)!=0){
	  alloclen1*=2;
	  seq1=(char*)realloc(seq1,alloclen1*sizeof(char));
	  qul1=(char*)realloc(qul1,alloclen1*sizeof(char));
	  toprint1=(char*)realloc(toprint1,alloclen1*sizeof(char));
	}
	strcpy(toprint1,seq1+clr_bgn1);
	toprint1[clr_end1-clr_bgn1]='\0';
      }


      if(frag->numLinks==0 || ! mate->flags.bits.isChaff){
	//	printf("%d should print by itself\n",frag->iid);


	printf(">" F_S64 " /type=singleton\n%s\n",
	       gkpFrag.readUID,toprint1);
      } else {
	int rv2;

	//	printf("%d and %d should print as mini-scaffold\n",frag->iid,mate->iid);

	rv2 = getGateKeeperFragmentStore(ScaffoldGraph->gkpStore.frgStore,mate->iid,&gkpMate);
	assert(rv2==0);
	
	if(getFragStore(ScaffoldGraph->fragStore,mate->iid,FRAG_S_ALL,fsmate)!=0){
	  fprintf(stderr,"Couldn't get fragment from frgStore for iid %d\n",mate->iid);
	  assert(0);
	} else {
	  getClearRegion_ReadStruct(fsmate, &clr_bgn2,&clr_end2, READSTRUCT_LATEST);
	  while(getSequence_ReadStruct(fsmate,seq2,qul2,alloclen2)!=0){
	    alloclen2*=2;
	    seq2=(char*)realloc(seq2,alloclen2*sizeof(char));
	    qul2=(char*)realloc(qul2,alloclen2*sizeof(char));
	    toprint2=(char*)realloc(toprint2,alloclen2*sizeof(char));
	  }
	  strcpy(toprint2,seq2+clr_bgn2);
	  toprint2[clr_end2-clr_bgn2]='\0';
	  //	  printf(" before rc, 2nd frg is:\n%s\n",toprint2);
	  Complement(toprint2,strlen(toprint2));
	}
	{
	  uint64 blockSize = 300;
	  CDS_UID_t uid;
	  int32  uidStatus;
	  CDS_UID_t interval_UID[4];
	  if(firstUID){
	    firstUID=0;
	    set_start_uid(UIDstart); /* used if readUID == FALSE */
	    get_uids(blockSize,interval_UID,realUID);
	  }

	  uidStatus = get_next_uid(&uid,realUID);
	  if( uidStatus != UID_CODE_OK )
	    {
	      get_uids(blockSize,interval_UID,realUID);
	      uidStatus = get_next_uid(&uid,realUID);
	    }	  
	  if( UID_CODE_OK != uidStatus )
	    { 
	      fprintf(stderr, "Could not get UID \n");
              assert(0);
	    }

	  // make sure the following chain of Ns is divisible by three; the exact
	  // length is arbitrary but Doug Rusch points out that by making it
	  // divisible by 3, we can get lucky and maintain the phase of a protein ...
	  // which helps in the auto-annotation of environmental samples
	  printf(">" F_S64 " /type=mini_scaffold /frgs=(" F_S64 "," F_S64 ")\n"
		 "%sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN%s\n",
		 uid,
		 gkpFrag.readUID,gkpMate.readUID,
		 toprint1,toprint2);
	}

      }

    } else {
      continue; // non-chaff fragment -- do nothing
    }
  }
  exit(0);
}
