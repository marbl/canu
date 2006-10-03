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


static char CM_ID[] = "$Id: asmInfoForFrg.c,v 1.6 2006-10-03 21:49:53 brianwalenz Exp $";


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
#include "Instrument_CGW.h"
int USE_SDB;
int USE_SDB_PART;


int main( int argc, char *argv[])
{
  int32 restartFromCheckpoint = NULLINDEX;
  Global_CGW *data;
  char *inputPath;
  char *prefix;
  MesgReader reader;
  MesgWriter writer;
  FILE *myerr = stderr; 
  FILE *myout = stdout; 
  char *outputPath = NULL;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int setSingleSid = FALSE, singleSid;
  int ckptNum = NULLINDEX;
  int frgIID;
  int useIndexRatherThanIID=0; 

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "c:f:g:n:i")) != EOF)){
      switch(ch) {
        case 'i':
          useIndexRatherThanIID=1;
          break;
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
	fprintf (stderr, "USAGE:  %s -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum>\n",argv[0]);
	exit (EXIT_FAILURE);
      }
  }

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);

  while(scanf("%d\n",&frgIID)==1){
    InfoByIID * info;
    CIFragT * frag;
    if(useIndexRatherThanIID){
      frgIID =  GetCIFragT(ScaffoldGraph->CIFrags, frgIID)->iid;
    }       
    // Don't use the convenience function, since we want to print the index
    info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, frgIID);
    if(info==NULL){
      printf("WARNING: couldn't identify fragment %d in assembly\n",frgIID);
      continue;
    }
    frag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
    assert(frag!=NULL);
    assert(GetGraphNode(ScaffoldGraph->ContigGraph,frag->contigID)!=NULL);
    printf("scaffold %d ",GetGraphNode(ScaffoldGraph->ContigGraph,frag->contigID)->scaffoldID);
    PrintFragment(frag, info->fragIndex, stdout);
  }

  exit(0);
}



