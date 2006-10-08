
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
#include "Instrument_CGW.h"
#define CMDIR "CloneMiddles"

static ScaffoldInstrumenter *si;

extern int do_draw_frags_in_CelamyScaffold;
extern int do_compute_missing_overlaps;
extern do_surrogate_tracking;
extern printMateUIDs;

void dumpCloneMiddle(int sid){
  char camname[1000];
  FILE *camfile=NULL;
  CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,sid);
  if ((isDeadCIScaffoldT(scaffold)) ||
      (scaffold->type != REAL_SCAFFOLD))
    return;


  sprintf(camname,CMDIR "/scf%d_cm.cam",scaffold->id);
  camfile = fopen(camname,"w");
  assert(camfile!=NULL);
  DumpCelamyColors(camfile);
  DumpCelamyMateColors(camfile);
  if(do_draw_frags_in_CelamyScaffold)
    DumpCelamyFragColors(camfile);

  //  fprintf(stderr,"pre CelamyScaffold");

  CelamyScaffold(camfile,scaffold,0,scaffold->bpLength.mean);

  //  fprintf(stderr,"past CelamyScaffold");

  InstrumentScaffold(ScaffoldGraph,
		     scaffold,
		     si,
		     InstrumenterVerbose2,
		     GlobalData->stderrc);

  //  fprintf(stderr,"past InstrumentScaffold");

  PrintScaffoldInstrumenterMateDetails(si,camfile,PRINTCELAMY);
  PrintExternalMateDetailsAndDists(ScaffoldGraph,si->bookkeeping.wExtMates,"\t",camfile,PRINTCELAMY);
  PrintUnmatedDetails(si,camfile,PRINTCELAMY);
    
  fclose(camfile);
}

void setup_ovlStore(void){

  ScaffoldGraph->frgOvlStore = New_OVL_Store ();
  Open_OVL_Store (ScaffoldGraph->frgOvlStore, GlobalData->OVL_Store_Name);

}

void finished_with_ovlStore(void){

  Free_OVL_Store (ScaffoldGraph->frgOvlStore);

}

void usage(char *pgm){
  fprintf(stderr, "usage: %s -f <frgStore> -g <gkpStore> -o <ovlStore> -c <ckpName> -n <ckpNum> [other options]\n", pgm);
  fprintf(stderr, "  META OPTION\n");
  fprintf(stderr, "    -p <prefix>          -- attempt to guess all the required options, if your assembly\n");
  fprintf(stderr, "                            follows runCA-OBT naming conventions.\n");
  fprintf(stderr, "  REQUIRED OPTIONS\n");
  fprintf(stderr, "    -f <FragStoreName>\n");
  fprintf(stderr, "    -g <GatekeeperStoreName>\n");
  fprintf(stderr, "    -o <OVLStoreName>\n");
  fprintf(stderr, "    -c <CkptFileName>\n");
  fprintf(stderr, "    -n <CkpPtNum>\n");
  fprintf(stderr, "  OPTIONAL OPTIONS\n");
  fprintf(stderr, "    -s <single scfIID>   -- generate a single scaffold\n");
  fprintf(stderr, "    -l <min length>      -- generate only scaffolds larger than min length\n");
  fprintf(stderr, "    -U                   -- name clones with the UID of their read\n");

}

int main (int argc , char * argv[] ) {

  Global_CGW *data;
  char *prefix;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setOvlStore = FALSE;
  int setPrefixName = FALSE;
  int ckptNum = NULLINDEX;
  int i, index;
  char subset_map[1000];
  char full_map[1000];
  char setSubsetMap=0;
  char setFullMap=0;
  char ovlPath[1000];
  int setFullOvl=0;
  int specificScf = NULLINDEX;
  int minLen=0;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,"c:f:g:n:s:o:p:l:SU")) != EOF)){
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
        case 'l':
          minLen=atoi(optarg);
          assert(minLen>0);
          break;
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case 'o':
          strcpy( data->OVL_Store_Name, argv[optind - 1]);
          setOvlStore = TRUE;
          break;	  
        case 'p':
          prefix = argv[optind - 1];
          sprintf(data->Frag_Store_Name, "%s.frgStore", prefix);
          sprintf(data->Gatekeeper_Store_Name, "%s.gkpStore", prefix);
          sprintf(data->OVL_Store_Name, "%s.ovlStore", prefix);
          sprintf(data->File_Name_Prefix, "7-CGW/%s", prefix);

          setFragStore = TRUE;
          setGatekeeperStore = TRUE;
          setOvlStore = TRUE;
          setPrefixName = TRUE;

          //  Find the checkpoint number by testing what files open.  We
          //  assume checkpoints are numbered contiguously.

          {
            int  foundFirst = 0;
            int  i = 0;

            ckptNum = -1;

            for (i=0; i<256; i++) {
              char         testname[1024];
              struct stat  teststat;

              sprintf(testname, "%s.ckp.%d", data->File_Name_Prefix, i);
              fprintf(stderr, "Testing '%s'\n", testname);
              if (stat(testname, &teststat) == 0) {
                foundFirst++;
              } else {
                if (foundFirst) {
                  //  Found the checkpoint number!  It's the one before this!
                  fprintf(stderr, "Checkpoint number %d found!\n", i-1);
                  ckptNum = i - 1;
                  break;
                }
              }
            }
          }

          if (ckptNum < 1) {
            fprintf(stderr, "ERROR:  I couldn't find the checkpoints.\n");
            exit(1);
          }

          break;
        case 's':
          specificScf = atoi(argv[optind - 1]);
          break;
        case 'S':
          do_surrogate_tracking=0;
          break;
        case 'U':
          printMateUIDs=1;
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0) || ( setOvlStore == 0)){
      fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d\n",
              argc, optind, setFragStore,setGatekeeperStore);

      usage(argv[0]);
      exit (-1);
    }
  }

  ScaffoldGraph = 
    LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);
  setup_ovlStore();

  {
    static DIR *camdir=NULL;
    /* global instrumenter */
    si = CreateScaffoldInstrumenter(ScaffoldGraph, INST_OPT_ALL);
    assert(si != NULL);
  
    camdir=opendir(CMDIR);
    if(camdir==NULL){
      system("mkdir " CMDIR);
      camdir=opendir(CMDIR);
      assert(camdir!=NULL);
    }
    closedir(camdir);

    do_draw_frags_in_CelamyScaffold=1;
    do_compute_missing_overlaps=1;

    // over all scfs in graph
    if(specificScf!=NULLINDEX){
      dumpCloneMiddle(specificScf);
    } else {
      int sid;
      for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++){
	if(GetGraphNode(ScaffoldGraph->ScaffoldGraph,sid)->bpLength.mean>=minLen){
	  dumpCloneMiddle(sid);
	}
      }
    }
    DestroyScaffoldInstrumenter(si);
  }

  finished_with_ovlStore ();

}

