
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
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
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
//extern int do_compute_missing_overlaps;

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

void print_olap(Long_Olap_Data_t olap){
  printf ("    %8d %8d %c %5d %5d %4.1f %4.1f\n",
          olap . a_iid,
          olap . b_iid,
          olap . flipped ? 'I' : 'N',
          olap . a_hang, olap . b_hang,
          olap . orig_erate / 10.0, olap . corr_erate / 10.0);
}

void usage(char *pgm){
  fprintf (stderr, "USAGE:  %s -g <GatekeeperStoreName> -o <OVLStoreName> -c <CkptFileName> -n <CkpPtNum>\n",
           pgm);
}

int main (int argc , char * argv[] ) {

  Global_CGW *data;
  char *prefix;
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

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,"c:g:n:1:2:o:")) != EOF)){
      switch(ch) {
        case 'c':
          strcpy( data->File_Name_Prefix, argv[optind - 1]);
          setPrefixName = TRUE;		  
          break;
        case 'g':
          strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
          setGatekeeperStore = TRUE;
          break;	  
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case 'o':
          strcpy( data->OVL_Store_Name, argv[optind - 1]);
          setOvlStore = TRUE;
          break;	  
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if((setPrefixName == FALSE) || (setGatekeeperStore == 0) || ( setOvlStore == 0)){
      fprintf(stderr,"* argc = %d optind = %d setGatekeeperStore = %d\n",
              argc, optind, setGatekeeperStore);

      usage(argv[0]);
      exit (-1);
    }
  }

  ScaffoldGraph = 
    LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);
  setup_ovlStore();

  {
    int frgiid;
    while ( scanf("%d",&frgiid) == 1 ){

      Long_Olap_Data_t  olap;
      static OVL_Stream_t  * my_stream = NULL;

      InfoByIID * info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, frgiid);
      if(info==NULL){
        fprintf(stderr,"Could not get info for fragment %d\n",frgiid);
        continue;
      }


      if(my_stream == NULL){
	my_stream = New_OVL_Stream ();
      } 

      Init_OVL_Stream (my_stream, frgiid, frgiid, ScaffoldGraph->frgOvlStore);

      while  (Next_From_OVL_Stream (& olap, my_stream)){
	char *infoString;
	int ovlIID = olap . b_iid;
	InfoByIID * info = GetInfoByIID(ScaffoldGraph->iidToFragIndex, ovlIID);
	CIFragT * ovlFrag = GetCIFragT(ScaffoldGraph->CIFrags, info->fragIndex);
	int utgID = ovlFrag->cid;
	ChunkInstanceT *unitig = GetGraphNode(ScaffoldGraph->CIGraph,utgID);
	int numInst = unitig->info.CI.numInstances;
	int ctgID,scfID;

	static char *locs=NULL;
	static int lenloc=0;
	int lenUsed = 0;
	if(locs==NULL){
	  lenloc = 1000;
	  locs = (char *) malloc(lenloc*sizeof(char));
	  assert(locs!=NULL);
	}
	locs[0]='\0';

	if(unitig->info.CI.numInstances>0){
	  if(numInst<=2){
	    safelyAppendInstInfo(&locs,unitig->info.CI.instances.in_line.instance1,&lenloc,&lenUsed);
	    if(numInst==2){
	      safelyAppendInstInfo(&locs,unitig->info.CI.instances.in_line.instance2,&lenloc,&lenUsed);
	    }
	  } else {
	    int i,n;
	    int32 *inst_list;
	    n=unitig->info.CI.numInstances;
	    assert(n == GetNumint32s(unitig->info.CI.instances.va));
	    inst_list = Getint32(unitig->info.CI.instances.va,0);
	    for(i=0;i<n;i++){
	      safelyAppendInstInfo(&locs,inst_list[i],&lenloc,&lenUsed);
	    }
	  }
	} else {
	  ctgID = unitig->info.CI.contigID;
	  scfID = unitig->scaffoldID;
	  sprintf(locs," utg:%d ctg:%d scf:%d",unitig->id,ctgID,scfID);
	}

	if  (olap . a_hang < 0){
	  printf("\tAEnd\t%d\t%s ",ovlIID,locs);
	  print_olap(olap);
	}
	if  (olap . b_hang > 0){
	  printf("\tBEnd\t%d\t%s ",ovlIID,locs);
	  print_olap(olap);
	}

      }
    }
  }
  finished_with_ovlStore ();

}
