
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
static char CM_ID[] = "$Id: AS_CGW_LoadCheckpoint.c,v 1.4 2005-03-22 19:48:33 jason_miller Exp $";


/*********************************************************************
 * Module:  AS_CGW_LoadCheckpoint
 * Description:
 *    For use with debugger to query values in a checkpoint
 * 
 *    Reference: 
 *
 *    Command Line Interface:
 *        $ loadcgw checkpointPath
 *
 *       CGBInputFiles: The file with new IUM,OUM, etc records to process. 
 *
 *       Checkpoint File: File named <outputName>.ckp.n
 *
 * 
 *********************************************************************/
//#define DEBUG 1
//#define DEBUG_BUCIS 1
//#define DEBUG_MERGE_SCAF 1

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
#include "Stats_CGW.h"
#include "Instrument_CGW.h"

FILE *  File_Open (const char * Filename, const char * Mode, int exitOnFailure);



int main(int argc, char *argv[]){
  Global_CGW *data;
  char *outputPath = NULL;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int ckptNum = NULLINDEX;
  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->stderro = stderr;
  data->stderrfp = fopen("loadcgw.stderr","w");
  data->logfp = fopen("loadcgw.log","w");

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "f:g:n:c:")) != EOF)){
      switch(ch) {
      case 'n':
	ckptNum = atoi(argv[optind - 1]);
	break;
      case 'c':
	{
	  strcpy( data->File_Name_Prefix, argv[optind - 1]);
	  setPrefixName = 1;

	}
	break;
      case 'f':
	{
	  strcpy( data->Frag_Store_Name, argv[optind - 1]);
	  setFragStore = 1;
	}
	break;
      case 'g':
	{
	  strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
	  setGatekeeperStore = 1;
	}
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
	fprintf (stderr, "USAGE:  loadcgw -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum>\n");
	exit (EXIT_FAILURE);
      }
  }

    ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix,ckptNum, FALSE);    
    fprintf(stderr,"* Hit Interrupt in Debugger to Proceed!!!! *\n");
    fflush(stderr);
    while(1){
      // wait for interrupt in debugger
    }
  }
