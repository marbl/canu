
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
static char CM_ID[] = "$Id: examineFragStore.c,v 1.3 2005-03-22 19:04:04 jason_miller Exp $";


/*********************************************************************
 * Module:  based on AS_CGW_LoadCheckpoint
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
#define FRAG_POSS_SIZE 200000

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
#include "AS_PER_SafeIO.h"

int printClearRanges( int fragIid );

int debug = 0;

static ReadStructp fsread = NULL;
FragStoreHandle fragStore;

int examineClearRanges( int fragIid )
{
  uint32 clr_bgn_orig, clr_bgn_latest, clr_end_orig, clr_end_latest;
  int setStatus = 0;
  
  if( fsread == NULL ) fsread = new_ReadStruct();

  if ( fragIid != -1)
  {
	getFragStore( fragStore, fragIid, FRAG_S_ALL, fsread);	

	getClearRegion_ReadStruct( fsread, &clr_bgn_orig, &clr_end_orig, READSTRUCT_ORIGINAL);
	getClearRegion_ReadStruct( fsread, &clr_bgn_latest, &clr_end_latest, READSTRUCT_LATEST);

	if (( clr_bgn_orig != clr_bgn_latest ) || ( clr_end_orig != clr_end_latest ))
	  printClearRanges( fragIid );	
  }
  return ( setStatus ); 
}

int printClearRanges( CDS_CID_t fragIid )
{
  uint32 clr_bgn, clr_end;
  int setStatus = 0;
  
  if( fsread == NULL ) fsread = new_ReadStruct();

  if ( fragIid != -1)
  {
	getFragStore( fragStore, fragIid, FRAG_S_ALL, fsread);	

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_ORIGINAL);
	fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_ORG clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_OVL);
	fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_OVL clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
	fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_CNS clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CGW);
	fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_CGW clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

	getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_LATEST);
	fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_LAT clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
			 fragIid, clr_bgn, clr_end, clr_end - clr_bgn);
  }

  return ( setStatus ); 
}

int main( int argc, char *argv[])
{
  Global_CGW *data;
  int setFragStoreName = FALSE;
  int ch, errflg=0;
  time_t t1;
  int64 i;
  char Frag_Store_Name[1024];
  
  if (0)
  {
	GlobalData  = data = CreateGlobal_CGW();
	data->stderrc = stderr;
	data->stderro = stderr;
	data->stderrfp = stderr;
	data->timefp = stderr;
	data->logfp = stderr;
  }
  
  /* Parse the argument list using "man 3 getopt". */ 
  optarg = NULL;
  while (!errflg && ((ch = getopt(argc, argv, "c:f:g:n:")) != EOF))
  {
	switch(ch) 
	{
	  case 'f':
	  {
		strcpy( Frag_Store_Name, argv[optind - 1]);
		setFragStoreName = TRUE;
	  }
	  break;
	  case '?':
		fprintf(stderr,"Unrecognized option -%c",optopt);
	  default :
		errflg++;
	}
  }
  if (setFragStoreName == 0)
  {
	fprintf(stderr, "* argc = %d optind = %d setFragStoreName = %d\n",
			argc, optind, setFragStoreName); // , setGatekeeperStore, setCkptNum, setPrefixName);
	fprintf (stderr, "USAGE:  %s -f <FragStoreName> \n", argv[0]);
	exit (EXIT_FAILURE);
  }

  t1 = time(0);
  fprintf( stderr, "====> Starting at %s\n", ctime(&t1));

  fragStore = openFragStore( Frag_Store_Name, "r");
  if ( fragStore == NULLSTOREHANDLE )
	assert(0);

  fprintf( stderr, "examining frags from " F_S64 " to " F_S64 "\n",
		   getFirstElemFragStore(fragStore), getLastElemFragStore(fragStore));
  
  for ( i = getFirstElemFragStore(fragStore); i <= getLastElemFragStore(fragStore); i++)
  {
	examineClearRanges( i );
	if ( (i + 1) % 10000 == 0)
	  fprintf( stderr, "examined frag " F_S64 "\n", i);
  }
  
  closeFragStore( fragStore );
  return 0;
}
