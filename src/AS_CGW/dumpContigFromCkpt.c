
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


static char CM_ID[] = "$Id: dumpContigFromCkpt.c,v 1.3 2005-03-22 19:48:37 jason_miller Exp $";


/*********************************************************************/

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

#include "cds.h"
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

#define NUM_STDDEV_CUTOFF 3.0

#define COMPARE_ARGS char *aseq, char *bseq, int beg, int end, int opposite, \
                    double erate, double thresh, int minlen, \
				 CompareOptions what



void usage(char *pgmname){
  fprintf (stderr, 
	   "USAGE:  %s -c <CkptFileName> -n <CkpPtNum> -f <frgStore> -g <gkpStore> [-r] -i ID\n"
	   "\t-r option causes reverse complement\n",
	   pgmname);
  exit (EXIT_FAILURE);
}

void SequenceComplement(char *sequence, char *quality);

void printContigViaCheckpoint( int contigID, int reverse){

  VA_TYPE(char) *ContigConsensus = NULL;
  VA_TYPE(char) *ContigQuality = NULL;
  char *Sequence;
  int len;

  ContigConsensus = CreateVA_char(256*1024);
  ContigQuality = CreateVA_char(256*1024);

  GetConsensus( ScaffoldGraph->ContigGraph, contigID, ContigConsensus, ContigQuality);
  Sequence = Getchar( ContigConsensus, 0);

  if(reverse){
    SequenceComplement( Sequence, NULL);
  }

  len = strlen(Sequence);
  fprintf(stdout,">contig %d len %d\n",contigID,len);
  {int i;
  for(i=0;i<len;i+=60){
    fprintf(stdout,"%.60s\n",Sequence+i);
  } }

  return;
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
  int i, index;
  int sid, startingGap, setStartingGap = FALSE;
  int numGaps = 0, gapNumber = 0;
  int numGapsVarTooSmall = 0;
  time_t t1;
  int overlapFound = 0, noOverlapFound = 0;
  int dumpSeqs=0;
  int reverse=0;
  int contigID=NULLINDEX;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->stderro = stderr;
  data->stderrfp = fopen("findMissedOverlaps.stderr","w");
  data->timefp = stderr;
  data->logfp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "c:f:g:n:s:ri:")) != EOF)){
      switch(ch) {
      case 'f':
	strcpy( data->Frag_Store_Name, argv[optind - 1]);
	setFragStore = TRUE;
	break;
      case 'g':
	strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
	setGatekeeperStore = TRUE;
	break;	  
      case 'c':
	strcpy( data->File_Name_Prefix, argv[optind - 1]);
	setPrefixName = TRUE;		  
	break;
      case 'n':
	ckptNum = atoi(argv[optind - 1]);
	break;
      case 'i':
	contigID = atoi(argv[optind - 1]);
	break;
      case 'r': 
	reverse=1;
      case '?':
      case 'h':
      default :
	errflg++;
	usage(argv[0]);
      }
    }
  }

  if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0) || (contigID == NULLINDEX))
    {
      usage(argv[0]);
    }
  

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, FALSE);

  printContigViaCheckpoint( contigID, reverse);
  return 0;
}
