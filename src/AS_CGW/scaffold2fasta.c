
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
static char CM_ID[] = "$Id: scaffold2fasta.c,v 1.1.1.1 2004-04-14 13:51:00 catmandew Exp $";


/*********************************************************************
 * Module:  scaffold2fasta
 * Description:
 *    Dumps a scaffold in fasta format
 * 
 *    Reference: 
 *
 *    Command Line Interface:
 *        $ scaffold2fasta -f <frgStore> -g <gkpStore> -n <ckptNum> -c <basename> -s <scaffNum>
 *
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
#include "GapWalkerREZ.h"
FILE *  File_Open (const char * Filename, const char * Mode, int exitOnFailure);

int main(int argc, char *argv[]){
  Global_CGW *data;
  char *outputPath = NULL;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int setStartScaff = FALSE;
  int ckptNum = NULLINDEX;
  CDS_CID_t startScaff = NULLINDEX;
  NodeCGW_T* scaff;
  size_t i;
  char fastaOutputFilename[1024];
  FILE *fastaOutputFile;

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->stderro = stderr;
  data->stderrfp = fopen("loadcgw.stderr","w");
  
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,	"f:g:n:c:s:")) != EOF)){
      switch(ch) 
      {
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case 's':
        {
          startScaff = atoi(argv[optind - 1]);
          setStartScaff = 1;
        }
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
    if((setPrefixName == FALSE) || (setFragStore == 0) ||
       (setGatekeeperStore == 0) || (setStartScaff == 0))
    {
      fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d outputPath = %s\n",
              argc, optind, setFragStore,setGatekeeperStore, outputPath);
      fprintf (stderr, 
               "USAGE:  scaffold2fasts -f <FragStore> -g <GkpStore> -c <CkptFileName> -n <CkpPtNum> -s <scaffNum>\n");
      exit (EXIT_FAILURE);
    }
  }
  
  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix,
                                                  ckptNum, FALSE);
  
  scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, startScaff);
  
  // make sure the scaffold is there
  assert(scaff != NULL);
  
  sprintf( fastaOutputFilename, "%s.scaff_" F_CID ".fa",
           data->File_Name_Prefix, scaff->id );
  fastaOutputFile = File_Open( fastaOutputFilename, "w", TRUE );
  
  // not interested in dead, not real, or singleton scaffolds
  if ((isDeadCIScaffoldT(scaff)) || (scaff->type != REAL_SCAFFOLD))
  {
    fprintf( stderr, "scaffold " F_CID " is dead (%d) or not real (%d)\n", 
             scaff->id, isDeadCIScaffoldT(scaff),
             (scaff->type != REAL_SCAFFOLD));
    exit(1);
  }
  else
  {
    CIScaffoldTIterator CIsTemp;
    ChunkInstanceT *contig = NULL, *nextContig = NULL;
    static VA_TYPE(char) *contigConsensus = NULL;
    static VA_TYPE(char) *contigQuality = NULL;
    char *sequence;
    
    if ( contigConsensus == NULL )
    {
      contigConsensus = CreateVA_char(1024);
      contigQuality = CreateVA_char(1024);
    }
    
    fprintf( fastaOutputFile,
             ">scaffold " F_CID ", scaffold length %.f\n", 
             scaff->id, scaff->bpLength.mean);
    
    InitCIScaffoldTIterator(ScaffoldGraph, scaff, TRUE, FALSE, &CIsTemp);
    while (NextCIScaffoldTIterator(&CIsTemp))
    {
      char temp;
      LengthT gapSize;
      
      contig = GetGraphNode( ScaffoldGraph->ContigGraph, CIsTemp.curr);      
      assert( contig != NULL);
      if ( CIsTemp.next != NULLINDEX )
      {
        nextContig = GetGraphNode( ScaffoldGraph->ContigGraph, CIsTemp.next);
        assert( nextContig != NULL);
      }
      
      GetConsensus( ScaffoldGraph->ContigGraph, contig->id, contigConsensus, contigQuality);
      sequence = Getchar( contigConsensus, 0);
      
      for ( i = 0; i < strlen( sequence ); i += 60)
      {
        if ( i + 60 < strlen( sequence ))
        {
          temp = sequence[i + 60];
          sequence[i + 60] = '\0';
          fprintf( fastaOutputFile, "%s\n", &sequence[i]);
          sequence[i + 60] = temp;
        }
        else
        {
          fprintf( fastaOutputFile, "%s", &sequence[i]);		  
        }
      }
      
      if ( CIsTemp.next != NULLINDEX )
      {
        CDS_COORD_t pos = 0;
        gapSize = FindGapLength( contig, nextContig, FALSE);
        if ( (int) gapSize.mean < 20 ) gapSize.mean = 20;
        for ( i = 0; i < (int) gapSize.mean; i++)
        {
          if ( (pos % 60) == 0)
            fprintf( fastaOutputFile, "\n");
          fprintf( fastaOutputFile, "N");
          pos++;
        }
        fprintf( fastaOutputFile, "\n");
      }
    }
  }
  return 0;
}

	  

