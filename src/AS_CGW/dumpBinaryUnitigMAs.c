
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
static char CM_ID[] = "$Id: dumpBinaryUnitigMAs.c,v 1.9 2007-02-08 02:46:00 brianwalenz Exp $";


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
#include "AS_PER_ReadStruct.h"
#include "PublicAPI_CNS.h"
#include "MultiAlignment_CNS.h"
#include "fixZLFContigs.h"

FILE *  File_Open (const char * Filename, const char * Mode, int exitOnFailure);

ReadStructp myRead;

void FindPattern(void)
{
  GraphCGW_T *graph = ScaffoldGraph->ContigGraph;
  GraphNodeIterator     nodes;
  ContigT		*ctg;
  int numMatches = 0;
  
  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  while((ctg = NextGraphNodeIterator(&nodes)) != NULL)
    {
      NodeCGW_T * utg;
    

      if(ctg->flags.bits.isDead) continue;
      if(ctg->flags.bits.isChaff) continue;
      if(ctg->info.Contig.numCI > 1) continue;
      if(ctg->scaffoldID != NULLINDEX) continue;

      utg = GetGraphNode(ScaffoldGraph->CIGraph, ctg->info.Contig.AEndCI);
      if(utg->type != UNRESOLVEDCHUNK_CGW ||
         utg->scaffoldID == NULLINDEX ||
         utg->info.CI.numInstances > 0) continue;

      /*
        fprintf(stderr, "utg " F_CID ", scaffold " F_CID ", type %s, rock: %c, potential rock: %c, stone: %c, potential stone: %c, numInstances: %d\n",
        utg->id, utg->scaffoldID,
        (utg->type == UNRESOLVEDCHUNK_CGW) ? "unresolved" :
        ((utg->type == DISCRIMINATORUNIQUECHUNK_CGW) ? "Dunique" :
        ((utg->type == UNIQUECHUNK_CGW) ? "Unique" :
        ((utg->type == RESOLVEDREPEATCHUNK_CGW) ? "resolved" : "null"))),
        utg->flags.bits.isRock ? 'Y' : 'N',
        utg->flags.bits.isPotentialRock ? 'Y' : 'N',
        utg->flags.bits.isStone ? 'Y' : 'N',
        utg->flags.bits.isPotentialStone ? 'Y' : 'N',
        utg->info.CI.numInstances);
      */

      fprintf(stdout, F_CID "\n", utg->id);
      numMatches++;
    }

  fprintf(stderr, "%d matches.\n", numMatches);
}

void PrintFasta(FILE * fp, char * header, char * seq,
                uint32 start, uint32 end)
{
  uint32 i;
  fprintf( fp, "%s", header );

  for( i = start; i < end; i++ )
    {
      if( i != start && (i-start) % 60 == 0 )
        fprintf( fp, "\n" );
      fprintf( fp, "%c", seq[i] );
    }
  fprintf( fp, "\n" );
}

void PrintRSFasta( FILE * fp, CDS_CID_t iid, ReadStructp rs)
{
  uint32 start, end;
  char seq[2048], qua[2048];
  char header[1024];

  getClearRegion_ReadStruct(rs, &start, &end, READSTRUCT_CGW);
  getSequence_ReadStruct(rs, seq, qua, 2048);
  sprintf(header, ">" F_CID "|Fragment|len=" F_SIZE_T "|clr=(%u,%u)\n",
          iid, strlen(seq), start, end);
  PrintFasta(fp, header, seq, start, end);
}

void PrintMAFasta(FILE * fp, MultiAlignT * ma)
{
  char header[1024];
  size_t length = strlen(Getchar(ma->consensus, 0));

  sprintf(header, ">" F_CID "|%s|len=" F_SIZE_T "\n",
          ma->id,
          (ma->u_list == NULL ||
           GetNumIntUnitigPoss(ma->u_list) == 0) ? "Unitig" : "Contig",
          length);
  PrintFasta(fp, header, Getchar(ma->consensus, 0), 0, length);
}

void DumpUnitigMultiAlignInfo ( CDS_CID_t unitigID )
{
  MultiAlignT *uma = CreateEmptyMultiAlignT();
  int i;
  FILE * fp;
  char fname[1024];

  uma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB,
                                       unitigID, TRUE);

  {
    sprintf(fname, "U" F_CID ".fa", unitigID);
    fp = fopen(fname, "w");
    PrintMAFasta(fp, uma);
    fclose(fp);
  }
  
  fprintf( stderr, "\tunitig %8" F_CIDP ", strlen( consensus ): %9" F_SIZE_TP "\n",
           unitigID, strlen( Getchar( uma->consensus, 0 )));

  for ( i = 0; i < GetNumIntMultiPoss( uma->f_list ); i++)
    {
      IntMultiPos *pos = GetIntMultiPos( uma->f_list, i);
    
      getFragStore(ScaffoldGraph->fragStore, pos->ident, FRAG_S_ALL, myRead);
      sprintf(fname,"U" F_CID "_%d_F" F_IID ".fa",
              unitigID, i, pos->ident);
      fp = fopen(fname, "w");
      PrintRSFasta(fp, pos->ident, myRead);
      fclose(fp);
    
      fprintf( stderr, "\t\tfragment %8" F_IIDP ", bgn: %10" F_COORDP ", "
               "end: %10" F_COORDP ", length: %10" F_COORDP ", source: %d\n", 
               pos->ident,
               pos->position.bgn, pos->position.end,
               abs(pos->position.bgn - pos->position.end),
               pos->sourceInt);
    }
}


void DumpContigMultiAlignInfo ( CDS_CID_t contigID )
{
  MultiAlignT *cma = CreateEmptyMultiAlignT();
  int i;
  FILE * fp;
  char fname[1024];

  fprintf( stderr,
           "------------------------------------------------------------\n");
  
  cma = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB,
                                       contigID, FALSE);

  {
    sprintf(fname, "C" F_CID ".fa", contigID);
    fp = fopen(fname, "w");
    PrintMAFasta(fp, cma);
    fclose(fp);
  }
  
  fprintf( stderr, "contig %8" F_CIDP ", strlen( consensus ): %9" F_SIZE_TP "\n",
           contigID, strlen( Getchar( cma->consensus, 0 )));

  for ( i = 0; i < GetNumIntMultiPoss( cma->f_list ); i++)
    {
      IntMultiPos *pos = GetIntMultiPos(cma->f_list,i);
    
      getFragStore(ScaffoldGraph->fragStore, pos->ident, FRAG_S_ALL, myRead);
      sprintf(fname,"C" F_CID "_%d_F" F_IID ".fa",
              contigID, i, pos->ident);
      fp = fopen(fname, "w");
      PrintRSFasta(fp, pos->ident, myRead);
      fclose(fp);
    
      fprintf( stderr, "\t\tfragment %8" F_IIDP ", "
               "bgn: %10" F_COORDP ", end: %10" F_COORDP ", length: %10" F_COORDP "\n", 
               pos->ident,
               pos->position.bgn, pos->position.end,
               abs(pos->position.bgn - pos->position.end));
    }
  for ( i = 0; i < GetNumIntUnitigPoss( cma->u_list ); i++)
    {
      IntUnitigPos *pos = GetIntUnitigPos( cma->u_list, i);
      NodeCGW_T *unitig = GetGraphNode( ScaffoldGraph->CIGraph, pos->ident);
    
      fprintf( stderr, "\tunitig %8" F_CIDP ", "
               "bgn: %10" F_COORDP ", end: %10" F_COORDP ", length: %10" F_COORDP "\n", 
               unitig->id,
               pos->position.bgn, pos->position.end,
               abs(pos->position.bgn - pos->position.end));
    
      DumpUnitigMultiAlignInfo( unitig->id );	
    }
  fprintf( stderr, "\n");
}

int main(int argc, char *argv[]){
  Global_CGW *data;
  char *outputPath = NULL;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int ckptNum = NULLINDEX;
  VA_TYPE(CDS_CID_t) * clist = CreateVA_CDS_CID_t(100);
  char * unitigIDsFile = NULL;
  
  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "i:f:g:n:c:u:")) != EOF)){
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
        case 'i':
          {
            CDS_CID_t iid = atoi(argv[optind-1]);
            AppendVA_CDS_CID_t(clist, &iid);
          }
          break;
        case 'u':
          unitigIDsFile = argv[optind-1];
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }
    if((setPrefixName == FALSE) || (setFragStore == 0) || (setGatekeeperStore == 0) || unitigIDsFile == NULL)
      {
	fprintf(stderr,"* argc = %d optind = %d setFragStore = %d setGatekeeperStore = %d outputPath = %s\n",
		argc, optind, setFragStore,setGatekeeperStore, outputPath);
	fprintf (stderr, "USAGE:  loadcgw -f <FragStoreName> -g <GatekeeperStoreName> -c <CkptFileName> -n <CkpPtNum>\n");
	exit (EXIT_FAILURE);
      }
  }

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix,ckptNum, FALSE);

  // initialize other important variables
  GlobalData->aligner=Local_Overlap_AS_forCNS;
  
  // initialize globals required by consensus
  cnslog = stderr;
  USE_SDB=1;
  RALPH_INIT = InitializeAlphTable();
  sequenceDB = ScaffoldGraph->sequenceDB;
  global_fragStore = ScaffoldGraph->fragStore;


  {
    /*
      Open file listing unitig IIDs
      Loop over them
      get the MA
      convert to protoIO
      write out
    */
    FILE * fp = fopen(unitigIDsFile, "r");
    FILE * fpout = fopen("outFile.mas", "w");
    VA_TYPE(CDS_CID_t) * iids = CreateVA_CDS_CID_t(100);
    char line[1024];
    int i;

    while(fgets(line, 1024, fp))
      {
        CDS_CID_t id = atoi(line);
        AppendVA_CDS_CID_t(iids, &id);
      }

    for(i = 0; i < GetNumVA_CDS_CID_t(iids); i++)
      {
        CDS_CID_t * idp = GetVA_CDS_CID_t(iids, i);
        MultiAlignT * ma =
          LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, *idp, TRUE);

        WriteMAToFile(ma, fpout);
      }
    fclose(fp);
    fclose(fpout);
  }
#ifdef NEVER
  /*
    Experiment to redo unitig 298448 in contig 517900 in scaffold 4000
    in dros checkpoint 25 in /prod/IR03/dros5_20011107/workbox/mike/dewim/
    Get multialignT of unitig
    Set up ZLFScaffold data structure
    Call FixZLFContigs()
  */
  {
    VA_TYPE(ZLFScaffold) * zlfScaffolds;
    ZLFScaffold zlfScaffold;
    ZLFContig zlfContig;
    MultiAlignT * ma;
    char * consensus;
    int i;

    ma =
      CloneMultiAlignT(LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
                                                     298448, TRUE));

    ma->id = 298448;
    // muck up part of the consensus that doesn't overlap other unitigs
    // to check that fix 'took'
    consensus = GetVA_char(ma->consensus, 0);
    for(i = 450; i < 1000; i++)
      {
        consensus[i] = 'C';
      }
          
    zlfContig.id = 517900;
    zlfContig.zlfUMAs = CreateVA_MultiAlignT(1);
    AppendVA_MultiAlignT(zlfContig.zlfUMAs, ma);

    zlfScaffold.id = 4000;
    zlfScaffold.zlfContigs = CreateVA_ZLFContig(1);
    AppendVA_ZLFContig(zlfScaffold.zlfContigs, &zlfContig);
    
    zlfScaffolds = CreateVA_ZLFScaffold(1);
    AppendVA_ZLFScaffold(zlfScaffolds, &zlfScaffold);

    FixZLFContigs(zlfScaffolds, TRUE, TRUE);
    
    // now dump a checkpoint, reload it & check again
    CheckpointScaffoldGraph(ScaffoldGraph, 1);
    DestroyScaffoldGraph(ScaffoldGraph);
    ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix,
                                                    ckptNum + 1, FALSE);
    
    ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
                                       298448, TRUE);

    consensus = GetVA_char(ma->consensus, 0);
    for(i = 450; i < 1000; i++)
      {
        if(consensus[i] != 'C')
          {
            fprintf(stderr, "Fix didn't 'take'!\n");
            break;
          }
      }
    if(i == 1000)
      fprintf(stderr, "Fix 'took'!\n");
  }
#endif

  /*
    fprintf(stderr,"* Hit Interrupt in Debugger to Proceed!!!! *\n");
    fflush(stderr);
    while(1){
    // wait for interrupt in debugger
    }
  */

  /* 
     myRead =  new_ReadStruct();

     {
     int i;
     for(i = 0; i < GetNumVA_CDS_CID_t(clist); i++)
     DumpContigMultiAlignInfo(*(GetVA_CDS_CID_t(clist, i)));
     }
  */
  /*
    {
    int i;
    for(i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++)
    {
    DistT *dptr;
    dptr = GetDistT(ScaffoldGraph->Dists,i);
        
    fprintf(stderr, "\tiid:%d, mean:%f, stddev:%f\n", i,
    dptr->mean, dptr->stddev);
    }
    }
  */
    
  /*
    {
    FILE * fp;
    char line[1024];
      
    fp = fopen("/prod/IR04/RAT_ASSEMBLY/workbox/dewim/reestimatedDists.txt", "r");
    assert(fp != NULL);
    while(fgets(line, 1023, fp) != NULL)
    {
    DistT *dptr;
    CDS_CID_t iid;
    float32 mean;
    float32 stddev;

    fprintf(stderr, "line is : %s", line);
        
    sscanf(line, "%d %f %f", &iid, &mean, &stddev);
    if(iid <= 0 || iid >= GetNumDistTs(ScaffoldGraph->Dists))
    continue;
        
    fprintf(stderr, "\tiid:%d, mean:%f, stddev:%f\n", iid, mean, stddev);
        
    dptr = GetDistT(ScaffoldGraph->Dists,iid);

    fprintf(stderr, "got dist pointer\n");

    fprintf(stdout, "%d: was (%f,%f) is (%f,%f)\n",
    iid, dptr->mean, dptr->stddev, mean, stddev);
    fflush(stdout);
        
    dptr->mean = mean;
    dptr->stddev = stddev;
    }
    fclose(fp);
    }
    fprintf(stderr,"Checkpoint %d written after distance re-estimation\n",
    ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph, 2);
  */
  return 0;
}

