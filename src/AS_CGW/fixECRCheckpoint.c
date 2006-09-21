
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
static char CM_ID[] = "$Id: fixECRCheckpoint.c,v 1.7 2006-09-21 21:34:01 brianwalenz Exp $";


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
#include "UtilsREZ.h"
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

typedef struct
{
  CDS_CID_t unitigID;
  CDS_CID_t contigID;
  CDS_CID_t scaffoldID;
} ThreeIDs;
VA_DEF(ThreeIDs);

static int ThreeIDsCompare(const ThreeIDs * a, const ThreeIDs * b)
{
  if(a->scaffoldID == b->scaffoldID)
    {
      if(a->contigID == b->contigID)
        return 0;
      else
        return (int) (a->contigID - b->contigID);
    }
  else
    return (int) (a->scaffoldID - b->scaffoldID);
}


MultiAlignT * EmptyMA = NULL;
void AppendNewMAToContig(ZLFContig * zlfContig, ThreeIDs * tidp)
{
  if(EmptyMA == NULL)
    EmptyMA = safe_calloc(1, sizeof(MultiAlignT));
  
  EmptyMA->id = tidp->unitigID;
  AppendVA_MultiAlignT(zlfContig->zlfUMAs, EmptyMA);
}


void AppendNewContigToScaffold(ZLFScaffold * zlfScaffold,
                               ZLFContig * zlfContig,
                               ThreeIDs * tidp)
{
  memset(zlfContig, 0, sizeof(ZLFContig));
  zlfContig->id = tidp->contigID;
  zlfContig->zlfUMAs = CreateVA_MultiAlignT(1);

  AppendNewMAToContig(zlfContig, tidp);
  
  AppendVA_ZLFContig(zlfScaffold->zlfContigs, zlfContig);
}


void AppendNewScaffold(VA_TYPE(ZLFScaffold) * zlfScaffolds,
                       ZLFScaffold * zlfScaffold,
                       ZLFContig * zlfContig,
                       ThreeIDs * tidp)
{
  memset(zlfScaffold, 0, sizeof(ZLFScaffold));
  zlfScaffold->id = tidp->scaffoldID;
  zlfScaffold->zlfContigs = CreateVA_ZLFContig(1);
  
  AppendNewContigToScaffold(zlfScaffold, zlfContig, tidp);

  AppendVA_ZLFScaffold(zlfScaffolds, zlfScaffold);
}


int main(int argc, char *argv[]){
  Global_CGW *data;
  char *outputPath = NULL;
  int setFragStore = FALSE;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int ckptNum = NULLINDEX;
  VA_TYPE(CDS_CID_t) * clist = CreateVA_CDS_CID_t(100);
  char * unitigIDFile = NULL;
  char * maFile = NULL;

  
  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "i:f:g:n:c:u:m:")) != EOF)){
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
          unitigIDFile = argv[optind-1];
          break;
        case 'm':
          maFile = argv[optind-1];
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if((setPrefixName == FALSE) || (setFragStore == 0) ||
       (setGatekeeperStore == 0) || unitigIDFile == NULL || maFile == NULL)
      {
	fprintf(stderr,
                "* argc = %d optind = %d setFragStore = %d "
                "setGatekeeperStore = %d outputPath = %s\n",
		argc, optind, setFragStore,setGatekeeperStore, outputPath);
	fprintf (stderr,
                 "USAGE:  loadcgw\n"
                 "\t-f <FragStoreName>\n"
                 "\t-g <GatekeeperStoreName>\n"
                 "\t-c <CkptFileName>\n"
                 "\t-n <CkpPtNum>\n"
                 "\t-u <unitigIIDFile>\n"
                 "\t-m <MultiAlignTStore>\n");
	exit (EXIT_FAILURE);
      }
  }

  ScaffoldGraph =
    LoadScaffoldGraphFromCheckpoint(data->File_Name_Prefix, ckptNum, TRUE);

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
    FILE * fp = fopen(unitigIDFile, "r");
    char line[1024];
    VA_TYPE(ThreeIDs) * threeIDs = CreateVA_ThreeIDs(100);
    HashTable_AS * tidHT;
    VA_TYPE(ZLFScaffold) * zlfScaffolds = CreateVA_ZLFScaffold(10);
    ZLFScaffold zlfScaffold;
    ZLFContig zlfContig;
    int i, j, k;

    while(fgets(line, 1024, fp))
      {
        ThreeIDs tid;
        ChunkInstanceT * ci;
      
        tid.unitigID = atoi(line);
        ci = GetGraphNode(ScaffoldGraph->CIGraph, tid.unitigID);
        tid.contigID = ci->info.CI.contigID;
        tid.scaffoldID = ci->scaffoldID;

        AppendVA_ThreeIDs(threeIDs, &tid);
      }
    fclose(fp);
    qsort(GetVA_ThreeIDs(threeIDs, 0),
          GetNumVA_ThreeIDs(threeIDs),
          sizeof(ThreeIDs),
          (int (*) (const void *, const void *)) ThreeIDsCompare);
          
    zlfScaffold.id = -1;
    zlfContig.id = -1;
    
    for(i = 0; i < GetNumVA_ThreeIDs(threeIDs); i++)
      {
        ThreeIDs * tidp = GetVA_ThreeIDs(threeIDs, i);

        if(tidp->scaffoldID == zlfScaffold.id)
          {
            // continue populating zlfScaffold
            if(tidp->contigID == zlfContig.id)
              {
                // continue populating the contig
                AppendNewMAToContig(&zlfContig, tidp);
              }
            else
              {
                // same scaffold, different contig
                AppendNewContigToScaffold(&zlfScaffold, &zlfContig, tidp);
              }
          }
        else
          {
            // new scaffold (& new contig)
            AppendNewScaffold(zlfScaffolds, &zlfScaffold, &zlfContig, tidp);
          }
      }

    // append MA pointers to hashtable for quick loading from protoIO file
    tidHT = CreateHashTable_int32_AS(GetNumVA_ThreeIDs(threeIDs));
    for(i = 0; i < GetNumVA_ZLFScaffold(zlfScaffolds); i++)
      {
        ZLFScaffold * zlfsp = GetVA_ZLFScaffold(zlfScaffolds, i);
      
        for(j = 0; j < GetNumVA_ZLFContig(zlfsp->zlfContigs); j++)
          {
            ZLFContig * zlfcp = GetVA_ZLFContig(zlfsp->zlfContigs, j);

            for(k = 0; k < GetNumVA_MultiAlignT(zlfcp->zlfUMAs); k++)
              {
                MultiAlignT * map = GetVA_MultiAlignT(zlfcp->zlfUMAs, k);

                if(LookupInHashTable_AS(tidHT,
                                        (void *) &(map->id),
                                        sizeof(map->id)))
                  {
                    fprintf(stderr, "ERROR: unitig " F_CID " listed multiple times!\n",
                            map->id);
                    exit(-1);
                  }
                else
                  {
                    InsertInHashTable_AS(tidHT,
                                         (void *) &(map->id),
                                         sizeof(map->id),
                                         (void *) map);
                  }
              }
          }
      }

    // now read MultiAlignT store
    fp = fopen(maFile, "r");
    {
      MultiAlignT ma;

      while(ReadMAFromFile(&ma, fp))
        {
          MultiAlignT * map = LookupInHashTable_AS(tidHT,
                                                   (void *) &(ma.id),
                                                   sizeof(ma.id));
          if(map != NULL)
            {
              assert(map->id == ma.id);

              // copy the fields - the variable arrays will persist
          
              CopyMultiAlignT(map, &ma);
            }
        }
      fclose(fp);

      // fix zlf contigs
      FixZLFContigs(zlfScaffolds, TRUE, FALSE);
      CheckpointScaffoldGraph(ScaffoldGraph, -1);
    }
  }
  return 0;
}

