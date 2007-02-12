
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
static char CM_ID[] = "$Id: TraceNodes.c,v 1.6 2007-02-12 22:16:56 brianwalenz Exp $";

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
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "Instrument_CGW.h"


typedef struct
{
  CDS_CID_t iid;
  VA_TYPE(CDS_CID_t) * utg_iids;
  VA_TYPE(CDS_CID_t) * scf_iids;
} Tracer;

VA_DEF(Tracer);

void Usage(char * message, char * prog_name)
{
  fprintf(stderr, "\n**** %s ****\n\n", message);
  fprintf(stderr, "Usage: %s\n", prog_name);
  fprintf(stderr,
          "required:\n"
          "<-a name>     assembly name prefix\n"
          "<-f frgStore> frgStore\n"
          "<-g gkpStore> gkpStore\n"
          "<-l #>        latest checkpoint number\n"
          "optional:\n"
          "<-e #>        earliest checkpoint number (default = 4)\n"
          "<-F>          go forward from earliest to latest checkpoints\n"
          "<-i>          instrument scaffolds\n"
          "<-n>          don't instrument intra-contig mates\n"
          "<-x>          don't instrument inter-contig mates\n"
          "<-b>          don't instrument breakpoints\n");
  
  fprintf(stderr,
          "<-v verbosity>         verbosity level:\n"
          "\tverbosity levels (default is 3):\n"
          "\t1) Scaffold graph & scaffold summary (if doing entire graph),\n"
          "\t2)  + contig summary (if doing entire graph),\n"
          "\t3) (+) stats for each scaffold,\n"
          "\t4)  + stats for each contig,\n"
          "\t5)  + details of each unhappy fragment (not implemented yet)\n");
  
  fprintf(stderr,
          "optional (may be specified multiple times):\n"
          "<-u IID>      unitig IID to trace containing scaffold(s)\n"
          "<-U file>     file listing unitig IIDs\n"
          "<-c IID>      contig IID to trace containing scaffold(s)\n"
          "<-C file>     file listing contig IIDs\n"
          "<-s IID>      scaffold IID to trace\n"
          "<-S file>     file listing scaffold IIDs\n"
          );
  exit(1);
}


int AddFileIIDsToList(char * filename, VA_TYPE(Tracer) * list)
{
  FILE * fp;
  char line[1024];
  Tracer tracer;

  // open the file that lists IIDs, one per line
  if((fp = fopen(filename, "r")) == NULL)
    {
      fprintf(stderr, "Failed to open file %s for reading\n", filename);
      return 1;
    }

  tracer.utg_iids = NULL;
  tracer.scf_iids = NULL;
  // read each line
  while(fgets(line, 1024, fp))
    {
      // get the IID and add it to the list
      tracer.iid = atoi(line);
      AppendVA_Tracer(list, &tracer);
    }

  // close the file & return success
  fclose(fp);
  return 0;
}


static int tracer_iid_compare(const Tracer * a, const Tracer * b)
{
  if( a->iid > b->iid ) return 1;
  if( a->iid == b->iid) return 0;
  return -1;
}


static int iid_compare(const CDS_CID_t * a, const CDS_CID_t * b)
{
  return (int) (*a - *b);
}


static int iid_pair_compare(const Tracer * a, const Tracer * b)
{
  CDS_CID_t * a_iid = GetVA_CDS_CID_t(a->scf_iids, 0);
  CDS_CID_t * b_iid = GetVA_CDS_CID_t(b->scf_iids, 0);
  return (int)(*a_iid - *b_iid);
}


VA_TYPE(Tracer) * PopulateScaffoldIIDs(ScaffoldGraphT * graph)
{
  VA_TYPE(Tracer) * scfs = NULL;
  Tracer tracer;

  scfs = CreateVA_Tracer(GetNumCIScaffoldTs(graph->CIScaffolds));

  if(scfs)
    {
      GraphNodeIterator scaffolds;
      CIScaffoldT * scaff;

      tracer.utg_iids = NULL;
      tracer.scf_iids = NULL;
      InitGraphNodeIterator(&scaffolds,
                            graph->ScaffoldGraph,
                            GRAPH_NODE_DEFAULT);
      while(NULL != (scaff = NextGraphNodeIterator(&scaffolds)))
        {
          if(scaff->flags.bits.isDead == FALSE && scaff->type == REAL_SCAFFOLD)
            {
              tracer.iid = scaff->id;
              AppendVA_Tracer(scfs, &tracer);
            }
        }
    }
  return scfs;
}


void AddContigUnitigs(ScaffoldGraphT * graph,
                      Tracer * scf,
                      CDS_CID_t contig_id)
{
  ContigTIterator unitig_iterator;
  ChunkInstanceT * unitig;
  
  InitContigTIterator(graph, contig_id, TRUE, FALSE, &unitig_iterator);
  while((unitig = NextContigTIterator(&unitig_iterator)) != NULL)
    {
      // ignore surrogates
      if(!unitig->flags.bits.isStoneSurrogate &&
         !unitig->flags.bits.isWalkSurrogate)
        {
          AppendVA_CDS_CID_t(scf->utg_iids, &(unitig->id));
        }
    }
}


int PopulateScaffoldUnitigIIDs(ScaffoldGraphT * graph, Tracer * scf)
{
  CIScaffoldTIterator CIsTemp;
  CIScaffoldT * scaffold = GetGraphNode(graph->ScaffoldGraph, scf->iid);
  int i = 0;

  if(scaffold == NULL)
    {
      fprintf(stderr, "Scaffold ID " F_CID " isn't in checkpoint!\n", scf->iid);
      return 1;
    }
  
  // add contig's unitigs to list
  scf->utg_iids = CreateVA_CDS_CID_t(100);
  if(scf->utg_iids == NULL)
    {
      fprintf(stderr, "Failed to create temp array of unitig iids\n");
      return 1;
    }
  
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIsTemp);
  for(i = 0; NextCIScaffoldTIterator(&CIsTemp); i++)
    {
      if(i == 0)
        AddContigUnitigs(graph, scf, CIsTemp.curr);
      if(CIsTemp.next != NULLINDEX && CIsTemp.next != CIsTemp.curr)
        AddContigUnitigs(graph, scf, CIsTemp.next);
    }
  return 0;
}


ScaffoldGraphT * LoadCheckpoint(int i, int read_write)
{
  ScaffoldGraphT * graph = NULL;
  
  fprintf(stderr, "Loading scaffold graph %s, checkpoint %d.\n",
          GlobalData->File_Name_Prefix, i);
  fprintf(stdout,
          "********************************************************\n");
  fprintf(stdout,
          "**************** Checkpoint %10d *****************\n", i);
  fprintf(stdout,
          "********************************************************\n");
  graph = LoadScaffoldGraphFromCheckpoint(GlobalData->File_Name_Prefix,
                                          i,
                                          FALSE);
  if(graph == NULL)
    {
      fprintf(stderr, "Failed to load graph from checkpoint %d!\n", i);
      return NULL;
    }
  return graph;
}


VA_TYPE(Tracer) * GetScaffoldIIDsFromNodeIIDs(ScaffoldGraphT * graph,
                                              VA_TYPE(Tracer) * nodes,
                                              int is_unitig)
{
  VA_TYPE(Tracer) * scfs1 = CreateVA_Tracer(GetNumVA_Tracer(nodes));
  VA_TYPE(Tracer) * scfs2 = NULL;

  if(scfs1 != NULL)
    {
      CDS_CID_t i;
      Tracer tracer;
      tracer.utg_iids = NULL;
      tracer.scf_iids = NULL;
      for(i = 0; i < GetNumVA_Tracer(nodes); i++)
        {
          Tracer * tracerp = GetVA_Tracer(nodes, i);
          ChunkInstanceT * node =
            GetGraphNode((is_unitig) ? graph->CIGraph : graph->ContigGraph,
                         tracerp->iid);

          if(node->scaffoldID != NULLINDEX)
            {
              tracer.iid = node->scaffoldID;
              AppendVA_Tracer(scfs1, &tracer);
            }
        }
    
      if(GetNumVA_Tracer(scfs1) > 1)
        {
          int j;
          CDS_CID_t last_iid = NULLINDEX;

          scfs2 = CreateVA_Tracer(GetNumVA_Tracer(nodes));
          if(scfs2 == NULL)
            {
              fprintf(stderr, "Failed to allocate scaffold IIDs array!\n");
              return NULL;
            }
          qsort(GetVA_Tracer(scfs1, 0),
                GetNumVA_Tracer(scfs1),
                sizeof(Tracer),
                (int (*) (const void *, const void *)) tracer_iid_compare);

          for(j = 0; j < GetNumVA_Tracer(scfs1); j++)
            {
              Tracer * tracerp = GetVA_Tracer(scfs1, j);
              // fprintf(stdout, F_CID "\n", tracerp->iid);
              if(tracerp->iid != last_iid)
                AppendVA_Tracer(scfs2, tracerp);
              last_iid = tracerp->iid;
            }
          DeleteVA_Tracer(scfs1);
        }
      else
        scfs2 = scfs1;
    }

  return scfs2;
}
      

int ProcessCheckpoint(int checkpoint,
                      ScaffoldGraphInstrumenter * sgi,
                      ScaffoldInstrumenter * si,
                      VA_TYPE(Tracer) * scfs,
                      int instrument,
                      InstrumenterVerbosity verbosity)
{
  int j, k;
  ScaffoldGraphT * graph = NULL;
  
  // load the checkpoint
  if((graph = LoadCheckpoint(checkpoint, FALSE)) == NULL)
    return 1;

  // loop over all scaffolds
  for(j = 0; j < GetNumVA_Tracer(scfs); j++)
    {
      Tracer * scf = GetVA_Tracer(scfs, j);
      CDS_CID_t last_iid = NULLINDEX;
    
      if(scf->scf_iids == NULL)
        {
          scf->scf_iids = CreateVA_CDS_CID_t(100);
          if(scf->scf_iids == NULL)
            {
              fprintf(stderr, "Failed to allocate temp array of scf iids\n");
              return 1;
            }
        }
      else
        ResetVA_CDS_CID_t(scf->scf_iids);
    
      // loop over all unitigs in scaffold & look up current scaffolds
      for(k = 0; k < GetNumVA_CDS_CID_t(scf->utg_iids); k++)
        {
          CDS_CID_t * utg_iid = GetVA_CDS_CID_t(scf->utg_iids, k);
          ChunkInstanceT * unitig = GetGraphNode(graph->CIGraph, *utg_iid);
      
          if(!unitig->flags.bits.isStoneSurrogate &&
             !unitig->flags.bits.isWalkSurrogate)
            AppendVA_CDS_CID_t(scf->scf_iids, &(unitig->scaffoldID));
        }
    
      // sort scf_iids to go in order & avoid duplicates
      if(GetNumVA_CDS_CID_t(scf->scf_iids) > 1)
        qsort(GetVA_CDS_CID_t(scf->scf_iids, 0),
              GetNumVA_CDS_CID_t(scf->scf_iids),
              sizeof(CDS_CID_t),
              (int (*) (const void *, const void *)) iid_compare);
    
      // print prior scaffold IDs
      fprintf(stdout, "Checkpoint %d, scaffold " F_CID " = ",
              checkpoint, scf->iid);
      for(k = 0; k < GetNumVA_CDS_CID_t(scf->scf_iids); k++)
        {
          CDS_CID_t * scf_iid = GetVA_CDS_CID_t(scf->scf_iids, k);
          if(*scf_iid != last_iid)
            {
              fprintf(stdout, "\t" F_CID, *scf_iid);
            }
          last_iid = *scf_iid;
        }
      fprintf(stdout, "\n");
    
      if(instrument && !sgi)
        {
          last_iid = NULLINDEX;
      
          for(k = 0; k < GetNumVA_CDS_CID_t(scf->scf_iids); k++)
            {
              CDS_CID_t * scf_iid = GetVA_CDS_CID_t(scf->scf_iids, k);
              if(*scf_iid != last_iid)
                {
                  CIScaffoldT * scaff = GetCIScaffoldT(graph->CIScaffolds,
                                                       *scf_iid);
                  InstrumentScaffold(graph, scaff, si, verbosity, stdout);
                }
              last_iid = *scf_iid;
            }
        }
    }
  
  if(instrument && sgi)
    {
      InstrumentScaffoldGraph(graph, sgi,
                              0, CDS_CID_MAX,
                              verbosity, stdout);
    }
  
  // free up for next round...
  DestroyScaffoldGraph(graph);

  return 0;
}

int main(int argc, char *argv[])
{
  int32  latestCheckPointNumber = -1;
  int32  earliestCheckPointNumber = 4;
  ScaffoldGraphT * graph = NULL;
  VA_TYPE(Tracer) * scfs = NULL;
  VA_TYPE(Tracer) * utgs = NULL;
  VA_TYPE(Tracer) * ctgs = NULL;
  int i, j;
  int instrument = 0;
  int go_forward = 0;
  InstrumenterVerbosity verbosity = InstrumenterVerbose3;
  cds_uint32 options = INST_OPT_ALL;
  ScaffoldInstrumenter * si = NULL;
  ScaffoldGraphInstrumenter * sgi = NULL;

  GlobalData = CreateGlobal_CGW();
  GlobalData->stderrc = stderr;
  GlobalData->File_Name_Prefix[0] = '\0';
  GlobalData->Gatekeeper_Store_Name[0] = '\0';
  
  // parse the command line
  {
    int ch, errflg = 0;
    
    optarg = NULL;
    /*
      "<-a name>     assembly name prefix\n"
      "<-g gkpStore> gkpStore\n"
      "<-l #>        latest checkpoint number\n"
      "optional:\n"
      "<-e #>        earliest checkpoint number (default = 4)\n"
      "<-F>          go forward from earliest to latest checkpoints\n"
      "<-i>          instrument scaffolds\n"
      "<-n>          don't instrument intra-contig mates\n"
      "<-x>          don't instrument inter-contig mates\n"
      "<-b>          don't instrument breakpoints\n");
  
      "\tverbosity levels (default is 3):\n"
      "\t1) Scaffold graph & scaffold summary (if doing entire graph),\n"
      "\t2)  + contig summary (if doing entire graph),\n"
      "\t3) (+) stats for each scaffold,\n"
      "\t4)  + stats for each contig,\n"
      "\t5)  + details of each unhappy fragment (not implemented yet)\n");
  
      "optional (may be specified multiple times):\n"
      "<-u IID>      unitig IID to trace containing scaffold(s)\n"
      "<-U file>     file listing unitig IIDs\n"
      "<-c IID>      contig IID to trace containing scaffold(s)\n"
      "<-C file>     file listing contig IIDs\n"
      "<-s IID>      scaffold IID to trace\n"
      "<-S file>     file listing scaffold IIDs\n"
    */
    while(!errflg &&
          ((ch = getopt( argc,
                         argv,
                         "a:f:g:l:e:Finxbv:u:U:c:C:s:S:" )) != EOF))
      {
        switch(ch)
          {
            case 'a':
              strcpy(GlobalData->File_Name_Prefix, optarg);
              break;
            case 'g':
              strcpy(GlobalData->Gatekeeper_Store_Name, optarg);
              break;
            case 'l':
              latestCheckPointNumber = atoi(optarg);
              if(latestCheckPointNumber < 4)
                Usage("No scaffolds prior to checkpoint 4!", argv[0]);
              break;
            case 'e':
              earliestCheckPointNumber = atoi(optarg);
              if(earliestCheckPointNumber < 4)
                Usage("No scaffolds prior to checkpoint 4!", argv[0]);
              break;
            case 'F':
              go_forward = 1;
              break;
            case 'i':
              instrument = 1;
              break;
            case 'n':
              options ^= INST_OPT_INTRA_MATES;
              break;
            case 'x':
              options ^= INST_OPT_INTER_MATES;
              break;
            case 'b':
              options ^= INST_OPT_BREAKPOINTS;
              break;
            case 'v':
              switch(atoi(optarg))
                {
                  case 1:
                    verbosity = InstrumenterVerbose1;
                    break;
                  case 2:
                    verbosity = InstrumenterVerbose2;
                    break;
                  case 3:
                    verbosity = InstrumenterVerbose3;
                    break;
                  case 4:
                    verbosity = InstrumenterVerbose4;
                    break;
                  case 5:
                    verbosity = InstrumenterVerbose5;
                    break;
                  default:
                    Usage("Only verbosity levels 1 - 5 supported", argv[0]);
                    break;
                }
              break;
            case 'u':
            case 'U':
              if(utgs == NULL)
                {
                  utgs = CreateVA_Tracer(100);
                  if(utgs == NULL)
                    {
                      fprintf(stderr,
                              "Failed to allocate unitig IID variable array\n");
                      return 1;
                    }
                }
              if(ch == 'u')
                {
                  Tracer tracer;
                  tracer.iid = atoi(optarg);  // unitig IID!
                  tracer.utg_iids = NULL;
                  tracer.scf_iids = NULL;
                  AppendVA_Tracer(utgs, &tracer);
                }
              else
                {
                  if(AddFileIIDsToList(optarg, utgs))
                    return 1;
                }
              break;
            case 'c':
            case 'C':
              if(ctgs == NULL)
                {
                  ctgs = CreateVA_Tracer(100);
                  if(ctgs == NULL)
                    {
                      fprintf(stderr,
                              "Failed to allocate contig IID variable array\n");
                      return 1;
                    }
                }
              if(ch == 'c')
                {
                  Tracer tracer;
                  tracer.iid = atoi(optarg);  // contig IID!
                  tracer.utg_iids = NULL;
                  tracer.scf_iids = NULL;
                  AppendVA_Tracer(ctgs, &tracer);
                }
              else
                {
                  if(AddFileIIDsToList(optarg, ctgs))
                    return 1;
                }
              break;
            case 's':
            case 'S':
              if(scfs == NULL)
                {
                  scfs = CreateVA_Tracer(100);
                  if(scfs == NULL)
                    {
                      fprintf(stderr,
                              "Failed to allocate scaffold IID variable array\n");
                      return 1;
                    }
                }
              if(ch == 's')
                {
                  Tracer tracer;
                  tracer.iid = atoi(optarg);
                  tracer.utg_iids = NULL;
                  tracer.scf_iids = NULL;
                  AppendVA_Tracer(scfs, &tracer);
                }
              else
                {
                  if(AddFileIIDsToList(optarg, scfs))
                    return 1;
                }
              break;
            default:
              Usage("Flag not recognized", argv[0]);
              break;
          }
      }
  }
  
  if(latestCheckPointNumber < earliestCheckPointNumber)
    Usage("Latest checkpoint must be greater than earliest checkpoint",
          argv[0]);
  if(strlen(GlobalData->File_Name_Prefix) == 0)
    Usage("Please specify an assembly checkpoint name", argv[0]);
  if(strlen(GlobalData->Gatekeeper_Store_Name) == 0)
    Usage("Please specify a gatekeeper store name", argv[0]);

  // load the checkpoint
  if(go_forward)
    {
      if((graph = LoadCheckpoint(earliestCheckPointNumber, FALSE)) == NULL)
        return 1;
    }
  else
    {
      if((graph = LoadCheckpoint(latestCheckPointNumber, FALSE)) == NULL)
        return 1;
    }
  
  // if user supplied unitig IIDs, convert to scaffold IIDs
  if(utgs)
    {
      VA_TYPE(Tracer) * scfs2 = GetScaffoldIIDsFromNodeIIDs(graph, utgs, 1);
      DeleteVA_Tracer(utgs);
      if(scfs2 == NULL)
        {
          fprintf(stderr, "Failed to convert unitig IIDs to scaffold IIDs!\n");
          return 1;
        }
      if(scfs)
        {
          int k;
          // append
          for(k = 0; k < GetNumVA_Tracer(scfs2); k++)
            {
              Tracer * tracerp = GetVA_Tracer(scfs2, k);
              AppendVA_Tracer(scfs, tracerp);
            }
          DeleteVA_Tracer(scfs2);
        }
      else
        scfs = scfs2;
    }

  // if user supplied contig IIDs, convert to scaffold IIDs
  if(ctgs)
    {
      VA_TYPE(Tracer) * scfs2 = GetScaffoldIIDsFromNodeIIDs(graph, ctgs, 0);
      DeleteVA_Tracer(ctgs);
      if(scfs2 == NULL)
        {
          fprintf(stderr, "Failed to convert contig IIDs to scaffold IIDs!\n");
          return 1;
        }
      if(scfs)
        {
          int k;
          // append
          for(k = 0; k < GetNumVA_Tracer(scfs2); k++)
            {
              Tracer * tracerp = GetVA_Tracer(scfs2, k);
              AppendVA_Tracer(scfs, tracerp);
            }
          DeleteVA_Tracer(scfs2);
        }
      else
        scfs = scfs2;
    }

  if(!scfs)
    {
      sgi = CreateScaffoldGraphInstrumenter(graph, options);
      if(sgi == NULL)
        {
          fprintf(stderr, "Failed to create scaffold graph instrumenter\n");
          return 1;
        }
    
      // do ALL scaffold IIDs
      fprintf(stderr, "Tracing all scaffold IIDs\n");
    
      if((scfs = PopulateScaffoldIIDs(graph)) == NULL)
        {
          fprintf(stderr, "Failed to get list of all scaffold iids\n");
          return 1;
        }
    }
  
  // at this point, scfs will be non-null
  // sort the scaffold IIDs - may make things go faster
  if(GetNumVA_Tracer(scfs) > 1)
    qsort(GetVA_Tracer(scfs, 0),
          GetNumVA_Tracer(scfs),
          sizeof(Tracer),
          (int (*) (const void *, const void *)) tracer_iid_compare);

  // get the unitig IIDs that are part of each scaffold
  for(i = 0; i < GetNumVA_Tracer(scfs); i++)
    {
      PopulateScaffoldUnitigIIDs(graph, GetVA_Tracer(scfs, i));
    }

  // print scaffolds etc
  fprintf(stdout, "Tracing %sfrom ckp %d to %d\n",
          (instrument) ? "and instrumenting " : "",
          latestCheckPointNumber, earliestCheckPointNumber);
  for(i = 0; i < GetNumVA_Tracer(scfs); i++)
    {
      Tracer * scf = GetVA_Tracer(scfs, i);
      fprintf(stdout, F_CID "\n", scf->iid);
    }

  //  if(0)
  {
    // print scaffold->uid lists
    fprintf(stdout, "SCF_IID\tUTG_IIDs\n");
    for(i = 0; i < GetNumVA_Tracer(scfs); i++)
      {
        Tracer * scf = GetVA_Tracer(scfs, i);
        fprintf(stdout, F_CID, scf->iid);
        for(j = 0;
            scf->utg_iids != NULL && j < GetNumVA_CDS_CID_t(scf->utg_iids);
            j++)
          {
            CDS_CID_t * id = GetVA_CDS_CID_t(scf->utg_iids, j);
            fprintf(stdout, "\t" F_CID, *id);
          }
        fprintf(stdout, "\n");
      }
  }

  // if instrumenting, allocate instrumenters & instrument...
  if(instrument)
    {
      if((si = CreateScaffoldInstrumenter(graph, options)) == NULL)
        {
          fprintf(stderr, "Failed to create scaffold instrumenter\n");
          return 1;
        }
      for(j = 0; j < GetNumVA_Tracer(scfs); j++)
        {
          Tracer * scf = GetVA_Tracer(scfs, j);
          CIScaffoldT * scaff = GetCIScaffoldT(graph->CIScaffolds, scf->iid);
          InstrumentScaffold(graph, scaff, si, verbosity, stdout);
        }
    }
  DestroyScaffoldGraph(graph);
  
  // loop over scaffold IIDs and instrument/print
  if(go_forward)
    {
      for(i = earliestCheckPointNumber - 1; i >= latestCheckPointNumber; i++)
        {
          ProcessCheckpoint(i, sgi, si, scfs, instrument, verbosity);
        }
    }
  else
    {
      for(i = latestCheckPointNumber - 1; i >= earliestCheckPointNumber; i--)
        {
          ProcessCheckpoint(i, sgi, si, scfs, instrument, verbosity);
        }
    }
  
  // clean up
  if(si)
    DestroyScaffoldInstrumenter(si);
  if(sgi)
    DestroyScaffoldGraphInstrumenter(sgi);
  DeleteVA_Tracer(scfs);

  return 0;
}
