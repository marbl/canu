
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
static char CM_ID[] = "$Id: InstrumentCheckpoint.c,v 1.9 2007-05-14 09:27:11 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "Instrument_CGW.h"

void Usage(char * prog_name)
{
  fprintf(stderr, "Usage: %s\n", prog_name);
  fprintf(stderr,
          "required:\n"
          "<-a checkPointName>    assembly prefix\n"
          "<-n checkPointNumber>  checkpoint number\n"
          "<-g gkpStore>          gatekeeper store name\n"
          "optional:\n"
          "<-i>                   don't examine intra-contig mate pairs\n"
          "<-o>                   don't examine inter-contig mate pairs\n"
          "<-b>                   don't identify breakpoints\n"
          "<-j>                   just dump scaffold gap sizes\n"
          "<-l #>                 lower limit on scaffold size to instrument\n"
          "<-u #>                 upper limit on scaffold size to instrument\n"
          "<-v verbosity>         verbosity level:\n");
  fprintf(stderr, "\tverbosity levels (default is 3):\n"
          "\t1) Scaffold graph & scaffold summary,\n"
          "\t2) + contig summary,\n"
          "\t3) + stats for each scaffold,\n"
          "\t4) + stats for each contig,\n"
          "\t5) + details of each unhappy fragment (not implemented yet)\n");
  fprintf(stderr,
          "optional (may be specified multiple times):\n"
          "<-s IID>               scaffold IID to instrument\n"
          "<-S file>              file listing scaffold IIDs\n"
          "<-c IID>               contig IID to instrument\n"
          "<-C file>              file listing contig IIDs\n"
          );
  
  exit(1);
}


int AddFileIIDsToList(char * filename, VA_TYPE(CDS_CID_t) * list)
{
  FILE * fp;
  char line[1024];

  // open the file that lists IIDs, one per line
  if((fp = fopen(filename, "r")) == NULL)
    {
      fprintf(stderr, "Failed to opten file %s for reading\n", filename);
      return 1;
    }

  // read each line
  while(fgets(line, 1024, fp))
    {
      // get the IID and add it to the list
      CDS_CID_t iid = atoi(line);
      AppendVA_CDS_CID_t(list, &iid);
    }

  // close the file & return success
  fclose(fp);
  return 0;
}


static int iid_compare( const CDS_CID_t * a, const CDS_CID_t * b )
{
  if( *a > *b ) return 1;
  if( *a == * b) return 0;
  return -1;
}


int main(int argc, char *argv[])
{
  int32  checkPointNumber = -1;
  ScaffoldGraphT * graph = NULL;
  InstrumenterVerbosity verbosity = InstrumenterVerbose3;
  VA_TYPE(CDS_CID_t) * scf_iids = NULL;
  VA_TYPE(CDS_CID_t) * ctg_iids = NULL;
  uint32 options = INST_OPT_ALL;
  CDS_CID_t lower_limit = 0;
  CDS_CID_t upper_limit = CDS_CID_MAX;
  int justGapSizes = 0;

  GlobalData = CreateGlobal_CGW();
  GlobalData->stderrc = stderr;
  GlobalData->File_Name_Prefix[0] = '\0';
  GlobalData->Gatekeeper_Store_Name[0] = '\0';

  // parse the command line
  {
    int ch, errflg = 0;
    
    optarg = NULL;
    /*
      a: assembly name
      n: checkpoint number
      g: gatekeeper store name
      v: verbosidty
      s: scaffold number(s)
      S: file with scaffold number(s)
      c: contig number(s)
      C: file with contig number(s)
      u: upper limit on scaffold size
      l: lower limit on scaffold size
    */
    while( !errflg &&
           ((ch = getopt( argc,
                          argv,
                          "a:n:g:iobjv:s:S:c:C:u:l:" )) != EOF) )
      {
        switch(ch)
          {
            case 'a':
              strcpy(GlobalData->File_Name_Prefix, optarg);
              break;
            case 'n':
              checkPointNumber = atoi(optarg);
              break;
            case 'g':
              strcpy(GlobalData->Gatekeeper_Store_Name, optarg);
              break;
            case 'i':
              options ^= INST_OPT_INTRA_MATES;
              break;
            case 'o':
              options ^= INST_OPT_INTER_MATES;
              break;
            case 'b':
              options ^= INST_OPT_BREAKPOINTS;
              break;
            case 'j':
              justGapSizes = 1;
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
                    fprintf(stderr, "Only verbosity levels 1 - 5 supported\n");
                    Usage(argv[0]);
                    break;
                }
              break;
            case 's':
            case 'S':
              {
                if(scf_iids == NULL)
                  {
                    scf_iids = CreateVA_CDS_CID_t(100);
                    if(scf_iids == NULL)
                      {
                        fprintf(stderr,
                                "Failed to allocate scaffold IID variable array\n");
                        return 1;
                      }
                  }
                if(ch == 's')
                  {
                    CDS_CID_t s_iid = atoi(optarg);
                    AppendVA_CDS_CID_t(scf_iids, &s_iid);
                  }
                else
                  {
                    if(AddFileIIDsToList(optarg, scf_iids))
                      return 1;
                  }
              }
              break;
            case 'c':
            case 'C':
              {
                if(ctg_iids == NULL)
                  {
                    ctg_iids = CreateVA_CDS_CID_t(100);
                    if(ctg_iids == NULL)
                      {
                        fprintf(stderr,
                                "Failed to allocate scaffold IID variable array\n");
                        return 1;
                      }
                  }
                if(ch == 'c')
                  {
                    CDS_CID_t c_iid = atoi(optarg);
                    AppendVA_CDS_CID_t(ctg_iids, &c_iid);
                  }
                else
                  {
                    if(AddFileIIDsToList(optarg, ctg_iids))
                      return 1;
                  }
              }
              break;
            case 'l':
              lower_limit = atoi(optarg);
              break;
            case 'u':
              upper_limit = atoi(optarg);
              break;
            default:
              Usage(argv[0]);
              break;
          }
      }
  }

  if(checkPointNumber < 0 ||
     strlen(GlobalData->File_Name_Prefix) == 0 ||
     strlen(GlobalData->Gatekeeper_Store_Name) == 0)
    Usage(argv[0]);

  if(lower_limit > upper_limit)
    {
      fprintf(stderr, "lower limit may not exceed upper limit!\n");
      Usage(argv[0]);
    }

  if(justGapSizes)
    verbosity = InstrumenterSilent;

  // load the checkpoint
  fprintf( stderr, "Loading scaffold graph %s, checkpoint %d.\n",
           GlobalData->File_Name_Prefix, checkPointNumber);
  graph = LoadScaffoldGraphFromCheckpoint(GlobalData->File_Name_Prefix,
                                          checkPointNumber,
                                          FALSE );
  if(graph == NULL) {
    fprintf(stderr, "Failed to load graph from checkpoint!\n");
    return 1;
  }

  if(!scf_iids && !ctg_iids)
    {
#if 0
        ContigInstrumenter * ci = CreateContigInstrumenter(graph, options);
        ScaffoldInstrumenter * si = CreateScaffoldInstrumenter(graph, options);
        GraphNodeIterator scaffolds;
        CIScaffoldT * scaff;
        CDS_CID_t i = 0;

        // make sure verbosity is high enough to generate output
        if(verbosity < InstrumenterVerbose3)
        verbosity = InstrumenterVerbose3;

        // loop over all scaffolds in the graph
        InitGraphNodeIterator(&scaffolds,
        graph->ScaffoldGraph,
        GRAPH_NODE_DEFAULT);
        while(NULL != (scaff = NextGraphNodeIterator(&scaffolds)))
        {
        if(scaff->flags.bits.isDead == FALSE && scaff->type == REAL_SCAFFOLD)
        {
        i++;
        fprintf(stderr, "\r" F_CID "\t" F_CID "\t%15.0fbp",
        i, scaff->id, scaff->bpLength.mean);
        if(InstrumentScaffold(graph, NULL, scaff, si, ci, verbosity, stdout))
        {
        fprintf(stderr,
        "Failed to instrument scaffold " F_CID "\n",scaff->scaffoldID);
        return 1;
        }
        }
        }
        DestroyScaffoldInstrumenter(si);
        DestroyContigInstrumenter(ci);
#endif

      ScaffoldGraphInstrumenter * sgi =
        CreateScaffoldGraphInstrumenter(graph, options);
      // no scaffold or contig IIDs specified, instrument the whole graph
      if(InstrumentScaffoldGraph(graph, sgi, lower_limit, upper_limit,
                                 verbosity, stdout))
        {
          fprintf(GlobalData->stderrc, "Failed to instrument scaffold graph\n");
          return 1;
        }
      DestroyScaffoldGraphInstrumenter(sgi);

    }
  else
    {
      int i;
      ContigInstrumenter * ci = CreateContigInstrumenter(graph, options);
    
      if(ctg_iids)
        {
          // sort the contig IIDs - may make things go faster
          if(GetNumVA_CDS_CID_t(ctg_iids) > 1)
            qsort(ctg_iids->Elements,
                  GetNumVA_CDS_CID_t(ctg_iids),
                  sizeof(CDS_CID_t),
                  (int (*) (const void *, const void *)) iid_compare);

          // loop over contig IIDs & instrument & print
          for(i = 0; i < GetNumVA_CDS_CID_t(ctg_iids); i++)
            {
              CDS_CID_t * c_iid = GetVA_CDS_CID_t(ctg_iids, i);
              ChunkInstanceT * contig = GetGraphNode(graph->RezGraph, *c_iid);
              if(contig != NULL)
                {
                  InstrumentContig(graph, NULL, NULL, contig, ci,
                                   0.0f, contig->bpLength.mean);
                  PrintContigInstrumenter(graph, ci, verbosity, "", stdout);
                }
              else
                fprintf(stderr, "Contig " F_CID " doesn't exist.\n", *c_iid);
            }
          DeleteVA_CDS_CID_t(ctg_iids);
        }
    
      if(scf_iids)
        {
          ScaffoldInstrumenter * si = CreateScaffoldInstrumenter(graph, options);
      
          // sort the scaffold IIDs - may make things go faster
          if(GetNumVA_CDS_CID_t(scf_iids) > 1)
            qsort(scf_iids->Elements,
                  GetNumVA_CDS_CID_t(scf_iids),
                  sizeof(CDS_CID_t),
                  (int (*) (const void *, const void *)) iid_compare);

          // make sure verbosity is high enough to generate output
          if(verbosity < InstrumenterVerbose3)
            verbosity = InstrumenterVerbose3;

          // loop over scaffold IIDs and instrument/print
          for(i = 0; i < GetNumVA_CDS_CID_t(scf_iids); i++)
            {
              CDS_CID_t * s_iid = GetVA_CDS_CID_t(scf_iids, i);
              CIScaffoldT * scaffold = GetGraphNode(graph->ScaffoldGraph, *s_iid);
              if(scaffold != NULL)
                InstrumentScaffold(graph, scaffold, si, verbosity, stdout);
              else
                fprintf(stderr, "Scaffold " F_CID " doesn't exist.\n", *s_iid);
            }
          DestroyScaffoldInstrumenter(si);
          DeleteVA_CDS_CID_t(scf_iids);
        }
      DestroyContigInstrumenter(ci);
    }

  return 0;
}
