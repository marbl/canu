
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
#include  <stdlib.h>
#include  <stdio.h>
#include  <unistd.h>
#include  <assert.h>
#include  <sys/types.h>
#include  <string.h>

#include "AS_global.h"
#include "ASMData.h"


void Usage(char * prog, char * message)
{
  if(message)
    fprintf(stderr, "ERROR: %s\n", message);
  fprintf(stderr, "Usage: %s [-s asmStore] [-u unsatisfied type] [-n stddevs] [-f fragment type] [-S scaffoldUID] [-g] [-d]\n",
          prog);
  fprintf(stderr,
          "\t-s     assembly store\n"
          
          "\t-u     unsatisfied mate pair type\n"
          "\t\t     multiple types are allowed, and default is all\n"
          "\t\t     S = stretched mate pairs   (----->                   <-----)\n"
          "\t\t     C = compressed mate pairs  (-----> <-----)\n"
          "\t\t     N = normal mate pairs      (----->     ----->)\n"
          "\t\t     A = antinormal mate pairs  (<-----     <-----)\n"
          "\t\t     O = outtie mate pair left & right (<-----     ----->)\n"
          "\t\t     o = outtie mate pair middle       (      -----      )\n"
          
          "\t-n     number of stddevs from mean to use for slop\n"
          "\t\t     (default is 3.0)\n"
          
          "\t-f     type of fragments to dump:\n"
          "\t\t     multiple types are allowed, and default is all\n"
          "\t\t     R = Celera Read\n"
          "\t\t     X = External WGS read\n"
          "\t\t     T = Transposon library read\n"
          "\t\t     E = BAC end\n"
          "\t\t     L = Lightly shotgunned BAC\n"
          "\t\t     U = Unfinished BAC\n"
          "\t\t     F = Finished BAC (incomplete)\n"
          "\t\t     C = Full BAC (complete)\n"
          "\t\t     S = Sts\n"
          "\t\t     u = Unitig\n"
          "\t\t     c = Contig\n"
          "\t\t     B = BacTig\n"
          "\t\t     G = BGLII read\n"
          
          "\t-S     scaffold UID for which to dump stretched mates\n"
          "\t\t     multiple UIDs are allowed, and default is all\n"

          "\t-g     write output for gnu-plot\n"
          
          "\t-d     process degenerate scaffolds\n"
          "\t\t     default is to process only real scaffolds\n");
  exit(-1);
}


int main(int argc, char ** argv)
{
  char * inputStoreName = NULL;
  char fragTypes[256];
  int doDegenerates = FALSE;
  float32 numStddevs = 3.0f;
  int errflg=FALSE;
  char bpTypes[NumBreakpointTypes];
  int32 i;
  BreakpointType bpt;
  CDS_UID_t uid;
  VA_TYPE(CDS_UID_t) * scaffolds = CreateVA_CDS_UID_t(10);
  char message[256];
  ASM_OutputType ot = ASM_CliqueFinderOutput;
  
  memset(fragTypes, 0, 256);
  memset(bpTypes, 0, NumBreakpointTypes * sizeof(char));
  
  {
    int ch;
    while (!errflg && ((ch = getopt(argc, argv, "s:n:u:f:S:gd")) != EOF))
    {
      switch(ch) 
      {
	case 's':
          inputStoreName = optarg;
          break;
        case 'g':
          ot = ASM_GnuplotOutput;
          break;
        case 'n':
          numStddevs = atof(optarg);
          break;
        case 'S':
          uid = STR_TO_UID(optarg, NULL, 10);
          AppendVA_CDS_UID_t(scaffolds, &uid);
          break;
        case 'u':
          if(strlen(optarg) > 1)
            Usage(argv[0], "Unsatisfied type must be one character");
          if(optarg[0] != 'N' &&
             optarg[0] != 'S' &&
             optarg[0] != 'C' &&
             optarg[0] != 'O' &&
             optarg[0] != 'A' &&
             optarg[0] != 'o')
          {
            sprintf(message, "%s is not an unsatisfied type!", optarg);
            Usage(argv[0], message);
          }
          else
          {
            bpTypes[(int) optarg[0]] = 1;
          }
          break;
        case 'f':
          if(strlen(optarg) > 1)
            Usage(argv[0], "Fragment types must be specified separately");
          if(optarg[0] != 'R' &&
             optarg[0] != 'X' &&
             optarg[0] != 'T' &&
             optarg[0] != 'E' &&
             optarg[0] != 'L' &&
             optarg[0] != 'U' &&
             optarg[0] != 'F' &&
             optarg[0] != 'C' &&
             optarg[0] != 'S' &&
             optarg[0] != 'u' &&
             optarg[0] != 'c' &&
             optarg[0] != 'B' &&
             optarg[0] != 'G')
          {
            sprintf(message, "%s is not a fragment type!", optarg);
            Usage(argv[0], message);
          }
          else
          {
            fragTypes[(int) optarg[0]] = 1;
          }
          break;
        case 'd':
          doDegenerates = TRUE;
          break;
        default:
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	  errflg++;
          break;
      }
    }
  }

  for(bpt = 0; bpt < NumBreakpointTypes; bpt++)
    if(bpTypes[bpt] != 0) break;
  if(bpt == NumBreakpointTypes)
    memset(bpTypes, 1, NumBreakpointTypes * sizeof(char));
  
  for(i = 0; i < 256; i++)
    if(fragTypes[i] != 0) break;
  if(i == 256)
    memset(fragTypes, 1, 256);

  if(errflg != FALSE || inputStoreName == NULL || numStddevs <= 0.f)
  {
    Usage(argv[0], "Missing or invalid parameter");
  }

  {
    AssemblyStore * asmStore;
    VA_TYPE(ASM_Quad) * quads;
    int32 maxI;
    ASM_SCFRecord scf;
    ASM_DSCRecord dsc;
    CloneData * cd;
    FILE * fout;
    char fname[1024];
    
    asmStore = OpenReadOnlyAssemblyStore(inputStoreName);
    
    if(GetNumVA_CDS_UID_t(scaffolds) == 0)
    {
      // do all scaffolds
      maxI = getNumASM_SCFs(asmStore->scfStore);
      for(i = 1; i <= maxI; i++)
      {
        cd = GetScaffoldCloneData(asmStore, i, TRUE, FALSE);
        getASM_SCFStore(asmStore->scfStore, i, &scf);
        for(bpt = 0; bpt < NumBreakpointTypes; bpt++)
        {
          if(bpTypes[bpt] != 0)
          {
            quads = IdentifyBadMateQuads(asmStore, cd, fragTypes,
                                         bpt, numStddevs);
            if(GetNumVA_ASM_Quad(quads) > 0)
            {
              CreateBMPFilename(fname, scf.uid, fragTypes, bpt, numStddevs);
              fout = fopen(fname, "w");
              PrintQuadrilaterals(quads, ot, fout);
              fclose(fout);
            }
            DeleteVA_ASM_Quad(quads);
          }
        }
        DeleteCloneData(cd);
      }
      if(doDegenerates)
      {
        maxI = getNumASM_DSCs(asmStore->dscStore);
        for(i = 1; i <= maxI; i++)
        {
          cd = GetScaffoldCloneData(asmStore, i, TRUE, TRUE);
          getASM_DSCStore(asmStore->dscStore, i, &dsc);
          for(bpt = 0; bpt < NumBreakpointTypes; bpt++)
          {
            if(bpTypes[bpt] != 0)
            {
              quads = IdentifyBadMateQuads(asmStore, cd, fragTypes,
                                           bpt, numStddevs);
              if(GetNumVA_ASM_Quad(quads) > 0)
              {
                CreateBMPFilename(fname, dsc.uid, fragTypes, bpt, numStddevs);
                fout = fopen(fname, "w");
                PrintQuadrilaterals(quads, ot, fout);
                fclose(fout);
              }
              DeleteVA_ASM_Quad(quads);
            }
          }
          DeleteCloneData(cd);
        }
      }
    }
    else
    {
      for(i = 0; i < GetNumVA_CDS_UID_t(scaffolds); i++)
      {
        CDS_UID_t * myUID = GetVA_CDS_UID_t(scaffolds, i);
        PHashValue_AS value;
        if(HASH_FAILURE != LookupInPHashTable_AS(asmStore->hashTable,
                                                 ASM_UID_NAMESPACE,
                                                 *myUID, &value))
        {
          if(value.type == AS_IID_SCF)
          {
            cd = GetScaffoldCloneData(asmStore, value.IID, TRUE, FALSE);
          }
          else if(value.type == AS_IID_DSC)
          {
            cd = GetScaffoldCloneData(asmStore, value.IID, TRUE, TRUE);
          }
          else
          {
            fprintf(stderr, F_UID " is a UID, but not for a scaffold\n",
                    *myUID);
            continue;
          }
          for(bpt = 0; bpt < NumBreakpointTypes; bpt++)
          {
            if(bpTypes[bpt] != 0)
            {
              quads = IdentifyBadMateQuads(asmStore, cd, fragTypes,
                                           bpt, numStddevs);
              if(GetNumVA_ASM_Quad(quads) > 0)
              {
                CreateBMPFilename(fname, *myUID, fragTypes, bpt, numStddevs);
                fout = fopen(fname, "w");
                PrintQuadrilaterals(quads, ot, fout);
                fclose(fout);
              }
              DeleteVA_ASM_Quad(quads);
            }
          }
          DeleteCloneData(cd);
        }
        else
        {
          fprintf(stderr, "Failed to find scaffold uid " F_UID " in store\n",
                  *myUID);
        }
      }
    }
  }
  return 0;
}
