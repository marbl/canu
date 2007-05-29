
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
#include  <string.h>

#include "AS_global.h"
#include "ASMData.h"

VA_DEF(CDS_UID_t);

int main(int argc, char ** argv)
{
  char * inputStoreName = NULL;
  CDS_UID_t uid;
  int doGnuplotOutput = FALSE;
  int doAssemblyCoords = FALSE;
  VA_TYPE(CDS_UID_t) * uids = CreateVA_CDS_UID_t(10);
  {
    int ch,errflg=FALSE;
    while (!errflg && ((ch = getopt(argc, argv, "s:S:ga")) != EOF))
    {
      switch(ch) 
      {
	case 's':
          inputStoreName = optarg;
          break;
        case 'S':
          uid = STR_TO_UID(optarg, NULL, 10);
          AppendVA_CDS_UID_t(uids, &uid);
          break;
        case 'g':
          doGnuplotOutput = TRUE;
          break;
        case 'a':
          doAssemblyCoords = TRUE;
          break;
        default:
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	  errflg++;
          break;
      }
    }
  }
  
  if(inputStoreName == NULL)
  {
    fprintf(stderr, "Usage: %s [-s asmStore] [-S scaffold UID] [-ga]\n"
            "\t-S        specify as many as you like\n"
            "\t\t        default is all\n"
            "\t-g        write output for plotting in gnuplot\n"
            "\t-a        write assembly coordinates\n"
            "\t\t        default is fasta file coordinates\n",
            argv[0]);
    exit(-1);
  }

  {
    AssemblyStore * asmStore;
    int32 i;
    int32 maxI;

    asmStore = OpenReadOnlyAssemblyStore(inputStoreName);
    if(GetNumVA_CDS_UID_t(uids) == 0)
    {
      maxI = getNumASM_SCFs(asmStore->scfStore);
      for(i = 1; i <= maxI; i++)
      {
        PrintScaffoldContigCoordinates(asmStore, i, doGnuplotOutput,
                                       doAssemblyCoords, stdout);
      }
    }
    else
    {
      for(i = 0; i < GetNumVA_CDS_UID_t(uids); i++)
      {
        CDS_UID_t * myUID = GetVA_CDS_UID_t(uids, i);
        uint64      iid;
        uint32      typ;

        if (HASH_FAILURE != LookupInHashTable_AS(asmStore->hashTable, *myUID, 0, &iid, &typ))
        {
          if(typ == AS_IID_SCF)
          {
            PrintScaffoldContigCoordinates(asmStore, iid,
                                           doGnuplotOutput,
                                           doAssemblyCoords, stdout);
          }
          else
          {
            fprintf(stderr, F_UID " is a UID, but not for a scaffold\n",
                    *myUID);
            continue;
          }
        }
        else
        {
          fprintf(stderr, F_UID " is not a UID in this assembly store\n",
                  *myUID);
          continue;
        }
      }
    }
    CloseAssemblyStore(asmStore);
  }
  
  return 0;
}
