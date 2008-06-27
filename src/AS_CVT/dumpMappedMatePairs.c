
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

/*************************************************************************/
/* Local include files */
/*************************************************************************/

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_Hash.h"
#include "AS_PER_gkpStore.h"

#include "ASMData.h"

VA_DEF(CDS_UID_t);

int main(int argc, char ** argv)
{
  char * asmStorePath = NULL;
  char * mapStorePath = NULL;
  CDS_UID_t uid = 0;
  VA_TYPE(CDS_UID_t) * uids = CreateVA_CDS_UID_t(10);

  // parse command line
  {
    int ch,errflg=FALSE;
    while (!errflg && ((ch = getopt(argc, argv, "s:m:u:")) != EOF))
    {
      switch(ch)
      {
        case 's':
          asmStorePath = optarg;
          break;
        case 'm':
          mapStorePath = optarg;
          break;
        case 'u':
          uid = STR_TO_UID(optarg, NULL, 10);
          AppendVA_CDS_UID_t(uids, &uid);
          break;
        default:
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	  errflg++;
          break;
      }
    }
  }

  // check command line requirements
  if(asmStorePath == NULL ||
     mapStorePath == NULL)
  {
    fprintf(stderr,
            "Usage: %s [-s assemblyStore] [-m mapStore] [-u chromUID]\n"
            "\t-u is optional. default is all chromosome/segments\n",
            argv[0]);
    exit(1);
  }

  {
    AssemblyStore * asmStore;
    MapStore * mapStore;
    CloneData * cd;
    ASM_CHRRecord chr;
    char flags[256];

    asmStore = OpenReadOnlyAssemblyStore(asmStorePath);
    mapStore = OpenReadOnlyMapStore(mapStorePath);

    memset(flags, 0, 256 * sizeof(char));
    flags[AS_INNIE] = 1;
    flags[AS_NORMAL] = 1;
    flags[AS_ANTI] = 1;
    flags[AS_OUTTIE] = 1;
    flags['E'] = 0;
    flags['M'] = 0;

    if(GetNumVA_CDS_UID_t(uids) == 0)
    {
      int32 i;
      for(i = 1; i <= getNumASM_CHRs(mapStore->chrStore); i++)
      {
        char fname[1024];
        FILE * fo;
        getASM_CHRStore(mapStore->chrStore, i, &chr);
        cd = GetChromosomeCloneData(asmStore, mapStore, chr.uid);

        sprintf(fname, "%03" F_UIDP ".mmps.txt", chr.uid);
        fo = fopen(fname, "w");
        if(fo == NULL)
        {
          fprintf(stderr, "Failed to open file %s for writing.\n", fname);
          exit(1);
        }
        fprintf(stderr,
                "Writing chromosome UID " F_UID " mapped mate pairs to %s\n",
                chr.uid, fname);
        PrintCloneData(chr.uid, cd, flags, fo);
        fclose(fo);
        DeleteCloneData(cd);
      }
    }
    else
    {
      int32 i;
      for(i = 0; i < GetNumVA_CDS_UID_t(uids); i++)
      {
        CDS_UID_t * uidp = GetVA_CDS_UID_t(uids, i);
        cd = GetChromosomeCloneData(asmStore, mapStore, *uidp);
        PrintCloneData(*uidp, cd, flags, stdout);
        DeleteCloneData(cd);
      }
    }

    CloseAssemblyStore(asmStore);
    CloseMapStore(mapStore);
  }
  return 0;
}
