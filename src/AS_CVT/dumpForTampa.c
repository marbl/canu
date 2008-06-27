
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

void PrintScaffoldCloneData(AssemblyStore * asmStore,
                            CDS_IID_t iid,
                            char * assemblyName)
{
  FILE * fp;
  char fname[2048];
  ASM_SCFRecord scf;
  CloneData * cd;

  cd = GetScaffoldCloneData(asmStore, iid, TRUE, FALSE);
  getASM_SCFStore(asmStore->scfStore, iid, &scf);
  sprintf(fname, "%s_" F_UID "_intra.txt", assemblyName, scf.uid);
  fp = fopen(fname, "w");
  PrintIntraCloneData(scf.uid, cd, fp);
  fclose(fp);

  sprintf(fname, "%s_" F_UID "_inter.txt", assemblyName, scf.uid);
  fp = fopen(fname, "w");
  PrintScaffoldElsewheres(asmStore, iid, 1, 0, fp);
  fclose(fp);

  DeleteCloneData(cd);
}


int main(int argc, char ** argv)
{
  char * inputStoreName = NULL;
  char * assemblyName = NULL;
  int errflg=FALSE;
  CDS_UID_t uid;
  VA_TYPE(CDS_UID_t) * uids = CreateVA_CDS_UID_t(10);

  {
    int ch;
    while (!errflg && ((ch = getopt(argc, argv, "a:s:u:")) != EOF))
    {
      switch(ch)
      {
        case 'a':
          assemblyName = optarg;
          break;
	case 's':
          inputStoreName = optarg;
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

  if(errflg != FALSE || inputStoreName == NULL || assemblyName == NULL)
  {
    fprintf(stderr, "Usage: %s [-a assembly] [-s store] [-u uid]\n", argv[0]);
    fprintf(stderr,
            "\t-a assembly assembly name prefix to use in generating output\n"
            "\t              filenames\n"
            "\t-s store    assembly store with all the data\n"
            "\t-u uid      limit to this scaffold. multiple -u's okay\n"
            "\t              (default is all)\n");
    exit(1);
  }

  {
    AssemblyStore * asmStore;
    int32 i;
    int32 maxI;

    asmStore = OpenReadOnlyAssemblyStore(inputStoreName);

    // print libraries
    {
      FILE * fp;
      char fname[2048];
      sprintf(fname, "%sLibs.txt", assemblyName);
      fp = fopen(fname, "w");
      PrintLibraries(asmStore, fp);
      fclose(fp);
    }

    if(GetNumVA_CDS_UID_t(uids) == 0)
    {
      maxI = getNumASM_SCFs(asmStore->scfStore);
      for(i = 1; i <= maxI; i++)
      {
        PrintScaffoldCloneData(asmStore, i, assemblyName);
      }
    }
    else
    {
      for(i = 0; i < GetNumVA_CDS_UID_t(uids); i++)
      {
        CDS_UID_t * uidp = GetVA_CDS_UID_t(uids, i);
        uint64      iid;
        uint32      typ;

        if(HASH_FAILURE == LookupInHashTable_AS(asmStore->hashTable, *uidp, 0, &iid, &typ))
        {
          fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n", *uidp);
          assert(0);
        }
        if(typ == AS_IID_SCF)
        {
          PrintScaffoldCloneData(asmStore, iid, assemblyName);
        }
        else
        {
          // skip degenerates (in DSCStore)
        }
      }
    }
    CloseAssemblyStore(asmStore);
  }

  return 0;
}
