
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

/*************************************************************************/
/* Local include files */
/*************************************************************************/

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_Hash.h"
#include "AS_PER_gkpStore.h"

#include "ASMData.h"

int main(int argc, char ** argv)
{
  char * asmStorePath = NULL;
  Scaffold_ID uid = 0;
  VA_TYPE(CDS_UID_t) * uids = CreateVA_CDS_UID_t(10);
  
  // parse command line
  {
    int ch,errflg=FALSE;
    while (!errflg && ((ch = getopt(argc, argv, "s:u:")) != EOF))
    {
      switch(ch) 
      {
        case 's':
          asmStorePath = optarg;
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
  if(asmStorePath == NULL)
  {
    fprintf(stderr,
            "Usage: %s [-s assemblyStore] [-u scaffoldUID]\n"
            "\t-u is optional. default is all scaffolds\n",
            argv[0]);
    exit(-1);
  }

  {
    AssemblyStore * asmStore;
    ASM_SCFRecord scf;
    
    asmStore = OpenReadOnlyAssemblyStore(asmStorePath);

    if(GetNumVA_CDS_UID_t(uids) == 0)
    {
      int32 i;
      for(i = 1; i <= getNumASM_SCFs(asmStore->scfStore); i++)
      {
        FILE * fp;
        char fname[2048];
        getASM_SCFStore(asmStore->scfStore, i, &scf);

        sprintf(fname, F_UID ".elsewhere.txt", scf.uid);
        fp = fopen(fname, "w");
        if(fp == NULL)
        {
          fprintf(stderr, "Failed to open file %s for writing.\n", fname);
          exit(1);
        }
        fprintf(stderr,
                "Writing scaffold " F_UID " inter-scaffold mate pairs to %s\n",
                scf.uid, fname);
        PrintScaffoldElsewheres(asmStore, i, 1, 0, fp);
        fclose(fp);
      }
    }
    else
    {
      // do nothing
    }

    CloseAssemblyStore(asmStore);
  }
  return 0;
}
