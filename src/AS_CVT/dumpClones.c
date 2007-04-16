
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
  char flags[256];
  int doDegenerates = FALSE;
  int errflg=FALSE;
  CDS_UID_t uid;
  VA_TYPE(CDS_UID_t) * uids = CreateVA_CDS_UID_t(10);
  
  memset(flags, 0, 256);
  {
    int ch;
    while (!errflg && ((ch = getopt(argc, argv, "s:dNIOAEDSCMu:")) != EOF))
    {
      switch(ch) 
      {
	case 's':
          inputStoreName = optarg;
          break;
        case 'N':
        case 'I':
        case 'O':
        case 'A':
        case 'E':
        case 'D':
        case 'S':
        case 'C':
        case 'M':
          flags[ch] = TRUE;
          break;
        case 'd':
          doDegenerates = TRUE;
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
  
  if(errflg != FALSE || inputStoreName == NULL)
  {
    fprintf(stderr, "Usage: %s [-s inputStore] [-dNIOAESC] [-u uid]\n", argv[0]);
    fprintf(stderr,
            "\t-d     process degenerate scaffolds, too\n"
            "\t-N     print normal mate pairs     (----->     ----->)\n"
            "\t-I     print innie mate pairs      (----->     <-----)\n"
            "\t-O     print outtie mate pairs     (<-----     ----->)\n"
            "\t-A     print antinormal mate pairs (<-----     <-----)\n"
            "\t-E     print if mate is in a different real scaffold\n"
            "\t-D     print if mate is in a degenerate scaffold\n"
            "\t-S     print if mate is in a multiply-placed surrogate\n"
            "\t-C     print if mate is chaff fragment\n"
            "\t-M     print if mate is missing from assembly\n"
            "\t-u uid limit to this scaffold. multiple -u's okay\n"
            "\t         (default is all)\n");
    exit(-1);
  }

  {
    AssemblyStore * asmStore;
    int32 i;
    int32 maxI;
    ASM_SCFRecord scf;
    ASM_DSCRecord dsc;
    CloneData * cd;

    asmStore = OpenReadOnlyAssemblyStore(inputStoreName);

    if(GetNumVA_CDS_UID_t(uids) == 0)
    {
      maxI = getNumASM_SCFs(asmStore->scfStore);
      for(i = 1; i <= maxI; i++)
      {
        FILE * fp;
        char fname[2048];
        cd = GetScaffoldCloneData(asmStore, i, TRUE, FALSE);
        getASM_SCFStore(asmStore->scfStore, i, &scf);
        sprintf(fname, F_UID ".mmps.txt", scf.uid);
        fp = fopen(fname, "w");
        PrintCloneData(scf.uid, cd, flags, fp);
        fclose(fp);
        DeleteCloneData(cd);
      }
      
      if(doDegenerates)
      {
        maxI = getNumASM_DSCs(asmStore->dscStore);
        for(i = 1; i <= maxI; i++)
        {
          cd = GetScaffoldCloneData(asmStore, i, TRUE, TRUE);
          getASM_DSCStore(asmStore->dscStore, i, &dsc);
          PrintCloneData(dsc.uid, cd, flags, stdout);
          DeleteCloneData(cd);
        }
      }
    }
    else
    {
      for(i = 0; i < GetNumVA_CDS_UID_t(uids); i++)
      {
        PHashValue_AS value;
        CDS_UID_t * uidp = GetVA_CDS_UID_t(uids, i);
        if(HASH_FAILURE == LookupInPHashTable_AS(asmStore->hashTable,
                                                 ASM_UID_NAMESPACE,
                                                 *uidp,
                                                 &value))
        {
          fprintf(stderr, "Failed to lookup " F_UID " in hashtable\n", *uidp);
          assert(0);
        }
        if(value.type == AS_IID_SCF)
        {
          cd = GetScaffoldCloneData(asmStore, value.IID, TRUE, FALSE);
          getASM_SCFStore(asmStore->scfStore, value.IID, &scf);
          PrintCloneData(scf.uid, cd, flags, stdout);
          DeleteCloneData(cd);
        }
        else
        {
          cd = GetScaffoldCloneData(asmStore, value.IID, TRUE, TRUE);
          getASM_DSCStore(asmStore->dscStore, value.IID, &dsc);
          PrintCloneData(dsc.uid, cd, flags, stdout);
          DeleteCloneData(cd);
        }
      }
    }
    CloseAssemblyStore(asmStore);
  }
  
  return 0;
}
