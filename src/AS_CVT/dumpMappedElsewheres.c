
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
  char * mapStorePath = NULL;
  
  // parse command line
  {
    int ch,errflg=FALSE;
    while (!errflg && ((ch = getopt(argc, argv, "s:m:")) != EOF))
    {
      switch(ch) 
      {
        case 's':
          asmStorePath = optarg;
          break;
        case 'm':
          mapStorePath = optarg;
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
            "Usage: %s [-s assemblyStore] [-m mapStore]\n",
            argv[0]);
    exit(-1);
  }

  {
    AssemblyStore * asmStore;
    MapStore * mapStore;
    int32 i;
    
    asmStore = OpenReadOnlyAssemblyStore(asmStorePath);
    mapStore = OpenReadOnlyMapStore(mapStorePath);

    for(i = 1; i <= getNumASM_CHRs(mapStore->chrStore); i++)
    {
      ASM_CHRRecord chr;
      char fname[1024];
      FILE * fo;
      getASM_CHRStore(mapStore->chrStore, i, &chr);
      
      sprintf(fname, "%03" F_UIDP ".txt", chr.uid);
      fo = fopen(fname, "w");
      if(fo == NULL)
      {
        fprintf(stderr, "Failed to open file %s for writing.\n", fname);
        exit(1);
      }
      fprintf(stderr, "Writing chromosome UID " F_UID " elsewhere mapped mate pairs to %s\n",
              chr.uid, fname);
      PrintChromosomeElsewheres(asmStore, mapStore, chr.uid, fo);
      fclose(fo);
    }

    CloseAssemblyStore(asmStore);
    CloseMapStore(mapStore);
  }
  return 0;
}
