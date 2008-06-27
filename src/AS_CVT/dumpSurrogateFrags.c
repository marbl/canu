
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

int main(int argc, char ** argv)
{
  char * inputStoreName = NULL;
  {
    int ch,errflg=FALSE;
    while (!errflg && ((ch = getopt(argc, argv, "s:")) != EOF))
    {
      switch(ch)
      {
	case 's':
          inputStoreName = optarg;
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
    fprintf(stderr, "Usage: %s [-s inputStore]\n", argv[0]);
    exit(1);
  }

  {
    AssemblyStore * asmStore;
    int32 i;
    int32 maxI;
    ASM_AFGRecord afg;
    ASM_UTGRecord utg;

    asmStore = OpenReadOnlyAssemblyStore(inputStoreName);
    maxI = getNumASM_AFGs(asmStore->afgStore);
    for(i = 1; i <= maxI; i++)
    {
      getASM_AFGStore(asmStore->afgStore, i, &afg);
      if(!afg.chaff)
      {
        getASM_UTGStore(asmStore->utgStore, afg.unitigIndex, &utg);
        if(utg.status == AS_SEP)
          fprintf(stdout, F_UID "\t" F_UID "\n", afg.uid, utg.uid);
      }
    }
    CloseAssemblyStore(asmStore);
  }

  return 0;
}
