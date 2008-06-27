
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
  char * deflineFilename = NULL;
  char * assemblyName = NULL;
  {
    int ch,errflg=FALSE;
    while (!errflg && ((ch = getopt(argc, argv, "s:d:n:")) != EOF))
    {
      switch(ch)
      {
	case 's':
          inputStoreName = optarg;
          break;
        case 'd':
          deflineFilename = optarg;
          break;
        case 'n':
          assemblyName = optarg;
          break;
        default:
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	  errflg++;
          break;
      }
    }
  }

  if(inputStoreName == NULL || deflineFilename == NULL || assemblyName == NULL)
  {
    fprintf(stderr, "Usage: %s [-s asmStore] [-d deflineFile] [-n assemblyName]\n", argv[0]);
    exit(1);
  }

  {
    AssemblyStore * asmStore;
    FILE * dfp;

    asmStore = OpenReadOnlyAssemblyStore(inputStoreName);

    InitializeATACFile(asmStore, stdout);

    dfp = fopen(deflineFilename, "r");
    PrintDeflineATACAxes(asmStore, dfp, assemblyName, stdout);
    fclose(dfp);

    PrintATACSurrogates(asmStore, assemblyName, stdout);

    CloseAssemblyStore(asmStore);
  }
  return 0;
}
