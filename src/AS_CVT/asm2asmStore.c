
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
  char * asmFilename = NULL;
  char * storePath = NULL;
  char * gkpStorePath = NULL;
  char * frgStorePath = NULL;
  AssemblyStore * asmStore;

  // parse command line
  {
    int ch,errflg=FALSE;
    while (!errflg && ((ch = getopt(argc, argv, "a:s:g:f:")) != EOF))
    {
      switch(ch) 
      {
	case 'a':
          asmFilename = optarg;
          break;
        case 's':
          storePath = optarg;
          break;
        case 'g':
          gkpStorePath = optarg;
          break;
        case 'f':
          frgStorePath = optarg;
          break;
        default:
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	  errflg++;
          break;
      }
    }
  }

  // check command line requirements
  if(asmFilename == NULL || storePath == NULL)
  {
    fprintf(stderr,
            "Usage: %s [-a asmFilename] [-s storePath] [-g gkpStore] [-f frgStore]\n"
            "\t-g & -f are optional, but strongly encouraged.\n",
            argv[0]);
    exit(-1);
  }

  // read asm file into data structures
  {
    FILE * fi = fopen(asmFilename, "r");

    assert(fi != NULL);
    
    asmStore = CreateAssemblyStoreFromASMFile(fi, storePath,
                                              gkpStorePath, frgStorePath);
    assert(asmStore != NULL);
  }

  CloseAssemblyStore(asmStore);
  return 0;
}
