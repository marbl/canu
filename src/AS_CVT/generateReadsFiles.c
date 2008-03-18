
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

#include "AS_global.h"
#include "ASMData.h"

int main(int argc, char ** argv)
{
  char * inputStoreName = NULL;
  int doInstances = FALSE;
  int doSingleSurrogates = FALSE;
  int doDegenerates = FALSE;
  int doChaff = FALSE;
  int doUnreferenced = FALSE;
  {
    int ch,errflg=FALSE;
    while (!errflg && ((ch = getopt(argc, argv, "s:iSdclu")) != EOF))
    {
      switch(ch) 
      {
	case 's':
          inputStoreName = optarg;
          break;
        case 'i':
          doInstances = TRUE;
          break;
        case 'S':
          doSingleSurrogates = TRUE;
          break;
        case 'd':
          doDegenerates = TRUE;
          break;
        case 'c':
          doChaff = TRUE;
          break;
        case 'u':
          doUnreferenced = TRUE;
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
    fprintf(stderr, "Usage: %s [-s inputStore] [-iSdclu]\n"
            "\t-i      print instances - fragments in surrogate unitigs\n"
            "\t          that were multiply placed\n"
            "\t-S      print fragments in uniquely placed surrogates\n"
            "\t-d      print fragments in degenerate scaffolds\n"
            "\t-c      print chaff fragments\n"
            "\t-u      print unreferenced fragments that are in scaffolds\n\n",
            argv[0]);
    fprintf(stderr,
            "NOTE: This program generates a file similar to the standard\n"
            "      reads.placed file. Coordinates are 1-offset\n"
            "      Stanard reads.placed fields are as follows:\n"
            "\t(a) NCBI ti number for read (or *, if none known)\n"
            "\t(b) read name\n"
            "\t(c) start of trimmed read on original read\n"
            "\t(d) number of bases in trimmed read\n"
            "\t(e) orientation on contig (0 = forward, 1 = reverse)\n"
            "\t(f) contig name\n"
            "\t(g) supercontig name\n"
            "\t(h) approximate start of trimmed read on contig\n"
            "\t(i) approximate start of trimmed read on supercontig.\n\n");
    fprintf(stderr,
            "This program fudges (a) and (b) a as follows:\n"
            "(a) read UID\n"
            "(b) concatenation of the following:\n"
            "   fragment type (one character),\n"
            "   parent type (C = chaff, D = degenerate, S = inSurrogate,\n"
            "                U = unreferenced, R = none of the above), and\n"
            "   read UID\n\n");
    fprintf(stderr,
            "If run with -i, fragments in multiply placed surrogates will\n"
            "be listed multiple times, once for each placement of the surrogate\n\n");
    exit(1);
  }

  {
    AssemblyStore * asmStore;

    asmStore = OpenReadOnlyAssemblyStore(inputStoreName);
    PrintReadsPlaced(asmStore,
                     doInstances,
                     doSingleSurrogates,
                     doDegenerates,
                     doChaff,
                     doUnreferenced,
                     stdout);
    CloseAssemblyStore(asmStore);
  }
  
  return 0;
}
