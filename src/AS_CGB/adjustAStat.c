
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

const char *mainid = "$Id: adjustAStat.c,v 1.4 2008-10-08 22:02:54 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_UTL_Hash.h"
#include "AS_UTL_fileIO.h"
#include "AS_MSG_pmesg.h"

FILE            *errorFP   = NULL;

static
void
usage(char *filename, int longhelp) {
fprintf(stdout, "In file %s, i am about to print the usage\n", filename);
  fprintf(stdout, "usage: %s -u <unitig list> -c <cgb file>\n", filename);
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "The file will edit the supplied CGB file using the list of unitig IDs and a-stats in unitig list\n");
  fprintf(stdout, "\n");
}

int
main(int argc, char **argv) {

  //  Options for everyone.
  //
  char            *inputFileName    = NULL;
  char            *cgbFileName      = NULL;
  char            *errorFile        = NULL;

  inputFileName   = NULL;
  cgbFileName     = NULL;
  errorFile       = NULL;
  errorFP         = stderr;

  int arg = 1;
  int err = 0;
  int hlp = 0;

  argc = AS_configure(argc, argv);

  while (arg < argc) {
    if        (strcmp(argv[arg], "-u") == 0) {
      inputFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-c") == 0) {
      cgbFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-E") == 0) {
      errorFile = argv[++arg];
    } else if (strcmp(argv[arg], "-h") == 0) {
      hlp++;
      err++;
    } else {
      err++;
    }
    arg++;
  }

  if ((err) || (inputFileName == NULL) || (cgbFileName == NULL)) {
    usage(argv[0], hlp);
    exit(1);
  }

   if (errorFile) {
      errno = 0;
      errorFP = fopen(errorFile, "w");
      if (errno) {
         fprintf(stderr, "%s: Couldn't open -E file '%s': %s\n", argv[0], errorFile, strerror(errno));
         exit(1);
      }
   }

   //  Construct an IID map to a-stat of objects we are supplied to be used for updating
   //
   HashTable_AS* IIDtoASTAT = CreateScalarHashTable_AS(32 * 1024);
   FILE            *F        = NULL;
   char             L[1024];

   errno = 0;
   F = fopen(inputFileName, "r");
   if (errno) {
      fprintf(errorFP, "%s: Couldn't open -u file '%s': %s\n", argv[0], inputFileName, strerror(errno));
      exit(1);
   }

   fgets(L, 1024, F);
   while (!feof(F)) {
      char *endPtr;
      AS_IID    iid = AS_IID_fromString(L, &endPtr);

      char in[2];
      char ch;
      sscanf(endPtr, " %1[XUR]", in);
      ch = in[0];

      UnitigFUR status = (UnitigFUR) ch;

      InsertInHashTable_AS(IIDtoASTAT, iid, 0, (uint64)status, 0);
      fgets(L, 1024, F);
   }
   fclose(F);

   // now read in the cgb file and modify the appropriate records astat values
   GenericMesg  *pmesg;
   errno = 0;
   F = fopen(cgbFileName, "r");
   if (errno) {
      fprintf(errorFP, "%s: Couldn't open -c file '%s': %s\n", argv[0], cgbFileName, strerror(errno));
      exit(1);
   }

   while( EOF != ReadProtoMesg_AS(F, &pmesg)) {
      const MessageType imesgtype = pmesg->t;

      // we only care about the IUM messages which contain the cov: stat
      if (pmesg->t == MESG_IUM) {
         IntUnitigMesg *o = (IntUnitigMesg *)pmesg->m;

      if (ExistsInHashTable_AS(IIDtoASTAT, o->iaccession, 0)) {
         o->unique_rept = (UnitigFUR)(LookupValueInHashTable_AS(IIDtoASTAT, o->iaccession, 0));
      }
    }
    WriteProtoMesg_AS(stdout, pmesg);
  }
  fclose(F);

  if (errorFile)
    fclose(errorFP);

  return(0);
}
