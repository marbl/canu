
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

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_fileIO.h"

void
usage(char *name) {
  fprintf(stderr, "usage: %s [-x] [-i] [-m message type] [-o outputfile] < <input file>\n", name);
  fprintf(stderr, "       -i      include the following messages in the next output\n");
  fprintf(stderr, "       -x      exclude the following messages from the next output\n");
  fprintf(stderr, "       -m      message\n");
  fprintf(stderr, "       -o      write output here\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "extractmessages attempts to construct a filter chain to put every message\n");
  fprintf(stderr, "into a specific file.  Using the -i and -x switches, you can specify messages\n");
  fprintf(stderr, "to include in the next file or to exclude from the next file.\n");
  fprintf(stderr, "For example:\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  extractmessages -i -m ICM -m IDS -o icm-and-ids -x -m IAF -o everythingelse > /dev/null\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "includes only ICM and IDS messages in the first file, then excludes IAF messages\n");
  fprintf(stderr, "from the second file, and everything else (here, just IAF messages) are written\n");
  fprintf(stderr, "to stdout.\n");
  fprintf(stderr, "\n");
}

  
int
main(int argc, char **argv) {
  int            msglist[NUM_OF_REC_TYPES + 1];
  FILE          *outfile[NUM_OF_REC_TYPES + 1];
  off_t          count[NUM_OF_REC_TYPES + 1];
  off_t          size[NUM_OF_REC_TYPES + 1];
  int            i;

  for (i=0; i<=NUM_OF_REC_TYPES; i++) {
    msglist[i] = 0;
    outfile[i] = 0L;
    count[i]   = 0;
    size[i]    = 0;
  }

  int arg = 1;
  int inc = 0;
  int err = 0;
  int msg = 0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-i") == 0) {
      inc = 1;
    } else if (strcmp(argv[arg], "-x") == 0) {
      inc = 0;
    } else if (strcmp(argv[arg], "-o") == 0) {
      errno = 0;
      FILE *F = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "%s: failed to open output file '%s': %s\n", argv[0], argv[arg], strerror(errno)), exit(1);

      //  Depending on the include flag, we either write all messages
      //  listed in our msglist (or write all message not in the
      //  msglist) to the freshly opened file.
      //
      if (inc) {
        //  Include message i in the output if it was listed
        for (i=1; i<=NUM_OF_REC_TYPES; i++)
          if ((outfile[i] == NULL) && (msglist[i] > 0))
            outfile[i] = F;
      } else {
        //  Include message i in the output if it was not listed
        for (i=1; i<=NUM_OF_REC_TYPES; i++)
          if ((outfile[i] == NULL) && (msglist[i] == 0))
            outfile[i] = F;
      }

      for (i=0; i<=NUM_OF_REC_TYPES; i++)
        msglist[i] = 0;
    } else if (strcmp(argv[arg], "-m") == 0) {
      int type = GetMessageType(argv[++arg]);
      if ((type >= 1) && (type <= NUM_OF_REC_TYPES)) {
        msglist[type]++;
        msg++;
      } else {
        fprintf(stderr, "%s: invalid message type '%s'.\n", argv[0], argv[arg]);
        err = 1;
      }
    } else if (strcmp(argv[arg], "-h") == 0) {
      err = 1;
    } else {
      int type = GetMessageType(argv[arg]);
      if ((type >= 1) && (type <= NUM_OF_REC_TYPES)) {
        msglist[type]++;
        msg++;
      } else {
        fprintf(stderr, "%s: invalid option '%s'.\n", argv[0], argv[arg]);
        err = 1;
      }
    }
    arg++;
  }

  if (err)
    usage(argv[0]), exit(1);

  //  Assume everything else goes to stdout.  We need to obey the inc
  //  flag, still, though.
  //
  if (inc) {
    //  Include message i in the output if it was listed
    for (i=1; i<=NUM_OF_REC_TYPES; i++)
      if ((outfile[i] == NULL) && (msglist[i] > 0))
        outfile[i] = stdout;
  } else {
    //  Include message i in the output if it was not listed
    for (i=1; i<=NUM_OF_REC_TYPES; i++)
      if ((outfile[i] == NULL) && (msglist[i] == 0))
        outfile[i] = stdout;
  }

  GenericMesg   *pmesg;
  off_t          currPos = 0;
  off_t          prevPos = 0;

  while (ReadProtoMesg_AS(stdin, &pmesg) != EOF) {
    assert(pmesg->t <= NUM_OF_REC_TYPES);

    currPos = AS_UTL_ftell(stdin);

    if (outfile[pmesg->t] != NULL) {
      count[pmesg->t]++;

      size[pmesg->t] += currPos - prevPos;

      WriteProtoMesg_AS(outfile[pmesg->t], pmesg);
    }

    prevPos = currPos;
  }

  for (i=0; i<=NUM_OF_REC_TYPES; i++)
    if (count[i] > 0)
      fprintf(stderr, "%s num "F_OFF_T" size "F_OFF_T" avg %f\n",
              MessageTypeName[i], count[i], size[i], (double)size[i] / count[i]);

  exit(0);
}

