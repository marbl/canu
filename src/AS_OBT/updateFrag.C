/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

const char *mainid = "$Id: updateFrag.C,v 1.4 2008-10-08 22:02:57 brianwalenz Exp $";

#include "constants.H"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
}

#include "util++.H"


uint32         lineMax = 128 * 1024;
char          *line    = 0L;

void
readLine(FILE *F) {
  fgets(line, lineMax, F);
  line[lineMax-1] = 0;
  assert(strlen(line) < (lineMax-1));
}

int
main(int argc, char **argv) {
  FILE    *iidFile           = 0L;
  char    *frgStore          = 0L;
  bool     doModify          = true;  //  Make this false for testing

  line = new char [lineMax];

  argc = AS_configure(argc, argv);

  if (argc < 4) {
    fprintf(stderr, "usage: %s [-iid iid] -frg frgStore\n", argv[0]);
    fprintf(stderr, "  -iid x   The iids of fragment that we need to update\n");
    fprintf(stderr, "  -frg f   'f' is our frag store\n");
    exit(1);
  }

  int arg=1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-frg", 2) == 0) {
      frgStore = argv[++arg];
    } else if (strncmp(argv[arg], "-iid", 2) == 0) {
      errno=0;
      iidFile = fopen(argv[++arg], "r");
      if (errno)
        fprintf(stderr, "Failed to open %s for reading the iid list: %s\n", argv[arg], strerror(errno)), exit(1);
      }
    arg++;
  }

  //  Open the frgStore, prepare for reading fragments
  //
  GateKeeperStore *gkp = openGateKeeperStore(frgStore, doModify);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open fragStore %s!\n", frgStore);
    exit(1);
  }

  gkpStore->frg = convertStoreToMemoryStore(gkpStore->frg);

  fragRecord fr;

  uint64  iid   = 0;
  uint32  left  = 0;
  uint32  right = 0;

  readLine(iidFile);
  while (!feof(iidFile)) {
    splitToWords  W(line);
    iid    = atoi(W[0]);

    //  Read the fragment from the store, compute the adjustment
    //  points.  All the values from the mapping are off by the original
    //  clear range, we add it back in as we decode the string.

    getFrag(gkp, iid, &fr, FRAG_S_INF | FRAG_S_QLT);

    uint32 qltLQ1 = getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_OBT);
    uint32 qltRQ1 = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_OBT);

    GateKeeperLibraryRecord  *gklr = getGateKeeperLibrary(gkp, getFragRecordLibraryIID(fr));

    //  Only proceed if we're mutable.
    //
    if ((gklr) && (gklr->doNotOverlapTrim)) {
    }
    else {
      left  = atoi(W[1]) + qltLQ1;
      right = atoi(W[2]) + qltLQ1;

      if (doModify) {
        setFragRecordClearRegion(fr, left, right, AS_READ_CLEAR_OBT);
        setFrag(gkp, iid, fr);
      }
    }
    readLine(iidFile);
  }

  closeGateKeeperStore(gkp);
  fclose(iidFile);

  return(0);
}
