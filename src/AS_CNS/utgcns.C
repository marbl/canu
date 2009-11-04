
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

const char *mainid = "$Id: utgcns.C,v 1.6 2009-11-04 17:21:08 brianwalenz Exp $";

#include "AS_global.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"

int
main (int argc, char **argv) {
  char   tmpName[FILENAME_MAX] = {0};

  char  *gkpName = NULL;

  char  *tigName = NULL;
  int32  tigVers = -1;
  int32  tigPart = -1;

  int32  utgTest = -1;

  int32  numFailures = 0;

  CNS_PrintKey printwhat=CNS_STATS_ONLY;

  CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                          CNS_OPTIONS_MIN_ANCHOR_DEFAULT };

  //  Comminucate to MultiAlignment_CNS.c that we are doing consensus and not cgw.
  thisIsConsensus = 1;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);
      tigPart = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-u") == 0) {
      utgTest = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-v") == 0) {
      printwhat = CNS_VIEW_UNITIG;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE_MULTIALIGN_OUTPUT = 1;

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == NULL) || (tigName == NULL)) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version partition [opts]\n", argv[0]);
    fprintf(stderr, "    -u id        Compute only unitig 'id' (must be in the correct partition!)\n");
    fprintf(stderr, "    -v           Show multialigns.\n");
    fprintf(stderr, "    -V           Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  tigStore = new MultiAlignStore(tigName, tigVers, tigPart, 0, TRUE, FALSE);

  gkpStore->gkStore_loadPartition(tigPart);

  uint32  b = 0;
  uint32  e = tigStore->numUnitigs();

  if (utgTest != -1) {
    b = utgTest;
    e = utgTest + 1;
  }

  for (uint32 i=b; i<e; i++) {
    MultiAlignT  *ma = tigStore->loadMultiAlign(i, TRUE);

    if (ma == NULL)
      //  Not in our partition, or deleted.
      continue;

    if (ma->data.num_frags > 1)
      fprintf(stderr, "Working on unitig %d (%d unitigs and %d fragments)\n",
              ma->maID, ma->data.num_unitigs, ma->data.num_frags);

    if (MultiAlignUnitig(ma, gkpStore, printwhat, &options)) {
      tigStore->insertMultiAlign(ma, TRUE, FALSE);
      DeleteMultiAlignT(ma);
    } else {
      fprintf(stderr, "MultiAlignUnitig()-- unitig %d failed.\n", i);
      numFailures++;
    }
  }

 finish:
  delete tigStore;

  fprintf(stderr, "\n");
  fprintf(stderr, "NumColumnsInUnitigs             = %d\n", NumColumnsInUnitigs);
  fprintf(stderr, "NumGapsInUnitigs                = %d\n", NumGapsInUnitigs);
  fprintf(stderr, "NumRunsOfGapsInUnitigReads      = %d\n", NumRunsOfGapsInUnitigReads);
  fprintf(stderr, "NumColumnsInContigs             = %d\n", NumColumnsInContigs);
  fprintf(stderr, "NumGapsInContigs                = %d\n", NumGapsInContigs);
  fprintf(stderr, "NumRunsOfGapsInContigReads      = %d\n", NumRunsOfGapsInContigReads);
  fprintf(stderr, "NumAAMismatches                 = %d\n", NumAAMismatches);
  fprintf(stderr, "NumVARRecords                   = %d\n", NumVARRecords);
  fprintf(stderr, "NumVARStringsWithFlankingGaps   = %d\n", NumVARStringsWithFlankingGaps);
  fprintf(stderr, "NumUnitigRetrySuccess           = %d\n", NumUnitigRetrySuccess);
  fprintf(stderr, "\n");

  if (numFailures) {
    fprintf(stderr, "WARNING:  Total number of unitig failures = %d\n", numFailures);
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus did NOT finish successfully.\n");
    return(1);
  }

  fprintf(stderr, "Consensus finished successfully.  Bye.\n");
  return(0);
}
