
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

const char *mainid = "$Id: utgcns.C,v 1.12 2011-11-29 11:50:00 brianwalenz Exp $";

#include "AS_global.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"


int
main (int argc, char **argv) {
  char  *gkpName = NULL;

  char  *tigName = NULL;
  int32  tigVers = -1;
  int32  tigPart = -1;

  int32  utgTest = -1;
  char  *utgFile = NULL;

  bool   forceCompute = false;

  int32  numFailures = 0;
  int32  numSkipped  = 0;

  CNS_PrintKey printwhat=CNS_STATS_ONLY;

  CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                          CNS_OPTIONS_MIN_ANCHOR_DEFAULT,
                          CNS_OPTIONS_DO_PHASING_DEFAULT };

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

      if (tigVers <= 0)
        fprintf(stderr, "invalid tigStore version (-t store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);
      if ((tigPart <= 0) && (argv[arg][0] != '.'))
        fprintf(stderr, "invalid tigStore partition (-t store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);

    } else if (strcmp(argv[arg], "-u") == 0) {
      utgTest = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-T") == 0) {
      utgFile = argv[++arg];

    } else if (strcmp(argv[arg], "-f") == 0) {
      forceCompute = true;

    } else if (strcmp(argv[arg], "-v") == 0) {
      printwhat = CNS_VIEW_UNITIG;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE_MULTIALIGN_OUTPUT++;

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == NULL) || (tigName == NULL)) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version partition [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    -u id        Compute only unitig 'id' (must be in the correct partition!)\n");
    fprintf(stderr, "    -T file      Test the computation of the unitig layout in 'file'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -f           Recompute unitigs that already have a multialignment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -v           Show multialigns.\n");
    fprintf(stderr, "    -V           Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  //  Open both stores for read only.
  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  tigStore = new MultiAlignStore(tigName, tigVers, tigPart, 0, FALSE, FALSE, FALSE);

  gkpStore->gkStore_loadPartition(tigPart);

  //  Decide on what to compute.  Either all unitigs, or a single unitig, or a special case test.

  uint32  b = 0;
  uint32  e = tigStore->numUnitigs();

  if (utgTest != -1) {
    b = utgTest;
    e = utgTest + 1;
  }

  if (utgFile != NULL) {
    errno = 0;
    FILE         *F = fopen(utgFile, "r");
    if (errno)
      fprintf(stderr, "Failed to open input unitig file '%s': %s\n", utgFile, strerror(errno)), exit(1);

    MultiAlignT  *ma       = CreateEmptyMultiAlignT();
    bool          isUnitig = false;

    while (LoadMultiAlignFromHuman(ma, isUnitig, F) == true) {
      if (ma->maID < 0)
        ma->maID = (isUnitig) ? tigStore->numUnitigs() : tigStore->numContigs();

      int32 firstFailed = 0;

      if (MultiAlignUnitig(ma, gkpStore, printwhat, &options, firstFailed)) {
      } else {
        fprintf(stderr, "MultiAlignUnitig()-- unitig %d failed at fragment index %d.\n",
                ma->maID, firstFailed);
        numFailures++;
      }
    }

    DeleteMultiAlignT(ma);

    b = e = 0;
  }

  //  Reopen for writing, if we have work to do.
  if (b < e) {
    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, tigPart, 0, TRUE, FALSE, TRUE);
  }

  fprintf(stderr, "Computing unitig consensus for b="F_U32" to e="F_U32"\n", b, e);

  //  Now the usual case.  Iterate over all unitigs, compute and update.
  for (uint32 i=b; i<e; i++) {
    MultiAlignT  *ma = tigStore->loadMultiAlign(i, TRUE);

    if (ma == NULL) {
      //  Not in our partition, or deleted.
      continue;
    }

    bool exists = (ma->consensus != NULL) && (GetNumchars(ma->consensus) > 1);

    if ((forceCompute == false) && (exists == true)) {
      //  Already finished unitig consensus.
      if (ma->data.num_frags > 1)
        fprintf(stderr, "Working on unitig %d (%d unitigs and %d fragments) - already computed, skipped\n",
                ma->maID, ma->data.num_unitigs, ma->data.num_frags);
      numSkipped++;
      continue;
    }

    if (ma->data.num_frags > 1)
      fprintf(stderr, "Working on unitig %d (%d unitigs and %d fragments)%s\n",
              ma->maID, ma->data.num_unitigs, ma->data.num_frags,
              (exists) ? " - already computed, recomputing" : "");

    int32 firstFailed = 0;

    if (MultiAlignUnitig(ma, gkpStore, printwhat, &options, firstFailed)) {
      tigStore->insertMultiAlign(ma, TRUE, FALSE);
      DeleteMultiAlignT(ma);
    } else {
      fprintf(stderr, "MultiAlignUnitig()-- unitig %d failed at fragment index %d.\n",
              ma->maID, firstFailed);
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
