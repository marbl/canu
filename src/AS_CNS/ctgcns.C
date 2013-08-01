
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

const char *mainid = "$Id: ctgcns.C,v 1.12 2012-01-05 18:29:43 brianwalenz Exp $";

#include "AS_global.H"
#include "MultiAlign.H"
#include "MultiAlignStore.H"
#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"

int
main (int argc, char **argv) {
  char   tmpName[FILENAME_MAX] = {0};

  char  *gkpName = NULL;

  char  *tigName = NULL;
  int32  tigVers = -1;
  int32  tigPart = -1;

  int32  ctgTest = -1;
  char  *ctgFile = NULL;

  bool   forceCompute = false;

  int32  numFailures = 0;
  int32  numSkipped  = 0;

  bool   showResult = false;

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

    } else if (strcmp(argv[arg], "-c") == 0) {
      ctgTest = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-T") == 0) {
      ctgFile = argv[++arg];

    } else if (strcmp(argv[arg], "-f") == 0) {
      forceCompute = true;

    } else if (strcmp(argv[arg], "-v") == 0) {
      showResult = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE_MULTIALIGN_OUTPUT++;

    } else if (strcmp(argv[arg], "-w") == 0) {
      options.smooth_win = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-P") == 0) {
      options.do_phasing = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == NULL) || (tigName == NULL)) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version partition [opts]\n", argv[0]);
    fprintf(stderr, "    -c id        Compute only contig 'id' (must be in the correct partition!)\n");
    fprintf(stderr, "    -T file      Test the computation of the contig layout in 'file'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -f           Recompute contigs that already have a multialignment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -v           Show multialigns.\n");
    fprintf(stderr, "    -V           Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -w ws        Smoothing window size\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  //  Open both stores for read only.
  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  tigStore = new MultiAlignStore(tigName, tigVers, 0, tigPart, FALSE, FALSE, FALSE);

  gkpStore->gkStore_loadPartition(tigPart);

  //  Decide on what to compute.  Either all contigs, or a single contig, or a special case test.
  uint32 b = 0;
  uint32 e = tigStore->numContigs();

  if (ctgTest != -1) {
    b = ctgTest;
    e = ctgTest + 1;
  }

  FORCE_UNITIG_ABUT = 1;

  if (ctgFile != NULL) {
    errno = 0;
    FILE         *F = fopen(ctgFile, "r");
    if (errno)
      fprintf(stderr, "Failed to open input contig file '%s': %s\n", ctgFile, strerror(errno)), exit(1);

    MultiAlignT  *ma       = CreateEmptyMultiAlignT();
    bool          isUnitig = false;

    while (LoadMultiAlignFromHuman(ma, isUnitig, F) == true) {
      if (ma->maID < 0)
        ma->maID = (isUnitig) ? tigStore->numUnitigs() : tigStore->numContigs();

      if (MultiAlignContig(ma, gkpStore, &options)) {
        if (showResult)
          PrintMultiAlignT(stdout, ma, gkpStore, false, false, AS_READ_CLEAR_LATEST);
      } else {
        fprintf(stderr, "MultiAlignContig()-- contig %d failed.\n", ma->maID);
        numFailures++;
      }
    }

    DeleteMultiAlignT(ma);

    b = e = 0;
  }

  //  Reopen for writing, if we have work to do.
  if (b < e) {
    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, 0, tigPart, TRUE, FALSE, TRUE);
  }

  //  Now the usual case.  Iterate over all unitigs, compute and update.
  for (uint32 i=b; i<e; i++) {
    MultiAlignT  *ma = tigStore->loadMultiAlign(i, FALSE);

    if (ma == NULL) {
      //  Not in our partition, or deleted.
      continue;
    }

    bool  exists = (ma->consensus != NULL) && (GetNumchars(ma->consensus) > 1);

    if ((forceCompute == false) && (exists == true)) {
      //  Already finished contig consensus.
      fprintf(stderr, "Working on contig %d (%d unitigs and %d fragments) - already computed, skipped\n",
              ma->maID, ma->data.num_unitigs, ma->data.num_frags);
      numSkipped++;
      continue;
    }

    fprintf(stderr, "Working on contig %d (%d unitigs and %d fragments)%s\n",
            ma->maID, ma->data.num_unitigs, ma->data.num_frags,
            (exists) ? " - already computed, recomputing" : "");

    if (MultiAlignContig(ma, gkpStore, &options)) {
      tigStore->insertMultiAlign(ma, FALSE, FALSE);
      if (showResult)
        PrintMultiAlignT(stdout, ma, gkpStore, false, false, AS_READ_CLEAR_LATEST);
      DeleteMultiAlignT(ma);
    } else {
      fprintf(stderr, "MultiAlignContig()-- contig %d failed.\n", ma->maID);
      numFailures++;
    }
  }

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
    fprintf(stderr, "WARNING:  Total number of contig failures = %d\n", numFailures);
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus did NOT finish successfully.\n");
    return(1);
  }

  fprintf(stderr, "Consensus finished successfully.  Bye.\n");
  return(0);
}
