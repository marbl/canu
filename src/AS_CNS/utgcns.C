
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

const char *mainid = "$Id: utgcns.C,v 1.16 2012-08-29 06:01:22 brianwalenz Exp $";

#include "AS_global.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"

#include "AS_UTL_decodeRange.H"

//  Create a new f_list for the ma that has no contained reads.
//  The original f_list is returned.
//
VA_TYPE(IntMultiPos) *
stashContains(MultiAlignT *ma, double threshold) {
  VA_TYPE(IntMultiPos) *fl = ma->f_list;

  int32  nOrig = GetNumIntMultiPoss(fl);
  int32  nRedu = 0;

  for (uint32 fi=0; fi<GetNumIntMultiPoss(fl); fi++) {
    IntMultiPos  *imp = GetIntMultiPos(fl, fi);

    if (imp->contained == 0)
      nRedu++;
  }

  if (1.0 - (double)nRedu / nOrig < threshold)
    return(NULL);

  fprintf(stderr, "Removed "F_S32" contains from unitig; processing only "F_S32" reads.\n",
          nOrig - nRedu, nRedu);

  ma->f_list = CreateVA_IntMultiPos(0);

  for (uint32 fi=0; fi<GetNumIntMultiPoss(fl); fi++) {
    IntMultiPos  *imp = GetIntMultiPos(fl, fi);

    if (imp->contained == 0)
      AppendVA_IntMultiPos(ma->f_list, imp);
  }

  return(fl);
}


//  Restores the f_list, and updates the position of non-contained reads.
//
void
unstashContains(MultiAlignT *ma, VA_TYPE(IntMultiPos) *fl) {

  if (fl == NULL)
    return;

  uint32   oldMax = 0;
  uint32   newMax = 0;

  //  For fragments not involved in the consensus computation, we'll scale their position linearly
  //  from the old max to the new max.
  //
  //  We probably should do an alignment to the consensus sequence to find the true location, but
  //  that's (a) expensive and (b) likely overkill for these unitigs.
  //
  for (uint32 fi=0, ci=0; fi<GetNumIntMultiPoss(fl); fi++) {
    IntMultiPos  *imp = GetIntMultiPos(fl, fi);

    if (oldMax < imp->position.bgn)    oldMax = imp->position.bgn;
    if (oldMax < imp->position.end)    oldMax = imp->position.end;
  }

  newMax = GetMultiAlignLength(ma);

  double sf = (double)newMax / oldMax;

  for (uint32 fi=0, ci=0; fi<GetNumIntMultiPoss(fl); fi++) {
    IntMultiPos  *imp = GetIntMultiPos(fl, fi);

    if (imp->contained == 0) {
      //  Copy the location used by consensus back to the original list
      SetVA_IntMultiPos(fl, fi, GetVA_IntMultiPos(ma->f_list, ci));
      ci++;

    } else {
      //  Adjust old position
      imp->position.bgn = sf * imp->position.bgn;
      imp->position.end = sf * imp->position.end;

      if (imp->position.bgn > newMax)  imp->position.bgn = newMax;
      if (imp->position.end > newMax)  imp->position.end = newMax;
    }
  }

  DeleteVA_IntMultiPos(ma->f_list);

  ma->f_list = fl;
}



int
main (int argc, char **argv) {
  char  *gkpName = NULL;

  char  *tigName = NULL;
  int32  tigVers = -1;
  int32  tigPart = -1;

  uint64 utgBgn = -1;
  uint64 utgEnd = -1;
  char  *utgFile = NULL;

  bool   forceCompute = false;

  int32  numFailures = 0;
  int32  numSkipped  = 0;

  bool   showResult = false;

  bool   ignoreContains = false;
  double ignoreContainT = 0.75;

  bool   inplace = false;
  bool   loadall = false;

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
      AS_UTL_decodeRange(argv[++arg], utgBgn, utgEnd);

    } else if (strcmp(argv[arg], "-T") == 0) {
      utgFile = argv[++arg];

    } else if (strcmp(argv[arg], "-f") == 0) {
      forceCompute = true;

    } else if (strcmp(argv[arg], "-v") == 0) {
      showResult = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE_MULTIALIGN_OUTPUT++;

    } else if (strcmp(argv[arg], "-nocontains") == 0) {
      ignoreContains = true;
      ignoreContainT = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-inplace") == 0) {
      inplace = true;

    } else if (strcmp(argv[arg], "-loadall") == 0) {
      loadall = true;

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if (gkpName == NULL)
    err++;
  if ((utgFile == NULL) && (tigName == NULL))
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version partition [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    -u b         Compute only unitig ID 'b' (must be in the correct partition!)\n");
    fprintf(stderr, "    -u b-e       Compute only unitigs from ID 'b' to ID 'e'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -T file      Test the computation of the unitig layout in 'file'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -f           Recompute unitigs that already have a multialignment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -v           Show multialigns.\n");
    fprintf(stderr, "    -V           Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  ADVANCED OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -nocontains F   Do not use contains, if the unitig is more than fraction F\n");
    fprintf(stderr, "                    contained reads.  The consensus sequence will be computed from\n");
    fprintf(stderr, "                    only the non-contained reads.  The layout will still retain the\n");
    fprintf(stderr, "                    contained reads, but their location will not be precise.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -inplace        Write the updated unitig to the same version it was read from.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -t S V P        If 'partition' is '.', use an unpartitioned tigStore/gkpStore.\n");
    fprintf(stderr, "    -loadall        If not partitioned, load ALL reads into memory.\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR:  No gkpStore (-g) supplied.\n");

    if ((utgFile == NULL) && (tigName == NULL))
      fprintf(stderr, "ERROR:  No tigStore (-t) OR no test unitig (-T) supplied.\n");

    exit(1);
  }

  //  Open gatekeeper for read only.

  gkpStore = new gkStore(gkpName, FALSE, FALSE);

  //  If we are testing a unitig, do that.

  if (utgFile != NULL) {
    errno = 0;
    FILE         *F = fopen(utgFile, "r");
    if (errno)
      fprintf(stderr, "Failed to open input unitig file '%s': %s\n", utgFile, strerror(errno)), exit(1);

    MultiAlignT  *ma       = CreateEmptyMultiAlignT();
    bool          isUnitig = false;

    while (LoadMultiAlignFromHuman(ma, isUnitig, F) == true) {
      //if (ma->maID < 0)
      //  ma->maID = (isUnitig) ? tigStore->numUnitigs() : tigStore->numContigs();

      if (MultiAlignUnitig(ma, gkpStore, &options, NULL)) {
        if (showResult)
          PrintMultiAlignT(stdout, ma, gkpStore, false, false, AS_READ_CLEAR_LATEST);
      } else {
        fprintf(stderr, "MultiAlignUnitig()-- unitig %d failed.\n", ma->maID);
        numFailures++;
      }
    }

    DeleteMultiAlignT(ma);

    exit(0);
  }

  //  Otherwise, we're computing unitigs from the store.  Open it for read only,
  //  and load the reads.

  tigStore = new MultiAlignStore(tigName, tigVers, tigPart, 0, FALSE, FALSE, FALSE);

  if (tigPart == 0) {
    if (loadall) {
      fprintf(stderr, "Loading all reads into memory.\n");
      gkpStore->gkStore_load(0, 0, GKFRAGMENT_QLT);
    }
  } else {
    gkpStore->gkStore_loadPartition(tigPart);
  }

  //  Decide on what to compute.  Either all unitigs, or a single unitig, or a special case test.

  uint32  b = 0;
  uint32  e = tigStore->numUnitigs();

  if (utgBgn != -1) {
    b = utgBgn;
    e = utgEnd + 1;
  }

  //  Reopen for writing, if we have work to do.

  if (b < e) {
    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, tigPart, 0, TRUE, inplace, !inplace);
  }

  fprintf(stderr, "Computing unitig consensus for b="F_U32" to e="F_U32"\n", b, e);

  //  Now the usual case.  Iterate over all unitigs, compute and update.

  for (uint32 i=b; i<e; i++) {
    MultiAlignT              *ma = tigStore->loadMultiAlign(i, TRUE);
    VA_TYPE(IntMultiPos)     *fl = NULL;

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

    //  Build a new ma if we're ignoring contains.  We'll need to put back the reads we remove
    //  before we add it to the store.

    if (ignoreContains)
      fl = stashContains(ma, ignoreContainT);

    if (MultiAlignUnitig(ma, gkpStore, &options, NULL)) {
      if (showResult)
        PrintMultiAlignT(stdout, ma, gkpStore, false, false, AS_READ_CLEAR_LATEST);

      if (ignoreContains)
        unstashContains(ma, fl);

      tigStore->insertMultiAlign(ma, TRUE, FALSE);

      DeleteMultiAlignT(ma);
    } else {
      fprintf(stderr, "MultiAlignUnitig()-- unitig %d failed.\n", ma->maID);
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
