
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2011, J. Craig Venter Institute.
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

const char *mainid = "$Id: utgcnsfix.C,v 1.1 2011-11-29 11:50:00 brianwalenz Exp $";

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

      if (tigVers <= 0)
        fprintf(stderr, "invalid tigStore version (-t store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);

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
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    -v           Show multialigns.\n");
    fprintf(stderr, "    -V           Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, TRUE, FALSE, FALSE);

  uint32  b = 0;
  uint32  e = tigStore->numUnitigs();

  fprintf(stderr, "Checking unitig consensus for b="F_U32" to e="F_U32"\n", b, e);

  for (uint32 i=b; i<e; i++) {
    MultiAlignT  *ma = tigStore->loadMultiAlign(i, TRUE);

    if (ma == NULL) {
      //  Not in our partition, or deleted.
      fprintf(stderr, "NULL ma for unitig "F_U32"\n", i);
      continue;
    }

    if ((ma->consensus != NULL) &&
        (GetNumchars(ma->consensus) > 1))
      //  Has consensus sequence already.
      continue;

    fprintf(stderr, "MultiAlignUnitig()-- Fixing unitig %d (%d unitigs and %d fragments)\n",
            ma->maID, ma->data.num_unitigs, ma->data.num_frags);

    int32 firstFailed = 0;

    while (GetNumIntMultiPoss(ma->f_list) > 0) {

      //  Compute consensus.  We don't care about the result, just 'firstFailed'.  If
      //  consensus returns success, set 'firstFailed' to be all fragments in the unitig.
      //
      if (MultiAlignUnitig(ma, gkpStore, printwhat, &options, firstFailed) == true)
        firstFailed = GetNumIntMultiPoss(ma->f_list);

      //  Build the 'split' unitig with the 'firstFailed' fragments of 'ma'.

      MultiAlignT *split = CreateEmptyMultiAlignT();
      int32        maBgn = firstFailed;
      int32        maEnd = GetNumIntMultiPoss(ma->f_list);

      //  Loop until we find a valid unitig.  It is unlikely that we'll ever actually
      //  loop here, since we just verified that the first 'maBgn' fragments consensus
      //  together, but we're trying to be bullet proof here.

      while (firstFailed > 0) {
        assert(GetNumIntMultiPoss(ma->f_list) >= firstFailed);

        //  Copy fragments to the new unitig
        ResetVA_IntMultiPos(split->f_list);

        for (uint32 i=0; i<maBgn; i++)
          AppendVA_IntMultiPos(split->f_list, GetIntMultiPos(ma->f_list, i));

        //  Make sure it works, and if it does, we keep it.  If not, try again.

        if (MultiAlignUnitig(split, gkpStore, printwhat, &options, firstFailed) == true) {
          assert(firstFailed == 0);
        } else {
          fprintf(stderr, "MultiAlignUnitig()--  unitig %d failed at fragment index %d.\n",
                  ma->maID, firstFailed);

          maBgn = firstFailed;
        }
      }

      assert(maBgn > 0);
      assert(firstFailed == 0);
      assert(GetNumIntMultiPoss(split->f_list) > 0);

      //  Guaranteed to have a new unitig now.  Add it to the store.

      split->maID = tigStore->numUnitigs();

      tigStore->insertMultiAlign(split, TRUE, FALSE);

      fprintf(stderr, "MultiAlignUnitig()--  Added unitig %d with %d fragments.\n",
              split->maID, split->data.num_frags);

      DeleteMultiAlignT(split);

      //  Update the ma unitig and continue splitting.

      for (int32 i=maBgn; i<maEnd; i++)
        SetVA_IntMultiPos(ma->f_list, i - maBgn, GetIntMultiPos(ma->f_list, i));

      ResetToRange_VA(ma->f_list, maEnd - maBgn);
    }  //  Until the original unitig is empty

    //  Now mark 'ma' as deleted.

    tigStore->deleteMultiAlign(ma->maID, true);
  }

 finish:
  delete tigStore;

  if (numFailures) {
    fprintf(stderr, "WARNING:  Total number of unitig failures = %d\n", numFailures);
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus did NOT finish successfully.\n");
    return(1);
  }

  fprintf(stderr, "Consensus finished successfully.  Bye.\n");
  return(0);
}
