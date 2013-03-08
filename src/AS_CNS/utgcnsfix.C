
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

const char *mainid = "$Id: utgcnsfix.C,v 1.8 2013-03-08 19:25:43 brianwalenz Exp $";

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

  AS_IID bgnID = 0;
  AS_IID onlID = UINT32_MAX;
  AS_IID endID = UINT32_MAX;

  bool   showResult = false;
  bool   doUpdate   = true;
  bool   doNothing  = false;
  bool   doCache    = false;

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

    } else if (strcmp(argv[arg], "-u") == 0) {
      onlID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-n") == 0) {
      doUpdate = false;

    } else if (strcmp(argv[arg], "-N") == 0) {
      doNothing = true;

    } else if (strcmp(argv[arg], "-v") == 0) {
      showResult = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE_MULTIALIGN_OUTPUT++;

    } else if (strcmp(argv[arg], "-l") == 0) {
      doCache = 1;

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
    fprintf(stderr, "    -u iid       Only fix unitig 'iid'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -n           Don't update tigStore with any fixes.\n");
    fprintf(stderr, "    -N           Don't do anything, just report which unitigs are broken.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -l           Load the entire gkpStore into memory (faster, but more memory)\n");
    exit(1);
  }

  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, doUpdate, FALSE, FALSE);

  if (doCache) {
    fprintf(stderr, "Loading all reads into memory.\n");
    gkpStore->gkStore_load(0, 0, GKFRAGMENT_QLT);
  }

  bgnID = 0;
  endID = tigStore->numUnitigs();

  if (onlID != UINT32_MAX) {
    if (onlID < endID) {
      bgnID = onlID;
      endID = onlID + 1;
    } else {
      fprintf(stderr, "ERROR: Invalid unitig "F_IID" (-u); only "F_IID" unitigs in the store.\n", onlID, endID);
      exit(1);
    }
  }
 
  fprintf(stderr, "Checking unitig consensus for b="F_U32" to e="F_U32"\n", bgnID, endID);

  for (uint32 i=bgnID; i<endID; i++) {
    MultiAlignT  *maOrig = tigStore->loadMultiAlign(i, TRUE);

    if (maOrig == NULL) {
      //  Not in our partition, or deleted.
      fprintf(stderr, "NULL ma for unitig "F_U32" -- already deleted?\n", i);
      continue;
    }

    if ((maOrig->consensus != NULL) &&
        (GetNumchars(maOrig->consensus) > 1))
      //  Has consensus sequence already.
      continue;

    if (doNothing == true) {
      fprintf(stderr, "unitig %d of length %d with %d fragments is broken\n",
              maOrig->maID, GetMultiAlignLength(maOrig), maOrig->data.num_frags);
      continue;
    }

    fprintf(stderr, "\nEvaluating unitig %d of length %d with %d fragments\n",
            maOrig->maID, GetMultiAlignLength(maOrig), maOrig->data.num_frags);

    //  This is a complicated algorithm, but only to make it bullet-proof.  The basic idea
    //  is to run consensus on the unitig.  Whatever fragments are placed get moved to
    //  a new unitig, reconsensed, and inserted to the store.  Whatever fragments do not
    //  get placed are put into the next unitig we loop on.

    MultiAlignT *maTest = CreateEmptyMultiAlignT();
    MultiAlignT *maNext = CreateEmptyMultiAlignT();
    MultiAlignT *maFixd = CreateEmptyMultiAlignT();

    maTest->maID = maOrig->maID;  //  For diagnostic output only
    maNext->maID = maOrig->maID;
    maFixd->maID = maOrig->maID;

    ReuseClone_VA(maTest->f_list, maOrig->f_list);  //  Copy the fragments to our 'test' unitig

    int32 *failed    = new int32 [GetNumIntMultiPoss(maOrig->f_list)];
    int32  lastAdded = 0;

    while (GetNumIntMultiPoss(maTest->f_list) > 0) {

      //  Compute consensus.  We don't care about the result, just what fragments failed.
      //
      MultiAlignUnitig(maTest, gkpStore, &options, failed);

      //  Build the 'fixed' unitig using the fragments placed in 'test', and the 'next' unitig
      //  using those not placed.

    tryAgain:
      fprintf(stderr, "\nFixing unitig %d ("F_SIZE_T" fragments remain)\n",
              maOrig->maID, GetNumIntMultiPoss(maTest->f_list));

      ResetVA_IntMultiPos(maFixd->f_list);
      ResetVA_IntMultiPos(maNext->f_list);

      ResetVA_IntUnitigPos(maFixd->u_list);
      ResetVA_IntUnitigPos(maNext->u_list);

      for (uint32 i=0; i < GetNumIntMultiPoss(maTest->f_list); i++) {
        IntMultiPos *imp = GetIntMultiPos(maTest->f_list, i);

        if (failed[i]) {
          AppendVA_IntMultiPos(maNext->f_list, imp);
          fprintf(stderr, "  save fragment idx=%d ident=%d for next pass\n", i, imp->ident);
        } else {
          AppendVA_IntMultiPos(maFixd->f_list, imp);
          lastAdded = i;
        }
      }

      //  Test the fixed unitig.  If it succeeds, a new unitig (with ident maID) is created in the
      //  multialign, and we proceed to insert this into the store.
      //
      //  If it fails, mark, in the original maTest unitig, the last fragment added as failed.  We
      //  could (probably) figure out which fragments failed and mark just those, but it would (seem
      //  to) add a lot of complicated (buggy) code.  Plus, we just don't expect anything to
      //  actually fail.

      maFixd->maID = tigStore->numUnitigs();

      fprintf(stderr, "Testing new unitig %d with "F_SIZE_T" fragments ("F_SIZE_T" remain)\n",
              maFixd->maID, GetNumIntMultiPoss(maFixd->f_list), GetNumIntMultiPoss(maNext->f_list));

      assert(GetNumIntMultiPoss(maFixd->f_list) > 0);

      if (MultiAlignUnitig(maFixd, gkpStore, &options, NULL) == false) {
        fprintf(stderr, "Unitig %d failed, again.\n", maFixd->maID);

        failed[lastAdded] = 1;
        //assert(0);
        goto tryAgain;
      }

      //  Add the new unitig to the store (the asserts are duplicated in cgw too).

      assert(1 == GetNumIntUnitigPoss(maFixd->u_list));  //  One unitig in the list (data.num_unitigs is updated on insertion)
      assert(maFixd->maID == GetIntUnitigPos(maFixd->u_list, 0)->ident);  //  Unitig has correct ident

      if (doUpdate)
        tigStore->insertMultiAlign(maFixd, TRUE, FALSE);

      fprintf(stderr, "Added unitig %d with "F_SIZE_T" fragments.\n",
              maFixd->maID, GetNumIntMultiPoss(maFixd->f_list));

      if (showResult)
        PrintMultiAlignT(stdout, maFixd, gkpStore, false, false, AS_READ_CLEAR_LATEST);

      //  Update the test unitig and continue splitting.

      maFixd->maID = maOrig->maID;  //  For diagnostic output.

      ResetVA_IntMultiPos(maTest->f_list);
      ReuseClone_VA(maTest->f_list, maNext->f_list);
    }  //  Until the test unitig is empty

    //  Now mark the original unitig as deleted.

    if (doUpdate)
      tigStore->deleteMultiAlign(maOrig->maID, true);
  }

 finish:
  delete tigStore;

  fprintf(stderr, "Consensus finished successfully.  Bye.\n");
  return(0);
}
