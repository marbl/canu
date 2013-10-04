
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

const char *mainid = "$Id$";

#include "AS_global.H"
#include "MultiAlign.H"
#include "MultiAlignStore.H"
#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"


int
main (int argc, char **argv) {
  char  *gkpName = NULL;

  char  *tigName = NULL;
  int32  tigVers = 0;
  int32  tigPart = 0;

  AS_IID bgnID = 0;
  AS_IID onlID = UINT32_MAX;
  AS_IID endID = UINT32_MAX;

  bool   showResult = false;
  bool   doUpdate   = true;
  bool   doNothing  = false;
  bool   doCache    = false;

  char  *outputName = NULL;
  FILE  *outputFile = NULL;

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

      if (argv[++arg][0] != '.')
        tigPart = atoi(argv[arg]);

      if (tigVers <= 0)
        fprintf(stderr, "invalid tigStore version (-t store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);
      if ((tigPart <= 0) && (argv[arg][0] != '.'))
        fprintf(stderr, "invalid tigStore partition (-t store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);

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

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputName = argv[++arg];

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((tigPart > 0) && (outputName == NULL))
    err++;
  if ((err) || (gkpName == NULL) || (tigName == NULL)) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version partition [opts]\n", argv[0]);
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
    fprintf(stderr, "\n");
    fprintf(stderr, "    -o           Partitioned output file.  If 'partition' is not '.' or '0' this must\n");
    fprintf(stderr, "                 be supplied.\n");
    fprintf(stderr, "\n");
    if ((tigPart > 0) && (outputName == NULL))
      fprintf(stderr, "ERROR: Output name (-o) required if operating on a partition (-t store name part)\n");
    exit(1);
  }

  if (outputName) {
    fprintf(stderr, "Saving fixed unitigs to '%s'; store NOT updated.\n", outputName);

    doUpdate = false;

    errno = 0;
    outputFile = fopen(outputName, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s' for write: %s\n", outputName, strerror(errno)), exit(1);
  }


  gkpStore = new gkStore(gkpName, false, false);
  tigStore = new MultiAlignStore(tigName, tigVers, tigPart, 0, doUpdate, false, false);

  if (tigPart != 0) {
    fprintf(stderr, "Loading reads into memory.\n");
    gkpStore->gkStore_loadPartition(tigPart);
  }

  else if (doCache) {
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
    MultiAlignT  *maOrig = tigStore->loadMultiAlign(i, true);

    if (maOrig == NULL) {
      //  Not in our partition, or deleted.
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
    int32 *failmap   = new int32 [GetNumIntMultiPoss(maOrig->f_list)];
    int32 *failag    = new int32 [GetNumIntMultiPoss(maOrig->f_list)];

    memset(failed,  0, sizeof(int32) * GetNumIntMultiPoss(maOrig->f_list));
    memset(failmap, 0, sizeof(int32) * GetNumIntMultiPoss(maOrig->f_list));
    memset(failag,  0, sizeof(int32) * GetNumIntMultiPoss(maOrig->f_list));

    while (GetNumIntMultiPoss(maTest->f_list) > 0) {

      //  Compute consensus.  We don't care about the result, just what fragments failed.
      //
      MultiAlignUnitig(maTest, gkpStore, &options, failed);

      //  Build the 'fixed' unitig using the fragments placed in 'test', and the 'next' unitig
      //  using those not placed.

    tryAgain:
      fprintf(stderr, " - Fixing unitig %d with "F_SIZE_T" fragments\n",
              maOrig->maID, GetNumIntMultiPoss(maTest->f_list));

      ResetVA_IntMultiPos(maFixd->f_list);
      ResetVA_IntMultiPos(maNext->f_list);

      ResetVA_IntUnitigPos(maFixd->u_list);
      ResetVA_IntUnitigPos(maNext->u_list);

      for (uint32 i=0; i < GetNumIntMultiPoss(maTest->f_list); i++) {
        IntMultiPos *imp = GetIntMultiPos(maTest->f_list, i);

        if (failed[i]) {
          AppendVA_IntMultiPos(maNext->f_list, imp);
          fprintf(stderr, "   - save fragment idx=%d ident=%d for next pass\n", i, imp->ident);
        } else {
          failmap[GetNumIntMultiPoss(maFixd->f_list)] = i;  //  Map read failag[] in maFixd to read i in the test unitig
          AppendVA_IntMultiPos(maFixd->f_list, imp);
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

      if (GetNumIntMultiPoss(maNext->f_list) == 0)
        fprintf(stderr, " - Testing new unitig %d with "F_SIZE_T" fragments\n",
                maFixd->maID, GetNumIntMultiPoss(maFixd->f_list));
      else
        fprintf(stderr, " - Testing new unitig %d with "F_SIZE_T" fragments ("F_SIZE_T" remain for next pass)\n",
                maFixd->maID, GetNumIntMultiPoss(maFixd->f_list), GetNumIntMultiPoss(maNext->f_list));

      assert(GetNumIntMultiPoss(maFixd->f_list) > 0);

      if (MultiAlignUnitig(maFixd, gkpStore, &options, failag) == false) {
        fprintf(stderr, "   - Unitig %d failed, again!\n", maFixd->maID);

        for (uint32 i=0; i<GetNumIntMultiPoss(maFixd->f_list); i++) {
          if (failag[i]) {
            fprintf(stderr, "     - fragment ident=%d failed.\n",
                    GetIntMultiPos(maFixd->f_list, i)->ident);
            failed[failmap[i]]++;

            assert(GetIntMultiPos(maTest->f_list, failmap[i])->ident == GetIntMultiPos(maFixd->f_list, i)->ident);
          }
        }

        //failed[lastAdded] = 1;
        //assert(0);
        goto tryAgain;
      }

      //  Add the new unitig to the store (the asserts are duplicated in cgw too).

      assert(1            == GetNumIntUnitigPoss(maFixd->u_list));        //  One unitig in the list (data.num_unitigs is updated on insertion)
      assert(maFixd->maID == GetIntUnitigPos(maFixd->u_list, 0)->ident);  //  Unitig has correct ident

      if (doUpdate) {
        tigStore->insertMultiAlign(maFixd, true, false);
        fprintf(stderr, " - Added unitig %d with "F_SIZE_T" fragments.\n",
                maFixd->maID, GetNumIntMultiPoss(maFixd->f_list));
      } else {
        //  Usually set by insertMultiAlign, need to force it here
        maFixd->data.num_frags   = GetNumIntMultiPoss(maFixd->f_list);
        maFixd->data.num_unitigs = GetNumIntUnitigPoss(maFixd->u_list);
      }

      if (outputFile) {
        maFixd->maID = -1;
        DumpMultiAlignForHuman(outputFile, maFixd, true);
        fprintf(stderr, " - Saved new unitig with "F_SIZE_T" fragments.\n",
                GetNumIntMultiPoss(maFixd->f_list));
      }

      if (showResult)
        PrintMultiAlignT(stdout, maFixd, gkpStore, false, false, AS_READ_CLEAR_LATEST);

      //  Update the test unitig and continue splitting.

      maFixd->maID = maOrig->maID;  //  For diagnostic output.

      ResetVA_IntMultiPos(maTest->f_list);
      ReuseClone_VA(maTest->f_list, maNext->f_list);
    }  //  Until the test unitig is empty

    delete [] failed;
    delete [] failmap;
    delete [] failag;

    //  Now mark the original unitig as deleted.  If we're writing to human output, a little more work
    //  is needed to remove reads first.  The loader notices this and performs the delete.

    if (outputFile) {
      ResetVA_IntMultiPos(maOrig->f_list);
      ResetVA_IntUnitigPos(maOrig->u_list);

      maOrig->data.num_frags   = 0;  //  Usually set by insertMultiAlign, need to force it here
      maOrig->data.num_unitigs = 0;

      DumpMultiAlignForHuman(outputFile, maOrig, true);
    }

    if (doUpdate)
      tigStore->deleteMultiAlign(maOrig->maID, true);
    else
      tigStore->unloadMultiAlign(maOrig->maID, true);
  }

 finish:
  delete tigStore;

  if (outputFile)
    fclose(outputFile);

  fprintf(stderr, "\nConsensus finished successfully.  Bye.\n");
  return(0);
}
