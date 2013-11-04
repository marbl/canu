
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

const char *mainid = "$Id$";

#include "AS_global.H"
#include "MultiAlign.H"
#include "MultiAlignStore.H"
#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"

#include "AS_UTL_decodeRange.H"



void
writeToOutFile(char *outName, int32 tigPart, MultiAlignT *ma) {
  errno = 0;

  FILE *outFile = fopen(outName, "a");
  if (errno)
    fprintf(stderr, "ctgcns:  Failed to write contig %d to output file '%s': %s\n",
            ma->maID, outName, strerror(errno)), exit(1);

  fwrite(&ma->maID, sizeof(int32), 1, outFile);
  fwrite(&tigPart,  sizeof(int32), 1, outFile);

  SaveMultiAlignTToStream(ma, outFile);

  fclose(outFile);
}



void
importFromFile(char *inName, int32 tigPart) {
  errno = 0;

  FILE *inFile = fopen(inName, "r");
  if (errno)
    fprintf(stderr, "ctgcns:  Failed to open output file '%s' for input: %s\n",
            inName, strerror(errno)), exit(1);

  int32   maID  = 0;
  int32   pt    = 0;

  fread(&maID, sizeof(int32), 1, inFile);
  fread(&pt,   sizeof(int32), 1, inFile);

  while (!feof(inFile)) {
    MultiAlignT *ma = LoadMultiAlignTFromStream(inFile);

    assert(ma->maID == maID);
    assert(tigPart  == pt);

    fprintf(stderr, "Working on contig %d (%d unitigs and %d fragments) - importing from load file\n",
            ma->maID, ma->data.num_unitigs, ma->data.num_frags);
    tigStore->insertMultiAlign(ma, false, false);

    DeleteMultiAlignT(ma);

    fread(&maID, sizeof(int32), 1, inFile);
    fread(&pt,   sizeof(int32), 1, inFile);
  }
    
  fclose(inFile);
}



int
main (int argc, char **argv) {
  char   tmpName[FILENAME_MAX] = {0};

  char  *gkpName = NULL;

  char  *tigName = NULL;
  int32  tigVers = -1;
  int32  tigPart = -1;

  int64  ctgBgn = -1;
  int64  ctgEnd = -1;

  char  *ctgName = NULL;
  char  *outName = NULL;
  char  *inName  = NULL;

  bool   forceCompute = false;

  int32  numFailures = 0;
  int32  numSkipped  = 0;

  bool   useUnitig  = false;
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
      AS_UTL_decodeRange(argv[++arg], ctgBgn, ctgEnd);

    } else if (strcmp(argv[arg], "-T") == 0) {
      ctgName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      outName = argv[++arg];
    } else if (strcmp(argv[arg], "-I") == 0) {
      inName = argv[++arg];

    } else if (strcmp(argv[arg], "-f") == 0) {
      forceCompute = true;

    } else if (strcmp(argv[arg], "-U") == 0) {
      useUnitig = true;

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
    fprintf(stderr, "    -c b         Compute only contig ID 'b' (must be in the correct partition!)\n");
    fprintf(stderr, "    -c b-e       Compute only contigs from ID 'b' to ID 'e'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -T file      Test the computation of the contig layout in 'file'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -f           Recompute contigs that already have a multialignment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -U           Reuse the unitig consensus for contigs with only a single\n");
    fprintf(stderr, "                 unitig (EXPERIMENTAL!)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -O file      Don't update tigStore, dump a binary file instead.\n");
    fprintf(stderr, "    -I file      Import binary file into tigStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -v           Show multialigns.\n");
    fprintf(stderr, "    -V           Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -w ws        Smoothing window size\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  //  Open both stores for read only.
  gkpStore = new gkStore(gkpName, false, false);
  tigStore = new MultiAlignStore(tigName, tigVers, 0, tigPart, false, false, false);

  gkpStore->gkStore_loadPartition(tigPart);

  //  Decide on what to compute.  Either all contigs, or a single contig, or a special case test.
  uint32 b = 0;
  uint32 e = tigStore->numContigs();

  if (ctgBgn != -1) {
    b = ctgBgn;
    e = ctgEnd + 1;
  }

  FORCE_UNITIG_ABUT = 1;

  if (ctgName != NULL) {
    errno = 0;
    FILE         *F = fopen(ctgName, "r");
    if (errno)
      fprintf(stderr, "Failed to open input contig file '%s': %s\n", ctgName, strerror(errno)), exit(1);

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
  if (((inName) || (b < e)) && (outName == NULL)) {
    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, 0, tigPart, true, false, true);
  }

  if (inName) {
    importFromFile(inName, tigPart);

    b = e = 0;
  }

  //  Now the usual case.  Iterate over all contigs, compute and update.
  for (uint32 i=b; i<e; i++) {
    MultiAlignT  *cma = tigStore->loadMultiAlign(i, false);

    if (cma == NULL) {
      //  Not in our partition, or deleted.
      continue;
    }

    bool  exists = (cma->consensus != NULL) && (GetNumchars(cma->consensus) > 1);

    if ((forceCompute == false) && (exists == true)) {
      //  Already finished contig consensus.
      fprintf(stderr, "Working on contig %d (%d unitigs and %d fragments) - already computed, skipped\n",
              cma->maID, cma->data.num_unitigs, cma->data.num_frags);

      numSkipped++;

      tigStore->unloadMultiAlign(cma->maID, false);

      continue;
    }

    int32         uID = GetIntUnitigPos(cma->u_list, 0)->ident;

    //  If this is a surrogate, we CANNOT reuse the unitig.  We need to process the contig so that
    //  the unplaced reads are stripped out.  A surrogate should have different contig and unitig
    //  IDs; we could also check the contig status.

    if ((cma->data.num_unitigs == 1) &&
        (cma->maID == uID) &&
        (useUnitig == true)) {
      fprintf(stderr, "Working on contig %d (%d unitigs and %d fragments) - reusing unitig %d consensus\n",
              cma->maID, cma->data.num_unitigs, cma->data.num_frags, uID);

      MultiAlignT  *uma = tigStore->loadMultiAlign(uID, true);

      uma->data = cma->data;

      tigStore->unloadMultiAlign(cma->maID, false);

      if (outName)
        writeToOutFile(outName, tigPart, uma);
      else
        tigStore->insertMultiAlign(uma, false, false);

      tigStore->unloadMultiAlign(uma->maID, true);

      continue;
    }

    fprintf(stderr, "Working on contig %d (%d unitigs and %d fragments)%s\n",
            cma->maID, cma->data.num_unitigs, cma->data.num_frags,
            (exists) ? " - already computed, recomputing" : "");

    if (MultiAlignContig(cma, gkpStore, &options)) {
      if (outName)
        writeToOutFile(outName, tigPart, cma);
      else
        tigStore->insertMultiAlign(cma, false, true);

      if (showResult)
        PrintMultiAlignT(stdout, cma, gkpStore, false, false, AS_READ_CLEAR_LATEST);

      tigStore->unloadMultiAlign(cma->maID, false);
    } else {
      fprintf(stderr, "MultiAlignContig()-- contig %d failed.\n", cma->maID);
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
