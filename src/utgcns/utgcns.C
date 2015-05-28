
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
#include "gkStore.H"
#include "tgStore.H"
#include "abAbacus.H"

#include "AS_UTL_decodeRange.H"
#include "AS_UTL_fasta.C"

#include "stashContains.H"

#include <map>
#include <algorithm>


int
main (int argc, char **argv) {
  char  *gkpName = NULL;

  char   *tigName = NULL;
  uint32  tigVers = UINT32_MAX;
  uint32  tigPart = UINT32_MAX;

  uint32  utgBgn  = UINT32_MAX;
  uint32  utgEnd  = UINT32_MAX;
  char   *utgFile = NULL;

  char *outResultsName = NULL;
  char *outLayoutsName = NULL;
  char *outSeqName     = NULL;

  FILE *outResultsFile = NULL;
  FILE *outLayoutsFile = NULL;
  FILE *outSeqFile     = NULL;

  bool   forceCompute = false;

  double errorRate    = 0.06;
  double errorRateMax = 0.40;
  uint32 minOverlap   = 40;

  int32  numFailures = 0;
  int32  numSkipped  = 0;

  bool   showResult = false;

  double maxCov = 0.0;
  uint32 maxLen = UINT32_MAX;

  bool   inplace  = false;
  bool   loadall  = false;

  uint32 verbosity = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);
      tigPart = atoi(argv[++arg]);

      if (argv[arg][0] == '.')
        tigPart = UINT32_MAX;

      if (tigVers == 0)
        fprintf(stderr, "invalid tigStore version (-T store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);
      if (tigPart == 0)
        fprintf(stderr, "invalid tigStore partition (-T store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);

    } else if (strcmp(argv[arg], "-u") == 0) {
      AS_UTL_decodeRange(argv[++arg], utgBgn, utgEnd);

    } else if (strcmp(argv[arg], "-t") == 0) {
      utgFile = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      outResultsName = argv[++arg];

    } else if (strcmp(argv[arg], "-L") == 0) {
      outLayoutsName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      outSeqName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      errorRate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-em") == 0) {
      errorRateMax = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      minOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-f") == 0) {
      forceCompute = true;

    } else if (strcmp(argv[arg], "-v") == 0) {
      showResult = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      verbosity++;

    } else if (strcmp(argv[arg], "-maxcoverage") == 0) {
      maxCov   = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-maxlength") == 0) {
      maxLen   = atof(argv[++arg]);

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
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "    -G gkpStore\n");
    fprintf(stderr, "    -T tigStore version partition\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -O results      Write computed tigs to binary output file 'results'\n");
    fprintf(stderr, "    -L layouts      Write computed tigs to layout output file 'layouts'\n");
    fprintf(stderr, "    -F fastq        Write computed tigs to fastq  output file 'fastq'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -u b            Compute only unitig ID 'b' (must be in the correct partition!)\n");
    fprintf(stderr, "    -u b-e          Compute only unitigs from ID 'b' to ID 'e'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -t file         Test the computation of the unitig layout in 'file'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -f              Recompute unitigs that already have a multialignment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -e e            Expect alignments at up to fraction e error\n");
    fprintf(stderr, "    -em m           Don't ever allow alignments more than fraction m error\n");
    fprintf(stderr, "    -l l            Expect alignments of at least l bases\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -v              Show multialigns.\n");
    fprintf(stderr, "    -V              Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  ADVANCED OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -n              Do not update the store after computing consensus.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -maxcoverage c  Use non-contained reads and the longest contained reads, up to\n");
    fprintf(stderr, "                    C coverage, for consensus generation.  The default is 0, and will\n");
    fprintf(stderr, "                    use all reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -maxlength l    Do not compute consensus for unitigs longer than l bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -inplace        Write the updated unitig to the same version it was read from.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -t S V P        If 'partition' is '.', use an unpartitioned tigStore/gkpStore.\n");
    fprintf(stderr, "    -loadall        Load ALL reads into memory.  Ignores partition if it exists.\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR:  No gkpStore (-g) supplied.\n");

    if ((utgFile == NULL) && (tigName == NULL))
      fprintf(stderr, "ERROR:  No tigStore (-t) OR no test unitig (-T) supplied.\n");

    exit(1);
  }

  errno = 0;

  //  Open output files

  if (outResultsName)
    outResultsFile = fopen(outResultsName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output results file '%s': %s\n", outResultsName, strerror(errno)), exit(1);

  if (outLayoutsName)
    outLayoutsFile = fopen(outLayoutsName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output layout file '%s': %s\n", outLayoutsName, strerror(errno)), exit(1);

  if (outSeqName)
    outSeqFile = fopen(outSeqName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output FASTQ file '%s': %s\n", outSeqName, strerror(errno)), exit(1);

  //  Open gatekeeper for read only, and load the partitioned data if tigPart > 0.

  fprintf(stderr, "Opening gkpStore.\n");
  gkStore   *gkpStore = new gkStore(gkpName, gkStore_readOnly, tigPart);

  //  Create a consensus object.

  fprintf(stderr, "Creating abacus.\n");
  abAbacus  *abacus   = new abAbacus(gkpStore);

  //  If we are testing a unitig, do that.

  if (utgFile != NULL) {
    fprintf(stderr, "utgFile not supported.\n");
    exit(1);
#if 0
    errno = 0;
    FILE  *F = fopen(utgFile, "r");
    if (errno)
      fprintf(stderr, "Failed to open input unitig file '%s': %s\n", utgFile, strerror(errno)), exit(1);

    MultiAlignT  *ma       = CreateEmptyMultiAlignT();
    bool          isUnitig = false;

    while (LoadMultiAlignFromHuman(tig, isUnitig, F) == true) {
      if (generateMultiAlignment(tig, gkpStore, NULL)) {
        if (showResult)
          abacus->getMultiAlignment()->printAlignment(abacus, stdout);

      } else {
        fprintf(stderr, "tig %d failed.\n", tig->tigID());
        numFailures++;
      }
    }

    DeleteMultiAlignT(ma);
#endif

    delete abacus;
    delete gkpStore;

    exit(0);
  }

  //  Otherwise, we're computing unitigs from the store.  Open it for read only.
  //  Outputs get written to a single output file.

  fprintf(stderr, "Opening tigStore.\n");
  tgStore *tigStore = new tgStore(tigName, tigVers);

  //  Decide on what to compute.  Either all unitigs, or a single unitig, or a special case test.

  uint32  b = 0;
  uint32  e = tigStore->numTigs() - 1;

  if (utgEnd > e)
    utgEnd = e;

  if (utgBgn != UINT32_MAX) {
    b = utgBgn;
    e = utgEnd;
  }

  fprintf(stderr, "Computing unitig consensus for b="F_U32" to e="F_U32" with errorRate %0.4f (max %0.4f) and minimum overlap "F_U32"\n",
          b, e, errorRate, errorRateMax, minOverlap);

  //  Now the usual case.  Iterate over all unitigs, compute and update.

  for (uint32 ti=b; ti<=e; ti++) {
    tgTig  *tig = tigStore->loadTig(ti);  //  Store owns the tig

    //fprintf(stderr, "tig %u\n", ti);

    //  Deleted?

    if (tig == NULL)
      continue;

    //  Are we parittioned?  Is this tig in our partition?

    if (tigPart != UINT32_MAX) {
      uint32  missingReads = 0;

      for (uint32 ii=0; ii<tig->numberOfChildren(); ii++)
        if (gkpStore->gkStore_getReadInPartition(tig->getChild(ii)->ident()) == NULL)
          missingReads++;

      if (missingReads) {
        //fprintf(stderr, "tig %u has %u reads and %u not in this partition\n",
        //        ti, tig->numberOfChildren(), missingReads);
        continue;
      }

      fprintf(stderr, "tig %u has %u reads and %u not in this partition\n",
              ti, tig->numberOfChildren(), missingReads);
    }

    bool exists = (tig->gappedLength() > 0);

    if ((forceCompute == false) && (exists == true)) {
      //  Already finished unitig consensus.
      if (tig->numberOfChildren() > 1)
        fprintf(stderr, "Working on unitig %d of length %d (%d children) - already computed, skipped\n",
                tig->tigID(), tig->layoutLength(), tig->numberOfChildren());
      numSkipped++;
      continue;
    }

    if (tig->layoutLength() > maxLen) {
      fprintf(stderr, "SKIP unitig %d of length %d (%d children) - too long, skipped\n",
              tig->tigID(), tig->layoutLength(), tig->numberOfChildren());
      continue;
    }

    if (tig->numberOfChildren() == 0) {
      fprintf(stderr, "SKIP unitig %d of length %d (%d children) - no children, skipped\n",
              tig->tigID(), tig->layoutLength(), tig->numberOfChildren());
      continue;
    }

    if (tig->numberOfChildren() > 1)
      fprintf(stderr, "Working on unitig %d of length %d (%d children)%s\n",
              tig->tigID(), tig->layoutLength(), tig->numberOfChildren(),
              (exists) ? " - already computed, recomputing" : "");

    //  Build a new ma if we're ignoring contains.  We'll need to put back the reads we remove
    //  before we add it to the store.

    savedChildren *origChildren = stashContains(tig, maxCov);

    tig->_utgcns_verboseLevel = verbosity;
    tig->_utgcns_smoothWindow = 11;
    tig->_utgcns_splitAlleles = false;
    tig->_utgcns_doPhasing    = false;

    if (generateMultiAlignment(tig, gkpStore, NULL, errorRate, errorRateMax, minOverlap)) {
      if (showResult)
        tig->display(stdout, gkpStore, 200, 3);

      unstashContains(tig, origChildren);

      //tigStore->insertTig(tig, true);           //  Store owns the tig still.

      if (outResultsFile)
        tig->saveToStream(outResultsFile);

      if (outLayoutsFile)
        tig->dumpLayout(outLayoutsFile);

      if (outSeqFile)
        AS_UTL_writeFastQ(outSeqFile,
                          tig->ungappedBases(), tig->ungappedLength(),
                          tig->ungappedQuals(), tig->ungappedLength(),
                          "@utg%08u\n", tig->tigID());

      tigStore->unloadTig(tig->tigID(), true);  //  Tell the store we're done with it

    } else {
      fprintf(stderr, "unitigConsensus()-- unitig %d failed.\n", tig->tigID());
      numFailures++;
    }
  }

 finish:
  delete abacus;
  delete tigStore;
  delete gkpStore;

  if (outResultsFile)  fclose(outResultsFile);
  if (outLayoutsFile)  fclose(outLayoutsFile);

#if 0
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
#endif

  if (numFailures) {
    fprintf(stderr, "WARNING:  Total number of unitig failures = %d\n", numFailures);
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus did NOT finish successfully.\n");

  } else {
    fprintf(stderr, "Consensus finished successfully.  Bye.\n");
  }

  return(0);
}
