
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_CNS/utgcns.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2009-OCT-05 to 2014-MAR-31
 *      are Copyright 2009-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-30 to 2015-AUG-07
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-28
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "tgStore.H"

#include "AS_UTL_decodeRange.H"

#include "stashContains.H"

#include "unitigConsensus.H"

#include <omp.h>
#include <map>
#include <algorithm>


int
main (int argc, char **argv) {
  char    *gkpName         = NULL;

  char    *tigName         = NULL;
  uint32   tigVers         = UINT32_MAX;
  uint32   tigPart         = UINT32_MAX;

  char    *tigFileName     = NULL;

  uint32   utgBgn          = UINT32_MAX;
  uint32   utgEnd          = UINT32_MAX;

  char    *outResultsName  = NULL;
  char    *outLayoutsName  = NULL;
  char    *outSeqNameA     = NULL;
  char    *outSeqNameQ     = NULL;
  char    *outPackageName  = NULL;

  FILE     *outResultsFile = NULL;
  FILE     *outLayoutsFile = NULL;
  FILE     *outSeqFileA    = NULL;
  FILE     *outSeqFileQ    = NULL;
  FILE     *outPackageFile = NULL;

  char    *inPackageName   = NULL;

  char      algorithm      = 'P';
  uint32    numThreads	   = 0;

  bool      forceCompute   = false;

  double    errorRate      = 0.12;
  double    errorRateMax   = 0.40;
  uint32    minOverlap     = 40;

  int32     numFailures    = 0;

  bool      showResult     = false;

  double    maxCov         = 0.0;
  uint32    maxLen         = UINT32_MAX;

  uint32    verbosity      = 0;

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
      tigFileName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      outResultsName = argv[++arg];

    } else if (strcmp(argv[arg], "-L") == 0) {
      outLayoutsName = argv[++arg];

    } else if (strcmp(argv[arg], "-A") == 0) {
      outSeqNameA = argv[++arg];

    } else if (strcmp(argv[arg], "-Q") == 0) {
      outSeqNameQ = argv[++arg];

    } else if (strcmp(argv[arg], "-quick") == 0) {
      algorithm = 'Q';
    } else if (strcmp(argv[arg], "-pbdagcon") == 0) {
      algorithm = 'P';
    } else if (strcmp(argv[arg], "-utgcns") == 0) {
      algorithm = 'U';

    } else if (strcmp(argv[arg], "-threads") == 0) {
      numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-p") == 0) {
      inPackageName = argv[++arg];

    } else if (strcmp(argv[arg], "-P") == 0) {
      outPackageName = argv[++arg];

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

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }

  if ((gkpName == NULL) && (inPackageName == NULL))
    err++;

  if ((tigFileName == NULL) && (tigName == NULL) && (inPackageName == NULL))
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  INPUT\n");
    fprintf(stderr, "    -G g            Load reads from gkStore 'g'\n");
    fprintf(stderr, "    -T t v p        Load unitigs from tgStore 't', version 'v', partition 'p'.\n");
    fprintf(stderr, "                      Expects reads will be in gkStore partition 'p' as well\n");
    fprintf(stderr, "                      Use p='.' to specify no partition\n");
    fprintf(stderr, "    -t file         Test the computation of the unitig layout in 'file'\n");
    fprintf(stderr, "                      'file' can be from:\n");
    fprintf(stderr, "                        'tgStoreDump -d layout' (human readable layout format)\n");
    fprintf(stderr, "                        'utgcns -L'             (human readable layout format)\n");
    fprintf(stderr, "                        'utgcns -O'             (binary multialignment format)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -p package      Load unitig and read from 'package' created with -P.  This\n");
    fprintf(stderr, "                    is usually used by developers.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  ALGORITHM\n");
    fprintf(stderr, "    -quick          No alignments, just paste read sequence into the unitig positions.\n");
    fprintf(stderr, "                    This is very fast, but the consensus sequence is formed from a mosaic\n");
    fprintf(stderr, "                    of read sequences, and there can be large indel.  This is useful for\n");
    fprintf(stderr, "                    checking intermediate assembly structure by mapping to reference, or\n");
    fprintf(stderr, "                    possibly for use as input to a polishing step.\n");
    fprintf(stderr, "    -pbdagcon       Use pbdagcon (https://github.com/PacificBiosciences/pbdagcon).\n");
    fprintf(stderr, "                    This is fast and robust.  It is the default algorithm.  It does not\n");
    fprintf(stderr, "                    generate a final multialignment output (the -v option will not show\n");
    fprintf(stderr, "                    anything useful).\n");
    fprintf(stderr, "    -utgcns         Use utgcns (the original Celera Assembler consensus algorithm)\n");
    fprintf(stderr, "                    This isn't as fast, isn't as robust, but does generate a final multialign\n");
    fprintf(stderr, "                    output.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  OUTPUT\n");
    fprintf(stderr, "    -O results      Write computed tigs to binary output file 'results'\n");
    fprintf(stderr, "    -L layouts      Write computed tigs to layout output file 'layouts'\n");
    fprintf(stderr, "    -A fasta        Write computed tigs to fasta  output file 'fasta'\n");
    fprintf(stderr, "    -Q fastq        Write computed tigs to fastq  output file 'fastq'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -P package      Create a copy of the inputs needed to compute the unitigs.  This\n");
    fprintf(stderr, "                    file can then be sent to the developers for debugging.  The unitig(s)\n");
    fprintf(stderr, "                    are not processed and no other outputs are created.  Ideally,\n");
    fprintf(stderr, "                    only one unitig is selected (-u, below).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  TIG SELECTION (if -T input is used)\n");
    fprintf(stderr, "    -u b            Compute only unitig ID 'b' (must be in the correct partition!)\n");
    fprintf(stderr, "    -u b-e          Compute only unitigs from ID 'b' to ID 'e'\n");
    fprintf(stderr, "    -f              Recompute unitigs that already have a multialignment\n");
    fprintf(stderr, "    -maxlength l    Do not compute consensus for unitigs longer than l bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  PARAMETERS\n");
    fprintf(stderr, "    -e e            Expect alignments at up to fraction e error\n");
    fprintf(stderr, "    -em m           Don't ever allow alignments more than fraction m error\n");
    fprintf(stderr, "    -l l            Expect alignments of at least l bases\n");
    fprintf(stderr, "    -maxcoverage c  Use non-contained reads and the longest contained reads, up to\n");
    fprintf(stderr, "                    C coverage, for consensus generation.  The default is 0, and will\n");
    fprintf(stderr, "                    use all reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  LOGGING\n");
    fprintf(stderr, "    -v              Show multialigns.\n");
    fprintf(stderr, "    -V              Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");


    if ((gkpName == NULL) && (inPackageName == NULL))
      fprintf(stderr, "ERROR:  No gkpStore (-G) and no package (-p) supplied.\n");

    if ((tigFileName == NULL) && (tigName == NULL)  && (inPackageName == NULL))
      fprintf(stderr, "ERROR:  No tigStore (-T) OR no test unitig (-t) OR no package (-p)  supplied.\n");

    exit(1);
  }

  errno = 0;

  //  Open output files.  If we're creating a package, the usual output files are not opened.

  if (outPackageName)
    outPackageFile = fopen(outPackageName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output package file '%s': %s\n", outPackageName, strerror(errno)), exit(1);

  if ((outResultsName) && (outPackageName == NULL))
    outResultsFile = fopen(outResultsName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output results file '%s': %s\n", outResultsName, strerror(errno)), exit(1);

  if ((outLayoutsName) && (outPackageName == NULL))
    outLayoutsFile = fopen(outLayoutsName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output layout file '%s': %s\n", outLayoutsName, strerror(errno)), exit(1);

  if ((outSeqNameA) && (outPackageName == NULL))
    outSeqFileA = fopen(outSeqNameA, "w");
  if (errno)
    fprintf(stderr, "Failed to open output FASTA file '%s': %s\n", outSeqNameA, strerror(errno)), exit(1);

  if ((outSeqNameQ) && (outPackageName == NULL))
    outSeqFileQ = fopen(outSeqNameQ, "w");
  if (errno)
    fprintf(stderr, "Failed to open output FASTQ file '%s': %s\n", outSeqNameQ, strerror(errno)), exit(1);

  if (numThreads > 0) {
    omp_set_num_threads(numThreads);
    fprintf(stderr, "number of threads     = %d (command line)\n", numThreads);
    fprintf(stderr, "\n");
  } else {
    fprintf(stderr, "number of threads     = %d (OpenMP default)\n", omp_get_max_threads());
    fprintf(stderr, "\n");
  }

  //  Open gatekeeper for read only, and load the partitioned data if tigPart > 0.

  gkStore                   *gkpStore          = NULL;
  tgStore                   *tigStore          = NULL;
  FILE                      *tigFile           = NULL;
  FILE                      *inPackageFile     = NULL;
  map<uint32, gkRead *>     *inPackageRead     = NULL;
  map<uint32, gkReadData *> *inPackageReadData = NULL;

  if (gkpName) {
    fprintf(stderr, "-- Opening gkpStore '%s' partition %u.\n", gkpName, tigPart);
    gkpStore = gkStore::gkStore_open(gkpName, gkStore_readOnly, tigPart);
  }

  if (tigName) {
    fprintf(stderr, "-- Opening tigStore '%s' version %u.\n", tigName, tigVers);
    tigStore = new tgStore(tigName, tigVers);
  }

  if (tigFileName) {
    fprintf(stderr, "-- Opening tigFile '%s'.\n", tigFileName);

    errno = 0;
    tigFile = fopen(tigFileName, "r");
    if (errno)
      fprintf(stderr, "Failed to open input tig file '%s': %s\n", tigFileName, strerror(errno)), exit(1);
  }

  if (inPackageName) {
    fprintf(stderr, "-- Opening package file '%s'.\n", inPackageName);

    errno = 0;
    inPackageFile = fopen(inPackageName, "r");
    if (errno)
      fprintf(stderr, "Failed to open input package file '%s': %s\n", inPackageName, strerror(errno)), exit(1);
  }

  //  Report some sizes.

  fprintf(stderr, "sizeof(abBead)     "F_SIZE_T"\n", sizeof(abBead));
  fprintf(stderr, "sizeof(abColumn)   "F_SIZE_T"\n", sizeof(abColumn));
  fprintf(stderr, "sizeof(abAbacus)   "F_SIZE_T"\n", sizeof(abAbacus));
  fprintf(stderr, "sizeof(abSequence) "F_SIZE_T"\n", sizeof(abSequence));

  //  Decide on what to compute.  Either all unitigs, or a single unitig, or a special case test.

  uint32  b = 0;
  uint32  e = UINT32_MAX;

  if (tigStore) {
    if (utgEnd > tigStore->numTigs() - 1)
      utgEnd = tigStore->numTigs() - 1;

    if (utgBgn != UINT32_MAX) {
      b = utgBgn;
      e = utgEnd;

    } else {
      b = 0;
      e = utgEnd;
    }

    fprintf(stderr, "-- Computing unitig consensus for b="F_U32" to e="F_U32" with errorRate %0.4f (max %0.4f) and minimum overlap "F_U32"\n",
            b, e, errorRate, errorRateMax, minOverlap);
  }

  else {
    fprintf(stderr, "-- Computing unitig consensus with errorRate %0.4f (max %0.4f) and minimum overlap "F_U32"\n",
            errorRate, errorRateMax, minOverlap);
  }

  fprintf(stderr, "\n");

  //  I don't like this loop control.

  for (uint32 ti=b; (e == UINT32_MAX) || (ti <= e); ti++) {
    tgTig  *tig = NULL;

    //  If a tigStore, load the tig.  The tig is the owner; it cannot be deleted by us.
    if (tigStore)
      tig = tigStore->loadTig(ti);

    //  If a tigFile or a package, create a new tig and fill it.  Obviously, we own it.
    if (tigFile || inPackageFile) {
      tig = new tgTig();

      if (tig->loadFromStreamOrLayout((tigFile != NULL) ? tigFile : inPackageFile) == false) {
        delete tig;
        break;
      }
    }

    //  No tig loaded, keep going.

    if (tig == NULL)
      continue;

    //  If a package, populate the read and readData maps with data from the package.

    if (inPackageFile) {
      inPackageRead      = new map<uint32, gkRead *>;
      inPackageReadData  = new map<uint32, gkReadData *>;

      for (int32 ii=0; ii<tig->numberOfChildren(); ii++) {
        uint32       readID = tig->getChild(ii)->ident();
        gkRead      *read   = (*inPackageRead)[readID]     = new gkRead;
        gkReadData  *data   = (*inPackageReadData)[readID] = new gkReadData;

        gkStore::gkStore_loadReadFromStream(inPackageFile, read, data);

        if (read->gkRead_readID() != readID)
          fprintf(stderr, "ERROR: package not in sync with tig.  package readID = %u  tig readID = %u\n",
                  read->gkRead_readID(), readID);
        assert(read->gkRead_readID() == readID);
      }
    }

    //  More 'not liking' - set the verbosity level for logging.

    tig->_utgcns_verboseLevel = verbosity;

    //  Are we parittioned?  Is this tig in our partition?

    if (tigPart != UINT32_MAX) {
      uint32  missingReads = 0;

      for (uint32 ii=0; ii<tig->numberOfChildren(); ii++)
        if (gkpStore->gkStore_getReadInPartition(tig->getChild(ii)->ident()) == NULL)
          missingReads++;

      if (missingReads) {
        //fprintf(stderr, "SKIP unitig %u with %u reads found only %u reads in partition, skipped\n",
        //        tig->tigID(), tig->numberOfChildren(), tig->numberOfChildren() - missingReads);
        continue;
      }
    }

    if (tig->length(true) > maxLen) {
      fprintf(stderr, "SKIP unitig %d of length %d (%d children) - too long, skipped\n",
              tig->tigID(), tig->length(true), tig->numberOfChildren());
      continue;
    }

    if (tig->numberOfChildren() == 0) {
      fprintf(stderr, "SKIP unitig %d of length %d (%d children) - no children, skipped\n",
              tig->tigID(), tig->length(true), tig->numberOfChildren());
      continue;
    }

    bool exists   = tig->consensusExists();

    if (tig->numberOfChildren() > 1)
      fprintf(stderr, "Working on unitig %d of length %d (%d children)%s%s\n",
              tig->tigID(), tig->length(true), tig->numberOfChildren(),
              ((exists == true)  && (forceCompute == false)) ? " - already computed"              : "",
              ((exists == true)  && (forceCompute == true))  ? " - already computed, recomputing" : "");

    //  Process the tig.  Remove deep coverage, create a consensus object, process it, and report the results.
    //  before we add it to the store.

    unitigConsensus  *utgcns       = new unitigConsensus(gkpStore, errorRate, errorRateMax, minOverlap);
    savedChildren    *origChildren = NULL;
    bool              success      = exists;

    //  Save the tig in the package?
    //
    //  The original idea was to dump the tig and all the reads, then load the tig and process as normal.
    //  Sadly, stashContains() rearranges the order of the reads even if it doesn't remove any.  The rearranged
    //  tig couldn't be saved (otherwise it would be rearranged again).  So, we were in the position of
    //  needing to save the original tig and the rearranged reads.  Impossible.
    //
    //  Instead, we save the origianl tig and original reads -- including any that get stashed -- then
    //  load them all back into a map for use in consensus proper.  It's a bit of a pain, and could
    //  have way more reads saved than necessary.

    if (outPackageFile) {
      utgcns->savePackage(outPackageFile, tig);
      fprintf(stderr, "  Packaged unitig %u into '%s'\n", tig->tigID(), outPackageName);
    }

    //  Compute consensus if it doesn't exist, or if we're forcing a recompute.  But only if we
    //  didn't just package it.

    if ((outPackageFile == NULL) &&
        ((exists == false) || (forceCompute == true))) {
      origChildren = stashContains(tig, maxCov, true);

      switch (algorithm) {
        case 'Q':
          success = utgcns->generateQuick(tig, inPackageRead, inPackageReadData);
          break;
        case 'P':
        default:
          success = utgcns->generatePBDAG(tig, inPackageRead, inPackageReadData);
          break;
        case 'U':
          success = utgcns->generate(tig, inPackageRead, inPackageReadData);
          break;
      }
    }

    //  If it was successful (or existed already), output.  Success is always false if the unitig
    //  was packaged, regardless of if it existed already.

    if (success == true) {
      if ((showResult) && (gkpStore))  //  No gkpStore if we're from a package.  Dang.
        tig->display(stdout, gkpStore, 200, 3);

      unstashContains(tig, origChildren);

      if (outResultsFile)
        tig->saveToStream(outResultsFile);

      if (outLayoutsFile)
        tig->dumpLayout(outLayoutsFile);

      if (outSeqFileA)
        tig->dumpFASTA(outSeqFileA, true);

      if (outSeqFileQ)
        tig->dumpFASTQ(outSeqFileQ, true);
    }

    //  Report failures.

    if ((success == false) && (outPackageFile == NULL)) {
      fprintf(stderr, "unitigConsensus()-- unitig %d failed.\n", tig->tigID());
      numFailures++;
    }

    //  Clean up, unloading or deleting the tig.

    delete utgcns;        //  No real reason to keep this until here.
    delete origChildren;  //  Need to keep it until after we display() above.

    if (tigStore)
      tigStore->unloadTig(tig->tigID(), true);  //  Tell the store we're done with it

    if (tigFile)
      delete tig;
  }

 finish:
  delete tigStore;

  gkpStore->gkStore_close();

  if (tigFile)         fclose(tigFile);
  if (outResultsFile)  fclose(outResultsFile);
  if (outLayoutsFile)  fclose(outLayoutsFile);
  if (outPackageFile)  fclose(outPackageFile);
  if (inPackageFile)   fclose(inPackageFile);

  if (numFailures) {
    fprintf(stderr, "WARNING:  Total number of unitig failures = %d\n", numFailures);
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus did NOT finish successfully.\n");

  } else {
    fprintf(stderr, "Consensus finished successfully.  Bye.\n");
  }

  return(numFailures != 0);
}
