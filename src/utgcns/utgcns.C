
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
#include "strings.H"

#include "sqStore.H"
#include "tgStore.H"

#include "stashContains.H"

#include "unitigConsensus.H"

#ifndef BROKEN_CLANG_OpenMP
#include <omp.h>
#endif
#include <map>
#include <algorithm>


int
main (int argc, char **argv) {
  char    *seqName         = NULL;

  char    *tigName         = NULL;
  uint32   tigVers         = UINT32_MAX;
  uint32   tigPart         = UINT32_MAX;

  char    *tigFileName     = NULL;

  uint32   tigBgn          = 0;
  uint32   tigEnd          = UINT32_MAX;

  char    *outResultsName  = NULL;
  char    *outLayoutsName  = NULL;
  char    *outSeqNameA     = NULL;
  char    *outSeqNameQ     = NULL;

  char    *exportName      = NULL;
  char    *importName      = NULL;

  char      algorithm      = 'P';
  char      aligner        = 'E';

  uint32    numThreads	   = omp_get_max_threads();

  double    errorRate      = 0.12;
  double    errorRateMax   = 0.40;
  uint32    minOverlap     = 40;

  uint32    numFailures    = 0;

  bool      showResult     = false;

  double    maxCov         = 0.0;
  uint32    maxLen         = UINT32_MAX;

  bool      onlyUnassem    = false;
  bool      onlyBubble     = false;
  bool      onlyContig     = false;

  bool      noSingleton    = false;

  uint32    verbosity      = 0;

  sqStore  *seqStore = NULL;
  tgStore  *tigStore = NULL;
  FILE     *tigFile  = NULL;

  FILE     *importFile  = NULL;
  FILE     *exportFile = NULL;

  FILE     *outResultsFile = NULL;
  FILE     *outLayoutsFile = NULL;
  FILE     *outSeqFileA    = NULL;
  FILE     *outSeqFileQ    = NULL;


  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);
      tigPart = atoi(argv[++arg]);

      if (argv[arg][0] == '.')
        tigPart = UINT32_MAX;

      if (tigVers == 0) {
        char *s = new char [1024];
        snprintf(s, 1024, "Invalid tigStore version (-T store version partition) '-T %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]);
        err.push_back(s);
      }

      if (tigPart == 0) {
        char *s = new char [1024];
        snprintf(s, 1024, "Invalid tigStore partition (-T store version partition) '-T %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);
        err.push_back(s);
      }

    } else if ((strcmp(argv[arg], "-u") == 0) ||
               (strcmp(argv[arg], "-tig") == 0)) {
      decodeRange(argv[++arg], tigBgn, tigEnd);

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
    } else if (strcmp(argv[arg], "-norealign") == 0) {
      algorithm = 'p';

    } else if (strcmp(argv[arg], "-edlib") == 0) {
      aligner = 'E';

    } else if (strcmp(argv[arg], "-threads") == 0) {
      numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-export") == 0) {
      exportName = argv[++arg];
    } else if (strcmp(argv[arg], "-import") == 0) {
      importName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      errorRate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-em") == 0) {
      errorRateMax = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      minOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-v") == 0) {
      showResult = true;

    } else if (strncmp(argv[arg], "-V", 2) == 0) {
      verbosity += strlen(argv[arg]) - 1;

    } else if (strcmp(argv[arg], "-maxcoverage") == 0) {
      maxCov   = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-maxlength") == 0) {
      maxLen   = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-onlyunassem") == 0) {
      onlyUnassem = true;

    } else if (strcmp(argv[arg], "-onlybubble") == 0) {
      onlyBubble = true;

    } else if (strcmp(argv[arg], "-onlycontig") == 0) {
      onlyContig = true;

    } else if (strcmp(argv[arg], "-nosingleton") == 0) {
      noSingleton = true;

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((seqName == NULL) && (tigName != NULL)) {
    seqName = new char [FILENAME_MAX];
    snprintf(seqName, FILENAME_MAX, "%s/partitionedReads.seqStore", tigName);
  }


  if ((seqName == NULL) && (importName == NULL))
    err.push_back("ERROR:  No seqStore (-S) and no package (-p) supplied.\n");

  if ((tigFileName == NULL) && (tigName == NULL)  && (importName == NULL))
    err.push_back("ERROR:  No tigStore (-T) OR no test tig (-t) OR no package (-p)  supplied.\n");


  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  INPUT\n");
    fprintf(stderr, "    -S g            Load reads from sqStore 'g'\n");
    fprintf(stderr, "    -T t v p        Load tig from tgStore 't', version 'v', partition 'p'.\n");
    fprintf(stderr, "                      Expects reads will be in sqStore partition 'p' as well\n");
    fprintf(stderr, "                      Use p='.' to specify no partition\n");
    fprintf(stderr, "    -t file         Test the computation of the tig layout in 'file'\n");
    fprintf(stderr, "                      'file' can be from:\n");
    fprintf(stderr, "                        'tgStoreDump -d layout' (human readable layout format)\n");
    fprintf(stderr, "                        'utgcns -L'             (human readable layout format)\n");
    fprintf(stderr, "                        'utgcns -O'             (binary multialignment format)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -import name    Load tig and reads from file 'name' created with -export.  This\n");
    fprintf(stderr, "                    is usually used by developers.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  ALGORITHM\n");
    fprintf(stderr, "    -quick          Stitch reads together to cover the contig.  The bases in the contig\n");
    fprintf(stderr, "                    is formed from exactly one read; no consensus sequence is computed.\n");
    fprintf(stderr, "                    This is useful for checking intermediate assembly structure by mapping\n");
    fprintf(stderr, "                    to reference, or as input to a polishing step.  Read positions will be\n");
    fprintf(stderr, "                    incorrect, and no BAM output is possible.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -pbdagcon       Use pbdagcon (https://github.com/PacificBiosciences/pbdagcon).\n");
    fprintf(stderr, "                    This is fast and robust.  It is the default algorithm.  It does not\n");
    fprintf(stderr, "                    generate a final multialignment output (the -v option will not show\n");
    fprintf(stderr, "                    anything useful).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -norealign      Disable alignment of reads back to the final consensus sequence.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  ALIGNER\n");
    fprintf(stderr, "    -edlib          Myers' O(ND) algorithm from Edlib (https://github.com/Martinsos/edlib).\n");
    fprintf(stderr, "                    This is the default (and, yes, there is no non-default aligner).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  OUTPUT\n");
    fprintf(stderr, "    -O results      Write computed tigs to binary output file 'results'\n");
    fprintf(stderr, "    -L layouts      Write computed tigs to layout output file 'layouts'\n");
    fprintf(stderr, "    -A fasta        Write computed tigs to fasta  output file 'fasta'\n");
    fprintf(stderr, "    -Q fastq        Write computed tigs to fastq  output file 'fastq'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -export name    Create a copy of the inputs needed to compute the tigs.  This\n");
    fprintf(stderr, "                    file can then be sent to the developers for debugging.  The tig(s)\n");
    fprintf(stderr, "                    are not processed and no other outputs are created.  Ideally,\n");
    fprintf(stderr, "                    only one tig is selected (-u, below).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  TIG SELECTION (if -T input is used)\n");
    fprintf(stderr, "    -tig b          Compute only tig ID 'b' (must be in the correct partition!)\n");
    fprintf(stderr, "    -tig b-e        Compute only tigs from ID 'b' to ID 'e'\n");
    fprintf(stderr, "    -u              Alias for -tig\n");
    fprintf(stderr, "    -maxlength l    Do not compute consensus for tigs longer than l bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -onlyunassem    Only compute consensus for unassembled tigs.\n");
    fprintf(stderr, "    -onlybubble     Only compute consensus for bubble tigs (there are no bubbles).\n");
    fprintf(stderr, "    -onlycontig     Only compute consensus for real unitigs/contigs.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -nosingleton    Do not compute consensus for singleton (single-read) tigs.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  PARAMETERS\n");
    fprintf(stderr, "    -e e            Expect alignments at up to fraction e error\n");
    fprintf(stderr, "    -em m           Don't ever allow alignments more than fraction m error\n");
    fprintf(stderr, "    -l l            Expect alignments of at least l bases\n");
    fprintf(stderr, "    -maxcoverage c  Use non-contained reads and the longest contained reads, up to\n");
    fprintf(stderr, "                    C coverage, for consensus generation.  The default is 0, and will\n");
    fprintf(stderr, "                    use all reads.\n");
    fprintf(stderr, "    -threads t      Use 't' compute threads; default 1.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  LOGGING\n");
    fprintf(stderr, "    -v              Show multialigns.\n");
    fprintf(stderr, "    -V              Enable debugging option 'verbosemultialign'.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }


  omp_set_num_threads(numThreads);


  //  Open inputs.

  if (seqName) {
    fprintf(stderr, "-- Opening seqStore '%s' partition %u.\n", seqName, tigPart);
    seqStore = sqStore::sqStore_open(seqName, sqStore_readOnly, tigPart);
  }

  if (tigName) {
    fprintf(stderr, "-- Opening tigStore '%s' version %u.\n", tigName, tigVers);
    tigStore = new tgStore(tigName, tigVers);
  }

  if (tigFileName) {
    fprintf(stderr, "-- Opening tigFile '%s'.\n", tigFileName);
    tigFile = AS_UTL_openInputFile(tigFileName);
  }

  //  Decide which reads to use.  Usually, we want the trimmed reads, but for HiFi reads, we want to
  //  use the original non-homopolymer compressed reads.

  if (seqStore) {
    if (seqStore->sqStore_getLibrary(1)->sqLibrary_readType() == SQ_READTYPE_PACBIO_HIFI) {
      fprintf(stderr, "HIFI\n");
      sqRead_setDefaultVersion(sqRead_raw);
    } else {
      sqRead_setDefaultVersion(sqRead_trimmed);
    }
  }

  //  Open the import/export files.

  if (exportName) {
    fprintf(stderr, "-- Opening output package '%s'.\n", exportName);
    exportFile = AS_UTL_openOutputFile(exportName);
  }

  if (importName) {
    fprintf(stderr, "-- Opening input package '%s'.\n", importName);
    importFile = AS_UTL_openInputFile(importName);
  }

  //  Open output files.  If we're creating a package, the usual output files are not opened.

  if ((exportName == NULL) && (outResultsName)) {
    fprintf(stderr, "-- Opening output results file '%s'.\n", outResultsName);
    outResultsFile = AS_UTL_openOutputFile(outResultsName);
  }

  if ((exportName == NULL) && (outLayoutsName)) {
    fprintf(stderr, "-- Opening output layouts file '%s'.\n", outLayoutsName);
    outLayoutsFile = AS_UTL_openOutputFile(outLayoutsName);
  }

  if ((exportName == NULL) && (outSeqNameA)) {
    fprintf(stderr, "-- Opening output FASTA file '%s'.\n", outSeqNameA);
    outSeqFileA    = AS_UTL_openOutputFile(outSeqNameA);
  }

  if ((exportName == NULL) && (outSeqNameQ)) {
    fprintf(stderr, "-- Opening output FASTQ file '%s'.\n", outSeqNameQ);
    outSeqFileQ    = AS_UTL_openOutputFile(outSeqNameQ);
  }

  //  Open sequence store for read only, and load the partitioned data if tigPart > 0.
  //  Decide on what to compute.  Either all tigs, or a single tig, or a special case test.

  map<uint32, sqRead *>     *importRead     = NULL;
  map<uint32, sqReadData *> *importReadData = NULL;

  if ((tigStore) &&
      (tigEnd > tigStore->numTigs() - 1))
    tigEnd = tigStore->numTigs() - 1;

  //  Write a nice header if we're computing stuff.

  uint32  nTigs       = 0;   //  For reporting at the end.
  uint32  nSingletons = 0;

  if (exportFile == NULL) {
    fprintf(stderr, "--\n");
    fprintf(stderr, "-- Computing consensus for b=" F_U32 " to e=" F_U32 " with errorRate %0.4f (max %0.4f) and minimum overlap " F_U32 "\n",
            tigBgn, tigEnd, errorRate, errorRateMax, minOverlap);
    fprintf(stderr, "--\n");
    fprintf(stdout, "                           ----------CONTAINED READS----------  -DOVETAIL  READS-\n");
    fprintf(stdout, "  tigID    length   reads      used coverage  ignored coverage      used coverage\n");
    fprintf(stdout, "------- --------- -------  -------- -------- -------- --------  -------- --------\n");
  }

  //
  //  If input from a file, either a package or a layout, load and process data until there isn't any more.
  //

  if (importFile) {
    tgTig                     *tig = new tgTig();
    map<uint32, sqRead *>      reads;
    map<uint32, sqReadData *>  datas;

    FILE  *importedLayouts = AS_UTL_openOutputFile(importName, '.', "layout", (importName != NULL));
    FILE  *importedReads   = AS_UTL_openOutputFile(importName, '.', "fasta",  (importName != NULL));

    while (((importFile) && (tig->importData(importFile, reads, datas, importedLayouts, importedReads) == true)) ||
           ((tigFile)    && (tig->loadFromStreamOrLayout(tigFile) == true))) {

      //  Stash excess coverage.

      savedChildren *origChildren = stashContains(tig, maxCov, true);

      //  Compute!

      tig->_utgcns_verboseLevel = verbosity;

      unitigConsensus  *utgcns  = new unitigConsensus(seqStore, errorRate, errorRateMax, minOverlap);
      bool              success = utgcns->generate(tig, algorithm, aligner, &reads, &datas);

      //  Show the result, if requested.

      if (showResult)
        tig->display(stdout, seqStore, 200, 3);

      //  Unstash.

      unstashContains(tig, origChildren);

      //  Save the result.

      if (outResultsFile)   tig->saveToStream(outResultsFile);
      if (outLayoutsFile)   tig->dumpLayout(outLayoutsFile);
      if (outSeqFileA)      tig->dumpFASTA(outSeqFileA);
      if (outSeqFileQ)      tig->dumpFASTQ(outSeqFileQ);

      //  Tidy up for the next tig.

      delete tig;
      tig = new tgTig();    //  Next loop needs an existing empty layout.
    }

    AS_UTL_closeFile(importedReads);
    AS_UTL_closeFile(importedLayouts);
  }

  //
  //  If output to a package file, load and dump data.  No filtering, everything is dumped.
  //

  else if (exportFile) {
    for (uint32 ti=tigBgn; ti<=tigEnd; ti++) {
      tgTig *tig = tigStore->loadTig(ti);

      if (tig) {
        nTigs++;
        tig->exportData(exportFile, seqStore, false);
      }
    }
  }

  //
  //  Otherwise, input is from a tigStore, process all tigs requested.

  else {
    for (uint32 ti=tigBgn; ti<=tigEnd; ti++) {
      tgTig *tig = tigStore->loadTig(ti);

      if ((tig == NULL) ||                  //  Ignore non-existent and
          (tig->numberOfChildren() == 0))   //  empty tigs.
        continue;

      //  Skip stuff we want to skip.

      if (((onlyUnassem == true) && (tig->_class != tgTig_unassembled)) ||
          ((onlyContig  == true) && (tig->_class != tgTig_contig)) ||
          ((onlyBubble  == true) && (tig->_class != tgTig_bubble)) ||
          ((noSingleton == true) && (tig->numberOfChildren() == 1)) ||
          (tig->length() > maxLen))
        continue;

      //  If partitioned, skip this tig if all the reads aren't in this partition.

      if (tigPart != UINT32_MAX) {
        uint32  missingReads = 0;

        for (uint32 ii=0; ii<tig->numberOfChildren(); ii++)
          if (seqStore->sqStore_readInPartition(tig->getChild(ii)->ident()) == false)
            missingReads++;

        if (missingReads)
          continue;
      }

      //  Log that we're processing.

      if (tig->numberOfChildren() > 1) {
        fprintf(stdout, "%7u %9u %7u", tig->tigID(), tig->length(), tig->numberOfChildren());
      }

      //  Stash excess coverage.

      savedChildren *origChildren = stashContains(tig, maxCov, true);

      if (origChildren != NULL) {
        nTigs++;
        fprintf(stdout, "  %8u %7.2fx %8u %7.2fx  %8u %7.2fx\n",
                origChildren->numContainsSaved,    origChildren->covContainsSaved,
                origChildren->numContainsRemoved,  origChildren->covContainsRemoved,
                origChildren->numDovetails,        origChildren->covDovetail);
      } else {
        nSingletons++;
      }

      //  Compute!

      tig->_utgcns_verboseLevel = verbosity;

      unitigConsensus  *utgcns  = new unitigConsensus(seqStore, errorRate, errorRateMax, minOverlap);
      bool              success = utgcns->generate(tig, algorithm, aligner);

      //  Show the result, if requested.

      if (showResult)
        tig->display(stdout, seqStore, 200, 3);

      //  Unstash.

      unstashContains(tig, origChildren);

      //  Save the result.

      if (outResultsFile)   tig->saveToStream(outResultsFile);
      if (outLayoutsFile)   tig->dumpLayout(outLayoutsFile);
      if (outSeqFileA)      tig->dumpFASTA(outSeqFileA);
      if (outSeqFileQ)      tig->dumpFASTQ(outSeqFileQ);

      //  Count failure.

      if (success == false) {
        fprintf(stderr, "unitigConsensus()-- tig %d failed.\n", tig->tigID());
        numFailures++;
      }

      //  Tidy up for the next tig.

      delete utgcns;        //  No real reason to keep this until here.
      delete origChildren;  //  Need to keep it until after we display() above.

      tigStore->unloadTig(tig->tigID(), true);  //  Tell the store we're done with it
    }
  }

  delete tigStore;

  seqStore->sqStore_close();

  AS_UTL_closeFile(tigFile, tigFileName);

  AS_UTL_closeFile(outResultsFile, outResultsName);
  AS_UTL_closeFile(outLayoutsFile, outLayoutsName);

  AS_UTL_closeFile(outSeqFileA, outSeqNameA);
  AS_UTL_closeFile(outSeqFileQ, outSeqNameQ);

  AS_UTL_closeFile(exportFile, exportName);
  AS_UTL_closeFile(importFile, importName);

  if (exportName != NULL) {
    fprintf(stdout, "\n");
    fprintf(stderr, "Exported %u tig%s to file '%s'.\n", nTigs, (nTigs == 1) ? "" : "s", exportName);
  }

  if (exportName == NULL) {
    fprintf(stdout, "\n");
    fprintf(stdout, "Processed %u tig%s and %u singleton%s.\n",
            nTigs, (nTigs == 1)             ? "" : "s",
            nSingletons, (nSingletons == 1) ? "" : "s");
    fprintf(stdout, "\n");

    if (numFailures) {
      fprintf(stderr, "WARNING:  %u tig%s failed.\n", numFailures, (numFailures == 1) ? "" : "s");
      fprintf(stderr, "\n");
      fprintf(stderr, "Consensus did NOT finish successfully.\n");
    } else {
      fprintf(stderr, "Consensus finished successfully.\n");
    }
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(numFailures != 0);
}
