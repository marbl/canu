
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"
#include "system.H"
#include "strings.H"

#include "merlin-globals.H"

#include "sqStore.H"
#include "tgStore.H"

#include "unitigConsensus.H"
#include "unitigPartition.H"


class cnsParameters {
public:
  void                    closeAndCleanup(void) {
    delete seqStore;   seqStore = NULL;
    delete tigStore;   tigStore = NULL;

    if (seqReads)
      for (auto it=seqReads->begin(); it != seqReads->end(); it++)
        delete it->second;

    delete seqReads;

    AS_UTL_closeFile(outResultsFile, outResultsName);
    AS_UTL_closeFile(outLayoutsFile, outLayoutsName);

    AS_UTL_closeFile(outSeqFileA, outSeqNameA);
    AS_UTL_closeFile(outSeqFileQ, outSeqNameQ);
  };

  char                   *seqName = nullptr;
  char                   *seqFile = nullptr;

  char                   *tigName = nullptr;
  uint32                  tigVers = UINT32_MAX;
  uint32                  tigPart = 0;

  uint32                  tigBgn  = 0;
  uint32                  tigEnd  = UINT32_MAX;

  char                   *outResultsName = nullptr;
  char                   *outLayoutsName = nullptr;
  char                   *outSeqNameA    = nullptr;
  char                   *outSeqNameQ    = nullptr;

  char                   *exportName     = nullptr;
  char                   *importName     = nullptr;
  bool                    dumpImport     = false;

  char                    algorithm = 'P';
  char                    aligner   = 'E';

  bool                    createPartitions = false;
  double                  partitionSize    = 1.00;   //  Size partitions to be 100% of the largest tig.
  double                  partitionScaling = 1.00;   //  Estimated tig length is 100% of actual tig length.
  double                  partitionReads   = 0.05;   //  5% of all reads can end up in a single partition.

  double                  errorRate      = 0.12;
  double                  errorRateMax   = 0.40;
  uint32                  errorRateMaxID = 0;

  uint32                  minOverlap   = 500;
  uint32                  minCoverage  = 0;

  uint32                  numFailures = 0;

  bool                    showResult = false;

  double                  maxCov = 0.0;
  uint32                  minLen = 0;
  uint32                  maxLen = UINT32_MAX;

  bool                    onlyUnassem = false;
  bool                    onlyContig  = false;

  bool                    noBubble    = false;
  bool                    noRepeat    = false;
  bool                    noSingleton = false;

  uint32                  verbosity = 0;

  sqStore                *seqStore = nullptr;
  std::map<uint32, sqRead *>  *seqReads = nullptr;
  tgStore                *tigStore = nullptr;

  FILE                   *outResultsFile = nullptr;
  FILE                   *outLayoutsFile = nullptr;
  FILE                   *outSeqFileA    = nullptr;
  FILE                   *outSeqFileQ    = nullptr;
  merlinGlobal           *merlinGlobal_  = nullptr;
  char                   *markerDB       = nullptr;
};



//  Scan the tigs to compute expected consensus effort.
//
//  The memory estimate is _very_ simple, just 1 GB memory for each 1 Mbp
//  of sequence (which is, of course, 1 KB memory for every base).
void
createPartitions(cnsParameters  &params) {
  tigPartitioning  tp;

  tp.loadTigInfo(params.tigStore, params.verbosity > 1);

  tp.greedilyPartition(params.partitionSize,
                       params.partitionScaling,
                       params.partitionReads,
                       params.verbosity > 0);

  tp.outputPartitions(params.seqStore, params.tigStore, params.tigName);

  FILE   *partFile = AS_UTL_openOutputFile(params.tigName, '/', "partitioning");
  tp.reportPartitioning(partFile);
  AS_UTL_closeFile(partFile);
}



void
printHeader(cnsParameters  &params) {
  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Computing consensus for b=" F_U32 " to e=" F_U32 " with errorRate %0.4f (max %0.4f) and minimum overlap " F_U32 "\n",
          params.tigBgn, params.tigEnd, params.errorRate, params.errorRateMax, params.minOverlap);
  fprintf(stderr, "--\n");
  fprintf(stdout, "                           ----------CONTAINED READS----------  -DOVETAIL  READS-\n");
  fprintf(stdout, "  tigID    length   reads      used coverage  ignored coverage      used coverage\n");
  fprintf(stdout, "------- --------- -------  -------- -------- -------- --------  -------- --------\n");
}


void
processImportedTigs(cnsParameters  &params) {
  readBuffer *importFile      = new readBuffer(params.importName);
  FILE       *importedLayouts = nullptr;
  FILE       *importedReads   = nullptr;

  if (params.dumpImport) {
    importedLayouts = AS_UTL_openOutputFile(params.importName, '.', "layout");
    importedReads   = AS_UTL_openOutputFile(params.importName, '.', "fasta");
  }

  tgTig                       *tig = new tgTig();
  std::map<uint32, sqRead *>   reads;

  sqRead_defaultVersion = sqRead_raw | sqRead_normal;

  while (tig->importData(importFile, reads, importedLayouts, importedReads) == true) {
    if ((params.tigBgn <= tig->tigID()) &&
        (tig->tigID()  <= params.tigEnd)) {
      fprintf(stdout, "%7u %9u %7u%s", tig->tigID(), tig->length(), tig->numberOfChildren(),
              (tig->numberOfChildren() == 1) ? "\n" : "");

      //  No longer stash excess coverage, flag reads to be ignored instead.  Singletons report no logging.

      // tgTigStashed S;

      // tig->stashContains(params.maxCov, S);

      // if (S.nBack > 0)
      //   fprintf(stdout, "  %8u %7.2fx %8u %7.2fx  %8u %7.2fx\n",
      //           S.nCont, (double)S.bCont / tig->length(),
      //           S.nStsh, (double)S.bStsh / tig->length(),
      //           S.nBack, (double)S.bBack / tig->length());

      //  Compute!

      tig->_utgcns_verboseLevel = params.verbosity;

      params.merlinGlobal_->markerDBname = params.markerDB;
      params.merlinGlobal_->load_Kmers(params.markerDB);
      unitigConsensus  *utgcns  = new unitigConsensus(params.seqStore, params.errorRate, params.errorRateMax, params.errorRateMaxID, params.minOverlap, params.minCoverage, params.merlinGlobal_);
      bool              success = utgcns->generate(tig, params.algorithm, params.aligner, params.maxCov, NULL, NULL, &reads);

      //  Show the result, if requested.

      if (params.showResult)
        tig->display(stdout, params.seqStore, 200, 3);

      //  Unstash.
      // tig->unstashContains();

      //  Save the result.

      if (params.outResultsFile)   tig->saveToStream(params.outResultsFile);
      if (params.outLayoutsFile)   tig->dumpLayout(params.outLayoutsFile);
      if (params.outSeqFileA)      tig->dumpFASTA(params.outSeqFileA);
      if (params.outSeqFileQ)      tig->dumpFASTQ(params.outSeqFileQ);

      delete utgcns;
      utgcns = nullptr;
    }
    else {
      fprintf(stdout, "%7u %9u %7u  skipped.\n", tig->tigID(), tig->length(), tig->numberOfChildren());
    }

    //  Tidy up for the next tig.

    for (auto it=reads.begin(); it != reads.end(); ++it) {
      delete it->second;
    }
    reads.clear();

    delete tig;
    tig = new tgTig();    //  Next loop needs an existing empty layout.
  }

  delete tig;

  AS_UTL_closeFile(importedReads);
  AS_UTL_closeFile(importedLayouts);

  delete importFile;
}



void
exportTigs(cnsParameters  &params) {

  fprintf(stderr, "-- Opening output package '%s'.\n", params.exportName);

  writeBuffer *exportFile = new writeBuffer(params.exportName, "w");
  uint32       nTigs      = 0;

  for (uint32 ti=params.tigBgn; ti<=params.tigEnd; ti++) {
    tgTig *tig = params.tigStore->loadTig(ti);

    if (tig) {
      nTigs++;
      tig->exportData(exportFile, params.seqStore, false);
    }
  }

  delete exportFile;

  fprintf(stdout, "\n");
  fprintf(stderr, "Exported %u tig%s to file '%s'.\n", nTigs, (nTigs == 1) ? "" : "s", params.exportName);
}



std::set<uint32>
loadProcessList(char *prefix, uint32 tigPart) {
  std::set<uint32>   processList;
  uint32             Lmax = 1024;
  uint32             Llen = 0;
  char              *L    = new char [Lmax];
  char              *N    = new char [FILENAME_MAX + 1];

  snprintf(N, FILENAME_MAX, "%s/partitioning", prefix);

  if ((tigPart > 0) &&             //  Partitioning requested, and
      (fileExists(N) == true)) {   //  partitioning file exists, load it.
    FILE *F = AS_UTL_openInputFile(N);

    while (AS_UTL_readLine(L, Llen, Lmax, F)) {
      splitToWords S(L);

      if (S.touint32(5) == tigPart)
        processList.insert(S.touint32(0));
    }

    AS_UTL_closeFile(F, N);
  }

  delete [] N;
  delete [] L;

  return(processList);
}



std::map<uint32, sqRead *> *
loadPartitionedReads(char *seqFile) {

  if (seqFile == NULL)
    return(NULL);

  //  Allocate space for the reads, and buffers to load them.

  std::map<uint32, sqRead *>  *reads = new std::map<uint32, sqRead *>;
  readBuffer                  *rb    = new readBuffer(seqFile);
  sqRead                      *rd    = new sqRead;

  uint64 magc;
  uint64 vers;
  uint64 defv;

  //  Read the header.

  if (rb->readIFFobject("MAGC", magc) == false)
    fprintf(stderr, "File '%s' isn't a utgcns seqFile: no magic number found.\n", seqFile), exit(1);

  if (magc != 0x5f5f656c69467173llu)
    fprintf(stderr, "File '%s' isn't a utgcns seqFile: found magic 0x%016lx.\n", seqFile, magc), exit(1);

  if (rb->readIFFobject("VERS", vers) == false)
    fprintf(stderr, "File '%s' isn't a utgcns seqFile: no file version found.\n", seqFile), exit(1);

  if (vers != 0x0000000000000001llu)
    fprintf(stderr, "File '%s' is a utgcns seqFile, but an unsupported version %lu.\n", seqFile, vers), exit(1);

  if (rb->readIFFobject("DEFV", defv) == false)
    fprintf(stderr, "File '%s' isn't a utgcns seqFile: no default version found.\n", seqFile), exit(1);

  sqRead_defaultVersion = (sqRead_which)defv;

  fprintf(stderr, "Loading %s reads from seqFile '%s'\n", toString(sqRead_defaultVersion), seqFile);

  //  Read the reads.

  while (sqStore::sqStore_loadReadFromBuffer(rb, rd) == true) {
    (*reads)[rd->sqRead_readID()] = rd;

    rd = new sqRead;
  }

  delete rd;
  delete rb;

  //  Return the reads.

  return(reads);
}



void
processTigs(cnsParameters  &params) {
  uint32   nTigs       = 0;
  uint32   nSingletons = 0;
  uint32   numFailures = 0;

  //  Load the partition file, if it exists.

  std::set<uint32>   processList = loadProcessList(params.tigName, params.tigPart);

  //  Load the partitioned reads, if they exist.

  params.seqReads = loadPartitionedReads(params.seqFile);

  //  Loop over all tigs, loading each one and processing if requested.

  for (uint32 ti=params.tigBgn; ti<=params.tigEnd; ti++) {

    if ((processList.size() > 0) &&       //  Ignore tigs not in our partition.
        (processList.count(ti) == 0))     //  (if a partition exists)
      continue;

    tgTig *tig = params.tigStore->loadTig(ti);

    if ((tig == NULL) ||                  //  Ignore non-existent and
        (tig->numberOfChildren() == 0))   //  empty tigs.
      continue;

    //  Skip stuff we want to skip.

    if (((params.onlyUnassem == true) && (tig->_class != tgTig_unassembled)) ||
        ((params.onlyContig  == true) && (tig->_class != tgTig_contig)) ||
        ((params.noSingleton == true) && (tig->numberOfChildren() == 1)) ||
        (tig->length() < params.minLen) ||
        (tig->length() > params.maxLen))
      continue;

    //  Skip repeats and bubbles.

    if (((params.noRepeat == true) && (tig->_suggestRepeat == true)) ||
        ((params.noBubble == true) && (tig->_suggestBubble == true)))
      continue;

    //  Log that we're processing.

    if (tig->numberOfChildren() > 1) {
      fprintf(stdout, "%7u %9u %7u", tig->tigID(), tig->length(), tig->numberOfChildren());
    }

    //  No longer stash excess coverage, flag reads to be ignored instead.  Singletons report no logging.

    // tgTigStashed S;

    // tig->stashContains(params.maxCov, S);

    // if (S.nBack > 0) {
    //   nTigs++;
    //   fprintf(stdout, "  %8u %7.2fx %8u %7.2fx  %8u %7.2fx\n",
    //           S.nCont, (double)S.bCont / tig->length(),
    //           S.nStsh, (double)S.bStsh / tig->length(),
    //           S.nBack, (double)S.bBack / tig->length());
    // } else {
    //   nSingletons++;
    // }

    //  Compute!

    tig->_utgcns_verboseLevel = params.verbosity;

    params.merlinGlobal_->markerDBname = params.markerDB;
    params.merlinGlobal_->load_Kmers(params.markerDB);
    unitigConsensus  *utgcns  = new unitigConsensus(params.seqStore, params.errorRate, params.errorRateMax, params.errorRateMaxID, params.minOverlap, params.minCoverage, params.merlinGlobal_);
    bool              success = utgcns->generate(tig, params.algorithm, params.aligner, params.maxCov, &nTigs, &nSingletons, params.seqReads);

    //  Show the result, if requested.

    if (params.showResult)
      tig->display(stdout, params.seqStore, 200, 3);

    //  Unstash.

    // tig->unstashContains();

    //  Save the result.

    if (params.outResultsFile)   tig->saveToStream(params.outResultsFile);
    if (params.outLayoutsFile)   tig->dumpLayout(params.outLayoutsFile);
    if (params.outSeqFileA)      tig->dumpFASTA(params.outSeqFileA);
    if (params.outSeqFileQ)      tig->dumpFASTQ(params.outSeqFileQ);

    //  Count failure.

    if (success == false) {
      fprintf(stderr, "unitigConsensus()-- tig %d failed.\n", tig->tigID());
      numFailures++;
    }

    //  Tidy up for the next tig.

    delete utgcns;        //  No real reason to keep this until here.

    params.tigStore->unloadTig(tig->tigID(), true);  //  Tell the store we're done with it
  }

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



int
main(int argc, char **argv) {
  cnsParameters  params;

  params.merlinGlobal_ = new merlinGlobal(argv[0]);

  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if      (strcmp(argv[arg], "-S") == 0) {
      params.seqName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-C") == 0){
       params.minCoverage = strtouint32(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-R") == 0) {
      params.seqFile = argv[++arg];
    }

    else if (strcmp(argv[arg], "-T") == 0) {
      params.tigName = argv[++arg];
      params.tigVers = strtouint32(argv[++arg]);

      if (params.tigVers == 0) {
        char *s = new char [1024];
        snprintf(s, 1024, "Invalid tigStore version (-T store v) '-T %s %s'.\n", argv[arg-1], argv[arg]);
        err.push_back(s);
      }
    }

    else if (strcmp(argv[arg], "-P") == 0) {
      params.tigPart = strtouint32(argv[++arg]);
    }


    else if ((strcmp(argv[arg], "-u") == 0) ||
             (strcmp(argv[arg], "-tig") == 0)) {
      decodeRange(argv[++arg], params.tigBgn, params.tigEnd);
    }

    else if (strcmp(argv[arg], "-O") == 0) {
      params.outResultsName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-L") == 0) {
      params.outLayoutsName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-A") == 0) {
      params.outSeqNameA = argv[++arg];
    }

    else if (strcmp(argv[arg], "-Q") == 0) {
      params.outSeqNameQ = argv[++arg];
    }

    //  Partition options

    else if (strcmp(argv[arg], "-partition") == 0) {
      params.createPartitions = true;
      params.partitionSize    = strtodouble(argv[++arg]);
      params.partitionScaling = strtodouble(argv[++arg]);
      params.partitionReads   = strtodouble(argv[++arg]);
    }

    //  Algorithm options

    else if (strcmp(argv[arg], "-quick") == 0) {
      params.algorithm = 'Q';
    }

    else if (strcmp(argv[arg], "-pbdagcon") == 0) {
      params.algorithm = 'P';
    }

    else if (strcmp(argv[arg], "-norealign") == 0) {
      params.algorithm = 'p';
    }

    else if (strcmp(argv[arg], "-edlib") == 0) {
      params.aligner = 'E';
    }

    else if (strcmp(argv[arg], "-threads") == 0) {
      setNumThreads(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-export") == 0) {
      params.exportName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-import") == 0) {
      params.importName = argv[++arg];
    }

    else if (strcmp(argv[arg], "-dumpimport") == 0) {
      params.dumpImport = true;
    }

    else if (strcmp(argv[arg], "-e") == 0) {
      params.errorRate = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-em") == 0) {
      params.errorRateMax = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-EM") == 0) {
      params.errorRateMaxID = strtouint32(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-l") == 0) {
      params.minOverlap = strtouint32(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-v") == 0) {
      params.showResult = true;
    }

    else if (strncmp(argv[arg], "-V", 2) == 0) {
      params.verbosity += strlen(argv[arg]) - 1;
    }

    else if (strcmp(argv[arg], "-maxcoverage") == 0) {
      params.maxCov   = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-minlength") == 0) {
      params.minLen   = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-maxlength") == 0) {
      params.maxLen   = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-onlyunassem") == 0) {
      params.onlyUnassem = true;
    }

    else if (strcmp(argv[arg], "-onlycontig") == 0) {
      params.onlyContig = true;
    }

    else if (strcmp(argv[arg], "-norepeat") == 0) {
      params.noSingleton = true;
    }
    else if (strcmp(argv[arg], "-nobubble") == 0) {
      params.noSingleton = true;
    }
    else if (strcmp(argv[arg], "-nosingleton") == 0) {
      params.noSingleton = true;
    }
    else if (strcmp(argv[arg], "-markers") == 0) {
      params.markerDB = argv[++arg];
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }


  if ((params.seqName == NULL) && (params.importName == NULL) && (params.seqFile == NULL))
    err.push_back("ERROR:  No sequence data!  Need one of seqStore (-S), read file (-R) or package (-p).\n");

  if ((params.tigName == NULL)  && (params.importName == NULL))
    err.push_back("ERROR:  No tigStore (-T) OR no test tig (-t) OR no package (-p) supplied.\n");


  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  INPUT\n");
    fprintf(stderr, "    -S g            Load reads from seqStore 'g'\n");
    fprintf(stderr, "    -R f            Load reads from partition file 'f'\n");
    fprintf(stderr, "    -T t v          Load tig from tigStore 't'.\n");
    fprintf(stderr, "    -t file         Test the computation of the tig layout in 'file'\n");
    fprintf(stderr, "                      'file' can be from:\n");
    fprintf(stderr, "                        'tgStoreDump -d layout' (human readable layout format)\n");
    fprintf(stderr, "                        'utgcns -L'             (human readable layout format)\n");
    fprintf(stderr, "                        'utgcns -O'             (binary multialignment format)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -import name    Load tig and reads from file 'name' created with -export.  This\n");
    fprintf(stderr, "                    is usually used by developers (and Verkko).\n");
    fprintf(stderr, "    -dumpimport     Write the layout and reads from the import file to files\n");
    fprintf(stderr, "                    'name.layout' and 'name.fasta'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -partition a b c\n");
    fprintf(stderr, "                    Create partitions in the tigStore.  Canu uses a=0.8 b=1.0 c=0.1.\n");
    fprintf(stderr, "                      a - Set partition size to be 'a * largest_tig'.  Any tig larger\n");
    fprintf(stderr, "                          than this size is placed entirely in one partition; it is not\n");
    fprintf(stderr, "                          split between partitions.\n");
    fprintf(stderr, "                      b - Scale each tig by 'b' when computing its size.  Only really useful\n");
    fprintf(stderr, "                          for adjusting for homopolymer compression; b=1.5 suggested.\n");
    fprintf(stderr, "                      c - Allow up to 'c * NR' reads per partition, where NR is the number\n");
    fprintf(stderr, "                          of reads in the assembly.\n");
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
    fprintf(stderr, "  TIG SELECTION (if -T input is used)\n");
    fprintf(stderr, "    -tig b          Compute only tig ID 'b' (must be in the correct partition!)\n");
    fprintf(stderr, "    -tig b-e        Compute only tigs from ID 'b' to ID 'e'\n");
    fprintf(stderr, "    -u              Alias for -tig\n");
    fprintf(stderr, "    -minlength l    Do not compute consensus for tigs shorter than l bases.\n");
    fprintf(stderr, "    -maxlength l    Do not compute consensus for tigs longer than l bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -onlyunassem    Only compute consensus for unassembled tigs.\n");
    fprintf(stderr, "    -onlycontig     Only compute consensus for real unitigs/contigs.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -norepeat       Do not compute consensus for repeat tigs.\n");
    fprintf(stderr, "    -nobubble       Do not compute consensus for bubble tigs.\n");
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


  //  Open inputs.
  //
  //  We want to use corrected and trimmed reads, but definitely not
  //  homopolymer compressed reads, regardless of what the store says is the
  //  default.

  if (params.seqName) {
    fprintf(stderr, "-- Opening seqStore '%s'.\n", params.seqName);
    params.seqStore = new sqStore(params.seqName, sqStore_readOnly);

    sqRead_setDefaultVersion(sqRead_defaultVersion & ~sqRead_compressed);
    fprintf(stderr, "-- Using %s reads.\n", toString(sqRead_defaultVersion));
  }

  if (params.seqFile) {
    fprintf(stderr, "-- Using seqFile '%s'.\n", params.seqFile);

    delete params.seqStore;    //  That was a lot of work just to get sqRead_defaultVersion!
    params.seqStore = NULL;
  }

  if (params.tigName) {
    fprintf(stderr, "-- Opening tigStore '%s' version %u.\n", params.tigName, params.tigVers);
    params.tigStore = new tgStore(params.tigName, params.tigVers);

    if (params.tigEnd > params.tigStore->numTigs() - 1)
      params.tigEnd = params.tigStore->numTigs() - 1;
  }

  //  Open output files.  If we're creating a package, the usual output files are not opened.

  if ((params.exportName == NULL) && (params.outResultsName)) {
    fprintf(stderr, "-- Opening output results file '%s'.\n", params.outResultsName);
    params.outResultsFile = AS_UTL_openOutputFile(params.outResultsName);
  }

  if ((params.exportName == NULL) && (params.outLayoutsName)) {
    fprintf(stderr, "-- Opening output layouts file '%s'.\n", params.outLayoutsName);
    params.outLayoutsFile = AS_UTL_openOutputFile(params.outLayoutsName);
  }

  if ((params.exportName == NULL) && (params.outSeqNameA)) {
    fprintf(stderr, "-- Opening output FASTA file '%s'.\n", params.outSeqNameA);
    params.outSeqFileA    = AS_UTL_openOutputFile(params.outSeqNameA);
  }

  if ((params.exportName == NULL) && (params.outSeqNameQ)) {
    fprintf(stderr, "-- Opening output FASTQ file '%s'.\n", params.outSeqNameQ);
    params.outSeqFileQ    = AS_UTL_openOutputFile(params.outSeqNameQ);
  }

  //
  //  Process!
  //

  if      (params.createPartitions) {
    createPartitions(params);
  }

  else if (params.importName) {
    printHeader(params);
    processImportedTigs(params);
  }

  else if (params.exportName) {
    exportTigs(params);
  }

  else if ((params.seqFile) ||
           (params.seqName)) {
    printHeader(params);
    processTigs(params);
  }

  else {
    fprintf(stderr, "How'd you do this?  I don't know what to do.  Oops.\n");
    exit(1);
  }

  //
  //  Cleanup!
  //

  params.closeAndCleanup();

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(0);
}
