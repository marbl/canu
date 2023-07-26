
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

#include "system.H"
#include "files.H"
#include "kmers.H"
#include "strings.H"
#include "sequence.H"

using namespace merylutil::kmers::v1;

#include "sqStore.H"

#include <vector>
#include <queue>

#define BATCH_SIZE      100
#define IN_QUEUE_LENGTH 3
#define OT_QUEUE_LENGTH 3


class hapData {
public:
  hapData(char *merylname, char *histoname, char *fastaname) {
    strncpy(merylName,  merylname, FILENAME_MAX);
    strncpy(histoName,  histoname, FILENAME_MAX);
    strncpy(outputName, fastaname, FILENAME_MAX);
  }

  ~hapData() {
    delete lookup;
    delete outputWriter;
  }

public:
  void   initializeKmerTable(uint32 maxMemory);

  void   initializeOutput(void) {
    outputWriter = new compressedFileWriter(outputName);
    outputFile   = outputWriter->file();
  };

public:
  char                  merylName[FILENAME_MAX+1]  = {0};
  char                  histoName[FILENAME_MAX+1]  = {0};
  char                  outputName[FILENAME_MAX+1] = {0};

  merylExactLookup     *lookup   = nullptr;
  uint32                minCount = 0;
  uint32                maxCount = uint32max;
  uint64                nKmers   = 0;

  uint32                nReads = 0;
  uint64                nBases = 0;

  compressedFileWriter *outputWriter = nullptr;
  FILE                 *outputFile   = nullptr;
};



class allData {
public:
  allData() {
    _numThreads = getMaxThreadsAllowed();
  };

  ~allData() {
    delete _seqStore;

    for (uint32 ii=0; ii<_haps.size(); ii++)
      delete _haps[ii];

    delete _ambiguousWriter;
  };

public:
  void      openInputs(void);
  void      openOutputs(void);
  void      loadHaplotypeData(void);

public:
  char                     *_seqName   = nullptr;   //  Input from a Canu seqStore
  uint32                    _idMin     = 1;
  uint32                    _idCur     = 1;
  uint32                    _idMax     = uint32max;
  sqStore                  *_seqStore  = nullptr;
  sqRead                    _read;
  uint32                    _numReads  = 0;

  std::queue<char const *>  _inNames;
  dnaSeqFile               *_seqs      = nullptr;   //  Input from FASTA/FASTQ files.
  uint32                    _seqCounts = 0;         // read counts for current file

  std::vector<hapData *>    _haps;

  double                    _minRatio        = 1.0;
  uint32                    _minOutputLength = 1000;

  char                     *_ambiguousName   = nullptr;
  compressedFileWriter     *_ambiguousWriter = nullptr;
  FILE                     *_ambiguous       = nullptr;
  uint32                    _ambiguousReads  = 0;
  uint64                    _ambiguousBases  = 0;

  uint32                    _filteredReads   = 0;
  uint64                    _filteredBases   = 0;

  uint32                    _numThreads      = 0;   //  Reset in constructor.
  uint32                    _maxMemory       = 0;
};



class thrData {
public:
  thrData()   {                    };
  ~thrData()  { delete [] matches; };

public:
  void          clearMatches(uint32 nHaps) {
    if (matches == nullptr)
      matches = new uint32 [nHaps];

    for (uint32 hh=0; hh<nHaps; hh++)
      matches[hh] = 0;
  };

public:
  uint32       *matches = nullptr;
};



class simpleString {
public:
  simpleString()   {                 };
  ~simpleString()  { delete [] _str; };

  void    clear(void) {
    _strLen = 0;

    if (_str)
      _str[0] = 0;
  };

  //  Allocate psuedo-page size chunks of characters.
  void    set(char const *ins, uint32 insLen=0) {

    if (insLen == 0)
      insLen = strlen(ins);

    resizeArray(_str, _strLen, _strMax, 8192 * ((insLen + 8192 + 1) / 8192), _raAct::doNothing);

    assert(insLen + 1 <= _strMax);

    memcpy(_str, ins, sizeof(char) * (insLen + 1));

    _strLen = insLen;
  }

  uint32  length(void)  { return(_strLen); };
  char   *string(void)  { return(_str);    };

  uint32  _strLen = 0;
  uint32  _strMax = 0;
  char   *_str    = nullptr;
};



class readBatch {
public:
  readBatch(uint32 batchSize) {
    _numReads = 0;
    _maxReads = batchSize;

    _names = new simpleString [_maxReads];
    _bases = new simpleString [_maxReads];
    _files = new uint32       [_maxReads];
  };

  ~readBatch() {
    delete [] _names;
    delete [] _bases;
    delete [] _files;  //  Closed elsewhere!
  };

  uint32         _maxReads = 0;     //  Maximum number of reads we can store here.
  uint32         _numReads = 0;     //  Actual number of reads stored here.

  simpleString  *_names = nullptr;  //  Name of each sequence.
  simpleString  *_bases = nullptr;  //  Bases in each sequence.
  uint32        *_files = nullptr;  //  File ID where each sequence should be output.
};





uint32
getMinFreqFromHistogram(char *histoName) {

  //  If the file doesn't exist, assume it's a number and return that.

  if (fileExists(histoName) == false)
    return(strtouint32(histoName));

  //  Otherwise, open the histogram file and load into memory.  This handles
  //  both 'meryl statistics' and 'meryl histogram' outputs.

  splitToWords  S;

  uint32        Llen     = 0;
  uint32        Lmax     = 0;
  char         *L        = new char [Lmax];

  uint32        histoLen = 0;
  uint32        histoMax = 1024;
  uint64       *histo    = new uint64 [histoMax];

  FILE *H = merylutil::openInputFile(histoName);
  while (merylutil::readLine(L, Llen, Lmax, H)) {
    S.split(L);

    uint32  f = 0;   //  Frequency
    uint64  v = 0;   //  Number of kmers at that frequency.

    if (S.numWords() >= 2) {  //  First word must be the frequency,
      f = S.touint32(0);     //  second word must be the number of
      v = S.touint64(1);     //  kmer that occur exactly f times.
    }

    if (f == 0)              //  If zero, assume it's a header line
      continue;              //  or some other such crud.

    if (f >= histoMax)       //  If big, we're done, we don't
      break;                 //  care about highly common kmers.

    histo[f] = v;
    histoLen = f + 1;
  }
  merylutil::closeFile(H, histoName);

  delete [] L;

  //  Use the histogram to decide on min and max frequency thresholds.
  //  Pick the frequency:
  //    For min, at the bottom of the trough between the noise and real.
  //    For max, after the peak, and that has count 75% of the min.
  //
  //  *
  //   *       **
  //    *     *  *
  //     *   *    *
  //      ***      *  +-max
  //       ^        * v
  //       +-min     ***
  //                    **
  //
  //  Both are picked using a 9-window rolling average.
  //
  //  The maximum cutoff is disabled, based on limited testing and the
  //  argument that even if they are repeat kmers, they are still unique to
  //  the haplotype, so are still useful.

  uint32  aveSize = 5;
  uint64  thisSum = 0;
  uint32  thisLen = 0;

  uint32  f = 1;

  uint32  minFreq = 1;
  uint64  minAve  = histo[f++];

  //for (uint32 ii=0; ii<aveSize; ii++) {
  //  thisSum += histo[f++];
  //  thisLen += 1;
  //}

  //  Scan forward until the rolling average starts to increase, declare the middle
  //  of the average the minimum.  Keep searching ahead until our current average is
  //  twice that of the minimum found.

  for (; f < histoLen; f++) {
    if (thisLen == aveSize) {
      thisSum += histo[f];
      thisSum -= histo[f - aveSize];
    }

    else {
      thisSum += histo[f];
      thisLen++;
    }

    if (thisSum / thisLen < minAve) {
      minFreq = f - thisLen/2;
      minAve  = thisSum / thisLen;
    }

    if (2 * minAve * aveSize < thisSum)   //  Over estimates the minimum sum when thisLen < aveSize - i.e., for
      break;                              //  frequencies 1, 2, 3 and 4.  Probably not an issue.
  }

  delete [] histo;

  return(minFreq);
}






void
hapData::initializeKmerTable(uint32 maxMemory) {

  //  Decide on a threshold below which we consider the kmers as useless noise.

  uint32 minFreq = getMinFreqFromHistogram(histoName);

  fprintf(stdout, "--  Haplotype '%s':\n", merylName);
  fprintf(stdout, "--   use kmers with frequency at least %u.\n", minFreq);

  //  Construct an exact lookup table.
  //
  //  If there is not valid merylName, do not load data.  This is only useful
  //  for testing getMinFreqFromHistogram() above.
  //
  //  Get this behavior with option '-H "" histo out.fasta',

  if (merylName[0]) {
    merylFileReader  *reader = new merylFileReader(merylName);

    lookup = new merylExactLookup();
    lookup->load(reader, maxMemory, 0, minFreq, UINT32_MAX);
    
    nKmers = lookup->nKmers();

    delete reader;
  }

  //  And report what we loaded.

  fprintf(stderr, "--   loaded %lu kmers.\n", nKmers);
};



//  Open inputs and check the range of reads to operate on.
void
allData::openInputs(void) {

  //  Open sqStore, remember the number of reads.

  if (_seqName) {
    _seqStore = new sqStore(_seqName);
    _numReads = _seqStore->sqStore_lastReadID();

    if (_numReads < _idMax)
      _idMax = _numReads;
  }

  //  Nothing to do for dnaSeqFiles, they're opened when the command
  //  line is parsed.
}



//  Open outputs for each haplotype, and one for ambiguous reads.
void
allData::openOutputs(void) {

  for (uint32 ii=0; ii<_haps.size(); ii++)
    _haps[ii]->initializeOutput();

  if (_ambiguousName) {
    _ambiguousWriter = new compressedFileWriter(_ambiguousName);
    _ambiguous       = _ambiguousWriter->file();
  }
}



//  Create meryl exact lookup structures for all the haplotypes.
void
allData::loadHaplotypeData(void) {
  uint32 memPerHap = _maxMemory / _haps.size();

  if (memPerHap == 0)   //  If zero, it would be allowed
    memPerHap = 1;      //  to use all available memory!

  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Loading haplotype data, using up to %u GB memory for each.\n", memPerHap);
  fprintf(stderr, "--\n");

  for (uint32 ii=0; ii<_haps.size(); ii++)
    _haps[ii]->initializeKmerTable(memPerHap);

  fprintf(stderr, "-- Data loaded.\n");
  fprintf(stderr, "--\n");
}







void *
loadReadBatch(void *G) {
  allData     *g = (allData   *)G;
  readBatch   *s = NULL;

  s = new readBatch(BATCH_SIZE);   //  We should be using recycled ones.
  //fprintf(stderr, "Alloc  readBatch s %p\n", s);

  s->_numReads = 0;

  dnaSeq  seq;

  while (s->_numReads < s->_maxReads) {
    uint32 rr = s->_numReads;   //  Where to put the read we're loading.

    //  Try to load a sequence from the seqStore.

    if ((g->_seqStore) &&
        (g->_idCur <= g->_idMax)) {
      uint32  readLen = g->_seqStore->sqStore_getReadLength(g->_idCur);

      if (readLen >= g->_minOutputLength) {
        g->_seqStore->sqStore_getRead(g->_idCur, &g->_read);

        s->_names[rr].set(g->_read.sqRead_name());
        s->_bases[rr].set(g->_read.sqRead_sequence(), readLen);
        s->_files[rr] = UINT32_MAX;

        s->_numReads++;
      }
      else {
        g->_filteredReads++;
        g->_filteredBases += readLen;
      }

      g->_idCur++;   //  Loaded a sequence!  (or tried to and failed)
      continue;      //  Continue on to loading the next one.
    }

    //  Nope, try to load from one of the sequence files.

    if (g->_inNames.empty() == false) {
      if (g->_seqs == nullptr) {
        fprintf(stdout, "-- Begin    processing file %s\n", g->_inNames.front());
        g->_seqs = new dnaSeqFile(g->_inNames.front());
      }

      if (g->_seqs->loadSequence(seq) == false) {   //  Failed to load a sequence, hit EOF.
        if (g->_seqCounts == 0) {
            fprintf(stdout, "--\n");
            fprintf(stdout, "-- ERROR:   loaded no reads from file %s, are you sure it is a valid fastq/fasta file?\n", g->_inNames.front());
            fprintf(stdout, "--\n");
            fflush(stdout);
            exit(1);
        }

        fprintf(stdout, "-- Finished processing file %s with %d records\n", g->_inNames.front(), g->_seqCounts);
        fprintf(stdout, "--\n");

        delete g->_seqs;                            //  Discard the file and try the next.

        g->_seqs = nullptr;
        g->_seqCounts = 0;
        g->_inNames.pop();

        continue;
      }

      if (seq.length() >= g->_minOutputLength) {            //  Loaded something.  If it's long
        s->_names[rr].set(seq.ident());                     //  enough, save it to our list.
        s->_bases[rr].set(seq.bases(), seq.length());
        s->_files[rr] = UINT32_MAX;

        s->_numReads++;
      }
      else {
        g->_filteredReads++;
        g->_filteredBases += seq.length();
      }

      g->_seqCounts++;
      continue;      //  Loaded (or skipped) a sequence.  Thank you, may I have another?
    }

    //  Nope, guess there's no more reads.

    break;
  }

  if (s->_numReads == 0) {
    delete s;
    s = NULL;
  }

  //fprintf(stderr, "Return readBatch s %p with %u/%u reads %p %p %p\n", s, s->_numReads, s->_maxReads, s->_names, s->_bases, s->_files);

  return(s);
};



void
processReadBatch(void *G, void *T, void *S) {
  allData     *g = (allData   *)G;
  thrData     *t = (thrData   *)T;
  readBatch   *s = (readBatch *)S;

  //fprintf(stderr, "Proces readBatch s %p with %u/%u reads %p %p %p\n", s, s->_numReads, s->_maxReads, s->_names, s->_bases, s->_files);

  uint32       nHaps   = g->_haps.size();
  uint32      *matches = new uint32 [nHaps];

  for (uint32 ii=0; ii<s->_numReads; ii++) {

    //  Count the number of matching kmers for each haplotype.
    //
    //  The kmer iteration came from merylOp-count.C and merylOp-countSimple.C.

    for (uint32 hh=0; hh<nHaps; hh++)
      matches[hh] = 0;

    kmerIterator  kiter(s->_bases[ii].string(),
                        s->_bases[ii].length());

    while (kiter.nextMer())
      for (uint32 hh=0; hh<nHaps; hh++)
        if ((g->_haps[hh]->lookup->value(kiter.fmer()) > 0) ||
            (g->_haps[hh]->lookup->value(kiter.rmer()) > 0))
          matches[hh]++;

    //  Find the haplotype with the most and second most matching kmers.

    uint32  hap1st = UINT32_MAX;   //  Index of the best haplotype.
    double  sco1st = 0.0;          //  Score of the best haplotype.

    uint32  hap2nd = UINT32_MAX;   //  Index of the second best haplotype.
    double  sco2nd = 0.0;          //  Score of the second best haplotype.

    for (uint32 hh=0; hh<nHaps; hh++) {
      double  sco = (double)matches[hh] / g->_haps[hh]->nKmers;

      if      (sco1st <= sco) {
        hap2nd = hap1st;    sco2nd = sco1st;
        hap1st = hh;        sco1st = sco;
      }
      else if (sco2nd <= sco) {
        hap2nd = hh;        sco2nd = sco;
      }

      assert(sco2nd <= sco1st);
    }

    if (0) {
      if (((sco2nd < DBL_MIN) && (sco1st > DBL_MIN)) ||
          ((sco2nd > DBL_MIN) && (sco1st / sco2nd > g->_minRatio)))
        fprintf(stdout, "hap1st %1u sco1st %9.7f matches %6u   hap2 %1u sco2nd %9.7f matches %6u -> %1u\n",
                hap1st, sco1st, matches[hap1st],
                hap2nd, sco2nd, matches[hap2nd], hap1st);
      else
        fprintf(stdout, "hap1st %1u sco1st %9.7f matches %6u   hap2 %1u sco2nd %9.7f matches %6u -> AMBIGUOUS\n",
                hap1st, sco1st, matches[hap1st],
                hap2nd, sco2nd, matches[hap2nd]);
    }

    //  Write the read to the 'best' haplotype, unless it's an ambiguous assignment.
    //
    //  By default, write to the ambiguous file.
    //
    //  Write to the best file only if
    //   - there is a non-zero best score and the second best is zero
    //   - the ratio of best to second best is bigger than some threshold

    s->_files[ii] = UINT32_MAX;

    if (((sco2nd < DBL_MIN) && (sco1st > DBL_MIN)) ||
        ((sco2nd > DBL_MIN) && (sco1st / sco2nd > g->_minRatio)))
      s->_files[ii] = hap1st;
  }

  delete [] matches;
}


void
outputReadBatch(void *G, void *S) {
  allData     *g = (allData   *)G;
  readBatch   *s = (readBatch *)S;
  FILE        *F = NULL;

  for (uint32 ii=0; ii<s->_numReads; ii++) {
    uint32 ff = s->_files[ii];

    if (ff == UINT32_MAX) {
      F = g->_ambiguous;

      g->_ambiguousReads += 1;
      g->_ambiguousBases += s->_bases[ii].length();

    } else {
      F = g->_haps[ff]->outputFile;

      g->_haps[ff]->nReads += 1;
      g->_haps[ff]->nBases += s->_bases[ii].length();
    }

    outputFASTA(F, s->_bases[ii].string(), s->_bases[ii].length(), 0,
                "%s", s->_names[ii].string());
  }

  delete s;    //  We should recycle this, but hard to do.
}






int
main(int argc, char **argv) {
  allData      *G          = new allData;
  bool          beVerbose  = false;

  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-S") == 0) {   //  INPUT READS and RANGE TO PROCESS
      G->_seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], G->_idMin, G->_idMax);

    } else if (strcmp(argv[arg], "-R") == 0) {
      while ((arg < argc) && (fileExists(argv[arg+1])))
        G->_inNames.push(argv[++arg]);

    } else if (strcmp(argv[arg], "-H") == 0) {   //  HAPLOTYPE SPECIFICATION
      G->_haps.push_back(new hapData(argv[arg+1], argv[arg+2], argv[arg+3]));
      arg += 3;

    } else if (strcmp(argv[arg], "-A") == 0) {
      G->_ambiguousName = argv[++arg];

    } else if (strcmp(argv[arg], "-cr") == 0) {  //  PARAMETERS
      G->_minRatio = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-cl") == 0) {
      G->_minOutputLength = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      G->_numThreads = setNumThreads(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory") == 0) {
      G->_maxMemory  = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;


    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  if ((G->_seqName == NULL) && (G->_inNames.size() == 0))
    err.push_back("No input sequences supplied with either (-S) or (-R).\n");
  if ((G->_seqName != NULL) && (G->_inNames.size() != 0))
    err.push_back("Only one type of input reads (-S or -R) supported.\n");
  if (G->_haps.size() < 2)
    err.push_back("Not enough haplotypes (-H) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "READ INPUTS:\n");
    fprintf(stderr, "  Expects PacBio or Nanopore reads in one or more Canu seqStore, FASTA or FASTQ\n");
    fprintf(stderr, "  inputs.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S seqStore                      path to input seqStore of reads to classify.\n");
    fprintf(stderr, "  -r bgn[-end]                     range of reads to operate on.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -R reads.fasta                   path to input FASTA or FASTQ of reads to classify.\n");
    fprintf(stderr, "                                   these may be uncompressed, gzip, bzip2 or xz compressed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "HAPLOTYPE INPUTS AND OUTPUTS\n");
    fprintf(stderr, "  Each -H option specifies a haplotype using three parameters.\n");
    fprintf(stderr, "    haplo-kmers.meryl       - haplotype specific kmers contained in a meryl database.\n");
    fprintf(stderr, "    parent-kmers.histogram  - a histogram of all parent kmers.\n");
    fprintf(stderr, "    haplo-output.fasta.gz   - output reads assigned to this haplotype.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -H haplo-kmers.meryl parent-kmers.histogram haplo-output.fasta.gz\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  The 'parent-kmers.histgram' is used to determine a noise threshold.  kmers\n");
    fprintf(stderr, "  that occur fewer than that many times are ignored as being likely noise kmers.\n");
    fprintf(stderr, "  Instead of a full histogram, a single integer can be supplied to directly\n");
    fprintf(stderr, "  set the threshold.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUTS:\n");
    fprintf(stderr, "  Haplotype-specific reads are written to 'haplo.fasta.gz' as specified in each -H\n");
    fprintf(stderr, "  option.  Reads not assigned to any haplotype are written to the file specified\n");
    fprintf(stderr, "  in the -A option, if not specified, they are silently discarded.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Output fasta files are 'gzip -1' compressed if they end in '.gz'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -A ambiguous.fasta.gz\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "PARAMETERS\n");
    fprintf(stderr, "  -cr ratio        minimum ratio between best and second best to classify\n");
    fprintf(stderr, "  -cl length       minimum length of output read\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -v               report how many batches per second are being processed\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  G->openInputs();
  G->openOutputs();

  G->loadHaplotypeData();

  thrData   *TD = new thrData [G->_numThreads];
  sweatShop *SS = new sweatShop(loadReadBatch, processReadBatch, outputReadBatch);

  SS->setNumberOfWorkers(G->_numThreads);

  for (uint32 ii=0; ii<G->_numThreads; ii++)
    SS->setThreadData(ii, TD + ii);

  SS->setLoaderBatchSize(1);
  SS->setLoaderQueueSize(G->_numThreads * IN_QUEUE_LENGTH);
  SS->setWorkerBatchSize(1);
  SS->setWriterQueueSize(G->_numThreads * OT_QUEUE_LENGTH);

  fprintf(stderr, "-- Processing reads in batches of %u reads each.\n", BATCH_SIZE);
  fprintf(stderr, "--\n");

  SS->run(G, beVerbose);

  //  Write the log to stdout.

  for (uint32 ii=0; ii<G->_haps.size(); ii++)
    fprintf(stdout, "-- %8u reads %12lu bases written to haplotype file %s.\n", G->_haps[ii]->nReads, G->_haps[ii]->nBases, G->_haps[ii]->outputName);
  fprintf(stdout, "-- %8u reads %12lu bases written to haplotype file %s.\n", G->_ambiguousReads, G->_ambiguousBases, G->_ambiguousName);
  fprintf(stdout, "--\n");
  fprintf(stdout, "-- %8u reads %12lu bases filtered for being too short.\n", G->_filteredReads, G->_filteredBases);
  fprintf(stdout, "--\n");

  delete    SS;
  delete [] TD;
  delete    G;

  fprintf(stderr, "-- Bye.\n");
  exit(0);
}
