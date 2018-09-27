
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-08
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "sqStore.H"

#include "files.H"
#include "kmers.H"
#include "strings.H"
#include "sequence.H"

#include "sweatShop.H"

#include <vector>
#include <queue>

using namespace std;


#define BATCH_SIZE      1024
#define IN_QUEUE_LENGTH 3
#define OT_QUEUE_LENGTH 3


class hapData {
public:
  hapData(char *name);
  ~hapData();

public:
  void   initializeKmerTable(uint32 minFrequency, uint32 maxFrequency);

  void   initializeOutput(void) {
    outputWriter = new compressedFileWriter(outputName);
    outputFile   = outputWriter->file();
  };

public:
  char                    merylName[FILENAME_MAX+1];
  char                    outputName[FILENAME_MAX+1];

  kmerCountExactLookup   *lookup;
  uint32                  minCount;
  uint32                  maxCount;
  uint64                  nKmers;

  compressedFileWriter   *outputWriter;
  FILE                   *outputFile;
};



class allData {
public:
  allData() {
    _seqName          = NULL;
    _idMin            = 1;
    _idCur            = 1;
    _idMax            = UINT32_MAX;
    _seqStore         = NULL;
    _numReads         = 0;

    //  _seqs and _haps are assumed to be clear already.

    _minF             = 0;
    _maxF             = UINT32_MAX;

    _minRatio         = 1.0;
    _minOutputLength  = 1000;

    _ambiguousName   = NULL;
    _ambiguousWriter = NULL;
    _ambiguous       = NULL;
  };

  ~allData() {

    if (_seqStore)
      _seqStore->sqStore_close();

    for (uint32 ii=0; ii<_haps.size(); ii++)
      delete _haps[ii];

    delete _ambiguousWriter;
  };

public:
  void      openInputs(void);
  void      openOutputs(void);
  void      loadHaplotypeData(void);
  

public:
  char                  *_seqName;   //  Input from a Canu seqStore
  uint32                 _idMin;
  uint32                 _idCur;
  uint32                 _idMax;
  sqStore               *_seqStore;
  sqReadData             _readData;
  uint32                 _numReads;

  queue<dnaSeqFile *>    _seqs;      //  Input from FASTA/FASTQ files.

  vector<hapData *>      _haps;

  uint32                 _minF;      //  Not actually set by
  uint32                 _maxF;      //  command line parameters

  double                 _minRatio;
  uint32                 _minOutputLength;

  char                  *_ambiguousName;
  compressedFileWriter  *_ambiguousWriter;
  FILE                  *_ambiguous;
};



class thrData {
public:
  thrData() {
    matches = NULL;
  };

  ~thrData() {
    delete [] matches;
  };

public:
  void          clearMatches(uint32 nHaps) {
    if (matches == NULL)
      matches = new uint32 [nHaps];

    for (uint32 hh=0; hh<nHaps; hh++)
      matches[hh] = 0;
  };


public:
  uint32       *matches;
};



class simpleString {
public:
  simpleString() {
    _strLen = 0;
    _strMax = 0;
    _str    = NULL;
  };
  ~simpleString() {
    delete [] _str;
  };

  void    clear(void) {
    _strLen = 0;

    if (_str)
      _str[0] = 0;
  };

  //  Allocate psuedo-page size chunks of characters.
  void    set(char *ins, uint32 insLen=0) {

    if (insLen == 0)
      insLen = strlen(ins);

    resizeArray(_str, _strLen, _strMax, 8192 * ((insLen + 8192 + 1) / 8192), resizeArray_doNothing);

    assert(insLen + 1 <= _strMax);

    memcpy(_str, ins, sizeof(char) * (insLen + 1));

    _strLen = insLen;
  }

  uint32  length(void) {
    return(_strLen);
  };

  char   *string(void)  {
    return(_str);
  };

  uint32  _strLen;
  uint32  _strMax;
  char   *_str;
};



class readBatch {
public:
  readBatch(uint32 batchSize) {
    _numReads = 0;
    _maxReads = batchSize;

    _names = new simpleString [_maxReads];
    _bases = new simpleString [_maxReads];
    _files = new FILE *       [_maxReads];

    //fprintf(stderr, "Alloc names %p\n", _names);
    //fprintf(stderr, "Alloc bases %p\n", _bases);
    //fprintf(stderr, "Alloc files %p\n", _files);
  };

  ~readBatch() {
    //fprintf(stderr, "Delet names %p\n", _names);
    //fprintf(stderr, "Delet bases %p\n", _bases);
    //fprintf(stderr, "Delet files %p\n", _files);

    delete [] _names;
    delete [] _bases;
    delete [] _files;  //  Closed elsewhere!
  };

  uint32         _maxReads;    //  Maximum number of reads we can store here.
  uint32         _numReads;    //  Actual number of reads stored here.

  simpleString  *_names;       //  Name of each sequence.
  simpleString  *_bases;       //  Bases in each sequence.
  FILE         **_files;       //  File where each sequence should be output.
};










hapData::hapData(char *name) {
  strncpy(merylName,  name, FILENAME_MAX);
  strncpy(outputName, name, FILENAME_MAX);

  uint32 outputNameLen = strlen(outputName);

  if ((outputNameLen > 6) &&
      (strcmp(outputName + outputNameLen - 6, ".meryl") == 0))
    outputName[outputNameLen-6] = 0;

  strcat(outputName, ".fasta.gz");

  lookup       = NULL;
  minCount     = 0;
  maxCount     = UINT32_MAX;
  nKmers       = 0;

  outputWriter = NULL;
  outputFile   = NULL;
};


hapData::~hapData() {
  delete lookup;
  delete outputWriter;
};


void
hapData::initializeKmerTable(uint32 minFrequency, uint32 maxFrequency) {
  kmerCountFileReader  *reader = new kmerCountFileReader(merylName);
  kmerCountStatistics  *stats  = reader->stats();

  //
  //  Use the reader histogram to decide on min and max frequency thresholds.
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

  uint32  aveSize = 9;
  uint64  thisSum = 0;
  uint32  f = 2;

  for (uint32 ii=0; ii<aveSize; ii++)
    thisSum += stats->numKmersAtFrequency(f++);

  uint32  minFreq = f-1 - aveSize/2, maxFreq = UINT32_MAX;
  uint64  minSum  = thisSum,         maxSum  = 0;

  //  Scan forward until the rolling average starts to increase, declare the middle
  //  of the average the minimum.  Keep searching ahead until our current average is
  //  twice that of the minimum found.

  for (; f < stats->numFrequencies(); f++) {
    thisSum = (stats->numKmersAtFrequency(f) + thisSum - stats->numKmersAtFrequency(f - aveSize));

    if (thisSum < minSum) {
      minFreq = f - aveSize/2;
      minSum  = thisSum;
    }

    if (2 * minSum < thisSum)
      break;
  }

  //  Continue scanning until we find the number of kmers falls below 0.75 * min.
#if 0
  for (; f < stats->numFrequencies(); f++) {
    thisSum = (stats->numKmersAtFrequency(f) + thisSum - stats->numKmersAtFrequency(f - aveSize));

    if (thisSum < 0.75 * minSum) {
      maxFreq = f - aveSize/2;
      maxSum  = thisSum;

      break;
    }
  }
#endif

  fprintf(stderr, "Use min freq %u sum %lu -- max freq %u sum %lu\n",
          minFreq, minSum,
          maxFreq, maxSum);

  //
  //  With those set, construct an exact lookup table.
  //

  lookup = new kmerCountExactLookup(reader, minFreq, maxFreq);
  nKmers = lookup->nKmers();

  //

  delete reader;
};



//  Open inputs and check the range of reads to operate on.
void
allData::openInputs(void) {

  //  Open sqStore, remember the number of reads.

  if (_seqName) {
    _seqStore = sqStore::sqStore_open(_seqName);
    _numReads = _seqStore->sqStore_getNumReads();

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

  for (uint32 ii=0; ii<_haps.size(); ii++)
    _haps[ii]->initializeKmerTable(_minF, _maxF);
}







void *
loadReadBatch(void *G) {
  allData     *g = (allData   *)G;
  readBatch   *s = NULL;

  s = new readBatch(BATCH_SIZE);   //  We should be using recycled ones.
  //fprintf(stderr, "Alloc  readBatch s %p\n", s);

  //char  *name  = new char [AS_MAX_READLEN + 1];
  //char  *bases = new char [AS_MAX_READLEN + 1];
  //char  *quals = new char [AS_MAX_READLEN + 1];

  s->_numReads = 0;

  dnaSeq  seq;

  while (s->_numReads < s->_maxReads) {
    uint32 rr = s->_numReads;   //  Where to put the read we're loading.

    //  Try to load a sequence from the seqStore.

    if ((g->_seqStore) && 
        (g->_idCur <= g->_idMax)) {
      uint32  readLen = g->_seqStore->sqStore_getRead(g->_idCur)->sqRead_sequenceLength();

      if (readLen >= g->_minOutputLength) {
        g->_seqStore->sqStore_loadReadData(g->_idCur, &g->_readData);

        s->_names[rr].set(g->_readData.sqReadData_getName());
        s->_bases[rr].set(g->_readData.sqReadData_getSequence(), readLen);
        s->_files[rr] = NULL;

        s->_numReads++;
      }

      g->_idCur++;   //  Loaded a sequence!  (or tried to and failed)
      continue;      //  Continue on to loading the next one.
    }

    //  Nope, try to load from one of the sequence files.

    if (g->_seqs.empty() == false) {
      if (g->_seqs.front()->loadSequence(seq) == false) {   //  Failed to load a sequence, hit EOF.
        delete g->_seqs.front();                            //  Discard the file and try the next.
        g->_seqs.pop();
        continue;
      }

      if (seq.length() >= g->_minOutputLength) {
        s->_names[rr].set(seq.name());
        s->_bases[rr].set(seq.bases(), seq.length());
        s->_files[rr] = NULL;

        s->_numReads++;
      }

      continue;      //  Loaded (or skipped) a sequence.  Thank you, may I have another?
    }

    //  Nope, guess there's no more reads.

    break;
  }

  if (s->_numReads == 0) {
    delete s;
    s = NULL;
  }

  //delete [] name;
  //delete [] bases;
  //delete [] quals;

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

    kmerTiny   fmer;    //  merSize() is set when the meryl databases are loaded.
    kmerTiny   rmer;

    uint32     kmerLoad   = 0;
    uint32     kmerValid  = fmer.merSize() - 1;
    uint32     kmerSize   = fmer.merSize();

    char      *name       = s->_names[ii].string();
    char      *bases      = s->_bases[ii].string();
    uint32     basesLen   = s->_bases[ii].length();

    for (uint32 ss=0; ss<basesLen; ss++) {
      if ((bases[ss] != 'A') && (bases[ss] != 'a') &&   //  If not valid DNA, don't
          (bases[ss] != 'C') && (bases[ss] != 'c') &&   //  make a kmer, and reset
          (bases[ss] != 'G') && (bases[ss] != 'g') &&   //  the count until the next
          (bases[ss] != 'T') && (bases[ss] != 't')) {   //  valid kmer is available.
        kmerLoad = 0;
        continue;
      }

      fmer.addR(bases[ss]);
      rmer.addL(bases[ss]);

      if (kmerLoad < kmerValid) {   //  If not a full kmer, increase the length we've
        kmerLoad++;                 //  got loaded, and keep going.
        continue;
      }

      for (uint32 hh=0; hh<nHaps; hh++) {
        if ((g->_haps[hh]->lookup->value(fmer) > 0) ||
            (g->_haps[hh]->lookup->value(rmer) > 0))
          matches[hh]++;
      }
    }

    //  Find the haplotype with the most and second most matching kmers.

    uint32  hap1 = UINT32_MAX,   hap2 = UINT32_MAX;
    double  sco1 = 0.0,          sco2 = 0.0;

    for (uint32 hh=0; hh<nHaps; hh++) {
      double  sco = (double)matches[hh] / g->_haps[hh]->nKmers;

      if      (sco1 <= sco) {
        hap2 = hap1;    sco2 = sco1;
        hap1 = hh;      sco1 = sco;
      }
      else if (sco2 <= sco) {
        hap2 = hh;      sco2 = sco;
      }

      assert(sco2 <= sco1);
    }

#if 0
    if (((sco2 < DBL_MIN) && (sco1 > DBL_MIN)) ||
        ((sco2 > DBL_MIN) && (sco1 / sco2 > minRatio)))
      fprintf(stdout, "hap1 %1u sco1 %9.7f matches %6u   hap2 %1u sco2 %9.7f matches %6u -> %1u\n",
              hap1, sco1, matches[hap1],
              hap2, sco2, matches[hap2], hap1);
    else
      fprintf(stdout, "hap1 %1u sco1 %9.7f matches %6u   hap2 %1u sco2 %9.7f matches %6u -> AMBIGUOUS\n",
              hap1, sco1, matches[hap1],
              hap2, sco2, matches[hap2]);
#endif

    //  Write the read to the 'best' haplotype, unless it's an ambiguous assignment.
    //
    //  By default, write to the ambiguous file.
    //
    //  Write to the best file only if
    //   - there is a non-zero best score and the second best is zero
    //   - the ratio of best to second best is bigger than some threshold
     
    s->_files[ii] = g->_ambiguous;

    if (((sco2 < DBL_MIN) && (sco1 > DBL_MIN)) ||
        ((sco2 > DBL_MIN) && (sco1 / sco2 > g->_minRatio)))
      s->_files[ii] = g->_haps[hap1]->outputFile;
  }

  delete [] matches;
}


void
outputReadBatch(void *G, void *S) {
  allData     *g = (allData   *)G;
  readBatch   *s = (readBatch *)S;

  for (uint32 ii=0; ii<s->_numReads; ii++) {
    if (s->_files[ii])
      AS_UTL_writeFastA(s->_files[ii],
                        s->_bases[ii].string(), s->_bases[ii].length(), 0,
                        ">%s\n", s->_names[ii].string());
  }

  delete s;    //  We should recycle this, but hard to do.
}






int
main(int argc, char **argv) {
  allData      *G          = new allData;
  uint32        numThreads = 1;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {   //  INPUT READS and RANGE TO PROCESS
      G->_seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], G->_idMin, G->_idMax);

    } else if (strcmp(argv[arg], "-R") == 0) {
      while ((arg < argc) && (fileExists(argv[arg+1])))
        G->_seqs.push(new dnaSeqFile(argv[++arg]));

    } else if (strcmp(argv[arg], "-H") == 0) {   //  HAPLOTYPE SPECIFICATION
      G->_haps.push_back(new hapData(argv[++arg]));

    } else if (strcmp(argv[arg], "-A") == 0) {
      G->_ambiguousName = argv[++arg];

    } else if (strcmp(argv[arg], "-cr") == 0) {  //  PARAMETERS
      G->_minRatio = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-cl") == 0) {
      G->_minOutputLength = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      numThreads = strtouint32(argv[++arg]);


    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((G->_seqName == NULL) && (G->_seqs.size() == 0))
    err.push_back("No input sequences supplied with either (-S) or (-R).\n");
  if ((G->_seqName != NULL) && (G->_seqs.size() != 0))
    err.push_back("Only one type of input reads (-S or -R) supported.\n");
  if (G->_haps.size() < 2)
    err.push_back("Not enough haplotypes (-H) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS:\n");
    fprintf(stderr, "  Expects PacBio or Nanopore reads in one or more Canu seqStore, FASTA or FASTQ\n");
    fprintf(stderr, "  inputs.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Expects meryl kmer databases for two or more haplotypes.  The databases should\n");
    fprintf(stderr, "  contain kmers that are only present in that haplotype.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S seqStore      path to input seqStore of reads to classify.\n");
    fprintf(stderr, "  -r bgn[-end]     range of reads to operate on.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -R seq.fasta     path to input FASTA or FASTQ of reads to classify.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -H haplo.meryl   path to input haplotype-specific kmer database.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUTS:\n");
    fprintf(stderr, "  Haplotype-specific reads are written to 'haplo.fasta.gz' corresponding to each\n");
    fprintf(stderr, "  'haplo.meryl' listed in -H options.  Reads not assigned to any haplotype are\n");
    fprintf(stderr, "  written to the file specified with the -A option, if not specified, they are\n");
    fprintf(stderr, "  silently discarded.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Output fasta files are 'gzip -1' compressed if they end in '.gz'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -A ambiguous.fasta.gz\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "PARAMETERS\n");
    fprintf(stderr, "  -cr ratio        minimum ratio between best and second best to classify\n");
    fprintf(stderr, "  -cl length       minimum length of output read\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  omp_set_num_threads(numThreads);  //  Lets the kmer data be loaded with threads.

  G->openInputs();
  G->openOutputs();

  G->loadHaplotypeData();
  
  thrData   *TD = new thrData [numThreads];
  sweatShop *SS = new sweatShop(loadReadBatch, processReadBatch, outputReadBatch);

  SS->setNumberOfWorkers(numThreads);

  for (uint32 ii=0; ii<numThreads; ii++)
    SS->setThreadData(ii, TD + ii);

  SS->setLoaderBatchSize(1);
  SS->setLoaderQueueSize(numThreads * IN_QUEUE_LENGTH);
  SS->setWorkerBatchSize(1);
  SS->setWriterQueueSize(numThreads * OT_QUEUE_LENGTH);

  SS->run(G, false);

  delete    SS;
  delete [] TD;
  delete    G;

  fprintf(stderr, "Bye.\n");
  exit(0);
}
