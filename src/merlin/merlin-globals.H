
/******************************************************************************
 *
 *  This is a k-mer based variant evaluation tool for polishing assemblies.
 *
 *  This software is based on:
 *    'Meryl'                  (https://github.com/marbl/meryl)
 *
 *  This is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */


#include "files.H"
#include "kmers.H"
#include "strings.H"
#include "sequence.H"
#include "system.H"
using namespace merylutil::kmers::v1;
#include "sqStore.H"

#include <string>
#include <vector>

#include "graph.H"

#ifndef MERFIN_GLOBAL
#define MERFIN_GLOBAL

#define OP_NONE       0
#define OP_BUILD      1

////////////////////////////////////////
//
//  This object holds the input sequence and output results.  A new object is
//  created for each input, passed to the compute function, then passed to
//  the output function and deleted there.
//
class merlinInput {
public:
public:
  merlinInput() {
  };
  ~merlinInput() {
  };

  //  Input data.

  dnaSeq         seq;
  kmerIterator   kiter;

  //  Some results.
  uint64         totalK = 0;   //  number of kmers in the sequence processed

  //  merlin-graph
  vector<kmer>   mersFound;   //  kmer
  vector<uint64> mersPos;     //  position
  vector<char>   mersOrient;  //  orientation
  vector<kmvalu> mersValue;

};



////////////////////////////////////////
//
//  This object holds per-thread data for -vmer work.
//   - a handful of arrays that we can recycle for each input
//   - a per-thread debug logging file
//
class
merlinThrData {
public:
  merlinThrData() {
  }
  ~merlinThrData() {
    delete oDebug;
  }

  //  Our ident.

  uint32                  threadID;

  //  For debug logging.

  compressedFileWriter   *oDebug   = nullptr;
};



////////////////////////////////////////
//
//  This object holds command line parameters and input data.  The input data
//  should be treated as read-only; multiple threads will access it with no
//  locking.
//
//  This object also holds any summary data for output.  This data should be
//  modified only in the sweatshop 'output' function - no locking is needed as
//  the sweatshop guarantees only one thread is active at a time here.
//
class merlinGlobal {
public:
  merlinGlobal(char *const execname) {
    execName = execname;
  }

  ~merlinGlobal() {
    //  input seq handle
    delete seqFile;

    //  lookup db
    delete markerLookup;
    delete pLookupTable;

    //  output handle
    delete graphGfaFile;
  }

public:
  void      load_Kmetric(void); // we don't need this yet...
  void      load_Sequence(void);
  void      load_Kmers(char*);

  void      open_Inputs(void);
  kmvalu    lookupMarker(kmer     fmer,
                         kmer     rmer);

  void      updateGraph(void *G, void *S);
  void      outputGraph(void);

  void      colorGraph(const char*);
  void      outputCsv(void);

  //  Parameters for build mode.
public:
  compressedFileWriter   *graphGfaFile  = nullptr;
  compressedFileWriter   *graphCsvFile  = nullptr;

  //  Input data.
public:
  char                   *seqName       = nullptr;
  dnaSeqFile             *seqFile       = nullptr;

  char                   *outName       = nullptr;

  char                   *markerDBname  = nullptr;
  char                   *hapADBname    = nullptr;
  char                   *hapBDBname    = nullptr;
  merylExactLookup       *markerLookup;

  uint64                  minV          = 0;           //  Kmer value limits for the
  uint64                  maxV          = uint64max;   //  above two lookup tables.

  char                   *pLookupTable  = nullptr;     //  Probability lookup table
  vector<uint32>          copyKmerK;
  vector<double>          copyKmerP;

  //  Run-time parameters.
public:
  uint32                  threads       = getMaxThreadsAllowed();
  double                  maxMemory     = getMaxMemoryAllowed() / 1024.0 / 1024.0 / 1024.0;

  uint32                  reportType    = OP_NONE;

  bool                    debug         = false;

public:
  double                  peak          = 0;

  char                   *execName      = nullptr;

  //  Graph
public:
  map<kmer, node*>        nodes;
  vector<edge*>           edges;
};

#endif  //  MERFIN_GLOBAL
