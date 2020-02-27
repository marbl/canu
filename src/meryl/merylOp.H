
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef MERYLOP_H
#define MERYLOP_H

#include "meryl.H"

//  These would (maybe) be handy, but they change the order of kmers, so we need
//  to sort again.
//
//  opReverse,             //  Reverse the kmer.     Accepts a single input only.
//  opComplement,          //  Complement the kmer.  Accepts a single input only.
//  opReverseComplement,   //  Rev-Compl the kmer.   Accepts a single input only.


enum merylOp {
  opCount,
  opCountForward,
  opCountReverse,
  opPassThrough,         //  Internal operation; behaves just like opUnion, but only a single input.

  opLessThan,
  opGreaterThan,
  opAtLeast,
  opAtMost,
  opEqualTo,
  opNotEqualTo,

  opIncrease,
  opDecrease,
  opMultiply,
  opDivide,
  opModulo,

  opUnion,
  opUnionMin,
  opUnionMax,
  opUnionSum,

  opIntersect,
  opIntersectMin,
  opIntersectMax,
  opIntersectSum,

  opDifference,
  opSymmetricDifference,

  opHistogram,
  opStatistics,

  opCompare,

  opNothing
};


enum merylVerbosity {
  sayNothing,
  sayStandard,
  sayConstruction,
  sayDetails,
  sayEverything
};


//  For the simple counting algorithm, how bit of a word to use for initial counts.
typedef uint16 lowBits_t;


class merylOperation {
public:
  merylOperation(merylOp op=opNothing, uint32 ff=UINT32_MAX, uint32 threads=1, uint64 memory=0);
  ~merylOperation();

private:
  void    clearInputs(void);
  void    checkInputs(const char *name);

  uint64  guesstimateNumberOfkmersInInput_dnaSeqFile(dnaSeqFile *sequence);
#ifdef CANU
  uint64  guesstimateNumberOfkmersInInput_sqStore(sqStore *store, uint32 bgnID, uint32 endID);
#endif
  uint64  guesstimateNumberOfkmersInInput(void);

  //  Returns true
  void    configureCounting(uint64   memoryAllowed,      //  Input:  Maximum allowed memory in bytes
                            bool    &useSimple_,         //  Output: algorithm to use
                            uint32  &wPrefix_,           //  Output: Number of bits in the prefix (== bucket address)
                            uint64  &nPrefix_,           //  Output: Number of prefixes there are (== number of buckets)
                            uint32  &wData_,             //  Output: Number of bits in kmer data
                            kmdata  &wDataMask_);        //  Output: A mask to return just the data of the mer);

public:
  void    addInput(merylOperation *operation);
  void    addInput(merylFileReader *reader);
  void    addInput(dnaSeqFile *sequence);
#ifdef CANU
  void    addInput(sqStore *store, uint32 segment, uint32 segmentMax);
#endif

  void    addOutput(merylFileWriter *writer);
  char   *getOutputName(void);

  void    finalize(void);

  void    addPrinter(FILE *printer, bool ACGTorder);
  char   *getPrinterName(void);

  void    setOperation(merylOp op) { _operation = op;    };
  merylOp getOperation(void)       { return(_operation); };

  void    setConstant(uint64 p)                { _mathConstant      = p; };
  void    setThreshold(uint64 p)               { _threshold         = p; };
  void    setFractionDistinct(double p)        { _fracDist          = p; };
  void    setWordFrequency(double p)           { _wordFreq          = p; };

  void    setExpectedNumberOfKmers(uint64 n)   { _expNumKmers       = n; };

  void    setCountSuffix(char *s)              {
    _countSuffixLength = strlen(s);

    for (uint32 ii=0; ii<_countSuffixLength; ii++) {
      _countSuffixString[ii] = s[ii];

      _countSuffix.addR(s[ii]);
    }
  };

  void    setMemoryLimit(uint64 m)             { _maxMemory         = m; };
  void    setThreadLimit(uint32 t)             { _maxThreads        = t; };

  bool    isCounting(void) {
    return((_operation == opCount)        ||
           (_operation == opCountForward) ||
           (_operation == opCountReverse));
  };

  bool    isNormal(void) {
    return(isCounting() == false);
  };

  bool    needsParameter(void) {
    return((_operation == opLessThan)     ||
           (_operation == opGreaterThan)  ||
           (_operation == opAtLeast)      ||
           (_operation == opAtMost)       ||
           (_operation == opEqualTo)      ||
           (_operation == opNotEqualTo)   ||
           (_operation == opIncrease)     ||
           (_operation == opDecrease)     ||
           (_operation == opMultiply)     ||
           (_operation == opDivide)       ||
           (_operation == opModulo));
  };

  //bool    needsRealParameter(void) {
  //  return((_operation == opMostFrequent));
  //};

  bool    isMultiSet(void)         { return(_isMultiSet); };

  kmer   &theFMer(void)            { return(_kmer);       };
  kmvalu  theValue(void)           { return(_value);      };

public:
  void    doCounting(void);
  void    convertToPassThrough(char *inputName, uint32 threadFile);

public:
  void    initializeThreshold(void);
  bool    initialize(void);

private:
  void    nextMer_findSmallestNormal(void);
  void    nextMer_findSmallestMultiSet(void);
  bool    nextMer_finish(void);

public:
  bool    nextMer(void);

  bool    validMer(void)           { return(_valid);  };

  void    countSimple(void);
  void    countThreads(uint32  wPrefix,
                       uint64  nPrefix,
                       uint32  wData,
                       kmdata  wDataMask);
  void    count(uint32  wPrefix,
                uint64  nPrefix,
                uint32  wData,
                kmdata  wDataMask);

  void    reportHistogram(void);
  void    reportStatistics(void);

  static
  void    onlyConfigure(void)               { _onlyConfig   = true;       };
  static
  void    showProgress(void)                { _showProgress = true;       };
  static
  void    increaseVerbosity(void) {
    if      (_verbosity == sayNothing)        _verbosity    = sayStandard;
    else if (_verbosity == sayStandard)       _verbosity    = sayConstruction;
    else if (_verbosity == sayConstruction)   _verbosity    = sayDetails;
    else if (_verbosity == sayDetails)        _verbosity    = sayEverything;
    else                                      _verbosity    = sayEverything;
  };
  static
  void    beQuiet(void)                     { _verbosity    = sayNothing; };


private:
  void    findMinCount(void);
  void    findMaxCount(void);
  void    findSumCount(void);

  vector<merylInput *>           _inputs;
  bool                           _isMultiSet;

  merylOp                        _operation;

  uint64                         _mathConstant;
  uint64                         _threshold;
  double                         _fracDist;
  double                         _wordFreq;

  char                           _countSuffixString[65];
  uint32                         _countSuffixLength;
  kmer                           _countSuffix;

  uint64                         _expNumKmers;

  uint32                         _maxThreads;
  uint64                         _maxMemory;

  merylHistogram                *_stats;

  merylFileWriter               *_output;    //  This is the main output object, but for streaming
  merylStreamWriter             *_writer;    //  operations, _writer is used.
  FILE                          *_printer;

  bool                           _printACGTorder;

  char                           _kmerString[65];

  uint32                         _fileNumber;

  uint32                         _actLen;
  kmvalu                        *_actCount;
  uint32                        *_actIndex;

  kmer                           _kmer;
  kmvalu                         _value;
  bool                           _valid;

  static bool                    _onlyConfig;
  static bool                    _showProgress;
  static merylVerbosity          _verbosity;
};


char const *toString(merylOp op);


#endif  //  MERYLOP_H
