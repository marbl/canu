
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#ifndef MERYLOP_H
#define MERYLOP_H

#include "meryl.H"

#include "sqStore.H"

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
  uint64  guesstimateNumberOfkmersInInput_sqStore(sqStore *store, uint32 bgnID, uint32 endID);
  uint64  guesstimateNumberOfkmersInInput(void);

  //  Returns true 
  void    configureCounting(uint64   memoryAllowed,      //  Input:  Maximum allowed memory in bytes
                            bool    &useSimple_,         //  Output: algorithm to use
                            uint32  &wPrefix_,           //  Output: Number of bits in the prefix (== bucket address)
                            uint64  &nPrefix_,           //  Output: Number of prefixes there are (== number of buckets)
                            uint32  &wData_,             //  Output: Number of bits in kmer data
                            uint64  &wDataMask_);        //  Output: A mask to return just the data of the mer);

public:
  void    addInput(merylOperation *operation);
  void    addInput(kmerCountFileReader *reader);
  void    addInput(dnaSeqFile *sequence);
  void    addInput(sqStore *store, uint32 segment, uint32 segmentMax);

  void    addOutput(kmerCountFileWriter *writer);
  char   *getOutputName(void);

  void    addPrinter(FILE *printer);
  char   *getPrinterName(void);

  void    setOperation(merylOp op) { _operation = op;    };
  merylOp getOperation(void)       { return(_operation); };

  void    setParameter(uint64 p)               { _parameter   = p;  };
  void    setExpectedNumberOfKmers(uint64 n)   { _expNumKmers = n;  };

  void    setMemoryLimit(uint64 m)             { _maxMemory = m;    };
  void    setThreadLimit(uint32 t)             { _maxThreads = t;   };

  bool    isCounting(void) {
    return((_operation == opCount)        ||
           (_operation == opCountForward) ||
           (_operation == opCountReverse));
  };

  bool    isNormal(void) {
    return(isCounting() == false);
  };

  bool    needsParameter(void) {
    return((_operation == opLessThan)    ||
           (_operation == opGreaterThan) ||
           (_operation == opAtLeast)     ||
           (_operation == opAtMost)      ||
           (_operation == opEqualTo)     ||
           (_operation == opNotEqualTo)  ||
           (_operation == opIncrease)    ||
           (_operation == opDecrease)    ||
           (_operation == opMultiply)    ||
           (_operation == opDivide)      ||
           (_operation == opModulo));
  };

  kmer   &theFMer(void)            { return(_kmer);   };
  uint64  theCount(void)           { return(_count);  };

public:
  void    doCounting(void);
  void    convertToPassThrough(char *inputName);

public:
  bool    initialize(void);
  bool    nextMer(void);

  bool    validMer(void)           { return(_valid);  };

  void    countSimple(void);
  void    count(uint32  wPrefix,
                uint64  nPrefix,
                uint32  wData,
                uint64  wDataMask);

  void    reportHistogram(void);

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

  merylOp                        _operation;

  uint64                         _parameter;
  uint64                         _expNumKmers;

  uint32                         _maxThreads;
  uint64                         _maxMemory;

  kmerCountStatistics           *_stats;

  kmerCountFileWriter           *_output;
  kmerCountStreamWriter         *_writer;
  FILE                          *_printer;

  uint32                         _fileNumber;

  uint32                         _actLen;
  uint64                        *_actCount;
  uint32                        *_actIndex;

  kmer                           _kmer;
  uint64                         _count;
  bool                           _valid;

  static bool                    _onlyConfig;
  static bool                    _showProgress;
  static merylVerbosity          _verbosity;
};


char const *toString(merylOp op);


#endif  //  MERYLOP_H
