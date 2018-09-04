
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

#include "meryl.H"



void
merylOperation::findMinCount(void) {
  _count = _actCount[0];
  for (uint32 ii=1; ii<_actLen; ii++)
    if (_actCount[ii] < _count)
      _count = _actCount[ii];
}



void
merylOperation::findMaxCount(void) {
  _count = _actCount[0];
  for (uint32 ii=1; ii<_actLen; ii++)
    if (_count < _actCount[ii])
      _count = _actCount[ii];
}



void
merylOperation::findSumCount(void) {
  _count = 0;
  for (uint32 ii=0; ii<_actLen; ii++)
    _count += _actCount[ii];
}



void
merylOperation::initializeThreshold(void) {

  //  If no thresholds to set, nothing to do.

  if ((_fracDist == DBL_MAX) &&
      (_wordFreq == DBL_MAX))
    return;

  //  The problem with using more than one database is that the number of
  //  distinct kmers is not known.

  if (_inputs.size() != 1) {
    fprintf(stderr, "ERROR: operation most-frequent can work with only one meryl database.\n");
    exit(1);
  }

  //  Later, we could allow streaming operations, and construction of statistics on
  //  the fly.

  bool    allDatabase = true;

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if (_inputs[ii]->isFromDatabase() == false) {
      fprintf(stderr, "ERROR: input '%s' to operation most-frequent is not a meryl database.\n",
              _inputs[ii]->_name);
      allDatabase = false;
    }
  }

  if (allDatabase == false)
    exit(1);

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    _inputs[ii]->_stream->loadStatistics();

  kmerCountStatistics  *stats = _inputs[0]->_stream->stats();


#warning "need to decide whic direction we're going?"

  if (_fracDist < DBL_MAX) {
    uint64  nKmers       = 0;
    uint64  nKmersTarget = _fracDist * stats->numDistinct();

    for (uint32 ii=0; ii<stats->numFrequencies(); ii++) {
      nKmers += stats->numKmersAtFrequency(ii);

      if (nKmers >= nKmersTarget) {
        _threshold = ii;
        break;
      }
    }

    fprintf(stderr, "For fraction-distinct %f, found threshold %lu\n", _fracDist, _threshold);
  }

#warning "rounding issues in word-frequency!"

  if (_wordFreq < DBL_MAX) {
    _threshold = _wordFreq * stats->numTotal();

    fprintf(stderr, "For word-frequency %f, found threshold %lu\n", _wordFreq, _threshold);
  }

  //  Cleanup.

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    _inputs[ii]->_stream->dropStatistics();
}



bool
merylOperation::initialize(void) {
  bool  proceed = true;

  fprintf(stderr, "INITIALIZE\n");

  //  Initialize all the inputs this operation might have.

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    _inputs[ii]->initialize();

  //  Set up the output for the specific kmer database file we're processing.
  //  Note that if this _was_ a counting operation, nextMer_doCounting() just
  //  above will have already created an output for the count, written the
  //  data, and removed the _output pointer.

  if (_output) {
    _output->initialize();
    _writer = _output->getStreamWriter(_fileNumber);
  }

  //  The threshold operations need to decide on a threshold based on the histogram.

  initializeThreshold();

  //  If configuring, or if the operation is pass-through with no output,
  //  don't stream the mers.  This only matters for the root node; the return
  //  value for all other nodes is ignored (those are called above, when
  //  initializing the inputs to this node).

  if (_onlyConfig == true)
    proceed = false;

  if ((_operation == opPassThrough) &&   //  'meryl print DATABASE' uses opPassThrough.
      (_printer   == NULL))              //  but has _printer set.
    proceed = false;                     //  Counting operations do not set _printer.

  return(proceed);
}



//  Perform the counting operation, then close the output.
//
void
merylOperation::doCounting(void) {

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    _inputs[ii]->initialize();

  bool    doSimple  = false;
  uint32  wPrefix   = 0;
  uint64  nPrefix   = 0;
  uint32  wData     = 0;
  uint64  wDataMask = 0;

  configureCounting(_maxMemory,
                    doSimple,
                    wPrefix,
                    nPrefix,
                    wData,
                    wDataMask);

  //if (_operation == opCountSimple)
  //  doSimple = true;

  if (_kmer.merSize() <= 16)
    doSimple = true;

  if (_maxMemory < (uint64)10 * 1024 * 1024 * 1024)
    doSimple = false;

  omp_set_num_threads(_maxThreads);

  if (doSimple)
    countSimple();
  else
    count(wPrefix, nPrefix, wData, wDataMask);

  clearInputs();

  delete _output;
  _output = NULL;
}



//  Convert the presumed counting operation into a pass-through operation.
//  The merylOpStack (meryl.C) will assign inputs/outputs only to the
//  first file, and that is handled in doCounting() above.
//
//  All we need to do here is reset the operation and add an input to the
//  freshly constructed meryl database.
//
//  If there is no inputName, we're just attempting to configure for Canu.
//
void
merylOperation::convertToPassThrough(char *inputName) {

  //clearInputs();
  //clearOutput();

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "merylOp::nextMer()-- CONVERTING '%s' to '%s'.\n",
            toString(opCount), toString(opPassThrough));

  _operation = opPassThrough;

  if ((inputName != NULL) && (inputName[0] != 0) && (_onlyConfig == false))
    addInput(new kmerCountFileReader(inputName));
}



bool
merylOperation::nextMer(void) {

  char  kmerString[256];

  //  Get some logging out of the way.

  if (_verbosity >= sayDetails) {
    fprintf(stderr, "\n");
    fprintf(stderr, "merylOp::nextMer()-- STARTING for operation %s\n",
            toString(_operation));

    if (_verbosity >= sayEverything)
      for (uint32 ii=0; ii<_inputs.size(); ii++)
        fprintf(stderr, "merylOp::nextMer()--   CURRENT STATE: input %s kmer %s count " F_U64 " %s\n",
                _inputs[ii]->_name,
                _inputs[ii]->_kmer.toString(kmerString),
                _inputs[ii]->_count,
                _inputs[ii]->_valid ? "valid" : "INVALID");
  }

  //  Grab the next mer for every input that was active in the last iteration.
  //  (on the first call, all inputs were 'active' last time)
  //
  for (uint32 ii=0; ii<_actLen; ii++) {
    if (_verbosity >= sayDetails)
      fprintf(stderr, "merylOp::nextMer()-- CALL NEXTMER on input actIndex " F_U32 "\n", _actIndex[ii]);
    _inputs[_actIndex[ii]]->nextMer();
  }

  _actLen = 0;

  //  Find the smallest kmer in the _inputs, and save their counts in _actCount.
  //  Mark which input was used in _actIndex.

#if 0
  if (_verbosity >= sayEverything)
    for (uint32 ii=0; ii<_inputs.size(); ii++)
      fprintf(stderr, "merylOp::nextMer()--   BEFORE OPERATION: input %s kmer %s count " F_U64 " %s\n",
              _inputs[ii]->_name,
              _inputs[ii]->_kmer.toString(kmerString),
              _inputs[ii]->_count,
              _inputs[ii]->_valid ? "valid" : "INVALID");
#endif

  //  Build a list of the inputs that have the smallest kmer.

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if (_inputs[ii]->_valid == false)
      continue;

    //  If we have no active kmer, or the input kmer is smaller than the one we
    //  have, reset the list.

    if ((_actLen == 0) ||
        (_inputs[ii]->_kmer < _kmer)) {
      _actLen = 0;
      _kmer              = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_count;
      _actIndex[_actLen] = ii;
      _actLen++;

      if (_verbosity >= sayDetails)
        fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s. reset\n", _kmer.toString(kmerString), _inputs[ii]->_name);
    }

    //  Otherwise, if the input kmer is the one we have, save the count to the list.

    else if (_inputs[ii]->_kmer == _kmer) {
      //_kmer             = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_count;
      _actIndex[_actLen] = ii;
      _actLen++;

      if (_verbosity >= sayDetails)
        fprintf(stderr, "merylOp::nextMer()-- Active kmer %s from input %s\n", _kmer.toString(kmerString), _inputs[ii]->_name);
    }

    //  Otherwise, the input kmer comes after the one we're examining, ignore it.

    else {
    }
  }

  //  If no active kmers, we're done.  Several bits of housekeeping need to be done:
  //
  //  Histogram operations need to finish up and report the histogram now.
  //  Alternatively, it could be done in the destructor.
  //
  //  Any outputs need to call finishIteration() to rename and/or merge their intermediate outputs.

  if (_actLen == 0) {
    if (_verbosity >= sayDetails) {
      fprintf(stderr, "merylOp::nextMer()-- No inputs found, all done here.\n");
      fprintf(stderr, "\n");
    }

    _valid = false;

    if (_operation == opHistogram)
      reportHistogram();

    delete _writer;
    _writer = NULL;

    return(false);
  }

  //  Otherwise, active kmers!  Figure out what the count should be.

  if (_verbosity >= sayDetails)
    fprintf(stderr, "merylOp::nextMer()-- op %s activeLen " F_U32 " kmer %s\n", toString(_operation), _actLen, _kmer.toString(kmerString));

  //  If math-subtract gets implemented, use negative-zero to mean "don't output" and positive-zero
  //  to mean zero.  For now, count=0 means don't output.

  //  Set the count to zero, meaning "don't output the kmer".  Intersect depends on this,
  //  skipping most of it's work if all files don't have the kmer.
  _count = 0;

  switch (_operation) {
    case opCount:
    case opCountForward:
    case opCountReverse:
      fprintf(stderr, "ERROR: got %s, but shouldn't have.\n", toString(_operation));
      assert(0);
      break;

    case opPassThrough:                     //  Result of counting kmers.  Guaranteed to have
      _count = _actCount[0];                //  exactly one input file.
      break;

    case opLessThan:
      _count = (_actCount[0]  < _threshold) ? _actCount[0] : 0;
      break;

    case opGreaterThan:
      _count = (_actCount[0]  > _threshold) ? _actCount[0] : 0;
      break;

    case opAtLeast:
      _count = (_actCount[0] >= _threshold) ? _actCount[0] : 0;
      break;

    case opAtMost:
      _count = (_actCount[0] <= _threshold) ? _actCount[0] : 0;
      break;

    case opEqualTo:
      _count = (_actCount[0] == _threshold) ? _actCount[0] : 0;
      break;

    case opNotEqualTo:
      _count = (_actCount[0] != _threshold) ? _actCount[0] : 0;
      break;

    case opIncrease:
      if (UINT64_MAX - _actCount[0] < _mathConstant)
        _count = UINT64_MAX;    //  OVERFLOW!
      else
        _count = _actCount[0] + _mathConstant;
      break;

    case opDecrease:
      if (_actCount[0] < _mathConstant)
        _count = 0;             //  UNDERFLOW!
      else
        _count = _actCount[0] - _mathConstant;
      break;

    case opMultiply:
      if (UINT64_MAX / _actCount[0] < _mathConstant)
        _count = UINT64_MAX;    //  OVERFLOW!
      else
        _count = _actCount[0] * _mathConstant;
      break;

    case opDivide:
      if (_mathConstant == 0)
        _count = 0;             //  DIVIDE BY ZERO!
      else
        _count = _actCount[0] / _mathConstant;
      break;

    case opModulo:
      if (_mathConstant == 0)
        _count = 0;             //  DIVIDE BY ZERO!
      else
        _count = _actCount[0] % _mathConstant;
      break;

    case opUnion:                           //  Union
      _count = _actLen;
      break;

    case opUnionMin:                        //  Union, retain smallest count
      findMinCount();
      break;

    case opUnionMax:                        //  Union, retain largest count
      findMaxCount();
      break;

    case opUnionSum:                        //  Union, sum all counts
      findSumCount();
      break;

    case opIntersect:                       //  Intersect
      if (_actLen == _inputs.size())
        _count = _actCount[0];
      break;

    case opIntersectMin:                    //  Intersect, retain smallest count
      if (_actLen == _inputs.size())
        findMinCount();
      break;

    case opIntersectMax:                    //  Intersect, retain largest count
      if (_actLen == _inputs.size())
        findMaxCount();
      break;

    case opIntersectSum:                    //  Intersect, sum all counts
      if (_actLen == _inputs.size())
        findSumCount();
      break;

    case opDifference:
      if ((_actLen == 1) && (_actIndex[0] == 0))
        _count = _actCount[0];
      break;

    case opSymmetricDifference:
      if (_actLen == 1)
        _count = _actCount[0];
      break;

    case opCompare:
      if       (_actLen == 1) {
        char  str[33];

        fprintf(stdout, "kmer %s only in input %u\n",
                _kmer.toString(str), _actIndex[0]);
      }
      else if ((_actLen == 2) && (_actCount[0] != _actCount[1])) {
        char  str[33];

        fprintf(stdout, "kmer %s has value %lu in input 1 != value %lu in input 2\n",
                _kmer.toString(str), _actCount[0], _actCount[1]);
      }
      else {
      }

    case opHistogram:
      break;

    case opNothing:
      break;
  }  

  //  If flagged for output, output!

  if ((_output != NULL) &&
      (_count  > 0)) {
    _writer->addMer(_kmer, _count);
  }

  //  If flagged for printing, print!

  if ((_printer != NULL) &&
      (_count > 0)) {
    char  flags[4] = { 0 };  //  Default, no flags (and no space) printed.

    if (_kmer.isCanonical()) {
      flags[0] = '\t';
      flags[1] = 'C';
    }

    if (_kmer.isPalindrome()) {
      flags[0] = '\t';
      flags[1] = 'P';
    }

    //fprintf(_printer, "%s\t" F_U64 "%s\n", _kmer.toString(kmerString), _actCount[0], flags);
    fprintf(_printer, "%s\t" F_U64 "\n", _kmer.toString(kmerString), _actCount[0]);
  }

  //  

  if (_verbosity >= sayDetails) {
    fprintf(stderr, "merylOp::nextMer()-- FINISHED for operation %s with kmer %s count " F_U64 "%s\n",
            toString(_operation), _kmer.toString(kmerString), _count, ((_output != NULL) && (_count != 0)) ? " OUTPUT" : "");
    fprintf(stderr, "\n");
  }

  return(true);
}
