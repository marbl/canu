
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

#include "meryl.H"
#ifdef CANU
#include "sqStore.H"
#endif

bool            merylOperation::_onlyConfig   = false;
bool            merylOperation::_showProgress = false;
merylVerbosity  merylOperation::_verbosity    = sayStandard;


merylOperation::merylOperation(merylOp op, uint32 ff, uint32 threads, uint64 memory) {

  _isMultiSet    = false;  //  set in initialize().

  _operation     = op;

  _mathConstant  = 0;
  _threshold     = UINT64_MAX;
  _fracDist      = DBL_MAX;
  _wordFreq      = DBL_MAX;

  memset(_countSuffixString, 0, sizeof(char) * 65);
  _countSuffixLength = 0;
  //_countSuffix.clear();   //  No clear?

  _expNumKmers   = 0;

  _maxThreads    = threads;
  _maxMemory     = memory;

  _stats         = NULL;

  _output        = NULL;
  _writer        = NULL;

  _printer       = NULL;

  _fileNumber    = ff;

  _actLen        = 0;
  _actCount      = new kmvalu [1024];
  _actIndex      = new uint32 [1024];

  _value         = 0;
  _valid         = true;
}


merylOperation::~merylOperation() {

  clearInputs();

  delete    _stats;

  assert(_writer == NULL);

  if (_fileNumber == 0)   //  The output is shared among all operations for this set of files.
    delete  _output;      //  Only one operation can delete it.

  if (_printer != stdout)
    AS_UTL_closeFile(_printer);

  delete [] _actCount;
  delete [] _actIndex;
}




void
merylOperation::clearInputs(void) {

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];

  _inputs.clear();

  _actLen = 0;
}



void
merylOperation::checkInputs(const char *name) {

  if ((_actLen > 1) && ((_operation == opPassThrough) ||
                        (_operation == opLessThan)    ||
                        (_operation == opGreaterThan) ||
                        (_operation == opEqualTo)     ||
                        (_operation == opHistogram))) {
    fprintf(stderr, "merylOp::addInput()-- ERROR: can't add input '%s' to operation '%s': only one input supported.\n",
            name, toString(_operation));
    exit(1);
  }
}



void
merylOperation::addInput(merylOperation *operation) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from operation '%s' to operation '%s'\n",
            toString(operation->_operation), toString(_operation));

  _inputs.push_back(new merylInput(operation));
  _actIndex[_actLen++] = _inputs.size() - 1;

  if (operation->_operation == opHistogram)
    fprintf(stderr, "ERROR: operation '%s' can't be used as an input: it doesn't supply kmers.\n", toString(operation->_operation)), exit(1);

  if (_operation == opHistogram)
    fprintf(stderr, "ERROR: operation '%s' can't take input from '%s': it can only accept databases.\n", toString(_operation), toString(operation->_operation)), exit(1);

  checkInputs(toString(operation->getOperation()));
}


void
merylOperation::addInput(merylFileReader *reader) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
            reader->filename(), toString(_operation));

  _inputs.push_back(new merylInput(reader->filename(), reader));
  _actIndex[_actLen++] = _inputs.size() - 1;

  checkInputs(reader->filename());
}


void
merylOperation::addInput(dnaSeqFile *sequence) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
            sequence->filename(), toString(_operation));

  _inputs.push_back(new merylInput(sequence->filename(), sequence));
  _actIndex[_actLen++] = _inputs.size() - 1;

  if (isCounting() == false)
    fprintf(stderr, "ERROR: operation '%s' cannot use sequence files as inputs.\n", toString(_operation)), exit(1);
}



#ifdef CANU
void
merylOperation::addInput(sqStore *store, uint32 segment, uint32 segmentMax) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from sqStore '%s' to operation '%s'\n",
            store->sqStore_path(), toString(_operation));

  _inputs.push_back(new merylInput(store->sqStore_path(), store, segment, segmentMax));
  _actIndex[_actLen++] = _inputs.size() - 1;

  if (isCounting() == false)
    fprintf(stderr, "ERROR: operation '%s' cannot use sqStore as inputs.\n", toString(_operation)), exit(1);
}
#endif



void
merylOperation::addOutput(merylFileWriter *writer) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding output to file '%s' from operation '%s'\n",
            writer->filename(), toString(_operation));

  if (_output)
    fprintf(stderr, "ERROR: already have an output set!\n"), exit(1);

  if (_operation == opHistogram)
    fprintf(stderr, "ERROR: operation '%s' can't use 'output' modifier.\n", toString(_operation));

  _output = writer;
}



char *
merylOperation::getOutputName(void) {
  if (_output)
    return(_output->filename());
  else
    return(NULL);
}



//  We're all done processing this operation.  Clean up what we can.
//  The _output CANNOT be deleted until all operations are done with it.
//  Yes, I should be pointer counting or something smart like that.
void
merylOperation::finalize(void) {

  clearInputs();

  delete [] _actCount;   _actCount = NULL;
  delete [] _actIndex;   _actIndex = NULL;
}



void
merylOperation::addPrinter(FILE *printer, bool ACGTorder) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding printer to %s from operation '%s'\n",
            (printer == stdout) ? "(stdout)" : "(a file)",
            toString(_operation));

  if (_printer)
    fprintf(stderr, "ERROR: already have a printer set!\n"), exit(1);

  if (_operation == opHistogram)
    fprintf(stderr, "ERROR: operation '%s' can't use 'output' modifier.\n", toString(_operation));

  _printer        = printer;
  _printACGTorder = ACGTorder;
}




char const *
toString(merylOp op) {
  switch (op) {
    case opCount:                return("opCount");                break;
    case opCountForward:         return("opCountForward");         break;
    case opCountReverse:         return("opCountReverse");         break;
    case opPassThrough:          return("opPassThrough");          break;

    case opLessThan:             return("opLessThan");             break;
    case opGreaterThan:          return("opGreaterThan");          break;
    case opAtLeast:              return("opAtLeast");              break;
    case opAtMost:               return("opAtMost");               break;
    case opEqualTo:              return("opEqualTo");              break;
    case opNotEqualTo:           return("opNotEqualTo");           break;

    case opIncrease:             return("opIncrease");             break;
    case opDecrease:             return("opDecrease");             break;
    case opMultiply:             return("opMultiply");             break;
    case opDivide:               return("opDivide");               break;
    case opModulo:               return("opModulo");               break;

    case opUnion:                return("opUnion");                break;
    case opUnionMin:             return("opUnionMin");             break;
    case opUnionMax:             return("opUnionMax");             break;
    case opUnionSum:             return("opUnionSum");             break;

    case opIntersect:            return("opIntersect");            break;
    case opIntersectMin:         return("opIntersectMin");         break;
    case opIntersectMax:         return("opIntersectMax");         break;
    case opIntersectSum:         return("opIntersectSum");         break;

    case opDifference:           return("opDifference");           break;
    case opSymmetricDifference:  return("opSymmetricDifference");  break;

    case opHistogram:            return("opHistogram");            break;
    case opStatistics:           return("opStatistics");           break;

    case opCompare:              return("opCompare");              break;

    case opNothing:              return("opNothing");              break;
  }

  assert(0);
  return(NULL);
}

