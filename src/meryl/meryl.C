
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
 *    kmer/meryl/meryl.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2004-APR-07
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-MAR-25 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAY-23 to 2009-AUG-07
 *      are Copyright 2005,2007-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#undef  NEW_MERYL   //  Loses last mer sometimes


#ifndef NEW_MERYL

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "meryl.H"

int
main(int argc, char **argv) {
  merylArgs   *args = new merylArgs(argc, argv);

  switch (args->personality) {
    case 'P':
      estimate(args);
      break;

    case 'B':
      build(args);
      break;

    case 'd':
      dumpDistanceBetweenMers(args);
      break;
    case 't':
      dumpThreshold(args);
      break;
    case 'p':
      dumpPositions(args);
      break;
    case 'c':
      countUnique(args);
      break;
    case 'h':
      plotHistogram(args);
      break;

    case PERSONALITY_MIN:
    case PERSONALITY_MINEXIST:
    case PERSONALITY_MAX:
    case PERSONALITY_MAXEXIST:
    case PERSONALITY_ADD:
    case PERSONALITY_AND:
    case PERSONALITY_NAND:
    case PERSONALITY_OR:
    case PERSONALITY_XOR:
      multipleOperations(args);
      break;

    case PERSONALITY_SUB:
    case PERSONALITY_DIFFERENCE:
    case PERSONALITY_ABS:
    case PERSONALITY_DIVIDE:
      binaryOperations(args);
      break;

    case PERSONALITY_LEQ:
    case PERSONALITY_GEQ:
    case PERSONALITY_EQ:
      unaryOperations(args);
      break;

    default:
      args->usage();
      fprintf(stderr, "%s: unknown personality.  Specify -P, -B, -S or -M!\n", args->execName);
      exit(1);
      break;
  }

  delete args;

  return(0);
}

#else

#include "AS_global.H"
#include "AS_UTL_fileIO.H"

#include "libmeryl.H"
#include "kMer.H"

#include <vector>
#include <stack>
using namespace std;


enum merylOp {
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
  opComplement,
  opNothing
};


char const *
toString(merylOp op) {
  switch (op) {
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
    case opComplement:           return("opComplement");           break;
    case opNothing:              return("opNothing");              break; 
  }

  assert(0);
  return(NULL);
}




class merylOperation;


class merylInput {
public:
  merylInput(char *n, merylStreamReader *s);
  merylInput(merylOperation *o);
  ~merylInput();

  void   nextMer(void);

  merylStreamReader  *_stream;
  merylOperation     *_operation;

  char                _name[FILENAME_MAX];

  kMer                _kmer;
  uint64              _count;
  bool                _valid;
};




class merylOperation {
public:
  merylOperation(merylOp op=opNothing);
  ~merylOperation();

  void    addInput(char *name, merylStreamReader *reader);
  void    addInput(merylOperation *operation);

  void    addOutput(char *name, merylStreamWriter *writer);

  void    setOperation(merylOp op) { _operation = op;    };
  merylOp getOperation(void)       { return(_operation); };

  kMer   &theFMer(void)            { return(_kmer);   };
  uint64  theCount(void)           { return(_count);  };

  bool    nextMer(void);
  bool    validMer(void)           { return(_valid);  };

private:
  void    findMinCount(void);
  void    findMaxCount(void);
  void    findSumCount(void);

  vector<merylInput *>           _inputs;

  merylOp                        _operation;

  char                           _outputName[FILENAME_MAX];
  merylStreamWriter             *_output;

  kMer                           _smallest;

  uint32                         _actLen;
  uint64                        *_actCount;
  uint32                        *_actIndex;

  kMer                           _kmer;
  uint64                         _count;
  bool                           _valid;
};







merylInput::merylInput(char *n, merylStreamReader *s) {
  _stream      = s;
  _operation   = NULL;
  _count       = 0;
  _valid       = false;

  strncpy(_name, n, FILENAME_MAX);
}


merylInput::merylInput(merylOperation *o) {
  _stream      = NULL;
  _operation   = o;
  _count       = 0;
  _valid       = false;

  strncpy(_name, toString(_operation->getOperation()), FILENAME_MAX);
}


merylInput::~merylInput() {
  fprintf(stderr, "Destroy input %s\n", _name);
  delete _stream;
  delete _operation;
}


void
merylInput::nextMer(void) {
  char kmerString[256];

  if (_stream) {
    _valid = _stream->nextMer();
    _kmer  = _stream->theFMer();
    _count = _stream->theCount();
  }

  if (_operation) {
    _valid = _operation->nextMer();
    _kmer  = _operation->theFMer();
    _count = _operation->theCount();
  }
}








merylOperation::merylOperation(merylOp op) {
  fprintf(stderr, "Create operation '%s'\n", toString(op));

  _operation     = op;

  _outputName[0] = 0;
  _output        = NULL;

  _actLen        = 0;
  _actCount      = new uint64 [1024];
  _actIndex      = new uint32 [1024];

  _count         = 0;
  _valid         = true;
}


merylOperation::~merylOperation() {
  fprintf(stderr, "Destroy op %s\n", toString(_operation));

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];

  delete    _output;
  delete [] _actCount;
  delete [] _actIndex;
}




void
merylOperation::addInput(char *name, merylStreamReader *reader) {
  fprintf(stderr, "Adding input file '%s' to operation '%s'\n",
          name, toString(_operation));

  _inputs.push_back(new merylInput(name, reader));

  _actIndex[_actLen++] = _inputs.size() - 1;
}


void
merylOperation::addInput(merylOperation *operation) {
  fprintf(stderr, "Adding input from operation '%s' to operation '%s'\n",
          toString(operation->_operation), toString(_operation));

  _inputs.push_back(new merylInput(operation));

  _actIndex[_actLen++] = _inputs.size() - 1;
}


void
merylOperation::addOutput(char *name, merylStreamWriter *writer) {
  fprintf(stderr, "Adding output to file '%s' from operation '%s'\n",
          name, toString(_operation));

  strncpy(_outputName, name, FILENAME_MAX);
  _output = writer;
}




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
  for (uint32 ii=0; ii<_actLen; ii++)
    _count += _actCount[ii];
}


bool
merylOperation::nextMer(void) {

  char  kmerString[256];

  //  Find the smallest kmer in the _inputs, and save their counts in _actCount.
  //  Mark which input was used in _actIndex.

  fprintf(stderr, "\n");

  //  Grab the next mer for every input that was active in the last iteration.
  //  (on the first call, all inputs were 'active' last time)
  //
  for (uint32 ii=0; ii<_actLen; ii++)
    _inputs[_actIndex[ii]]->nextMer();

  _actLen = 0;

  //  Log.

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "input %s kmer %s count %lu%s\n",
            _inputs[ii]->_name,
            _inputs[ii]->_kmer.merToString(kmerString),
            _inputs[ii]->_count,
            _inputs[ii]->_valid ? " valid" : "");
  }

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if (_inputs[ii]->_valid == false)
      continue;

    //  Of we have no active kmer, or the input kmer is smaller than the one we
    //  have, reset the list.

    if ((_actLen == 0) ||
        (_inputs[ii]->_kmer < _kmer)) {
      _actLen = 0;
      _kmer              = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_count;
      _actIndex[_actLen] = ii;
      _actLen++;

      fprintf(stderr, "Active kmer %s from input %s. reset\n", _kmer.merToString(kmerString), _inputs[ii]->_name);
    }

    //  Otherwise, if the input kmer is the one we have, save the count to the list.

    else if (_inputs[ii]->_kmer == _kmer) {
      //_kmer             = _inputs[ii]->_kmer;
      _actCount[_actLen] = _inputs[ii]->_count;
      _actIndex[_actLen] = ii;
      _actLen++;

      fprintf(stderr, "Active kmer %s from input %s\n", _kmer.merToString(kmerString), _inputs[ii]->_name);
    }

    //  Otherwise, the input kmer comes after the one we're examining, ignore it.

    else {
    }
  }

  //  If no active kmers, we're done.

  fprintf(stderr, "op %s activeLen %u kmer %s\n", toString(_operation), _actLen, _kmer.merToString(kmerString));

  if (_actLen == 0) {
    _valid = false;
    return(false);
  }

  //  Otherwise, active kmers!  Figure out what the count should be.

  //  If math-subtract gets implemented, use negative-zero to mean "don't output" and positive-zero
  //  to mean zero.  For now, count=0 means don't output.

  //  Set the count to zero, meaning "don't output the kmer".  Intersect depends on this,
  //  skipping most of it's work if all files don't have the kmer.
  _count = 0;

  switch (_operation) {
    case opUnion:
      _count = 1;
      break;

    case opUnionMin:
      findMinCount();
      break;

    case opUnionMax:
      findMaxCount();
      break;

    case opUnionSum:
      findSumCount();
      break;

    case opIntersect:                       //  Intersect, retain count of first file
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
      if ((_actLen == 1) &&                 //  If only found in one set, and that
          (_actIndex[0] == 0))              //  set is the first one, save it.
        _count = _actCount[0];
      break;

    case opSymmetricDifference:
      break;

    case opComplement:
      break;

    case opNothing:
      break;
  }  

  //  If flagged for output, output!

  if ((_output != NULL) &&
      (_count  != 0)) {
    fprintf(stderr, "OUTPUT %s count %lu\n", _kmer.merToString(kmerString), _count);
    _output->addMer(_kmer, _count);
  }

  return(true);
}




int
main(int argc, char **argv) {
  stack<merylOperation *>   opStack;

#warning "need to scan the input files to make sure they're all the same mersize"
  uint32                    merSize = 22;

  uint32                    outputArg = UINT32_MAX;

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    char    mcidx[FILENAME_MAX];
    char    mcdat[FILENAME_MAX];
    char    opt[FILENAME_MAX];
    uint32  optLen = strlen(argv[arg]);

    uint32   creating    = 0;
    uint32   terminating = 0;

    strncpy(opt, argv[arg], FILENAME_MAX);

    //  If we have [ as the first character, make a new operation (with no
    //  operation set yet) and add it to the inputs, then push on the stack.
    //
    //  We can get 0 or 1 open bracket at a time.  Seeing two in a row is an error,
    //  but it isn't caught.

    if (opStack.empty() == true) {
      creating = true;
    }

    if (opt[0] == '[') {
      strncpy(opt, argv[arg]+1, FILENAME_MAX);
      optLen--;

      creating = true;
    }

    //  If we have a ] as the last character, strip it off and remember.
    //
    //  We can get any number of closing brackets.

    while (opt[optLen-1] == ']') {
      opt[optLen-1] = 0;
      optLen--;

      terminating++;
    }

    //  Now that brackets are stripped, make meryl database names for the arg.

    snprintf(mcidx, FILENAME_MAX, "%s.mcidx", opt);
    snprintf(mcdat, FILENAME_MAX, "%s.mcdat", opt);

    merylOp             op     = opNothing;
    merylStreamWriter  *writer = NULL;
    merylStreamReader  *reader = NULL;

    if      (opt[0] == 0)
      ;  //  Got a single bracket, nothing to do here except make it not be an error.

    else if (strcmp(opt, "union") == 0)                  op = opUnion;
    else if (strcmp(opt, "union-min") == 0)              op = opUnionMin;
    else if (strcmp(opt, "union-max") == 0)              op = opUnionMax;
    else if (strcmp(opt, "union-sum") == 0)              op = opUnionSum;
    else if (strcmp(opt, "intersect") == 0)              op = opIntersect;
    else if (strcmp(opt, "intersect-min") == 0)          op = opIntersectMin;
    else if (strcmp(opt, "intersect-max") == 0)          op = opIntersectMax;
    else if (strcmp(opt, "intersect-sum") == 0)          op = opIntersectSum;
    else if (strcmp(opt, "difference") == 0)             op = opDifference;
    else if (strcmp(opt, "symmetric-difference") == 0)   op = opSymmetricDifference;
    else if (strcmp(opt, "complement") == 0)             op = opComplement;

    //  If we see 'output', flag the next arg as being the output name.
    //  If this arg is flagged as output, add an output using the bracket-stripped name.
    //  If this arg is a valid meryl file, make it an input.

    else if (strcmp(opt, "output") == 0)       outputArg = arg+1;

    else if (arg == outputArg)                 writer = new merylStreamWriter(opt, merSize, 0, merSize/2, false);

    else if (AS_UTL_fileExists(mcidx) &&
             AS_UTL_fileExists(mcdat))         reader = new merylStreamReader(opt);

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option or meryl database files not found for '%s'.\n", opt);
      err.push_back(s);
    }


    //  Create a new operation or set inputs or output.  Or do nothing.

    if (op != opNothing) {
      merylOperation *newOp = new merylOperation(op);

      if (opStack.empty() == false)
        opStack.top()->addInput(newOp);

      opStack.push(newOp);
      op = opNothing;
    }

    if ((writer != NULL) &&              //  If nothing on the stack, wait for an operation
        (opStack.size() > 0)) {          //  to show up.
      opStack.top()->addOutput(opt, writer);
      writer = NULL;
    }

    if (reader != NULL) {
      opStack.top()->addInput(opt, reader);
      reader = NULL;
    }

    //  If we're done adding inputs to the operation, load the first mer from
    //  the inputs, and pop it off the stack.

    for (; terminating > 0; terminating--) {
      opStack.top()->nextMer();
      opStack.pop();
    }

    arg++;
  }

  //  If any errors, fail.

  if (err.size() > 0) {
    exit(1);
  }

  //  Now just walk through the kmers until nothing is left.

  merylOperation *op = opStack.top();

  fprintf(stderr, "START operation %s\n", toString(op->getOperation()));

  while (op->nextMer() == true)
    ;

  //  Done!

  fprintf(stderr, "DONE operation %s\n", toString(op->getOperation()));

  delete op;


  return(0);
}


#endif
