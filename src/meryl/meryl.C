
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
 *    Brian P. Walenz beginning on 2017-SEP-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "meryl.H"
#include "strings.H"
#include "system.H"


//  In meryOp-count.C
uint64
findMaxInputSizeForMemorySize(uint32 kMerSize, uint64 memorySize);


bool
isDigit(char c) {
  return(('0' <= c) && (c <= '9'));
}

bool
isNumber(char *s, char dot='.') {

  if ((s    == NULL) ||
      (s[0] == 0))
    return(false);

  for (uint32 ii=0; s[ii] != 0; ii++)
    if ((isDigit(s[ii]) == false) &&
        (s[ii] != dot))
      return(false);

  return(true);
}



class merylOpStack {
public:
  merylOpStack(uint32 nFiles) {
    _nFiles     = nFiles;
    _stacks     = new stack <merylOperation *> [nFiles];
    _operations = new vector<merylOperation *> [nFiles];
  };

  ~merylOpStack() {
    delete [] _stacks;
    delete [] _operations;
  };

  uint32    numberOfOperations(void) {
    return(_operations[0].size());
  };

  uint32    numberOfFiles(void) {
    return(_nFiles);
  };

  merylOperation *getOp(uint32 file) {
    return(_stacks[file].top());
  };

  merylOperation *getOp(uint32 opNum, uint32 fileNum) {
    assert(fileNum < _nFiles);
    assert(opNum   < _operations[fileNum].size());
    return(_operations[fileNum][opNum]);
  };

  void            pushOp(merylOp opName, uint32 threads, uint64 memory) {
    for (uint32 ff=0; ff<_nFiles; ff++) {
      merylOperation *newOp = new merylOperation(opName, ff, threads, memory);

      if (_stacks[ff].empty() == false)        //  If a command exists, the new command
        _stacks[ff].top()->addInput(newOp);    //  supplies input to the existing command.

      _stacks[ff].push(newOp);                 //  Make the new command the current command.
      _operations[ff].push_back(newOp);
    }
  };

  void            popOp(void) {
    for (uint32 ff=0; ff<_nFiles; ff++)
      _stacks[ff].pop();
  };

  //  For input from a database, we need to create new reader objects for
  //  each thread - for simplicity, we just make a new object for each input
  //  file.
  //
  void    addInput(merylFileReader *reader) {

    //#pragma omp parallel for schedule(dynamic, 1)
    for (uint32 ff=0; ff<_nFiles; ff++)
      _stacks[ff].top()->addInput(new merylFileReader(reader->filename(), ff));

    delete reader;
  };

  //  dnaSeqFile and sqStore inputs are only used by counting, and that
  //  doesn't use the threaded merylOperation scheme here.
  //
  void    addInput(dnaSeqFile *sequence) {
    _stacks[0].top()->addInput(sequence);
  };

#ifdef CANU
  void    addInput(sqStore *store, uint32 segment, uint32 segmentMax) {
    _stacks[0].top()->addInput(store, segment, segmentMax);
  };
#endif

  //  Counting operations only need the output associated with the first
  //  file, and associating with other files just makes life difficult and/or
  //  dangerous, so don't.
  //
  void    addOutput(char *writerName) {
    merylFileWriter   *writer = new merylFileWriter(writerName);

    if (isCounting())
      _stacks[0].top()->addOutput(writer);
    else
      for (uint32 ff=0; ff<_nFiles; ff++)
        _stacks[ff].top()->addOutput(writer);
  };

  void    addPrinter(char *printerName, bool ACGTorder) {
    char  T[FILENAME_MAX+1] = { 0 };
    char  N[FILENAME_MAX+1] = { 0 };

    if ((printerName == NULL) ||
        (strcmp(printerName, "-") == 0)) {
      for (uint32 ff=0; ff<_nFiles; ff++)
        _stacks[ff].top()->addPrinter(stdout, ACGTorder);
      return;
    }

    strncpy(T, printerName, FILENAME_MAX);

    char   *pre = T;
    char   *suf = strchr(T, '#');
    uint32  len = 0;

    while ((suf) && (*suf == '#')) {
      *suf = 0;
      len++;
      suf++;
    }

    for (uint32 ff=0; ff<_nFiles; ff++) {
      if (len == 0)
        snprintf(N, FILENAME_MAX, "%s.%d", printerName, ff);
      else
        snprintf(N, FILENAME_MAX, "%s%0*d%s", pre, len, ff, suf);

      _stacks[ff].top()->addPrinter(AS_UTL_openOutputFile(N), ACGTorder);
   }
  };



  void       setThreshold(uint64 p) {
    for (uint32 ff=0; ff<_nFiles; ff++) {
      _stacks[ff].top()->setThreshold(p);
      _stacks[ff].top()->setConstant(p);
    }
  };

  void       setFractionDistinct(double p) {
    for (uint32 ff=0; ff<_nFiles; ff++)
      _stacks[ff].top()->setFractionDistinct(p);
  };

  void       setWordFrequency(double p) {
    for (uint32 ff=0; ff<_nFiles; ff++)
      _stacks[ff].top()->setWordFrequency(p);
  };

  void       setExpectedNumberOfKmers(uint64 n) {
    for (uint32 ff=0; ff<_nFiles; ff++)
      _stacks[ff].top()->setExpectedNumberOfKmers(n);
  };

  void       setCountSuffix(char *s) {
    for (uint32 ff=0; ff<_nFiles; ff++)
      _stacks[ff].top()->setCountSuffix(s);
  };




  void       setMemoryLimit(uint64 m) {
    for (uint32 ff=0; ff<_nFiles; ff++)
      _stacks[ff].top()->setMemoryLimit(m);
  };

  void       setThreadLimit(uint32 t) {
    for (uint32 ff=0; ff<_nFiles; ff++)
      _stacks[ff].top()->setThreadLimit(t);
  };

  uint64     size(void)                    { return(_stacks[0].size());                  };
  bool       empty(void)                   { return(_stacks[0].size() == 0);             };

  bool       isCounting(void)              { return(_stacks[0].top()->isCounting());     };
  bool       isNormal(void)                { return(_stacks[0].top()->isNormal());       };
  bool       needsParameter(void)          { return(_stacks[0].top()->needsParameter()); };

private:
  uint32                     _nFiles;
  stack<merylOperation *>   *_stacks;
  vector<merylOperation *>  *_operations;
};




int
main(int argc, char **argv) {
  uint32                    optStringLen = 0;
  char                      optString[FILENAME_MAX+1];
  char                      inoutName[FILENAME_MAX+1];
  char                      indexName[FILENAME_MAX+1];
  char                      sqInfName[FILENAME_MAX+1];
  char                      sqRdsName[FILENAME_MAX+1];

#warning HARD CODED NUMBER OF FILES!
  merylOpStack              opStack(64);

  merylOp                   opName         = opNothing;

  uint32                    outputArg      = UINT32_MAX;
  uint32                    printerArg     = UINT32_MAX;

  char                     *writerName     = NULL;
  char                     *printerName    = NULL;

  bool                      printACGTorder = false;

  merylFileReader          *reader         = NULL;
  dnaSeqFile               *sequence       = NULL;
#ifdef CANU
  sqStore                  *store          = NULL;
#endif

  uint32                    terminating    = 0;

  uint32                    physThreads    = omp_get_max_threads();     //  Absolute maximum limits on
  uint64                    physMemory     = getPhysicalMemorySize();   //  memory= and threads= values.

  uint32                    allowedThreads = 0;                         //  Global limits, if memory= or
  uint64                    allowedMemory  = physMemory;                //  threads= is set before any operation.

  uint32                    usedThreads    = 0;                         //  The number of threads we need to request.

  uint32                    segment        = 1;
  uint32                    segmentMax     = 1;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  for (int32 arg=1; arg < argc; arg++) {

    //  Save a copy of the string.  Ugly code here, but simplifies handling
    //  of the brackets later (and argv might be read-only).

    optStringLen = strlen(argv[arg]);

    strncpy(optString, argv[arg], FILENAME_MAX);

    //  Ignore '[' at the start of the string.  Their purpose is in the matching ']' which
    //  tells us to stop adding inputs to the current command.
    //
    //  There should only be one opening bracket.

    if (optString[0] == '[') {
      for (uint32 ii=0; ii<optStringLen; ii++)
        optString[ii] = optString[ii+1];

      optStringLen--;
    }

    //  If we have a ] as the last character, strip it off and remember that we need to
    //  close the command on the stack after we process this arg.
    //
    //  We can get any number of closing brackets.

    while ((optStringLen > 0) &&
           (optString[optStringLen-1] == ']')) {
      optString[optStringLen-1] = 0;
      optStringLen--;

      terminating++;
    }

    //  Save a few copies of the command line word.

    strncpy(inoutName, optString, FILENAME_MAX);
    snprintf(indexName, FILENAME_MAX, "%s/merylIndex", optString);
    snprintf(sqInfName, FILENAME_MAX, "%s/info",       optString);
    snprintf(sqRdsName, FILENAME_MAX, "%s/reads",      optString);

    //  Scan for debug options.

    if (strcmp(optString, "dumpIndex") == 0) {               //  Report the index for the dataset.
      arg++;                                                 //  It's just the parameters used for encoding.
      delete new merylFileReader(argv[arg++], true);
      continue;
    }

    if (strcmp(optString, "dumpFile") == 0) {                //  Dump the index for a single data file.
      arg++;
      dumpMerylDataFile(argv[arg++]);
      continue;
    }

    //  Scan for options.  If any trigger, set the option and move on to the next word on the command line.

    if     ((strcmp(optString, "-h")   == 0) ||
            (strcmp(optString, "help") == 0)) {
      err.push_back(NULL);
      continue;
    }

    //  Kmer size.
    else if ((optStringLen > 2) &&
             (strncmp(optString, "k=", 2) == 0) &&
             (isNumber(optString + 2) == true)) {
      kmerTiny::setSize(strtouint32(optString + 2));
      continue;
    }

    //  Number of kmers expected for counting.
    else if ((opStack.size() > 0) &&
             (optStringLen > 2) &&
             (strncmp(optString, "n=", 2) == 0) &&
             (isNumber(optString + 2) == true)) {
      opStack.setExpectedNumberOfKmers(strtouint64(optString + 2));
      continue;
    }

    //  A suffix to filter kmers by when counting.
    else if ((optStringLen > 13) &&
             (strncmp(optString, "count-suffix=", 13) == 0)) {
      fprintf(stderr, "COUNT SUFFIX '%s'\n", optString + 12);
      opStack.setCountSuffix(optString + 13);
      continue;
    }


    //  Threshold values for less-than, etc, specifed as a fraction
    //  of the total distinct kmers, or as a word-frequency, or as
    //  an absolute count.

    else if ((opStack.size() > 0) &&
             (opStack.needsParameter() == true) &&
             (strncmp(optString, "d=", 2) == 0) &&
             (isNumber(optString + 2))) {
      opStack.setFractionDistinct(strtodouble(optString + 2));
      continue;
    }
    else if ((opStack.size() > 0) &&
             (opStack.needsParameter() == true) &&
             (strncmp(optString, "distinct=", 9) == 0) &&
             (isNumber(optString + 9))) {
      opStack.setFractionDistinct(strtodouble(optString + 9));
      continue;
    }

    else if ((opStack.size() > 0) &&
             (opStack.needsParameter() == true) &&
             (strncmp(optString, "f=", 2) == 0) &&
             (isNumber(optString + 2))) {
      opStack.setWordFrequency(strtodouble(optString + 2));
      continue;
    }
    else if ((opStack.size() > 0) &&
             (opStack.needsParameter() == true) &&
             (strncmp(optString, "word-frequency=", 15) == 0) &&
             (isNumber(optString + 15))) {
      opStack.setWordFrequency(strtodouble(optString + 15));
      continue;
    }

    else if ((opStack.size() > 0) &&
             (opStack.needsParameter() == true) &&
             (strncmp(optString, "t=", 15) == 0) &&
             (isNumber(optString + 2))) {
      opStack.setThreshold(strtouint64(optString + 2));
      continue;
    }
    else if ((opStack.size() > 0) &&
             (opStack.needsParameter() == true) &&
             (strncmp(optString, "threshold=", 10) == 0) &&
             (isNumber(optString + 10))) {
      opStack.setThreshold(strtouint64(optString + 10));
      continue;
    }
    else if ((opStack.size() > 0) &&
             (opStack.needsParameter() == true) &&
             (isNumber(optString))) {
      opStack.setThreshold(strtouint64(optString));
      continue;
    }




    //  Memory limit, either global or per-task.
    //
    //  We used to fail if the memory requested was larger than the physical memory
    //  on the machine.  I removed the check so that 'configureOnly' operations will
    //  work on any size machine.
    //
    else if ((optStringLen > 7) &&
             (strncmp(optString, "memory=", 7) == 0) &&
             (isNumber(optString + 7) == true)) {
      uint64 memory = (uint64)(strtodouble(optString + 7) * 1024 * 1024 * 1024);

      if (opStack.size() == 0)
        allowedMemory = memory;
      else
        opStack.setMemoryLimit(memory);

      continue;
    }

    //  Thread limit, either global or per-task.
    else if ((optStringLen > 8) &&
             (strncmp(optString, "threads=", 8) == 0) &&
             (isNumber(optString + 8) == true)) {
      uint32 threads = strtouint32(optString + 8);

      if (threads > usedThreads)          //  Remember the highest thread count requested.
        usedThreads = threads;

      if (opStack.size() == 0)            //  If nothing on the stack, set the global threads allowed.
        allowedThreads = threads;         //  Otherwise, set the threads allowed for just this piece.
      else
        opStack.setThreadLimit(threads);

      continue;
    }

    //  Segment of input, for counting from seqStore.
#ifdef CANU
    else if ((optStringLen > 8) &&
             (strncmp(optString, "segment=", 8) == 0) &&
             (isNumber(optString + 8, '/') == true)) {
      decodeRange(optString + 8, segment, segmentMax);
      continue;
    }
#endif

    else if (strncmp(optString, "-V", 2) == 0) {      //  Anything that starts with -V
      for (uint32 vv=1; vv<strlen(optString); vv++)   //  increases verbosity by the
        merylOperation::increaseVerbosity();          //  number of letters.
      continue;
    }

    else if (strcmp(optString, "-Q") == 0) {
      merylOperation::beQuiet();
      continue;
    }

    else if (strcmp(optString, "-P") == 0) {
      merylOperation::showProgress();
      continue;
    }

    else if (strcmp(optString, "-C") == 0) {
      merylOperation::onlyConfigure();
      continue;
    }

    else if (strcmp(optString, "-E") == 0) {
      findMaxInputSizeForMemorySize(strtouint32(argv[arg+1]), (uint64)(1000000000 * strtodouble(argv[arg+2])));
      continue;
    }



    //
    //  Parse this word.  Decide if it's a new operation, or an output name, or an input file.
    //

    if      (0 == optStringLen)
      ;  //  Got a single bracket, nothing to do here except make it not be an error.

    else if (0 == strcmp(optString, "count"))                  opName = opCount;
    else if (0 == strcmp(optString, "count-forward"))          opName = opCountForward;
    else if (0 == strcmp(optString, "count-reverse"))          opName = opCountReverse;

    else if (0 == strcmp(optString, "less-than"))              opName = opLessThan;
    else if (0 == strcmp(optString, "greater-than"))           opName = opGreaterThan;
    else if (0 == strcmp(optString, "at-least"))               opName = opAtLeast;
    else if (0 == strcmp(optString, "at-most"))                opName = opAtMost;
    else if (0 == strcmp(optString, "equal-to"))               opName = opEqualTo;
    else if (0 == strcmp(optString, "not-equal-to"))           opName = opNotEqualTo;

    else if (0 == strcmp(optString, "increase"))               opName = opIncrease;
    else if (0 == strcmp(optString, "decrease"))               opName = opDecrease;
    else if (0 == strcmp(optString, "multiply"))               opName = opMultiply;
    else if (0 == strcmp(optString, "divide"))                 opName = opDivide;
    else if (0 == strcmp(optString, "modulo"))                 opName = opModulo;

    else if (0 == strcmp(optString, "union"))                  opName = opUnion;
    else if (0 == strcmp(optString, "union-min"))              opName = opUnionMin;
    else if (0 == strcmp(optString, "union-max"))              opName = opUnionMax;
    else if (0 == strcmp(optString, "union-sum"))              opName = opUnionSum;

    else if (0 == strcmp(optString, "intersect"))              opName = opIntersect;
    else if (0 == strcmp(optString, "intersect-min"))          opName = opIntersectMin;
    else if (0 == strcmp(optString, "intersect-max"))          opName = opIntersectMax;
    else if (0 == strcmp(optString, "intersect-sum"))          opName = opIntersectSum;

    else if (0 == strcmp(optString, "difference"))             opName = opDifference;
    else if (0 == strcmp(optString, "symmetric-difference"))   opName = opSymmetricDifference;

    else if (0 == strcmp(optString, "histogram"))              opName = opHistogram;
    else if (0 == strcmp(optString, "statistics"))             opName = opStatistics;

    else if (0 == strcmp(optString, "compare"))                opName = opCompare;

    //  Handle output names.

    else if (0 == strcmp(optString, "output")) {          //  Flag the next arg as the output name for a database
      outputArg = arg + 1;                                //  if we see 'output'.
    }

    else if (arg == outputArg) {                          //  If this is the output name, make a new
      writerName = duplicateString(inoutName);            //  output writer.
    }

    //  Handle printer names.

    else if (0 == strcmp(optString, "print")) {           //  Flag the next arg as the output name for printing
      printerArg     = arg + 1;                           //  if we see 'print'.
      printACGTorder = false;
    }
    else if (0 == strcmp(optString, "printACGT")) {       //  Flag the next arg as the output name for printing
      printerArg     = arg + 1;                           //  if we see 'print'.
      printACGTorder = true;
    }

    else if ((arg == printerArg) &&                       //  If this is the printer name, and not a meryl database, make
             (fileExists(indexName) == false)) {          //  a new file to print to.  Note: this isn't triggered if the
      printerName = duplicateString(inoutName);           //  arg is an op, the if-cascade stops when the op is parsed.
    }

    //  Handle inputs.

    else if (fileExists(indexName) == true) {             //  Make a reader if the arg is a meryl database.
      reader = new merylFileReader(inoutName);
    }

    else if ((opStack.size() > 0) &&                      //  If a counting command exists, add a sequence file.
             (opStack.isCounting()   == true) &&
             (fileExists(inoutName)  == true)) {
      sequence = new dnaSeqFile(inoutName);
    }

#ifdef CANU
    else if ((opStack.size() > 0) &&                      //  If a counting command exists, add a Canu seqStore.
             (opStack.isCounting()   == true) &&
             (fileExists(sqInfName)  == true) &&
             (fileExists(sqRdsName)  == true)) {
      store = new sqStore(inoutName);
    }
#endif

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Don't know what to do with '%s'.", optString);
      err.push_back(s);
    }

    //
    //  With the argument parsed, do something with it.
    //  Order is quite important here, as they're all independent tests.
    //

    if ((arg == printerArg) &&                            //  We wanted to find a printer name here, but found
        (printerName == NULL)) {                          //  something else; make the print go to stdout.
      printerName = duplicateString("-");
    }

    if ((printerName != NULL) &&                          //  If a printer and a reader exist, but no operation,
        (reader      != NULL) &&                          //  the user asked to just print a database.
        (opName      == opNothing)) {                     //  Add a pass-through operation to give it
      opName = opPassThrough;                             //  something to hang on to.
    }

    if (opName != opNothing) {                            //  Add any just-parsed command to the stack.
      opStack.pushOp(opName,
                     (allowedThreads == 0) ? physThreads : allowedThreads,
                     allowedMemory);
      opName = opNothing;
    }

    //
    //  Attach outputs and inputs to the top operation on the stack.  If the stack is empty,
    //  do nothing with the output/input this time, wait until the next argument.
    //

    if ((writerName != NULL) &&
        (opStack.size() > 0)) {
      opStack.addOutput(writerName);
      delete [] writerName;
      writerName = NULL;
    }

    if ((printerName != NULL) &&
        (opStack.size() > 0)) {
      opStack.addPrinter(printerName, printACGTorder);
      delete [] printerName;
      printerName = NULL;
    }

    if ((reader != NULL) &&
        (opStack.size() > 0)) {
      opStack.addInput(reader);
      reader = NULL;
    }

    if ((sequence != NULL) &&
        (opStack.size() > 0)) {
      opStack.addInput(sequence);
      sequence = NULL;
    }

#ifdef CANU
    if ((store != NULL) &&
        (opStack.size() > 0)) {
      opStack.addInput(store, segment, segmentMax);
      store      = NULL;
      segment    = 1;
      segmentMax = 1;
    }
#endif

    //  Finally, if we've been told to terminate the command, do so.

    for (; terminating > 0; terminating--)
      opStack.popOp();
  }

  //  If any errors, fail.

  if ((argc == 1) ||        //  No commands
      (err.size() > 0)) {   //  Errors
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  A meryl command line is formed as a series of commands and files, possibly\n");
    fprintf(stderr, "  grouped using square brackets.  Each command operates on the file(s) that\n");
    fprintf(stderr, "  are listed after it.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  COMMANDS:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    print                display kmers on the screen as 'kmer<tab>count'.  accepts exactly one input.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    count                Count the occurrences of canonical kmers in the input.  must have 'output' specified.\n");
    fprintf(stderr, "    count-forward        Count the occurrences of forward kmers in the input.  must have 'output' specified.\n");
    fprintf(stderr, "    count-reverse        Count the occurrences of reverse kmers in the input.  must have 'output' specified.\n");
    fprintf(stderr, "      k=<K>              create mers of size K bases (mandatory).\n");
    fprintf(stderr, "      n=<N>              expect N mers in the input (optional; for precise memory sizing).\n");
    fprintf(stderr, "      memory=M           use no more than (about) M GB memory.\n");
    fprintf(stderr, "      threads=T          use no more than T threads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    less-than N          return kmers that occur fewer than N times in the input.  accepts exactly one input.\n");
    fprintf(stderr, "    greater-than N       return kmers that occur more than N times in the input.  accepts exactly one input.\n");
    fprintf(stderr, "    equal-to N           return kmers that occur exactly N times in the input.  accepts exactly one input.\n");
    fprintf(stderr, "    not-equal-to N       return kmers that do not occur exactly N times in the input.  accepts exactly one input.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    increase X           add X to the count of each kmer.\n");
    fprintf(stderr, "    decrease X           subtract X from the count of each kmer.\n");
    fprintf(stderr, "    multiply X           multiply the count of each kmer by X.\n");
    fprintf(stderr, "    divide X             divide the count of each kmer by X.\n");
    fprintf(stderr, "    modulo X             set the count of each kmer to the remainder of the count divided by X.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    union                return kmers that occur in any input, set the count to the number of inputs with this kmer.\n");
    fprintf(stderr, "    union-min            return kmers that occur in any input, set the count to the minimum count\n");
    fprintf(stderr, "    union-max            return kmers that occur in any input, set the count to the maximum count\n");
    fprintf(stderr, "    union-sum            return kmers that occur in any input, set the count to the sum of the counts\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    intersect            return kmers that occur in all inputs, set the count to the count in the first input.\n");
    fprintf(stderr, "    intersect-min        return kmers that occur in all inputs, set the count to the minimum count.\n");
    fprintf(stderr, "    intersect-max        return kmers that occur in all inputs, set the count to the maximum count.\n");
    fprintf(stderr, "    intersect-sum        return kmers that occur in all inputs, set the count to the sum of the counts.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    difference           return kmers that occur in the first input, but none of the other inputs\n");
    fprintf(stderr, "    symmetric-difference return kmers that occur in exactly one input\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  MODIFIERS:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    output O             write kmers generated by the present command to an output  meryl database O\n");
    fprintf(stderr, "                         mandatory for count operations.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  EXAMPLES:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Example:  Report 22-mers present in at least one of input1.fasta and input2.fasta.\n");
    fprintf(stderr, "            Kmers from each input are saved in meryl databases 'input1' and 'input2',\n");
    fprintf(stderr, "            but the kmers in the union are only reported to the screen.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "            meryl print \\\n");
    fprintf(stderr, "                    union \\\n");
    fprintf(stderr, "                      [count k=22 input1.fasta output input1] \\\n");
    fprintf(stderr, "                      [count k=22 input2.fasta output input2]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Example:  Find the highest count of each kmer present in both files, save the kmers to\n");
    fprintf(stderr, "            database 'maxCount'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "            meryl intersect-max input1 input2 output maxCount\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Example:  Find unique kmers common to both files.  Brackets are necessary\n");
    fprintf(stderr, "            on the first 'equal-to' command to prevent the second 'equal-to' from\n");
    fprintf(stderr, "            being used as an input to the first 'equal-to'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "            meryl intersect [equal-to 1 input1] equal-to 1 input2\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii] != NULL)
        fprintf(stderr, "%s\n", err[ii]);

    exit(1);
  }

  //  If nothing on the stack, nothing to do.

  if (opStack.size() == 0)
    exit(1);

  //  Pop the stack until we get back to the root operation.

  while (opStack.size() > 1)
    opStack.popOp();

  //  Enable threads.  If usedThreads is non-zero, the user supplied a threads= argument, and
  //  usedThreads is set to the maximum of those.  Otherwise, allow one thread per CPU.

  {
    uint32  threads = (usedThreads > 0) ? usedThreads : physThreads;

    fprintf(stderr, "Enabling %u threads.\n", threads);
    omp_set_num_threads(threads);
  }

  //  opHistogram is limited to showing only histograms already stored in a database.
  //  opHistogram cannot take input from anything but a database.
  //  opHistogram does not generate kmer outputs.
  //  So, if the top op is histogram, all we can do is load the histogram and dump it.
  //
  //  Eventully, maybe, opHistogram will pass through mers (but then we need to figure out
  //  where to report the histogram).
  //
  //  Eventually, maybe, opHistogram will allow input from a kmer stream.

  if (opStack.getOp(0)->getOperation() == opHistogram) {
    opStack.getOp(0)->reportHistogram();
    exit(0);
  }

  if (opStack.getOp(0)->getOperation() == opStatistics) {
    opStack.getOp(0)->reportStatistics();
    exit(0);
  }

  //  Counting operations are a big headache.  They don't fit into the
  //  tree nicely:
  //   - they do their own threading, so only one thread can start the operation
  //   - when done, they transform themselves into a pass-through operation that
  //     simply reads the (just counted) input and passes kmers through.
  //
  //  So, we special case them here.  Process in order, counting, writing the
  //  output, and converting to a pass-through operation.

  for (uint32 opNum=0; opNum<opStack.numberOfOperations(); opNum++) {
    merylOperation *op = opStack.getOp(opNum, 0);
    char            name[FILENAME_MAX + 1] = { 0 };

    if (op->isCounting() == false)                        //  If not a counting operation,
      continue;                                           //  skip it.

    if (op->getOutputName())                              //  Save the output name, so we
      strncpy(name, op->getOutputName(), FILENAME_MAX);   //  know which input to open later.

    op->doCounting();                                     //  Do the counting.

    for (uint32 fn=0;                                     //  Convert all the files to pass
         fn < opStack.numberOfFiles();                    //  through operations.
         fn++)
      opStack.getOp(opNum, fn)->convertToPassThrough(name, fn);
  }

  //  If there is an operation (debug operations and -h have no operations)
  //  keep calling nextMer() on that top operation until there are no more mers.
  //
  //  op->initialize() returns false if we're only configuring, or if the top
  //  operation on the stack is a counting operation.

  uint32  nf = opStack.numberOfFiles();

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<nf; ff++) {
    merylOperation *op = opStack.getOp(ff);

    if (op->initialize() == true)
      while (op->nextMer() == true)
        ;

    op->finalize();
  }

  //  Now that everything is done, delete!
  //  Output presents a problem, in that everyone has a copy
  //  of it, but only one can delete it.  This is hardcoded
  //  in merylOperation::~merylOperation() so that only
  //  the first one (ff==0 here) deletes the output object.

  for (uint32 ff=0; ff<nf; ff++) {
    merylOperation *op = opStack.getOp(ff);

    delete op;
  }

  fprintf(stderr, "Bye.\n");

  return(0);
}
