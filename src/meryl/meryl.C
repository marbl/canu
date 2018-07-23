
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
#include "system.H"


bool
isDigit(char c) {
  return(('0' <= c) && (c <= '9'));
}

bool
isNumber(char *s) {

  if (s == NULL)
    return(false);

  for (uint32 ii=0; s[ii] != 0; ii++)
    if ((isDigit(s[ii]) == false) &&
        (s[ii] != '.'))
      return(false);

  return(true);
}




int
main(int argc, char **argv) {
  uint32                    optStringLen = 0;
  char                      optString[FILENAME_MAX+1];
  char                      inoutName[FILENAME_MAX+1];
  char                      indexName[FILENAME_MAX+1];

  stack<merylOperation *>   opStack;
  merylOp                   opName         = opNothing;
  uint32                    outputArg      = UINT32_MAX;
  uint32                    printArg       = UINT32_MAX;
  kmerCountFileWriter      *writer         = NULL;
  FILE                     *printer        = NULL;
  kmerCountFileReader      *reader         = NULL;
  dnaSeqFile               *sequence       = NULL;

  uint32                    terminating    = 0;

  uint32                    physThreads    = omp_get_max_threads();     //  Absolute maximum limits on
  uint64                    physMemory     = getPhysicalMemorySize();   //  memory= and threads= values.

  uint32                    allowedThreads = physThreads;               //  Global limits, if memory= or
  uint64                    allowedMemory  = physMemory;                //  threads= is set before any operation.


  vector<char *>  err;
  for (int32 arg=1; arg < argc; arg++) {

    //  Save a few copies of the command line word.

    optStringLen = strlen(argv[arg]);

    strncpy(optString, argv[arg], FILENAME_MAX);
    strncpy(inoutName, argv[arg], FILENAME_MAX);
    strncpy(indexName, argv[arg], FILENAME_MAX);
    strncat(indexName, "/merylIndex", FILENAME_MAX - optStringLen - 1);

    //  Scan for debug options.

    if (strcmp(optString, "dumpIndex") == 0) {
      arg++;
      delete new kmerCountFileReader(argv[arg++], true, true);
      continue;
    }

    if (strcmp(optString, "dumpFile") == 0) {
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
      opStack.top()->setExpectedNumberOfKmers(strtouint64(optString + 2));
      continue;
    }

    //  Threshold values for less-than, greater-than and equal-to are just a number.
    else if ((opStack.size() > 0) &&
             (opStack.top()->needsParameter() == true) &&
             (isNumber(optString))) {
      opStack.top()->setParameter(strtouint64(optString));
      continue;
    }

    //  Memory limit, either global or per-task.
    else if ((optStringLen > 7) &&
             (strncmp(optString, "memory=", 7) == 0) &&
             (isNumber(optString + 7) == true)) {
      uint64 memory = (uint64)(strtodouble(optString + 7) * 1024 * 1024 * 1024);

      if (memory > physMemory) {
        char *s = new char [1024];
        snprintf(s, 1024, "Requested memory '%s' (GB) is more than physical memory %.2f GB.",
                 optString, physMemory / 1024.0 / 1024.0 / 1024.0);
        err.push_back(s);
      }

      if (opStack.size() == 0)
        allowedMemory = memory;
      else
        opStack.top()->setMemoryLimit(memory);
             
      continue;
    }

    //  Thread limit, either global or per-task.
    else if ((optStringLen > 8) &&
             (strncmp(optString, "threads=", 8) == 0) &&
             (isNumber(optString + 8) == true)) {
      uint32 threads = strtouint32(optString + 8);

      if (opStack.size() == 0) {
        allowedThreads = threads;
        omp_set_num_threads(allowedThreads);
      }
      else {
        opStack.top()->setThreadLimit(threads);
      }

      continue;
    }


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

    //  Now, parse this word.  Decide if it's a new operation, or an output name, or an input file.

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

    else if (0 == strcmp(optString, "output"))            //  Flag the next arg as the output name for a database
      outputArg = arg + 1;                                //  if we see 'output'.

    else if (0 == strcmp(optString, "print"))             //  Flag the next arg as the output name for printing
      printArg = arg + 1;                                 //  if we see 'print'.

    else if (arg == outputArg)                            //  If this is the output name, make a new
      writer = new kmerCountFileWriter(inoutName);        //  output writer.

    else if ((arg == printArg) &&                         //  This this _should_ have been an output name,
             (directoryExists(inoutName))) {              //  but is instead a directory; user is asking for
      printer = stdout;                                   //  the database in the directory to be printed.
      reader  = new kmerCountFileReader(inoutName);
    }

    else if (arg == printArg)                             //  If this is the printer name, make a new file to print
      printer = AS_UTL_openOutputFile(inoutName);         //  to.  Note: this isn't triggered if the arg is an op.

    else if ((opStack.size() > 0) &&                      //  If a command exists,
             (opStack.top()->isCounting()  == false) &&   //  and it isn't for counting,
             (fileExists(indexName) == true))             //  and the meryl index file exists,
      reader = new kmerCountFileReader(inoutName);        //  add a kmerCountFile as input to the current command.

    else if ((opStack.size() > 0) &&                      //  If a command exists,
             (opStack.top()->isCounting()  == true) &&    //  and it IS for counting,
             (fileExists(inoutName)  == true))            //  and the file exists,
      sequence = new dnaSeqFile(inoutName);               //  add a sequence file as input to the current command.

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Don't know what to do with '%s'.", optString);
      err.push_back(s);
    }

    //  A couple of special cases for printing that can't be handled above.

    if ((printer != NULL) &&       //  A printer and a reader exist but there isn't anything
        (reader  != NULL) &&       //  on the stack, we've been told to print a meryl database.
        (opStack.size() == 0))     //  Use the NoOp op opPassThrough just to give the print
      opName = opPassThrough;      //  something to hang on.

    if ((arg == printArg) &&       //  This this _should_ have been an output name,
        (opName != opNothing))     //  but is instead an operation, send the printer
      printer = stdout;            //  to stdout.

    //  Now, do something with the parsed word.
    //    If 'op' is set, make a new command.
    //    If 'writer' exists, set the output of the top most command to that.
    //    If 'printer' exists, set the printer of the top most command to that.
    //    If 'reader' or 'sequence' exist, add it to the inputs of the top most command.

    if (opName != opNothing) {
      merylOperation *newOp = new merylOperation(opName, allowedThreads, allowedMemory);

      if (opStack.empty() == false)        //  If a command exists, the new command
        opStack.top()->addInput(newOp);    //  supplies input to the existing command.

      opStack.push(newOp);                 //  Make the new command the current command.
      opName = opNothing;
    }

    if ((writer != NULL) &&                //  Add the writer to the top most command, if
        (opStack.size() > 0)) {            //  one actually exists.  If not, wait until the
      opStack.top()->addOutput(writer);    //  command is created.
      writer = NULL;
    }

    if ((printer != NULL) &&               //  Same as for writer, except for the printer.
        (opStack.size() > 0)) {
      opStack.top()->addPrinter(printer);
      printer = NULL;
    }

    if (reader != NULL) {                  //  Add the reader to the top most command.  The
      opStack.top()->addInput(reader);     //  top most command always exists (else we'd error out
      reader = NULL;                       //  when creating the reader object above).
    }

    if (sequence != NULL) {                //  Same story, different object.
      opStack.top()->addInput(sequence);
      sequence = NULL;
    }

    //  Finally, if we've been told to terminate the command, do so.

    for (; terminating > 0; terminating--)
      opStack.pop();
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
    opStack.pop();

  //  opHistogram is limited to showing only histograms already stored in a database.
  //  opHistogram cannot take input from anything but a database.
  //  opHistogram does not generate kmer outputs.
  //  So, if the top op is histogram, all we can do is load the histogram and dump it.
  //
  //  Eventully, maybe, opHistogram will pass through mers (but then we need to figure out
  //  where to report the histogram).
  //
  //  Eventually, maybe, opHistogram will allow input from a kmer stream.

  if (opStack.top()->getOperation() == opHistogram) {
    opStack.top()->nextMer();          //  To load the file.
    opStack.top()->reportHistogram();
    exit(0);
  }
   
  //  If there is an operation (debug operations and -h have no operations)
  //  keep calling nextMer() on that top operation until there are no more mers.

  merylOperation *op = opStack.top();

  //fprintf(stderr, "Detected %u available threads and %.3f GB memory.\n",
  //        physThreads, physMemory / 1024.0 / 1024.0 / 1024.0);

  //  Initialization and counting are special cases that we don't want
  //  to test for over and over and over on every kmer.  If initialize
  //  returns false, then the first operation was a counting operation,
  //  and we don't need to run through the kmers.

  if (op->initialize(true) == true)
    while (op->nextMer() == true)
      ;

  delete op;  //  Deletes all the child operations too.

  fprintf(stderr, "Bye.\n");

  return(0);
}
