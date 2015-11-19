
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
 *    kmer/meryl/args.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2004-MAR-31 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-05 to 2004-JUL-01
 *      are Copyright 2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAR-20 to 2014-APR-11
 *      are Copyright 2005-2011,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2015-MAY-29
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "meryl.H"


//  Some string handling utilities.
//
bool
writeString(const char *str, FILE *F) {
  errno = 0;

  uint32 len = 0;
  if (str) {
    len = (uint32)strlen(str) + 1;
    fwrite(&len, sizeof(uint32), 1, F);
    fwrite( str, sizeof(char), len, F);
  } else {
    fwrite(&len, sizeof(uint32), 1, F);
  }

  if (errno) {
    fprintf(stderr, "writeString()-- Failed to write string of length "F_U32": %s\n", len, strerror(errno));
    fprintf(stderr, "writeString()-- First 80 bytes of string is:\n");
    fprintf(stderr, "%80.80s\n", str);
    return(false);
  }

  return(true);
}

char*
readString(FILE *F) {
  errno = 0;

  uint32 len = 0;
  fread(&len, sizeof(uint32), 1, F);
  if (errno) {
    fprintf(stderr, "readString()-- Failed to read string: %s\n", strerror(errno));
    exit(1);
  }

  char *str = 0L;

  if (len > 0) {
    str = new char [len];
    fread(str, sizeof(char), len, F);
    if (errno) {
      fprintf(stderr, "readString()-- Failed to read string: %s\n", strerror(errno));
      exit(1);
    }
  }

  return(str);
}

char*
duplString(char *str) {
  char   *dupstr = 0L;
  if (str) {
    uint32  len = (uint32)strlen(str);
    dupstr = new char [len+1];
    strcpy(dupstr, str);
  }
  return(dupstr);
}



void
merylArgs::usage(void) {
  fprintf(stderr, "usage: %s [personality] [global options] [options]\n", execName);
  fprintf(stderr, "\n");
  fprintf(stderr, "where personality is:\n");
  fprintf(stderr, "        -P -- compute parameters\n");
  fprintf(stderr, "        -B -- build table\n");
  fprintf(stderr, "        -S -- scan table\n");
  fprintf(stderr, "        -M -- \"math\" operations\n");
  fprintf(stderr, "        -D -- dump table\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "-P:  Given a sequence file (-s) or an upper limit on the\n");
  fprintf(stderr, "     number of mers in the file (-n), compute the table size\n");
  fprintf(stderr, "     (-t in build) to minimize the memory usage.\n");
  fprintf(stderr, "        -m #          (size of a mer; required)\n");
  fprintf(stderr, "        -c #          (homopolymer compression; optional)\n");
  fprintf(stderr, "        -p            (enable positions)\n");
  fprintf(stderr, "        -s seq.fasta  (seq.fasta is scanned to determine the number of mers)\n");
  fprintf(stderr, "        -n #          (compute params assuming file with this many mers in it)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Only one of -s, -n need to be specified.  If both are given\n");
  fprintf(stderr, "     -s takes priority.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "-B:  Given a sequence file (-s) and lots of parameters, compute\n");
  fprintf(stderr, "     the mer-count tables.  By default, both strands are processed.\n");
  fprintf(stderr, "        -f            (only build for the forward strand)\n");
  fprintf(stderr, "        -r            (only build for the reverse strand)\n");
  fprintf(stderr, "        -C            (use canonical mers, assumes both strands)\n");
  fprintf(stderr, "        -L #          (DON'T save mers that occur less than # times)\n");
  fprintf(stderr, "        -U #          (DON'T save mers that occur more than # times)\n");
  fprintf(stderr, "        -m #          (size of a mer; required)\n");
  fprintf(stderr, "        -c #          (homopolymer compression; optional)\n");
  fprintf(stderr, "        -p            (enable positions)\n");
  fprintf(stderr, "        -s seq.fasta  (sequence to build the table for)\n");
  fprintf(stderr, "        -o tblprefix  (output table prefix)\n");
  fprintf(stderr, "        -v            (entertain the user)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     By default, the computation is done as one large sequential process.\n");
  fprintf(stderr, "     Multi-threaded operation is possible, at additional memory expense, as\n");
  fprintf(stderr, "     is segmented operation, at additional I/O expense.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Threaded operation: Split the counting in to n almost-equally sized\n");
  fprintf(stderr, "     pieces.  This uses an extra h MB (from -P) per thread.\n");
  fprintf(stderr, "        -threads n    (use n threads to build)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Segmented, sequential operation: Split the counting into pieces that\n");
  fprintf(stderr, "     will fit into no more than m MB of memory, or into n equal sized pieces.\n");
  fprintf(stderr, "     Each piece is computed sequentially, and the results are merged at the end.\n");
  fprintf(stderr, "     Only one of -memory and -segments is needed.\n");
  fprintf(stderr, "        -memory mMB   (use at most m MB of memory per segment)\n");
  fprintf(stderr, "        -segments n   (use n segments)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Segmented, batched operation: Same as sequential, except this allows\n");
  fprintf(stderr, "     each segment to be manually executed in parallel.\n");
  fprintf(stderr, "     Only one of -memory and -segments is needed.\n");
  fprintf(stderr, "        -memory mMB     (use at most m MB of memory per segment)\n");
  fprintf(stderr, "        -segments n     (use n segments)\n");
  fprintf(stderr, "        -configbatch    (create the batches)\n");
  fprintf(stderr, "        -countbatch n   (run batch number n)\n");
  fprintf(stderr, "        -mergebatch     (merge the batches)\n");
  fprintf(stderr, "     Initialize the compute with -configbatch, which needs all the build options.\n");
  fprintf(stderr, "     Execute all -countbatch jobs, then -mergebatch to complete.\n");
  fprintf(stderr, "       meryl -configbatch -B [options] -o file\n");
  fprintf(stderr, "       meryl -countbatch 0 -o file\n");
  fprintf(stderr, "       meryl -countbatch 1 -o file\n");
  fprintf(stderr, "       ...\n");
  fprintf(stderr, "       meryl -countbatch N -o file\n");
  fprintf(stderr, "       meryl -mergebatch N -o file\n");
  fprintf(stderr, "     Batched mode can run on the grid.\n");
  fprintf(stderr, "        -sge        jobname      unique job name for this execution.  Meryl will submit\n");
  fprintf(stderr, "                                 jobs with name mpjobname, ncjobname, nmjobname, for\n");
  fprintf(stderr, "                                 phases prepare, count and merge.\n");
  fprintf(stderr, "        -sgebuild \"options\"    any additional options to sge, e.g.,\n");
  fprintf(stderr, "        -sgemerge \"options\"    \"-p -153 -pe thread 2 -A merylaccount\"\n");
  fprintf(stderr, "                                 N.B. - -N will be ignored\n");
  fprintf(stderr, "                                 N.B. - be sure to quote the options\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "-M:  Given a list of tables, perform a math, logical or threshold operation.\n");
  fprintf(stderr, "     Unless specified, all operations take any number of databases.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Math operations are:\n");
  fprintf(stderr, "        min       count is the minimum count for all databases.  If the mer\n");
  fprintf(stderr, "                  does NOT exist in all databases, the mer has a zero count, and\n");
  fprintf(stderr, "                  is NOT in the output.\n");
  fprintf(stderr, "        minexist  count is the minimum count for all databases that contain the mer\n");
  fprintf(stderr, "        max       count is the maximum count for all databases\n");
  fprintf(stderr, "        add       count is sum of the counts for all databases\n");
  fprintf(stderr, "        sub       count is the first minus the second (binary only)\n");
  fprintf(stderr, "        abs       count is the absolute value of the first minus the second (binary only)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Logical operations are:\n");
  fprintf(stderr, "        and       outputs mer iff it exists in all databases\n");
  fprintf(stderr, "        nand      outputs mer iff it exists in at least one, but not all, databases\n");
  fprintf(stderr, "        or        outputs mer iff it exists in at least one database\n");
  fprintf(stderr, "        xor       outputs mer iff it exists in an odd number of databases\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Threshold operations are:\n");
  fprintf(stderr, "        lessthan x            outputs mer iff it has count <  x\n");
  fprintf(stderr, "        lessthanorequal x     outputs mer iff it has count <= x\n");
  fprintf(stderr, "        greaterthan x         outputs mer iff it has count >  x\n");
  fprintf(stderr, "        greaterthanorequal x  outputs mer iff it has count >= x\n");
  fprintf(stderr, "        equal x               outputs mer iff it has count == x\n");
  fprintf(stderr, "     Threshold operations work on exactly one database.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "        -s tblprefix  (use tblprefix as a database)\n");
  fprintf(stderr, "        -o tblprefix  (create this output)\n");
  fprintf(stderr, "        -v            (entertain the user)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     NOTE:  Multiple tables are specified with multiple -s switches; e.g.:\n");
  fprintf(stderr, "              %s -M add -s 1 -s 2 -s 3 -s 4 -o all\n", execName);
  fprintf(stderr, "     NOTE:  It is NOT possible to specify more than one operation:\n");
  fprintf(stderr, "              %s -M add -s 1 -s 2 -sub -s 3\n", execName);
  fprintf(stderr, "            will NOT work.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "-D:  Dump the table (not all of these work).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     -Dd        Dump a histogram of the distance between the same mers.\n");
  fprintf(stderr, "     -Dt        Dump mers >= a threshold.  Use -n to specify the threshold.\n");
  fprintf(stderr, "     -Dc        Count the number of mers, distinct mers and unique mers.\n");
  fprintf(stderr, "     -Dh        Dump (to stdout) a histogram of mer counts.\n");
  fprintf(stderr, "     -s         Read the count table from here (leave off the .mcdat or .mcidx).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
}



void
merylArgs::clear(void) {

  execName           = 0L;
  options            = 0L;

  beVerbose          = false;
  doForward          = true;
  doReverse          = false;
  doCanonical        = false;

  inputFile          = 0L;
  outputFile         = 0L;
  queryFile          = 0L;

  merSize            = 20;
  merComp            = 0;
  positionsEnabled   = false;

  numMersEstimated   = 0;
  numMersActual      = 0;

  numBasesActual     = 0;

  mersPerBatch       = 0;
  basesPerBatch      = 0;

  numBuckets         = 0;
  numBuckets_log2    = 0;
  merDataWidth       = 0;
  merDataMask        = uint64ZERO;
  bucketPointerWidth = 0;

  numThreads         = 0;
  memoryLimit        = 0;
  segmentLimit       = 0;
  configBatch        = false;
  countBatch         = false;
  mergeBatch         = false;
  batchNumber        = 0;

  sgeJobName         = 0L;
  sgeBuildOpt        = 0L;
  sgeMergeOpt        = 0L;
  isOnGrid           = false;

  lowCount           = 0;
  highCount          = ~lowCount;
  desiredCount       = 0;

  outputCount        = 0;
  outputAll          = 0;
  outputPosition     = 0;

  mergeFilesMax      = 0;
  mergeFilesLen      = 0;
  mergeFiles         = 0L;

  personality        = 0;
}




merylArgs::merylArgs(int argc, char **argv) {

  clear();

  execName           = duplString(argv[0]);

  if (argc == 1) {
    usage();
    exit(1);
  }

  //  Count how many '-s' switches there are, then allocate space
  //  for them in mergeFiles.  We also sum the length of all options,
  //  so we can copy them into an 'options' string used when we
  //  resubmit to the grid.
  //
  uint32  optionsLen = 0;
  for (int arg=1; arg < argc; arg++) {
    optionsLen += strlen(argv[arg]) + 1;
    if (strcmp(argv[arg], "-s") == 0)
      mergeFilesMax++;
  }

  mergeFiles   = new char * [mergeFilesMax];
  options      = new char   [2 * optionsLen + 1];
  options[0]   = 0;

  bool fail = false;

  char *optptr = options;

  for (int arg=1; arg < argc; arg++) {
    if (arg > 1)
      *optptr++ = ' ';

    //  Arg!  If the arg has spaces or other stuff that the shell
    //  needs escaped we need to escape them again.  So, we copy byte
    //  by byte and insert escapes at the right points.

    for (char *op=argv[arg]; *op; op++, optptr++) {
      if (isspace(*op) || !isalnum(*op))
        if ((*op != '-') && (*op != '_') && (*op != '.') && (*op != '/'))
          *optptr++ = '\\';
      *optptr = *op;
    }

    //strcat(options, argv[arg]);
  }


  //  Parse the options
  //
  for (int arg=1; arg < argc; arg++) {
    if        (strncmp(argv[arg], "-V", 2) == 0) {
      fprintf(stdout, "meryl the Mighty Mer Counter version (no version)\n");
      exit(0);
    } else if (strcmp(argv[arg], "-m") == 0) {
      arg++;
      merSize = strtouint32(argv[arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      arg++;
      merComp = strtouint32(argv[arg]);
    } else if (strcmp(argv[arg], "-p") == 0) {
      positionsEnabled = true;
    } else if (strcmp(argv[arg], "-s") == 0) {
      arg++;
      delete [] inputFile;
      inputFile                   = duplString(argv[arg]);
      mergeFiles[mergeFilesLen++] = duplString(argv[arg]);
    } else if (strcmp(argv[arg], "-n") == 0) {
      arg++;
      numMersEstimated = strtouint64(argv[arg]);
    } else if (strcmp(argv[arg], "-f") == 0) {
      doForward   = true;
      doReverse   = false;
      doCanonical = false;
    } else if (strcmp(argv[arg], "-r") == 0) {
      doForward   = false;
      doReverse   = true;
      doCanonical = false;
    } else if (strcmp(argv[arg], "-C") == 0) {
      doForward   = false;
      doReverse   = false;
      doCanonical = true;
    } else if (strcmp(argv[arg], "-L") == 0) {
      arg++;
      lowCount = strtouint32(argv[arg]);
    } else if (strcmp(argv[arg], "-U") == 0) {
      arg++;
      highCount = strtouint32(argv[arg]);
    } else if (strcmp(argv[arg], "-o") == 0) {
      arg++;
      delete [] outputFile;
      outputFile = duplString(argv[arg]);
    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;

    } else if (strcmp(argv[arg], "-P") == 0) {
      personality = 'P';
    } else if (strcmp(argv[arg], "-B") == 0) {
      personality = 'B';
    } else if (strcmp(argv[arg], "-S") == 0) {
      personality = 'S';
    } else if (strcmp(argv[arg], "-M") == 0) {
      arg++;
      if        (strcmp(argv[arg], "merge") == 0) {
        personality = PERSONALITY_MERGE;
      } else if (strcmp(argv[arg], "min") == 0) {
        personality = PERSONALITY_MIN;
      } else if (strcmp(argv[arg], "minexist") == 0) {
        personality = PERSONALITY_MINEXIST;
      } else if (strcmp(argv[arg], "max") == 0) {
        personality = PERSONALITY_MAX;
      } else if (strcmp(argv[arg], "maxexist") == 0) {
        personality = PERSONALITY_MAXEXIST;
      } else if (strcmp(argv[arg], "add") == 0) {
        personality = PERSONALITY_ADD;
      } else if (strcmp(argv[arg], "sub") == 0) {
        personality = PERSONALITY_SUB;
      } else if (strcmp(argv[arg], "abs") == 0) {
        personality = PERSONALITY_ABS;
      } else if (strcmp(argv[arg], "divide") == 0) {
        personality = PERSONALITY_DIVIDE;
      } else if (strcmp(argv[arg], "and") == 0) {
        personality = PERSONALITY_AND;
      } else if (strcmp(argv[arg], "nand") == 0) {
        personality = PERSONALITY_NAND;
      } else if (strcmp(argv[arg], "or") == 0) {
        personality = PERSONALITY_OR;
      } else if (strcmp(argv[arg], "xor") == 0) {
        personality = PERSONALITY_XOR;
      } else if (strcmp(argv[arg], "lessthan") == 0) {
        personality = PERSONALITY_LEQ;
        arg++;
        desiredCount = strtouint32(argv[arg]) - 1;
      } else if (strcmp(argv[arg], "lessthanorequal") == 0) {
        personality = PERSONALITY_LEQ;
        arg++;
        desiredCount = strtouint32(argv[arg]);
      } else if (strcmp(argv[arg], "greaterthan") == 0) {
        personality = PERSONALITY_GEQ;
        arg++;
        desiredCount = strtouint32(argv[arg]) + 1;
      } else if (strcmp(argv[arg], "greaterthanorequal") == 0) {
        personality = PERSONALITY_GEQ;
        arg++;
        desiredCount = strtouint32(argv[arg]);
      } else if (strcmp(argv[arg], "equal") == 0) {
        personality = PERSONALITY_EQ;
        arg++;
        desiredCount = strtouint32(argv[arg]);
      } else {
        fprintf(stderr, "ERROR: unknown math personality %s\n", argv[arg]);
        exit(1);
      }
    } else if (strcmp(argv[arg], "-Dd") == 0) {
      personality = 'd';
    } else if (strcmp(argv[arg], "-Dt") == 0) {
      personality = 't';
    } else if (strcmp(argv[arg], "-Dp") == 0) {
      personality = 'p';
    } else if (strcmp(argv[arg], "-Dc") == 0) {
      personality = 'c';
    } else if (strcmp(argv[arg], "-Dh") == 0) {
      personality = 'h';
    } else if (strcmp(argv[arg], "-memory") == 0) {
      arg++;
      memoryLimit = strtouint64(argv[arg]) * 1024 * 1024;
    } else if (strcmp(argv[arg], "-segments") == 0) {
      arg++;
      segmentLimit = strtouint64(argv[arg]);
    } else if (strcmp(argv[arg], "-threads") == 0) {
      arg++;
      numThreads   = strtouint32(argv[arg]);
    } else if (strcmp(argv[arg], "-configbatch") == 0) {
      personality = 'B';
      configBatch = true;
      countBatch  = false;
      mergeBatch  = false;
      batchNumber = uint32ZERO;
    } else if (strcmp(argv[arg], "-countbatch") == 0) {
      arg++;
      personality = 'B';
      configBatch = false;
      countBatch  = true;
      mergeBatch  = false;
      batchNumber = strtouint32(argv[arg]);
    } else if (strcmp(argv[arg], "-mergebatch") == 0) {
      personality = 'B';
      configBatch = false;
      countBatch  = false;
      mergeBatch  = true;
      batchNumber = uint32ZERO;
    } else if (strcmp(argv[arg], "-sge") == 0) {
      sgeJobName = argv[++arg];
    } else if (strcmp(argv[arg], "-sgebuild") == 0) {
      sgeBuildOpt = argv[++arg];
    } else if (strcmp(argv[arg], "-sgemerge") == 0) {
      sgeMergeOpt = argv[++arg];
    } else if (strcmp(argv[arg], "-forcebuild") == 0) {
      isOnGrid = true;
    } else {
      fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
      fail = true;
    }
  }

  //  Using threads is only useful if we are not a batch.
  //
  if ((numThreads > 0) && (configBatch || countBatch || mergeBatch)) {
    if (configBatch)
      fprintf(stderr, "WARNING: -threads has no effect with -configbatch, disabled.\n");
    if (countBatch)
      fprintf(stderr, "WARNING: -threads has no effect with -countbatch, disabled.\n");
    if (mergeBatch)
      fprintf(stderr, "WARNING: -threads has no effect with -mergebatch, disabled.\n");
    numThreads = 0;
  }

  //  SGE is not useful unless we are in batch mode.
  //
  if (sgeJobName && !configBatch && !countBatch && !mergeBatch) {
    fprintf(stderr, "ERROR: -sge not useful unless in batch mode (replace -B with -configbatch)\n");
    exit(1);
  }

  if (fail)
    exit(1);
}



merylArgs::merylArgs(const char *prefix) {

  clear();

  char *filename = new char [strlen(prefix) + 17];
  sprintf(filename, "%s.merylArgs", prefix);

  errno = 0;
  FILE *F = fopen(filename, "rb");
  if (errno) {
    fprintf(stderr, "merylArgs::readConfig()-- Failed to open '%s': %s\n", filename, strerror(errno));
    exit(1);
  }

  char  magic[17] = {0};
  fread(magic, sizeof(char), 16, F);
  if (strncmp(magic, "merylBatcherv02", 16) != 0) {
    fprintf(stderr, "merylArgs::readConfig()-- '%s' doesn't appear to be a merylArgs file.\n", filename);
    exit(1);
  }

  //  Load the config, then reset the pointers.

  fread(this, sizeof(merylArgs), 1, F);

  execName    = readString(F);
  options     = 0L;
  inputFile   = readString(F);
  outputFile  = readString(F);
  queryFile   = 0L;
  sgeJobName  = readString(F);
  sgeBuildOpt = readString(F);
  sgeMergeOpt = readString(F);

  mergeFiles  = new char* [mergeFilesLen];

  for (uint32 i=0; i<mergeFilesLen; i++)
    mergeFiles[i] = readString(F);

  fclose(F);

  delete [] filename;
}



merylArgs::~merylArgs() {
  delete [] execName;
  delete [] options;
  delete [] inputFile;
  delete [] outputFile;

  for (uint32 i=0; i<mergeFilesLen; i++)
    delete [] mergeFiles[i];

  delete [] mergeFiles;
}



bool
merylArgs::writeConfig(void) {
  char *filename;

  filename = new char [strlen(outputFile) + 17];
  sprintf(filename, "%s.merylArgs", outputFile);

  errno = 0;
  FILE *F = fopen(filename, "wb");
  if (errno) {
    fprintf(stderr, "merylArgs::writeConfig()-- Failed to open '%s': %s\n", filename, strerror(errno));
    exit(1);
  }

  fwrite("merylBatcherv02", sizeof(char), 16, F);

  fwrite(this, sizeof(merylArgs), 1, F);

  writeString(execName,    F);
  writeString(inputFile,   F);
  writeString(outputFile,  F);
  writeString(sgeJobName,  F);
  writeString(sgeBuildOpt, F);
  writeString(sgeMergeOpt, F);

  for (uint32 i=0; i<mergeFilesLen; i++)
    writeString(mergeFiles[i], F);

  fclose(F);

  return(true);
}
