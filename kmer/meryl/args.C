#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bio++.H"
#include "meryl.H"


//  Some string handling utilities.
//
bool
writeString(const char *str, FILE *F) {
  errno = 0;

  u32bit len = 0;
  if (str) {
    len = (u32bit)strlen(str) + 1;
    fwrite(&len, sizeof(u32bit), 1, F);
    fwrite( str, sizeof(char), len, F);
  } else {
    fwrite(&len, sizeof(u32bit), 1, F);
  }

  if (errno) {
    fprintf(stderr, "writeString()-- Failed to write string of length "u32bitFMT": %s\n", len, strerror(errno));
    fprintf(stderr, "writeString()-- First 80 bytes of string is:\n");
    fprintf(stderr, "%80.80s\n", str);
    return(false);
  }

  return(true);
}

char*
readString(FILE *F) {
  errno = 0;

  u32bit len = 0;
  fread(&len, sizeof(u32bit), 1, F);
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
    u32bit  len = (u32bit)strlen(str);
    dupstr = new char [len+1];
    strcpy(dupstr, str);
  }
  return(dupstr);
}




const char *usagestring =
"usage: %s [personality] [global options] [options]\n"
"\n"
"where personality is:\n"
"        -P -- compute parameters\n"
"        -B -- build table\n"
"        -S -- scan table\n"
"        -M -- \"math\" operations\n"
"        -D -- dump table\n"
"\n"
"global options\n"
"        -stats file   (write run time statistics to file)\n"
"\n"
"-P:  Given a sequence file (-s) or an upper limit on the\n"
"     number of mers in the file (-n), compute the table size\n"
"     (-t in build) to minimize the memory usage.\n"
"        -m #          (size of a mer; required)\n"
"        -c #          (homopolymer compression; optional)\n"
"        -s seq.fasta  (seq.fasta is scanned to determine the number of mers)\n"
"        -n #          (compute params assuming file with this many mers in it)\n"
"\n"
"     Only one of -s, -n need to be specified.  If both are given\n"
"     -s takes priority.\n"
"\n"
"\n"
"-B:  Given a sequence file (-s) and lots of parameters, compute\n"
"     the mer-count tables.  By default, both strands are processed.\n"
"        -f            (only build for the forward strand)\n"
"        -r            (only build for the reverse strand)\n"
"        -C            (use canonical mers, assumes both strands)\n"
"        -L #          (DON'T save mers that occur less than # times)\n"
"        -U #          (DON'T save mers that occur more than # times)\n"
"        -m #          (size of a mer; required)\n"
"        -c #          (homopolymer compression; optional)\n"
"        -s seq.fasta  (sequence to build the table for)\n"
"        -o tblprefix  (output table prefix)\n"
"        -v            (entertain the user)\n"
"\n"
"     By default, the computation is done as one large sequential process.\n"
"     Multi-threaded operation is possible, at additional memory expense, as\n"
"     is segmented operation, at additional I/O expense.\n"
"\n"
"     Threaded operation: Split the counting in to n almost-equally sized\n"
"     pieces.  This uses an extra h MB (from -P) per thread.\n"
"        -threads n    (use n threads to build)\n"
"\n"
"     Segmented, sequential operation: Split the counting into pieces that\n"
"     will fit into no more than m MB of memory, or into n equal sized pieces.\n"
"     Each piece is computed sequentially, and the results are merged at the end.\n"
"     Only one of -memory and -segments is needed.\n"
"        -memory mMB   (use at most m MB of memory per segment)\n"
"        -segments n   (use n segments)\n"
"\n"
"     Segmented, batched operation: Same as sequential, except this allows\n"
"     each segment to be manually executed in parallel.\n"
"     Only one of -memory and -segments is needed.\n"
"        -memory mMB     (use at most m MB of memory per segment)\n"
"        -segments n     (use n segments)\n"
"        -configbatch    (create the batches)\n"
"        -countbatch n   (run batch number n)\n"
"        -mergebatch     (merge the batches)\n"
"     Initialize the compute with -configbatch, which needs all the build options.\n"
"     Execute all -countbatch jobs, then -mergebatch to complete.\n"
"       meryl -configbatch -B [options] -o file\n"
"       meryl -countbatch 0 -o file\n"
"       meryl -countbatch 1 -o file\n"
"       ...\n"
"       meryl -countbatch N -o file\n"
"       meryl -mergebatch N -o file\n"
"     Batched mode can run on the grid.\n"
"        -sge        jobname      unique job name for this execution.  Meryl will submit\n"
"                                 jobs with name mpjobname, ncjobname, nmjobname, for\n"
"                                 phases prepare, count and merge.\n"
"        -sgebuild \"options\"    any additional options to sge, e.g.,\n"
"        -sgemerge \"options\"    \"-p -153 -pe thread 2 -A merylaccount\"\n"
"                                 N.B. - -N will be ignored\n"
"                                 N.B. - be sure to quote the options\n"
"\n"
"-M:  Given a list of tables, perform a math, logical or threshold operation.\n"
"     Unless specified, all operations take any number of databases.\n"
"\n"
"     Math operations are:\n"
"        min       count is the minimum count for all databases.  If the mer\n"
"                  does NOT exist in all databases, the mer has a zero count, and\n"
"                  is NOT in the output.\n"
"        minexist  count is the minimum count for all databases that contain the mer\n"
"        max       count is the maximum count for all databases\n"
"        add       count is sum of the counts for all databases\n"
"        sub       count is the first minus the second (binary only)\n"
"        abs       count is the absolute value of the first minus the second (binary only)\n"
"\n"
"     Logical operations are:\n"
"        and       outputs mer iff it exists in all databases\n"
"        nand      outputs mer iff it exists in at least one, but not all, databases\n"
"        or        outputs mer iff it exists in at least one database\n"
"        xor       outputs mer iff it exists in an odd number of databases\n"
"\n"
"     Threshold operations are:\n"
"        lessthan x            outputs mer iff it has count <  x\n"
"        lessthanorequal x     outputs mer iff it has count <= x\n"
"        greaterthan x         outputs mer iff it has count >  x\n"
"        greaterthanorequal x  outputs mer iff it has count >= x\n"
"        equal x               outputs mer iff it has count == x\n"
"     Threshold operations work on exactly one database.\n"
"\n"
"        -s tblprefix  (use tblprefix as a database)\n"
"        -o tblprefix  (create this output)\n"
"        -v            (entertain the user)\n"
"\n"
"     NOTE:  Multiple tables are specified with multiple -s switches; e.g.:\n"
"              %s -M add -s 1 -s 2 -s 3 -s 4 -o all\n"
"     NOTE:  It is NOT possible to specify more than one operation:\n"
"              %s -M add -s 1 -s 2 -sub -s 3\n"
"            will NOT work.\n"
"\n"
"\n"
"-D:  Dump the table (not all of these work).\n"
"\n"
"     -Dt        Dump mers >= a threshold.  Use -n to specify the threshold.\n"
"     -Dc        Count the number of mers, distinct mers and unique mers.\n"
"     -Dh        Dump (to stdout) a histogram of mer counts.\n"
"     -s         Read the count table from here (leave off the .mcdat or .mcidx).\n"
"\n"
"\n";



void
merylArgs::usage(void) {
  fprintf(stderr, usagestring, execName, execName, execName);
}







merylArgs::merylArgs(int argc, char **argv) {
  execName           = duplString(argv[0]);
  options            = 0L;

  beVerbose          = false;

  inputFile          = 0L;
  outputFile         = 0L;

  merSize            = 20;
  merComp            = 0;

  doForward          = true;
  doReverse          = false;
  doCanonical        = false;

  numMersEstimated   = 0;
  numMersActual      = 0;

  numBuckets         = 0;
  numBuckets_log2    = 0;
  mersPerBatch       = 0;
  merDataWidth       = 0;
  merDataMask        = u64bitZERO;
  bucketPointerWidth = 0;
  //bucketPointerMask  = u64bitZERO;

  memoryLimit        = 0;
  segmentLimit       = 0;
  numThreads         = 0;
  configBatch        = false;
  countBatch         = false;
  mergeBatch         = false;
  batchNumber        = 0;

  sgeJobName         = 0L;
  sgeBuildOpt        = 0L;
  sgeMergeOpt        = 0L;

  //  We could possibly do (getenv("SGE_TASK_ID") != 0L), but then we
  //  break on non-SGE grids.
  isOnGrid           = false;

  lowCount           = 0;
  highCount          = ~lowCount;
  desiredCount       = 0;

  mergeFilesMax      = 0;
  mergeFilesLen      = 0;
  mergeFiles         = 0L;

  personality        = 0;

  statsFile          = 0L;

  if (argc == 1) {
    usage();
    exit(1);
  }

  //  Count how many '-s' switches there are, then allocate space
  //  for them in mergeFiles.  We also sum the length of all options,
  //  so we can copy them into an 'options' string used when we
  //  resubmit to the grid.
  //
  u32bit  optionsLen = 0;
  for (int arg=1; arg < argc; arg++) {
    optionsLen += strlen(argv[arg]) + 1;
    if (strcmp(argv[arg], "-s") == 0)
      mergeFilesMax++;
  }

  mergeFiles   = new char * [mergeFilesMax];
  options      = new char   [optionsLen + 1];
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
      merSize = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-c") == 0) {
      arg++;
      merComp = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-s") == 0) {
      arg++;
      delete [] inputFile;
      inputFile                   = duplString(argv[arg]);
      mergeFiles[mergeFilesLen++] = duplString(argv[arg]);
    } else if (strcmp(argv[arg], "-stats") == 0) {
      arg++;
      delete [] statsFile;
      statsFile = duplString(argv[arg]);
    } else if (strcmp(argv[arg], "-n") == 0) {
      arg++;
      numMersEstimated = strtou64bit(argv[arg], 0L);
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
      lowCount = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-U") == 0) {
      arg++;
      highCount = strtou32bit(argv[arg], 0L);
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
      } else if (strcmp(argv[arg], "add") == 0) {
        personality = PERSONALITY_ADD;
      } else if (strcmp(argv[arg], "sub") == 0) {
        personality = PERSONALITY_SUB;
      } else if (strcmp(argv[arg], "abs") == 0) {
        personality = PERSONALITY_ABS;
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
        desiredCount = strtou32bit(argv[arg], 0L) - 1;
      } else if (strcmp(argv[arg], "lessthanorequal") == 0) {
        personality = PERSONALITY_LEQ;
        arg++;
        desiredCount = strtou32bit(argv[arg], 0L);
      } else if (strcmp(argv[arg], "greaterthan") == 0) {
        personality = PERSONALITY_GEQ;
        arg++;
        desiredCount = strtou32bit(argv[arg], 0L) + 1;
      } else if (strcmp(argv[arg], "greaterthanorequal") == 0) {
        personality = PERSONALITY_GEQ;
        arg++;
        desiredCount = strtou32bit(argv[arg], 0L);
      } else if (strcmp(argv[arg], "equal") == 0) {
        personality = PERSONALITY_EQ;
        arg++;
        desiredCount = strtou32bit(argv[arg], 0L);
      } else {
        fprintf(stderr, "ERROR: unknown math personality %s\n", argv[arg]);
        exit(1);
      }
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
      memoryLimit = strtou64bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-segments") == 0) {
      arg++;
      segmentLimit = strtou64bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-threads") == 0) {
      arg++;
      numThreads   = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-configbatch") == 0) {
      personality = 'B';
      configBatch = true;
      countBatch  = false;
      mergeBatch  = false;
      batchNumber = u32bitZERO;
    } else if (strcmp(argv[arg], "-countbatch") == 0) {
      arg++;
      personality = 'B';
      configBatch = false;
      countBatch  = true;
      mergeBatch  = false;
      batchNumber = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-mergebatch") == 0) {
      personality = 'B';
      configBatch = false;
      countBatch  = false;
      mergeBatch  = true;
      batchNumber = u32bitZERO;
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
  char *filename;
  char  magic[17];

  filename = new char [strlen(prefix) + 17];
  sprintf(filename, "%s.merylArgs", prefix);

  FILE *F = fopen(filename, "rb");
  if (errno) {
    fprintf(stderr, "merylArgs::readConfig()-- Failed to open '%s': %s\n", filename, strerror(errno));
    exit(1);
  }

  fread(magic, sizeof(char), 16, F);
  if (strncmp(magic, "merylBatcherv02", 16) != 0) {
    fprintf(stderr, "merylArgs::readConfig()-- '%s' doesn't appear to be a merylArgs file.\n", filename);
    exit(1);
  }

  fread(this, sizeof(merylArgs), 1, F);

  execName    = readString(F);
  inputFile   = readString(F);
  outputFile  = readString(F);
  statsFile   = readString(F);
  sgeJobName  = readString(F);
  sgeBuildOpt = readString(F);
  sgeMergeOpt = readString(F);

  mergeFiles = new char* [mergeFilesLen];
  for (u32bit i=0; i<mergeFilesLen; i++)
    mergeFiles[i] = readString(F);

  fclose(F);

  delete [] filename;
}



merylArgs::~merylArgs() {
  delete [] execName;
  delete [] inputFile;
  delete [] outputFile;

  for (u32bit i=0; i<mergeFilesLen; i++)
    delete [] mergeFiles[i];

  delete [] mergeFiles;
  delete [] statsFile;
}



bool
merylArgs::writeConfig(void) {
  char *filename;

  filename = new char [strlen(outputFile) + 17];
  sprintf(filename, "%s.merylArgs", outputFile);

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
  writeString(statsFile,   F);
  writeString(sgeJobName,  F);
  writeString(sgeBuildOpt, F);
  writeString(sgeMergeOpt, F);

  for (u32bit i=0; i<mergeFilesLen; i++)
    writeString(mergeFiles[i], F);

  fclose(F);

  return(true);
}
