#include "posix.H"
#include "searchGENOME.H"
#include "buildinfo-seagen.h"
#include "buildinfo-libkmer.h"
#include "buildinfo-libbio.h"
#include "buildinfo-libutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

configuration::configuration(void) {

  _beVerbose               = false;

  _merSize                 = 20;
  _merSkip                 = 0;
  _numSearchThreads        = 4;

  _doReverse               = true;
  _doForward               = true;

  _maxDiagonal             = 25;
  _maxGap                  = 0;
  _qsOverlap               = 15;
  _dsOverlap               = 15;

  //  Alternate match extension scheme
  _extendWeight            = 0;
  _extendMinimum           = 0;
  _extendAlternate         = false;

  _maxIntronLength         = 1000000000;

  _smallSequenceCutoff     = 0;

  _minLengthSingle         = 0;
  _minCoverageSingle       = 0.0;
  _minLengthMultiple       = 0;
  _minCoverageMultiple     = 0.0;

  _dbFileName              = 0L;
  _qsFileName              = 0L;
  _maskFileName            = 0L;
  _onlyFileName            = 0L;
  _outputFileName          = 0L;
  _queryMatchFileName      = 0L;
  _statsFileName           = 0L;

  _outputFile              = STDOUT_FILENO;
  _matchCountsFile         = -1;

  _tableTemporaryFileName  = 0L;
  _tableFileName           = 0L;
  _tableBuildOnly          = false;

  _binaryOutput            = false;

  _qsFASTA                 = 0L;
  _maskDB                  = 0L;
  _onlyDB                  = 0L;
  _positions               = 0L;

  _numberOfQueries         = 0;

  _startTime               = getTime();
  _initTime                = _startTime;
  _buildTime               = _startTime;
  _searchTime              = _startTime;

  _loaderQueue             = 16 * 1024;
  _loaderSleep.tv_sec      = 1;
  _loaderSleep.tv_nsec     = 0;
  _loaderWarnings          = false;

  _searchSleep.tv_sec      = 0;
  _searchSleep.tv_nsec     = 10000000;

  _writerQueue             = 32 * 1024;
  _writerSleep.tv_sec      = 1;
  _writerSleep.tv_nsec     = 0;
  _writerWarnings          = false;
}

configuration::~configuration() {

  if (_beVerbose) {
    u32bit nq = _qsFASTA->getNumberOfSequences();
    double tm = _searchTime - _buildTime;
    fprintf(stderr, "\n"u32bitFMTW(7)" sequences in %5.2f seconds, %8.3f per second.\n", nq, tm, nq/tm);
  }

  errno = 0;
  close(_outputFile);
  close(_matchCountsFile);
  if (errno)
    fprintf(stderr, "Couldn't close to the output file '%s': %s\n", config._outputFileName, strerror(errno));

  delete _qsFASTA;
  delete _maskDB;
  delete _onlyDB;
  delete _positions;

  dumpStats();
}

static char const *usageString =
"usage: %s [--buildinfo] [options]\n"
"\n"
"Algorithm Options:\n"
"    -mersize k              Use k-mers\n"
"    -merskip j              Skip j mers between each mer inserted into table\n"
"    -forward                Search only the normal query sequences\n"
"    -reverse                Search only the reverse-complemented query sequences\n"
"    -maxdiagonal d\n"
"    -maxgap g\n"
"    -qoverlap q\n"
"    -doverlap d\n"
"    -maxintron m\n"
"    -smallsequence\n"
"    -singlelength l\n"
"    -singlecoverage c\n"
"    -multiplelength l\n"
"    -multiplecoverage c\n"
"    -extendweight w\n"
"    -extendminimum m\n"
"\n"
"Process Options\n"
"    -numthreads n           Use n search threads\n"
"\n"
"    -loaderqueue h          Size of the loader queue\n"
"    -loadersleep t          Time the loader will sleep when its output queue is full\n"
"    -loaderwarnings         Enable warning messages for the loader\n"
"\n"
"    -searchsleep t          Time the searcher will sleep when it has no input\n"
"\n"
"    -writerqueue h          Size of the output queue\n"
"    -writersleep t          Time the writer will sleep when it has nothing to write\n"
"    -writerwarnings         Enable warning messages for the writer\n"
"\n"
"\n"
"    -buildtables datfile    If 'datfile' doesn't exist, build the tables, write\n"
"                            them to 'datfile' and exit.\n"
"    -usetables datfile      Load the tables from 'datfile' file and do the compute.\n"
"                            If 'datfile' doesn't exist, an implicit -buildtables is\n"
"                            performed.\n"
"Input Options:\n"
"    -mask f                 Ignore all mers listed in file f\n"
"    -only f                 Use only the mers listed in file f\n"
"    -cdna c.fasta           Query sequences (the cDNA, the stream)\n"
"    -stream                 An alias for -cdna\n"
"    -genomic g.fasta        Database sequences (the genome, the table)\n"
"    -table                  An alias for -genomic)\n"
"    -use #,#,#,#            using only those sequences specified\n"
"    -use file               using only those sequences listed in the file\n"
"\n"
"Output Options\n"
"    -verbose                Entertain the user\n"
"    -binary                 Write the hits in a binary format\n"
"    -output f               Write output to file f\n"
"    -count f                Write counts of hits to file f\n"
"    -stats f                Write resource statistics to f\n";

void
configuration::usage(char *name) {
  fprintf(stderr, usageString, name);
}



void
configuration::read(int argc, char **argv) {
  int arg = 1;

  if (argc < 2) {
    usage(argv[0]);
    exit(1);
  }

  while (arg < argc) {
    if        (strcmp(argv[arg], "-mersize") == 0) {
      arg++;
      _merSize = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-merskip") == 0) {
      arg++;
      _merSkip = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-numthreads") == 0) {
      arg++;
      _numSearchThreads = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-mask") == 0) {
      arg++;
      _maskFileName = argv[arg];
    } else if (strcmp(argv[arg], "-only") == 0) {
      arg++;
      _onlyFileName = argv[arg];
    } else if (strcmp(argv[arg], "-cdna") == 0) {
      arg++;
      _qsFileName = argv[arg];
    } else if (strcmp(argv[arg], "-stream") == 0) {
      arg++;
      _qsFileName = argv[arg];
    } else if (strcmp(argv[arg], "-genomic") == 0) {
      arg++;
      _dbFileName = argv[arg];
    } else if (strcmp(argv[arg], "-table") == 0) {
      arg++;
      _dbFileName = argv[arg];
    } else if (strcmp(argv[arg], "-buildtemporary") == 0) {
      arg++;
      _tableTemporaryFileName = argv[arg];
    } else if (strcmp(argv[arg], "-buildtables") == 0) {
      arg++;
      _tableFileName  = argv[arg];
      _tableBuildOnly = true;
    } else if (strcmp(argv[arg], "-usetables") == 0) {
      arg++;
      _tableFileName  = argv[arg];
      _tableBuildOnly = false;
    } else if (strcmp(argv[arg], "-positions") == 0) {
      arg++;
      _tableFileName  = argv[arg];
      _tableBuildOnly = false;
    } else if (strcmp(argv[arg], "-use") == 0) {
      arg++;
      _useList.parse(argv[arg]);
    } else if (strcmp(argv[arg], "-forward") == 0) {
      _doForward = true;
      _doReverse = false;
    } else if (strcmp(argv[arg], "-reverse") == 0) {
      _doReverse = true;
      _doForward = false;
    } else if (strcmp(argv[arg], "-verbose") == 0) {
      _beVerbose = true;
    } else if (strcmp(argv[arg], "-binary") == 0) {
      _binaryOutput = true;
    } else if (strcmp(argv[arg], "-output") == 0) {
      arg++;
      _outputFileName = argv[arg];
    } else if (strcmp(argv[arg], "-count") == 0) {
      arg++;
      _queryMatchFileName = argv[arg];
    } else if (strcmp(argv[arg], "-stats") == 0) {
      arg++;
      _statsFileName = argv[arg];
    } else if (strcmp(argv[arg], "-maxdiagonal") == 0) {
      arg++;
      _maxDiagonal = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-maxgap") == 0) {
      arg++;
      _maxGap = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-qoverlap") == 0) {
      arg++;
      _qsOverlap = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-doverlap") == 0) {
      arg++;
      _dsOverlap = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-maxintron") == 0) {
      arg++;
      _maxIntronLength = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-smallsequence") == 0) {
      arg++;
      _smallSequenceCutoff = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-singlelength") == 0) {
      arg++;
      _minLengthSingle = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-multiplelength") == 0) {
      arg++;
      _minLengthMultiple = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-singlecoverage") == 0) {
      arg++;
      _minCoverageSingle = atof(argv[arg]);
    } else if (strcmp(argv[arg], "-multiplecoverage") == 0) {
      arg++;
      _minCoverageMultiple = atof(argv[arg]);
    } else if (strncmp(argv[arg], "-extendweight", 7) == 0) {
      arg++;
      _extendWeight = atoi(argv[arg]);
      _extendAlternate = true;
    } else if (strncmp(argv[arg], "-extendminimum", 7) == 0) {
      arg++;
      _extendMinimum = atoi(argv[arg]);
      _extendAlternate = true;

    } else if (strncmp(argv[arg], "-loaderqueue", 8) == 0) {
      _loaderQueue = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-loadersleep",         8) == 0) {
      setTime(&_loaderSleep, atof(argv[++arg]));
    } else if (strncmp(argv[arg], "-loaderwarnings",      8) == 0) {
      _loaderWarnings = true;

    } else if (strncmp(argv[arg], "-searchsleep",         8) == 0) {
      setTime(&_searchSleep, atof(argv[++arg]));

    } else if (strncmp(argv[arg], "-writerqueue", 8) == 0) {
      _writerQueue = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-writersleep",         8) == 0) {
      setTime(&_writerSleep, atof(argv[++arg]));
    } else if (strncmp(argv[arg], "-writerwarnings",      8) == 0) {
      _writerWarnings = true;

    } else if (strncmp(argv[arg], "--buildinfo", 3) == 0) {
      buildinfo_seagen(stderr);
      buildinfo_libkmer(stderr);
      buildinfo_libbio(stderr);
      buildinfo_libutil(stderr);
      exit(1);
    } else {
      fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  //
  //  Make sure some constraints are met
  //

  if (_maskFileName && _onlyFileName) {
    fprintf(stderr, "ERROR:  At most one of -mask and -only may be used.\n");
    exit(-1);
  }

  //
  //  Check that the mers are at least adjacent
  //
  if (_merSkip >= _merSize) {
    fprintf(stderr, "ERROR:  Mers are not adjacent; make sure merskip <= mersize.\n");
    exit(-1);
  }

  //  Fail if we don't get reasonable signal criteria
  //
  if (((_minLengthSingle   == 0) && (_minCoverageSingle   == 0.0)) ||
      ((_minLengthMultiple == 0) && (_minCoverageMultiple == 0.0)))
    fprintf(stderr, "WARNING:  Minimum match lengths not specified.  All matches will be reported.\n");


  //  Open output file
  //
  if (_outputFileName) {
    errno = 0;
    _outputFile = open(_outputFileName,
                      O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                      S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno) {
      fprintf(stderr, "Couldn't open the output file '%s'?\n%s\n", _outputFileName, strerror(errno));
      exit(1);
    }
  }


  if (_queryMatchFileName) {
    errno = 0;
    _matchCountsFile = open(_queryMatchFileName,
                            O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                            S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno) {
      fprintf(stderr, "Couldn't open the match counts file '%s'?\n%s\n", _queryMatchFileName, strerror(errno));
      exit(1);
    }
  }


  //  Gotta go somewhere!
  //
  _startTime = getTime();
}




void
configuration::display(FILE *out) {
  if ((out == stdout) && (_beVerbose)) {
    fprintf(out, "--Using these Options--\n");
    fprintf(out, "numSearchThreads    = "u32bitFMT"\n",   _numSearchThreads);

    fprintf(out, "\n");

    fprintf(out, "loaderQueue         = "u32bitFMT"\n", _loaderQueue);
    fprintf(out, "loaderSleep         = %f\n", (double)_loaderSleep.tv_sec + (double)_loaderSleep.tv_nsec * 1e-9);
    fprintf(out, "loaderWarnings      = %s\n", _loaderWarnings ? "true" : "false");
    fprintf(out, "searchSleep         = %f\n", (double)_searchSleep.tv_sec + (double)_searchSleep.tv_nsec * 1e-9);
    fprintf(out, "writerQueue         = "u32bitFMT"\n", _writerQueue);
    fprintf(out, "writerSleep         = %f\n", (double)_writerSleep.tv_sec + (double)_writerSleep.tv_nsec * 1e-9);
    fprintf(out, "writerWarnings      = %s\n", _writerWarnings ? "true" : "false");

    fprintf(out, "\n");


    fprintf(out, "merSize             = "u32bitFMT"\n",   _merSize);
    fprintf(out, "merSkip             = "u32bitFMT"\n",   _merSkip);
    fprintf(out, "doReverse           = %s\n",   _doReverse ? "true" : "false");
    fprintf(out, "doForward           = %s\n",   _doForward ? "true" : "false");
    fprintf(out, "\n");


    fprintf(out, "--Using these Parameters--\n");
    fprintf(out, "maxDiagonal         = "u32bitFMT"\n",   _maxDiagonal);
    fprintf(out, "maxGap              = "u32bitFMT"\n",   _maxGap);
    fprintf(out, "qsOverlap           = "u32bitFMT"\n",   _qsOverlap);
    fprintf(out, "dsOverlap           = "u32bitFMT"\n",   _dsOverlap);
    fprintf(out, "maxIntron           = "u32bitFMT"\n",   _maxIntronLength);
    fprintf(out, "smallSeqCutoff      = "u32bitFMT"\n",   _smallSequenceCutoff);
    fprintf(out, "minLengthSingle     = "u32bitFMT"\n",   _minLengthSingle   + _merSize);
    fprintf(out, "minCoverageSingle   = %f\n",   _minCoverageSingle);
    fprintf(out, "minLengthMultiple   = "u32bitFMT"\n",   _minLengthMultiple + _merSize);
    fprintf(out, "minCoverageMultiple = %f\n",   _minCoverageMultiple);
    fprintf(out, "\n");
    fprintf(out, "--Using these Files--\n");
    if (_dbFileName)
      fprintf(out, "dbFile              = '%s'\n", _dbFileName);
    else
      fprintf(out, "dbFile              = None Specified.\n");
      
    if (_outputFileName)
      fprintf(out, "outputFile          = '%s'\n", _outputFileName);
    else
      fprintf(out, "outputFile          = None Specified.\n");

    if (_queryMatchFileName)
      fprintf(out, "countFile           = '%s'\n", _queryMatchFileName);
    else
      fprintf(out, "countFile           = None Specified.\n");

    if (_statsFileName)
      fprintf(out, "statsFile           = '%s'\n", _statsFileName);
    else
      fprintf(out, "statsFile           = None Specified.\n");

    if (_qsFileName)
      fprintf(out, "qsFile              = '%s'\n", _qsFileName);
    else
      fprintf(out, "qsFile              = None Specified.\n");

    if (_maskFileName)
      fprintf(out, "maskFile            = '%s'\n", _maskFileName);
    else
      fprintf(out, "maskFile            = None Specified.\n");

    if (_onlyFileName)
      fprintf(out, "onlyFile            = '%s'\n", _onlyFileName);
    else
      fprintf(out, "onlyFile            = None Specified.\n");

    fprintf(out, "\n");
  }
}




void
configuration::dumpStats(void) {

  if (_statsFileName) {
    errno = 0;
    FILE *F = fopen(_statsFileName, "wb");
    if (errno) {
      fprintf(stderr, "Couldn't open the stats file '%s': %s\nStats going to stderr.", _statsFileName, strerror(errno));
      _statsFileName = 0L;
      F = stderr;
    }

    display(F);
    write_rusage(F);
      
    fprintf(F, "wallClockTimes--------------------------\n");
    fprintf(F, "init:     %9.5f\n", _initTime   - _startTime);
    fprintf(F, "build:    %9.5f\n", _buildTime  - _initTime);
    fprintf(F, "search:   %9.5f\n", _searchTime - _buildTime);
    fprintf(F, "total:    %9.5f\n", getTime()   - _startTime);

#if 0
    fprintf(F, "searchThreadInfo------------------------\n");
    for (u64bit i=0; i<_numSearchThreads; i++)
      if (threadStats[i])
        fprintf(F, threadStats[i]);
#endif

    if (_statsFileName)
      fclose(F);
  }


}
