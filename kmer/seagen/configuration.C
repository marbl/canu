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

  _tableTemporaryFileName  = 0L;
  _tableFileName           = 0L;
  _tableBuildOnly          = false;

  _binaryOutput            = false;

  _startTime               = 0.0;
  _initTime                = 0.0;
  _buildTime               = 0.0;
  _searchTime              = 0.0;
  _totalTime               = 0.0;

  _loaderHighWaterMark     = 16 * 1024;
  _loaderSleep.tv_sec      = 1;
  _loaderSleep.tv_nsec     = 0;
  _loaderWarnings          = false;

  _searchSleep.tv_sec      = 0;
  _searchSleep.tv_nsec     = 10000000;

  _writerHighWaterMark     = 32 * 1024;
  _writerSleep.tv_sec      = 1;
  _writerSleep.tv_nsec     = 0;
  _writerWarnings          = false;
}

configuration::~configuration() {
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
"    -loaderhighwatermark h  Size of the loader queue\n"
"    -loadersleep t          Time the loader will sleep when its output queue is full\n"
"    -loaderwarnings         Enable warning messages for the loader\n"
"    -searchsleep t          Time the searcher will sleep when it has no input\n"
"    -writerhighwatermark h  Size of the output queue\n"
"    -writersleep t          Time the writer will sleep when it has nothing to write\n"
"    -writerwarnings         Enable warning messages for the writer\n"
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
    } else if (strcmp(argv[arg], "-use") == 0) {
      arg++;
      config._useList.parse(argv[arg]);
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

    } else if (strncmp(argv[arg], "-loaderhighwatermark", 8) == 0) {
      _loaderHighWaterMark = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-loadersleep",         8) == 0) {
      setTime(&_loaderSleep, atof(argv[++arg]));
    } else if (strncmp(argv[arg], "-loaderwarnings",      8) == 0) {
      _loaderWarnings = true;
    } else if (strncmp(argv[arg], "-searchsleep",         8) == 0) {
      setTime(&_searchSleep, atof(argv[++arg]));
    } else if (strncmp(argv[arg], "-writerhighwatermark", 8) == 0) {
      _writerHighWaterMark = atoi(argv[++arg]);
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

  if (_numSearchThreads > MAX_THREADS) {
    fprintf(stderr, "ERROR:  Threads are limited to %d.\n", MAX_THREADS);
    exit(1);
  }

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
}

void
configuration::display(FILE *out) {
  if ((out == stdout) && (_beVerbose)) {
    fprintf(out, "--Using these Options--\n");
#ifdef TRUE64BIT
    fprintf(out, "numSearchThreads    = %u\n",   _numSearchThreads);
#else
    fprintf(out, "numSearchThreads    = %lu\n",   _numSearchThreads);
#endif

    fprintf(out, "\n");

#ifdef TRUE64BIT
    fprintf(out, "loaderHighWaterMark = %u\n", _loaderHighWaterMark);
#else
    fprintf(out, "loaderHighWaterMark = %lu\n", _loaderHighWaterMark);
#endif
    fprintf(out, "loaderSleep         = %f\n", (double)_loaderSleep.tv_sec + (double)_loaderSleep.tv_nsec * 1e-9);
    fprintf(out, "loaderWarnings      = %s\n", _loaderWarnings ? "true" : "false");
    fprintf(out, "searchSleep         = %f\n", (double)_searchSleep.tv_sec + (double)_searchSleep.tv_nsec * 1e-9);
#ifdef TRUE64BIT
    fprintf(out, "writerHighWaterMark = %u\n", _writerHighWaterMark);
#else
    fprintf(out, "writerHighWaterMark = %lu\n", _writerHighWaterMark);
#endif
    fprintf(out, "writerSleep         = %f\n", (double)_writerSleep.tv_sec + (double)_writerSleep.tv_nsec * 1e-9);
    fprintf(out, "writerWarnings      = %s\n", _writerWarnings ? "true" : "false");

    fprintf(out, "\n");


#ifdef TRUE64BIT
    fprintf(out, "merSize             = %u\n",   _merSize);
    fprintf(out, "merSkip             = %u\n",   _merSkip);
#else
    fprintf(out, "merSize             = %lu\n",   _merSize);
    fprintf(out, "merSkip             = %lu\n",   _merSkip);
#endif
    fprintf(out, "doReverse           = %s\n",   _doReverse ? "true" : "false");
    fprintf(out, "doForward           = %s\n",   _doForward ? "true" : "false");
    fprintf(out, "\n");


    fprintf(out, "--Using these Parameters--\n");
#ifdef TRUE64BIT
    fprintf(out, "maxDiagonal         = %u\n",   _maxDiagonal);
    fprintf(out, "maxGap              = %u\n",   _maxGap);
    fprintf(out, "qsOverlap           = %u\n",   _qsOverlap);
    fprintf(out, "dsOverlap           = %u\n",   _dsOverlap);
    fprintf(out, "maxIntron           = %u\n",   _maxIntronLength);
    fprintf(out, "smallSeqCutoff      = %u\n",   _smallSequenceCutoff);
    fprintf(out, "minLengthSingle     = %u\n",   _minLengthSingle   + _merSize);
    fprintf(out, "minCoverageSingle   = %lf\n",   _minCoverageSingle);
    fprintf(out, "minLengthMultiple   = %u\n",   _minLengthMultiple + _merSize);
    fprintf(out, "minCoverageMultiple = %lf\n",   _minCoverageMultiple);
#else
    fprintf(out, "maxDiagonal         = %lu\n",   _maxDiagonal);
    fprintf(out, "maxGap              = %lu\n",   _maxGap);
    fprintf(out, "qsOverlap           = %lu\n",   _qsOverlap);
    fprintf(out, "dsOverlap           = %lu\n",   _dsOverlap);
    fprintf(out, "maxIntron           = %lu\n",   _maxIntronLength);
    fprintf(out, "smallSeqCutoff      = %lu\n",   _smallSequenceCutoff);
    fprintf(out, "minLengthSingle     = %lu\n",   _minLengthSingle   + _merSize);
    fprintf(out, "minCoverageSingle   = %lf\n",   _minCoverageSingle);
    fprintf(out, "minLengthMultiple   = %lu\n",   _minLengthMultiple + _merSize);
    fprintf(out, "minCoverageMultiple = %lf\n",   _minCoverageMultiple);
#endif
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
