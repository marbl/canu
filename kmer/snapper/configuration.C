#include "posix.H"
#include "snapper2.H"
#include "buildinfo-snapper2.h"
#include "buildinfo-existDB.h"
#include "buildinfo-positionDB.h"
#include "buildinfo-libbri.h"
#include "buildinfo-libsim4.h"
#include "buildinfo-libsim4polish.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

configuration::configuration(void) {

  _beVerbose            = false;

  _merSize              = 20;
  _merSkip              = 0;
  _numSearchThreads     = 4;

  _doReverse            = true;
  _doForward            = true;
  _doValidation         = false;
  _doValidationFileName = 0L;

  _Lo                   = 0.5;
  _Hi                   = 1.0;
  _Va                   = 0.6;

  _maxDiagonal          = 25;

  //  Alternate match extension scheme
  _extendWeight         = 2.0;
  _extendMinimum        = 100;

  _minHitLength         = 0;
  _minHitCoverage       = 0.2;

  _minMatchIdentity     = 98;
  _minMatchCoverage     = 96;

  _dbFileName           = 0L;
  _psFileName           = 0L;
  _qsFileName           = 0L;
  _maskFileName         = 0L;
  _onlyFileName         = 0L;
  _outputFileName       = 0L;
  _statsFileName        = 0L;

  _buildOnly            = false;

  _startTime            = getTime();
  _initTime             = 0.0;
  _buildTime            = 0.0;
  _searchTime           = 0.0;

  _loaderHighWaterMark  = 16 * 1024;
  _loaderSleep.tv_sec   = 1;
  _loaderSleep.tv_nsec  = 0;
  _loaderWarnings       = false;

  _searchSleep.tv_sec   = 0;
  _searchSleep.tv_nsec  = 10000000;

  _writerHighWaterMark  = 32 * 1024;
  _writerSleep.tv_sec   = 1;
  _writerSleep.tv_nsec  = 0;
  _writerWarnings       = false;
}

configuration::~configuration() {
}

static char const *usageString =
"usage: %s [--buildinfo] [options]\n"
"\n"
"Algorithm Options:\n"
"    -forward                Search only the normal cDNA\n"
"    -reverse                Search only the reverse-complement cDNA\n"
"\n"
"    -mersize k              Use k-mers\n"
"    -merskip l              Skip l mers between\n"
"    -maxdiagonal d          Maximum diagonal gap within a hit (25)\n"
"    -minhitlength l         Minimum length for a hit to be polished (0)\n"
"    -minhitcoverage c       Minimum coverage for a hit to be polished (0.2, 0.0 to 1.0)\n"
"    -minmatchidentity i     Minimum percent identity for matches (98, integer)\n"
"    -minmatchcoverage c     Minimum coverage for matches (96, integer)\n"
"    -extendweight w         For each unhit base, extend by this much (2)\n"
"    -extendminimum e        Always extend hits by at least this much (100)\n"
"\n"
"  Filter and Filter Validation\n"
"    -validate               Enable tuning of the filter (expensive!)\n"
"    -setfilter L H V        Use { L,H,V } as the filter parameters\n"
"\n"
"Input Options:\n"
"    -mask f                 Ignore (only use) all mers listed in file f\n"
"    -only f\n"
"    -maskn f n              Ignore (only use) the mers listed in meryl prefix f\n"
"    -onlyn f n              For mask, mers with count >= n are masked.\n"
"                            For only, mers with count <= n are used.\n"
"    -queries c.fasta        Query sequences\n"
"    -genomic g.fasta        Database sequences\n"
"    -positions p.positionDB Build and save / use positionDB.  Assumes you aren't using -use\n"
"    -buildonly              Only do the build and save.\n"
"    -use [...]\n"
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
"Output Options\n"
"    -verbose                Entertain the user\n"
"    -output f               Write output to file f\n"
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
    } else if (strcmp(argv[arg], "-maskn") == 0) {
      arg++;
      _maskPrefix    = argv[arg];
      arg++;
      _maskThreshold = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-onlyn") == 0) {
      arg++;
      _onlyPrefix    = argv[arg];
      arg++;
      _onlyThreshold = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-queries") == 0) {
      arg++;
      _qsFileName = argv[arg];
    } else if (strcmp(argv[arg], "-genomic") == 0) {
      arg++;
      _dbFileName = argv[arg];
    } else if (strcmp(argv[arg], "-positions") == 0) {
      arg++;
      _psFileName = argv[arg];
    } else if (strcmp(argv[arg], "-buildonly") == 0) {
      _buildOnly = argv[arg];
    } else if (strcmp(argv[arg], "-use") == 0) {
      arg++;
      _useList.parse(argv[arg]);
    } else if (strcmp(argv[arg], "-forward") == 0) {
      _doForward = true;
      _doReverse = false;
    } else if (strcmp(argv[arg], "-reverse") == 0) {
      _doReverse = true;
      _doForward = false;
    } else if (strcmp(argv[arg], "-validate") == 0) {
      arg++;
      _doValidation         = true;
      _doValidationFileName = argv[arg];
    } else if ((strcmp(argv[arg], "-setfilter") == 0) ||
               (strcmp(argv[arg], "-lhv") == 0) ||
               (strcmp(argv[arg], "-LHV") == 0)) {
      arg++;  _Lo = atof(argv[arg]);
      arg++;  _Hi = atof(argv[arg]);
      arg++;  _Va = atof(argv[arg]);
    } else if (strcmp(argv[arg], "-verbose") == 0) {
      _beVerbose = true;
    } else if (strcmp(argv[arg], "-output") == 0) {
      arg++;
      _outputFileName = argv[arg];
    } else if (strcmp(argv[arg], "-stats") == 0) {
      arg++;
      _statsFileName = argv[arg];
    } else if (strcmp(argv[arg], "-maxdiagonal") == 0) {
      arg++;
      _maxDiagonal = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-minhitlength") == 0) {
      arg++;
      _minHitLength = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-minhitcoverage") == 0) {
      arg++;
      _minHitCoverage = atof(argv[arg]);
    } else if (strcmp(argv[arg], "-minmatchidentity") == 0) {
      arg++;
      _minMatchIdentity = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-minmatchcoverage") == 0) {
      arg++;
      _minMatchCoverage = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-extendweight", 7) == 0) {
      arg++;
      _extendWeight = atof(argv[arg]);
    } else if (strncmp(argv[arg], "-extendminimum", 7) == 0) {
      arg++;
      _extendMinimum = atoi(argv[arg]);
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
      buildinfo_snapper2(stderr);
#if 0
      buildinfo_existDB(stderr);
      buildinfo_positionDB(stderr);
      buildinfo_libbri(stderr);
      buildinfo_libsim4(stderr);
      buildinfo_libsim4polish(stderr);
#endif
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
    fprintf(stderr, "ERROR:  Number of threads is limited to %d.\n", MAX_THREADS);
    exit(1);
  }

  if (_maskFileName && _onlyFileName) {
    fprintf(stderr, "ERROR:  At most one of -mask and -only may be used.\n");
    exit(-1);
  }

  if (_merSkip >= _merSize) {
    fprintf(stderr, "ERROR:  Mers are not adjacent; make sure merskip <= mersize.\n");
    exit(-1);
  }

  //  Fail if no query sequences
  //
  if ((_qsFileName == 0L) && (_buildOnly == false)) {
    fprintf(stderr, "ERROR: No query file supplied.\n");
    exit(1);
  }

  if (_dbFileName == 0L) {
    fprintf(stderr, "ERROR: No genome file supplied.\n");
    exit(1);
  }
}

void
configuration::display(FILE *out) {
  if ((out == stdout) && (_beVerbose)) {
    fprintf(out, "--Using these Algorithm Options--\n");

    fprintf(out, "merSize             = "u32bitFMT"\n",   _merSize);
    fprintf(out, "merSkip             = "u32bitFMT"\n",   _merSkip);
    fprintf(out, "doReverse           = %s\n",   _doReverse ? "true" : "false");
    fprintf(out, "doForward           = %s\n",   _doForward ? "true" : "false");
    fprintf(out, "\n");
    fprintf(out, "maxDiagonal         = "u32bitFMT"\n",   _maxDiagonal);
    fprintf(out, "minHitLength        = "u32bitFMT"\n",   _minHitLength + _merSize);
    fprintf(out, "minHitCoverage      = %lf\n",  _minHitCoverage);
    fprintf(out, "minMatchIdentity    = "u32bitFMT"\n",   _minMatchIdentity);
    fprintf(out, "minMatchCoverage    = "u32bitFMT"\n",   _minMatchCoverage);
    fprintf(out, "\n");

    if (_doValidation)
      fprintf(out, "--VALIDATION ENABLED--\n\n");

    fprintf(out, "--Using these Process Options--\n");
    fprintf(out, "\n");
    fprintf(out, "numSearchThreads    = "u32bitFMT"\n",   _numSearchThreads);
    fprintf(out, "loaderHighWaterMark = "u32bitFMT"\n", _loaderHighWaterMark);
    fprintf(out, "loaderSleep         = %f\n", (double)_loaderSleep.tv_sec + (double)_loaderSleep.tv_nsec * 1e-9);
    fprintf(out, "loaderWarnings      = %s\n", _loaderWarnings ? "true" : "false");
    fprintf(out, "searchSleep         = %f\n", (double)_searchSleep.tv_sec + (double)_searchSleep.tv_nsec * 1e-9);
    fprintf(out, "writerHighWaterMark = "u32bitFMT"\n", _writerHighWaterMark);
    fprintf(out, "writerSleep         = %f\n", (double)_writerSleep.tv_sec + (double)_writerSleep.tv_nsec * 1e-9);
    fprintf(out, "writerWarnings      = %s\n", _writerWarnings ? "true" : "false");
    fprintf(out, "\n");
    fprintf(out, "--Using the use-list:--\n");
    fprintf(out, "\n");
    fprintf(out, "XXXX:  Need this!\n");
    fprintf(out, "\n");
    fprintf(out, "--Using these Files--\n");
    fprintf(out, "Genomic File          = '%s'\n", _dbFileName);
    fprintf(out, "Positions File        = '%s'\n", _psFileName);
    fprintf(out, "Query File            = '%s'\n", _qsFileName);

    if (_maskFileName)
      fprintf(out, "maskFile            = '%s'\n", _maskFileName);
    else
      fprintf(out, "maskFile            = None Specified.\n");

    if (_onlyFileName)
      fprintf(out, "onlyFile            = '%s'\n", _onlyFileName);
    else
      fprintf(out, "onlyFile            = None Specified.\n");

    if (_outputFileName)
      fprintf(out, "outputFile          = '%s'\n", _outputFileName);
    else
      fprintf(out, "outputFile          = None Specified.\n");

    if (_statsFileName)
      fprintf(out, "statsFile           = '%s'\n", _statsFileName);
    else
      fprintf(out, "statsFile           = None Specified.\n");

    fprintf(out, "\n");
  }
}
