#include "posix.H"
#include "snapper2.H"
#include "buildinfo-snapper2.h"
#include "buildinfo-existDB.h"
#include "buildinfo-positionDB.h"
#include "buildinfo-libbri.h"
#include "buildinfo-libfasta.h"
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

  _L                    = 0.5;
  _H                    = 1.0;
  _V                    = 0.6;

  _maxDiagonal          = 25;

  //  Alternate match extension scheme
  _extendWeight         = 2.0;
  _extendMinimum        = 100;

  _minHitLength         = 0;
  _minHitCoverage       = 0.2;

  _minMatchIdentity     = 98;
  _minMatchCoverage     = 96;

  _dbFileName           = 0L;
  _qsFileName           = 0L;
  _maskFileName         = 0L;
  _onlyFileName         = 0L;
  _outputFileName       = 0L;
  _statsFileName        = 0L;

  _useList              = 0L;
  _useListLen           = 0;
  _useListMax           = 0;

  _startTime            = 0.0;
  _initTime             = 0.0;
  _buildTime            = 0.0;
  _searchTime           = 0.0;
  _totalTime            = 0.0;

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
"    -maxdiagonal d\n"
"    -minhitlength l         Minimum length for a hit to be polished\n"
"    -minhitcoverage c       Minimum coverage for a hit to be polished (0.0 to 1.0)\n"
"    -minmatchidentity i     Minimum percent identity for matches (integer)\n"
"    -minmatchcoverage c     Minimum coverage for matches (integer)\n"
"    -extendweight w         (default is 2)\n"
"    -extendminimum e        (default is 100)\n"
"\n"
"  Filter and Filter Validation\n"
"    -validate               Enable tuning of the filter (expensive!)\n"
"    -setfilter L H V        Use { L,H,V } as the filter parameters\n"
"\n"
"Input Options:\n"
"    -mask f                 Ignore (only use) all mers listed in file f\n"
"    -only f\n"
"    -maskn n f              Ignore (only use) the mers listed in meryl prefix f\n"
"    -onlyn n f              For mask, mers with count >= n are masked.\n"
"                            For only, mers with count <= n are used.\n"
"    -queries c.fasta        Query sequences\n"
"    -genomic g.fasta        Database sequences\n"
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
configuration::addToUse(u32bit v) {

  if (_useListLen >= _useListMax) {
    _useListMax <<= 1;
    use_s *u = new use_s [_useListMax];
    for (u32bit i=0; i<_useListLen; i++) {
      u[i].seq  = _useList[i].seq;
      u[i].size = _useList[i].size;
    }
    delete [] _useList;
    _useList = u;
  }

  _useList[_useListLen].seq  = v;
  _useList[_useListLen].size = 0;
  _useListLen++;
}


void
configuration::parseUseLine(char *line) {
  u32bit  v=0, u=0;
  char   *rest = line;

  _useListLen = 0;
  _useListMax = 1024;
  _useList    = new use_s [_useListMax];


  //  line can be either a list of numbers, or a file.  See if "line" opens
  //  as a file.  If it does, read in the numbers.
  //
  FILE *F = fopen(line, "r");
  if (F) {
    if (_beVerbose)
      fprintf(stderr, "Reading use list from '%s'\n", line);
    while (!feof(F)) {
#ifdef TRUE64BIT
      fscanf(F, " %u ", &v);
#else
      fscanf(F, " %lu ", &v);
#endif
      addToUse(v);
    }
    fclose(F);
  } else {
    while (*line) {

      //  We are at a number.  Get it.
      //
      v = (u32bit)strtoul(line, &rest, 10);
      line = rest;
      
      //  If we get ',' or 0, add the single number to the list.
      //  If we get '-', decode the rest of the range.
      //
      switch (*line) {
        case 0:
        case ',':
          addToUse(v);
          break;
        case '-':
          line++;
          u = (u32bit)strtoul(line, &rest, 10);
          line = rest;

          if (v > u)
            for (; u <= v; u++)
              addToUse(u);
          else
            for (; v <= u; v++)
              addToUse(v);
          break;
        default:
          fprintf(stderr, "Invalid -use specification -- got number, but no separator.\n");
          fprintf(stderr, "Trouble starts at '%s'\n", line);
          exit(1);
          break;
      }

      //  We should be at a ',' (or end of line).  Anything else is
      //  an error.  If ',', move to the next number.
      //
      switch (*line) {
        case 0:
          break;
        case ',':
          line++;
          break;
        default:
          fprintf(stderr, "Invalid -use specification -- didn't get separator.\n");
          fprintf(stderr, "Trouble starts at '%s'\n", line);
          exit(1);
          break;
      }
    }
  }
}


static
int
useListSortHelper(const void *a, const void *b) {
  use_s const *A = (use_s const *)a;
  use_s const *B = (use_s const *)b;

  if (A->seq < B->seq)
    return(-1);
  if (A->seq > B->seq)
    return(1);

  return(0);
}


#ifdef TRUE64BIT
char const *useListMessage       = "Completing use list: (%u)\n";
char const *useAllGenomicMessage = "Using all sequences in the genomic (%u)\n";
char const *notInMessage         = "WARNING: Sequence %u in -use list is not in '%s'\n";
#else
char const *useListMessage       = "Completing use list: (%lu)\n";
char const *useAllGenomicMessage = "Using all sequences in the genomic (%lu)\n";
char const *notInMessage         = "WARNING: Sequence %lu in -use list is not in '%s'\n";
#endif

//  Removes duplicate entries, sorts
//
void
configuration::completeUseList(u32bit numSeqs) {

  if (_beVerbose)
    fprintf(stderr, useListMessage, _useListLen);

  //  If no use list given, create one using the db.  Otherwise,
  //  fix the existing one.
  //
  if (_useListLen == 0) {
    if (_beVerbose)
      fprintf(stderr, useAllGenomicMessage, numSeqs);
    _useListLen = 0;
    _useListMax = numSeqs;
    _useList    = new use_s [_useListMax];

    for (u32bit i=numSeqs; i--; ) {
      _useList[_useListLen].seq   = i;
      _useList[_useListLen].size  = 0;
      _useList[_useListLen].start = 0;
      _useListLen++;
    }
  } else {
    char    *seen    = new char [numSeqs];
    u32bit   seenLen = 0;
    for (u32bit i=numSeqs; i--; )
      seen[i] = 0;

    for (u32bit i=0; i<_useListLen; i++) {
      if (_useList[i].seq >= numSeqs) {
        fprintf(stderr, notInMessage, i, _dbFileName);
      } else {
        if (seen[_useList[i].seq] == 0) {
          seen[_useList[i].seq] = 1;
          seenLen++;
        }
      }
    }

    delete [] _useList;

    _useListLen = 0;
    _useListMax = seenLen;
    _useList    = new use_s [_useListMax];

    for (u32bit i=numSeqs; i--; ) {
      if (seen[i]) {
        _useList[_useListLen].seq   = i;
        _useList[_useListLen].size  = 0;
        _useList[_useListLen].start = 0;

        _useListLen++;
      }
    }
  }

  qsort(_useList, _useListLen, sizeof(use_s), useListSortHelper);
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
    } else if (strcmp(argv[arg], "-use") == 0) {
      arg++;
      parseUseLine(argv[arg]);
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
      arg++;  _L = atof(argv[arg]);
      arg++;  _H = atof(argv[arg]);
      arg++;  _V = atof(argv[arg]);
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
      buildinfo_libfasta(stderr);
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
  if (_qsFileName == 0L) {
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
#ifdef TRUE64BIT
    fprintf(out, "maxDiagonal         = %u\n",   _maxDiagonal);
    fprintf(out, "minHitLength        = %u\n",   _minHitLength + _merSize);
    fprintf(out, "minHitCoverage      = %lf\n",  _minHitCoverage);
    fprintf(out, "minMatchIdentity    = %u\n",   _minMatchIdentity);
    fprintf(out, "minMatchCoverage    = %u\n",   _minMatchCoverage);
#else
    fprintf(out, "maxDiagonal         = %lu\n",   _maxDiagonal);
    fprintf(out, "minHitLength        = %lu\n",   _minHitLength + _merSize);
    fprintf(out, "minHitCoverage      = %lf\n",   _minHitCoverage);
    fprintf(out, "minMatchIdentity    = %lu\n",   _minMatchIdentity);
    fprintf(out, "minMatchCoverage    = %lu\n",   _minMatchCoverage);
#endif

    fprintf(out, "\n");

    if (_doValidation) {
      fprintf(out, "--VALIDATION ENABLED--\n");
      fprintf(out, "\n");
    }

    fprintf(out, "--Using these Process Options--\n");
    fprintf(out, "\n");
#ifdef TRUE64BIT
    fprintf(out, "numSearchThreads    = %u\n",   _numSearchThreads);
#else
    fprintf(out, "numSearchThreads    = %lu\n",   _numSearchThreads);
#endif

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
    fprintf(out, "--Using the use-list:--\n");
    fprintf(out, "\n");
    fprintf(out, "XXXX:  Need this!\n");
    fprintf(out, "\n");

    fprintf(out, "\n");
    fprintf(out, "--Using these Files--\n");
    fprintf(out, "Genomic File          = '%s'\n", _dbFileName);
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
