#include "posix.H"
#include "seatac.H"
#include "buildinfo-seatac.h"
#include "buildinfo-libkmer.h"
#include "buildinfo-libbio.h"
#include "buildinfo-libutil.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sharedObj.H"

configuration::configuration(void) {

  _beVerbose           = false;

  _merSize             = 20;
  _merSkip             = 0;
  _numSearchThreads    = 4;

  _doReverse           = true;
  _doForward           = true;

  _maxDiagonal         = 25;
  _maxGap              = 0;
  _qsOverlap           = 15;
  _dsOverlap           = 15;

  _minLength           = 20;

  _dbFileName          = 0L;
  _qsFileName          = 0L;
  _maskFileName        = 0L;
  _onlyFileName        = 0L;
  _outputFileName      = 0L;
  _statsFileName       = 0L;

  _tableFileName       = 0L;

  _useListLen          = 0;
  _useListMax          = 0;
  _useList             = 0L;
  _useListString       = 0L;

  _filtername          = 0L;
  _filteropts          = 0L;
  _filterObj           = 0L;

  _startTime           = 0.0;
  _initTime            = 0.0;
  _buildTime           = 0.0;
  _searchTime          = 0.0;
  _totalTime           = 0.0;

  _loaderHighWaterMark = 2;
  _loaderSleep.tv_sec  = 1;
  _loaderSleep.tv_nsec = 0;
  _loaderWarnings      = false;

  _searchSleep.tv_sec  = 0;
  _searchSleep.tv_nsec = 10000000;

  _writerHighWaterMark = 256;
  _writerSleep.tv_sec  = 1;
  _writerSleep.tv_nsec = 0;
  _writerWarnings      = false;
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
"    -minelength l\n"
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
"    -buildonly datfile      If 'datfile' doesn't exist, build the tables, write\n"
"                            them to 'datfile' and exit.  If 'datfile' exists\n"
"                            AND is a complete and valid file, load the tables\n"
"                            from the file and do the compute.\n"
"\n"
"Filtering Options\n"
"    -filtername x.so        Use the shared object x.so as a filter method.\n"
"    -filteropts opts        The string 'opts' is passed to the filter on creation.\n"
"\n"
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

  _useListLen    = 0;
  _useListMax    = 1024;
  _useList       = new use_s [_useListMax];
  _useListString = strdup(line);

  //  line can be either a list of numbers, or a file.  See if "line" opens
  //  as a file.  If it does, read in the numbers.
  //
  FILE *F = fopen(line, "r");
  if (F) {
    if (_beVerbose)
      fprintf(stderr, "Reading use list from '%s'\n", line);
    while (!feof(F)) {
      fscanf(F, " "u32bitFMT" ", &v);
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



//  Removes duplicate entries, sorts
//
void
configuration::completeUseList(FastAWrapper *db) {

  if (_beVerbose)
    fprintf(stderr, "Completing use list: ("u32bitFMT")\n", _useListLen);

  //  If no use list given, create one using the db.  Otherwise,
  //  fix the existing one.
  //
  if (_useListLen == 0) {
    if (_beVerbose)
      fprintf(stderr, "Using all sequences in the genomic ("u32bitFMT")\n", db->getNumberOfSequences());
    _useListLen = 0;
    _useListMax = db->getNumberOfSequences();
    _useList    = new use_s [_useListMax];

    for (u32bit i=db->getNumberOfSequences(); i--; ) {
      _useList[_useListLen].seq   = i;
      _useList[_useListLen].size  = 0;
      _useList[_useListLen].start = 0;
      _useListLen++;
    }
  } else {
    char    *seen    = new char [db->getNumberOfSequences()];
    u32bit   seenLen = 0;
    for (u32bit i=db->getNumberOfSequences(); i--; )
      seen[i] = 0;

    for (u32bit i=0; i<_useListLen; i++) {
      if (_useList[i].seq >= db->getNumberOfSequences()) {
        fprintf(stderr, "WARNING: Sequence "u32bitFMT" in -use list is not in '%s'\n", i, _dbFileName);
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

    for (u32bit i=db->getNumberOfSequences(); i--; ) {
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
    } else if (strcmp(argv[arg], "-buildonly") == 0) {
      arg++;
      _tableFileName = argv[arg];
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
    } else if (strcmp(argv[arg], "-use") == 0) {
      arg++;
      parseUseLine(argv[arg]);
    } else if (strcmp(argv[arg], "-forward") == 0) {
      _doForward = true;
      _doReverse = false;
    } else if (strcmp(argv[arg], "-reverse") == 0) {
      _doReverse = true;
      _doForward = false;
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
    } else if (strcmp(argv[arg], "-maxgap") == 0) {
      arg++;
      _maxGap = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-qoverlap") == 0) {
      arg++;
      _qsOverlap = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-doverlap") == 0) {
      arg++;
      _dsOverlap = atoi(argv[arg]);
    } else if (strcmp(argv[arg], "-minlength") == 0) {
      arg++;
      _minLength = atoi(argv[arg]);
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
    } else if (strcmp(argv[arg], "-filtername") == 0) {
      arg++;
      _filtername = argv[arg];
      _filterObj  = new sharedObj(argv[arg]);
    } else if (strcmp(argv[arg], "-filteropts") == 0) {
      arg++;
      _filteropts = argv[arg];
    } else if (strncmp(argv[arg], "--buildinfo", 3) == 0) {
      buildinfo_seatac(stderr);
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
    exit(1);
  }

  //
  //  Check that the mers are at least adjacent
  //
  if (_merSkip >= _merSize) {
    fprintf(stderr, "ERROR:  Mers are not adjacent; make sure merskip <= mersize.\n");
    exit(1);
  }

  //
  //  Test that we can build filter and stat objects
  //
  if (_filtername) {
    filterObj *testf = new filterObj(_filterObj, _filteropts);
    delete testf;

    statObj *tests = new statObj(_filterObj, _filteropts);
    delete tests;
  }
}


void
configuration::display(FILE *out) {
  if ((out == stdout) && (_beVerbose)) {
    fprintf(out, "--Using these Options--\n");
    fprintf(out, "beVerbose           = %s\n", _beVerbose ? "enabled" : "disabled");
    fprintf(out, "numSearchThreads    = "u32bitFMT"\n",   _numSearchThreads);
    fprintf(out, "\n");
    fprintf(out, "loaderHighWaterMark = "u32bitFMT"\n", _loaderHighWaterMark);
    fprintf(out, "loaderSleep         = %f\n", (double)_loaderSleep.tv_sec + (double)_loaderSleep.tv_nsec * 1e-9);
    fprintf(out, "loaderWarnings      = %s\n", _loaderWarnings ? "true" : "false");
    fprintf(out, "searchSleep         = %f\n", (double)_searchSleep.tv_sec + (double)_searchSleep.tv_nsec * 1e-9);
    fprintf(out, "writerHighWaterMark = "u32bitFMT"\n", _writerHighWaterMark);
    fprintf(out, "writerSleep         = %f\n", (double)_writerSleep.tv_sec + (double)_writerSleep.tv_nsec * 1e-9);
    fprintf(out, "writerWarnings      = %s\n", _writerWarnings ? "true" : "false");
    fprintf(out, "\n");
    fprintf(out, "--Using these Parameters--\n");
    fprintf(out, "maxDiagonal         = "u32bitFMT"\n",   _maxDiagonal);
    fprintf(out, "maxGap              = "u32bitFMT"\n",   _maxGap);
    fprintf(out, "qsOverlap           = "u32bitFMT"\n",   _qsOverlap);
    fprintf(out, "dsOverlap           = "u32bitFMT"\n",   _dsOverlap);
    fprintf(out, "minLength           = "u32bitFMT"\n",   _minLength + _merSize);
    fprintf(out, "\n");
    fprintf(out, "merSize             = "u32bitFMT"\n",   _merSize);
    fprintf(out, "merSkip             = "u32bitFMT"\n",   _merSkip);
    fprintf(out, "doReverse           = %s\n",   (_doReverse) ? "true" : "false");
    fprintf(out, "doForward           = %s\n",   (_doForward) ? "true" : "false");
    fprintf(out, "\n");
    fprintf(out, "filterName          = %s\n", (_filtername) ? _filtername : "None Specified.");
    fprintf(out, "filterOpts          = %s\n", (_filteropts) ? _filteropts : "None Specified.");
    fprintf(out, "\n");
    fprintf(out, "--Using these Sequences--\n");
    fprintf(out, "useList             = %s\n", (_useListString) ? _useListString : "Every sequence.");
    fprintf(out, "\n");
    fprintf(out, "--Using these Files--\n");
    fprintf(out, "dbFile              = %s\n", (_dbFileName) ? _dbFileName : "None Specified.");
    fprintf(out, "qsFile              = %s\n", (_qsFileName) ? _qsFileName : "None Specified.");
    fprintf(out, "maskFile            = %s\n", (_maskFileName) ? _maskFileName : "None Specified.");
    fprintf(out, "onlyFile            = %s\n", (_onlyFileName) ? _onlyFileName : "None Specified.");
    fprintf(out, "outputFile          = %s\n", (_outputFileName) ? _outputFileName : "None Specified.");
    fprintf(out, "statsFile           = %s\n", (_statsFileName) ? _statsFileName : "None Specified.");
    fprintf(out, "tableFile           = %s\n", (_tableFileName) ? _tableFileName : "None Specified.");
    fprintf(out, "\n");
  }
}



void
configuration::writeATACheader(FILE *out) {
  fprintf(out, "! format atac 1.0\n");
  fprintf(out, "/seatacBeVerbose=%s\n", _beVerbose ? "enabled" : "disabled");
  fprintf(out, "/seatacNumSearchThreads="u32bitFMT"\n", _numSearchThreads);
  fprintf(out, "/seatacLoaderHighWaterMark="u32bitFMT"\n", _loaderHighWaterMark);
  fprintf(out, "/seatacLoaderSleep=%f\n", (double)_loaderSleep.tv_sec + (double)_loaderSleep.tv_nsec * 1e-9);
  fprintf(out, "/seatacLoaderWarnings=%s\n", _loaderWarnings ? "true" : "false");
  fprintf(out, "/seatacSearchSleep=%f\n", (double)_searchSleep.tv_sec + (double)_searchSleep.tv_nsec * 1e-9);
  fprintf(out, "/seatacWriterHighWaterMark="u32bitFMT"\n", _writerHighWaterMark);
  fprintf(out, "/seatacWriterSleep=%f\n", (double)_writerSleep.tv_sec + (double)_writerSleep.tv_nsec * 1e-9);
  fprintf(out, "/seatacWriterWarnings=%s\n", _writerWarnings ? "true" : "false");
  fprintf(out, "/seatacMaxDiagonal="u32bitFMT"\n", _maxDiagonal);
  fprintf(out, "/seatacMaxGap="u32bitFMT"\n", _maxGap);
  fprintf(out, "/seatacQsOverlap="u32bitFMT"\n", _qsOverlap);
  fprintf(out, "/seatacDsOverlap="u32bitFMT"\n", _dsOverlap);
  fprintf(out, "/seatacMinLength="u32bitFMT"\n", _minLength + _merSize);
  fprintf(out, "/seatacMerSize="u32bitFMT"\n", _merSize);
  fprintf(out, "/seatacMerSkip="u32bitFMT"\n", _merSkip);
  fprintf(out, "/seatacDoReverse=%s\n", (_doReverse) ? "true" : "false");
  fprintf(out, "/seatacDoForward=%s\n", (_doForward) ? "true" : "false");
  fprintf(out, "/seatacFilterName=%s\n", (_filtername) ? _filtername : "None Specified.");
  fprintf(out, "/seatacFilterOpts=%s\n", (_filteropts) ? _filteropts : "None Specified.");
  fprintf(out, "/seatacUseList=%s\n", (_useListString) ? _useListString : "Every sequence.");
  fprintf(out, "/seatacDbFile =%s\n", (_dbFileName) ? _dbFileName : "None Specified.");
  fprintf(out, "/seatacQsFile =%s\n", (_qsFileName) ? _qsFileName : "None Specified.");
  fprintf(out, "/seatacMaskFile=%s\n", (_maskFileName) ? _maskFileName : "None Specified.");
  fprintf(out, "/seatacOnlyFile=%s\n", (_onlyFileName) ? _onlyFileName : "None Specified.");
  fprintf(out, "/seatacOutputFile=%s\n", (_outputFileName) ? _outputFileName : "None Specified.");
  fprintf(out, "/seatacStatsFile=%s\n", (_statsFileName) ? _statsFileName : "None Specified.");
  fprintf(out, "/seatacTableFile=%s\n", (_tableFileName) ? _tableFileName : "None Specified.");
}
