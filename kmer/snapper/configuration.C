#include "snapper2.H"

configuration::configuration(void) {

  _beVerbose            = false;

  _merSize              = 20;
  _merSkip              = 0;
  _numSearchThreads     = 4;

  _doReverse            = true;
  _doForward            = true;
  _doValidation         = false;
  _doValidationFileName = 0L;

  _doAlignments         = true;

  _Lo                   = 0.5;
  _Hi                   = 1.0;
  _Va                   = 0.6;

  _maxDiagonal          = 25;

  //  Alternate match extension scheme
  _extendWeight         = 2.0;
  _extendMinimum        = 100;

  _repeatThreshold      = 3;

  _minHitLength         = 0;
  _minHitCoverage       = 0.2;

  _minMatchIdentity     = 98;
  _minMatchCoverage     = 96;

  _afEnabled            = false;
  _afThreshold          = 0.25;
  _afLength             = 64;
  _afInit               = 5;

  _discardExonLength    = 64;
  _discardExonQuality   = 90;
  _splitMatches         = true;

  _dbFileName           = 0L;
  _psFileName           = 0L;
  _qsFileName           = 0L;

  _maskFileName         = 0L;
  _onlyFileName         = 0L;

  _ignoreThreshold      = 0;

  _maskPrefix           = 0L;
  _maskThreshold        = 0;
  _onlyPrefix           = 0L;
  _onlyThreshold        = 0;

  _outputFileName       = 0L;
  _logmsgFileName       = 0L;
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
"usage: %s [options]\n"
"\n"
"Algorithm Options:\n"
"    -forward                Search only the normal cDNA.\n"
"    -reverse                Search only the reverse-complement cDNA.\n"
"\n"
"    -mersize k              Use k-mers.\n"
"    -merskip l              Skip l mers between.\n"
"    -maxdiagonal d          Maximum diagonal gap within a hit (25).\n"
"    -minhitlength l         Minimum length for a hit to be polished (0).\n"
"    -minhitcoverage c       Minimum coverage for a hit to be polished (0.2, 0.0 to 1.0).\n"
"    -minmatchidentity i     Minimum percent identity for matches (98, integer).\n"
"    -minmatchcoverage c     Minimum coverage for matches (96, integer).\n"
"    -discardexonlength l    Discard exons less than l bp long (64).\n"
"    -discardexonquality p   Discard exons less than p percent identity (90).\n"
"    -extendweight w         For each unhit base, extend by this much (2).\n"
"    -extendminimum e        Always extend hits by at least this much (100).\n"
"    -repeatthreshold t      Tune hits to expect t local repeat count (3).\n"
"\n"
"Filter and Filter Validation:\n"
"    -setfilter L H V        Use { L,H,V } as the filter parameters.\n"
"    -validate               Enable tuning of the filter (expensive!).\n"
"\n"
"Masking Options:\n"
"    -ignore n               Ignore mers with count more than n.\n"
"    -mask f                 Ignore (only use) all mers listed in file f.\n"
"    -only f\n"
"    -maskn f n              Ignore (only use) the mers listed in meryl prefix f.\n"
"    -onlyn f n              For mask, mers with count >= n are masked.\n"
"                            For only, mers with count <= n are used.\n"
"\n"
"Input Options:\n"
"    -queries c.fasta        Query sequences.\n"
"    -genomic g.fasta        Database sequences.\n"
"    -positions p.positionDB Build and save / use positionDB.  Assumes you aren't using -use.\n"
"    -buildonly              Only do the build and save.\n"
"    -use [...]\n"
"\n"
"Process Options:\n"
"    -numthreads n           Use n search threads.\n"
"    -loaderhighwatermark h  Size of the loader queue.\n"
"    -loadersleep t          Time the loader will sleep when its output queue is full.\n"
"    -loaderwarnings         Enable warning messages for the loader.\n"
"    -searchsleep t          Time the searcher will sleep when it has no input.\n"
"    -writerhighwatermark h  Size of the output queue.\n"
"    -writersleep t          Time the writer will sleep when it has nothing to write.\n"
"    -writerwarnings         Enable warning messages for the writer.\n"
"\n"
"Output Options:\n"
"    -verbose                Entertain the user with useless statistics.\n"
"    -output f               Write output to file f.\n"
"    -{no}aligns             Enable/Disable full alignments.  Enabled by default.\n"
"    -log f                  Write some debugging/logging information to file f.  This\n"
"                            is mostly for developers, and does NOT provide useful\n"
"                            information unless you know the guts of snapper.\n"
"    -stats f                Write resource usage statistics to f.\n";

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
      _merSize = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-merskip") == 0) {
      arg++;
      _merSkip = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-numthreads") == 0) {
      arg++;
      _numSearchThreads = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-ignore") == 0) {
      ++arg;
      _ignoreThreshold = strtou32bit(argv[arg], 0L);
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
      _maskThreshold = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-onlyn") == 0) {
      arg++;
      _onlyPrefix    = argv[arg];
      arg++;
      _onlyThreshold = strtou32bit(argv[arg], 0L);
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
    } else if (strcmp(argv[arg], "-aligns") == 0) {
      _doAlignments = true;
    } else if (strcmp(argv[arg], "-noaligns") == 0) {
      _doAlignments = false;
    } else if (strcmp(argv[arg], "-log") == 0) {
      arg++;
      _logmsgFileName = argv[arg];
    } else if (strcmp(argv[arg], "-stats") == 0) {
      arg++;
      _statsFileName = argv[arg];
    } else if (strcmp(argv[arg], "-maxdiagonal") == 0) {
      arg++;
      _maxDiagonal = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-minhitlength") == 0) {
      arg++;
      _minHitLength = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-minhitcoverage") == 0) {
      arg++;
      _minHitCoverage = atof(argv[arg]);
    } else if (strcmp(argv[arg], "-minmatchidentity") == 0) {
      arg++;
      _minMatchIdentity = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-minmatchcoverage") == 0) {
      arg++;
      _minMatchCoverage = strtou32bit(argv[arg], 0L);

    } else if (strcmp(argv[arg], "-af") == 0) {
      _afEnabled   = true;
    } else if (strcmp(argv[arg], "-afthreshold") == 0) {
      arg++;
      _afThreshold = atof(argv[arg]);
      _afEnabled   = true;
    } else if (strcmp(argv[arg], "-aflength") == 0) {
      arg++;
      _afLength    = strtou32bit(argv[arg], 0L);
      _afEnabled   = true;
    } else if (strcmp(argv[arg], "-afinit") == 0) {
      arg++;
      _afInit      = strtou32bit(argv[arg], 0L);
      _afEnabled   = true;


    } else if (strcmp(argv[arg], "-discardexonlength") == 0) {
      arg++;
      _discardExonLength = strtou32bit(argv[arg], 0L);
    } else if (strcmp(argv[arg], "-discardexonquality") == 0) {
      arg++;
      _discardExonQuality = strtou32bit(argv[arg], 0L);
    } else if (strncmp(argv[arg], "-extendweight", 8) == 0) {
      arg++;
      _extendWeight = atof(argv[arg]);
    } else if (strncmp(argv[arg], "-extendminimum", 8) == 0) {
      arg++;
      _extendMinimum = strtou32bit(argv[arg], 0L);
    } else if (strncmp(argv[arg], "-repeatthreshold", 8) == 0) {
      arg++;
      _repeatThreshold = strtou32bit(argv[arg], 0L);
    } else if (strncmp(argv[arg], "-loaderhighwatermark", 8) == 0) {
      _loaderHighWaterMark = strtou32bit(argv[++arg], 0L);
    } else if (strncmp(argv[arg], "-loadersleep",         8) == 0) {
      setTime(&_loaderSleep, atof(argv[++arg]));
    } else if (strncmp(argv[arg], "-loaderwarnings",      8) == 0) {
      _loaderWarnings = true;
    } else if (strncmp(argv[arg], "-searchsleep",         8) == 0) {
      setTime(&_searchSleep, atof(argv[++arg]));
    } else if (strncmp(argv[arg], "-writerhighwatermark", 8) == 0) {
      _writerHighWaterMark = strtou32bit(argv[++arg], 0L);
    } else if (strncmp(argv[arg], "-writersleep",         8) == 0) {
      setTime(&_writerSleep, atof(argv[++arg]));
    } else if (strncmp(argv[arg], "-writerwarnings",      8) == 0) {
      _writerWarnings = true;
    } else {
      fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
      exit(1);
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

  if ((_afThreshold < 0) || (_afThreshold > 1.0)) {
    fprintf(stderr, "ERROR: Invalid afThreshold %f, should be 0.0 <= t <= 1.0\n", _afThreshold);
    exit(1);
  }
  if (64 < _afLength) {
    fprintf(stderr, "ERROR: Invalid afLength "u32bitFMT", should be < 64.\n", _afLength);
    exit(1);
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
    fprintf(out, "doAlignments        = %s\n",   _doAlignments ? "true" : "false");
    fprintf(out, "\n");
    fprintf(out, "maxDiagonal         = "u32bitFMT"\n",   _maxDiagonal);
    fprintf(out, "minHitLength        = "u32bitFMT"\n",   _minHitLength + _merSize);
    fprintf(out, "minHitCoverage      = %f\n",  _minHitCoverage);
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
    fprintf(out, "Genomic File        = '%s'\n", _dbFileName);
    fprintf(out, "Positions File      = '%s'\n", _psFileName);
    fprintf(out, "Query File          = '%s'\n", _qsFileName);
    fprintf(out, "maskFile            = '%s'\n", (_maskFileName) ? _maskFileName : "None Specified.");
    fprintf(out, "onlyFile            = '%s'\n", (_onlyFileName) ? _onlyFileName : "None Specified.");
    fprintf(out, "maskPrefix          = '%s'\n", (_maskPrefix) ? _maskPrefix : "None Speficied.");
    fprintf(out, "maskThreshold       = "u32bitFMT"\n", _maskThreshold);
    fprintf(out, "onlyPrefix          = '%s'\n", (_onlyPrefix) ? _onlyPrefix : "None Specified.");
    fprintf(out, "onlyThreshold       = "u32bitFMT"\n", _onlyThreshold);
    fprintf(out, "outputFile          = '%s'\n", (_outputFileName) ? _outputFileName : "None Specified.");
    fprintf(out, "logmsgFile          = '%s'\n", (_logmsgFileName) ? _logmsgFileName : "None Speficied.");
    fprintf(out, "statsFile           = '%s'\n", (_statsFileName) ? _statsFileName : "NoneSpecified.");
    fprintf(out, "\n");
  }
}
