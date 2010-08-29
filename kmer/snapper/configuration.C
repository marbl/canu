#include "snapper2.H"



configuration::configuration() {
  _beVerbose            = false;

  _KBmerSize            = 20;
  _KBcompression        = 0;
  _KBspacingTemplate    = 0L;

  _merSkip              = 0;

  _numSearchThreads     = 4;

  _doReverse            = true;
  _doForward            = true;
  _doValidation         = false;
  _doValidationFileName = 0L;

  _doAlignments         = false;

  _Lo                   = 0.5;
  _Hi                   = 1.0;
  _Va                   = 0.6;

  _maxDiagonal          = 25;

  //  Alternate match extension scheme
  _extendWeight         = 2.0;
  _extendMinimum        = 100;
  _extendMaximum        = 2000;

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
  _polishOptimally      = false;

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
}

configuration::~configuration() {
}




void
configuration::read(int argc, char **argv) {

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-mersize") == 0) {
      _KBmerSize = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-merskip") == 0) {
      _merSkip = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-compression") == 0) {
      _KBcompression = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-template") == 0) {
      _KBspacingTemplate = argv[++arg];
    } else if (strcmp(argv[arg], "-numthreads") == 0) {
      _numSearchThreads = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-ignore") == 0) {
      _ignoreThreshold = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-mask") == 0) {
      _maskFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-only") == 0) {
      _onlyFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-maskn") == 0) {
      _maskPrefix    = argv[++arg];
      _maskThreshold = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-onlyn") == 0) {
      _onlyPrefix    = argv[++arg];
      _onlyThreshold = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-queries") == 0) {
      _qsFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-genomic") == 0) {
      _dbFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-positions") == 0) {
      _psFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-buildonly") == 0) {
      _buildOnly = argv[arg];
    } else if (strcmp(argv[arg], "-forward") == 0) {
      _doForward = true;
      _doReverse = false;
    } else if (strcmp(argv[arg], "-reverse") == 0) {
      _doReverse = true;
      _doForward = false;
    } else if (strcmp(argv[arg], "-validate") == 0) {
      _doValidation         = true;
      _doValidationFileName = argv[++arg];
    } else if ((strcmp(argv[arg], "-setfilter") == 0) ||
               (strcmp(argv[arg], "-lhv") == 0) ||
               (strcmp(argv[arg], "-LHV") == 0)) {
      _Lo = atof(argv[++arg]);
      _Hi = atof(argv[++arg]);
      _Va = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-verbose") == 0) {
      _beVerbose = true;
    } else if (strcmp(argv[arg], "-output") == 0) {
      _outputFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-aligns") == 0) {
      _doAlignments = true;
    } else if (strcmp(argv[arg], "-noaligns") == 0) {
      _doAlignments = false;
    } else if (strcmp(argv[arg], "-log") == 0) {
      _logmsgFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-stats") == 0) {
      _statsFileName = argv[++arg];
    } else if (strcmp(argv[arg], "-dp") == 0) {
      _polishOptimally = true;
    } else if (strcmp(argv[arg], "-maxdiagonal") == 0) {
      _maxDiagonal = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-minhitlength") == 0) {
      _minHitLength = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-minhitcoverage") == 0) {
      _minHitCoverage = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-minmatchidentity") == 0) {
      _minMatchIdentity = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-minmatchcoverage") == 0) {
      _minMatchCoverage = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-af") == 0) {
      _afEnabled   = true;
    } else if (strcmp(argv[arg], "-afthreshold") == 0) {
      _afThreshold = atof(argv[++arg]);
      _afEnabled   = true;
    } else if (strcmp(argv[arg], "-aflength") == 0) {
      _afLength    = strtou32bit(argv[++arg], 0L);
      _afEnabled   = true;
    } else if (strcmp(argv[arg], "-afinit") == 0) {
      _afInit      = strtou32bit(argv[++arg], 0L);
      _afEnabled   = true;

    } else if (strcmp(argv[arg], "-discardexonlength") == 0) {
      _discardExonLength = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-discardexonquality") == 0) {
      _discardExonQuality = strtou32bit(argv[++arg], 0L);
    } else if (strncmp(argv[arg], "-extendweight", 8) == 0) {
      _extendWeight = atof(argv[++arg]);
    } else if (strncmp(argv[arg], "-extendminimum", 8) == 0) {
      _extendMinimum = strtou32bit(argv[++arg], 0L);
    } else if (strncmp(argv[arg], "-extendmaximum", 8) == 0) {
      _extendMaximum = strtou32bit(argv[++arg], 0L);
    } else if (strncmp(argv[arg], "-repeatthreshold", 8) == 0) {
      _repeatThreshold = strtou32bit(argv[++arg], 0L);

    } else {
      fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  //
  //  Make sure some constraints are met
  //

  if (_maskFileName && _onlyFileName)
    fprintf(stderr, "ERROR:  At most one of -mask and -only may be used.\n"), err++;

  if (_merSkip >= _KBmerSize)
    fprintf(stderr, "ERROR:  Mers are not adjacent; make sure merskip <= mersize.\n"), err++;

  if ((_KBcompression) || (_KBspacingTemplate))
    fprintf(stderr, "ERROR:  Mer compression and spacing not supported right now.  :-(\n"), err++;

  if ((_afThreshold < 0) || (_afThreshold > 1.0))
    fprintf(stderr, "ERROR: Invalid afThreshold %f, should be 0.0 <= t <= 1.0\n", _afThreshold), err++;

  if (64 < _afLength)
    fprintf(stderr, "ERROR: Invalid afLength "u32bitFMT", should be < 64.\n", _afLength), err++;

  if ((_qsFileName == 0L) && (_buildOnly == false))
    fprintf(stderr, "ERROR: No query file supplied.\n"), err++;

  if (_dbFileName == 0L)
    fprintf(stderr, "ERROR: No genome file supplied.\n"), err++;

  //
  //  Be helpful.
  //

  if (err) {
    fprintf(stderr, "usage: %s [options]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Algorithm Options:\n");
    fprintf(stderr, "    -forward                Search only the normal cDNA.\n");
    fprintf(stderr, "    -reverse                Search only the reverse-complement cDNA.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -mersize k              Use k-mers.\n");
    fprintf(stderr, "    -merskip l              Skip l mers between.\n");
    fprintf(stderr, "    -compression c          Compress homopolymer runs to c letters.\n");
    fprintf(stderr, "    -template t             Use spaced seed template t.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -dp                     Optimially polish (broken)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -maxdiagonal d          Maximum diagonal gap within a hit (25).\n");
    fprintf(stderr, "    -minhitlength l         Minimum length for a hit to be polished (0).\n");
    fprintf(stderr, "    -minhitcoverage c       Minimum coverage for a hit to be polished (0.2, 0.0 to 1.0).\n");
    fprintf(stderr, "    -minmatchidentity i     Minimum percent identity for matches (98, integer).\n");
    fprintf(stderr, "    -minmatchcoverage c     Minimum coverage for matches (96, integer).\n");
    fprintf(stderr, "    -discardexonlength l    Discard exons less than l bp long (64).\n");
    fprintf(stderr, "    -discardexonquality p   Discard exons less than p percent identity (90).\n");
    fprintf(stderr, "    -extendweight w         For each unhit base, extend by this much (2).\n");
    fprintf(stderr, "    -extendminimum e        Extend hits by at least this much (100).\n");
    fprintf(stderr, "    -extendmaximum e        Extend hits by at most this much (2000).\n");
    fprintf(stderr, "    -repeatthreshold t      Tune hits to expect t local repeat count (3).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Filter and Filter Validation:\n");
    fprintf(stderr, "    -setfilter L H V        Use { L,H,V } as the filter parameters.\n");
    fprintf(stderr, "    -validate               Enable tuning of the filter (expensive!).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Masking Options:\n");
    fprintf(stderr, "    -ignore n               Ignore mers with count more than n.\n");
    fprintf(stderr, "    -mask f                 Ignore (only use) all mers listed in file f.\n");
    fprintf(stderr, "    -only f\n");
    fprintf(stderr, "    -maskn f n              Ignore (only use) the mers listed in meryl prefix f.\n");
    fprintf(stderr, "    -onlyn f n              For mask, mers with count >= n are masked.\n");
    fprintf(stderr, "                            For only, mers with count <= n are used.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input Options:\n");
    fprintf(stderr, "    -queries c.fasta        Query sequences.\n");
    fprintf(stderr, "    -genomic g.fasta        Database sequences.\n");
    fprintf(stderr, "    -positions p.positionDB Build and save / use positionDB.  Assumes you aren't using -use.\n");
    fprintf(stderr, "    -buildonly              Only do the build and save.\n");
    fprintf(stderr, "    -use [...]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Process Options:\n");
    fprintf(stderr, "    -numthreads n           Use n search threads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output Options:\n");
    fprintf(stderr, "    -verbose                Entertain the user with useless statistics.\n");
    fprintf(stderr, "    -output f               Write output to file f.\n");
    fprintf(stderr, "    -{no}aligns             Enable/Disable full alignments.  Enabled by default.\n");
    fprintf(stderr, "    -log f                  Write some debugging/logging information to file f.  This\n");
    fprintf(stderr, "                            is mostly for developers, and does NOT provide useful\n");
    fprintf(stderr, "                            information unless you know the guts of snapper.\n");
    fprintf(stderr, "    -stats f                Write resource usage statistics to f.\n");
    exit(1);
  }
}
