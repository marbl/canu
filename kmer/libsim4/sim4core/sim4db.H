#ifndef SIM4DBPARAMETERS_H
#define SIM4DBPARAMETERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "libbri.H"
#include "buildnumber.H"

#if 0
class sim4dbTiming {
public:
  double        _startTime;
  double        _buildTime;
  double        _searchTime;
  double        _sortTime;
  double        _linkTime;
  double        _exoncoreTime;
  double        _dbReadTime;
  double        _scReadTime;
  double        _qsReadTime;
  double        _paramTime;
  double        _extendTime;

  sim4dbTiming() {
    _startTime    = 0.0;
    _buildTime    = 0.0;
    _searchTime   = 0.0;
    _sortTime     = 0.0;
    _linkTime     = 0.0;
    _exoncoreTime = 0.0;
    _dbReadTime   = 0.0;
    _scReadTime   = 0.0;
    _qsReadTime   = 0.0;
    _paramTime    = 0.0;
    _extendTime   = 0.0;
  };
  ~sim4dbTiming() {
  };
};
#endif


class sim4dbParameters {
public:
  unsigned int  _numThreads;

  char         *_cdnaFileName;
  char         *_scriptFileName;
  char         *_databaseFileName;
  char         *_outputFileName;
  char         *_statsFileName;
  char         *_touchFileName;

#if 0
  bool          _useDefaultCforK;
#endif
  bool          _printAlignments;
  bool          _beVerbose;        //  Print progress
  bool          _beExplicit;       //  Print each script line as we process

  bool          _findAllExons;
  double        _minCoverage;
  int           _minCoverageLength;
  int           _minPercentExonIdentity;
  bool          _includeDefLine;

  int           _alwaysReport;

  bool          _ignorePolyTails;
  double        _polyTailPercent;

  int           _mspThresh1;
  int           _mspThresh2;
#if ABORT_EXPENSIVE
  unsigned int  _mspLimitAbsolute;  //  Number of MSPs allowed per hit
  double        _mspLimitPercent;   //  Same, as percentage of length of cDNA
#endif
  int           _relinkWeight;

  bool          _dontForceCanonicalSplicing;
  bool          _forceStrandPrediction;

  sim4dbParameters() {
    _numThreads       = 0;
    _cdnaFileName     = 0L;
    _scriptFileName   = 0L;
    _databaseFileName = 0L;
    _outputFileName   = 0L;
    _statsFileName    = 0L;
    _touchFileName    = 0L;

    _printAlignments        = false;
    _beVerbose              = false;
    _beExplicit             = false;

    _findAllExons           = false;
    _minCoverage            = 0.0;
    _minCoverageLength      = 0;
    _minPercentExonIdentity = 0;
    _includeDefLine         = true;

    _alwaysReport           = 0;

    _ignorePolyTails        = true;
    _polyTailPercent        = 0.75;

    _mspThresh1             = 0;
    _mspThresh2             = 0;
#if ABORT_EXPENSIVE
    _mspLimitAbsolute       = 0;
    _mspLimitPercent        = 0.0;
#endif

    _relinkWeight           = DEFAULT_RELINK_WEIGHT;

    _dontForceCanonicalSplicing = false;
    _forceStrandPrediction      = false;
  };
  ~sim4dbParameters() {
  };

  void usage(char *name) {
    fprintf(stderr, "usage: %s [-v] [-V]\n", name);
    //fprintf(stderr "        [-t numThreads]\n");
    fprintf(stderr, "       [-cdna file] [-genomic file] [-script file] [-output file] [-stats file] [-touch file]\n");
    fprintf(stderr, "       [-mincoverage percent] [-minidentity percent, 0.0 <-> 1.0] [-minlength length-in-bp]\n");
    fprintf(stderr, "       [-nodefline] [-alignments]\n");
    fprintf(stderr, "       [-polytails] [-cut x]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -v            print status while running\n");
    fprintf(stderr, "       -V            print script lines as they are processed\n");
    //fprintf(stderr, "       -t            use numThreads threads (default = adaptive)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -cdna         use these cDNA sequences\n");
    fprintf(stderr, "       -genomic      use these genomic sequences\n");
    fprintf(stderr, "       -script       use this script file\n");
    fprintf(stderr, "       -output       write output to this file\n");
    fprintf(stderr, "       -stats        write execution statistics to this file\n");
    fprintf(stderr, "       -touch        create this file when the program finishes execution\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -mincoverage  iteratively find all exon models with the specified\n");
    fprintf(stderr, "                     minimum PERCENT COVERAGE\n");
    fprintf(stderr, "       -minidentity  iteratively find all exon models with the specified\n");
    fprintf(stderr, "                     minimum PERCENT EXON IDENTITY\n");
    fprintf(stderr, "       -minlength    iteratively find all exon models with the specified\n");
    fprintf(stderr, "                     minimum ABSOLUTE COVERAGE (number of bp matched)\n");
    fprintf(stderr, "       -alwaysreport always report <number> exon models, even if they\n");
    fprintf(stderr, "                     are below the quality thresholds\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "         If no mincoverage or minidentity or minlength is given, only\n");
    fprintf(stderr, "         the best exon model is returned.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "         You will probably want to specify ALL THREE of mincoverage,\n");
    fprintf(stderr, "         minidentity and minlength!  Don't assume the default values\n");
    fprintf(stderr, "         are what you want!\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "         You will DEFINITELY want to specify at least one of mincoverage,\n");
    fprintf(stderr, "         minidentity and minlength with alwaysreport!  If you don't, mincoverage\n");
    fprintf(stderr, "         will be set to 90 and minidentity to 95 -- to reduce the number of\n");
    fprintf(stderr, "         spurious matches when a good match is found.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -nodeflines   don't include the defline in the output\n");
    fprintf(stderr, "       -alignments   print alignments\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -polytails    DON'T mask poly-A and poly-T tails.\n");
    fprintf(stderr, "       -cut          Trim marginal exons if A/T %% > x (poly-AT tails)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -noncanonical Don't force canonical splice sites\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -forcestrand  Force the strand prediction to always be\n");
    fprintf(stderr, "                     'forward' or 'reverse'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  The following are for use only by immortals.\n");
    fprintf(stderr, "       -H            set the relink weight factor\n");
    fprintf(stderr, "       -K            set the first MSP threshold\n");
    fprintf(stderr, "       -C            set the second MSP threshold\n");
    fprintf(stderr, "       -Ma           set the limit of the number of MSPs allowed\n");
    fprintf(stderr, "       -Mp           same, as percentage of bases in cDNA\n");
    fprintf(stderr, "       -W            ?\n");
    fprintf(stderr, "       -X            ?\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "This is sim4dbseq, built on %s %s (build number %d)\n", __DATE__, __TIME__, buildNumber());
    fprintf(stderr, "\n");
#ifdef INTERSPECIES
#warn INTERSPECIES MODE ENABLED IN sim4db.H
    fprintf(stderr, "!!! NOTE !!!  This is a modified version tuned for interspecies comparisons\n");
    fprintf(stderr, "Direct questions to liliana.florea@celera.com, brian.walenz@celera.com\n");
    fprintf(stderr, "\n");
#endif
  };

  void read(int argc, char **argv) {
    int arg = 1;

    while (arg < argc) {
      switch (argv[arg][1]) {
      case 'a':
        if (strncmp(argv[arg], "-ali", 4) == 0) {
          _printAlignments = true;
        } else if (strncmp(argv[arg], "-alw", 4) == 0) {
          _findAllExons = true;
          arg++;
          _alwaysReport = atoi(argv[arg]);
        } else {
          //fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
        }
        break;
      case 'c':
        if (strncmp(argv[arg], "-cd", 3) == 0) {
          arg++;
          _cdnaFileName = argv[arg];
        } else if (strncmp(argv[arg], "-cu", 3) == 0) {
          arg++;
          _polyTailPercent = atof(argv[arg]);
          if (_polyTailPercent < 0.0) {
            fprintf(stderr, "WARNING:  -cut adjusted to 0.0 (you gave %f)!\n", _polyTailPercent);
            _polyTailPercent = 0.0;
          }
          if (_polyTailPercent > 1.0) {
            fprintf(stderr, "WARNING:  -cut adjusted to 1.0 (you gave %f)!\n", _polyTailPercent);
            _polyTailPercent = 1.0;
          }
        } else {
          //fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
        }
        break;
      case 'g':
        arg++;
        _databaseFileName = argv[arg];
        break;
      case 'm':
        if (strncmp(argv[arg], "-minc", 5) == 0) {
          _findAllExons = true;
          arg++;
          _minCoverage  = atoi(argv[arg]) / 100.0;
        } else if (strncmp(argv[arg], "-mini", 5) == 0) {
          _findAllExons = true;
          arg++;
          _minPercentExonIdentity = atoi(argv[arg]);
        } else if (strncmp(argv[arg], "-minl", 5) == 0) {
          _findAllExons = true;
          arg++;
          _minCoverageLength = atoi(argv[arg]);
        } else {
          //fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
        }
        break;
      case 'n':
        if (strncmp(argv[arg], "-nod", 4) == 0) {
          _includeDefLine = false;
        } else if (strncmp(argv[arg], "-non", 4) == 0) {
          _dontForceCanonicalSplicing = true;
        }
        break;
      case 'f':
        if (strncmp(argv[arg], "-f", 2) == 0)
          _forceStrandPrediction = true;
        break;
      case 'o':
        arg++;
        _outputFileName = argv[arg];
        break;
      case 'p':
        _ignorePolyTails = false;
        break;
      case 's':
        if (strncmp(argv[arg], "-sc", 3) == 0) {
          arg++;
          _scriptFileName = argv[arg];
        } else if (strncmp(argv[arg], "-st", 3) == 0) {
          arg++;
          _statsFileName = argv[arg];
        } else {
          //fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
        }
        break;
      case 't':
        if (strncmp(argv[arg], "-th", 3) == 0) {
          arg++;
          _numThreads = (unsigned int)atoi(argv[arg]);
        } else if (strncmp(argv[arg], "-to", 3) == 0) {
          arg++;
          _touchFileName = argv[arg];
        } else {
          //fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
        }
        break;
      case 'v':
        _beVerbose = true;
        break;
      case 'V':
        _beExplicit = true;
        break;
      case 'H':
        arg++;
        _relinkWeight = atoi(argv[arg]);
        break;
      case 'K':
        arg++;
        _mspThresh1 = atoi(argv[arg]);
        break;
      case 'C':
        arg++;
        _mspThresh2 = atoi(argv[arg]);
        break;
#if ABORT_EXPENSIVE
        //
        //  Both must be supplied!
        //
      case 'M':
        switch (argv[arg][2]) {
          case 'a':
            arg++;
            _mspLimitAbsolute = atoi(argv[arg]);
            break;
          case 'p':
            arg++;
            _mspLimitPercent = atof(argv[arg]);
            break;
          default:
            break;
        }
        break;
#endif
      default:
        //fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
        break;
      }
      arg++;
    }

    //fprintf(stderr, "mspLimits: a=%d p=%f\n", _mspLimitAbsolute, _mspLimitPercent);

    if ((_findAllExons) &&
        (_minCoverage == 0.0) &&
        (_minPercentExonIdentity == 0) &&
        (_minCoverageLength == 0)) {
      _minCoverage            = 0.90;
      _minCoverageLength      = 0;
      _minPercentExonIdentity = 90;
    }
  };

  void display(void) {
  };
};

extern sim4dbParameters dbParams;

void sim4db(char           **scriptLines,
            unsigned int     scriptLinesNum,
            SequenceManager *DBs,
            SequenceManager *ESTs,
            int              fOutput);

void sim4dball(SequenceManager *DBs,
               SequenceManager *ESTs,
               int              fOutput);


#endif  //  SIM4DBPARAMETERS_H
