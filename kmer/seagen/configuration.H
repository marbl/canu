#ifndef CONFIGURATION_H
#define CONFIGURATION_H

//
//  A singleton for working with the command line parameters.
//

#include "bio++.H"
#include "util++.H"

class configuration {
public:
  bool             _beVerbose;

  u32bit           _merSize;
  u32bit           _numSearchThreads;

  u32bit           _merSkip;

  bool             _doReverse;
  bool             _doForward;

  u32bit           _maxDiagonal;
  u32bit           _maxGap;
  u32bit           _qsOverlap;
  u32bit           _dsOverlap;

  u32bit           _extendWeight;
  u32bit           _extendMinimum;
  bool             _extendAlternate;

  u32bit           _maxIntronLength;

  u32bit           _smallSequenceCutoff;

  u32bit           _minLengthSingle;
  double           _minCoverageSingle;

  double           _minCoverageMultiple;
  u32bit           _minLengthMultiple;

  char            *_dbFileName;
  char            *_qsFileName;
  char            *_maskFileName;
  char            *_onlyFileName;
  char            *_outputFileName;
  char            *_queryMatchFileName;
  char            *_statsFileName;

  int              _outputFile;
  int              _matchCountsFile;

  char            *_tableTemporaryFileName;

  char            *_tableFileName;
  bool             _tableBuildOnly;

  bool             _binaryOutput;

  seqStream        _useList;

  //  Running state
  //
  seqFile          *_qsFASTA;
  existDB          *_maskDB;
  existDB          *_onlyDB;
  positionDB       *_positions;

  u32bit            _numberOfQueries;

  //  Wall clock times
  //
  double           _startTime;
  double           _initTime;
  double           _buildTime;
  double           _searchTime;

  //  sweatShop queue sizes; maximum number of sequences pre-loaded,
  //  maximum number of results waiting for output.
  //
  u32bit           _loaderQueue;
  struct timespec  _loaderSleep;
  bool             _loaderWarnings;

  struct timespec  _searchSleep;

  u32bit           _writerQueue;
  struct timespec  _writerSleep;
  bool             _writerWarnings;

  configuration();
  ~configuration();

  void  usage(char *name);
  void  read(int argc, char **argv);
  void  display(FILE *out=stdout);

  void  dumpStats(void);

  void  setTime(struct timespec *ts, double t) {
    ts->tv_sec  = (time_t)floor(t);
    ts->tv_nsec = (long)((t - ts->tv_sec) * 1e9);
  };
private:
  void  addToUse(u32bit v);
  void  parseUseLine(char *line);
};

extern configuration          config;

#endif  //  CONFIGURATION_H
