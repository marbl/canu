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

  uint32           _merSize;
  uint32           _numSearchThreads;

  uint32           _merSkip;

  bool             _doReverse;
  bool             _doForward;

  uint32           _maxDiagonal;
  uint32           _maxGap;
  uint32           _qsOverlap;
  uint32           _dsOverlap;

  uint32           _extendWeight;
  uint32           _extendMinimum;
  bool             _extendAlternate;

  uint32           _maxIntronLength;

  uint32           _smallSequenceCutoff;

  uint32           _minLengthSingle;
  double           _minCoverageSingle;

  double           _minCoverageMultiple;
  uint32           _minLengthMultiple;

  char            *_dbFileName;
  char            *_qsFileName;
  char            *_maskFileName;
  char            *_onlyFileName;
  char            *_outputFileName;
  char            *_queryMatchFileName;

  int              _outputFile;
  int              _matchCountsFile;

  char            *_tableTemporaryFileName;

  char            *_tableFileName;
  bool             _tableBuildOnly;

  bool             _binaryOutput;

  seqStream       *_dbSTREAM;

  //  Running state
  //
  seqCache         *_qsFASTA;
  existDB          *_maskDB;
  existDB          *_onlyDB;
  positionDB       *_positions;

  uint32            _numberOfQueries;

  //  Wall clock times
  //
  double           _startTime;
  double           _initTime;
  double           _buildTime;
  double           _searchTime;

  //  sweatShop queue sizes; maximum number of sequences pre-loaded,
  //  maximum number of results waiting for output.
  //
  uint32           _loaderQueue;
  struct timespec  _loaderSleep;
  bool             _loaderWarnings;

  struct timespec  _searchSleep;

  uint32           _writerQueue;
  struct timespec  _writerSleep;
  bool             _writerWarnings;

  configuration();
  ~configuration();

  void  usage(char *name);
  void  read(int argc, char **argv);

  void  setTime(struct timespec *ts, double t) {
    ts->tv_sec  = (time_t)floor(t);
    ts->tv_nsec = (long)((t - ts->tv_sec) * 1e9);
  };
private:
  void  addToUse(uint32 v);
  void  parseUseLine(char *line);
};

extern configuration          config;

#endif  //  CONFIGURATION_H
