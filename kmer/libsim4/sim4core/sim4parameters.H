#ifndef SIM4_PARAMETERS_H
#define SIM4_PARAMETERS_H

#include "mspManager.H"
#include "sim4defines.H"
#include "ckalloc.H"

class sim4parameters {
public:
  sim4parameters();
  ~sim4parameters();

  void  setFindAllExons(bool x=true) {
    _findAllExons = x;
  };
  void  setMinCoverage(double x) {
    _minCoverage = x;
  };
  void  setMinCoverageLength(int l) {
    _minCoverageLength = l;
  };
  void  setMinPercentExonIdentity(int l) {
    _minPercentExonIdentity = l;
  };
  void  setIncludeDefLine(bool x=true) {
    _includeDefLine = x;
  };
  void  setPrintAlignments(bool x=true) {
    _printAlignments = x;
  };

  void  setAlwaysReport(int n) {
    _alwaysReport = n;
  };

  void  setIgnorePolyTails(bool x) {
    _ignorePolyTails = x;
  };
  void  setPolyTailPercent(double x) {
    _polyTailPercent = x;
  };

  void  setMSPThreshold1(int t) {
    _mspThresh1 = t;
  };
  void  setMSPThreshold2(int t) {
    _mspThresh2 = t;
  };

  void  setMSPLimitAbsolute(u32bit t) {
    _mspLimitAbsolute = t;
  };
  void  setMSPLimitPercent(double p) {
    _mspLimitPercent = p;
  };

  void  setRelinkWeight(int x) {
    _relinkWeight = x;
  };

  void  setWordSize(int w) {
    _wordSize = w;
  };
  void  setWordSizeInt(int w) {
    _wordSizeInt = w;
  };
  void  setWordSizeExt(int w) {
    _wordSizeExt = w;
  };

  void  setSpacedSeed(char *z) {
    _spacedSeed = strsave(z);
    _isSetSpacedSeed = true;
  };

  void  setSpliceModel(int j) { 
    _spliceModel = j;
    _isSetSpliceModel = true;
  };

  void  setDontForceCanonicalSplicing(bool x=true) {
    _dontForceCanonicalSplicing = x;
  };

  void  setForceStrandPrediction(bool x=true) {
    _forceStrandPrediction = x;
  };

  void  setSlideIntrons(bool x=true) {
    _slideIntrons = x;
  };

  void  setInterspecies(bool x=true) {
    _interspecies = x;

    if (_interspecies) {
      _percentError =  0.45;
      _match        =  1;
      _imismatch    = -1;
      _vmismatch    = -3;

      if (_isSetSpliceModel == false)
         _spliceModel = SPLICE_GENESPLICER;

      if (_isSetSpacedSeed == false)
         _spacedSeed = strsave(DEFAULT_SPACED_SEED);
      _spacedSeedInt = strsave(DEFAULT_SPACED_SEED_INT);
      _spacedSeedExt = strsave(DEFAULT_SPACED_SEED_EXT);

    } else {
      _percentError =  0.20;
      _match        =  1;
      _imismatch    = -5;
      _vmismatch    = -5;
    }
  }

  int   setSpliceMutex(void) {
    return pthread_mutex_init(&_splice_mutex,NULL);
  }


private:
  double        _minCoverage;
  int           _minCoverageLength;
  int           _minPercentExonIdentity;

  int           _alwaysReport;

  bool          _findAllExons;
  bool          _includeDefLine;
  bool          _printAlignments;

  bool          _dontForceCanonicalSplicing;
  bool          _forceStrandPrediction;

#ifdef __APPLE__
  int           _pad1;
  bool          _pad2;
  bool          _pad3;
#endif

  bool          _ignorePolyTails;
  double        _polyTailPercent;

  int           _mspThresh1;
  int           _mspThresh2;

  //  For aborting expensive polishes
  //
  double        _mspLimitPercent;   //  Same as below, as percentage of length of cDNA
  u32bit        _mspLimitAbsolute;  //  Number of MSPs allowed per hit

  int           _relinkWeight;

  int           _wordSize;
  int           _wordSizeInt;
  int           _wordSizeExt;

  char         *_spacedSeed;
  char         *_spacedSeedInt;
  char         *_spacedSeedExt;
  bool          _isSetSpacedSeed;

  int           _spliceModel;
  bool          _isSetSpliceModel;

  pthread_mutex_t
                _splice_mutex;

  // Interspecies comparison options.
  //    _percentError is the former #defined P
  //    _match is the former MATCH
  //    _misMatch is the former MISMATCH
  //

  bool          _slideIntrons;

  // Interspecies comparison options.
  //    _percentError is the former #defined P
  //    _match is the former MATCH
  //    _misMatch is the former MISMATCH
  //
  bool          _interspecies;
#ifdef __APPLE__
  bool          _pad4[3];
  u32bit        _pad5;
#endif

  double        _percentError;
  int           _match;
  int           _imismatch;
  int           _vmismatch;
  int           _imatchdiff;
  int           _vmatchdiff;
  
  friend class Sim4;
};


#endif  //  SIM4_PARAMETERS_H
