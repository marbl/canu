#ifndef SIM4_PARAMETERS_H
#define SIM4_PARAMETERS_H

#include "mspManager.H"

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
  void  setDontForceCanonicalSplicing(bool x=true) {
    _dontForceCanonicalSplicing = x;
  };

  void  setForceStrandPrediction(bool x=true) {
    _forceStrandPrediction = x;
  };


  void  setInterspecies(bool x=true) {
    _interspecies = x;

    if (_interspecies) {
      _percentError =  0.25;
      _match        =  1;
      _mismatch     = -3;
    } else {
      _percentError =  0.20;
      _match        =  1;
      _mismatch     = -5;
    }
  }


private:
  bool          _findAllExons;
  double        _minCoverage;
  int           _minCoverageLength;
  int           _minPercentExonIdentity;

  bool          _includeDefLine;
  bool          _printAlignments;

  int           _alwaysReport;

  bool          _ignorePolyTails;
  double        _polyTailPercent;

  int           _mspThresh1;
  int           _mspThresh2;

  //  For aborting expensive polishes
  //
  u32bit        _mspLimitAbsolute;  //  Number of MSPs allowed per hit
  double        _mspLimitPercent;   //  Same, as percentage of length of cDNA

  int           _relinkWeight;

  bool          _dontForceCanonicalSplicing;
  bool          _forceStrandPrediction;

  // Interspecies comparison options.
  //    _percentError is the former #defined P
  //    _match is the former MATCH
  //    _misMatch is the former MISMATCH
  //
  bool          _interspecies;
  double        _percentError;
  int           _match;
  int           _mismatch;

  friend class Sim4;
};


#endif  //  SIM4_PARAMETERS_H
