#ifndef SIM4_PARAMETERS_H
#define SIM4_PARAMETERS_H

#include "mspManager.H"

class sim4parameters {
public:
  sim4parameters();
  ~sim4parameters();

  void  setFindAllExons(bool x) {
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
  void  setIncludeDefLine(bool x) {
    _includeDefLine = x;
  };
  void  setPrintAlignments(bool x) {
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

#if ABORT_EXPENSIVE
  void  setMSPLimitAbsolute(u32bit t) {
    _mspLimitAbsolute = t;
  };
  void  setMSPLimitPercent(double p) {
    _mspLimitPercent = p;
  };
#endif

  void  setRelinkWeight(int x) {
    _relinkWeight = x;
  };
  void  setDontForceCanonicalSplicing(bool x) {
    _dontForceCanonicalSplicing = x;
  };

  void  setForceStrandPrediction(bool x) {
    _forceStrandPrediction = x;
  };


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
#if ABORT_EXPENSIVE
  u32bit        _mspLimitAbsolute;  //  Number of MSPs allowed per hit
  double        _mspLimitPercent;   //  Same, as percentage of length of cDNA
#endif
  int           _relinkWeight;

  bool          _dontForceCanonicalSplicing;
  bool          _forceStrandPrediction;

  friend class Sim4;
};


#endif  //  SIM4_PARAMETERS_H
