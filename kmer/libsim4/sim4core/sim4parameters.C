#include "sim4parameters.H"
#include "sim4defines.H"

sim4parameters::sim4parameters() {
    _findAllExons               = false;
    _minCoverage                = 0.0;
    _minCoverageLength          = 0;
    _minPercentExonIdentity     = 0;

    _includeDefLine             = true;
    _printAlignments            = false;

    _alwaysReport               = 0;

    _ignorePolyTails            = true;
    _polyTailPercent            = 0.75;

    _mspThresh1                 = 0;
    _mspThresh2                 = 0;

    _mspLimitAbsolute           = 0;
    _mspLimitPercent            = 0.0;

    _relinkWeight               = DEFAULT_RELINK_WEIGHT;

    _wordSize                   = 12;
    _wordSizeInt                = 0;
    _wordSizeExt                = 0;

    _dontForceCanonicalSplicing = false;
    _forceStrandPrediction      = false;

    _slideIntrons               = true;

    _interspecies               = false;
    _percentError               = 0.20;
    _match                      =  1;
    _mismatch                   = -5;
}

sim4parameters::~sim4parameters() {
}
