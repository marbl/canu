
#ifndef MULTIALIGNUNITIG_H
#define MULTIALIGNUNITIG_H

static const char *rcsid_MULTIALIGNUNITIG_H = "$Id:  $";

#include "AS_global.H"

#include "tgStore.H"
#include "abAbacus.H"

class ALNoverlap;
class NDalign;

class unitigConsensus {
public:
  unitigConsensus(gkStore  *gkpStore_,
                  double    errorRate_,
                  double    errorRateMax_,
                  uint32    minOverlap_);
  ~unitigConsensus();

  bool   generate(tgTig     *tig,
                  uint32    *failed);

  int32  initialize(uint32 *failed); 

  void   setErrorRate(double errorRate_)   { errorRate  = errorRate_;  };
  void   setMinOverlap(uint32 minOverlap_) { minOverlap = minOverlap_; };

  bool   showProgress(void)         { return(tig->_utgcns_verboseLevel >= 1); };  //  One -V, displays which reads are processing
  bool   showAlgorithm(void)        { return(tig->_utgcns_verboseLevel >= 2); };  //  
  bool   showPlacementBefore(void)  { return(tig->_utgcns_verboseLevel >= 3); };
  bool   showPlacement(void)        { return(tig->_utgcns_verboseLevel >= 3); };
  bool   showAlignments(void)       { return(tig->_utgcns_verboseLevel >= 4); };
  bool   showMultiAlignments(void)  { return(tig->_utgcns_verboseLevel >= 4); };

  void   reportStartingWork(void);
  void   reportFailure(uint32 *failed);
  void   reportSuccess(uint32 *failed);

  int32  moreFragments(void)  { tiid++;  return (tiid < numfrags); };

  int32  computePositionFromAnchor(void);
  int32  computePositionFromLayout(void);
  int32  computePositionFromAlignment(void);

  void   rebuild(bool recomputeFullConsensus, bool display);

  bool   rejectAlignment(bool allowBhang, bool allowAhang, ALNoverlap *O);

  int32  alignFragment(void);
  void   applyAlignment(void);

  void   generateConsensus(void);
  void   exportToTig(void);

private:
  gkStore        *gkpStore;

  tgTig          *tig;
  uint32          numfrags;    //  == tig->numberOfChildren()

  uint32          traceLen;
  uint32          traceMax;
  int32          *trace;

  int32           traceABgn;  //  Start of the aligned regions, discovered in alignFragment(),
  int32           traceBBgn;  //  used in applyAlignment().

  abAbacus       *abacus;
  abMultiAlignID  multialign;  //  Should really be storing the abMultiAlignID

  //  The two positions below are storing the low/high coords for the read.
  //  They do not encode the orientation in the coordinates.
  //
  tgPosition     *utgpos;      //  Original unitigger location.
  tgPosition     *cnspos;      //  Actual location in frankenstein.

  int32           tiid;        //  This frag IID
  int32           piid;        //  Anchor frag IID - if -1, not valid

  int32           frankensteinLen;
  int32           frankensteinMax;
  char           *frankenstein;
  abBeadID       *frankensteinBof;

  uint32          minOverlap;
  double          errorRate;
  double          errorRateMax;

  NDalign        *oaPartial;
  NDalign        *oaFull;
};


#endif
