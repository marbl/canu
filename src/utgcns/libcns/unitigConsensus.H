
#ifndef MULTIALIGNUNITIG_H
#define MULTIALIGNUNITIG_H

static const char *rcsid_MULTIALIGNUNITIG_H = "$Id:  $";

#include "AS_global.H"

//#include "gkStore.H"
#include "tgStore.H"
#include "aligners.H"
#include "abAbacus.H"

//  uint32  printWidth;
//  uint32  printSpacing;

//extern uint32 MULTIALIGN_PRINT_WIDTH;       //  In MultiAlignPrint.c
//extern uint32 MULTIALIGN_PRINT_SPACING;

//typedef enum {
//  CNS_SMOOTH = 1, // only eliminate pairwise construction artifacts
//  CNS_POLYX  = 2, // align poly-X regions
//  CNS_INDEL  = 4  // push apart mushed block indels
//}  CNS_RefineLevel;


//  Options to things in MultiAligment_CNS.c
//extern int32 allow_neg_hang;



class unitigConsensus {
public:
  unitigConsensus(tgTig *tig_) {
    tig          = tig_;
    numfrags     = tig->numberOfChildren();
    //fraglist = tig->_children;
    //fragback = NULL;
    //opp      = opp_;
    trace        = NULL;
    abacus       = NULL;
    multialign   = abMultiAlignID();
    utgpos       = NULL;
    cnspos       = NULL;
    tiid         = 0;
    piid         = -1;

    frankensteinLen = 0;
    frankensteinMax = 0;
    frankenstein    = NULL;
    frankensteinBof = NULL;
  };

  ~unitigConsensus() {
    delete [] trace;
    //delete    multialign;
    delete    abacus;
    //if (manode)
    //  DeleteMANode(manode->lid);

    delete [] utgpos;
    delete [] cnspos;
    delete [] frankenstein;
    delete [] frankensteinBof;
  };

  int32  initialize(gkStore *gkpStore, uint32 *failed); 

  bool   showProgress(void)         { return(tig->_utgcns_verboseLevel >= 1); };
  bool   showAlgorithm(void)        { return(tig->_utgcns_verboseLevel >= 2); };
  bool   showPlacementBefore(void)  { return(tig->_utgcns_verboseLevel >= 3); };
  bool   showPlacement(void)        { return(tig->_utgcns_verboseLevel >= 3); };
  bool   showAlignments(void)       { return(tig->_utgcns_verboseLevel >= 4); };

  void   reportStartingWork(void);
  void   reportFailure(uint32 *failed);
  void   reportSuccess(uint32 *failed);

  int32  moreFragments(void)  { tiid++;  return (tiid < numfrags); };

  int32  computePositionFromAnchor(void);
  int32  computePositionFromLayout(void);
  int32  computePositionFromAlignment(void);

  void   rebuild(bool recomputeFullConsensus);

  bool   rejectAlignment(bool allowBhang, bool allowAhang, ALNoverlap *O);

  int32  alignFragment(void);
  void   applyAlignment(void);

  void   generateConsensus(void);

private:
  tgTig          *tig;
  uint32          numfrags;    //  == tig->numberOfChildren()

  //IntMultiPos    *fraglist;
  //IntMultiPos    *fragback;

  uint32          traceLen;
  uint32          traceMax;
  int32          *trace;

  int32           traceBgn;

  abAbacus       *abacus;
  abMultiAlignID  multialign;  //  Should really be storing the abMultiAlignID

  //  The two positions below are storing the low/high coords for the read.
  //  They do not encode the orientation in the coordinates.
  //
  tgPosition     *utgpos;      //  Original unitigger location.
  tgPosition     *cnspos;      //  Actual location in frankenstein.

  //int32           ovl;    //  Expected overlap in bases to the frankenstein
  //int32           ahang;  //  Expected hangs to the frankenstein
  //int32           bhang;

  int32           tiid;        //  This frag IID
  int32           piid;        //  Anchor frag IID - if -1, not valid

  int32           frankensteinLen;
  int32           frankensteinMax;
  char           *frankenstein;
  abBeadID       *frankensteinBof;
};


#endif
