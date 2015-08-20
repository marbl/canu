
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_CNS/MultiAlignment_CNS.h
 *    src/AS_CNS/MultiAlignment_CNS_private.H
 *    src/AS_CNS/MultiAlignment_CNS_private.h
 *    src/utgcns/libcns/MultiAlignment_CNS_private.H
 *
 *  Modifications by:
 *
 *    Gennady Denisov from 2005-MAY-23 to 2007-OCT-25
 *      are Copyright 2005-2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUL-08 to 2013-AUG-01
 *      are Copyright 2005-2009,2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2006-FEB-13 to 2008-FEB-13
 *      are Copyright 2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-JAN-28 to 2009-SEP-25
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2011-SEP-21
 *      are Copyright 2011 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz beginning on 2014-NOV-17
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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

  int32  alignFragment(bool forceAlignment=false);
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
