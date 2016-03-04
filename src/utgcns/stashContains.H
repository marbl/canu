
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-APR-09 to 2015-MAY-09
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-06
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef STASH_CONTAINS_H
#define STASH_CONTAINS_H

#include "AS_global.H"
#include "gkStore.H"
#include "tgStore.H"


#include <map>
#include <algorithm>


class readLength {
public:
  uint32    idx;
  int32     len;

  bool operator<(const readLength &that) const {
    return(len < that.len);
  };
  bool operator>(const readLength &that) const {
    return(len > that.len);
  };
};


class savedChildren {
public:
  savedChildren(tgTig *tig) {
    childrenLen = tig->_childrenLen;
    childrenMax = tig->_childrenMax;
    children    = tig->_children;

    numContains = 0;
    covContain  = 0.0;
    percContain = 0.0;

    numContainsRemoved = 0;
    covContainsRemoved = 0.0;

    numContainsSaved   = 0;
    covContainsSaved   = 0.0;

    numDovetails = 0;
    covDovetail  = 0.0;
    percDovetail = 0.0;
  };

  void   reportDetected(FILE *out, uint32 id) {
    fprintf(out, "  unitig %d detected "F_S32" contains (%.2fx, %.2f%%) "F_S32" dovetail (%.2fx, %.2f%%)\n",
            id,
            numContains,  covContain,  percContain,
            numDovetails, covDovetail, percDovetail);
  };

  void   reportRemoved(FILE *out, uint32 id) {
      fprintf(out, "    unitig %d removing "F_S32" (%.2fx) contained reads; processing only "F_S32" contained (%.2fx) and "F_S32" dovetail (%.2fx) reads\n",
              id,
              numContainsRemoved, covContainsRemoved,
              numContainsSaved,   covContainsSaved,
              numDovetails,       covDovetail);
  };


  //  The saved children proper.

  uint32      childrenLen;
  uint32      childrenMax;
  tgPosition *children;

  //  Stats on the filtering, for logging by the caller

  uint32      numContains;
  double      covContain;
  double      percContain;

  uint32      numContainsRemoved;
  double      covContainsRemoved;

  uint32      numContainsSaved;
  double      covContainsSaved;

  uint32      numDovetails;
  double      covDovetail;
  double      percDovetail;
};


savedChildren *
stashContains(tgTig  *tig,
              double  maxCov,
              bool    beVerbose = false);


void
unstashContains(tgTig         *tig,
                savedChildren *saved);


#endif // STASH_CONTAINS_H
