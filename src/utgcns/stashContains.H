
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef STASH_CONTAINS_H
#define STASH_CONTAINS_H

#include "runtime.H"
#include "sqStore.H"
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
  savedChildren() {
    childrenLen        = 0;
    childrenMax        = 0;
    children           = NULL;

    numContains        = 0;
    covContain         = 0.0;
    percContain        = 0.0;

    numDovetails       = 0;
    covDovetail        = 0.0;
    percDovetail       = 0.0;

    numContainsSaved   = 0;
    covContainsSaved   = 0.0;

    numContainsRemoved = 0;
    covContainsRemoved = 0.0;
  };

  ~savedChildren() {
    delete [] children;
  };


  //  The saved children proper.

  uint32      childrenLen;
  uint32      childrenMax;
  tgPosition *children;

  //  Stats on the filtering, for logging by the caller

  uint32      numContains;
  double      covContain;
  double      percContain;

  uint32      numDovetails;
  double      covDovetail;
  double      percDovetail;

  uint32      numContainsSaved;
  double      covContainsSaved;

  uint32      numContainsRemoved;
  double      covContainsRemoved;
};


savedChildren *
stashContains(tgTig  *tig,
              double  maxCov,
              bool    beVerbose = false);


void
unstashContains(tgTig         *tig,
                savedChildren *saved);


#endif // STASH_CONTAINS_H
