
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
 *    src/overlapInCore/overlapPair-readCache.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUN-16 to 2015-JUN-23
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"
#include "tgStore.H"

class overlapReadCache {
public:
  overlapReadCache(gkStore *gkpStore_, uint64 memLimit);
  ~overlapReadCache();

private:
  void         loadRead(uint32 id);
  void         loadReads(set<uint32> reads);
  void         markForLoading(set<uint32> &reads, uint32 id);

public:
  void         loadReads(ovOverlap *ovl, uint32 nOvl);
  void         loadReads(tgTig *tig);

  void         purgeReads(void);

  char        *getRead(uint32 id) {
    assert(readLen[id] > 0);
    return(readSeqFwd[id]);
  };

  uint32       getLength(uint32 id) {
    assert(readLen[id] > 0);
    return(readLen[id]);
  };

private:
  gkStore     *gkpStore;
  uint32       nReads;

  uint32      *readAge;
  uint32      *readLen;
  char       **readSeqFwd;
  //char       **readSeqRev;  //  Save it, or recompute?

  gkReadData   readdata;

  uint64       memoryLimit;
};



