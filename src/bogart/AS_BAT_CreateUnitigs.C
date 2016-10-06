
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
 *    Brian P. Walenz beginning on 2016-OCT-03
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_AssemblyGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_TigVector.H"





//  Break on at a specific position.  In converting to unitigs, the position
//  is the end of a read with an intersection.
//
//  _bgn == true  -> reads that begin at after position are in the region
//  _end == false -> reads that end before position are in the region

class breakPointEnd {
public:
  breakPointEnd(uint32 tigID, uint32 pos, bool bgn) {
    _tigID    = tigID;
    _pos      = pos;
    _bgn      = bgn;
  };
  ~breakPointEnd() {
  };

  bool     operator<(breakPointEnd const &that) const {
    uint64  a =      _tigID;  a <<= 32;  a |=      _pos;  a <<= 1;  a |=      _bgn;  //  Because _tigID is 32-bit
    uint64  b = that._tigID;  b <<= 32;  b |= that._pos;  b <<= 1;  b |= that._bgn;

    return(a < b);
  };

  uint32  _tigID;
  uint32  _pos;
  bool    _bgn;
};







void
copyTig(TigVector    &tigs,
        Unitig       *oldtig) {
  Unitig  *newtig = tigs.newUnitig(true);

  newtig->_isUnassembled = oldtig->_isUnassembled;
  newtig->_isBubble      = oldtig->_isBubble;
  newtig->_isRepeat      = oldtig->_isRepeat;
  newtig->_isCircular    = oldtig->_isCircular;

  for (uint32 fi=0; fi<oldtig->ufpath.size(); fi++)
    newtig->addRead(oldtig->ufpath[fi], 0, false);
}





//  Split a tig based on read ends.

uint32
splitTig(TigVector                &tigs,
         Unitig                   *tig,
         vector<breakPointEnd>    &BP,
         Unitig                  **newTigs,
         int32                    *lowCoord,
         uint32                   *nMoved,
         bool                      doMove) {

  if (doMove == true) {
    memset(newTigs,  0, sizeof(Unitig *) * BP.size());
    memset(lowCoord, 0, sizeof(int32)    * BP.size());
  } else {
    memset(nMoved,   0, sizeof(uint32)   * BP.size());
  }

  if (doMove)
    for (uint32 tt=0; tt < BP.size() - 1; tt++)
      writeLog("splitTig()-- piece %2u from %8u %c to %8u %c\n",
               tt,
               BP[tt  ]._pos, BP[tt  ]._bgn ? 't' : 'f',
               BP[tt+1]._pos, BP[tt+1]._bgn ? 't' : 'f');

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode     &read   = tig->ufpath[fi];
    uint32      lo     = read.position.min();
    uint32      hi     = read.position.max();

    //  Find the intervals the end points of the read fall into.  Suppose we're trying to place
    //  the long read.  It begins in piece 1 and ends in piece 6.
    //
    //
    //   [----1---][----3----]---4---[--5---]------6-----]   Piece and boundary condition
    //   ------
    //      --------------------------------------
    //        -----
    //             ------
    //                  ------
    //                               ----
    //                                  -----
    //                                          ----------
    //
    //  The long read can not go in piece 1, as it would span the end boundary.  Piece 2 is
    //  of size zero between pieces 1 and 3, and we can place the read there.  Or, we can place
    //  it in piece 6 (we prefer piece 6).

    uint32 bgnBP = UINT32_MAX;
    uint32 endBP = UINT32_MAX;
    uint32 finBP = UINT32_MAX;

    //  Find the pieces the end points are in.

    for (uint32 tt=0; tt < BP.size()-1; tt++) {
      uint32  p = BP[tt  ]._pos;   bool  pb = BP[tt  ]._bgn;
      uint32  n = BP[tt+1]._pos;   bool  nb = BP[tt+1]._bgn;

      if ((p <= lo) && (lo < n))      //  If bgn == true  -- p == lo is in this region
        bgnBP = tt;

      if ((p < hi) && (hi <= n))      //  If bgn == false -- hi == n is in this region
        endBP = tt;
    }

    //  If both pieces are the same, we're done.

    if (bgnBP == endBP) {
      finBP = bgnBP;
    }

    //  If the next BP is a bgn boundary, we can still place the read in this piece.  It'll extend
    //  off the end, but we don't care.

    else if (BP[bgnBP+1]._bgn == true) {
      finBP = bgnBP;
    }

    //  If not, the next boundary is an end point, and we cannot place the read in this piece.
    //  If the endBP piece doesn't have restrictions on the begin, we can place the read there.

    else if (BP[endBP]._bgn == false) {
      finBP = endBP;
    }

    //  Well, shucks.  No place to put the read.  Search for an unbounded region between bgnBP and
    //  endBP.  There must be one, because bgnBP ends with a bgn=false boundary, and endBP begins
    //  with a bgn=true boundary.  If there are no intermediate boundaries, we can place the read in
    //  the middle.  If there are intermediate boundaries, we'll still have some piece that is
    //  unbounded.

    else {
      for (finBP=bgnBP+1; finBP < endBP; finBP++) {
        if ((BP[finBP  ]._bgn == false) &&
            (BP[finBP+1]._bgn == true))
          break;
      }

      if (finBP == endBP)
        writeLog("splitTig()-- failed to place read %u %u-%u in a region.  found bgn=%u and end=%u\n",
                 read.ident, read.position.bgn, read.position.end, bgnBP, endBP);
      assert(finBP != endBP);
    }

    //  Make a new tig, if needed

    if ((doMove == true) && (newTigs[finBP] == NULL)) {
      writeLog("splitTig()-- Make new tig at read %u\n", read.ident);
      lowCoord[finBP] = read.position.min();
      newTigs[finBP]  = tigs.newUnitig(true);
    }

    //  Now move the read, or account for moving it.

    if (doMove) {
      writeLog("splitTig()-- Move read %8u %8u-%-8u to piece %2u tig %6u\n",
               read.ident, read.position.bgn, read.position.end, finBP, newTigs[finBP]->id());
      newTigs[finBP]->addRead(read, -lowCoord[finBP], false);
    }
    else {
      //writeLog("splitTig()-- Move read %u %u-%u to piece %u (pos=%u)\n", read.ident, read.position.bgn, read.position.end, bp, BP[finBP]._pos);
      nMoved[finBP]++;
    }
  }

  //  Return the number of tigs created.

  uint32  nTigsCreated = 0;

  for (uint32 ii=0; ii<BP.size(); ii++)
    if (nMoved[ii] > 0)
      nTigsCreated++;

  return(nTigsCreated);
}






void
createUnitigs(AssemblyGraph  *AG,
              TigVector      &contigs,
              TigVector      &unitigs) {

  vector<breakPointEnd>   breaks;

  //  Check the reads at the end of every tig for intersections to other tigs.

  writeLog("\n");
  writeLog("Finding breakpoints.\n");
  writeLog("\n");

  for (uint32 ti=0; ti<contigs.size(); ti++) {
    Unitig    *tig = contigs[ti];

    if (tig == NULL)
      continue;
 
    if (tig->_isUnassembled == true)    //  Edge is FROM an unassembled thing, ignore it.
      continue;

    uint32  fi = tig->firstRead()->ident;
    uint32  li = tig->lastRead() ->ident;

    //  Give this tig a bogus breakpoint, just to get it in the list.  If there are no break points,
    //  it won't be split.  These also serve as sentinels during splitting.

    breaks.push_back(breakPointEnd(ti, 0,                true));    //  Add one at the start of the tig
    breaks.push_back(breakPointEnd(ti, tig->getLength(), false));   //  And one at the end

    //  Check the first read.

    if (AG->getForward(fi).size() + AG->getForward(li).size() > 0)
      writeLog("createUnitigs()-- tig %u len %u first read %u with %lu edges - last read %u with %lu edges\n",
               ti, tig->getLength(),
               fi, li,
               AG->getForward(fi).size(),
               AG->getForward(li).size());

    for (uint32 pp=0; pp<AG->getForward(fi).size(); pp++) {
      BestPlacement  &pf = AG->getForward(fi)[pp];

      //writeLog("createUnitigs()-- read %u pp %u - edge to position %d-%d isContig %u isRepeat %u FIRST\n",
      //        fi, pp, pf.placedBgn, pf.placedEnd, pf.isContig, pf.isRepeat);

      if ((pf.isRepeat == true) ||   //  isRepeat means olap is on the same end as the tig, but to a diff tig.
          (pf.isContig == true))     //  isContig means olap is on the same end as the tig, and to the same tig.
        continue;

      if (pf.bestC.b_iid > 0)        //  And if C, the read is fully contained in the other read, and can't break the tig.
        continue;

      BAToverlap     best  = (pf.best5.b_iid > 0) ? pf.best5 : pf.best3;
      Unitig        *btig  = contigs[ contigs.inUnitig(best.b_iid) ];
      ufNode        *read  = &btig->ufpath[ contigs.ufpathIdx(best.b_iid) ];
      bool           isL   = (read->position.isForward()) ? best.BEndIs5prime() : best.BEndIs3prime();
      uint32         coord = (isL == true) ? read->position.min() : read->position.max();

      if (btig->_isUnassembled == true)    //  Edge is TO an unassembled thing, ignore it.
        continue;

      breaks.push_back(breakPointEnd(pf.tigID, coord, isL));

      writeLog("splitThinEdge()-- read %6u splits tig %5u at coordinate %8u via intersection with read %6u isL %u\n",
              fi, pf.tigID, coord, best.b_iid, isL);
    }

    //  Check the last read

    for (uint32 pp=0; pp<AG->getForward(li).size(); pp++) {
      BestPlacement  &pf = AG->getForward(li)[pp];

      //writeLog("createUnitigs()-- read %u pp %u - edge to position %d-%d isContig %u isRepeat %u LAST\n",
      //        li, pp, pf.placedBgn, pf.placedEnd, pf.isContig, pf.isRepeat);

      if ((pf.isRepeat == true) ||
          (pf.isContig == true))
        continue;

      if (pf.bestC.b_iid > 0)
        continue;

      BAToverlap     best  = (pf.best5.b_iid > 0) ? pf.best5 : pf.best3;
      Unitig        *btig  = contigs[ contigs.inUnitig(best.b_iid) ];
      ufNode        *read  = &btig->ufpath[ contigs.ufpathIdx(best.b_iid) ];
      bool           isL   = (read->position.isForward()) ? best.BEndIs5prime() : best.BEndIs3prime();
      uint32         coord = (isL == true) ? read->position.min() : read->position.max();

      if (btig->_isUnassembled == true)    //  Edge is TO an unassembled thing, ignore it.
        continue;

      breaks.push_back(breakPointEnd(pf.tigID, coord, isL));

      writeLog("splitThinEdge()-- read %6u splits tig %5u at coordinate %8u via intersection with read %6u isL %u\n",
              li, pf.tigID, coord, best.b_iid, isL);
    }
  }


  //  The splitTigs function operates only on a single tig.  Sort the break points
  //  by tig id to find all the break points for each tig.

  sort(breaks.begin(), breaks.end());

  writeLog("\n");
  writeLog("createUnitigs()-- Found %u breakpoints.\n", breaks.size());
  writeLog("\n");

  //  Allocate space for breaking tigs.  These are _vastly_ too big, but guaranteed.

  vector<breakPointEnd>  BP;

  Unitig **newTigs   = new Unitig * [breaks.size() + 2];  //  Plus two, because we add an extra
  int32   *lowCoord  = new int32    [breaks.size() + 2];  //  break at the start and end
  uint32  *nMoved    = new uint32   [breaks.size() + 2];  //  of each set.

  //  Walk through the breaks, making a new vector of breaks for each tig.

  uint32  ss = 0;
  uint32  ee = 0;
  
  while (ss < breaks.size()) {
    Unitig  *tig = contigs[breaks[ss]._tigID];

    //  Find the last break point for this tig.  (Technically, the one after the last, but...)

    while ((ee < breaks.size()) && (breaks[ss]._tigID == breaks[ee]._tigID))
      ee++;

    //  Make a new vector for those break points.

    BP.clear();

    writeLog("createUnitigs()-- Process BPs from ss=%u to ee=%u\n", ss, ee);

    for (uint32 bb=ss; bb<ee; bb++) {
      //writeLog("createUnitigs()--   BP[%3u] pos %u bgn %c RAW\n", bb, breaks[bb]._pos, (breaks[bb]._bgn) ? 't' : 'f');
       if ((BP.size() == 0) ||
           (BP.back()._pos != breaks[bb]._pos) ||
          (BP.back()._bgn != breaks[bb]._bgn))
        BP.push_back(breaks[bb]);
    }

    writeLog("createUnitigs()-- tig %u found %u breakpoints (%u-%u)\n",
             tig->id(), BP.size(), ss, ee);
    for (uint32 bb=0; bb<BP.size(); bb++)
      writeLog("createUnitigs()--   BP[%2u] pos %8u %c\n", bb, BP[bb]._pos, (BP[bb]._bgn) ? 't' : 'f');

    //  Split the tig.  Copy it into the unitigs TigVector too.

    uint32  nTigs = splitTig(contigs, tig, BP, newTigs, lowCoord, nMoved, false);

    if (nTigs > 1)
      splitTig(unitigs, tig, BP, newTigs, lowCoord, nMoved, true);
    else
      copyTig(unitigs, tig);

    if (nTigs > 1)
      writeLog("createUnitigs()-- tig %u created %u new tigs.\n",
               tig->id(), nTigs);
    else
      writeLog("createUnitigs()-- tig %u copied.\n",
               tig->id(), nTigs);


    //reportTigsCreated(tig, BP, nTigs, newTigs, nMoved);

    //  Delete the old tig, if we didn't copy it to unitigs.
    //
    //if (nTigs > 1) {
    //  contigs[tig->id()] = NULL;
    //  delete tig;
    //}

    ss = ee;   //  Reset for the next iteration.
  }

  //  Cleanup.

  delete [] newTigs;
  delete [] lowCoord;
  delete [] nMoved;
}

