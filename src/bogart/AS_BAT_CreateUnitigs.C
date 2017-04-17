
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

#include "AS_BAT_CreateUnitigs.H"



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







Unitig *
copyTig(TigVector    &tigs,
        Unitig       *oldtig) {
  Unitig  *newtig = tigs.newUnitig(false);

  newtig->_isUnassembled = oldtig->_isUnassembled;
  newtig->_isBubble      = oldtig->_isBubble;
  newtig->_isRepeat      = oldtig->_isRepeat;
  newtig->_isCircular    = oldtig->_isCircular;

  for (uint32 fi=0; fi<oldtig->ufpath.size(); fi++)
    newtig->addRead(oldtig->ufpath[fi], 0, false);

  return(newtig);
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
      writeLog("splitTig()-- new tig %u (id=%u) at read %u %u-%u\n", tigs.size(), finBP, read.ident, read.position.min(), read.position.max());
      lowCoord[finBP] = read.position.min();
      newTigs[finBP]  = tigs.newUnitig(false);
    }

    //  Now move the read, or account for moving it.

    if (doMove) {
      //writeLog("splitTig()-- Move read %8u %8u-%-8u to piece %2u tig %6u\n",
      //         read.ident, read.position.bgn, read.position.end, finBP, newTigs[finBP]->id());
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



static
void
checkRead(AssemblyGraph *AG,
          TigVector &contigs,
          vector<breakPointEnd> &breaks,
          Unitig  *tgA, ufNode *rdA,
          bool isFirst) {

  for (uint32 pp=0; pp<AG->getForward(rdA->ident).size(); pp++) {
    BestPlacement  &pf = AG->getForward(rdA->ident)[pp];

    //  If a contained edge, we cannot split the other tig; it is correct (this read is contained in the other read).

    if (pf.bestC.b_iid > 0) {
      writeLog("createUnitigs()-- read %6u edgeTo tig %5u read %6u position %d-%d CONTAINED\n",
               rdA->ident, contigs.inUnitig(pf.bestC.b_iid), pf.bestC.b_iid, pf.placedBgn, pf.placedEnd);
      continue;
    }

    //  Decide which overlap we want to be using, based on the orientation of the read in the tig,
    //  and if it is the first or last read.
    //
    //           first == true                  first == false
    //  best5    fwd   == true  --------->      fwd   == false  <---------
    //  best3    fwd   == false <----------     fwd   == true   --------->

    BAToverlap     best = (isFirst == rdA->position.isForward()) ? pf.best5 : pf.best3;

    //  If there is no overlap on the expected end, well, that's it, nothing we can do but give up.
    //  Don't bother logging if it is the internal edge (which it shouldn't ever be, because those shouldn't
    //  be in the graph, right?)

    if (best.b_iid == 0) {
      uint32  rdC = (isFirst == rdA->position.isForward()) ? pf.best3.b_iid : pf.best5.b_iid;   //  Grab the other edge
      uint32  tgC = contigs.inUnitig(rdC);

      if (tgC != tgA->id())
        writeLog("createUnitigs()-- read %6u edgeTo tig %5u read %6u position %d-%d WRONG_END\n",
                 rdA->ident, tgC, rdC, pf.placedBgn, pf.placedEnd);
      continue;
    }

    //  Grab the tig and read we overlap to.

    Unitig   *tgB  = contigs[ contigs.inUnitig(best.b_iid) ];
    ufNode   *rdB = &tgB->ufpath[ contigs.ufpathIdx(best.b_iid) ];

    //  And find the coordinate of the break based on the orientation of the rdB and the overlap.
    //  isLow is true if the read is forward and the overlap is off of its 5' end, or
    //                if the read is reverse and the overlap is off of its 3' end

    bool      isLow = (rdB->position.isForward()) ? best.BEndIs5prime() : best.BEndIs3prime();
    uint32    coord = (isLow == true) ? rdB->position.min() : rdB->position.max();

    //  With all that done, throw out the edge if the overlap was used to form the contig itself.
    //
    //  We used to also throw out edges to validated repeats (pf.isRepeat == true), but those are
    //  indistinguishable from bubbles.

    if (pf.isContig == true) {
      writeLog("createUnitigs()-- read %6u edgeTo tig %5u at coordinate %8u via intersection with read %6u IS_%s\n",
               rdA->ident, tgB->id(), coord, rdB->ident, (pf.isContig == true) ? "CONTIG" : "REPEAT");
      continue;
    }

    //  Also chuck it out if it is to garbage.

    if (tgB->_isUnassembled == true) {
      writeLog("createUnitigs()-- read %6u edgeTo tig %5u read %6u UNASSEMBLED\n",
               rdA->ident, tgB->id(), rdB->ident);
      continue;
    }

    //  If here, we're all golden!

    writeLog("splitThinEdge()-- read %6u splits tig %5u at coordinate %8u via intersection with read %6u isLow %u\n",
             rdA->ident, pf.tigID, coord, rdB->ident, isLow);
    breaks.push_back(breakPointEnd(pf.tigID, coord, isLow));
  }
}



void
stripNonBackboneFromStart(TigVector &unitigs, Unitig *tig, bool isFirst) {
  vector<ufNode>   ufpath;
  uint32           ii = 0;

  while (RI->isBackbone(tig->ufpath[ii].ident) == false) { //  Find the first backbone read,
    unitigs.registerRead(tig->ufpath[ii].ident);
    writeLog("WARNING: unitig %u %s read %u is not backbone, removing.\n",
             tig->id(),
             isFirst ? "first" : "last ",
             tig->ufpath[ii].ident);
    ii++;
  }

  while (ii < tig->ufpath.size()) {               //  and copy to a new vector.
    ufpath.push_back(tig->ufpath[ii]);
    writeLog("SAVE     unitig %u %s read %u IS     backbone.\n",
             tig->id(),
             isFirst ? "first" : "last ",
             tig->ufpath[ii].ident);
    ii++;
  }

  tig->ufpath.swap(ufpath);                       //  assign the new vector to the tig
  tig->cleanUp();                                 //  adjust zero, find new length
  tig->reverseComplement();                       //  rebuild the idx mappings, and reverse for the next phase
}



void
createUnitigs(AssemblyGraph   *AG,
              TigVector       &contigs,
              TigVector       &unitigs,
              vector<tigLoc>  &unitigSource) {

  vector<breakPointEnd>   breaks;

  //  Check the reads at the end of every tig for intersections to other tigs.  If the read has a
  //  compatible overlap to the middle of some other tig, split the other tig into multiple unitigs.

  writeLog("\n");
  writeLog("----------------------------------------\n");
  writeLog("Finding contig-end to contig-middle intersections.\n");

  for (uint32 ti=0; ti<contigs.size(); ti++) {
    Unitig    *tig = contigs[ti];

    if (tig == NULL)
      continue;

    if (tig->_isUnassembled == true)    //  Edge is FROM an unassembled thing, ignore it.
      continue;

    //  Give this tig a pair of bogus breakpoints at the ends, just to get it in the list.  If there
    //  are no break points, it won't be split.  These also serve as sentinels during splitting.

    breaks.push_back(breakPointEnd(ti, 0,                true));    //  Add one at the start of the tig
    breaks.push_back(breakPointEnd(ti, tig->getLength(), false));   //  And one at the end

    //  Find break points in other tigs using the first and last reads.

    ufNode *fi = tig->firstRead();
    ufNode *li = tig->lastRead();

    if (AG->getForward(fi->ident).size() + AG->getForward(li->ident).size() > 0)
      writeLog("\ncreateUnitigs()-- tig %u len %u first read %u with %lu edges - last read %u with %lu edges\n",
               ti, tig->getLength(),
               fi->ident, AG->getForward(fi->ident).size(),
               li->ident, AG->getForward(li->ident).size());

    checkRead(AG, contigs, breaks, tig, fi, true);
    checkRead(AG, contigs, breaks, tig, li, false);
  }

  //  The splitTigs function operates only on a single tig.  Sort the break points
  //  by tig id to find all the break points for each tig.

  sort(breaks.begin(), breaks.end());

  writeLog("\n");
  writeLog("createUnitigs()-- Found %u breakpoints.\n", breaks.size());

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

    for (uint32 bb=ss; bb<ee; bb++)
      if ((BP.size() == 0) ||
          (BP.back()._pos != breaks[bb]._pos) ||
          (BP.back()._bgn != breaks[bb]._bgn))
        BP.push_back(breaks[bb]);

    writeLog("\n");

    if (BP.size() > 2)
      writeLog("createUnitigs()-- contig %u found %u breakpoint%s\n",
               tig->id(), BP.size()-2, (BP.size()-2 != 1) ? "s" : "");

    //  Split the tig.  Copy it into the unitigs TigVector too.

    uint32  nTigs = splitTig(contigs, tig, BP, newTigs, lowCoord, nMoved, false);

    if (nTigs > 1) {
      splitTig(unitigs, tig, BP, newTigs, lowCoord, nMoved, true);
      writeLog("createUnitigs()-- contig %u was split into %u unitigs, %u through %u.\n",  //  Can't use newTigs, because
               tig->id(), nTigs, unitigs.size() - nTigs, unitigs.size() - 1);              //  there are holes in it
    }

    else {
      newTigs[0] = copyTig(unitigs, tig);
      writeLog("createUnitigs()-- contig %u copied into unitig %u.\n", tig->id(), newTigs[0]->id());
    }

    //  Remember where these unitigs came from.

    unitigSource.resize(unitigs.size() + 1);

    for (uint32 tt=0; tt<nTigs; tt++) {
      if (newTigs[tt]) {
        uint32  id = newTigs[tt]->id();

        writeLog("createUnitigs()-- piece %3u -> tig %u from contig %u %u-%u\n",
                 tt, id, tig->id(), lowCoord[tt], lowCoord[tt] + newTigs[tt]->getLength());

        unitigSource[id].cID  = tig->id();
        unitigSource[id].cBgn = lowCoord[tt];
        unitigSource[id].cEnd = lowCoord[tt] + newTigs[tt]->getLength();
        unitigSource[id].uID  = id;
      }
    }

    //  Reset for the next iteration.

    ss = ee;
  }

  //  Remove non-backbone reads from the ends of unitigs.  These confound graph building because
  //  they can be missing overlaps.
  //
  //  If the last read in the tig is not a backbone read, we can remove it and all reads that come
  //  after it (because those reads are contained).

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig    *tig = unitigs[ti];

    if (tig == NULL)
      continue;

    //  First, check if we have any backbone reads.  If we have none, leave it as is.

    uint32  bbReads = 0;
    uint32  nbReads = 0;

    for (uint32 li=0; li<tig->ufpath.size(); li++) {
      if (RI->isBackbone(tig->ufpath[li].ident) == true)
        bbReads++;
      else
        nbReads++;
    }

    if (bbReads == 0)
      continue;

    //  Now remove non-backbone reads from the start of the tig.

    writeLog("unitig %u with %u reads, %u backbone and %u unplaced.\n",
             tig->id(), tig->ufpath.size(), bbReads, nbReads);

    stripNonBackboneFromStart(unitigs, tig, true);
    stripNonBackboneFromStart(unitigs, tig, false);
  }

  //  Cleanup.

  delete [] newTigs;
  delete [] lowCoord;
  delete [] nMoved;
}

