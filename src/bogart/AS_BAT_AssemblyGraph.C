
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

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_AssemblyGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_PlaceReadUsingOverlaps.H"

#include "intervalList.H"
#include "stddev.H"

#undef  FILTER_DENSE_BUBBLES_FROM_GRAPH
#define FILTER_DENSE_BUBBLES_THRESHOLD    3   //  Retain bubbles if they have fewer than this number of edges to other tigs

#undef  LOG_GRAPH
#undef  LOG_GRAPH_ALL



static
void
logAGbuild(uint32                     fi,
           uint32                     pp,
           vector<overlapPlacement>  &placements,
           const char                *message) {

  if (logFileFlagSet(LOG_PLACE_UNPLACED) == false)
    return;

  writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.10f %s\n",
           fi, pp,
           placements[pp].tigID,
           placements[pp].position.bgn, placements[pp].position.end,
           placements[pp].verified.bgn, placements[pp].verified.end,
           placements[pp].fCoverage,
           placements[pp].erate(),
           message);
}



static
void
logAGbuild(uint32                     fi,
           uint32                     pp,
           vector<overlapPlacement>  &placements,
           bool                       is5,
           bool                       is3,
           bool                       onLeft,
           bool                       onRight,
           const char                *message) {

  if (logFileFlagSet(LOG_PLACE_UNPLACED) == false)
    return;

  writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f Fidx %6u Lidx %6u is5 %d is3 %d onLeft %d onRight %d %s\n",
           fi, pp,
           placements[pp].tigID,
           placements[pp].position.bgn, placements[pp].position.end,
           placements[pp].verified.bgn, placements[pp].verified.end,
           placements[pp].fCoverage,
           placements[pp].erate(),
           placements[pp].tigFidx, placements[pp].tigLidx,
           is5, is3, onLeft, onRight,
           message);
}



static
void
logAGbuild(uint32                     fi,
           uint32                     pp,
           vector<overlapPlacement>  &placements,
           uint32                     idC,
           uint32                     id5,
           uint32                     id3,
           const char                *message) {

  if (logFileFlagSet(LOG_PLACE_UNPLACED) == false)
    return;

  writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f %s C=%8d 5=%8d 3=%8d\n",
           fi, pp,
           placements[pp].tigID,
           placements[pp].position.bgn, placements[pp].position.end,
           placements[pp].verified.bgn, placements[pp].verified.end,
           placements[pp].fCoverage,
           placements[pp].erate(),
           message,
           idC,
           id5,
           id3);
}



void
AssemblyGraph::buildGraph(const char   *UNUSED(prefix),
                          double        deviationRepeat,
                          double        repeatLimit,
                          TigVector    &tigs) {
  uint32  fiLimit    = RI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  //  Just some logging.  Count the number of reads we try to place.

  uint32   nToPlaceContained = 0;
  uint32   nToPlace          = 0;
  uint32   nPlacedContained  = 0;
  uint32   nPlaced           = 0;
  uint32   nFailedContained  = 0;
  uint32   nFailed           = 0;

  for (uint32 fid=1; fid<RI->numReads()+1; fid++) {
    if (tigs.inUnitig(fid) == 0)   //  Unplaced, don't care.
      continue;

    if (OG->isContained(fid))
      nToPlaceContained++;
    else
      nToPlace++;
  }

  writeStatus("\n");
  writeStatus("AssemblyGraph()-- allocating vectors for placements, %.3fMB\n",   //  vector<> is 24 bytes, pretty tiny.
              (sizeof(vector<BestPlacement>) + sizeof(vector<BestReverse>)) * (fiLimit + 1) / 1048576.0);

  _pForward = new vector<BestPlacement> [fiLimit + 1];
  _pReverse = new vector<BestReverse>   [fiLimit + 1];

  writeStatus("AssemblyGraph()-- finding edges for %u reads (%u contained), ignoring %u unplaced reads, with %d thread%s.\n",
              nToPlaceContained + nToPlace,
              nToPlaceContained,
              RI->numReads() - nToPlaceContained - nToPlace,
              numThreads, (numThreads == 1) ? "" : "s");

  //  Do the placing!

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi<RI->numReads()+1; fi++) {
    uint32   fiTigID = tigs.inUnitig(fi);

    if ((fiTigID == 0) ||                          //  Unplaced, don't care.
        (tigs[fiTigID]->_isUnassembled == true))   //  Unassembled, don't care.
      continue;

    //  Also in AS_BAT_MarkRepeatReads.C
    if (OG->isBubble(fi))                          //  Ignore bubble reads.
      continue;

    if (OG->isSpur(fi))                            //  Ignore spur reads.
      continue;

    //  Find ALL potential placements, regardless of error rate.

    uint32                     fiLen  = RI->readLength(fi);
    ufNode                    *fiRead = &tigs[fiTigID]->ufpath[ tigs.ufpathIdx(fi) ];
    int32                      fiMin  = fiRead->position.min();
    int32                      fiMax  = fiRead->position.max();
    vector<overlapPlacement>   placements;

    placeReadUsingOverlaps(tigs, NULL, fi, placements, placeRead_all, repeatLimit);

    //  For each placement decide if the overlap is compatible with the tig.

    for (uint32 pp=0; pp<placements.size(); pp++) {
      Unitig  *tig     = tigs[placements[pp].tigID];          //  Tig we're placed in.

      int32    utgmin  = placements[pp].position.min();       //  Position of placement in a unitig.
      int32    utgmax  = placements[pp].position.max();
      bool     utgfwd  = placements[pp].position.isForward();

      int32    ovlmin  = placements[pp].verified.min();       //  Placement in unitig, verified by overlaps.
      int32    ovlmax  = placements[pp].verified.max();

      bool     is5     = (placements[pp].covered.bgn == 0)     ? true : false;   //  Placement covers the 5' end of the read
      bool     is3     = (placements[pp].covered.end == fiLen) ? true : false;   //  Placement covers the 3' end of the read

      assert(placements[pp].covered.bgn < placements[pp].covered.end);     //  Coverage is always forward.


      //  Ignore placements in singletons, and placements that aren't
      //  overlaps (occurs when this read is a container for some small
      //  repeat read).
      //
      if ((tig->ufpath.size() <= 1) ||
          ((is5 == false) && (is3 == false))) {
        continue;
      }

      //  Ignore the placement if it is in the same place as the read came from.
      //
      if ((placements[pp].tigID == fiTigID) &&   //  Read placed in the same tig as it came from, and
          (utgmin <= fiMax) &&                   //  placement overlaps with the source position, therefore,
          (fiMin <= utgmax)) {                   //  this placement is where the read came from.
        logAGbuild(fi, pp, placements, "IDENTITY_PLACEMENT");
        continue;
      }

      //  Ignore the placement if it isn't compatible with the reads at the
      //  same location.
      //
      if (tig->overlapConsistentWithTig(deviationRepeat, ovlmin, ovlmax, placements[pp].erate()) < 0.5) {
        logAGbuild(fi, pp, placements, "HIGH_ERROR");
        continue;
      }

     if (placements[pp].fCoverage < 0.01) {
        logAGbuild(fi, pp, placements, "LOW_COVERAGE");
        continue;
     }
      //  A valid placement!
      //
      //  Decide if the overlap is to the:
      //    left  (towards 0) or
      //    right (towards infinity) of us on the tig.

      bool  onLeft  = (((utgfwd == true)  && (is5 == true)) ||
                       ((utgfwd == false) && (is3 == true)));

      bool  onRight = (((utgfwd == true)  && (is3 == true)) ||
                       ((utgfwd == false) && (is5 == true)));

      logAGbuild(fi, pp, placements, is5, is3, onLeft, onRight, "VALID_PLACEMENT");


      //  Find the reads we have overlaps to.
      //
      //  The range of reads here is the first and last read in the tig
      //  layout that overlaps with ourself.  We don't need to check that the
      //  reads overlap in the layout: the only false case I can think of
      //  involves contained reads.
      //
      //  READ:                         -----------------------------------
      //  TIG: Fidx  -----------------------------
      //  TIG: (1)      ------
      //  TIG:                 --------------------------------------
      //  TIG: Lidx                        -----------------------------------------
      //       (2)                                ------
      //
      //  The short read is placed at (1), but also has an overlap to us at (2).

      set<uint32>  tigReads;

      for (uint32 rr=placements[pp].tigFidx; rr <= placements[pp].tigLidx; rr++)
        tigReads.insert(tig->ufpath[rr].ident);

      //  Find the thickest overlap on each end.
      //
      //  Scan all overlaps.  Decide if the overlap is to the L or R of the
      //  _placed_ read, and save the thickest overlap on the 5' or 3' end of
      //  the read.

      uint32       no  = 0;
      BAToverlap  *ovl = OC->getOverlaps(fi, no);

      uint32  thickestC = UINT32_MAX, thickestCident = 0;
      uint32  thickest5 = UINT32_MAX, thickest5len   = 0;
      uint32  thickest3 = UINT32_MAX, thickest3len   = 0;

      for (uint32 oo=0; oo<no; oo++) {
        if (tigReads.count(ovl[oo].b_iid) == 0)   //  Don't care about overlaps to reads not in the set.
          continue;

        uint32  olapLen = RI->overlapLength(ovl[oo].a_iid, ovl[oo].b_iid, ovl[oo].a_hang, ovl[oo].b_hang);

        if      (ovl[oo].AisContainer() == true) {
          continue;
        }

        else if ((ovl[oo].AisContained() == true) && (is5 == true) && (is3 == true)) {
          if (thickestCident < ovl[oo].evalue) {
            thickestC      = oo;
            thickestCident = ovl[oo].evalue;
          }
        }

        else if ((ovl[oo].AEndIs5prime() == true) && (is5 == true)) {
          if (thickest5len < olapLen) {
            thickest5      = oo;
            thickest5len   = olapLen;
          }
        }

        else if ((ovl[oo].AEndIs3prime() == true) && (is3 == true)) {
          if (thickest3len < olapLen) {
            thickest3      = oo;
            thickest3len   = olapLen;
          }
        }
      }

      //  If we have both 5' and 3' edges, delete the containment edge.

      if ((thickest5 < no) &&
          (thickest3 < no)) {
        thickestC = UINT32_MAX;
      }

      //  If we have a containment edge, delete the 5' and 3' edges.

      if (thickestC < no) {
        thickest5 = UINT32_MAX;
        thickest3 = UINT32_MAX;
      }

      //  If no edges, log failure.

      if ((thickestC == UINT32_MAX) &&
          (thickest5 == UINT32_MAX) &&
          (thickest3 == UINT32_MAX)) {
        logAGbuild(fi, pp, placements, is5, is3, onLeft, onRight, "NO_EDGES");
        continue;
      }

      //  Create a new BestPlacement edge and save it on the list of placements for this read.

      BestPlacement  bp;

      bp.tigID     = placements[pp].tigID;

      bp.olapMin   = placements[pp].verified.min();
      bp.olapMax   = placements[pp].verified.max();

      if (thickestC < no)   bp.bestC = ovl[thickestC];
      if (thickest5 < no)   bp.best5 = ovl[thickest5];
      if (thickest3 < no)   bp.best3 = ovl[thickest3];

      if (bp.bestC.b_iid != 0) {                 //  Simple sanity check.  Ensure that contained
        assert(bp.best5.b_iid == 0);             //  edges have no dovetail edges.  This screws up
        assert(bp.best3.b_iid == 0);             //  the logic when outputting the graph.
      }

      assert((bp.bestC.a_hang <= 0) && (bp.bestC.b_hang >= 0));  //  ALL contained edges should be this.
      assert((bp.best5.a_hang <= 0) && (bp.best5.b_hang <= 0));  //  ALL 5' edges should be this.
      assert((bp.best3.a_hang >= 0) && (bp.best3.b_hang >= 0));  //  ALL 3' edges should be this.

      _pForward[fi].push_back(bp);

      //  Now just some logging of success.

      if (thickestC != UINT32_MAX)
        logAGbuild(fi, pp, placements, bp.bestC.b_iid, bp.best5.b_iid, bp.best3.b_iid, "CONTAINED");
      else
        logAGbuild(fi, pp, placements, bp.bestC.b_iid, bp.best5.b_iid, bp.best3.b_iid, "DOVETAIL");
    }  //  Over all placements
  }  //  Over all reads


  //  Make an index into BestPlacement.  Each read has a list of the read a
  //  BestPlacement comes from and the index of that placement.
  //
  //  Using this, you can get a list of all incoming edges to a given read.
  //    for (ii) {
  //      uint32  sourcei = _pReverse[fi].readID;
  //      uint32  sourcep = _pReverse[fi].placeID;
  //
  //      BestPlacement bp = _pForward[sourcei][sourcep];
  //    }
  //
  writeStatus("AssemblyGraph()-- building reverse edges.\n");

  for (uint32 fi=1; fi<RI->numReads()+1; fi++) {
    for (uint32 pp=0; pp<_pForward[fi].size(); pp++) {
      BestPlacement &bp = _pForward[fi][pp];
      BestReverse    br(fi, pp);

      if (bp.bestC.b_iid != 0)   _pReverse[bp.bestC.b_iid].push_back(br);
      if (bp.best5.b_iid != 0)   _pReverse[bp.best5.b_iid].push_back(br);
      if (bp.best3.b_iid != 0)   _pReverse[bp.best3.b_iid].push_back(br);
    }
  }

  writeStatus("AssemblyGraph()-- build complete.\n");
}
