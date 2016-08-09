
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
 *    Brian P. Walenz beginning on 2016-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_AssemblyGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_TigVector.H"

#include "intervalList.H"
#include "stddev.H"

#undef  LOG_GRAPH      //  Report the (forward) graph edges in the log file


AssemblyGraph::AssemblyGraph(const char   *prefix,
                             double        deviationGraph,
                             double        deviationBubble,
                             double        deviationRepeat,
                             TigVector &unitigs) {
  uint32  fiLimit    = FI->numFragments();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  //  Just some logging.  Count the number of reads we try to place.

  uint32   nToPlaceContained = 0;
  uint32   nToPlace          = 0;
  uint32   nPlacedContained  = 0;
  uint32   nPlaced           = 0;
  uint32   nFailedContained  = 0;
  uint32   nFailed           = 0;

  for (uint32 fid=1; fid<FI->numFragments()+1; fid++) {
    if (Unitig::fragIn(fid) == 0)   //  Unplaced, don't care.  These didn't assemble, and aren't contained.
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

  writeStatus("AssemblyGraph()-- finding edges for %u reads (%u contained), ignoring %u unplaced reads, with %d threads.\n",
              nToPlaceContained + nToPlace,
              nToPlaceContained,
              FI->numFragments() - nToPlaceContained - nToPlace,
              numThreads);

  //  Do the placing!

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi<FI->numFragments()+1; fi++) {
    bool  enableLog = true;

    uint32   fiTigID = Unitig::fragIn(fi);

    if (fiTigID == 0)  //  Unplaced, don't care.
      continue;

    //  Grab a bit about this read.

    uint32   fiLen  = FI->fragmentLength(fi);
    ufNode  *fiRead = &unitigs[fiTigID]->ufpath[ Unitig::pathPosition(fi) ];
    uint32   fiMin  = fiRead->position.min();
    uint32   fiMax  = fiRead->position.max();

    BestPlacement  bp;  //  Result of the placement, saved in the master vector at the end.

    //  Find ALL potential placements, regardless of error rate.

    vector<overlapPlacement>   placements;

    placeFragUsingOverlaps(unitigs, NULL, fi, placements);

#ifdef LOG_GRAPH
    writeLog("AG()-- working on frag %u with %u placements\n", fi, placements.size());
#endif

    //  For each placement decide if the overlap is compatible with the tig.
    //  Filter out placements that are too small.

    //uint32   minOverlap = 500;

    for (uint32 pp=0; pp<placements.size(); pp++) {
      Unitig *tig = unitigs[placements[pp].tigID];

      double  erate = placements[pp].errors / placements[pp].aligned;

      //  Ignore placements in singletons.
      if (tig->ufpath.size() <= 1) {
#ifdef LOG_GRAPH
        writeLog("AG()-- frag %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f SINGLETON\n",
                 fi, pp,
                 placements[pp].tigID,
                 placements[pp].position.bgn, placements[pp].position.end,
                 placements[pp].verified.bgn, placements[pp].verified.end,
                 placements[pp].fCoverage, erate);
#endif
        continue;
      }

      uint32   utgmin  = placements[pp].position.min();  //  Placement in unitig.
      uint32   utgmax  = placements[pp].position.max();
      bool     utgfwd  = placements[pp].position.isForward();

      uint32   ovlmin  = placements[pp].verified.min();  //  Placement in unitig, verified by overlaps.
      uint32   ovlmax  = placements[pp].verified.max();

      //  Ignore placements that are from thin overlaps.
      //if (ovlmax - ovlmin < minOverlap)
      //  continue;

      //  Decide if the overlap is on our 5' or 3' end.  If both, we've been aligned fully.
      bool  is5  = (placements[pp].covered.min() == 0)     ? true : false;
      bool  is3  = (placements[pp].covered.max() == fiLen) ? true : false;

      //  Ignore placements that aren't overlaps (contained reads placed inside this read will do this).
      if ((is5 == false) && (is3 == false)) {
#ifdef LOG_GRAPH
        writeLog("AG()-- frag %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f SPANNED_REPEAT\n",
                 fi, pp,
                 placements[pp].tigID,
                 placements[pp].position.bgn, placements[pp].position.end,
                 placements[pp].verified.bgn, placements[pp].verified.end,
                 placements[pp].fCoverage, erate);
#endif
        continue;
      }

      //  Decide if the overlap is to the left (towards 0) or right (towards infinity) of us on the tig.
      bool  onLeft  = (((utgfwd == true)  && (is5 == true)) ||
                       ((utgfwd == false) && (is3 == true))) ? true : false;

      bool  onRight = (((utgfwd == true)  && (is3 == true)) ||
                       ((utgfwd == false) && (is5 == true))) ? true : false;

      assert((is5 == true) || (is3 == true));


      //  Decide if this is already captured in a tig.  If so, we'll emit to GFA, but omit from our
      //  internal graph.
      bool  isTig = false;

      if ((placements[pp].tigID == fiTigID) && (utgmin <= fiMax) && (fiMin <= utgmax))
        isTig = true;

      //  Decide if the placement is complatible with the other reads in the tig.

#define REPEAT_FRACTION   0.5

      if ((isTig == false) &&
          (tig->overlapConsistentWithTig(deviationRepeat, ovlmin, ovlmax, erate) < REPEAT_FRACTION)) {
#ifdef LOG_GRAPH
        if ((enableLog == true) && (logFileFlagSet(LOG_PLACE_UNPLACED)))
          writeLog("AG()-- frag %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f HIGH_ERROR\n",
                   fi, pp,
                   placements[pp].tigID,
                   placements[pp].position.bgn, placements[pp].position.end,
                   placements[pp].verified.bgn, placements[pp].verified.end,
                   placements[pp].fCoverage, erate);
#endif
        continue;
      }

      //  A valid placement!

#if 0
      writeLog("AG()-- frag %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f Fidx %6u Lidx %6u is5 %d is3 %d onLeft %d onRight %d  VALID_PLACEMENT\n",
               fi, pp,
               placements[pp].tigID,
               placements[pp].position.bgn, placements[pp].position.end,
               placements[pp].verified.bgn, placements[pp].verified.end,
               placements[pp].fCoverage, erate,
               placements[pp].tigFidx, placements[pp].tigLidx,
               is5, is3, onLeft, onRight);
#endif


      //  Find the overlaps for this read.  We'll use this to find the thickest overlaps to this
      //  tig.  We'd like to use a map, to avoid searching later, but because of duplicates, the
      //  best we can (easily) do is just remember the ID's we have an overlap to.

      set<uint32>  olapIDs;
      uint32       no  = 0;
      BAToverlap  *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

      for (uint32 oo=0; oo<no; oo++)
        olapIDs.insert(ovl[oo].b_iid);

      //  Find the thickest overlaps on either side of the read.
      //
      //  The logic here is to first classify the type of overlap (contained, dovetail to the left,
      //  etc) then to decide if our read cares about these types of overlaps, and if so, save the
      //  best.

      uint32  thickest5 = 0, thickest5len   = 0;
      uint32  thickest3 = 0, thickest3len   = 0;
      uint32  thickestC = 0, thickestCident = 0;

      for (uint32 rr=placements[pp].tigFidx; rr <= placements[pp].tigLidx; rr++) {
        ufNode  *read    = &tig->ufpath[rr];
        uint32   readmin = read->position.min();
        uint32   readmax = read->position.max();

        //  No oevrlap?
        if (olapIDs.count(read->ident) == 0)
          continue;

        //  Can't use ourself as the thickest edge!  (useless test; it's captured above)
        if (read->ident == fi)
          continue;

        //  Is the read contained in us?
        if        ((ovlmin <= readmin) && (readmax <= ovlmax))
          continue;

        //  Are we contained in the read?  We should save the 'best' overlap, but we don't know identities here.
        else if ((readmin <= ovlmin) && (ovlmax <= readmax)) {
          if ((onLeft == true) && (onRight == true)) {
            //writeLog("AG()-- thickestC %u %d-%d\n", read->ident, read->position.bgn, read->position.end);
            thickestC      = read->ident;
            thickestCident = 0;
          }
        }

        //  Is the read to our left?
        else if ((readmin <= ovlmin) && (ovlmin <= readmax)) {
          if ((onLeft == true) && (readmax - ovlmin > thickest5len)) {
            //writeLog("AG()-- thickest5 %u len %u (was %u len %u) OVL %d %d READ %d %d\n", read->ident, readmax - ovlmin, thickest5, thickest5len, ovlmin, ovlmax, readmin, readmax);
            thickest5    = read->ident;
            thickest5len = readmax - ovlmin;
          }
        }

        //  Is the read to our right?
        else if ((readmin <= ovlmax) && (ovlmax <= readmax)) {
          if ((onRight == true) && (ovlmax - readmin > thickest3len)) {
            //writeLog("AG()-- thickest3 %u len %u (was %u len %u) OVL %d %d READ %d %d\n", read->ident, readmax - ovlmin, thickest3, thickest3len, ovlmin, ovlmax, readmin, readmax);
            thickest3    = read->ident;
            thickest3len = ovlmax - readmin;
          }
        }
      }

      //  Save the edge.

      bp.tigID     = placements[pp].tigID;

      bp.placedBgn = placements[pp].position.bgn;
      bp.placedEnd = placements[pp].position.end;

      bp.olapBgn   = placements[pp].verified.bgn;
      bp.olapEnd   = placements[pp].verified.end;

      bp.isRepeat  = false;
      bp.isBubble  = false;
      bp.isUnitig  = isTig;

      //  If there are best edges off the 5' or 3' end, grab all the overlaps, find the particular
      //  overlap, and generate new BestEdgeOverlaps for them.

      if ((thickestC == 0) &&
          (thickest5 == 0) &&
          (thickest3 == 0)) {
#ifdef LOG_GRAPH
          writeLog("AG()-- frag %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f NO_EDGES Fidx %6u Lidx %6u is5 %d is3 %d onLeft %d onRight %d\n",
                   fi, pp,
                   placements[pp].tigID,
                   placements[pp].position.bgn, placements[pp].position.end,
                   placements[pp].verified.bgn, placements[pp].verified.end,
                   placements[pp].fCoverage, erate,
                   placements[pp].tigFidx, placements[pp].tigLidx,
                   is5, is3, onLeft, onRight);
#endif
        continue;
      }
      assert((thickestC != 0) ||
             (thickest5 != 0) ||
             (thickest3 != 0));

      //  Find the actual overlap for each of the thickest

      bp.bestC = BAToverlapInt();
      bp.best5 = BAToverlapInt();
      bp.best3 = BAToverlapInt();

      for (uint32 oo=0; oo<no; oo++) {
        if ((ovl[oo].b_iid == thickestC) && (ovl[oo].AisContained()))
          bp.bestC.set(ovl[oo]);

        if ((ovl[oo].b_iid == thickest5) && (ovl[oo].AEndIs5prime()))
          bp.best5.set(ovl[oo]);

        if ((ovl[oo].b_iid == thickest3) && (ovl[oo].AEndIs3prime()))
          bp.best3.set(ovl[oo]);
      }

      //  Save the BestPlacement

      uint32       ff = _pForward[fi].size();

      _pForward[fi].push_back(bp);

      //  And now just log.

#ifdef LOG_GRAPH
      if (thickestC != 0) {
        writeLog("AG()-- frag %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f CONTAINED %8d (%8d %8d)%s\n",
                 fi, pp,
                 placements[pp].tigID,
                 placements[pp].position.bgn, placements[pp].position.end,
                 placements[pp].verified.bgn, placements[pp].verified.end,
                 placements[pp].fCoverage, erate,
                 thickestC, thickest5, thickest3,
                 (isTig == true) ? " IN_UNITIG" : "");
      } else {
        writeLog("AG()-- frag %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f DOVETAIL (%8d) %8d %8d%s\n",
                 fi, pp,
                 placements[pp].tigID,
                 placements[pp].position.bgn, placements[pp].position.end,
                 placements[pp].verified.bgn, placements[pp].verified.end,
                 placements[pp].fCoverage, erate,
                 thickestC, thickest5, thickest3,
                 (isTig == true) ? " IN_UNITIG" : "");
      } 
#endif
   }  //  Over all placements
  }  //  Over all reads


  //  Create the reverse links.  This can't be done in the parallel loop above without synchronization.

  writeStatus("AssemblyGraph()-- building reverse edges.\n",
              nToPlaceContained + nToPlace,
              nToPlaceContained,
              FI->numFragments() - nToPlaceContained - nToPlace,
              numThreads);

  for (uint32 fi=1; fi<FI->numFragments()+1; fi++) {
    for (uint32 ff=0; ff<_pForward[fi].size(); ff++) {
      BestPlacement &bp = _pForward[fi][ff];
      BestReverse    br(fi, ff);

      if (bp.bestC.b_iid != 0)
        _pReverse[bp.bestC.b_iid].push_back(br);

      if (bp.best5.b_iid != 0)
        _pReverse[bp.best5.b_iid].push_back(br);

      if (bp.best3.b_iid != 0)
        _pReverse[bp.best3.b_iid].push_back(br);
    }
  }

  writeStatus("AssemblyGraph()-- complete.\n");

  //  Cleanup.

  //writeStatus("placeContains()-- Placed %u contained reads and %u unplaced reads.\n", nPlacedContained, nPlaced);
  //writeStatus("placeContains()-- Failed to place %u contained reads (too high error suspected) and %u unplaced reads (lack of overlaps suspected).\n", nFailedContained, nFailed);
}





//  SWIPED FROM BestOverlapGraph::reportBestEdges

void
AssemblyGraph::reportGraph(const char *prefix, const char *label) {
 char   N[FILENAME_MAX];
  FILE *BEG = NULL;

  writeStatus("AssemblyGraph()-- generating '%s.%s.edges.gfa'.\n", prefix, label);

  sprintf(N, "%s.%s.assembly.gfa", prefix, label);

  BEG = fopen(N, "w");

  if (BEG == NULL)
    return;
     
  fprintf(BEG, "H\tVN:Z:bogart/edges\n");

  //  First, figure out what sequences are used.

  uint32   *used = new uint32 [FI->numFragments() + 1];

  memset(used, 0, sizeof(uint32) * (FI->numFragments() + 1));

  for (uint32 fi=1; fi<FI->numFragments() + 1; fi++) {

    for (uint32 pp=0; pp<_pForward[fi].size(); pp++)
      if (_pForward[fi][pp].isRepeat == false)
        used[fi] = 1;

    //if (_pForward[fi].size() > 0)
    //  used[fi] = 1;

    for (uint32 pp=0; pp<_pForward[fi].size(); pp++) {
      used[_pForward[fi][pp].bestC.b_iid] = 1;
      used[_pForward[fi][pp].best5.b_iid] = 1;
      used[_pForward[fi][pp].best3.b_iid] = 1;
    }
  }

  //  Then write those sequences.

  for (uint32 fi=1; fi<FI->numFragments() + 1; fi++)
    if (used[fi] == 1)
      fprintf(BEG, "S\tread%08u\t*\tLN:i:%u\n", fi, FI->fragmentLength(fi));

  //  Now, report edges.  GFA wants edges in exactly this format:
  //
  //       -------------
  //             -------------
  //
  //  with read orientation given by +/-.  Conveniently, this is what we've saved (for the edges).

  uint32  nConly = 0;  //  Number of placements with a contained overlap and no pair of edges
  uint32  nCdove = 0;  //  Number of placements with a cotnained overlap but a  pair of edges

  uint32  nRepeat = 0;

  uint32  nC = 0;
  uint32  n5 = 0;
  uint32  n3 = 0;

  for (uint32 fi=1; fi<FI->numFragments() + 1; fi++) {
    for (uint32 pp=0; pp<_pForward[fi].size(); pp++) {

      //  Skip eges to resolved repeats.
      if (_pForward[fi][pp].isRepeat == true) {
        nRepeat++;
        continue;
      }

      BAToverlapInt *bestedgeC = &_pForward[fi][pp].bestC;
      BAToverlapInt *bestedge5 = &_pForward[fi][pp].best5;
      BAToverlapInt *bestedge3 = &_pForward[fi][pp].best3;

      if (bestedgeC->b_iid != 0)  //  Not ambiguous.
        if ((bestedge5->b_iid == 0) && (bestedge3->b_iid == 0))
          nConly++;
        else
          nCdove++;

      if ((bestedgeC->b_iid != 0) && (bestedge5->b_iid == 0) && (bestedge3->b_iid == 0)) {
        int32  olaplen = FI->fragmentLength(fi);

        assert((bestedgeC->a_hang <= 0) && (bestedgeC->b_hang >= 0));  //  ALL contained edges should be this.

        nC++;
        fprintf(BEG, "C\tread%08u\t+\tread%08u\t%c\t%u\t%uM\n",
                fi,
                bestedgeC->b_iid, bestedgeC->flipped ? '-' : '+',
                -bestedgeC->a_hang,
                olaplen);
      }

      if (bestedge5->b_iid != 0) {
        int32  olaplen = FI->overlapLength(fi, bestedge5->b_iid, bestedge5->a_hang, bestedge5->b_hang);

        assert((bestedge5->a_hang <= 0) && (bestedge5->b_hang <= 0));  //  ALL 5' edges should be this.

        n5++;
        fprintf(BEG, "L\tread%08u\t-\tread%08u\t%c\t%uM 5p\n",
                fi,
                bestedge5->b_iid, bestedge5->BEndIs3prime() ? '-' : '+',
                olaplen);
      }

      if (bestedge3->b_iid != 0) {
        int32  olaplen = FI->overlapLength(fi, bestedge3->b_iid, bestedge3->a_hang, bestedge3->b_hang);

        assert((bestedge3->a_hang >= 0) && (bestedge3->b_hang >= 0));  //  ALL 3' edges should be this.

        n3++;
        fprintf(BEG, "L\tread%08u\t+\tread%08u\t%c\t%uM 3p\n",
                fi,
                bestedge3->b_iid, bestedge3->BEndIs3prime() ? '-' : '+',
                FI->overlapLength(fi, bestedge3->b_iid, bestedge3->a_hang, bestedge3->b_hang));
      }
    }
  }

  fclose(BEG);

  writeStatus("AssemblyGraph()-- output %u contained edges, %u with dovetail edges and %u contained only\n", nC, nCdove, nConly);
  writeStatus("AssemblyGraph()-- output %u 5' edges\n", n5);
  writeStatus("AssemblyGraph()-- output %u 3' edges\n", n3);
  writeStatus("AssemblyGraph()-- skipped %u reads to resolved repeat regions\n");
}
