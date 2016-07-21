
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
#include "AS_BAT_UnitigVector.H"

#include "intervalList.H"
#include "stddev.H"



AssemblyGraph::AssemblyGraph(const char   *prefix,
                             double        deviationGraph,
                             double        deviationBubble,
                             double        deviationRepeat,
                             UnitigVector &unitigs) {
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

  writeStatus("AssemblyGraph()-- finding edges for %u reads (%u contained), ignoring %u unplaced reads, with %d threads.\n",
              nToPlaceContained + nToPlace,
              nToPlaceContained,
              FI->numFragments() - nToPlaceContained - nToPlace,
              numThreads);
  writeStatus("AssemblyGraph()-- allocating vectors for placements, %.3fMB\n",   //  vector<> is 24 bytes, pretty tiny.
              sizeof(vector<BestPlacement>) * (fiLimit + 1) / 1048576.0);

  _placements = new vector<BestPlacement> [fiLimit + 1];

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

    placeFragUsingOverlaps(unitigs, AS_MAX_ERATE, NULL, fi, placements);

    //  For each placement decide if the overlap is compatible with the tig.
    //  Filter out placements that are too small.

    uint32   minOverlap = 500;

    for (uint32 i=0; i<placements.size(); i++) {
      Unitig *tig = unitigs[placements[i].tigID];

      //  Ignore placements in contains.
      if (tig->ufpath.size() == 1)
        continue;

      uint32  utgmin  = placements[i].position.min();  //  Placement in unitig.
      uint32  utgmax  = placements[i].position.max();

      uint32  ovlmin  = placements[i].verified.min();  //  Placement in unitig, verified by overlaps.
      uint32  ovlmax  = placements[i].verified.max();

      //  Ignore placements that are from thin overlaps.
      if (ovlmax - ovlmin < minOverlap)
        continue;

      //  Ignore placements that aren't overlaps (contained reads placed inside this read will do this).
      bool  is5 = (placements[i].covered.min() == 0)     ? true : false;
      bool  is3 = (placements[i].covered.max() == fiLen) ? true : false;

      if ((is5 == false) && (is3 == false)) {
        //writeLog("is5=false is3=false covered %d %d fiLen %d\n",
        //         placements[i].covered.min(), placements[i].covered.max(), fiLen);
        continue;
      }

      //  Decide if this is already captured in a tig.  If so, we'll emit to GFA, but omit from our
      //  internal graph.
      bool  isTig = false;

      if ((placements[i].tigID == fiTigID) && (utgmin <= fiMax) && (fiMin <= utgmax))
        isTig = true;

      //  Decide if the placement is complatible with the other reads in the tig.

      double  erate = placements[i].errors / placements[i].aligned;

      if ((isTig == false) &&
          (tig->overlapConsistentWithTig(5.0, ovlmin, ovlmax, erate) < 0.5)) {
        if ((enableLog == true) && (logFileFlagSet(LOG_PLACE_UNPLACED)))
          writeLog("frag %8u tested tig %6u (%6u reads) at %8u-%8u (cov %7.5f erate %6.4f) - HIGH ERROR\n",
                   fi, placements[i].tigID, tig->ufpath.size(), placements[i].position.bgn, placements[i].position.end, placements[i].fCoverage, erate);
        continue;
      }

      //  A valid placement.  Find the thickest overlaps on either side of the read.

      if (placements[i].tigFidx > placements[i].tigLidx)
        fprintf(stderr, "Invalid placement indices: tigFidx %u tigLidx %u\n", placements[i].tigFidx, placements[i].tigLidx);
      assert(placements[i].tigFidx <= placements[i].tigLidx);

      uint32  thickest5 = UINT32_MAX, thickest5len   = 0;
      uint32  thickest3 = UINT32_MAX, thickest3len   = 0;
      uint32  thickestC = UINT32_MAX, thickestCident = 0;

      uint32  n5 = 0;
      uint32  n3 = 0;

      for (uint32 rr=placements[i].tigFidx; rr <= placements[i].tigLidx; rr++) {
        ufNode  *read    = &tig->ufpath[rr];
        uint32   readmin = read->position.min();
        uint32   readmax = read->position.max();

        //  Can't use ourself as the thickest edge!
        if (read->ident == fi) {
          //writeLog("frag %u (len %u) utg %u-%u (%u) to read %u at %u-%u %u IDENT\n",
          //         fi, FI->fragmentLength(fi), utgmin, utgmax, utgmax-utgmin, read->ident, readmin, readmax, readmax-readmin);
          continue;
        }

        //  Is the read contained in us?
        if        ((utgmin <= readmin) && (readmax <= utgmax)) {
          //writeLog("frag %u (len %u) utg %u-%u (%u) to read %u at %u-%u %u contained\n",
          //         fi, FI->fragmentLength(fi), utgmin, utgmax, utgmax-utgmin, read->ident, readmin, readmax, readmax-readmin);
          continue;
        }

        //  Are we contained in the read?  We should save the 'best' overlap, but we don't know identities here.
        else if ((readmin <= utgmin) && (utgmax <= readmax)) {
          thickestC      = rr;
          thickestCident = 0;
        }

        //  Is the read to our right?
        else if ((readmin < utgmax) && (utgmax < readmax) && (utgmax - readmin > thickest3len)) {
          n3++;
          thickest3    = rr;
          thickest3len = utgmax - readmin;
        }

        //  Is the read to our left?
        else if ((readmin < utgmin) && (utgmin < readmax) && (readmax - utgmin > thickest5len)) {
          n5++;
          thickest5    = rr;
          thickest5len = readmax - utgmin;
        }

        //writeLog("frag %u (len %u) utg %u-%u (%u) to read %u at %u-%u %u\n",
        //         fi, FI->fragmentLength(fi), utgmin, utgmax, utgmax-utgmin, read->ident, readmin, readmax, readmax-readmin);

        //  Otherwise an inferior overlap.
      }

      //  Save the edge.

      bp.tigID     = placements[i].tigID;

      bp.placedBgn = placements[i].position.bgn;
      bp.placedEnd = placements[i].position.end;

      bp.olapBgn   = placements[i].verified.bgn;
      bp.olapEnd   = placements[i].verified.end;

      bp.thickestC = (thickestC != UINT32_MAX) ? tig->ufpath[thickestC].ident : 0;
      bp.thickest5 = (thickest5 != UINT32_MAX) ? tig->ufpath[thickest5].ident : 0;
      bp.thickest3 = (thickest3 != UINT32_MAX) ? tig->ufpath[thickest3].ident : 0;

      //  If there are best edges off the 5' or 3' end, grab all the overlaps, find the particular
      //  overlap, and generate new BestEdgeOverlaps for them.

      if ((bp.thickestC != 0) ||
          (bp.thickest5 != 0) ||
          (bp.thickest3 != 0)) {
        uint32       no  = 0;
        BAToverlap  *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

        for (uint32 oo=0; oo<no; oo++) {
          if ((ovl[oo].b_iid == bp.thickestC) && (ovl[oo].AisContained()))
            bp.bestC.set(ovl[oo]);

          if ((ovl[oo].b_iid == bp.thickest5) && (ovl[oo].AEndIs5prime()))
            bp.best5.set(ovl[oo]);

          if ((ovl[oo].b_iid == bp.thickest3) && (ovl[oo].AEndIs3prime()))
            bp.best3.set(ovl[oo]);
        }
      }

      if (isTig == true)
        _placements[fi].push_back(bp);

      //  And now just log.

#if 0
      if (bp.thickestC != 0) {
        writeLog("frag %8u -> tig %6u %6u reads  contained %6u     - at %8u-%-8u  cov %7.5f erate %6.4f olap %5u%s\n",
                 fi, bp.tigID, tig->ufpath.size(),
                 bp.thickestC,
                 bp.placedBgn, bp.placedEnd, placements[i].fCoverage, erate, ovlmax - ovlmin,
                 (isTig == true) ? " unitig" : "");
      } else {
        writeLog("frag %8u -> tig %6u %6u reads  dovetail %6u %6u at %8u-%-8u  cov %7.5f erate %6.4f olap %5u%s  n5=%u n3=%u\n",
                 fi, bp.tigID, tig->ufpath.size(),
                 bp.thickest5, bp.thickest3,
                 bp.placedBgn, bp.placedEnd, placements[i].fCoverage, erate, ovlmax - ovlmin,
                 (isTig == true) ? " unitig" : "", n5, n3);
      } 
#endif

   }  //  Over all placements
  }  //  Over all reads


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

  sprintf(N, "%s.%s.edges.gfa", prefix, label);

  BEG = fopen(N, "w");

  if (BEG == NULL)
    return;
     
  fprintf(BEG, "H\tVN:Z:bogart/edges\n");

  //  First, figure out what sequences are used.

  uint32   *used = new uint32 [FI->numFragments() + 1];

  memset(used, 0, sizeof(uint32) * (FI->numFragments() + 1));

  for (uint32 fi=1; fi<FI->numFragments() + 1; fi++) {
    if (_placements[fi].size() > 0)
      used[fi] = 1;

    for (uint32 pp=0; pp<_placements[fi].size(); pp++) {
      used[_placements[fi][pp].bestC.b_iid] = 1;
      used[_placements[fi][pp].best5.b_iid] = 1;
      used[_placements[fi][pp].best3.b_iid] = 1;
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

  uint32  nC = 0;
  uint32  n5 = 0;
  uint32  n3 = 0;

  for (uint32 fi=1; fi<FI->numFragments() + 1; fi++) {
    for (uint32 pp=0; pp<_placements[fi].size(); pp++) {
      BAToverlapInt *bestedgeC = &_placements[fi][pp].bestC;
      BAToverlapInt *bestedge5 = &_placements[fi][pp].best5;
      BAToverlapInt *bestedge3 = &_placements[fi][pp].best3;

      if (bestedgeC->b_iid != 0)
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

  writeStatus("AssemblyGraph()-- output %u contained edges, %u with dovetail edges\n", nC, nCdove);
  writeStatus("AssemblyGraph()-- output %u 5' edges\n", n5);
  writeStatus("AssemblyGraph()-- output %u 3' edges\n", n3);
}
