
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
                             TigVector    &tigs) {
  uint32  fiLimit    = FI->numReads();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize  = (fiLimit < 100 * numThreads) ? numThreads : fiLimit / 99;

  //  Just some logging.  Count the number of reads we try to place.

  uint32   nToPlaceContained = 0;
  uint32   nToPlace          = 0;
  uint32   nPlacedContained  = 0;
  uint32   nPlaced           = 0;
  uint32   nFailedContained  = 0;
  uint32   nFailed           = 0;

  for (uint32 fid=1; fid<FI->numReads()+1; fid++) {
    if (Unitig::readIn(fid) == 0)   //  Unplaced, don't care.  These didn't assemble, and aren't contained.
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
              FI->numReads() - nToPlaceContained - nToPlace,
              numThreads);

  //  Do the placing!

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 fi=1; fi<FI->numReads()+1; fi++) {
    bool  enableLog = true;

    uint32   fiTigID = Unitig::readIn(fi);

    if (fiTigID == 0)  //  Unplaced, don't care.
      continue;

    //  Grab a bit about this read.

    uint32   fiLen  = FI->readLength(fi);
    ufNode  *fiRead = &tigs[fiTigID]->ufpath[ Unitig::pathPosition(fi) ];
    uint32   fiMin  = fiRead->position.min();
    uint32   fiMax  = fiRead->position.max();

    //  Find ALL potential placements, regardless of error rate.

    vector<overlapPlacement>   placements;

    placeReadUsingOverlaps(tigs, NULL, fi, placements);

#ifdef LOG_GRAPH
    //writeLog("AG()-- working on read %u with %u placements\n", fi, placements.size());
#endif

    //  For each placement decide if the overlap is compatible with the tig.

    for (uint32 pp=0; pp<placements.size(); pp++) {
      Unitig *tig = tigs[placements[pp].tigID];

      double  erate = placements[pp].errors / placements[pp].aligned;

      //  Ignore placements in singletons.
      if (tig->ufpath.size() <= 1) {
#if 0
#ifdef LOG_GRAPH
        writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f SINGLETON\n",
                 fi, pp,
                 placements[pp].tigID,
                 placements[pp].position.bgn, placements[pp].position.end,
                 placements[pp].verified.bgn, placements[pp].verified.end,
                 placements[pp].fCoverage, erate);
#endif
#endif
        continue;
      }

      uint32   utgmin  = placements[pp].position.min();  //  Placement in unitig.
      uint32   utgmax  = placements[pp].position.max();
      bool     utgfwd  = placements[pp].position.isForward();

      uint32   ovlmin  = placements[pp].verified.min();  //  Placement in unitig, verified by overlaps.
      uint32   ovlmax  = placements[pp].verified.max();

      assert(placements[pp].covered.bgn < placements[pp].covered.end);     //  Coverage is always forward.

      bool  is5  = (placements[pp].covered.bgn == 0)     ? true : false;   //  Placement covers the 5' end of the read
      bool  is3  = (placements[pp].covered.end == fiLen) ? true : false;   //  Placement covers the 3' end of the read


      //  Ignore placements that aren't overlaps (contained reads placed inside this read will do this).
      if ((is5 == false) && (is3 == false)) {
#if 0
#ifdef LOG_GRAPH
        writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f SPANNED_REPEAT\n",
                 fi, pp,
                 placements[pp].tigID,
                 placements[pp].position.bgn, placements[pp].position.end,
                 placements[pp].verified.bgn, placements[pp].verified.end,
                 placements[pp].fCoverage, erate);
#endif
#endif
        continue;
      }

      //  Decide if the overlap is to the left (towards 0) or right (towards infinity) of us on the tig.
      bool  onLeft  = (((utgfwd == true)  && (is5 == true)) ||
                       ((utgfwd == false) && (is3 == true))) ? true : false;

      bool  onRight = (((utgfwd == true)  && (is3 == true)) ||
                       ((utgfwd == false) && (is5 == true))) ? true : false;


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
          writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f HIGH_ERROR\n",
                   fi, pp,
                   placements[pp].tigID,
                   placements[pp].position.bgn, placements[pp].position.end,
                   placements[pp].verified.bgn, placements[pp].verified.end,
                   placements[pp].fCoverage, erate);
#endif
        continue;
      }

      //  A valid placement!  Create a BestPlacement for it.

      BestPlacement  bp;

#if 0
      writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f Fidx %6u Lidx %6u is5 %d is3 %d onLeft %d onRight %d  VALID_PLACEMENT\n",
               fi, pp,
               placements[pp].tigID,
               placements[pp].position.bgn, placements[pp].position.end,
               placements[pp].verified.bgn, placements[pp].verified.end,
               placements[pp].fCoverage, erate,
               placements[pp].tigFidx, placements[pp].tigLidx,
               is5, is3, onLeft, onRight);
#endif

      //  Find the reads we have overlaps to.  The range of reads here is the first and last read in
      //  the tig layout that overlaps with ourself.  We don't need to check that the reads overlap in the
      //  layout: the only false case I can think of involves contained reads.
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

      //  Scan all overlaps.  Decide if the overlap is to the L or R of the _placed_ read, and save
      //  the thickest overlap on the 5' or 3' end of the read.

      uint32       no  = 0;
      BAToverlap  *ovl = OC->getOverlaps(fi, AS_MAX_EVALUE, no);

      uint32  thickestC = 0, thickestCident = 0;
      uint32  thickest5 = 0, thickest5len   = 0;
      uint32  thickest3 = 0, thickest3len   = 0;

      for (uint32 oo=0; oo<no; oo++) {
        if (tigReads.count(ovl[oo].b_iid) == 0)   //  Don't care about overlaps to reads not in the set.
          continue;

        uint32  olapLen = FI->overlapLength(ovl[oo].a_iid, ovl[oo].b_iid, ovl[oo].a_hang, ovl[oo].b_hang);

        if      (ovl[oo].AisContainer() == true) {
          continue;
        }

        else if (ovl[oo].AisContained() == true) {
          if (thickestCident < ovl[oo].evalue) {
            thickestC      = oo;
            thickestCident = ovl[oo].evalue;
            bp.bestC.set(ovl[oo]);
          }
        }

        else if (ovl[oo].AEndIs5prime() == true) {
          if (thickest5len < olapLen) {
            thickest5      = oo;
            thickest5len   = olapLen;
            bp.best5.set(ovl[oo]);
          }
        }

        else if (ovl[oo].AEndIs3prime() == true) {
          if (thickest3len < olapLen) {
            thickest3      = oo;
            thickest3len   = olapLen;
            bp.best3.set(ovl[oo]);
          }
        }
      }

      //  If we have both 5' and 3' edges, delete the containment edge.

      if ((bp.best5.b_iid != 0) && (bp.best3.b_iid != 0)) {
        thickestC      = 0;
        thickestCident = 0;
        bp.bestC       = BAToverlapInt();
      }

      //  Save the edge.

      bp.tigID     = placements[pp].tigID;

      bp.placedBgn = placements[pp].position.bgn;
      bp.placedEnd = placements[pp].position.end;

      bp.olapBgn   = placements[pp].verified.bgn;
      bp.olapEnd   = placements[pp].verified.end;

      bp.isContig  = isTig;  //  Currently, contigs == unitigs.
      bp.isUnitig  = isTig;
      bp.isBubble  = false;
      bp.isRepeat  = false;

      //  If there are best edges off the 5' or 3' end, grab all the overlaps, find the particular
      //  overlap, and generate new BestEdgeOverlaps for them.

      if ((thickestC == 0) &&
          (thickest5 == 0) &&
          (thickest3 == 0)) {
#ifdef LOG_GRAPH
          writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f NO_EDGES Fidx %6u Lidx %6u is5 %d is3 %d onLeft %d onRight %d\n",
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

      //  Save the BestPlacement

      uint32       ff = _pForward[fi].size();

      _pForward[fi].push_back(bp);

      //  And now just log.

#ifdef LOG_GRAPH
      if (isTig == false)
        if (thickestC != 0) {
          writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f CONTAINED %8d (%8d %8d)%s\n",
                   fi, pp,
                   placements[pp].tigID,
                   placements[pp].position.bgn, placements[pp].position.end,
                   placements[pp].verified.bgn, placements[pp].verified.end,
                   placements[pp].fCoverage, erate,
                   bp.bestC.b_iid, bp.best5.b_iid, bp.best3.b_iid,
                   (isTig == true) ? " IN_UNITIG" : "");
        } else {
          writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f DOVETAIL (%8d) %8d %8d%s\n",
                   fi, pp,
                   placements[pp].tigID,
                   placements[pp].position.bgn, placements[pp].position.end,
                   placements[pp].verified.bgn, placements[pp].verified.end,
                   placements[pp].fCoverage, erate,
                   bp.bestC.b_iid, bp.best5.b_iid, bp.best3.b_iid,
                   (isTig == true) ? " IN_UNITIG" : "");
        } 
#endif
    }  //  Over all placements
  }  //  Over all reads


  //  Create the reverse links.  This can't be done in the parallel loop above without synchronization.

  writeStatus("AssemblyGraph()-- building reverse edges.\n",
              nToPlaceContained + nToPlace,
              nToPlaceContained,
              FI->numReads() - nToPlaceContained - nToPlace,
              numThreads);

  for (uint32 fi=1; fi<FI->numReads()+1; fi++) {
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

  //  First, figure out what sequences are used.  A sequence is used if it has forward edges,
  //  or if it is referred to by a forward edge.

  uint32   *used = new uint32 [FI->numReads() + 1];

  memset(used, 0, sizeof(uint32) * (FI->numReads() + 1));

  for (uint32 fi=1; fi<FI->numReads() + 1; fi++) {
    if (_pForward[fi].size() > 0)
      used[fi] = 1;

    for (uint32 pp=0; pp<_pForward[fi].size(); pp++) {
      used[_pForward[fi][pp].bestC.b_iid] = 1;
      used[_pForward[fi][pp].best5.b_iid] = 1;
      used[_pForward[fi][pp].best3.b_iid] = 1;
    }
  }

  //  Then write those sequences.

  for (uint32 fi=1; fi<FI->numReads() + 1; fi++)
    if (used[fi] == 1)
      fprintf(BEG, "S\tread%08u\t*\tLN:i:%u\n", fi, FI->readLength(fi));

  //  Now, report edges.  GFA wants edges in exactly this format:
  //
  //       -------------
  //             -------------
  //
  //  with read orientation given by +/-.  Conveniently, this is what we've saved (for the edges).

  uint32  nCdove0 = 0;
  uint32  nCdove1 = 0;
  uint32  nCdove2 = 0;

  uint64  nTig[3] = {0,0,0};  //  Number of edges - both contig and unitig
  uint64  nCtg[3] = {0,0,0};  //  Number of edges - contig only
  uint64  nUtg[3] = {0,0,0};  //  Number of edges - unitig only (should be zero)
  uint64  nAsm[3] = {0,0,0};  //  Number of edges - between contigs

  uint64  nBubble = 0;
  uint64  nRepeat = 0;

  for (uint32 fi=1; fi<FI->numReads() + 1; fi++) {
    for (uint32 pp=0; pp<_pForward[fi].size(); pp++) {
      BAToverlapInt *bestedgeC = &_pForward[fi][pp].bestC;
      BAToverlapInt *bestedge5 = &_pForward[fi][pp].best5;
      BAToverlapInt *bestedge3 = &_pForward[fi][pp].best3;

      bool  outputC = (bestedgeC->b_iid != 0) && (bestedge5->b_iid == 0) && (bestedge3->b_iid == 0);
      bool  output5 = (bestedge5->b_iid != 0);
      bool  output3 = (bestedge3->b_iid != 0);

      //  Some statistics - number of contained reads with dovetail edges

      if (bestedgeC->b_iid != 0) {
        if      ((bestedge5->b_iid == 0) && (bestedge3->b_iid == 0))
          nCdove0++;  //  Contained edge with no dovetail edges

        else if ((bestedge5->b_iid != 0) && (bestedge3->b_iid != 0))
          nCdove2++;  //  Contained edge with both dovetail edges

        else
          nCdove1++;  //  Contained edge with one dovetail edge
      }

      //  Some statistics - number of edges of each type (in a contig, in a unitig, in both (tig), in neither (asm))

      if ((_pForward[fi][pp].isContig == true)  && (_pForward[fi][pp].isUnitig == true)) {
        if (outputC == true)   nTig[0]++;
        if (output5 == true)   nTig[1]++;
        if (output3 == true)   nTig[2]++;
      }

      if ((_pForward[fi][pp].isContig == true)  && (_pForward[fi][pp].isUnitig == false)) {
        if (outputC == true)   nCtg[0]++;
        if (output5 == true)   nCtg[1]++;
        if (output3 == true)   nCtg[2]++;
      }

      if ((_pForward[fi][pp].isContig == false) && (_pForward[fi][pp].isUnitig == true)) {
        if (outputC == true)   nUtg[0]++;
        if (output5 == true)   nUtg[1]++;
        if (output3 == true)   nUtg[2]++;
      }

      if ((_pForward[fi][pp].isContig == false) && (_pForward[fi][pp].isUnitig == false)) {
        if (outputC == true)   nAsm[0]++;
        if (output5 == true)   nAsm[1]++;
        if (output3 == true)   nAsm[2]++;
      }

      if (_pForward[fi][pp].isBubble == true) {
        nBubble++;
        continue;
      }

      if (_pForward[fi][pp].isRepeat == true) {
        nRepeat++;
        continue;
      }

      //  Check sanity.

      if (outputC)  assert((bestedgeC->a_hang <= 0) && (bestedgeC->b_hang >= 0));  //  ALL contained edges should be this.
      if (output5)  assert((bestedge5->a_hang <= 0) && (bestedge5->b_hang <= 0));  //  ALL 5' edges should be this.
      if (output3)  assert((bestedge3->a_hang >= 0) && (bestedge3->b_hang >= 0));  //  ALL 3' edges should be this.

      //  Finally, output the edge.

      if (outputC)
        fprintf(BEG, "C\tread%08u\t+\tread%08u\t%c\t%u\t%uM\tic:i:%d\tiu:i:%d\tib:i:%d\tir:i:%d\n",
                fi,
                bestedgeC->b_iid, bestedgeC->flipped ? '-' : '+',
                -bestedgeC->a_hang,
                FI->readLength(fi),
                _pForward[fi][pp].isContig,
                _pForward[fi][pp].isUnitig,
                _pForward[fi][pp].isBubble,
                _pForward[fi][pp].isRepeat);

      if (output5)
        fprintf(BEG, "L\tread%08u\t-\tread%08u\t%c\t%uM\tic:i:%d\tiu:i:%d\tib:i:%d\tir:i:%d\n",
                fi,
                bestedge5->b_iid, bestedge5->BEndIs3prime() ? '-' : '+',
                FI->overlapLength(fi, bestedge5->b_iid, bestedge5->a_hang, bestedge5->b_hang),
                _pForward[fi][pp].isContig,
                _pForward[fi][pp].isUnitig,
                _pForward[fi][pp].isBubble,
                _pForward[fi][pp].isRepeat);

      if (output3)
        fprintf(BEG, "L\tread%08u\t+\tread%08u\t%c\t%uM\tic:i:%d\tiu:i:%d\tib:i:%d\tir:i:%d\n",
                fi,
                bestedge3->b_iid, bestedge3->BEndIs3prime() ? '-' : '+',
                FI->overlapLength(fi, bestedge3->b_iid, bestedge3->a_hang, bestedge3->b_hang),
                _pForward[fi][pp].isContig,
                _pForward[fi][pp].isUnitig,
                _pForward[fi][pp].isBubble,
                _pForward[fi][pp].isRepeat);
    }
  }

  fclose(BEG);

  //  And report statistics.

  writeStatus("AssemblyGraph()-- %8"F_U64P" contained edges with zero dovetail edges (output as a contained edge)\n", nCdove0);
  writeStatus("AssemblyGraph()-- %8"F_U64P" contained edges with one  dovetail edge  (output as a single dovetail edge)\n",  nCdove1);
  writeStatus("AssemblyGraph()-- %8"F_U64P" contained edges with two  dovetail edges (output as a pair of dovetail edges)\n", nCdove2);
  writeStatus("\n");
  writeStatus("AssemblyGraph()-- %8"F_U64P" bubble placements\n", nBubble);
  writeStatus("AssemblyGraph()-- %8"F_U64P" repeat placements\n", nRepeat);
  writeStatus("\n");
  writeStatus("AssemblyGraph()-- Intratig edges:     %8"F_U64P" contained  %8"F_U64P" 5'  %8"F_U64P" 3' (in both contig and unitig)\n", nTig[0], nTig[1], nTig[2]);
  writeStatus("AssemblyGraph()-- Contig only edges:  %8"F_U64P" contained  %8"F_U64P" 5'  %8"F_U64P" 3'\n", nCtg[0], nCtg[1], nCtg[2]);
  writeStatus("AssemblyGraph()-- Unitig only edges:  %8"F_U64P" contained  %8"F_U64P" 5'  %8"F_U64P" 3'\n", nUtg[0], nUtg[1], nUtg[2]);
  writeStatus("AssemblyGraph()-- Intercontig edges:  %8"F_U64P" contained  %8"F_U64P" 5'  %8"F_U64P" 3' (in neither contig nor unitig)\n", nAsm[0], nAsm[1], nAsm[2]);
}
