
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


void
AssemblyGraph::buildReverseEdges(void) {

  writeStatus("AssemblyGraph()-- building reverse edges.\n");

  for (uint32 fi=1; fi<RI->numReads()+1; fi++)
    _pReverse[fi].clear();

  for (uint32 fi=1; fi<RI->numReads()+1; fi++) {
    for (uint32 ff=0; ff<_pForward[fi].size(); ff++) {
      BestPlacement &bp = _pForward[fi][ff];
      BestReverse    br(fi, ff);

      //  Ensure that contained edges have no dovetail edges.  This screws up the logic when
      //  rebuilding and outputting the graph.

      if (bp.bestC.b_iid != 0) {
        assert(bp.best5.b_iid == 0);
        assert(bp.best3.b_iid == 0);
      }

      //  Add reverse edges if the forward edge exists

      if (bp.bestC.b_iid != 0)   _pReverse[bp.bestC.b_iid].push_back(br);
      if (bp.best5.b_iid != 0)   _pReverse[bp.best5.b_iid].push_back(br);
      if (bp.best3.b_iid != 0)   _pReverse[bp.best3.b_iid].push_back(br);

      //  Check sanity.

      assert((bp.bestC.a_hang <= 0) && (bp.bestC.b_hang >= 0));  //  ALL contained edges should be this.
      assert((bp.best5.a_hang <= 0) && (bp.best5.b_hang <= 0));  //  ALL 5' edges should be this.
      assert((bp.best3.a_hang >= 0) && (bp.best3.b_hang >= 0));  //  ALL 3' edges should be this.
    }
  }
}



void
AssemblyGraph::buildGraph(const char   *UNUSED(prefix),
                          double        deviationRepeat,
                          TigVector    &tigs,
                          bool          tigEndsOnly) {
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
    if (tigs.inUnitig(fid) == 0)   //  Unplaced, don't care.  These didn't assemble, and aren't contained.
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
    bool  enableLog = true;

    uint32   fiTigID = tigs.inUnitig(fi);

    if (fiTigID == 0)  //  Unplaced, don't care.
      continue;

    if (tigEndsOnly == true) {
      uint32 f = tigs[fiTigID]->firstRead()->ident;
      uint32 l = tigs[fiTigID]->lastRead()->ident;

      if ((f != fi) && (l != fi))    //  Not the first read and not the last read,
        continue;                    //  Don't care.
    }

    //  Grab a bit about this read.

    uint32   fiLen  = RI->readLength(fi);
    ufNode  *fiRead = &tigs[fiTigID]->ufpath[ tigs.ufpathIdx(fi) ];
    int32    fiMin  = fiRead->position.min();
    int32    fiMax  = fiRead->position.max();

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
#ifdef LOG_GRAPH
        writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f SINGLETON\n",
                 fi, pp,
                 placements[pp].tigID,
                 placements[pp].position.bgn, placements[pp].position.end,
                 placements[pp].verified.bgn, placements[pp].verified.end,
                 placements[pp].fCoverage, erate);
#endif
        continue;
      }

      int32    utgmin  = placements[pp].position.min();  //  Placement in unitig.
      int32    utgmax  = placements[pp].position.max();
      bool     utgfwd  = placements[pp].position.isForward();

      int32    ovlmin  = placements[pp].verified.min();  //  Placement in unitig, verified by overlaps.
      int32    ovlmax  = placements[pp].verified.max();

      assert(placements[pp].covered.bgn < placements[pp].covered.end);     //  Coverage is always forward.

      bool  is5  = (placements[pp].covered.bgn == 0)     ? true : false;   //  Placement covers the 5' end of the read
      bool  is3  = (placements[pp].covered.end == fiLen) ? true : false;   //  Placement covers the 3' end of the read


      //  Ignore placements that aren't overlaps (contained reads placed inside this read will do this).
      if ((is5 == false) && (is3 == false)) {
#ifdef LOG_GRAPH_ALL
        writeLog("AG()-- read %8u placement %2u -> tig %7u placed %9d-%9d verified %9d-%9d cov %7.5f erate %6.4f SPANNED_REPEAT\n",
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


      //  Decide if this is already captured in a tig.  If so, we'll emit to GFA, but omit from our
      //  internal graph.
      bool  isTig = false;

      if ((placements[pp].tigID == fiTigID) && (utgmin <= fiMax) && (fiMin <= utgmax))
        isTig = true;

      //  Decide if the placement is complatible with the other reads in the tig.

#define REPEAT_FRACTION   0.5

      if ((isTig == false) &&
          (tig->overlapConsistentWithTig(deviationRepeat, ovlmin, ovlmax, erate) < REPEAT_FRACTION)) {
#ifdef LOG_GRAPH_ALL
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

#ifdef LOG_GRAPH
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
            bp.bestC       = ovl[oo];
          }
        }

        else if ((ovl[oo].AEndIs5prime() == true) && (is5 == true)) {
          if (thickest5len < olapLen) {
            thickest5      = oo;
            thickest5len   = olapLen;
            bp.best5       = ovl[oo];
          }
        }

        else if ((ovl[oo].AEndIs3prime() == true) && (is3 == true)) {
          if (thickest3len < olapLen) {
            thickest3      = oo;
            thickest3len   = olapLen;
            bp.best3       = ovl[oo];
          }
        }
      }

      //  If we have both 5' and 3' edges, delete the containment edge.

      if ((bp.best5.b_iid != 0) && (bp.best3.b_iid != 0)) {
        thickestC = UINT32_MAX;   thickestCident = 0;   bp.bestC = BAToverlap();
      }

      //  If we have a containment edge, delete the 5' and 3' edges.

      if (bp.bestC.b_iid != 0) {
        thickest5 = UINT32_MAX;   thickest5len = 0;     bp.best5 = BAToverlap();
        thickest3 = UINT32_MAX;   thickest3len = 0;     bp.best3 = BAToverlap();
      }


      //  Save the edge.

      bp.tigID     = placements[pp].tigID;

      bp.placedBgn = placements[pp].position.bgn;
      bp.placedEnd = placements[pp].position.end;

      bp.olapBgn   = placements[pp].verified.bgn;
      bp.olapEnd   = placements[pp].verified.end;

      bp.isContig  = isTig;
      bp.isUnitig  = false;
      bp.isBubble  = false;
      bp.isRepeat  = false;

      //  If there are best edges off the 5' or 3' end, grab all the overlaps, find the particular
      //  overlap, and generate new BestEdgeOverlaps for them.

      if ((thickestC == UINT32_MAX) &&
          (thickest5 == UINT32_MAX) &&
          (thickest3 == UINT32_MAX)) {
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
      if (thickestC != UINT32_MAX) {
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

  buildReverseEdges();

  writeStatus("AssemblyGraph()-- build complete.\n");
}






void
placeAsContained(TigVector     &tigs,
                 uint32         fi,
                 BestPlacement &bp) {
  BestEdgeOverlap   edge(bp.bestC);
  ufNode            read;
  Unitig           *tig = tigs[ tigs.inUnitig(edge.readId()) ];

  if (tig->placeRead(read, fi, bp.bestC.AEndIs3prime(), &edge) == false) {
    fprintf(stderr, "WARNING: placeAsContained failed for fi=%u\n", fi);
    assert(0);
  }

  bp.tigID     = tig->id();

  bp.placedBgn = read.position.bgn;
  bp.placedEnd = read.position.end;

  bp.olapBgn = INT32_MIN;  //  We don't know the overlapping region (without a lot
  bp.olapEnd = INT32_MAX;  //  of work) so make it invalid.

  bp.isContig = (tigs.inUnitig(fi) == tigs.inUnitig(edge.readId()));
}



//  This test is correct, but it isn't used correctly.  When rebuilding the graph, we don't know if
//  a read is fully covered.  If it isn't fully covered, it isn't 'inContig' even if the positions
//  overlap.
bool
areReadsOverlapping(TigVector  &tigs,
                    uint32      ai,
                    uint32      bi) {
  Unitig  *at = tigs[ tigs.inUnitig(ai) ];
  Unitig  *bt = tigs[ tigs.inUnitig(bi) ];

  if (at != bt)
    return(false);

  ufNode  &ar = at->ufpath[ tigs.ufpathIdx(ai) ];
  ufNode  &br = bt->ufpath[ tigs.ufpathIdx(bi) ];

  return((ar.position.min() < br.position.max()) &&
         (br.position.min() < ar.position.max()));
}



void
placeAsDovetail(TigVector     &tigs,
                uint32         fi,
                BestPlacement &bp) {
  BestEdgeOverlap   edge5(bp.best5),  edge3(bp.best3);
  ufNode            read5,            read3;

  if ((bp.best5.b_iid > 0) && (bp.best3.b_iid > 0)) {
    Unitig  *tig5 = tigs[ tigs.inUnitig(edge5.readId()) ];
    Unitig  *tig3 = tigs[ tigs.inUnitig(edge3.readId()) ];

    assert(tig5->id() == tig3->id());

    if ((tig5->placeRead(read5, fi, bp.best5.AEndIs3prime(), &edge5) == false) ||
        (tig3->placeRead(read3, fi, bp.best3.AEndIs3prime(), &edge3) == false)) {
      fprintf(stderr, "WARNING: placeAsDovetail 5' 3' failed for fi=%u\n", fi);
      assert(0);
    }

    bp.tigID     = tig5->id();
    bp.placedBgn = (read5.position.bgn + read3.position.bgn) / 2;
    bp.placedEnd = (read5.position.end + read3.position.end) / 2;

#if 0
    bp.isContig  = (areReadsOverlapping(tigs, fi, bp.best5.b_iid) &&
                    areReadsOverlapping(tigs, fi, bp.best3.b_iid));
#else
    if ((bp.isContig == true) &&              //  Remove the isContig mark if this read is now
        (tigs.inUnitig(fi) != bp.tigID))     //  in a different tig than the two edges (which is unlikely).
      bp.isContig = false;
#endif
  }

  else if (bp.best5.b_iid > 0) {
    Unitig  *tig5 = tigs[ tigs.inUnitig(edge5.readId()) ];

    if (tig5->placeRead(read5, fi, bp.best5.AEndIs3prime(), &edge5) == false) {
      fprintf(stderr, "WARNING: placeAsDovetail 5' failed for fi=%u\n", fi);
      assert(0);
    }

    bp.tigID     = tig5->id();
    bp.placedBgn = read5.position.bgn;
    bp.placedEnd = read5.position.end;

#if 0
    bp.isContig = areReadsOverlapping(tigs, fi, bp.best5.b_iid);
#else
    if ((bp.isContig == true) &&              //  Remove the isContig mark if this read is now
        (tigs.inUnitig(fi) != bp.tigID))     //  in a different tig than the edge.
      bp.isContig = false;
#endif
  }

  else if (bp.best3.b_iid > 0) {
    Unitig  *tig3 = tigs[ tigs.inUnitig(edge3.readId()) ];

    if (tig3->placeRead(read3, fi, bp.best3.AEndIs3prime(), &edge3) == false) {
      fprintf(stderr, "WARNING: placeAsDovetail 3' failed for fi=%u\n", fi);
      assert(0);
    }

    bp.tigID     = tig3->id();
    bp.placedBgn = read3.position.bgn;
    bp.placedEnd = read3.position.end;

#if 0
    bp.isContig  = areReadsOverlapping(tigs, fi, bp.best3.b_iid);
#else
    if ((bp.isContig == true) &&              //  Remove the isContig mark if this read is now
        (tigs.inUnitig(fi) != bp.tigID))     //  in a different tig than the edge.
      bp.isContig = false;
#endif
  }

  assert(tigs[bp.tigID] != NULL);

  bp.olapBgn = INT32_MIN;  //  We don't know the overlapping region (without a lot
  bp.olapEnd = INT32_MAX;  //  of work) so make it invalid.
}




void
AssemblyGraph::rebuildGraph(TigVector     &tigs) {

  writeStatus("AssemblyGraph()-- rebuilding\n");

  uint64   nContain = 0;
  uint64   nSame    = 0;
  uint64   nSplit   = 0;

  for (uint32 fi=1; fi<RI->numReads()+1; fi++) {
    for (uint32 ff=0; ff<_pForward[fi].size(); ff++) {
      BestPlacement   &bp = _pForward[fi][ff];

      //  Figure out which tig each of our three overlaps is in.

      uint32  t5 = (bp.best5.b_iid > 0) ? tigs.inUnitig(bp.best5.b_iid) : UINT32_MAX;
      uint32  t3 = (bp.best3.b_iid > 0) ? tigs.inUnitig(bp.best3.b_iid) : UINT32_MAX;

      //writeLog("AssemblyGraph()-- rebuilding read %u edge %u with overlaps %u %u %u\n",
      //         fi, ff, bp.bestC.b_iid, bp.best5.b_iid, bp.best3.b_iid);

      //  If a containment relationship, place it using the contain and update the placement.

      if (bp.bestC.b_iid > 0) {
        assert(bp.best5.b_iid == 0);
        assert(bp.best3.b_iid == 0);

        nContain++;
        placeAsContained(tigs, fi, bp);
      }

      //  Otherwise, dovetails.  If both overlapping reads are in the same tig, place it and update
      //  the placement.

      else if ((t5 == t3) ||           //  Both in the same tig
               (t5 == UINT32_MAX) ||   //  5' overlap isn't set
               (t3 == UINT32_MAX)) {   //  3' overlap isn't set
        nSame++;
        placeAsDovetail(tigs, fi, bp);
      }

      //  Otherwise, yikes, our overlapping reads are in different tigs!  We need to make new
      //  placements and delete the current one.

      else {
        BestPlacement   bp5 = bp;
        BestPlacement   bp3 = bp;

        bp5.best3 = BAToverlap();   //  Erase the 3' overlap
        bp3.best5 = BAToverlap();   //  Erase the 5' overlap

        assert(bp5.best5.b_iid != 0);  //  Overlap must exist!
        assert(bp3.best3.b_iid != 0);  //  Overlap must exist!

        nSplit++;
        placeAsDovetail(tigs, fi, bp5);
        placeAsDovetail(tigs, fi, bp3);

        //  Add the two placements to our list.  We let one placement overwrite the current
        //  placement, move the placement after that to the end of the list, and overwrite
        //  that placement with our other new one.

        uint32  ll = _pForward[fi].size();

        //  There's a nasty case when ff is the last currently on the list; there isn't an ff+1
        //  element to move to the end of the list.  So, we add a new element to the list -
        //  guaranteeing there is always an ff+1 element - then move, then replace.

        _pForward[fi].push_back(BestPlacement());

        _pForward[fi][ll] = _pForward[fi][ff+1];

        _pForward[fi][ff]   = bp5;
        _pForward[fi][ff+1] = bp3;

        //  Skip the edge we just added.

        ff++;
      }
    }
  }

  buildReverseEdges();

  writeStatus("AssemblyGraph()-- rebuild complete.\n");
}





//  Filter edges that originate from the middle of a tig.
//  Need to save interior edges as long as they are consistent with a boundary edge.

void
AssemblyGraph::filterEdges(TigVector     &tigs) {
  uint64  nUnitig = 0;
  uint64  nContig = 0;
  uint64  nBubble = 0;
  uint64  nRepeat = 0;

  uint64  nMiddleFiltered = 0, nMiddleReads = 0;
  uint64  nRepeatFiltered = 0, nRepeatReads = 0;

  uint64  nIntersecting = 0;

  uint64  nRepeatEdges = 0;
  uint64  nBubbleEdges = 0;

  writeStatus("AssemblyGraph()-- filtering edges\n");

  //  Mark edges that are from the interior of a tig as 'repeat'.

  for (uint32 fi=1; fi<RI->numReads()+1; fi++) {
    if (_pForward[fi].size() == 0)
      continue;

    uint32       tT     =  tigs.inUnitig(fi);
    Unitig      *tig    =  tigs[tT];
    ufNode      &read   =  tig->ufpath[tigs.ufpathIdx(fi)];

    bool         hadMiddle = false;

    for (uint32 ff=0; ff<_pForward[fi].size(); ff++) {
      BestPlacement   &bp = _pForward[fi][ff];

      //  Edges forming the tig are not repeats.

      if (bp.isUnitig == true)    continue;
      if (bp.isContig == true)    continue;

      //  Edges from the end of a tig are not repeats.

      if (((read.position.min() == 0)                && (read.position.isForward()) && (bp.best5.b_iid  > 0) && (bp.best3.b_iid == 0)) ||
          ((read.position.min() == 0)                && (read.position.isReverse()) && (bp.best5.b_iid == 0) && (bp.best3.b_iid  > 0)) ||
          ((read.position.max() == tig->getLength()) && (read.position.isForward()) && (bp.best5.b_iid == 0) && (bp.best3.b_iid  > 0)) ||
          ((read.position.max() == tig->getLength()) && (read.position.isReverse()) && (bp.best5.b_iid  > 0) && (bp.best3.b_iid == 0))) {
        nIntersecting++;
        continue;
      }

      nMiddleFiltered++;

      bp.isRepeat = true;
      hadMiddle   = true;
    }

    if (hadMiddle)
      nMiddleReads++;
  }

  //  Filter edges that hit too many tigs

  for (uint32 fi=1; fi<RI->numReads()+1; fi++) {
    if (_pForward[fi].size() == 0)
      continue;

    uint32       tT     =  tigs.inUnitig(fi);
    Unitig      *tig    =  tigs[tT];
    ufNode      &read   =  tig->ufpath[tigs.ufpathIdx(fi)];

    set<uint32>  hits;

    for (uint32 ff=0; ff<_pForward[fi].size(); ff++) {
      BestPlacement   &bp = _pForward[fi][ff];

      assert(bp.isUnitig == false);

      if (bp.isUnitig == true)   { continue; }   //  Skip edges that are in tigs
      if (bp.isContig == true)   { continue; }   //
      if (bp.isRepeat == true)   { continue; }   //  Skip edges that are already ignored

      hits.insert(bp.tigID);
    }

    //  If only a few other tigs are involved, keep all.

    if (hits.size() > 0)
      writeLog("AG()-- read %u in tig %u has edges to %u tigs\n", fi, tT, hits.size());


#ifdef FILTER_DENSE_BUBBLES_FROM_GRAPH
    if (hits.size() <= FILTER_DENSE_BUBBLES_THRESHOLD)
      continue;

    //  Otherwise, mark all edges as repeat.

    nRepeatReads++;

    for (uint32 ff=0; ff<_pForward[fi].size(); ff++) {
      BestPlacement   &bp = _pForward[fi][ff];

      assert(bp.isUnitig == false);

      if (bp.isUnitig == true)   { continue; }   //  Skip edges that are in tigs
      if (bp.isContig == true)   { continue; }   //
      if (bp.isRepeat == true)   { continue; }   //  Skip edges that are already ignored

      nRepeatFiltered++;

      bp.isRepeat = true;
    }
#endif
  }

  //  Generate statistics

  for (uint32 fi=1; fi<RI->numReads()+1; fi++) {
    for (uint32 ff=0; ff<_pForward[fi].size(); ff++) {
      BestPlacement   &bp = _pForward[fi][ff];

      if (bp.isUnitig == true)   { nUnitig++;  continue; }
      if (bp.isContig == true)   { nContig++;  continue; }
      if (bp.isRepeat == true)   { nRepeatEdges++;       }
      if (bp.isRepeat == false)  { nBubbleEdges++;       }
    }
  }

  //  Report

  writeStatus("AssemblyGraph()-- " F_U64 " contig edges and " F_U64 " unitig edges.\n", nContig, nUnitig);
  writeStatus("AssemblyGraph()-- " F_U64 " bubble edges and " F_U64 " repeat edges.\n", nBubble, nRepeat);
  writeStatus("AssemblyGraph()-- " F_U64 " middle contig edges filtered from " F_U64 " reads.\n", nMiddleFiltered, nMiddleReads);
  writeStatus("AssemblyGraph()-- " F_U64 " repeat end edges filtered from " F_U64 " reads.\n", nRepeatFiltered, nRepeatReads);
  writeStatus("AssemblyGraph()-- " F_U64 " repeat edges (not output).\n", nRepeatEdges);
  writeStatus("AssemblyGraph()-- " F_U64 " bubble edges.\n", nBubbleEdges);
  writeStatus("AssemblyGraph()-- " F_U64 " intersecting edges (from the end of a tig to somewhere else).\n", nIntersecting);
}






bool
reportReadGraph_reportEdge(TigVector      &tigs,
                           BestPlacement  &pf,
                           bool            skipBubble,
                           bool            skipRepeat,
                           bool           &reportC,
                           bool           &report5,
                           bool           &report3) {
  reportC = false;
  report5 = false;
  report3 = false;

  if ((skipBubble == true) && (pf.isBubble == true))
    return(false);

  if ((skipRepeat == true) && (pf.isRepeat == true))
    return(false);

  //  If the destination isunassembled, all edges are ignored.
  if ((tigs[pf.tigID] == NULL) || (tigs[pf.tigID]->_isUnassembled == true))
    return(false);

  reportC = (tigs.inUnitig(pf.bestC.b_iid) != 0) && (tigs[ tigs.inUnitig(pf.bestC.b_iid) ]->_isUnassembled == false);
  report5 = (tigs.inUnitig(pf.best5.b_iid) != 0) && (tigs[ tigs.inUnitig(pf.best5.b_iid) ]->_isUnassembled == false);
  report3 = (tigs.inUnitig(pf.best3.b_iid) != 0) && (tigs[ tigs.inUnitig(pf.best3.b_iid) ]->_isUnassembled == false);

  if ((reportC == false) &&
      (report5 == false) &&
      (report3 == false))
    return(false);

  return(true);
}


//  SWIPED FROM BestOverlapGraph::reportBestEdges

void
AssemblyGraph::reportReadGraph(TigVector &tigs, const char *prefix, const char *label) {
 char   N[FILENAME_MAX];
  FILE *BEG = NULL;

  bool  skipBubble      = true;
  bool  skipRepeat      = true;
  bool  skipUnassembled = true;

  uint64  nEdgeToUnasm = 0;

  writeStatus("AssemblyGraph()-- generating '%s.%s.assembly.gfa'.\n", prefix, label);

  snprintf(N, FILENAME_MAX, "%s.%s.assembly.gfa", prefix, label);

  BEG = fopen(N, "w");

  if (BEG == NULL)
    return;

  fprintf(BEG, "H\tVN:Z:bogart/edges\n");

  //  First, figure out what sequences are used.  A sequence is used if it has forward edges,
  //  or if it is referred to by a forward edge.

  uint32   *used = new uint32 [RI->numReads() + 1];

  memset(used, 0, sizeof(uint32) * (RI->numReads() + 1));

  for (uint32 fi=1; fi<RI->numReads() + 1; fi++) {
    for (uint32 pp=0; pp<_pForward[fi].size(); pp++) {
      BestPlacement  &pf = _pForward[fi][pp];
      bool            reportC=false, report5=false, report3=false;

      if ((tigs.inUnitig(pf.bestC.b_iid) != 0) && (tigs[ tigs.inUnitig(pf.bestC.b_iid) ]->_isUnassembled == true))
        nEdgeToUnasm++;
      if ((tigs.inUnitig(pf.best5.b_iid) != 0) && (tigs[ tigs.inUnitig(pf.best5.b_iid) ]->_isUnassembled == true))
        nEdgeToUnasm++;
      if ((tigs.inUnitig(pf.best3.b_iid) != 0) && (tigs[ tigs.inUnitig(pf.best3.b_iid) ]->_isUnassembled == true))
        nEdgeToUnasm++;

      if (reportReadGraph_reportEdge(tigs, pf, skipBubble, skipRepeat, reportC, report5, report3) == false)
        continue;

      used[fi] = 1;

      if (reportC)  used[pf.bestC.b_iid] = 1;
      if (report5)  used[pf.best5.b_iid] = 1;
      if (report3)  used[pf.best3.b_iid] = 1;
    }
  }

  writeStatus("AssemblyGraph()-- Found " F_U64 " edges to unassembled contigs.\n", nEdgeToUnasm);

  //  Then write those sequences.

  for (uint32 fi=1; fi<RI->numReads() + 1; fi++)
    if (used[fi] == 1)
      fprintf(BEG, "S\tread%08u\t*\tLN:i:%u\n", fi, RI->readLength(fi));

  delete [] used;


  //  Now, report edges.  GFA wants edges in exactly this format:
  //
  //       -------------
  //             -------------
  //
  //  with read orientation given by +/-.  Conveniently, this is what we've saved (for the edges).

  uint64  nTig[3] = {0,0,0};  //  Number of edges - both contig and unitig
  uint64  nCtg[3] = {0,0,0};  //  Number of edges - contig only
  uint64  nUtg[3] = {0,0,0};  //  Number of edges - unitig only (should be zero)
  uint64  nAsm[3] = {0,0,0};  //  Number of edges - between contigs

  uint64  nBubble = 0;
  uint64  nRepeat = 0;

  for (uint32 fi=1; fi<RI->numReads() + 1; fi++) {
    for (uint32 pp=0; pp<_pForward[fi].size(); pp++) {
      BestPlacement  &pf = _pForward[fi][pp];
      bool            reportC=false, report5=false, report3=false;

      if (reportReadGraph_reportEdge(tigs, pf, skipBubble, skipRepeat, reportC, report5, report3) == false)
        continue;

      //  Some statistics - number of edges of each type (in a contig, in a unitig, in both (tig), in neither (asm))

      if ((pf.isContig == true)  && (pf.isUnitig == true)) {
        if (reportC == true)   nTig[0]++;
        if (report5 == true)   nTig[1]++;
        if (report3 == true)   nTig[2]++;
      }

      if ((pf.isContig == true)  && (pf.isUnitig == false)) {
        if (reportC == true)   nCtg[0]++;
        if (report5 == true)   nCtg[1]++;
        if (report3 == true)   nCtg[2]++;
      }

      if ((pf.isContig == false) && (pf.isUnitig == true)) {
        if (reportC == true)   nUtg[0]++;
        if (report5 == true)   nUtg[1]++;
        if (report3 == true)   nUtg[2]++;
      }

      if ((pf.isContig == false) && (pf.isUnitig == false)) {
        if (reportC == true)   nAsm[0]++;
        if (report5 == true)   nAsm[1]++;
        if (report3 == true)   nAsm[2]++;
      }

      //  Finally, output the edge.

      if (reportC)
        fprintf(BEG, "C\tread%08u\t+\tread%08u\t%c\t%u\t%uM\tic:i:%d\tiu:i:%d\tib:i:%d\tir:i:%d\n",
                fi,
                pf.bestC.b_iid, pf.bestC.flipped ? '-' : '+',
                -pf.bestC.a_hang,
                RI->readLength(fi),
                pf.isContig,
                pf.isUnitig,
                pf.isBubble,
                pf.isRepeat);

      if (report5)
        fprintf(BEG, "L\tread%08u\t-\tread%08u\t%c\t%uM\tic:i:%d\tiu:i:%d\tib:i:%d\tir:i:%d\n",
                fi,
                pf.best5.b_iid, pf.best5.BEndIs3prime() ? '-' : '+',
                RI->overlapLength(fi, pf.best5.b_iid, pf.best5.a_hang, pf.best5.b_hang),
                pf.isContig,
                pf.isUnitig,
                pf.isBubble,
                pf.isRepeat);

      if (report3)
        fprintf(BEG, "L\tread%08u\t+\tread%08u\t%c\t%uM\tic:i:%d\tiu:i:%d\tib:i:%d\tir:i:%d\n",
                fi,
                pf.best3.b_iid, pf.best3.BEndIs3prime() ? '-' : '+',
                RI->overlapLength(fi, pf.best3.b_iid, pf.best3.a_hang, pf.best3.b_hang),
                pf.isContig,
                pf.isUnitig,
                pf.isBubble,
                pf.isRepeat);
    }
  }

  fclose(BEG);

  //  And report statistics.

  writeStatus("AssemblyGraph()-- %8" F_U64P " bubble placements\n", nBubble);
  writeStatus("AssemblyGraph()-- %8" F_U64P " repeat placements\n", nRepeat);
  writeStatus("\n");
  writeStatus("AssemblyGraph()-- Intratig edges:     %8" F_U64P " contained  %8" F_U64P " 5'  %8" F_U64P " 3' (in both contig and unitig)\n", nTig[0], nTig[1], nTig[2]);
  writeStatus("AssemblyGraph()-- Contig only edges:  %8" F_U64P " contained  %8" F_U64P " 5'  %8" F_U64P " 3'\n", nCtg[0], nCtg[1], nCtg[2]);
  writeStatus("AssemblyGraph()-- Unitig only edges:  %8" F_U64P " contained  %8" F_U64P " 5'  %8" F_U64P " 3'\n", nUtg[0], nUtg[1], nUtg[2]);
  writeStatus("AssemblyGraph()-- Intercontig edges:  %8" F_U64P " contained  %8" F_U64P " 5'  %8" F_U64P " 3' (in neither contig nor unitig)\n", nAsm[0], nAsm[1], nAsm[2]);
}

