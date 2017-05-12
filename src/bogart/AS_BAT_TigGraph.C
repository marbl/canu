
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
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_AssemblyGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_PlaceReadUsingOverlaps.H"

#include "AS_BAT_TigGraph.H"

#undef  SHOW_EDGES
#undef  SHOW_EDGES_UNPLACED   //  Generates a lot of noise
#undef  SHOW_EDGES_VERBOSE


class  grEdge {
public:
  grEdge() {
    tigID    = 0;
    bgn      = 0;
    end      = 0;
    fwd      = false;
    extended = false;
    deleted  = true;
  };

  grEdge(uint32 t, int32 b, int32 e, bool f) {
    tigID    = t;
    bgn      = b;
    end      = e;
    fwd      = f;
    extended = false;
    deleted  = false;
  };

  uint32  tigID;      //  Which tig we're placing this in
  int32   bgn;        //  Location of overlap
  int32   end;        //
  bool    fwd;        //  Overlap indicates tgB is forward (tgA is defined to be forward)

  bool    extended;
  bool    deleted;
};



void
emitEdges(TigVector      &tigs,
          Unitig         *tgA,
          bool            tgAflipped,
          FILE           *BEG,
          vector<tigLoc> &tigSource) {
  vector<overlapPlacement>   placements;
  vector<grEdge>             edges;

  //  Place the first read.

  ufNode   *rdA    = tgA->firstRead();
  uint32    rdAlen = RI->readLength(rdA->ident);

  placeReadUsingOverlaps(tigs, NULL, rdA->ident, placements, placeRead_all);

  //
  //  Somewhere we need to weed out the high error overlaps - Unitig::overlapConsistentWithTig() won't work
  //  because we're at the end of the tig and can have 1x of coverage.
  //

  //  Convert those placements into potential edges.
  //
  //  Overview: from this first placement, we'll try to extend the tig-tig alignment to generate the
  //  full edge.  In pictures:
  //
  //         <-----------------------------------------------  tgA
  //     rd1 ------------>
  //     rd2   -------------->
  //     rd3            <----------------
  //
  //  -------------------------------->    tgB (we don't care about its reads)
  //
  //  We'll place rd1 in tgB, then place rd2 and extend the alignment, then rd3 and notice that
  //  we've covered all of tgB, so an edge is emitted.  If, say, rd2 failed to align fully, we'd
  //  still extend the alignment, and let the total failure of rd3 kill the edge.

  for (uint32 pp=0; pp<placements.size(); pp++) {
    uint32  tgBid  = placements[pp].tigID;
    Unitig *tgB    = tigs[tgBid];
    uint32  tgBlen = tigs[tgBid]->getLength();

    int32    bgn  = placements[pp].verified.min();
    int32    end  = placements[pp].verified.max();

    if ((tgA->id() == tgBid) &&            //  If placed in the same tig and
        (bgn <= rdA->position.max()) &&    //  at the same location, skip it.
        (rdA->position.min() <= end))
      continue;

    if (tgB->_isUnassembled == true)       //  Ignore placements to unassembled crud.
      continue;

    //  For this to be a valid starting edge, the read must be placed from it's beginning.  In the
    //  picture above, rd1 must be placed fully to it's 5' end.  The 3' end can flop around; if the
    //  tig-tig alignment isn't true, then rd2 will fail to align.  Note thhat if the tig-tig
    //  alignment is fully captured by only rd1, its 3' end will flop around, tgB will be fully covered,
    //  and the edge will be emitted.

    if (((rdA->isForward() == true) && (placements[pp].covered.bgn > 0)) ||
        ((rdA->isReverse() == true) && (placements[pp].covered.end < rdAlen))) {
#ifdef SHOW_EDGES
      writeLog("emitEdges()-- edge --- - tig %6u read %8u %8u-%-8u placed bases %8u-%-8u in tig %6u %8u-%-8u - INCOMPLETELY PLACED outside\n",
               tgA->id(),
               rdA->ident, rdA->position.bgn, rdA->position.end,
               placements[pp].covered.bgn, placements[pp].covered.end,
               tgBid, bgn, end);
#endif
      continue;
    }

    //  Now, if the placed read didn't get placed to it's other end, and it's placed in the middle
    //  of the tig, reject the placement.

    if (((rdA->isForward() == true) && (placements[pp].covered.end < rdAlen) && (bgn > 100) && (end + 100 < tgBlen)) ||
        ((rdA->isReverse() == true) && (placements[pp].covered.bgn > 0)      && (bgn > 100) && (end + 100 < tgBlen))) {
#ifdef SHOW_EDGES
      writeLog("emitEdges()-- edge --- - tig %6u read %8u %8u-%-8u placed bases %8u-%-8u in tig %6u %8u-%-8u - INCOMPLETELY PLACED inside\n",
               tgA->id(),
               rdA->ident, rdA->position.bgn, rdA->position.end,
               placements[pp].covered.bgn, placements[pp].covered.end,
               tgBid, bgn, end, tgBlen);
#endif
      continue;
    }

#ifdef SHOW_EDGES
    writeLog("emitEdges()-- edge %3u - tig %6u read %8u %8u-%-8u placed bases %8u-%-8u in tig %6u %8u-%-8u quality %f\n",
             edges.size(),
             tgA->id(),
             rdA->ident, rdA->position.bgn, rdA->position.end,
             placements[pp].covered.bgn, placements[pp].covered.end,
             tgBid, bgn, end,
             (double)placements[pp].errors / placements[pp].aligned);
#endif

    //  Decide the orientation of the second tig based on the orientation of the read and its
    //  alignment.  If the orientations are the same, then the second tig doesn't need to be
    //  flipped.
    //
    //        <-------------------------------------
    //        <---  read in first tig
    //
    //        <---  alignment on second tig  -  so if not the same, the second tig needs to be
    //    ------------------->               -  flipped to make the alignment work

    bool fwd = false;

    if (((rdA->isForward() == true)  && (placements[pp].verified.isForward() == true)) ||
        ((rdA->isForward() == false) && (placements[pp].verified.isForward() == false)))
      fwd = true;

    edges.push_back(grEdge(tgBid, bgn, end, fwd));
  }

  //  Technically, we should run through the edges and emit those that are already satisfied.  But
  //  we can defer this until after the second read is processed.  Heck, we could defer until all
  //  reads are processed, but cleaning up the list makes us a little faster, and also lets us short
  //  circuit when we run out of potential edges before we run out of reads in the tig.

  //  While there are still placements to process, march down the reads in this tig, adding to the
  //  appropriate placement.

  for (uint32 fi=1; (fi<tgA->ufpath.size()) && (edges.size() > 0); fi++) {
    ufNode  *rdA    = &tgA->ufpath[fi];
    uint32   rdAlen = RI->readLength(rdA->ident);

    placeReadUsingOverlaps(tigs, NULL, rdA->ident, placements, placeRead_all);

    //  Mark every edge as being not extended.

    for (uint32 ee=0; ee<edges.size(); ee++)
      edges[ee].extended = false;

    //  Merge the new placements with the saved placements.

    for (uint32 pp=0; pp<placements.size(); pp++) {
      uint32   tgBid  = placements[pp].tigID;
      Unitig  *tgB    = tigs[tgBid];
      uint32   tgBlen = tigs[tgBid]->getLength();
      int32    bgn    = placements[pp].verified.min();
      int32    end    = placements[pp].verified.max();

      //  Ignore placements to unassembled crud.  Just an optimization.  We'd filter these out
      //  when trying to associate it with an existing overlap.

      if (tgB->_isUnassembled == true)
        continue;

      //  Accept the placement only if it is for the whole read, or if it is touching the end of the target tig.

      if (((placements[pp].covered.bgn > 0) ||
           (placements[pp].covered.end < rdAlen)) &&
          (bgn       > 100) &&
          (end + 100 < tgBlen)) {
#ifdef SHOW_EDGES_UNPLACED
        writeLog("emitEdges()-- read %5u incomplete placement covering %5u-%-5u in at %5u-%-5u in tig %4u\n",
                 rdA->ident, placements[pp].covered.bgn, placements[pp].covered.end, bgn, end, tgBid);
#endif
        continue;
      }

      for (uint32 ee=0; ee<edges.size(); ee++) {
        if (edges[ee].deleted == true)     //  Invalid or already finished edge.
          continue;

        if ((tgBid != edges[ee].tigID) ||    //  Wrong tig, keep looking.
            (end < edges[ee].bgn) ||       //  No intersection, keep looking.
            (edges[ee].end < bgn))
          continue;

        //  Otherwise, the right tig, and we intersect.  Extend the interval and mark it as extended.

        //  We're trusting that we don't find some bizarre repeat that would let us match ABC in
        //  tgA against CAB in the target tig.  If not, we'll need to keep count of which direction
        //  we extend things in.


        //  Fail if most of the extension is to the wrong side.  We always move to higher
        //  coordinates on tgA.  If tgB is forward, it should move to higher coordinates too.

        int32  nbgn = min(edges[ee].bgn, bgn);
        int32  nend = max(edges[ee].end, end);

        if ((edges[ee].fwd == true) &&
            (bgn - nbgn > nend - end)) {  //  If we decrease bgn more than we increased end, fail
#ifdef SHOW_EDGES_UNPLACED
        writeLog("emitEdges()-- edge %3u - extend from %5u-%-5u to %5u-%-5u -- placed read %5u at %5u-%-5u in tig %4u - wrong direction\n",
                 ee,
                 edges[ee].bgn, edges[ee].end,
                 nbgn, nend,
                 rdA->ident, bgn, end, tgBid);
#endif
          continue;
        }

        //  The reverse case is a bit tricky since we're tracking min/max posiiton on tgB.
        //  When we extend on tgA, we expect the bgn to decrease on tgB and the end to stay the same.

        if ((edges[ee].fwd == false) &&
            (nend - end > bgn - nbgn)) {  //  If we increase end more than we decreased bgn, fail
#ifdef SHOW_EDGES_UNPLACED
          writeLog("emitEdges()-- edge %3u - extend from %5u-%-5u to %5u-%-5u -- placed read %5u at %5u-%-5u in tig %4u - wrong direction\n",
                   ee,
                   edges[ee].bgn, edges[ee].end,
                   nbgn, nend,
                   rdA->ident, bgn, end, tgBid);
#endif
          continue;
        }

#ifdef SHOW_EDGES
        writeLog("emitEdges()-- edge %3u - extend from %5u-%-5u to %5u-%-5u -- placed read %5u at %5u-%-5u in tig %4u\n",
                 ee,
                 edges[ee].bgn, edges[ee].end,
                 nbgn, nend,
                 rdA->ident, bgn, end, tgBid);
#endif

        edges[ee].bgn      = nbgn;
        edges[ee].end      = nend;
        edges[ee].extended = true;
      }
    }

    //  Emit edges that are complete and mark them as done.
    //
    //  A better idea is to see if this read is overlapping with the first/last read
    //  in the other tig, and we're close enough to the end, instead of these silly 100bp thresholds.

    //  For edges making circles, when tgA == tgB, we need to flip tgB if tgA is flipped.

    for (uint32 ee=0; ee<edges.size(); ee++) {
      bool  tgBflipped = (edges[ee].tigID == tgA->id()) && (tgAflipped);

      bool  sameContig = false;

      if ((tigSource.size() > 0) && (tigSource[tgA->id()].cID == tigSource[edges[ee].tigID].cID))
        sameContig = true;

      if ((edges[ee].fwd == false) && (edges[ee].bgn <= 100)) {
#ifdef SHOW_EDGES_VERBOSE
        writeLog("emitEdges()-- edge %3u - tig %6u %s edgeTo tig %6u %s of length %6u (%6u-%6u)\n",
                 ee,
                 tgA->id(),       tgAflipped ? "<--" : "-->",
                 edges[ee].tigID, tgBflipped ? "-->" : "<--",
                 edges[ee].end - edges[ee].bgn, edges[ee].bgn, edges[ee].end);
#endif
        fprintf(BEG, "L\ttig%08u\t%c\ttig%08u\t%c\t%uM%s\n",
                edges[ee].tigID, tgBflipped ? '+' : '-',
                tgA->id(),       tgAflipped ? '-' : '+',
                edges[ee].end - edges[ee].bgn,
                (sameContig == true) ? "\tcv:A:T" : "\tcv:A:F");
        edges[ee].deleted = true;
      }

      if ((edges[ee].fwd == true) && (edges[ee].end + 100 >= tigs[edges[ee].tigID]->getLength())) {
#ifdef SHOW_EDGES_VERBOSE
        writeLog("emitEdges()-- edge %3u - tig %6u %s edgeTo tig %6u %s of length %6u (%6u-%6u)\n",
                 ee,
                 tgA->id(),       tgAflipped ? "<--" : "-->",
                 edges[ee].tigID, tgBflipped ? "<--" : "-->",
                 edges[ee].end - edges[ee].bgn, edges[ee].bgn, edges[ee].end);
#endif
        fprintf(BEG, "L\ttig%08u\t%c\ttig%08u\t%c\t%uM%s\n",
                edges[ee].tigID, tgBflipped ? '-' : '+',
                tgA->id(),       tgAflipped ? '-' : '+',
                edges[ee].end - edges[ee].bgn,
                (sameContig == true) ? "\tcv:A:T" : "\tcv:A:F");
        edges[ee].deleted = true;
      }
    }

    //  A bit of cleverness.  If we emit edges before dealing with deleted and non-extended edges, the first
    //  time we hit this code we'll emit edges for both the first read and the second read.

    for (uint32 ee=0; ee<edges.size(); ee++) {
      bool  tgBflipped = (edges[ee].tigID == tgA->id()) && (tgAflipped);

      if (edges[ee].fwd == false)
        tgBflipped = !tgBflipped;

      if (edges[ee].deleted == true)
        continue;

      if (edges[ee].extended == true)
        continue;

#ifdef SHOW_EDGES
      writeLog("emitEdges()-- tig %6u %s edgeTo tig %6u %s [0 %u-%u %u] UNSATISFIED at read %u #%u\n",
               tgA->id(),       tgAflipped ? "<--" : "-->",
               edges[ee].tigID, tgBflipped ? "<--" : "-->",
               edges[ee].bgn, edges[ee].end, tigs[edges[ee].tigID]->getLength(),
               rdA->ident, fi);
#endif

      edges[ee].deleted = true;
    }

    //  Compress the edges list (optional, but messes up logging if not done) to remove the deleted edges.

    uint32 oo = 0;

    for (uint32 ee=0; ee<edges.size(); ee++) {
      if (edges[ee].deleted == false) {   //  Not deleted, so copy it to the output vector
        if (ee != oo)                     //  at location oo.
          edges[oo] = edges[ee];
        oo++;
      }
    }

    edges.resize(oo);   //  Reset the vector to size we ended up with.

    //  And now place the next read in the source tig.
  }

  //  Any edges still on the list aren't edges, so we're all done without needing to check anything.

#ifdef SHOW_EDGES
  for (uint32 ee=0; ee<edges.size(); ee++) {
    bool  tgBflipped = (edges[ee].tigID == tgA->id()) && (tgBflipped);

    if (edges[ee].fwd == false)
      tgBflipped = !tgBflipped;

    if (edges[ee].extended == false)
      writeLog("emitEdges()-- tig %6u %s edgeTo tig %6u %s [0 %u-%u %u] UNSATISFIED after all reads\n",
               tgA->id(),       tgAflipped ? "<--" : "-->",
               edges[ee].tigID, tgBflipped ? "<--" : "-->",
               edges[ee].bgn, edges[ee].end, tigs[edges[ee].tigID]->getLength());
  }
#endif
}



//  Unlike placing bubbles and repeats, we don't have enough coverage to do any
//  fancy filtering based on the error profile.  We thus fall back to using
//  the filtering for best edges.

void
reportTigGraph(TigVector &tigs,
               vector<tigLoc> &tigSource,
               const char *prefix,
               const char *label) {
  char   BEGn[FILENAME_MAX];
  char   BEDn[FILENAME_MAX];

  writeLog("\n");
  writeLog("----------------------------------------\n");
  writeLog("Generating graph\n");

  writeStatus("AssemblyGraph()-- generating '%s.%s.gfa'.\n", prefix, label);

  snprintf(BEGn, FILENAME_MAX, "%s.%s.gfa", prefix, label);
  snprintf(BEDn, FILENAME_MAX, "%s.%s.bed", prefix, label);

  FILE *BEG =                          fopen(BEGn, "w");
  FILE *BED = (tigSource.size() > 0) ? fopen(BEDn, "w") : NULL;

  if (BEG == NULL)
    return;

  //  Write a header.  You've gotta start somewhere!

  fprintf(BEG, "H\tVN:Z:bogart/edges\n");

  //  Then write the sequences used in the graph.  Unlike the read and contig graphs, every sequence
  //  in our set is output.  By construction, only valid unitigs are in it.  Though we occasionally
  //  make a disconnected unitig and need to split it again.

  for (uint32 ti=1; ti<tigs.size(); ti++)
    if ((tigs[ti] != NULL) && (tigs[ti]->_isUnassembled == false))
      fprintf(BEG, "S\ttig%08u\t*\tLN:i:%u\n", ti, tigs[ti]->getLength());

  //  Run through all the tigs, emitting edges for the first and last read.

  for (uint32 ti=1; ti<tigs.size(); ti++) {
    Unitig  *tgA = tigs[ti];

    if ((tgA == NULL) || (tgA->_isUnassembled == true))
      continue;

    //if (ti == 4)
    //  logFileFlags |= LOG_PLACE_READ;

#ifdef SHOW_EDGES
    writeLog("\n");
    writeLog("reportTigGraph()-- tig %u len %u reads %u - firstRead %u\n",
             ti, tgA->getLength(), tgA->ufpath.size(), tgA->firstRead()->ident);
#endif

    emitEdges(tigs, tgA, false, BEG, tigSource);

#ifdef SHOW_EDGES
    writeLog("\n");
    writeLog("reportTigGraph()-- tig %u len %u reads %u - lastRead %u\n",
             ti, tgA->getLength(), tgA->ufpath.size(), tgA->lastRead()->ident);
#endif

    tgA->reverseComplement();
    emitEdges(tigs, tgA, true, BEG, tigSource);
    tgA->reverseComplement();

    if ((tigSource.size() > 0) && (tigSource[ti].cID != UINT32_MAX))
      fprintf(BED, "ctg%08u\t%u\t%u\tutg%08u\t%u\t%c\n",
              tigSource[ti].cID,
              tigSource[ti].cBgn,
              tigSource[ti].cEnd,
              ti,
              0,
              '+');

    //logFileFlags &= ~LOG_PLACE_READ;
  }

  if (BEG)   fclose(BEG);
  if (BED)   fclose(BED);

  //  And report statistics.

}
