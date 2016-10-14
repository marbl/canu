
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


#define SHOW_EDGES


class  grEdge {
public:
  grEdge() {
    tigID    = 0;
    bgn      = 0;
    end      = 0;
    extended = false;
    deleted  = true;
  };

  grEdge(uint32 t, int32 b, int32 e) {
    tigID    = t;
    bgn      = b;
    end      = e;
    extended = false;
    deleted  = false;
  };

  uint32  tigID;    //  Which tig we're placing this in
  int32   bgn;
  int32   end;

  bool    extended;
  bool    deleted;
};



void
emitEdges(TigVector &tigs,
          Unitig    *tgA,
          bool       isForward,
          FILE      *BEG) {
  vector<overlapPlacement>   placements;
  vector<grEdge>             edges;

  //  Place the first read.

  ufNode   *rdA = tgA->firstRead();

  placeReadUsingOverlaps(tigs, NULL, rdA->ident, placements, placeRead_all);

  //  Convert those placements into potential edges.

  for (uint32 pp=0; pp<placements.size(); pp++) {
    uint32   tid = placements[pp].tigID;
    int32    bgn = placements[pp].verified.min();
    int32    end = placements[pp].verified.max();

    if ((tgA->id() == tid) &&              //  If placed in the same tig and
        (bgn <= rdA->position.max()) &&    //  at the same location, skip it.
        (rdA->position.min() <= end))
      continue;

#ifdef SHOW_EDGES
    writeLog("emitEdges()-- tig %6u read %8u %8u-%-8u -> tig %6u %8u-%-8u\n",
             tgA->id(), 
             rdA->ident, rdA->position.min(), rdA->position.max(),
             tid, bgn, end);
#endif

    edges.push_back(grEdge(tid, bgn, end));
  }

  //  Technically, we should run through the edges and emit those that are already satisfied.
  //  But we can defer this until after the second read is processed.  Well, we could defer
  //  until all reads are processed, but cleaning up the list makes us a little faster, and
  //  also lets us short circuit when we run out of potential edges before we run out of reads
  //  in the tig.

  //  While there are still placements to process, march down the reads in this tig, adding to the
  //  appropriate placement.

  for (uint32 fi=1; (fi<tgA->ufpath.size()) && (edges.size() > 0); fi++) {
    ufNode  *rdA = &tgA->ufpath[fi];

    placeReadUsingOverlaps(tigs, NULL, rdA->ident, placements, placeRead_all);

    //  Mark every edge as being not extended.

    for (uint32 ee=0; ee<edges.size(); ee++)
      edges[ee].extended = false;

    //  Merge the new placements with the saved placements.

    for (uint32 pp=0; pp<placements.size(); pp++) {
      uint32   tid = placements[pp].tigID;
      int32    bgn = placements[pp].verified.min();
      int32    end = placements[pp].verified.max();

      for (uint32 ee=0; ee<edges.size(); ee++) {
        if (edges[ee].deleted == true)     //  Invalid or already finished edge.
          continue;

        if ((tid != edges[ee].tigID) ||    //  Wrong tig, keep looking.
            (end < edges[ee].bgn) ||       //  No intersection, keep looking.
            (edges[ee].end < bgn))
          continue;

        //  Otherwise, the right tig, and we intersect.  Extend the interval and mark it as extended.

        //  We're trusting that we don't find some bizarre repeat that would let us match ABC in
        //  tgA against CAB in the target tig.  If not, we'll need to keep count of which direction
        //  we extend things in.

        edges[ee].bgn      = min(edges[ee].bgn, bgn);
        edges[ee].end      = max(edges[ee].end, end);
        edges[ee].extended = true;
      }
    }

    //  Emit edges that are complete and mark them as done.
    //
    //  A better idea is to see if this read is overlapping with the first/last read
    //  in the other tig, and we're close enough to the end, instead of these silly 100bp thresholds.

    for (uint32 ee=0; ee<edges.size(); ee++) {
      if (edges[ee].bgn <= 100) {
#ifdef SHOW_EDGES_VERBOSE
        writeLog("emitEdges()-- tig %6u %s edgeTo tig %6u %s of length %6u\n",
                 tgA->id(), isForward ? "<--" : "-->",
                 edges[ee].tigID, "-->",
                 edges[ee].end - edges[ee].bgn);
#endif
        fprintf(BEG, "L\ttig%08u\t%c\ttig%08u\t%c\t%uM\n",
                tgA->id(), isForward ? '-' : '+',
                edges[ee].tigID, '+',
                edges[ee].end - edges[ee].bgn);
        edges[ee].deleted = true;
      }

      if (edges[ee].end + 100 >= tigs[edges[ee].tigID]->getLength()) {
#ifdef SHOW_EDGES_VERBOSE
        writeLog("emitEdges()-- tig %6u %s edgeTo tig %6u %s of length %6u\n",
                 tgA->id(), isForward ? "<--" : "-->",
                 edges[ee].tigID, "<--",
                 edges[ee].end - edges[ee].bgn);
#endif
        fprintf(BEG, "L\ttig%08u\t%c\ttig%08u\t%c\t%uM\n",
                tgA->id(), isForward ? '-' : '+',
                edges[ee].tigID, '-',
                edges[ee].end - edges[ee].bgn);
        edges[ee].deleted = true;
      }
    }

    //  A bit of cleverness.  If we emit edges before dealing with deleted and non-extended edges, the first
    //  time we hit this code we'll emit edges for both the first read and the second read.

    for (uint32 ee=0; ee<edges.size(); ee++) {
      if (edges[ee].extended == true)
        continue;

#ifdef SHOW_EDGES
      writeLog("emitEdges()-- tig %6u %s edgeTo tig %6u %s [0 %u-%u %u] UNSATISFIED at read %u #%u\n",
               tgA->id(), isForward ? "<--" : "-->",
               edges[ee].tigID, "-->",
               edges[ee].bgn, edges[ee].end, tigs[edges[ee].tigID]->getLength(),
               rdA->ident, fi);
#endif
                 
      edges[ee].deleted = true;
    }

    //  Compress the edges list (optional) to remove the deleted edges.

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
  for (uint32 ee=0; ee<edges.size(); ee++)
    if (edges[ee].extended == false)
      writeLog("emitEdges()-- tig %6u %s edgeTo tig %6u %s [0 %u-%u %u] UNSATISFIED after all reads\n",
               tgA->id(), isForward ? "<--" : "-->",
               edges[ee].tigID, "-->",
               edges[ee].bgn, edges[ee].end, tigs[edges[ee].tigID]->getLength());
#endif
}



//  Unlike placing bubbles and repeats, we don't have enough coverage to do any
//  fancy filtering based on the error profile.  We thus fall back to using
//  the filtering for best edges.

void
reportTigGraph(TigVector &tigs, const char *prefix, const char *label) {
 char   N[FILENAME_MAX];
  FILE *BEG = NULL;

  writeStatus("AssemblyGraph()-- generating '%s.unitigs.gfa'.\n", prefix);

  sprintf(N, "%s.unitigs.gfa", prefix);

  BEG = fopen(N, "w");

  if (BEG == NULL)
    return;

  //  Write a header.  You've gotta start somewhere!

  fprintf(BEG, "H\tVN:Z:bogart/edges\n");

  //  Then write the sequences used in the graph.  Unlike the read and contig graphs, every sequence
  //  in our set is output.  By construction, only valid unitigs are in it.  Though we occasionally
  //  make a disconnected unitig and need to split it again.

  for (uint32 ti=1; ti<tigs.size(); ti++)
    if (tigs[ti] != NULL)
      fprintf(BEG, "S\ttig%08u\t*\tLN:i:%u\n", ti, tigs[ti]->getLength());

  //  Run through all the tigs, emitting edges for the first and last read.

  for (uint32 ti=1; ti<tigs.size(); ti++) {
    Unitig  *tgA = tigs[ti];

    if (tgA == NULL)
      continue;

#ifdef SHOW_EDGES
    writeLog("\n");
    writeLog("reportTigGraph()-- tig %u len %u reads %u - firstRead %u lastRead %u\n",
             ti, tgA->getLength(), tgA->ufpath.size(), tgA->firstRead()->ident, tgA->lastRead()->ident);
#endif

    emitEdges(tigs, tgA, true,  BEG);

    tgA->reverseComplement();
    emitEdges(tigs, tgA, false, BEG);
    tgA->reverseComplement();
  }

  fclose(BEG);

  //  And report statistics.

}
