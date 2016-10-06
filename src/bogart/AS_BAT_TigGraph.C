
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


#undef  SHOW_EDGES
#undef  SHOW_EDGES_FILTERING

//  Grab overlaps for the first read, save any that are to the end of some other tig and of decent quality.
void
saveTigEdges(TigVector &tigs, uint32 fi, vector<BAToverlap> &tigEdges) {
  uint32               no  = 0;
  BAToverlap          *ovl = OC->getOverlaps(fi, no);

  for (uint32 oo=0; oo<no; oo++) {
    if (ovl[oo].filtered == true) {    //  Overlap is garbage.
#ifdef SHOW_EDGES_FILTERING
      writeLog("overlap A: %u tig %u  B: %u tig %u - is garbage.\n",
               ovl[oo].a_iid, tigs.inUnitig(ovl[oo].a_iid),
               ovl[oo].b_iid, tigs.inUnitig(ovl[oo].b_iid));
#endif
      continue;
    }

    if (ovl[oo].AisContainer() ||      //  Overlap can't be a graph edge.
        ovl[oo].AisContained()) {
#ifdef SHOW_EDGES_FILTERING
      writeLog("overlap A: %u tig %u  B: %u tig %u - is containment.\n",
               ovl[oo].a_iid, tigs.inUnitig(ovl[oo].a_iid),
               ovl[oo].b_iid, tigs.inUnitig(ovl[oo].b_iid));
#endif
      continue;
    }

    assert(fi == ovl[oo].a_iid);

    uint32  ai = tigs.inUnitig(ovl[oo].a_iid);
    uint32  bi = tigs.inUnitig(ovl[oo].b_iid);

    if ((ai == 0) || (bi == 0)) {      //  Overlap to something that isn't in a unitig.
#ifdef SHOW_EDGES_FILTERING
      writeLog("overlap A: %u tig %u  B: %u tig %u - is between non-tig things.\n",
               ovl[oo].a_iid, tigs.inUnitig(ovl[oo].a_iid),
               ovl[oo].b_iid, tigs.inUnitig(ovl[oo].b_iid));
#endif
      continue;
    }

    if ((tigs[bi]->firstRead()->ident != ovl[oo].b_iid) &&    //  Overlap isn't to the end of
        (tigs[bi]->lastRead() ->ident != ovl[oo].b_iid)) {    //  some other tig.
#ifdef SHOW_EDGES_FILTERING
      writeLog("overlap A: %u tig %u  B: %u tig %u - is not to an end read.\n",
               ovl[oo].a_iid, tigs.inUnitig(ovl[oo].a_iid),
               ovl[oo].b_iid, tigs.inUnitig(ovl[oo].b_iid));
#endif
      continue;
    }

    //  Success!

    tigEdges.push_back(ovl[oo]);
  }
}



//  Unlike placing bubbles and repeats, we don't have enough coverage to do any
//  fancy filtering based on the error profile.  We thus fall back to using
//  the filtering for best edges.

void
reportTigGraph(TigVector &tigs, const char *prefix, const char *label) {
 char   N[FILENAME_MAX];
  FILE *BEG = NULL;

  bool  skipBubble      = true;
  bool  skipRepeat      = true;
  bool  skipUnassembled = true;

  uint64  nEdgeToUnasm = 0;

  writeStatus("AssemblyGraph()-- generating '%s.unitigs.gfa'.\n", prefix);

  sprintf(N, "%s.unitigs.gfa", prefix);

  BEG = fopen(N, "w");

  if (BEG == NULL)
    return;

  //  Write a header.

  fprintf(BEG, "H\tVN:Z:bogart/edges\n");

  //  Then write the sequences used in the graph.  Unlike the read and contig graphs, every sequence
  //  in our set is output.  By construction, only valid unitigs are in it.

  for (uint32 ti=1; ti<tigs.size(); ti++)
    fprintf(BEG, "S\ttig%08u\t*\tLN:i:%u\n", ti, tigs[ti]->getLength());

  //  A list of the edges to output.  A list of <readID,readID> that we want to output.

  vector<BAToverlap>   tigEdges;

  for (uint32 ti=1; ti<tigs.size(); ti++) {
#ifdef SHOW_EDGES
    writeLog("reportTigGraph()-- tig %u len %u reads %u - firstRead %u lastRead %u\n",
             ti, tigs[ti]->getLength(), tigs[ti]->ufpath.size(), tigs[ti]->firstRead()->ident, tigs[ti]->lastRead()->ident);
#endif

    saveTigEdges(tigs, tigs[ti]->firstRead()->ident, tigEdges);

    if (tigs[ti]->ufpath.size() > 1)
      saveTigEdges(tigs, tigs[ti]->lastRead()->ident, tigEdges);
  }

  //  Now, report edges.  GFA wants edges in exactly this format:
  //
  //       -------------
  //             -------------
  //
  //  with read orientation given by +/-.  Conveniently, this is what we've saved (for the edges).

  for (uint32 te=0; te<tigEdges.size(); te++) {
    uint32          ra      = tigEdges[te].a_iid;
    uint32          rb      = tigEdges[te].b_iid;
    bool            flipped = tigEdges[te].flipped;

    Unitig *tgA     = tigs[ tigs.inUnitig(ra) ];
    Unitig *tgB     = tigs[ tigs.inUnitig(rb) ];

    ufNode *rdA     = &tgA->ufpath[ tigs.ufpathIdx(ra) ];
    ufNode *rdB     = &tgB->ufpath[ tigs.ufpathIdx(rb) ];

    bool    rdAspan = ((rdA->position.min() == 0) && (rdA->position.max() == tgA->getLength())) ? true : false;
    bool    rdBspan = ((rdB->position.min() == 0) && (rdB->position.max() == tgA->getLength())) ? true : false;

    bool    rdAbgn  = (rdA->position.min() == 0) ? true : false;    //  Read at the begin or end of the tig?
    bool    rdBbgn  = (rdB->position.min() == 0) ? true : false;

    bool    rdAfwd  = rdA->position.isForward();                    //  Read forward or flipped in the tig?
    bool    rdBfwd  = rdB->position.isForward();

    bool    rdA3p  = tigEdges[te].AEndIs3prime();                   //  Overlap on the 5' or 3' end of the read?
    bool    rdB3p  = tigEdges[te].BEndIs3prime();

    //  GFA wants edges in normal form:      ---------
    //                                            -----------
    //
    //  Our task here is to decide on an orientation for each read so that layout is formed.
    //  Mostly, it's decided by which end of the tig the read is on.

    bool    tgAfwd  = rdA->position.isForward();  //  Default to the single read case.
    bool    tgBfwd  = rdB->position.isForward();
    bool    invalid = false;

    //  Transfer the overlap to the tig.  We check for invalid cases too.

    //  If the read doesn't span the tig, decode the position and orientation
    //  of the read and use that (and the overlap) to figure out orientation of the tig.

    //  Some of these invalid orientations are from short unitigs with overlaps to
    //  both the first and last read:
    //
    //  tigA--------------->                Whatever read is at the end of this tig
    //      tigB-------------------->       has overlaps to both of these reads,
    //          11111111111111>             but only the overlap to #1 is a valid
    //                22222222222222>       graph edge.

    if (rdAspan == false) {
      if (rdAbgn == true) {
        if ((rdAfwd == true)  && (rdA3p == true))    invalid = true;
        if ((rdAfwd == true)  && (rdA3p == false))   tgAfwd  = false;
        if ((rdAfwd == false) && (rdA3p == true))    tgAfwd  = false;
        if ((rdAfwd == false) && (rdA3p == false))   invalid = true;
      }

      else {
        if ((rdAfwd == true)  && (rdA3p == true))    tgAfwd  = true;
        if ((rdAfwd == true)  && (rdA3p == false))   invalid = true;
        if ((rdAfwd == false) && (rdA3p == true))    invalid = true;
        if ((rdAfwd == false) && (rdA3p == false))   tgAfwd  = true;
      }
    }

    //  Otherwise, the tig is the read, and it's easier (but the same idea).

    else {
      if ((rdAfwd == true)  && (rdA3p == true))    tgAfwd = true;
      if ((rdAfwd == true)  && (rdA3p == false))   tgAfwd = false;
      if ((rdAfwd == false) && (rdA3p == true))    tgAfwd = false;
      if ((rdAfwd == false) && (rdA3p == false))   tgAfwd = true;
    }

    //  Do it all again, but negated, for tig B.

    if (rdBspan == false) {
      if (rdBbgn == true) {
        if ((rdBfwd == true)  && (rdB3p == true))    invalid = true;
        if ((rdBfwd == true)  && (rdB3p == false))   tgBfwd  = true;
        if ((rdBfwd == false) && (rdB3p == true))    tgBfwd  = true;
        if ((rdBfwd == false) && (rdB3p == false))   invalid = true;
      }

      else {
        if ((rdBfwd == true)  && (rdB3p == true))    tgBfwd  = false;
        if ((rdBfwd == true)  && (rdB3p == false))   invalid = true;
        if ((rdBfwd == false) && (rdB3p == true))    invalid = true;
        if ((rdBfwd == false) && (rdB3p == false))   tgBfwd  = false;
      }
    }
    else {
      if ((rdBfwd == true)  && (rdB3p == true))    tgBfwd = false;
      if ((rdBfwd == true)  && (rdB3p == false))   tgBfwd = true;
      if ((rdBfwd == false) && (rdB3p == true))    tgBfwd = true;
      if ((rdBfwd == false) && (rdB3p == false))   tgBfwd = false;
    }

    //  If we assume that the edge is valid - that the resulting contig overlap makes sense - we just
    //  need orient each tig based on which end the read is on.  Everything else should be constrained.

    uint32  olapLen = RI->overlapLength(ra, rb, tigEdges[te].a_hang, tigEdges[te].b_hang);

    if (invalid == false)
      fprintf(BEG, "L\ttig%08u\t%c\ttig%08u\t%c\t%uM\n",
              tgA->id(), (tgAfwd == true) ? '+' : '-',
              tgB->id(), (tgBfwd == true) ? '+' : '-',
              olapLen);

#ifdef SHOW_EDGES
    writeLog("edge read %6u %7u-%-7u (tig %7u len=%7u) -- read %7u %7u-%-7u (tig %7u len=%7u) -- olap rdA3p %d rdB3p %d len %7u -- %s\n",
             ra, rdA->position.bgn, rdA->position.end, tgA->id(), tgA->getLength(),
             rb, rdB->position.bgn, rdB->position.end, tgB->id(), tgB->getLength(),
             rdA3p, rdB3p, olapLen,
             (invalid == true) ? "BAD" : "");
#endif
  }

  fclose(BEG);

  //  And report statistics.

}
