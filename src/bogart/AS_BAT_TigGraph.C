
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



//  Grab overlaps for the first read, save any that are to the end of some other tig and of decent quality.
void
saveTigEdges(TigVector &tigs, uint32 fi, uint32 *used, vector<BAToverlap> &tigEdges) {
  uint32               no  = 0;
  BAToverlap          *ovl = OC->getOverlaps(fi, no);

  for (uint32 oo=0; oo<no; oo++) {
    if (ovl[oo].filtered == true) {    //  Overlap is garbage.
      //writeLog("overlap A: %u tig %u  B: %u tig %u - is garbage.\n",
      //         ovl[oo].a_iid, tigs.inUnitig(ovl[oo].a_iid),
      //         ovl[oo].b_iid, tigs.inUnitig(ovl[oo].b_iid));
      continue;
    }

    if (ovl[oo].AisContainer() ||      //  Overlap can't be a graph edge.
        ovl[oo].AisContained()) {
      //writeLog("overlap A: %u tig %u  B: %u tig %u - is containment.\n",
      //         ovl[oo].a_iid, tigs.inUnitig(ovl[oo].a_iid),
      //         ovl[oo].b_iid, tigs.inUnitig(ovl[oo].b_iid));
      continue;
    }

    assert(fi == ovl[oo].a_iid);

    uint32  ai = tigs.inUnitig(ovl[oo].a_iid);
    uint32  bi = tigs.inUnitig(ovl[oo].b_iid);

    if ((ai == 0) || (bi == 0)) {      //  Overlap to something that isn't in a unitig.
      //writeLog("overlap A: %u tig %u  B: %u tig %u - is between non-tig things.\n",
      //         ovl[oo].a_iid, tigs.inUnitig(ovl[oo].a_iid),
      //         ovl[oo].b_iid, tigs.inUnitig(ovl[oo].b_iid));
      continue;
    }

    if ((tigs[bi]->firstRead()->ident != ovl[oo].b_iid) &&    //  Overlap isn't to the end of
        (tigs[bi]->lastRead() ->ident != ovl[oo].b_iid)) {    //  some other tig.
      //writeLog("overlap A: %u tig %u  B: %u tig %u - is not to an end read.\n",
      //         ovl[oo].a_iid, tigs.inUnitig(ovl[oo].a_iid),
      //         ovl[oo].b_iid, tigs.inUnitig(ovl[oo].b_iid));
      continue;
    }

    //  Success!

    used[ai] = 1;
    used[bi] = 1;

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

  fprintf(BEG, "H\tVN:Z:bogart/edges\n");

  //  A list of the edges to output.  A list of <readID,readID> that we want to output.

  vector<BAToverlap>   tigEdges;

  //  First, figure out what sequences are used.  A tig is used if it has any (valid) edges from
  //  either read on the end.

  uint32   *used = new uint32 [tigs.size() + 1];

  memset(used, 0, sizeof(uint32) * (tigs.size() + 1));

  for (uint32 ti=0; ti<tigs.size(); ti++) {
    if (tigs[ti] == NULL)
      continue;

    writeLog("TIG %u firstRead %u lastRead %u\n", ti, tigs[ti]->firstRead()->ident, tigs[ti]->lastRead()->ident);

    used[ti] = 1;   //  Every tig is reported.  Even those with no edges.

    saveTigEdges(tigs, tigs[ti]->firstRead()->ident, used, tigEdges);
    saveTigEdges(tigs, tigs[ti]->lastRead() ->ident, used, tigEdges);
  }

  //  Then write the sequences used in the graph.

  for (uint32 ti=0; ti<tigs.size(); ti++)
    if (used[ti] == 1)
      fprintf(BEG, "S\ttig%08u\t*\tLN:i:%u\n", ti, tigs[ti]->getLength());



  //  Now, report edges.  GFA wants edges in exactly this format:
  //
  //       -------------
  //             -------------
  //
  //  with read orientation given by +/-.  Conveniently, this is what we've saved (for the edges).

  for (uint32 te=0; te<tigEdges.size(); te++) {
    uint32          ra      = tigEdges[te].a_iid;         //  I miss my buddy.
    uint32          rb      = tigEdges[te].b_iid;         //  So lonely by myself.
    uint32          olapLen = RI->overlapLength(ra, rb, tigEdges[te].a_hang, tigEdges[te].b_hang);
    bool            flipped = tigEdges[te].flipped;

    Unitig *tiga  = tigs[ tigs.inUnitig(ra) ];
    Unitig *tigb  = tigs[ tigs.inUnitig(rb) ];

    ufNode *reada = &tiga->ufpath[ tigs.ufpathIdx(ra) ];
    ufNode *readb = &tigb->ufpath[ tigs.ufpathIdx(rb) ];

    //  Which end of the tig is the read at?
    bool    rdAbgn  = (reada->position.min() == 0) ? true : false;
    bool    rdBbgn  = (readb->position.min() == 0) ? true : false;

    //  Orientation of the read in the tig?
    bool    rdAfwd  = reada->position.isForward();
    bool    rdBfwd  = readb->position.isForward();

    //  Orientation of the overlap?
    bool    rdA3p  = tigEdges[te].AEndIs3prime();
    bool    rdB3p  = tigEdges[te].BEndIs3prime();

    //  If we assume that the edge is valid - that the resulting contig overlap makes sense - we just
    //  need orient each tig based on which end the read is on.  Everything else should be constrained.

    fprintf(BEG, "L\ttig%08u\t%c\ttig%08u\t%c\t%uM\n",     //    ------------>
            tiga->id(), (rdAbgn == false) ? '+' : '-',     //              ---     rdAbgn == false
            tigb->id(), (rdBbgn == true)  ? '+' : '-',     //               <-----------
            olapLen);                                      //               ---    rdBbgn == false

    //  Check for a handful of invalid overlaps.  Assume that if there are invalid layouts generated,
    //  we'll get all possibilities.
    //
    //  If the read is at the same ends of the tig, and the reads are the same oriet, and the olap is not flipped, we conflict.
    //    ------------>    tigA
    //              -->    rdA
    //        --------->   tigB
    //               -->   rdB
    //
    if ((rdAbgn == rdBbgn) && (rdAfwd == rdBfwd) && (flipped == false)) {
      writeLog("bad edge read %u %u-%u tig %u (len=%u) %c -- read %u %u-%u tig %u (len=%u) %c -- len %u\n",
               ra, reada->position.bgn, reada->position.end, tiga->id(), tiga->getLength(), (rdAbgn == false) ? '+' : '-',
               rb, readb->position.bgn, readb->position.end, tigb->id(), tigb->getLength(), (rdBbgn == true)  ? '+' : '-',
               olapLen);
    }
    //assert(!((rdAbgn == rdBbgn) && (rdAfwd == rdBfwd) && (flipped == false)));
  }

  delete [] used;

  fclose(BEG);

  //  And report statistics.

}
