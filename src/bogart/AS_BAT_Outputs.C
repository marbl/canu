
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
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_Outputs.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2014-MAR-31
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUN-05
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Instrumentation.H"

#include "AS_BAT_Outputs.H"

//#include "AS_CGB_histo.H"
//#include "tgStore.H"




//  Massage the Unitig into a MultiAlignT (also used in SplitChunks_CGW.c)
void
unitigToTig(tgTig       *tig,
            uint32       tigid,
            Unitig      *utg) {

  tig->clear();

  tig->_tigID           = tigid;
  utg->_tigID           = tigid;

  tig->_coverageStat    = 1.0;  //  Default to just barely unique
  tig->_microhetProb    = 1.0;  //  Default to 100% probability of unique

  if      (utg->_isUnassembled == true)
    tig->_class = tgTig_unassembled;

  else if (utg->_isBubble == true)
    tig->_class = tgTig_bubble;

  else
    tig->_class = tgTig_contig;

  tig->_suggestRepeat   = (utg->_isRepeat   == true);
  tig->_suggestCircular = (utg->_isCircular == true);

  tig->_layoutLen       = utg->getLength();

  resizeArray(tig->_children, tig->_childrenLen, tig->_childrenMax, utg->ufpath.size(), resizeArray_doNothing);

  map<uint32,bool>  forward;
  map<uint32,bool>  allreads;

  //  Just for stats, build a map fo the reads in the unitig.

  for (uint32 ti=0; ti<utg->ufpath.size(); ti++)
    allreads[utg->ufpath[ti].ident] = true;

  //  Process all reads.

  for (uint32 ti=0; ti<utg->ufpath.size(); ti++) {
    ufNode        *frg   = &utg->ufpath[ti];

    //  Remember that we've placed this read, and if it was forward or reverse.
    forward[frg->ident] = (frg->position.bgn < frg->position.end);

    //  If the first read, just dump it in the unitig with no parent.
    if (ti == 0) {
      tig->addChild()->set(frg->ident, 0, 0, 0, frg->position.bgn, frg->position.end);
      continue;
    }

    //  Otherwise, find the thickest overlap to any read already placed in the unitig.

    uint32         olapsLen = 0;
    BAToverlap    *olaps = OC->getOverlaps(frg->ident, AS_MAX_EVALUE, olapsLen);

    uint32         tt     = UINT32_MAX;
    uint32         ttLen  = 0;
    double         ttErr  = DBL_MAX;

    int32          ah     = 0;
    int32          bh     = 0;

    uint32         notPresent = 0;  //  Potential parent isn't in the unitig
    uint32         notPlaced  = 0;  //  Potential parent isn't placed yet
    uint32         negHang    = 0;  //  Potential parent has a negative hang to a placed read
    uint32         goodOlap   = 0;

    for (uint32 oo=0; oo<olapsLen; oo++) {

      if (allreads.count(olaps[oo].b_iid) == 0) {
        notPresent++;
        continue;
      }

      if (forward.count(olaps[oo].b_iid) == 0) {       //  Potential parent not placed yet
        notPlaced++;
        continue;
      }

      uint32  l = FI->overlapLength(olaps[oo].a_iid, olaps[oo].b_iid, olaps[oo].a_hang, olaps[oo].b_hang);

      //  Compute the hangs, so we can ignore those that would place this read before the parent.
      //  This is a flaw somewhere in bogart, and should be caught and fixed earlier.

      //  Consensus is expecting the have the hangs for the parent read, not this read, and some
      //  fiddling is needed to flip the overlap for this:
      //    First, swap the reads so it's b-vs-a.
      //    Then, flip the overlap if the b read is in the unitig flipped.

      int32 ah = (olaps[oo].flipped == false) ? (-olaps[oo].a_hang) : (olaps[oo].b_hang);
      int32 bh = (olaps[oo].flipped == false) ? (-olaps[oo].b_hang) : (olaps[oo].a_hang);

      if (forward[olaps[oo].b_iid] == false) {
        swap(ah, bh);
        ah = -ah;
        bh = -bh;
      }

      //  If the ahang is negative, we flubbed up somewhere, and want to place this read before
      //  the parent (even though positions say to place it after, because we sorted by position).

      if (ah < 0) {
        //fprintf(stderr, "ERROR: read %u in tig %u has negative ahang from parent read %u, ejected.\n",
        //        frg->ident, ti, olaps[oo].b_iid);
        negHang++;
        continue;
      }

      //  The overlap is good.  Count it as such.

      goodOlap++;

      //  If the overlap is worse than the one we already have, we don't care.

      if ((l < ttLen) ||                  //  Too short
          (ttErr < olaps[oo].erate)) {    //  Too noisy
        continue;
      }

      tt    = oo;
      ttLen = l;
      ttErr = olaps[oo].erate;
    }

    //  If no thickest overlap, we screwed up somewhere.  Complain and eject the read.

    if (tt == UINT32_MAX) {
      fprintf(stderr, "ERROR: read %u in tig %u has no overlap to any previous read, ejected.  %u overlaps total.  %u negative hang.  %u to read not in tig.  %u to read later in tig.  %u good overlaps.\n",
              frg->ident, tig->tigID(), olapsLen, negHang, notPresent, notPlaced, goodOlap);
      continue;
    }

    tig->addChild()->set(frg->ident, olaps[tt].b_iid, ah, bh, frg->position.bgn, frg->position.end);
  }

  //fprintf(stderr, "unitigToTig()--  tig %u has %u children\n", tig->_tigID, tig->_childrenLen);
}



void
writeUnitigsToStore(UnitigVector  &unitigs,
                    char          *fileprefix,
                    char          *tigStorePath,
                    uint32         frg_count_target,
                    bool           isFinal) {
  uint32      utg_count              = 0;
  uint32      frg_count              = 0;
  uint32      prt_count              = 1;
  char        filename[FILENAME_MAX] = {0};
  uint32     *partmap                = new uint32 [unitigs.size()];

  //  This code closely follows that in AS_CGB_unitigger.c::output_the_chunks()

  if (isFinal)
    checkUnitigMembership(unitigs);

  // Open up the initial output file

  sprintf(filename, "%s.iidmap", fileprefix);
  FILE *iidm = fopen(filename, "w");
  assert(NULL != iidm);

  sprintf(filename, "%s.partitioning", fileprefix);
  FILE *part = fopen(filename, "w");
  assert(NULL != part);

  sprintf(filename, "%s.partitioningInfo", fileprefix);
  FILE *pari = fopen(filename, "w");
  assert(NULL != pari);

  //  Step through all the unitigs once to build the partition mapping and IID mapping.

  tgStore     *tigStore = new tgStore(tigStorePath);
  tgTig       *tig      = new tgTig;

  for (uint32 tigID=0, ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if ((utg == NULL) || (utg->getNumFrags() == 0))
      continue;

    assert(utg->getLength() > 0);

    //  Convert the bogart tig to a tgTig and save to the store.

    unitigToTig(tig, (isFinal) ? tigID : ti, utg);
    tigID++;

    tigStore->insertTig(tig, false);

    //  Increment the partition if the current one is too large.

    if ((frg_count + utg->getNumFrags() >= frg_count_target) &&
        (frg_count                      >  0)) {
      fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",
              prt_count, utg_count, frg_count);

      prt_count++;
      utg_count = 0;
      frg_count = 0;
    }

    //  Note that the tig is included in this partition.

    utg_count += 1;
    frg_count += utg->getNumFrags();

    //  Map the tig to a partition, and log both the tig-to-partition map and the partition-to-read map.

    fprintf(iidm, "bogart "F_U32" -> tig "F_U32" (in partition "F_U32" with "F_U32" frags)\n",
            utg->id(),
            utg->tigID(),
            prt_count,
            utg->getNumFrags());

    for (uint32 fragIdx=0; fragIdx<utg->getNumFrags(); fragIdx++)
      fprintf(part, "%d\t%d\n", prt_count, utg->ufpath[fragIdx].ident);
  }

  fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",   //  Don't forget to log the last partition!
          prt_count, utg_count, frg_count);

  fclose(pari);
  fclose(part);
  fclose(iidm);

  delete    tig;
  delete    tigStore;
}


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
//
//  Wow, this is ancient.
//
void
writeOverlapsUsed(UnitigVector &unitigs,
                  char         *prefix) {
  char   N[FILENAME_MAX];

  sprintf(N, "%s.unused.best.edges", prefix);

  FILE  *F = fopen(N, "w");

  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];
    Unitig  *ovl = NULL;
    char     tyt = 'C';

    if (tig == NULL)
      continue;

    if (tig->_isUnassembled)  tyt = 'U';
    if (tig->_isBubble)       tyt = 'B';
    if (tig->_isRepeat)       tyt = 'R';
    if (tig->_isCircular)     tyt = 'O';

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  *frg = &tig->ufpath[fi];
      ufNode  *oth = NULL;

      //  Report the unused best edge

      BestEdgeOverlap *be5 = OG->getBestEdgeOverlap(frg->ident, false);
      uint32   rd5 = (be5 == NULL) ?    0 : be5->fragId();
      Unitig  *tg5 = (be5 == NULL) ? NULL : unitigs[Unitig::fragIn(rd5)];
      char     ty5 = 'C';

      if ((tg5 != NULL) && (tg5->tigID() != tig->tigID())) {
        uint32  ord = Unitig::pathPosition(rd5);
        ufNode *oth = &tg5->ufpath[ord];

        if (tig->_isUnassembled)  ty5 = 'U';
        if (tig->_isBubble)       ty5 = 'B';
        if (tig->_isRepeat)       ty5 = 'R';
        if (tig->_isCircular)     ty5 = 'O';

        fprintf(F, "tig %7u %c read %8u at %9u %-9u %c' -- %8d %-8d -- tig %7u %c read %8u at %9u %-9u %c'\n",
                tig->tigID(), tyt, frg->ident, frg->position.bgn, frg->position.end, '5',
                be5->ahang(), be5->bhang(),
                tg5->tigID(), ty5, oth->ident, oth->position.bgn, oth->position.end, (be5->frag3p() == false) ? '5' : '3');
      }

      BestEdgeOverlap *be3 = OG->getBestEdgeOverlap(frg->ident, true);
      uint32   rd3 = (be3 == NULL) ?    0 : be3->fragId();
      Unitig  *tg3 = (be3 == NULL) ? NULL : unitigs[Unitig::fragIn(rd3)];
      char     ty3 = 'C';

      if ((tg3 != NULL) && (tg3->tigID() != tig->tigID())) {
        uint32  ord = Unitig::pathPosition(rd3);
        ufNode *oth = &tg3->ufpath[ord];

        if (tig->_isUnassembled)  ty3 = 'U';
        if (tig->_isBubble)       ty3 = 'B';
        if (tig->_isRepeat)       ty3 = 'R';
        if (tig->_isCircular)     ty3 = 'O';

        fprintf(F, "tig %7u %c read %8u at %9u %-9u %c' -- %8d %-8d -- tig %7u %c read %8u at %9u %-9u %c'\n",
                tig->tigID(), tyt, frg->ident, frg->position.bgn, frg->position.end, '3',
                be3->ahang(), be3->bhang(),
                tg3->tigID(), ty3, oth->ident, oth->position.bgn, oth->position.end, (be3->frag3p() == false) ? '5' : '3');
      }
    }
  }

  fclose(F);
}
