
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
 *    src/AS_CNS/MultiAlignUnitig.C
 *    src/AS_CNS/MultiAlignUnitig.c
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/MultiAlignUnitig.C
 *
 *  Modifications by:
 *
 *    Michael Schatz on 2004-SEP-23
 *      are Copyright 2004 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2005-MAR-22
 *      are Copyright 2005 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2005-MAR-30 to 2008-FEB-13
 *      are Copyright 2005-2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Gennady Denisov from 2005-MAY-09 to 2008-JUN-06
 *      are Copyright 2005-2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUN-16 to 2013-OCT-04
 *      are Copyright 2005-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2005-SEP-29 to 2006-OCT-03
 *      are Copyright 2005-2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-FEB-27 to 2009-JUN-05
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2011-OCT-27
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-AUG-11
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-17
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "unitigConsensus.H"

// for pbdagcon
#include "Alignment.H"
#include "AlnGraphBoost.H"
#include "SimpleAligner.H"

#include "NDalign.H"

#include <set>

using namespace std;



//  Define this.  Use the faster aligner from overlapper.  If not defined,
//  a full O(n^2) DP is computed.
//
#undef  WITH_NDALIGN
#define WITH_NDALIGN



unitigConsensus::unitigConsensus(gkStore  *gkpStore_,
                                 double    errorRate_,
                                 double    errorRateMax_,
                                 uint32    minOverlap_) {

  gkpStore        = gkpStore_;

  tig             = NULL;
  numfrags        = 0;
  trace           = NULL;
  abacus          = NULL;
  utgpos          = NULL;
  cnspos          = NULL;
  tiid            = 0;
  piid            = -1;

  minOverlap      = minOverlap_;
  errorRate       = errorRate_;
  errorRateMax    = errorRateMax_;

  oaPartial       = NULL;
  oaFull          = NULL;
}


unitigConsensus::~unitigConsensus() {
  delete [] trace;
  delete    abacus;

  delete [] utgpos;
  delete [] cnspos;

  delete    oaPartial;
  delete    oaFull;
}




void
unitigConsensus::reportStartingWork(void) {
  if (showProgress())
    fprintf(stderr, "unitigConsensus()-- processing read %u/%u id %d pos %d,%d anchor %d,%d,%d -- length %u\n",
            tiid+1, numfrags,
            utgpos[tiid].ident(),
            utgpos[tiid].min(),
            utgpos[tiid].max(),
            utgpos[tiid].anchor(),
            utgpos[tiid].aHang(),
            utgpos[tiid].bHang(),
            abacus->numberOfColumns());

  if (showPlacementBefore())
    for (int32 x=0; x<=tiid; x++)
      fprintf(stderr, "unitigConsensus()-- mid %10d  utgpos %7d,%7d  cnspos %7d,%7d  anchor %10d,%6d,%6d\n",
              utgpos[x].ident(),
              utgpos[x].min(), utgpos[x].max(),
              cnspos[x].min(), cnspos[x].max(),
              utgpos[x].anchor(), utgpos[x].aHang(), utgpos[x].bHang());
}


void
unitigConsensus::reportFailure(void) {
  fprintf(stderr, "unitigConsensus()-- failed to align fragment %d in unitig %d.\n",
          utgpos[tiid].ident(), tig->tigID());
}


void
unitigConsensus::reportSuccess() {
  //fprintf(stderr, "unitigConsensus()-- fragment %d aligned in unitig %d.\n",
  //        utgpos[tiid].ident(), tig->tigID());
}



//  Dump the unitig and reads to a single file.  We should also probably save any parameters,
//  but then it's not clear what to do with them.
//
bool
unitigConsensus::savePackage(FILE   *outPackageFile,
                             tgTig  *tig) {

  //  Saving the tig is easy, just use the standard dump.

  tig->saveToStream(outPackageFile);

  //  Saving the reads is also easy, but it's a non-standard dump.

  for (uint32 ii=0; ii<tig->numberOfChildren(); ii++)
    gkpStore->gkStore_saveReadToStream(outPackageFile, tig->getChild(ii)->ident());

  return(true);
}



bool
unitigConsensus::generate(tgTig                     *tig_,
                          map<uint32, gkRead *>     *inPackageRead_,
                          map<uint32, gkReadData *> *inPackageReadData_) {

  tig      = tig_;
  numfrags = tig->numberOfChildren();

  if (initialize(inPackageRead_, inPackageReadData_) == FALSE) {
    fprintf(stderr, "generate()--  Failed to initialize for tig %u with %u children\n", tig->tigID(), tig->numberOfChildren());
    goto returnFailure;
  }

  while (moreFragments()) {
    reportStartingWork();

    //  First attempt, all default parameters

    if (computePositionFromAnchor()    && alignFragment())  goto applyAlignment;
    if (computePositionFromLayout()    && alignFragment())  goto applyAlignment;
    if (computePositionFromAlignment() && alignFragment())  goto applyAlignment;

    //  Second attempt, default parameters after recomputing consensus sequence.

    if (showAlgorithm())
      fprintf(stderr, "generateMultiAlignment()-- recompute full consensus\n");

    recomputeConsensus(showMultiAlignments());

    if (computePositionFromAnchor()    && alignFragment())  goto applyAlignment;
    if (computePositionFromLayout()    && alignFragment())  goto applyAlignment;
    if (computePositionFromAlignment() && alignFragment())  goto applyAlignment;

    //  Third attempot, use whatever aligns.  (alignFragment(true) forced it to align, but that's breaking the consensus with garbage alignments)

    if (computePositionFromAlignment() && alignFragment(true))  goto applyAlignment;

    //  Nope, failed to align.

    reportFailure();
    continue;

  applyAlignment:
    setErrorRate(errorRate);
    setMinOverlap(minOverlap);

    reportSuccess();

    abacus->applyAlignment(tiid, traceABgn, traceBBgn, trace, traceLen);

    refreshPositions();
  }

  generateConsensus(tig);

  return(true);

 returnFailure:
  fprintf(stderr, "generateMultiAlignment()-- unitig %d FAILED.\n", tig->tigID());

  //  tgTig should have no changes.

  return(false);
}


bool
unitigConsensus::generatePBDAG(tgTig                     *tig_,
                               map<uint32, gkRead *>     *inPackageRead_,
                               map<uint32, gkReadData *> *inPackageReadData_) {
    tig      = tig_;
    numfrags = tig->numberOfChildren();

    if (initialize(inPackageRead_, inPackageReadData_) == FALSE) {
      fprintf(stderr, "generatePBDAG()--  Failed to initialize for tig %u with %u children\n", tig->tigID(), tig->numberOfChildren());
      return(false);
    }

    // first we need to load into Unitig data structure the quick cns
    Unitig utg;
    utg.id = tig->tigID();

    utg.seq = string(tig->_layoutLen, 'N');

    // build a quick consensus to align to, just smash together sequences.
    for (int i = 0; i < numfrags; i++) {
        gkRead  *read    = gkpStore->gkStore_getRead(utgpos[i].ident());
        uint32   readLen = read->gkRead_sequenceLength();

        uint32 start = utgpos[i].min();
        uint32 end = utgpos[i].max();

        if (start > utg.seq.length()) {
          start = utg.seq.length() - 1;
        }
        if (end - start > readLen) {
           end = start + readLen;
        }
        if (end > utg.seq.length()) {
           end = utg.seq.length() - 1;
        }

        abSequence  *seq      = abacus->getSequence(i);
        char        *fragment = seq->getBases();

        for (int j = start; j < end; j++) {
           if (utg.seq[j] == 'N') {
              utg.seq[j] = fragment[j - start];
           }
        }
    }
    AlnGraphBoost ag(utg.seq);

    // compute alignments of each sequence in parallel
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numfrags; i++) {
        bool placed = computePositionFromLayout();
        dagcon::Alignment aln;
        SimpleAligner align;

        // for each fragment align it
        abSequence  *seq      = abacus->getSequence(i);
        char        *fragment = seq->getBases();

        aln.start = utgpos[i].min();
        aln.end = utgpos[i].max();
        aln.frgid = utgpos[i].ident();
        aln.qstr = string(fragment);
        aln.tstr = utg.seq.substr(aln.start, aln.end-aln.start);

        align.align(aln, errorRate);
        if (aln.qstr.size() == 0) {
            cnspos[i].setMinMax(0, 0);
            continue;
        }
        cnspos[i].setMinMax(aln.start, aln.end);
        dagcon::Alignment norm = normalizeGaps(aln);

        // not thread safe to add to graph concurrently, so lock while adding
#pragma omp critical (graphAdd)
        ag.addAln(norm);
    }

    // merge the nodes and call consensus
    ag.mergeNodes();
    std::string cns = ag.consensus(1);

    // save consensus
    resizeArrayPair(tig->_gappedBases, tig->_gappedQuals, 0, tig->_gappedMax, (uint32) cns.length() + 1, resizeArray_doNothing);
    std::string::size_type len = 0;
    for(len = 0; len < cns.size(); len++) {
      tig->_gappedBases[len] = cns[len];
      tig->_gappedQuals[len] = CNS_MIN_QV;
    }
    //  Terminate the string.
    tig->_gappedBases[len] = 0;
    tig->_gappedQuals[len] = 0;
    tig->_gappedLen = len;
    tig->_layoutLen = len;

    assert(len < tig->_gappedMax);
    return true;
}


bool
unitigConsensus::generateQuick(tgTig                     *tig_,
                               map<uint32, gkRead *>     *inPackageRead_,
                               map<uint32, gkReadData *> *inPackageReadData_) {

  tig      = tig_;
  numfrags = tig->numberOfChildren();

  if (initialize(inPackageRead_, inPackageReadData_) == FALSE) {
    fprintf(stderr, "generateMultiAlignment()--  Failed to initialize for tig %u with %u children\n", tig->tigID(), tig->numberOfChildren());
    return(false);
  }

  //  The quick variety doesn't generate alignments, it just pastes read bases into consensus.  It
  //  still needs to find a placement for the read, and it uses the other placed reads for that.

  while (moreFragments()) {
    reportStartingWork();

    piid = -1;
    bool placed = computePositionFromLayout();

    gkRead  *read    = gkpStore->gkStore_getRead(utgpos[tiid].ident());
    uint32   readLen = read->gkRead_sequenceLength();

    uint32 start = cnspos[tiid].min();
    uint32 end = cnspos[tiid].max();

    // if we couldn't place the read, fall back to utg positions
    if (placed == false) {
      start = abacus->numberOfColumns() - 1;
      end   = start + readLen;
    }

    uint32   bHang   = end - abacus->numberOfColumns();
    if (bHang <= 0) {
       //  this read doesn't add anything, skip it
       continue;
    }

    //  check if our positions are wonky, adjust the end to match reality
    if (start > abacus->numberOfColumns()) {
      start = abacus->numberOfColumns();
    }
    if (end - start > readLen) {
      end = start + readLen;
    }
    cnspos[tiid].setMinMax(start, end);

    //  appendBases() will append only if new bases are needed.  Otherwise, it just
    //  sets the first/last bead position.

    abacus->appendBases(tiid,
                        cnspos[tiid].min(),
                        cnspos[tiid].max());

    //  I _think_ we need to rebuild iff bases are added.  This also resets positions for each read.
    //  Until someone complains this is too slow, it's left in.

    refreshPositions();
  }

  generateConsensus(tig);

  return(true);
}





int
unitigConsensus::initialize(map<uint32, gkRead *>     *inPackageRead,
                            map<uint32, gkReadData *> *inPackageReadData) {

  int32 num_columns = 0;
  //int32 num_bases   = 0;

  if (numfrags == 0) {
    fprintf(stderr, "utgCns::initialize()-- unitig has no children.\n");
    return(false);
  }

  utgpos = new tgPosition [numfrags];
  cnspos = new tgPosition [numfrags];

  memcpy(utgpos, tig->getChild(0), sizeof(tgPosition) * numfrags);
  memcpy(cnspos, tig->getChild(0), sizeof(tgPosition) * numfrags);

  traceLen   = 0;
  trace      = new int32 [2 * AS_MAX_READLEN];

  traceABgn  = 0;
  traceBBgn  = 0;

  memset(trace, 0, sizeof(int32) * 2 * AS_MAX_READLEN);

  abacus     = new abAbacus();

  //  Clear the cnspos position.  We use this to show it's been placed by consensus.
  //  Guess the number of columns we'll end up with.
  //  Initialize abacus with the reads.

  for (int32 i=0; i<numfrags; i++) {
    cnspos[i].setMinMax(0, 0);

    num_columns  = (utgpos[i].min() > num_columns) ? utgpos[i].min() : num_columns;
    num_columns  = (utgpos[i].max() > num_columns) ? utgpos[i].max() : num_columns;

    abacus->addRead(gkpStore,
                    utgpos[i].ident(),
                    utgpos[i]._askip, utgpos[i]._bskip,
                    utgpos[i].isReverse(),
                    inPackageRead,
                    inPackageReadData);
  }

  //  Check for duplicate reads

  {
    set<uint32>  dupFrag;

    for (uint32 i=0; i<numfrags; i++) {
      if (utgpos[i].isRead() == false) {
        fprintf(stderr, "unitigConsensus()-- Unitig %d FAILED.  Child %d is not a read.\n",
                tig->tigID(), utgpos[i].ident());
        return(false);
      }

      if (dupFrag.find(utgpos[i].ident()) != dupFrag.end()) {
        fprintf(stderr, "unitigConsensus()-- Unitig %d FAILED.  Child %d is a duplicate.\n",
                tig->tigID(), utgpos[i].ident());
        return(false);
      }

      dupFrag.insert(utgpos[i].ident());
    }
  }

  //  Initialize with the first read.

  abacus->applyAlignment(0, 0, 0, NULL, 0);

  //  And set the placement of the first read.

  cnspos[0].setMinMax(0, abacus->numberOfColumns());

  return(true);
}



int
unitigConsensus::computePositionFromAnchor(void) {

  assert(piid == -1);

  uint32 anchor = utgpos[tiid].anchor();

  if (anchor == 0)
    //  No anchor?!  Damn.
    goto computePositionFromAnchorFail;

  for (piid = tiid-1; piid >= 0; piid--) {
    abSequence *aseq = abacus->getSequence(piid);

    if (anchor != aseq->gkpIdent())
      //  Not the anchor.
      continue;

    if ((cnspos[piid].min() == 0) &&
        (cnspos[piid].max() == 0))
      //  Is the anchor, but that isn't placed.
      goto computePositionFromAnchorFail;

    if ((utgpos[piid].max() < utgpos[tiid].min()) ||
        (utgpos[tiid].max() < utgpos[piid].min())) {
      //  Is the anchor, and anchor is placed, but the anchor doesn't agree with the placement.
      if (showPlacement())
        fprintf(stderr, "computePositionFromAnchor()-- anchor %d at utg %d,%d doesn't agree with my utg %d,%d.  FAIL\n",
                anchor,
                utgpos[piid].min(), utgpos[piid].max(),
                utgpos[tiid].min(), utgpos[tiid].max());
      goto computePositionFromAnchorFail;
    }

    //  Scale the hangs by the change in the anchor size between bogart and consensus.

#if 0
    double   anchorScale = (double)(cnspos[piid].max() - cnspos[piid].min()) / (double)(utgpos[piid].max() - utgpos[piid].min());

    if (showPlacement())
      fprintf(stderr, "computePositionFromAnchor()--  frag %u in anchor %u -- hangs %d,%d -- scale %f -- final hangs %.0f,%.0f\n",
              utgpos[tiid].ident(),
              utgpos[piid].ident(),
              utgpos[tiid].aHang(),
              utgpos[tiid].bHang(),
              anchorScale,
              utgpos[tiid].aHang() * anchorScale,
              utgpos[tiid].bHang() * anchorScale);

    cnspos[tiid].setMinMax(cnspos[piid].min() + utgpos[tiid].aHang() * anchorScale,
                           cnspos[piid].max() + utgpos[tiid].bHang() * anchorScale);

    //  Hmmm, but if we shrank the read too much, add back in some of the length.  We want to end up
    //  with the read scaled by anchorScale, and centered on the hangs.

    int32   fragmentLength = utgpos[tiid].max() - utgpos[tiid].min();

    if ((cnspos[tiid].min() >= cnspos[tiid].max()) ||
        (cnspos[tiid].max() - cnspos[tiid].min() < 0.75 * fragmentLength)) {
      int32  center = (cnspos[tiid].min() + cnspos[tiid].max()) / 2;

      if (showPlacement()) {
        fprintf(stderr, "computePositionFromAnchor()--  frag %u in anchor %u -- too short.  reposition around center %d with adjusted length %.0f\n",
                utgpos[tiid].ident(),
                utgpos[piid].ident(),
                center, fragmentLength * anchorScale);
      }

      cnspos[tiid].setMinMax(center - fragmentLength * anchorScale / 2,
                             center + fragmentLength * anchorScale / 2);

      //  We seem immune to having a negative position.  We only use this to pull out a region from
      //  the partial consensus to align to.
      //
      //if (cnspos[tiid].min() < 0) {
      //  cnspos[tiid].min() = 0;
      //  cnspos[tiid].max() = fragmentLength * anchorScale;
      //}
    }
#else
    assert(0 <= utgpos[tiid].aHang());

    uint32  bgn = abacus->getColumn(piid, cnspos[piid].min() + utgpos[tiid].aHang() - cnspos[piid].min());
    uint32  end = abacus->getColumn(piid, cnspos[piid].max() + utgpos[tiid].bHang() - cnspos[piid].min());

    cnspos[tiid].setMinMax(bgn, end);
#endif

    assert(cnspos[tiid].min() < cnspos[tiid].max());

    if (showPlacement())
      fprintf(stderr, "computePositionFromAnchor()-- anchor %d at %d,%d --> beg,end %d,%d (tigLen %d)\n",
              anchor,
              cnspos[piid].min(), cnspos[piid].max(),
              cnspos[tiid].min(), cnspos[tiid].max(),
              abacus->numberOfColumns());
    return(true);
  }

 computePositionFromAnchorFail:
  cnspos[tiid].setMinMax(0, 0);

  piid = -1;

  return(false);
}



int
unitigConsensus::computePositionFromLayout(void) {
  int32   thickestLen = 0;

  assert(piid == -1);

  //  Find the thickest qiid overlap to any cnspos fragment
  for (int32 qiid = tiid-1; qiid >= 0; qiid--) {
    if ((utgpos[tiid].min() < utgpos[qiid].max()) &&
        (utgpos[tiid].max() > utgpos[qiid].min()) &&
        ((cnspos[qiid].min() != 0) ||
         (cnspos[qiid].max() != 0))) {
      cnspos[tiid].setMinMax(cnspos[qiid].min() + utgpos[tiid].min() - utgpos[qiid].min(),
                             cnspos[qiid].max() + utgpos[tiid].max() - utgpos[qiid].max());

      //  This assert triggers.  It results in 'ooo' below being negative, and we
      //  discard this overlap anyway.
      //
      //assert(cnspos[tiid].min() < cnspos[tiid].max());

      int32 ooo = MIN(cnspos[tiid].max(), abacus->numberOfColumns()) - cnspos[tiid].min();

#if 1
      if (showPlacement())
        fprintf(stderr, "computePositionFromLayout()-- layout %d at utg %d,%d cns %d,%d --> utg %d,%d cns %d,%d -- overlap %d\n",
                utgpos[qiid].ident(),
                utgpos[qiid].min(), utgpos[qiid].max(), cnspos[qiid].min(), cnspos[qiid].max(),
                utgpos[tiid].min(), utgpos[tiid].max(), cnspos[tiid].min(), cnspos[tiid].max(),
                ooo);
#endif

      //  Occasionally we see an overlap in the original placement (utgpos overlap) by after
      //  adjusting our fragment to the consensus position, we no longer have an overlap.  This
      //  seems to be caused by a bad original placement.
      //
      //  Example:
      //  utgpos[a] = 13480,14239    cnspos[a] = 13622,14279
      //  utgpos[b] = 14180,15062
      //
      //  Our placement is 200bp different at the start, but close at the end.  When we compute the
      //  new start placement, it starts after the end of the A read -- the utgpos say the B read
      //  starts 700bp after the A read, which is position 13622 + 700 = 14322....50bp after A ends.

      if ((cnspos[tiid].min() < abacus->numberOfColumns()) &&
          (thickestLen < ooo)) {
        thickestLen = ooo;

        assert(cnspos[tiid].min() < cnspos[tiid].max());  //  But we'll still assert cnspos is ordered correctly.

        int32 ovl   = ooo;
        int32 ahang = cnspos[tiid].min();
        int32 bhang = cnspos[tiid].max() - abacus->numberOfColumns();

        piid  = qiid;
      }
    }
  }

  //  If we have a VALID thickest placement, use that (recompute the placement that is likely
  //  overwritten -- ahang, bhang and piid are still correct).

  if (thickestLen >= minOverlap) {
    assert(piid != -1);

    cnspos[tiid].setMinMax(cnspos[piid].min() + utgpos[tiid].min() - utgpos[piid].min(),
                           cnspos[piid].max() + utgpos[tiid].max() - utgpos[piid].max());

    assert(cnspos[tiid].min() < cnspos[tiid].max());

    if (showPlacement())
      fprintf(stderr, "computePositionFromLayout()-- layout %d at %d,%d --> beg,end %d,%d (tigLen %d)\n",
              utgpos[piid].ident(),
              cnspos[piid].min(), cnspos[piid].max(),
              cnspos[tiid].min(), cnspos[tiid].max(),
              abacus->numberOfColumns());

    return(true);
  }

  cnspos[tiid].setMinMax(0, 0);

  piid = -1;

  return(false);
}



//  Occasionally we get a fragment that just refuses to go in the correct spot.  Search for the
//  correct placement in all of consensus, update ahang,bhang and retry.
//
//  We don't expect to have big negative ahangs, and so we don't allow them.  To unlimit this, use
//  "-fragmentLen" instead of the arbitrary cutoff below.
int
unitigConsensus::computePositionFromAlignment(void) {

  assert(piid == -1);

  int32        minlen      = minOverlap;
  int32        ahanglimit  = -10;

  abSequence  *seq         = abacus->getSequence(tiid);
  char        *fragment    = seq->getBases();
  int32        fragmentLen = seq->length();

  bool         foundAlign  = false;

  //
  //  Try NDalign.
  //

  if (foundAlign == false) {

    if (oaPartial == false)
      oaPartial = new NDalign(pedLocal, errorRate, 17);  //  partial allowed!

    oaPartial->initialize(0, abacus->bases(), abacus->numberOfColumns(), 0, abacus->numberOfColumns(),
                          1, fragment,        fragmentLen,               0, fragmentLen,
                          false);

    if ((oaPartial->findMinMaxDiagonal(minOverlap) == true) &&
        (oaPartial->findSeeds(false)               == true) &&
        (oaPartial->findHits()                     == true) &&
        (oaPartial->chainHits()                    == true) &&
        (oaPartial->processHits()                  == true)) {

      cnspos[tiid].setMinMax(oaPartial->abgn(), oaPartial->aend());

      //fprintf(stderr, "computePositionFromAlignment()-- cnspos[%3d] mid %d %d,%d (from NDalign)\n", tiid, utgpos[tiid].ident(), cnspos[tiid].min(), cnspos[tiid].max());

      foundAlign = true;
    }
  }

  //
  //  Fail.
  //

  if (foundAlign == false) {
    cnspos[tiid].setMinMax(0, 0);
    piid = -1;

    if (showAlgorithm())
      fprintf(stderr, "computePositionFromAlignment()-- Returns fail (no alignment).\n");
    return(false);
  }

  //  From the overlap and existing placements, find the thickest overlap, to set the piid and
  //  hangs, then reset the original placement based on that anchors original placement.
  //
  //  To work with fixFailures(), we need to scan the entire fragment list.  This isn't so bad,
  //  really, since before we were scanning (on average) half of it.

  assert(cnspos[tiid].min() < cnspos[tiid].max());

  int32   thickestLen = 0;

  for (int32 qiid = numfrags-1; qiid >= 0; qiid--) {
    if ((tiid != qiid) &&
        (cnspos[tiid].min() < cnspos[qiid].max()) &&
        (cnspos[tiid].max() > cnspos[qiid].min())) {
      int32 ooo = (MIN(cnspos[tiid].max(), cnspos[qiid].max()) -
                   MAX(cnspos[tiid].min(), cnspos[qiid].min()));

      if (thickestLen < ooo) {
        thickestLen = ooo;

        int32 ovl   = ooo;
        int32 ahang = cnspos[tiid].min();
        int32 bhang = cnspos[tiid].max() - abacus->numberOfColumns();

        piid  = qiid;
      }
    }
  }

  //  No thickest?  Dang.

  if (thickestLen == 0) {
    cnspos[tiid].setMinMax(0, 0);
    piid = -1;
    if (showAlgorithm())
      fprintf(stderr, "computePositionFromAlignment()-- Returns fail (no thickest).\n");
    return(false);
  }

  //  Success, yay!

  assert(piid != -1);

  if (showPlacement())
    fprintf(stderr, "computePositionFromAlignment()-- layout %d at %d,%d --> beg,end %d,%d (tigLen %d)\n",
            utgpos[piid].ident(),
            cnspos[piid].min(), cnspos[piid].max(),
            cnspos[tiid].min(), cnspos[tiid].max(),
            abacus->numberOfColumns());

  return(true);
}


void
unitigConsensus::generateConsensus(tgTig *tig) {

  abacus->recallBases(true);  //  Do one last base call, using the full works.

  abacus->refine(abAbacus_Smooth);
  abacus->mergeColumns(true);

  abacus->refine(abAbacus_Poly_X);
  abacus->mergeColumns(true);

  abacus->refine(abAbacus_Indel);
  abacus->mergeColumns(true);

  abacus->recallBases(true);  //  The bases are possibly all recalled, depending on the above refinements keeping things consistent.
  //abacus->refreshColumns();    //  Definitely needed, this copies base calls into _cnsBases and _cnsQuals.

  //  Copy the consensus and positions into the tig.

  abacus->getConsensus(tig);
  abacus->getPositions(tig);

  //  While we have fragments in memory, compute the microhet probability.  Ideally, this would be
  //  done in CGW when loading unitigs (the only place the probability is used) but the code wants
  //  to load sequence and quality for every fragment, and that's too expensive.
}



//  Update the position of each fragment in the consensus sequence.
//  Update the anchor/hang of the fragment we just placed.
void
unitigConsensus::refreshPositions(void) {

  for (int32 i=0; i<=tiid; i++) {
    if ((cnspos[i].min() == 0) &&
        (cnspos[i].max() == 0))
      //  Uh oh, not placed originally.
      continue;

    abColumn *fcol = abacus->readTofBead[i].column;
    abColumn *lcol = abacus->readTolBead[i].column;

    cnspos[i].setMinMax(fcol->position(),
                        lcol->position() + 1);

    assert(cnspos[i].min() >= 0);
    assert(cnspos[i].max() > cnspos[i].min());
  }

  if (piid >= 0)
    utgpos[tiid].setAnchor(utgpos[piid].ident(),
                           cnspos[tiid].min() - cnspos[piid].min(),
                           cnspos[tiid].max() - cnspos[piid].max());

  piid = -1;
}



//  Run abacus to rebuild the consensus sequence.  VERY expensive.
void
unitigConsensus::recomputeConsensus(bool display) {

  //abacus->recallBases(false);  //  Needed?  We should be up to date.

  abacus->refine(abAbacus_Smooth);
  abacus->mergeColumns(false);

  abacus->refine(abAbacus_Poly_X);
  abacus->mergeColumns(false);

  abacus->refine(abAbacus_Indel);
  abacus->mergeColumns(false);

  abacus->recallBases(false);  //  Possibly not needed.  If this is removed, the following refresh is definitely needed.
  //abacus->refreshColumns();    //  Definitely needed, this copies base calls into _cnsBases and _cnsQuals.

  refreshPositions();

  if (display)
    abacus->display(stderr);
}



//  This stub lets alignFragmnet() cleanup and return on alignment failures.  The original
//  implementation did the same thing with a goto to the end of the function.  Opening up the if
//  statements exposed variable declarations that prevented the goto from compiling.
//
bool
unitigConsensus::alignFragmentFailure(void) {
  cnspos[tiid].setMinMax(0, 0);
  piid = -1;

  if (showAlgorithm())
    fprintf(stderr, "alignFragment()-- No alignment found.\n");

  return(false);
}


//  Generates an alignment of the current read to the partial consensus.
//  The primary output is a trace stored in the object data.

bool
unitigConsensus::alignFragment(bool forceAlignment) {

  assert((cnspos[tiid].min() != 0) || (cnspos[tiid].max() != 0));
  assert(piid != -1);

  assert(cnspos[tiid].min() < cnspos[tiid].max());

  abSequence *bSEQ    = abacus->getSequence(tiid);
  char       *fragSeq = bSEQ->getBases();
  int32       fragLen = bSEQ->length();

  //  Decide on how much to align.  Pick too little of consensus, and we leave some of the read
  //  unaligned.  Pick too much, and the read aligns poorly.
  //
  //  endTrim is trimmed from the 3' of the read.  This is the stuff we don't expect to align to consensus.
  //  bgnExtra is trimmed from the 5' of consensus.  Same idea as endTrim.
  //  endExtra is trimmed from the 3' of consensus.  Only for contained reads.
  //
  //  These values are adjusted later, in trimStep increments, based on the alignments returned.
  //
  //  Of the two choices, making the two Extra's small at the start is probably safer.  That case is
  //  easy to detect, and easy to fix.
  //

  //  The expectedAlignLen is almost always an underestimate.  Any gaps inserted will make the real
  //  alignment length longer.  This used to be multiplied by the error rate.
  int32  expectedAlignLen = cnspos[tiid].max() - cnspos[tiid].min();

  //  If the read is contained, the full read is aligned.
  //  Otherwise, an extra 1/32 of the align length is added for padding.
  int32  fragBgn = 0;
  int32  fragEnd = (cnspos[tiid].max() < abacus->numberOfColumns()) ? (fragLen) : (33 * expectedAlignLen / 32);

  if (fragEnd > fragLen)
    fragEnd = fragLen;

  //  Given the usual case of using an actual overlap to a read in the multialign to find the region
  //  to align to, we expect that region to be nearly perfect.  Thus, we shouldn't need to extend it
  //  much.  If anything, we'll need to extend the 3' end of the read.

  int32  bgnExtra = 10;    //  Start with a small 'extra' allowance, easy to make bigger.
  int32  endExtra = 10;    //

  int32  trimStep = max(10, expectedAlignLen / 50);  //  Step by 10 bases or do at most 50 steps.

  //  Find an alignment!

  bool  allowedToTrim = true;

  assert(abacus->bases()[abacus->numberOfColumns()] == 0);  //  Consensus must be NUL terminated
  assert(fragSeq[fragLen]                           == 0);  //  The read must be NUL terminated

 alignFragmentAgain:

  //  Truncate consensus and the read to prevent false alignments.

  if (cnspos[tiid].max() + endExtra > abacus->numberOfColumns())
    endExtra = abacus->numberOfColumns() - cnspos[tiid].max();

  int32 cnsBgn     = MAX(0, cnspos[tiid].min() - bgnExtra);   //  Start position in consensus
  int32 cnsEnd     = cnspos[tiid].max() + endExtra;           //  Truncation of consensus
  int32 cnsEndBase = abacus->bases()[cnsEnd];                 //  Saved base (if not truncated, it's the NUL byte at the end)

  char *aseq         = abacus->bases() + cnsBgn;
  char *bseq         = fragSeq;

  char  fragEndBase  = bseq[fragEnd];                         //  Saved base

  abacus->bases()[cnsEnd] = 0;                                //  Do the truncations.
  bseq[fragEnd]           = 0;

  //  Report!

  if (showAlgorithm())
    fprintf(stderr, "alignFragment()-- Allow bgnExtra=%d and endExtra=%d (cnsBgn=%d cnsEnd=%d cnsLen=%d) (fragBgn=0 fragEnd=%d fragLen=%d)\n",
            bgnExtra, endExtra, cnsBgn, cnsEnd, abacus->numberOfColumns(), fragEnd, fragLen);

  //  Create new aligner object.  'Global' in this case just means to not stop early, not a true global alignment.

  if (oaFull == false)
    oaFull = new NDalign(pedGlobal, errorRate, 17);

  oaFull->initialize(0, aseq, cnsEnd  - cnsBgn,   0, cnsEnd  - cnsBgn,
                     1, bseq, fragEnd - fragBgn,  0, fragEnd - fragBgn,
                     false);

  //  Generate a null hit, then align it and then realign, from both endpoints, and save the better
  //  of the two.

  if ((oaFull->makeNullHit() == true) &&
      (oaFull->processHits() == true)) {
    if (showAlignments())
      oaFull->display("utgCns::alignFragment()--", true);

    oaFull->realignBackward(showAlgorithm(), showAlignments());
    oaFull->realignForward (showAlgorithm(), showAlignments());
  }

  //  Restore the bases we removed to end the strings early.

  if (cnsEndBase)   abacus->bases()[cnsEnd] = cnsEndBase;
  if (fragEndBase)  bseq[fragEnd]           = fragEndBase;

  //  If no alignment, bail.

  if (oaFull->length() == 0)
    return(alignFragmentFailure());

  //
  //  Check quality and fail if it sucks.
  //

  bool  isBad = oaFull->scanDeltaForBadness(showAlgorithm(), showAlignments());

  //  Check for bad (under) trimming of input sequences.
  //
  //  If the alignment is bad, and we hit the start of the consensus sequence (or the end of
  //  same), chances are good that the aligner returned a (higher scoring) global alignment instead
  //  of a (lower scoring) local alignment.  Trim off some of the extension and try again.

  if ((allowedToTrim == true) && (isBad == true) && (oaFull->ahg5() == 0) && (bgnExtra > 0)) {
    int32  adj = (bgnExtra < trimStep) ? 0 : bgnExtra - trimStep;

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- alignment is bad, hit the trimmed start of consensus, decrease bgnExtra from %u to %u\n", bgnExtra, adj);

    bgnExtra = adj;
    goto alignFragmentAgain;
  }

  if ((allowedToTrim == true) && (isBad == true) && (oaFull->ahg3() == 0) && (endExtra > 0)) {
    int32  adj = (endExtra < trimStep) ? 0 : endExtra - trimStep;

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- alignment is bad, hit the trimmed end of consensus, decrease endExtra from %u to %u\n", endExtra, adj);

    endExtra = adj;
    goto alignFragmentAgain;
  }

  //  Check for bad (over) trimming of input sequences.  Bad if:
  //     we don't hit the start of the read, and we chopped the start of consensus.
  //     we don't hit the end   of the read, and we chopped the end   of consensus.
  //     we do    hit the end   of the read, and we chopped the end   of the read.
  //     (we don't chop the start of the read, so the fourth possible case never happens)

  allowedToTrim = false;  //  No longer allowed to reduce bgnExtra or endExtra.  We'd hit infinite loops otherwise.

  if ((oaFull->bhg5() > 0) && (cnsBgn > 0)) {
    int32  adj = bgnExtra + 2 * oaFull->bhg5();

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- hit the trimmed start of consensus, increase bgnExtra from %u to %u\n", bgnExtra, adj);

    bgnExtra = adj;
    goto alignFragmentAgain;
  }

  if ((oaFull->bhg3() > 0) && (cnsEnd < abacus->numberOfColumns())) {
    int32  adj = endExtra + 2 * oaFull->bhg3();

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- hit the trimmed end of consensus, increase endExtra from %u to %u\n", endExtra, adj);

    endExtra = adj;
    goto alignFragmentAgain;
  }

  if ((oaFull->bhg3() == 0) && (fragEnd < fragLen)) {
    int32  adj = (fragEnd + trimStep < fragLen) ? fragEnd + trimStep : fragLen;

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- hit the trimmed end of the read, increase fragEnd from %d to %d\n", fragEnd, adj);

    fragEnd = adj;
    goto alignFragmentAgain;
  }

  //  If we get here, and it's still bad, well, not much we can do.

  if ((forceAlignment == false) && (isBad == true)) {
    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- alignment bad after realigning\n");
    return(alignFragmentFailure());
  }

  if ((forceAlignment == false) && (oaFull->erate() > errorRate)) {
    if (showAlgorithm()) {
      fprintf(stderr, "utgCns::alignFragment()-- alignment is low quality: %f > %f\n",
              oaFull->erate(), errorRate);
      oaFull->display("utgCns::alignFragment()-- ", true);
    }
    return(alignFragmentFailure());
  }

  //  Otherwise, its a good alignment.  Process the trace to the 'consensus-format' and return true.
  //
  //  Set the begin points of the trace.  We probably need at least one of abgn and bbgn to be
  //  zero.  If both are nonzero, then we have a branch in the alignment.
  //
  //  If traceABgn is negative, we insert gaps into A before it starts, but I'm not sure how that works.
  //
  //  Overlap encoding:
  //    Add N matches or mismatches.  If negative, insert a base in the first sequence.  If
  //    positive, delete a base in the first sequence.
  //
  //  Consensus encoding:
  //    If negative, align (-trace - apos) bases, then add a gap in A.
  //    If positive, align ( trace - bpos) bases, then add a gap in B.
  //

  if (oaFull->abgn() > 0)   assert(oaFull->bbgn() == 0);  //  read aligned fully if consensus isn't
  if (oaFull->bbgn() > 0)   assert(oaFull->abgn() == 0);  //  read extends past the begin, consensus aligned fully
  if (oaFull->bbgn() > 0)   assert(cnsBgn == 0);          //  read extends past the begin, consensus not trimmed at begin


  traceABgn = cnsBgn + oaFull->abgn() - oaFull->bbgn();
  traceBBgn =            oaFull->bbgn();

  int32   apos = oaFull->abgn();
  int32   bpos = oaFull->bbgn();

  traceLen = 0;

  for (uint32 ii=0; ii<oaFull->deltaLen(); ii++, traceLen++) {
    if (oaFull->delta()[ii] < 0) {
      apos += -oaFull->delta()[ii] - 1;
      bpos += -oaFull->delta()[ii];

      trace[traceLen] = -apos - cnsBgn - 1;

    } else {
      apos +=  oaFull->delta()[ii];
      bpos +=  oaFull->delta()[ii] - 1;  // critical

      trace[traceLen] = bpos + 1;
    }
  }

  trace[traceLen] = 0;

  return(true);
}
