
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_PopulateUnitig.H"


static
void
populateUnitig(Unitig           *unitig,
               BestEdgeOverlap  *bestnext) {

  assert(unitig->getLength() > 0);

  if ((bestnext == NULL) || (bestnext->readId() == 0)) {
    if (logFileFlagSet(LOG_BUILD_UNITIG)) {
      writeLog("                nothing\n");
      writeLog("\n");
    }
    return;
  }

  //  The ID of the last read in the unitig, and the end we should walk off of it.
  ufNode  read    = unitig->ufpath.back();
  int32   lastID  = read.ident;
  bool    last3p  = (read.position.bgn < read.position.end);

  uint32  nAdded  = 0;

  //  While there are reads to add AND those reads to add are not already in a unitig,
  //  construct a reverse-edge, and add the read.

  while ((bestnext->readId() != 0) &&
         (unitig->inUnitig(bestnext->readId()) == 0)) {
    BestEdgeOverlap  bestprev;

    //  Reverse nextedge (points from the unitig to the next read to add) so
    //  that it points from the next read to add back to something in the
    //  unitig.  If the reads are innie/outtie, we need to reverse the
    //  overlap to maintain that the A read is forward.

    if (last3p == bestnext->read3p())
      bestprev.set(lastID, last3p,  bestnext->bhang(),  bestnext->ahang(), bestnext->evalue());
    else
      bestprev.set(lastID, last3p, -bestnext->ahang(), -bestnext->bhang(), bestnext->evalue());

    //  'bestprev' points from read 'bestnext->readId()' end 'bestnext->read3p()'
    //                      to read 'lastID' end 'last3p'.
    //
    //  Compute the placement of the 'bestnext' read using the 'bestprev'
    //  edge (because that's how placeRead() wants to work).

    if (unitig->placeRead(read, bestnext->readId(), bestnext->read3p(), &bestprev)) {
      unitig->addRead(read);

      if (logFileFlagSet(LOG_BUILD_UNITIG))
        writeLog("                %s %7u; %c' ->\n",
                 (OG->isSpur(bestnext->readId()) == true) ? "spur" : "read",
                 bestnext->readId(),
                 bestnext->read3p() ? '5' : '3');


      nAdded++;
    }

    else {
      fprintf(stderr, "ERROR:  Failed to place read %d into BOG path.\n", read.ident);
      flushLog();
      assert(0);
    }

    //  Set up for the next read

    lastID  = read.ident;
    last3p  = (read.position.bgn < read.position.end);

    bestnext = OG->getBestEdgeOverlap(lastID, last3p);
  }

  if (logFileFlagSet(LOG_BUILD_UNITIG)) {
    if (bestnext->readId() == 0)
      writeLog("                nothing\n");
    else
      writeLog("                read %7u in tig %u\n",
               bestnext->readId(),
               unitig->inUnitig(bestnext->readId()));

    writeLog("tig %6u STOP after adding %u reads.\n", unitig->id(), nAdded);
    writeLog("\n");
  }
}



void
populateUnitig(TigVector &tigs,
               int32      fi) {

  //  Don't bother making tigs for deleted, contained, zombies, coverage gap,
  //  lopsided, et cetera, reads.

  if ((RI->readLength(fi) == 0) ||        //  Skip deleted
      (tigs.inUnitig(fi) != 0))           //  Skip placed
    return;

  if ((OG->isContained(fi)   == true) ||  //  Don't start a unitig if contained,
      (OG->isCoverageGap(fi) == true))    //  coverage gap, or lopsided.
    return;

  //  Grab the best edges for the candidate seed read.

  BestEdgeOverlap  *bestedge5 = OG->getBestEdgeOverlap(fi, false);
  BestEdgeOverlap  *bestedge3 = OG->getBestEdgeOverlap(fi, true);

  BestEdgeOverlap  *backedge5 = OG->getBestEdgeOverlap(bestedge5->readId(), bestedge5->read3p());
  BestEdgeOverlap  *backedge3 = OG->getBestEdgeOverlap(bestedge3->readId(), bestedge3->read3p());

  //  If this read has non-mutual best edges don't seed a tig with it.  There
  //  are (at least) four levels we could test:
  //    1) One end has a non-mutual best edge.
  //    2) Both ends have non-mutual best edge.
  //    3) One end of the read has no incoming best edge.
  //    4) Both ends have no incoming best edge.
  //
  //  On ecoli, the longest (and correct) path is seeded by a read that
  //  matches cases (1) and (3).
  //
  //  We've found one instance, in a human, where (4) occurred and led to a
  //  misassembly.

  //  This is case (1) and (2).  Case (2) is implemented.  Requiring case (1)
  //  didn't break ecoli, but it did change the layout slightly.  It 'feels
  //  like' case (1) might be too strict, but BPW hasn't tested on anything.
#if 1
  bool edgeTo5 = ((backedge5->readId() == fi) && (backedge5->read3p() == false));
  bool edgeTo3 = ((backedge3->readId() == fi) && (backedge3->read3p() ==  true));

  if ((edgeTo5 == false) ||
      (edgeTo3 == false)) {
    if (logFileFlagSet(LOG_BUILD_UNITIG))
      writeLog("tig ------ seed read %7u; non-mutual best edges, not using as a seed.  (edge to: 5' %s 3' %s)\n",
               fi,
               (edgeTo5) ? "yes" : "no",
               (edgeTo3) ? "yes" : "no");
    return;
  }
#endif

  //  This is case (3) and (4).  Case (3) is implemented.  The only argment
  //  against case (3) is a very very tiny DNA fragment with only two reads
  //  covering it.  In that case we'd never seed a tig (because both reads
  //  would have one end uncovered).
#if 0
  bool bestEdgeTo5 = false;
  bool bestEdgeTo3 = false;

  for (uint32 ii=1; ii<RI->numReads()+1; ii++) {
    BestEdgeOverlap  *e5 = OG->getBestEdgeOverlap(ii, false);
    BestEdgeOverlap  *e3 = OG->getBestEdgeOverlap(ii,  true);

    if (e5->readId() == fi)   ((e5->read3p() == false) ? bestEdgeTo5 : bestEdgeTo3) = true;
    if (e3->readId() == fi)   ((e3->read3p() == false) ? bestEdgeTo5 : bestEdgeTo3) = true;
  }

  if ((bestEdgeTo5 == false) ||
      (bestEdgeTo3 == false)) {
    if (logFileFlagSet(LOG_BUILD_UNITIG))
      writeLog("tig ------ seed read %7u; no best edge to me, not using as a seed.  (edge to: 5' %s 3' %s)\n",
               fi,
               (bestEdgeTo5) ? "yes" : "no",
               (bestEdgeTo3) ? "yes" : "no");
    return;
  }
#endif

  //  Add a first read -- to be 'compatable' with the old code, the first read is added
  //  reversed, we walk off of its 5' end, flip it, and add the 3' walk.

  Unitig *utg = tigs.newUnitig();

  utg->addRead(ufNode(fi, RI->readLength(fi), 0));

  //  Add reads as long as there is a path to follow...from the 3' end of the first read.

  if (bestedge3->readId()) {
    uint32  rid = utg->ufpath.back().ident;

    if (logFileFlagSet(LOG_BUILD_UNITIG))
      writeLog("tig %6u seed %s %7u; 5' ->\n",
               utg->id(), OG->isSpur(rid) ? "spur" : "read", rid);

    populateUnitig(utg, bestedge5);
  }

  //  Flip the tig around (and don't sort coordinates).

  utg->reverseComplement(false);

  //  Stick on reads on the beginning of the tig (that used to be the end).

  if (bestedge3->readId()) {
    uint32  rid = utg->ufpath.back().ident;

    if (logFileFlagSet(LOG_BUILD_UNITIG))
      writeLog("tig %6u seed %s %7u; 3' ->\n",
               utg->id(), OG->isSpur(rid) ? "spur" : "read", rid);

    populateUnitig(utg, bestedge3);
  }

  //  Enabling this reverse complement is known to degrade the assembly.  It is not known WHY it
  //  degrades the assembly.
  //
  //utg->reverseComplement(false);
}
