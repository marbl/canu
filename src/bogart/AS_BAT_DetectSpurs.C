
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



bool
detectSpur(TigVector &tigs,
           Unitig    *tig,
           uint32    &nEdges,
           uint32    &nPotential,
           uint32    &nVerified) {

  //  Grab the first read, and the edge out of the tig.  From that edge,
  //  figure out the tig that we intersect, and grab the read we intersect.
  //
  //  Note!  The other edge is NOT guaranteed to be back to the tig we're in
  //  (it usually will be), and this edge is also NOT guaranteed to be to a
  //  different tig (it usually will be).  So there really isn't anything
  //  useful we can test to make sure we've got the correct edge.

  ufNode           *fn      = tig->firstRead();
  BestEdgeOverlap  *edge    = OG->getBestEdgeOverlap(fn->ident, (fn->isForward() == true) ? false : true);

  if (edge->readId() == 0)
    return(false);

  uint32            isectID = tigs.inUnitig(edge->readId());
  Unitig           *isect   = tigs[isectID];

  if ((isectID == 0) ||
      (isect   == NULL))
    return(false);

  int32             ireadIdx =  tigs.ufpathIdx(edge->readId());
  ufNode           *iread    = &isect->ufpath[ireadIdx];

  writeLog("tig %6u read %7u %c' intersects tig %6u read %7u %c' (#%d).\n",
           tig->id(),   fn->ident, (fn->isForward() == true) ? '5' : '3',
           isect->id(), iread->ident, edge->read3p() ? '3' : '5', ireadIdx);
  nEdges++;

  //  edge->readId() -- the read we intersect.
  //  edge->read3p() -- true if we come into the 3' end of the read
  //
  //  There are four cases we that lead to a spur:
  //    the read is forward, we come into the 3' end, and it's near the start of the tig.
  //    the read is reverse, we come into the 5' end, and it's near the start of the tig.
  //
  //    the read is forward, we come into the 3' end, and it's near the end   of the tig.
  //    the read is reverse, we come into the 5' end, and it's near the end   of the tig.
  //
  uint32 nr    = isect->getNumReads();
  bool   isFwd = iread->isForward();
  bool   is3p  = edge->read3p();
  bool   isBgn = ((nr > 15) && (ireadIdx     < 5));
  bool   isEnd = ((nr > 15) && (ireadIdx + 5 > nr));

  bool   checkSpur = (((isFwd == true)  && (is3p == true)  && (isBgn == true )) ||
                      ((isFwd == false) && (is3p == false) && (isBgn == true )) ||
                      ((isFwd == true)  && (is3p == true)  && (isEnd == true)) ||
                      ((isFwd == false) && (is3p == false) && (isEnd == true)));

  writeLog("isFwd %u is3p %u isBgn %u isEnd %u nr %u\n", isFwd, is3p, isBgn, isEnd);

  if (checkSpur == false)
    return(false);

  writeLog("tig %6u read %7u %c' intersects tig %6u read %7u %c' and it's a potential spur.\n",
          tig->id(),   fn->ident, (fn->isForward() == true) ? '5' : '3',
          isect->id(), iread->ident, edge->read3p() ? '3' : '5');
  nPotential++;

  //  If the read at the end of the tig goes somewhere, it's not a spur.

  if (isBgn == true) {
    ufNode           *fread = isect->firstRead();
    BestEdgeOverlap  *e5    = OG->getBestEdgeOverlap(fread->ident, false);
    BestEdgeOverlap  *e3    = OG->getBestEdgeOverlap(fread->ident, true);

    if ((fread->isForward() == true) && (e5->readId() != 0))
      return(false);
    if ((fread->isReverse() == true) && (e3->readId() != 0))
      return(false);
  }

  if (isBgn == false) {
    ufNode           *lread = isect->lastRead();
    BestEdgeOverlap  *e5    = OG->getBestEdgeOverlap(lread->ident, false);
    BestEdgeOverlap  *e3    = OG->getBestEdgeOverlap(lread->ident, true);

    if ((lread->isForward() == true) && (e3->readId() != 0))
      return(false);
    if ((lread->isReverse() == true) && (e5->readId() != 0))
      return(false);
  }

  //  And so we have found a spur.

  writeLog("tig %6u read %7u %c' intersects tig %6u read %7u %c' and it's a VERIFIED spur.\n",
          tig->id(),   fn->ident, (fn->isForward() == true) ? '5' : '3',
          isect->id(), iread->ident, edge->read3p() ? '3' : '5');
  nVerified++;

  return(true);
}



void
detectSpurs(TigVector &tigs) {
  uint32   nTested       =   0;
  uint32   nEdges[2]     = { 0, 0 };
  uint32   nPotential[2] = { 0, 0 };
  uint32   nVerified[2]  = { 0, 0 };

  for (uint32 ti=1; ti<tigs.size(); ti++) {
    Unitig *tig = tigs[ti];

    //  If no tig, or small, or useless, don't let it chop off spurs.

    if ((tig                 == NULL) ||
        (tig->ufpath.size()  <= 10) ||
        (tig->_isUnassembled == true))
      continue;

    nTested++;

    detectSpur(tigs, tig, nEdges[0], nPotential[0], nVerified[0]);
    tig->reverseComplement();

    detectSpur(tigs, tig, nEdges[1], nPotential[1], nVerified[1]);
    tig->reverseComplement();
  }

  writeStatus("detectSpur() done.\n");
  writeStatus("tested         %4u\n", nTested);
  writeStatus("nEdges      5' %4u   3' %4u\n", nEdges[0],     nEdges[1]);
  writeStatus("nPotential     %4u      %4u\n", nPotential[0], nPotential[1]);
  writeStatus("nVerified      %4u      %4u\n", nVerified[0],  nVerified[1]);
}


