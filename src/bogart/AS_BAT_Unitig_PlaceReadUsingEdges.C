
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


#undef  DEBUG_PLACE_READ


ufNode
placeRead_contained(uint32           readId,
                    ufNode          &parent,
                    BestEdgeOverlap *edge) {

  bool   pFwd  = parent.position.isForward();
  int32  pMin  = parent.position.min();
  int32  pMax  = parent.position.max();

  assert(pMin < pMax);

  //  Reverse the overlap.  read3p here means the overlap is flipped.
  int32  ahang = (edge->read3p() == false) ? -edge->ahang() : edge->bhang();
  int32  bhang = (edge->read3p() == false) ? -edge->bhang() : edge->ahang();

  //  Depending on the parent orientation...
  //
  //     pMin         pMax        pMin         pMax
  //     ---------------->        <----------------
  //     ahang ----- bhang        bhang ----- ahang
  //     > 0         < 0          < 0         > 0

  int32 fMin = (pFwd == true) ? pMin + ahang : pMin - bhang;
  int32 fMax = (pFwd == true) ? pMax + bhang : pMax - ahang;

  //int32  fMin = pMin + ((read3p == false) ? -edge->ahang() : edge->bhang());   //  * intraScale
  //int32  fMax = pMax + ((read3p == false) ? -edge->bhang() : edge->ahang());   //  * interScale

  assert(fMin < fMax);

  //  We don't know the true length of the overlap, and our hang-based math tends to shrink reads.
  //  Reset the end coordinate using the actual length of the read.

#if 0
#warning NOT RESETTING fMax BASED ON READ LENGTH
  writeLog("placeCont()-- read %u %d-%d with hangs %d %d places read %u at %d-%d reset to %d\n",
           parent.ident,
           parent.position.min(), parent.position.max(),
           ahang, bhang,
           readId,
           fMin, fMax,
           fMin + RI->readLength(readId));
#endif
  fMax = fMin + RI->readLength(readId);

  //  Orientation is straightforward, based on the orient of the parent, and the flipped flag.

  bool fFwd = (((pFwd == true)  && (edge->read3p() == false)) ||  //  parent is fwd, olap is not flipped
               ((pFwd == false) && (edge->read3p() == true)));    //  parent is rev, olap is     flipped

  ufNode   read;

  read.ident        = readId;
  read.contained    = 0;
  read.parent       = edge->readId();       //  == parent->ident
  read.ahang        = 0;                    //  Not used in bogart, set on output
  read.bhang        = 0;                    //  Not used in bogart, set on output
  read.position.bgn = (fFwd) ? fMin : fMax;
  read.position.end = (fFwd) ? fMax : fMin;

#ifdef DEBUG_PLACE_READ
  writeLog("placeCont()-- parent %7d pos %7d,%7d -- edge to %7d %c' hangs %7d %7d -- read %7d C' -- placed %7d-%7d oriented %s %7d-%7d %f%% of length\n",
           parent.ident, parent.position.bgn, parent.position.end,
           edge->readId(), (edge->read3p()) ? '3' : '5', edge->ahang(), edge->bhang(),
           readId,
           fMin, fMax, (fFwd) ? "rev" : "fwd", read.position.bgn, read.position.end,
           100.0 * (read.position.max() - read.position.min()) / RI->readLength(readId));
#endif

  return(read);
}





ufNode
placeRead_dovetail(uint32           readId,
                   bool             read3p,
                   ufNode          &parent,
                   BestEdgeOverlap *edge) {

  //  We have an 'edge' from 'readId' end 'read3p' back to 'parent'.
  //  Use that to compute the placement of 'read'.

  bool   pFwd  = parent.position.isForward();
  int32  pMin  = parent.position.min();
  int32  pMax  = parent.position.max();

  assert(pMin < pMax);

  //  Scale the hangs based on the placed versus actual length of the parent read.

  //double  intraScale = (double)(pMax - pMin) / RI->readLength(parent.ident);  //  Within the parent read overlap
  //double  interScale = 1.0;                                                       //  Outside the parent read overlap

  //  We're given an edge from the read-to-place back to the parent.  Reverse the edge so it points
  //  from the parent to the read-to-place.
  //
  //  The canonical edge is from a forward parent to the child.
  //
  //      -P----\-->      +b
  //      +a  ---v--------C-
  //
  //  To reverse the edge:
  //
  //    If child is forward, swapping the order of the reads results in a canonical overlap.  The
  //    hangs become negative.
  //
  //      -P----\-->      +b    ---->   -a  ---/--------C>
  //      +a  ---v--------C>    ---->   -P----v-->      -b
  //
  //    If child is reverse, swapping the order of the reads results in a backwards canonical
  //    overlap, and we need to flip end-to-end also.  The hangs are swapped.
  //
  //      -P----\-->      +b    ---->   -C--------\-->  +a
  //      +a  <--v--------C-    ---->   +b      <--v----P-
  //
  int32  ahang = (read3p == false) ? -edge->ahang() : edge->bhang();
  int32  bhang = (read3p == false) ? -edge->bhang() : edge->ahang();

  //  The read is placed 'to the right' of the parent if
  //    pFwd == true  and edge points to 3' end
  //    pFwd == false and edge points to 5' end
  //
  bool  toRight = (pFwd == edge->read3p());

  //  If placing 'to the right', we add hangs.  Else, subtract the swapped hangs.

  int32  fMin = 0;
  int32  fMax = 0;

  if (toRight) {
    fMin = pMin + ahang;
    fMax = pMax + bhang;
  } else {
    fMin = pMin - bhang;
    fMax = pMax - ahang;
  }

  //int32  fMin = pMin + ((read3p == false) ? -edge->ahang() : edge->bhang());   //  * intraScale
  //int32  fMax = pMax + ((read3p == false) ? -edge->bhang() : edge->ahang());   //  * interScale

  assert(fMin < fMax);

  //  We don't know the true length of the overlap, and our hang-based math tends to shrink reads.
  //  Reset the end coordinate using the actual length of the read.

#if 0
#warning NOT RESETTING fMax BASED ON READ LENGTH
  writeLog("placeDovs()-- read %u %d-%d with hangs %d %d places read %u at %d-%d reset to %d\n",
           parent.ident,
           parent.position.min(), parent.position.max(),
           ahang, bhang,
           readId,
           fMin, fMax,
           fMin + RI->readLength(readId));
#endif
  fMax = fMin + RI->readLength(readId);


  //  Orientation is a bit more complicated, with eight cases (drawing pictures helps).
  //
  //    edge from read3p=true  to forward parent 3p -> reverse
  //    edge from read3p=false to reverse parent 3p -> reverse
  //    edge from read3p=false to forward parent 5p -> reverse
  //    edge from read3p=true  to reverse parent 5p -> reverse
  //
  //    edge from read3p=true  to reverse parent 3p -> forward
  //    edge from read3p=false to forward parent 3p -> forward
  //    edge from read3p=false to reverse parent 5p -> forward
  //    edge from read3p=true  to forward parent 5p -> forward
  //
  bool fFwd = (((read3p == true)  && (pFwd == true)  && (edge->read3p() == true))  ||
               ((read3p == false) && (pFwd == false) && (edge->read3p() == true))  ||
               ((read3p == false) && (pFwd == true)  && (edge->read3p() == false)) ||
               ((read3p == true)  && (pFwd == false) && (edge->read3p() == false))) ? false : true;

  ufNode   read;

  read.ident        = readId;
  read.contained    = 0;
  read.parent       = edge->readId();       //  == parent->ident
  read.ahang        = 0;                    //  Not used in bogart, set on output
  read.bhang        = 0;                    //  Not used in bogart, set on output
  read.position.bgn = (fFwd) ? fMin : fMax;
  read.position.end = (fFwd) ? fMax : fMin;

#ifdef DEBUG_PLACE_READ
  writeLog("placeDove()-- parent %7d pos %7d,%7d -- edge to %7d %c' hangs %7d %7d -- read %7d %c' -- placed %7d-%7d oriented %s %7d-%7d %f%% of length\n",
           parent.ident, parent.position.bgn, parent.position.end,
           edge->readId(), (edge->read3p()) ? '3' : '5', edge->ahang(), edge->bhang(),
           readId, (read3p) ? '3' : '5',
           fMin, fMax, (fFwd) ? "rev" : "fwd", read.position.bgn, read.position.end,
           100.0 * (read.position.max() - read.position.min()) / RI->readLength(readId));
#endif

  return(read);
}





//  Place a read into this tig using an edge from the read to some read in this tig.
//
bool
Unitig::placeRead(ufNode          &read,      //  output placement
                  uint32           readId,    //  id of read we want to place
                  bool             read3p,    //  end of read 'edge' is from, meaningless if contained
                  BestEdgeOverlap *edge) {    //  edge to read in this tig

  assert(readId > 0);
  assert(readId <= RI->numReads());

  read.ident             = readId;
  read.contained         = 0;
  read.parent            = 0;
  read.ahang             = 0;
  read.bhang             = 0;
  read.position.bgn      = 0;
  read.position.end      = 0;

  //  No best edge?  Hard to place without one.
  assert(edge != NULL);
  if (edge == NULL)
    return(false);

  //  Empty best edge?  Still hard to place.
  assert(edge->readId() != 0);
  if (edge->readId() == 0)
    return(false);

  //  Edge not pointing to a read in this tig?
  assert(inUnitig(edge->readId()) == id());
  if (inUnitig(edge->readId()) != id())
    return(false);

  //  Grab the index of the parent read.

  uint32 bidx = ufpathIdx(edge->readId());
  assert(edge->readId() == ufpath[bidx].ident);

  //  Now, just compute the placement and return success!

  if (((edge->ahang() >= 0) && (edge->bhang() <= 0)) ||
      ((edge->ahang() <= 0) && (edge->bhang() >= 0)))
    read = placeRead_contained(readId, ufpath[bidx], edge);
  else
    read = placeRead_dovetail(readId, read3p, ufpath[bidx], edge);

  return(true);
}
