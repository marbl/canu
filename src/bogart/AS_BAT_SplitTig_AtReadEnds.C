
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

#include "AS_BAT_SplitTig.H"


class breakDesc {
public:
  breakDesc() {
    bgnAtAfter = 0;
    bgnBefore  = UINT32_MAX;
    bgnSet     = false;

    endAtAfter = 0;
    endBefore  = UINT32_MAX;
    endSet     = false;

    lowCoord   = INT32_MAX;
    nRepeat    = 0;
    nUnique    = 0;

    tig        = nullptr;
  };

  uint32  bgnAtAfter;   //  read must begin  after this coord.
  uint32  bgnBefore;    //  read must begin before this coord.
  bool    bgnSet;
 
  uint32  endAtAfter;   //  read must end  after this coord.
  uint32  endBefore;    //  read must end before this coord.
  bool    endSet;

  int32   lowCoord;     //  lowest old-tig coord, to offset reads to zero in the new tig

  uint32  nRepeat;      //  number of reads moved to this new tig
  uint32  nUnique;      //  number of reads moved to this new tig

  Unitig *tig;          //  the new tig associated with this block
};


void
splitTigAtReadEnds(TigVector                &tigs,
                   Unitig                   *tig,
                   vector<breakReadEnd>     &BE,
                   intervalList<int32>      &tigMarksR) {

  //  Convert the list of breakReadEnd into a list of breakDesc.
  //
  //  --------------------
  //    X-------------------
  //      -------------------
  //        Y-------------------
  //           --------------------
  //
  //  For breakHigh == false (X and Y), a new segment is created at the break
  //  point and extending to the right with no limit.  The previous segment
  //  is limited to have only reads that begin before the break point.
  //    Point X will:
  //      -Create a segment with bgnBefore set to X (this is the first segment)
  //      -Create a segment with bgnAtAfter set to X
  //    Point Y will:
  //      -Set bgnBefore in the current segment to Y
  //      -Create a segment with bgnAtAfter set to Y
  //
  //  --------------------
  //    -------------------W
  //      -------------------
  //        -------------------X
  //           --------------------
  //             Y-------------------
  //               --------------------
  //                 -------------------Z
  //                   --------------------
  //
  //  for breakHigh == true (X and Z), the current segment must not have any
  //  reads extending past the break point.  If the current segment does not have
  //  endBefore set, the current segment is used, otherwise a new segment is created.
  //    Point W will:
  //      -Set endBefore to W in the current segment (this is the first segment)
  //      -Set endAtAfter to W in all future regions - this prevents reads from being
  //       assigned to both W and X.  The alternate is to pick the first matching
  //       region.
  //    Point X will:
  //      -Since the current segment already has an endBefore, a new segment will be
  //       created with endBefore set to X.
  //    Point Y will:
  //      -Set bgnBefore in the current segment to Y  (as above)
  //      -Since the current segment has an endBefore, a new segment needs to be created
  //       to catch the reads between X and Y.  endAtAfter=X and bgnBefore=Y are set. 
  //      -Create a segment with bgnAtAfter set to Y  (as above)
  //    Point Z will:
  //      -Set endBefore to Z in the current segment (as above for point W)
  //
  //  In both examples, a final segment needs to be added to catch the reads
  //  at the end of the tig.

  //
  //  IMPORTANT!  breakReadEnd MUST come in read order.
  //
  //    A  --------------------X
  //    B   --------------------
  //    C    -----
  //    D     Y1------------------
  //    E      -----
  //    F       --------------------
  //    G        Y2------------------
  //    H         ----
  //    I          -------------------Z
  //    J           --------------------
  //
  //  The desired outcome is:
  //    reads that end before X to go into         -> 1  A
  //    reads that begin before Y1                 -> 2  B C      (C is ambiguous, either 1 or 2)
  //    reads that begin after Y1 (but before Y2)  -> 3  D E F    (E is ambiguous, either 1, 2 or 3)
  //    reads that begin after Y2 and end before Z -> 4  G H I    (H is ambiguous, either 1, 2, 3 or 4)
  //    reads that end after Z                     -> 5  J
  //
  //  But if we sort by position:
  //    reads that begin before Y1                 -> 1  A B C
  //    reads that begin after Y1 but before Y2    -> 2  D E F
  //    reads that begin after Y2 and end before X -> 3  H
  //    reads that end before Z                    -> 4  G I
  //    reads that end after Z                     -> 5  J
  //

  uint32     nb = 0;
  uint32     mb = BE.size() * 2;

  breakDesc *breaks = new breakDesc [mb];

  for (uint32 bp=0; bp<BE.size(); bp++) {
    if (BE[bp].breakHigh == false) {
      breaks[nb].bgnBefore = BE[bp].splitCoord;          //  Restrict current region to reads                   ----------  <- this read!
      breaks[nb].bgnSet    = true;                       //  that begin before the break point.                    Y----------

      if (breaks[nb].endSet == true) {                   //  If there is an end point set for the current       ----------X
        nb++;                                            //  block, make a new block to catch all the reads        ------------- <- this read!
        breaks[nb].endAtAfter = breaks[nb-1].endBefore;  //  that end after that block but begin before the          Y-------------
        breaks[nb].bgnBefore  = BE[bp].splitCoord;       //  block we're adding next.
      }

      nb++;                                              //  Add a new block, requiring reads
      breaks[nb].bgnAtAfter = BE[bp].splitCoord;         //  to begin at or after the point.
      breaks[nb].bgnSet     = true;                      //

      breaks[nb].endAtAfter = BE[bp].splitCoord;         //  Remove any end restriction set below.
    }

    else {
      if (breaks[nb].endSet == true)                     //  If there is already an end point set, make a new block.
        nb++;                                            //  Otherwise, we just add to the existing breakLow block.

      //if (breaks[nb].bgnSet == false) {                  //  This is optional, since bgn < end anyway.
      //  breaks[nb].bgnBefore = BE[bp].splitCoord + 1;
      //  breaks[nb].bgnSet    = true;
      //}

      breaks[nb].endBefore  = BE[bp].splitCoord + 1;     //  Require reads in this region to end
      breaks[nb].endSet     = true;                      //  before the point.

      for (uint32 nn=nb+1; nn<mb; nn++)                  //  Only reads ending after this point can be used in future regions.
        breaks[nn].endAtAfter = BE[bp].splitCoord + 1;   //  This prevents multiple placements in later breakHigh regions.
    }
  }

  if (breaks[nb].endSet == true) {                       //  If the last block added is from a high break
    nb++;                                                //  point, add another block to catch all the reads
    breaks[nb].endAtAfter = breaks[nb-1].endBefore;      //  that end after it.
  }

  nb++;                                                  //  The last block is always left "open".

  //  Display the rules we'll use for moving reads.

  writeLog("\n");
  writeLog("Rules for splitting tig %u:\n", tig->id());
  for (uint32 bp=0; bp<nb; bp++)
    writeLog("  %10s <= bgn < %-10s  AND  %10s <= end < %-10s --> block %u\n",
             toDec(breaks[bp].bgnAtAfter), (breaks[bp].bgnBefore <= tig->getLength()) ? toDec(breaks[bp].bgnBefore) : "end-of-tig",
             toDec(breaks[bp].endAtAfter), (breaks[bp].endBefore <= tig->getLength()) ? toDec(breaks[bp].endBefore) : "end-of-tig",
             bp);
  writeLog("\n");

  //  Iterate through reads, placing the read in the first bin it fits into.
  //  Continue through the other bins, just to count for ambiguous placements.
  //
  //  As iterating through, remember the lowest coord for each block, so we can
  //  offset reads.  Though this can be skipped if we let the cleanup routines
  //  fix things up.

  uint32    *rPlace  = new uint32 [tig->ufpath.size()];
  uint32    *nPlaces = new uint32 [tig->ufpath.size()];

  bool fail=false;

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode  &read     = tig->ufpath[fi];
    uint32   lo       = read.position.min();
    uint32   hi       = read.position.max();
    bool     isRepeat = false;

    rPlace [fi] = 0;
    nPlaces[fi] = 0;

    //  Decide if this read is in a repeat region, or considered unique.
    //  This is only so we can label newly created contigs are 'repeat'.

    for (uint32 rr=0; rr<tigMarksR.numberOfIntervals(); rr++)
      if (isContained(lo, hi, tigMarksR.lo(rr), tigMarksR.hi(rr)) == true) {
        isRepeat = true;
        break;
      }

    //  Scan all the regions, counting the number that this read can get
    //  placed into.  Also set the minimum coord seen in each block, so we
    //  can offset the reads when we move them to the new tig.

    for (uint32 bp=0; bp<nb; bp++) {
      if ((breaks[bp].bgnAtAfter <= lo) && (lo < breaks[bp].bgnBefore) && 
          (breaks[bp].endAtAfter <= hi) && (hi < breaks[bp].endBefore)) {
        rPlace [fi]  = bp;
        nPlaces[fi] += 1;
      }
    }

    //  Set the low coord for the new tig if the read is placed uniquely.

    if (nPlaces[fi] == 1) {
      uint32  bp = rPlace[fi];

      breaks[bp].lowCoord = min(breaks[bp].lowCoord, read.position.min());

      if (isRepeat == true)   breaks[bp].nRepeat++;
      else                    breaks[bp].nUnique++;
    }
  }

  //  Log any failures to place reads, then blow up.

  bool  blowUp = false;

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode  &read     = tig->ufpath[fi];

    if (nPlaces[fi] == 0) {
      writeLog("ERROR: read %8u %8d,%-8d was not placed into a block\n", read.ident, read.position.min(), read.position.max());
      blowUp = true;
      flushLog();
    }
  }
  assert(blowUp == false);

  //  Make new tigs for any that will have reads moved to it.

  writeLog("Tig Creation:\n");

  for (uint32 bp=0; bp<nb; bp++) {
    if (breaks[bp].nRepeat + breaks[bp].nUnique == 0) {
      writeLog("  Block %2u has %6u repeat and %6u unique reads, tig   not created.\n", bp, 0, 0);
    }

    else {
      breaks[bp].tig = tigs.newUnitig(false);

      breaks[bp].tig->_isBubble = tig->_isBubble;                              //  Inherit it's previous bubble status.
      breaks[bp].tig->_isRepeat = (breaks[bp].nRepeat > breaks[bp].nUnique);   //  Assign a new repeat/unique status.

      writeLog("  Block %2u has %6u repeat and %6u unique reads, tig %5u created.\n", bp, breaks[bp].nRepeat, breaks[bp].nUnique, breaks[bp].tig->id());
    }
  }

  writeLog("\n");

  //  And move reads to their new home.

  for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
    ufNode  &read  = tig->ufpath[fi];
    uint32   lo    = read.position.min();
    uint32   hi    = read.position.max();
    uint32   bp    = rPlace[fi];

    //  Do the easy case first.

    if (nPlaces[fi] == 1) {
      //writeLog("read %8u %8d,%-8d into block %2u\n", read.ident, read.position.min(), read.position.max(), bp);

      breaks[bp].tig->addRead(read, -breaks[bp].lowCoord, false);
      continue;
    }

    //  Otherwise, it wasn't placed.  Eject it from the tig and log failure.

    tigs.registerRead(read.ident);

    if (nPlaces[fi] == 0)
      writeLog("  WARNING: read %8u %8d,%-8d into no block -- ejected\n", read.ident, read.position.min(), read.position.max());
    else
      writeLog("  WARNING: read %8u %8d,%-8d into multiple blocks -- ejected\n", read.ident, read.position.min(), read.position.max());
  }

  //  That's it.

  writeLog("\n");

  delete [] nPlaces;
  delete [] rPlace;
  delete [] breaks;
}
