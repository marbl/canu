
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

#include "trimReads.H"

#include "intervals.H"


bool
largestCovered(ovOverlap    *ovl,
               uint32        ovlLen,
               uint32        readID,
               uint32        readLen,
               uint32 UNUSED(ibgn),
               uint32        iend,
               uint32       &fbgn,
               uint32       &fend,
               char         *logMsg,
               uint32        errorValue,
               uint32        minOverlap,
               uint32        minCoverage,
               uint32        minReadLength) {

  bool verbose = false;

  logMsg[0] = 0;

  assert(readID == ovl[0].a_iid);
  assert(ovlLen > 0);

  intervalList<uint32>  IL;
  intervalList<uint32>  ID;
  int32                 iid = readID;

  uint32                nSkip = 0;
  uint32                nUsed = 0;

  for (uint32 i=0; i<ovlLen; i++) {
    uint32 tbgn = ovl[i].a_bgn();
    uint32 tend = ovl[i].a_end();

    assert(tbgn < tend);
    assert(iid == ovl[i].a_iid);

    if (ovl[i].evalue() > errorValue || (tend-tbgn < minOverlap)) {
      //  Overlap is crappy.
      //fprintf(stderr, "skip %2u\n", i);
      nSkip++;
      continue;
    }

//fprintf(stderr, "Processing interval %d - %d in read %d and currently I have %d stored\n", tbgn, tend, readID, IL.numberOfIntervals());
// check if this overlap ends at a point we already captured
bool doSkip = false;
for (uint32 j = 0; j < IL.numberOfIntervals(); j++) {
//fprintf(stderr, "Comparing interval from %d - %d in read %d to %d - %d with start bounds %d - %d\n", tbgn, tend, readID, IL.lo(j), IL.hi(j), (int32)(tbgn - floor(minOverlap/10)), (int32)(tbgn + floor(minOverlap/10)));
   if (tbgn > 15 && (int32)(tbgn - floor(minOverlap/50)) <= (int32)IL.lo(j) && (int32)(tbgn + floor(minOverlap/50)) >= (int32)IL.lo(j)) {
       doSkip = true;
//      fprintf(stderr, "Skip interval in read %ul from %ul - %ul because it is too close to %ul - %ul\n", readID, tbgn, tend, IL.lo(j), IL.hi(j));
      break;
   } 
//fprintf(stderr, "Comparing interval from %d - %d in read %d to %d - %d with end bounds %d - %d\n", tbgn, tend, readID, IL.lo(j), IL.hi(j), (int32)(tend - floor(minOverlap/10)), (int32)(tend + floor(minOverlap/10)));

   if (tend < readLen-15 && (int32)(tend - floor(minOverlap/50)) <= IL.hi(j) && (int32)(tend + floor(minOverlap/50)) >= IL.hi(j)) {
      doSkip = true;
//      fprintf(stderr, "Skip interval in read %ul from %ul - %ul because it is too close to %ul - %ul\n", readID, tbgn, tend, IL.lo(j), IL.hi(j));
      break;
   }
}
if (doSkip) {
      nSkip++;
continue;
}

    //fprintf(stderr, "save %2u from %d - %d\n", i, tbgn, tend);
    nUsed++;


uint32 trimLen = floor((tend-tbgn) * 0.05);
if (trimLen > 250) trimLen = 250;
//fprintf(stderr, "For interval from %d - %d the trim is %d\n", tbgn, tend, trimLen);

    if (tend + trimLen  >= readLen)
       IL.add(tbgn, tend - tbgn);
    else
       IL.add(tbgn, tend - tbgn - trimLen);
  }

  if (verbose)
    for (uint32 it=0; it<IL.numberOfIntervals(); it++)
      fprintf(stderr, "IL - %3u - %5u %5u\n", readID, IL.lo(it), IL.hi(it));

  //  I thought I'd allow low coverage at the end of the read, but not internally, but that is hard,
  //  and largely unnecessary.  We'll just not be assembling at low coverage joins, which is
  //  acceptable.

  if (minCoverage > 0) {
    intervalDepth<uint32>  DE(IL);

    uint32  it = 0;
    uint32  ib = 0;
    uint32  ie = 0;

    while (it < DE.numberOfIntervals()) {
      if (verbose)
        fprintf(stderr, "DE - %3u - %5u %5u depth %u\n", readID, DE.lo(it), DE.hi(it), DE.depth(it));

      if (DE.depth(it) < minCoverage) {
        //  Dropped below good coverage depth.  If we have an interval, save it.  Reset.
        if (ie > ib) {
          //fprintf(stderr, "AD1 %d-%d len %d\n", ib, ie, ie - ib);
          ID.add(ib, ie - ib);
        }
        ib = 0;
        ie = 0;

      } else if ((ib == 0) && (ie == 0)) {
        //  Depth is good.  If no current interval, make a new one.
        ib = DE.lo(it);
        ie = DE.hi(it);
        //fprintf(stderr, "NE1 %d-%d len %d\n", ib, ie, ie - ib);

      } else if (ie == DE.lo(it)) {
        //  Depth is good.  If this interval is adjacent to the current, extend.
        ie = DE.hi(it);
        //fprintf(stderr, "EXT %d-%d len %d\n", ib, ie, ie - ib);

      } else {
        //  Depth is good, but we just had a gap in coverage.  Save any current interval.  Reset.
        if (ie > ib) {
          //fprintf(stderr, "AD2 %d-%d len %d\n", ib, ie, ie - ib);
          ID.add(ib, ie - ib);
        }
        ib = DE.lo(it);
        ie = DE.hi(it);
        //fprintf(stderr, "NE2 %d-%d len %d\n", ib, ie, ie - ib);
      }

      it++;
    }

    if (ie > ib) {
      //fprintf(stderr, "AD3 %d-%d len %d\n", ib, ie, ie - ib);
      ID.add(ib, ie - ib);
    }
  }

  //  Now that we've created depth, merge the intervals.

  IL.merge(minOverlap);

  //  IL - covered interavls enforcing a minimum overlap size (these can overlap)
  //  ID - covered intervals enforcing a minimum depth (these cannot overlap)
  //
  //  Create new intervals from the intersection of IL and ID.
  //
  //  This catches one nasty case, where a thin overlap has more than minDepth coverage.
  //
  //         -------------               3x coverage
  //          -------------              all overlaps 1 or 2 dashes long
  //                     ---------
  //                      -----------

  if (minCoverage > 0) {
    intervalList<uint32> FI;

    uint32  li = 0;
    uint32  di = 0;

    while ((li < IL.numberOfIntervals()) &&
           (di < ID.numberOfIntervals())) {
      uint32   ll = IL.lo(li);
      uint32   lh = IL.hi(li);
      uint32   dl = ID.lo(di);
      uint32   dh = ID.hi(di);
      uint32   nl  = 0;
      uint32   nh  = 0;

      //  If they intersect, make a new region

      if ((ll <= dl) && (dl < lh)) {
        nl = dl;
        nh = (lh < dh) ? lh : dh;
      }

      if ((dl <= ll) && (ll < dh)) {
        nl = ll;
        nh = (lh < dh) ? lh : dh;
      }

      if (nl < nh)
        FI.add(nl, nh - nl);

      //  Advance the list with the earlier region.

      if (lh <= dh)
        //  IL ends at or before ID
        li++;

      if (dh <= lh) {
        //  ID ends at or before IL
        di++;
      }
    }

    //  Replace the intervals to use with the intersection.

    IL = FI;
  }

  ////////////////////////////////////////

  if (IL.numberOfIntervals() == 0) {
    if (verbose)
      fprintf(stderr, "no high quality overlaps\n");
    strcpy(logMsg, "\tno high quality overlaps");
    return(false);
  }

  sprintf(logMsg, "\tskipped %u overlaps; used %u overlaps", nSkip, nUsed);

  fbgn = IL.lo(0);
  fend = IL.hi(0);
  if (verbose)
    fprintf(stderr, "IL - %3u - %5u %5u  SAVE\n", readID, fbgn, fend);

  for (uint32 it=1; it<IL.numberOfIntervals(); it++) {
    if (IL.hi(it) - IL.lo(it) > fend - fbgn) {
      fbgn = IL.lo(it);
      fend = IL.hi(it);
      if (verbose)
        fprintf(stderr, "IL - %3u - %5u %5u  SAVE\n", readID, IL.lo(it), IL.hi(it));
    } else {
      if (verbose)
        fprintf(stderr, "IL - %3u - %5u %5u\n", readID, IL.lo(it), IL.hi(it));
    }
  }

  return(true);
}
