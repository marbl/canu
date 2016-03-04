
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
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "splitReads.H"


//  Examine overlaps for a specific pattern indicating a read that flips back on itself.
//  This requires multiple overlaps from the same read, in opposite orientation, to work.
//
//  If this read isn't suspected to contain sub reads, all we can do is mark some
//  overlaps as junk.
//
//  If this read is suspected to contain sub reads, we want to mark the junction region.




uint32
intervalOverlap(uint32 b1, uint32 e1, uint32 b2, uint32 e2) {
  uint32 minmax = MIN(e1, e2);
  uint32 maxmin = MAX(b1, b2);

  if (minmax > maxmin)
    return(minmax - maxmin);

  return(0);
}



bool
doCheckSubRead(gkStore *gkp, uint32 id) {
  gkRead     *read = gkp->gkStore_getRead(id);
  gkLibrary  *libr = gkp->gkStore_getLibrary(read->gkRead_libraryID());

  return(libr->gkLibrary_checkForSubReads() == true);
}





//  Populate w->blist with intervals where a suspected subread junction occurs.

void
detectSubReads(gkStore               *gkp,
               workUnit              *w,
               FILE                  *subreadFile,
               bool                   subreadFileVerbose) {

  assert(w->adjLen > 0);
  assert(doCheckSubRead(gkp, w->id) == true);

  map<uint32, uint32>  secondIdx;
  map<uint32, uint32>  numOlaps;

  bool                 largePalindrome = false;
  intervalList<int32>  BAD;
  intervalList<int32>  BADall;

  //  Count the number of overlaps for each b_iid, and remember the last index.  There are supposed to
  //  be at most two overlaps per ID pair, so if we remember the last, and iterate through, we can
  //  get both.

  for (uint32 ii=0; ii<w->adjLen; ii++) {
    secondIdx[w->adj[ii].b_iid] = ii;
    numOlaps [w->adj[ii].b_iid]++;
  }

  //  Scan overlaps.  For any pair of b_iid, with overlaps in opposite directions, compute a 'bad'
  //  interval where a suspected flip occurs.

  for (uint32 ii=0; ii<w->adjLen; ii++) {
    adjOverlap  *aii = w->adj + ii;

    if (numOlaps[w->adj[ii].b_iid] == 1) {
      //  Only one overlap, can't indicate sub read!
      //if ((subreadFile) && (subreadFileVerbose))
      //  fprintf(subreadFile, "oneOverlap                 %u (%u-%u) %u (%u-%u) -- can't indicate subreads\n",
      //          w->adj[ii].a_iid, w->adj[ii].aovlbgn, w->adj[ii].aovlend, w->adj[ii].b_iid, w->adj[ii].bovlbgn, w->adj[ii].bovlend);
      continue;
    }

    //  We should never get more than two overlaps per read pair.
    if (numOlaps[w->adj[ii].b_iid] > 2) {
      fprintf(stderr, "ERROR: more overlaps than expected for pair %u %u.\n",
              w->adj[ii].a_iid, w->adj[ii].b_iid);
      continue;
    }
    assert(numOlaps[w->adj[ii].b_iid] == 2);

    uint32         jj = secondIdx[w->adj[ii].b_iid];
    adjOverlap   *ajj = w->adj + jj;

    assert(jj < w->adjLen);

    if (ii == jj) {
      //  Already did this one!
      //if ((subreadFile) && (subreadFileVerbose))
      //  fprintf(subreadFile, "sameOverlap                %u (%u-%u) %u (%u-%u)\n",
      //          w->adj[ii].a_iid, w->adj[ii].aovlbgn, w->adj[ii].aovlend, w->adj[ii].b_iid, w->adj[ii].bovlbgn, w->adj[ii].bovlend);
      continue;
    }

    //  The two overlaps should be for the same reads.

    assert(w->adj[ii].a_iid == w->adj[jj].a_iid);
    assert(w->adj[ii].b_iid == w->adj[jj].b_iid);

    //  And opposite orientations.

    if (w->adj[ii].flipped == w->adj[jj].flipped) {
      fprintf(stderr, "ERROR: same orient duplicate overlaps for pair %u %u\n",
              w->adj[ii].a_iid, w->adj[ii].b_iid);
      continue;
    }
    assert(w->adj[ii].flipped != w->adj[jj].flipped);


    bool  AcheckSub = (doCheckSubRead(gkp, w->adj[ii].a_iid) == true);
    bool  BcheckSub = (doCheckSubRead(gkp, w->adj[ii].b_iid) == true);

    assert(AcheckSub == true);  //  Otherwise we wouldn't be in this function!

    //  Decide what type of duplicate we have.
    //    Overlap on the A read -=> B read is potentially sub read containing -=> don't use overlaps
    //    Overlap on the B read -=> A read is potentially sub read containing -=> split this read

    uint32 Aoverlap = intervalOverlap(w->adj[ii].aovlbgn, w->adj[ii].aovlend, w->adj[jj].aovlbgn, w->adj[jj].aovlend);
    uint32 Boverlap = intervalOverlap(w->adj[ii].bovlbgn, w->adj[ii].bovlend, w->adj[jj].bovlbgn, w->adj[jj].bovlend);

    //  If there is no overlap anywhere, we're not sure what is going on.  This could be a genomic
    //  repeat.  Leave the overlaps alone.
    //
    if ((Aoverlap == 0) &&
        (Boverlap == 0))
      continue;

    //  Remember if the overlapping ovelap is large - we'll later check if the bad region falls
    //  within here, and if there are enough spanning reads not trim.  We also use this as one more
    //  count of BAD.
    //
    if ((AcheckSub) && (Aoverlap > 1000) &&
        (BcheckSub) && (Boverlap > 1000)) {
      uint32  dist = (w->adj[ii].a_iid > w->adj[ii].b_iid) ? (w->adj[ii].a_iid - w->adj[ii].b_iid) : (w->adj[ii].b_iid - w->adj[ii].a_iid);

      if (subreadFile)
        fprintf(subreadFile, "  II %8u (%6u-%6u) %8u (%6u-%6u)  JJ %8u (%6u-%6u) %8u (%6u-%6u) %s\n",
                w->adj[ii].a_iid, w->adj[ii].aovlbgn, w->adj[ii].aovlend, w->adj[ii].b_iid, w->adj[ii].bovlbgn, w->adj[ii].bovlend,
                w->adj[jj].a_iid, w->adj[jj].aovlbgn, w->adj[jj].aovlend, w->adj[jj].b_iid, w->adj[jj].bovlbgn, w->adj[jj].bovlend,
                (dist > 5) ? " PALINDROME WARNING--FAR-IID--WARNING" : "PALINDROME");

      largePalindrome  = true;
    }

#if 0
    //  Otherwise, if the overlaps overlap on both reads by significant chunks, don't believe
    //  either.  These are possibly both chimeric reads, at least PacBio junction reads.
    //
    //  Or an inverted repeat.
    //
    if ((AcheckSub) && (Aoverlap > 50) &&
        (BcheckSub) && (Boverlap > 50)) {
      if (subreadFile)
        fprintf(subreadFile, "BothOv     %u (%u-%u) %u (%u-%u)  %u (%u-%u) %u (%u-%u)\n",
                w->adj[ii].a_iid, w->adj[ii].aovlbgn, w->adj[ii].aovlend, w->adj[ii].b_iid, w->adj[ii].bovlbgn, w->adj[ii].bovlend,
                w->adj[jj].a_iid, w->adj[jj].aovlbgn, w->adj[jj].aovlend, w->adj[jj].b_iid, w->adj[jj].bovlbgn, w->adj[jj].bovlend);
    }
#endif

#if 0
    //  Stronger overlap in the A reads.  The B read looks like it has subreads, which is perfectly fine
    //  evidence for us.  Unless they span a junction.
    //
    if ((BcheckSub) && (Boverlap < Aoverlap)) {
      if (subreadFile)
        fprintf(subreadFile, "BcheckSub  %u (%u-%u) %u (%u-%u)  %u (%u-%u) %u (%u-%u)\n",
                w->adj[ii].a_iid, w->adj[ii].aovlbgn, w->adj[ii].aovlend, w->adj[ii].b_iid, w->adj[ii].bovlbgn, w->adj[ii].bovlend,
                w->adj[jj].a_iid, w->adj[jj].aovlbgn, w->adj[jj].aovlend, w->adj[jj].b_iid, w->adj[jj].bovlbgn, w->adj[jj].bovlend);
    }
#endif


    //  It looks like A has sub reads if the B read has a strong overlap in overlaps, and the A read does not
    //  have a strong overlap.

    if ((Aoverlap > 250) ||
        (Boverlap < 250))
      //  A strong overlap in the A read, there isn't a sub read junction we can identifiy, OR
      //  A weak overlap in the B read, and we expected the B read to align to both of the A subreads.
      continue;


    //  Decide on a region in the read that is suspected to contain the chimer junction.
    //
    //  In the true case: ii overlap is first on the read; bad region from the end of this overlap
    //  to the start of the jj overlap.
    //
    //  Note that sometimes overlaps extend through the junction.  This will just flip the region
    //  around.  We're expecting to find non-overlapping overlaps, but if we find overlapping ones,
    //  the bad interval is still between the end points.
    //
    //    -------------->                 ------------>
    //                    <---------  vs             <---------
    //
    uint32  badbgn = (w->adj[ii].aovlbgn < w->adj[jj].aovlbgn) ? w->adj[ii].aovlend : w->adj[jj].aovlend;
    uint32  badend = (w->adj[ii].aovlbgn < w->adj[jj].aovlbgn) ? w->adj[jj].aovlbgn : w->adj[ii].aovlbgn;

    if (badbgn > badend) {
      uint32  a = badbgn;
      badbgn = badend;
      badend = a;
    }
    assert(badbgn <= badend);

    if (subreadFile)
      fprintf(subreadFile, "  II %8u (%6u-%6u) %8u (%6u-%6u)  JJ %8u (%6u-%6u) %8u (%6u-%6u)  BAD %6u-%6u size %6u %s\n",
              w->adj[ii].a_iid, w->adj[ii].aovlbgn, w->adj[ii].aovlend, w->adj[ii].b_iid, w->adj[ii].bovlbgn, w->adj[ii].bovlend,
              w->adj[jj].a_iid, w->adj[jj].aovlbgn, w->adj[jj].aovlend, w->adj[jj].b_iid, w->adj[jj].bovlbgn, w->adj[jj].bovlend,
              badbgn, badend, badend - badbgn,
              (badend - badbgn <= SUBREAD_LOOP_MAX_SIZE) ? "(EVIDENCE)" : "(too far)");


    //  A true subread signature will have a small bad interval (10 bases) and largely agree on the
    //  interval.  False signature will have a large size, and not agree.  We only check for size
    //  though.
    //
    if (badend - badbgn <= SUBREAD_LOOP_MAX_SIZE)
      BAD.add(badbgn, badend - badbgn);

    //  Save all plausible pairs.
    //
    if (badend - badbgn <= SUBREAD_LOOP_EXT_SIZE)
      BADall.add(badbgn, badend - badbgn);
  }

  //
  //  Merge all the 'bad' intervals.  Save the merged intervals for later use.
  //

  BAD.merge();
  BADall.merge();



  for (uint32 bb=0; bb<BAD.numberOfIntervals(); bb++) {
    uint32  numSpan = 0;
    uint32  allHits = 0;

    //  Find the BADall interval that corresponds to this one.  This BAD interval must be contained
    //  in a BADall (because it contains all bad intervals, while BAD is just the close stuff).
    //  Once we find it, remember the number of reads for later use.

    for (uint32 aa=0; aa<BADall.numberOfIntervals(); aa++)
      if ((BADall.lo(aa) <= BAD.lo(bb)) && (BAD.hi(bb) <= BADall.hi(aa)))
        allHits += BADall.count(aa);

    assert(allHits != 0);

    //  Count the number of reads that span this region.  If the spanning read is not from a library
    //  that might contain subreads, give it more weight.

    for (uint32 ii=0; ii<w->adjLen; ii++)
      if ((w->adj[ii].aovlbgn + 100 < BAD.lo(bb)) && (BAD.hi(bb) + 100 < w->adj[ii].aovlend))
        numSpan += (doCheckSubRead(gkp, w->adj[ii].a_iid)) ? 1 : 2;

    if (subreadFile)
      fprintf(subreadFile, "AcheckSub region %u ("F_S32"-"F_S32") with %u hits %u bighits - span %u largePalindrome %s\n",
              w->adj[0].a_iid, BAD.lo(bb), BAD.hi(bb), BAD.count(bb), allHits,
              numSpan, largePalindrome ? "true" : "false");

    if (numSpan > 9)
      //  If there are 10 or more spanning read (equivalents) this is not a subread junction.  There
      //  is plenty of evidence it is true.
      continue;

    if (BAD.count(bb) + allHits / 4 + largePalindrome < 3)
      //  If 2 or fewer reads claim this is a sub read junction, skip it.  Evidence is weak.
      continue;

    if (subreadFile)
      fprintf(subreadFile, "CONFIRMED BAD REGION %d-%d\n", BAD.lo(bb), BAD.hi(bb));

    w->blist.push_back(badRegion(w->id, badType_subread, BAD.lo(bb), BAD.hi(bb)));
  }
}
