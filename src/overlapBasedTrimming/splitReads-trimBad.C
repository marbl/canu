
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
 *    Brian P. Walenz beginning on 2016-JAN-19
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "splitReads.H"

//  If after chimer trimming a read had a bad interval in the clear, just delete the read.
//  Evidence said it was both good and bad.
//
//  But first, if the bad interval just touches the clear, trim it out.

//  Process the list of bad intervals in blist to find a single clear range.
//  Store the result in w->clrBgn and w->clrEnd.

void
trimBadInterval(gkStore               *gkp,
                workUnit              *w,
                uint32                 minReadLength,
                FILE                  *subreadFile,
                bool                   subreadFileVerbose) {

  if (w->blist.size() == 0)
    return;

  intervalList<int32>  goodRegions;

  //  Build an interval list of all the bad regions, and invert them into good regions.

  for (uint32 bb=0; bb<w->blist.size(); bb++)
    goodRegions.add(w->blist[bb].bgn, w->blist[bb].end - w->blist[bb].bgn);

  goodRegions.invert(w->clrBgn, w->clrEnd);

  //  Find the largest good region, save it in the output clear range.  If there are no
  //  regions (the whole read was marked bad?), default to a bougs clear range.

  w->clrBgn = UINT32_MAX;
  w->clrEnd = UINT32_MAX;

  for (uint32 rr=0; rr<goodRegions.numberOfIntervals(); rr++) {
    if ((w->clrEnd - w->clrBgn) < (goodRegions.hi(rr) - goodRegions.lo(rr))) {
      w->clrBgn = goodRegions.lo(rr);
      w->clrEnd = goodRegions.hi(rr);
    }
  }

  //  If the largest isn't big enough, remember that, and annotate the log appropriately

  if (w->clrEnd - w->clrBgn < minReadLength)
    w->isOK = false;

  //  For logging, find the two bordering bad regions

  if (subreadFile) {
    vector<uint32>   loBad;
    vector<uint32>   hiBad;

    uint32   spur5 = UINT32_MAX;
    uint32   spur3 = UINT32_MAX;

    for (uint32 rr=0; rr<w->blist.size(); rr++) {
      if (w->blist[rr].type == badType_5spur) {
        assert(spur5 == UINT32_MAX);  //  Should be at most one 5' spur region
        spur5 = rr;
      }

      if (w->blist[rr].type == badType_3spur) {
        assert(spur3 == UINT32_MAX);  //  Should be at most one 3' spur region
        spur3 = rr;
      }

      if (w->blist[rr].end == w->clrBgn)
        loBad.push_back(rr);

      if (w->clrEnd == w->blist[rr].bgn)
        hiBad.push_back(rr);
    }

    //  Can't really say much about the number of regions we find.  There could be none (if the
    //  biggest good region extends to the end of the read), exactly one (what we expect), or more
    //  than one (if two algorithms agree on a bad region).

    char  *logPtr = w->logMsg;

    sprintf(logPtr, "iid %6u trim %7u %7u",
            w->id, w->clrBgn, w->clrEnd);
    while (*logPtr)
      logPtr++;

    if (w->isOK == false)
      sprintf(logPtr, " TOO_SHORT");
    while (*logPtr)
      logPtr++;

    if (spur5 != UINT32_MAX)
      sprintf(logPtr, " (5'spur %7u %7u)", w->blist[spur5].bgn, w->blist[spur5].end);
    while (*logPtr)
      logPtr++;

    if (spur3 != UINT32_MAX)
      sprintf(logPtr, " (3'spur %7u %7u)", w->blist[spur3].bgn, w->blist[spur3].end);
    while (*logPtr)
      logPtr++;

    for (uint32 xx=0; xx<loBad.size(); xx++) {
      uint32 x = loBad[xx];

      sprintf(logPtr, " (%s %7u %7u)", w->blist[x].typeName(), w->blist[x].bgn, w->blist[x].end);
      while (*logPtr)
        logPtr++;
    }

    for (uint32 xx=0; xx<hiBad.size(); xx++) {
      uint32 x = hiBad[xx];

      sprintf(logPtr, " (%s %7u %7u)", w->blist[x].typeName(), w->blist[x].bgn, w->blist[x].end);
      while (*logPtr)
        logPtr++;
    }

    logPtr[0] = '\n';
    logPtr[1] = 0;
  }
}

