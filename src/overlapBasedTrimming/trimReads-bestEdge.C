
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

#include <vector>
#include <functional>
#include <algorithm>

using namespace std;

//  Generates plots for each trim.
#undef  GNUPLOT

//  Logging to the screen.
#undef  VERBOSE


bool
bestEdge(ovOverlap  *ovl,
         uint32       ovlLen,
         uint32       readID,
         uint32       readLen,
         uint32       ibgn,
         uint32       iend,
         uint32      &fbgn,
         uint32      &fend,
         char        *logMsg,
         uint32       errorValue,
         uint32       minOverlap,
         uint32       minCoverage,
         uint32       minReadLength) {

  fbgn      = ibgn;
  fend      = iend;
  logMsg[0] = 0;

  assert(readID == ovl[0].a_iid);
  assert(ovlLen > 0);

  //
  //  Trim once, using the largest covered rule.  This gets rid of any gaps in overlap coverage,
  //  etc.
  //

  uint32  lbgn = 0;
  uint32  lend = 0;

  if (largestCovered(ovl, ovlLen, readID, readLen, ibgn, iend, lbgn, lend, logMsg, errorValue, minOverlap, minCoverage, minReadLength) == false)
    return(false);

  //
  //  Trim again, to maximize overlap length.
  //

  vector<uint32>    trim5;
  vector<uint32>    trim3;

  uint32            nContained = 0;

#ifdef GNUPLOT
  vector<uint32>    trim5iid, trim5sco;
  vector<uint32>    trim3iid, trim3sco;
#endif

  //  For each overlap, add potential trim points where the overlap ends.

  for (uint32 i=0; i<ovlLen; i++) {
    uint32 tbgn = ovl[i].a_bgn();
    uint32 tend = ovl[i].a_end();

    if ((lend <= tbgn) ||
        (tend <= lbgn))
      //  Doesn't intersect the largest covered region.
      continue;

    assert(tbgn < tend);

    assert(ibgn <= tbgn);
    assert(tend <= iend);

    assert(readID == ovl[i].a_iid);

    if (ovl[i].evalue() > errorValue)
      //  Overlap is crappy.
      continue;

    //  Add trim points as long as they aren't outside the largest covered region.

    if (lbgn <= tbgn)
      trim5.push_back(tbgn);

    if (tend <= lend)
      trim3.push_back(tend);

    //  If overlap indicates this read is contained, we're done.  There isn't any
    //  trimming needed.  We can set the final trim and return, or let the rest
    //  of the algorithm/logging finish.

    if ((tbgn == ibgn) &&
        (tend == iend))
      nContained++;
  }

  //  If the read is contained in more than one other read, we're done.  No trimming needed.

  if (nContained >= 2) {
    //fbgn = ibgn;
    //fend = iend;
    //return(true);

    trim5.clear();
    trim3.clear();
  }

  //  Add trim points for the largest covered, i.e., no trimming past what largest covered did.

  trim5.push_back(lbgn);
  trim3.push_back(lend);

  //  Duplicate removal, and possibly the processing algorithm, need the trim points
  //  sorted from outside-in.

  sort(trim5.begin(), trim5.end(), std::less<uint32>());
  sort(trim3.begin(), trim3.end(), std::greater<uint32>());

  //  Remove duplicate points (easier here than in the loops below)

  {
    uint32 old = 0;
    uint32 sav = 0;

    for (old=0, sav=0; old<trim5.size(); old++)
      if (trim5[old] != trim5[sav])
        trim5[++sav] = trim5[old];
    trim5.resize(sav+1);

#ifdef GNUPLOT
    trim5iid.resize(sav+1);
    trim5sco.resize(sav+1);
#endif

    for (old=0, sav=0; old<trim3.size(); old++)
      if (trim3[old] != trim3[sav])
        trim3[++sav] = trim3[old];
    trim3.resize(sav+1);

#ifdef GNUPLOT
    trim3iid.resize(sav+1);
    trim3sco.resize(sav+1);
#endif
  }

  uint32   best5pt    = 0;  //  Index into trim5
  uint32   best5score = 0;  //  Score for the best index
  uint32   best5iid   = 0;  //  IID for the best index
  uint32   best5end   = 0;  //  Coord of the end point for the best index

  uint32   best3pt    = 0;
  uint32   best3score = 0;
  uint32   best3iid   = 0;
  uint32   best3bgn   = UINT32_MAX;

  //  Find the best 5' point.

  for (uint32 pt=0; pt < trim5.size(); pt++) {
    uint32   triml    = trim5[pt];
    uint32   score    = 0;
    uint32   sciid    = 0;
    uint32   scend    = 0;

    if (best5score >= readLen - triml)
      //  Not possible to get a higher score by trimming more.
      break;

    //fprintf(stderr, "trim5 pt %u out of %u\n", pt, trim5.size());

    for (uint32 i=0; i < ovlLen; i++) {
      uint32 tbgn = ibgn + ovl[i].a_bgn();
      uint32 tend = ibgn + ovl[i].a_end();

      if ((triml <  tbgn) ||
          (tend  <= triml))
        //  Alignment starts after of the trim point; not a valid overlap.
        //  or, alignment ends before the trim point (trimmed out).
        continue;

      //  Limit the overlap to the largest covered.
      if (tend > lend)
        tend = lend;

      assert(tend >= triml);

      uint32 tlen = tend - triml;

      assert(tlen <= readLen);

      //  Save the best score.  Break ties by favoring the current best.

      if (((score  < tlen)) ||
          ((score == tlen) && (ovl[i].b_iid == best5iid))) {
        score = tlen;
        sciid = ovl[i].b_iid;
        scend = tend;
#ifdef GNUPLOT
        trim5sco[pt] = score;
        trim5iid[pt] = sciid;
#endif
      }
    }

    //  Give up if we're not finding anything longer after 1/3 of the current read.  Previous
    //  attempts at this wanted to ensure that the current read (sciid) is the same as the best
    //  (best5iid) but ties and slightly different begin/end points make this unreliable.  For
    //  example, an overlap from 100-500 and an overlap from 200-501.  The first is longest up until
    //  trim point 200, then the second becomes longest.  The BEST longest is still the first
    //  overlap, at trim point 100.

    if (score < best5score * 0.66)
      break;

    //  Save a new best if the score is better, and the endpoint is more interior to the read.

    if ((best5score < score) &&
        (best5end   < scend)) {
      //fprintf(stderr, "RESET 5 end to pt %d score %d iid %d end %d\n",
      //        pt, score, sciid, scend);
      best5pt    = pt;
      best5score = score;
      best5iid   = sciid;
      best5end   = scend;
    }
  }

  //fprintf(stderr, "BEST at %u position %u pt %u\n", best5score, trim5[best5pt], best5pt);

  //  Find the best 3' point.

  for (uint32 pt=0; pt<trim3.size(); pt++) {
    uint32   trimr    = trim3[pt];;
    uint32   score    = 0;
    uint32   sciid    = 0;
    uint32   scbgn    = 0;

    if (best3score >= trimr - 0)
      //  Not possible to get a higher score by trimming more.
      break;

    //fprintf(stderr, "trim3 pt %u out of %u\n", pt, trim3.size());

    for (uint32 i=0; i < ovlLen; i++) {
      uint32 tbgn = ibgn + ovl[i].a_bgn();
      uint32 tend = ibgn + ovl[i].a_end();

      if ((tend < trimr) ||
          (trimr <= tbgn))
        //  Alignment ends before the trim point; not a valid overlap,
        //  or, alignment starts after the trim point (trimmed out)
        continue;

      //  Limit the overlap to the largest covered.
      if (tbgn < lbgn)
        tbgn = lbgn;

      assert(trimr >= tbgn);

      uint32 tlen = trimr - tbgn;

      assert(tlen <= readLen);

      if (((score  < tlen)) ||
          ((score == tlen) && (ovl[i].b_iid == best3iid))) {
        score = tlen;
        sciid = ovl[i].b_iid;
        scbgn = tbgn;
#ifdef GNUPLOT
        trim3sco[pt] = score;
        trim3iid[pt] = sciid;
#endif
      }
    }

    if (score < best3score * 0.66)
      break;

    if ((best3score < score) &&
        (scbgn      < best3bgn)) {
      //fprintf(stderr, "RESET 3 end to pt %d score %d iid %d bgn %d\n",
      //        pt, score, sciid, scbgn);
      best3pt    = pt;
      best3score = score;
      best3iid   = sciid;
      best3bgn   = scbgn;
    }
  }

#ifdef GNUPLOT
  {
    char  D[FILENAME_MAX];
    char  G[FILENAME_MAX];
    char  S[FILENAME_MAX];

    FILE *F;

    snprintf(D, FILENAME_MAX, "trim-%08d.dat", readID);
    snprintf(G, FILENAME_MAX, "trim-%08d.gp",  readID);
    snprintf(S, FILENAME_MAX, "gnuplot < trim-%08d.gp", readID);

    F = fopen(D, "w");
    for (uint32 i=0; i<MAX(trim5.size(), trim3.size()); i++) {
      if      (i < trim5.size() && i < trim3.size())
        fprintf(F, "trim[%03d] pt %4d len %4d iid %7d -- pt %4d len %4d iid %7d \n",
                i,
                trim5[i], trim5sco[i], trim5iid[i],
                trim3[i], trim3sco[i], trim3iid[i]);
      else if (i < trim5.size())
        fprintf(F, "trim[%03d] pt %4d len %4d iid %7d -- pt %4d len %4d iid %7d \n",
                i,
                trim5[i], trim5sco[i], trim5iid[i],
                0, 0, 0);
      else
        fprintf(F, "trim[%03d] pt %4d len %4d iid %7d -- pt %4d len %4d iid %7d \n",
                i,
                0, 0, 0,
                trim3[i], trim3sco[i], trim3iid[i]);
    }
    AS_UTL_closeFile(F, D);


    F = fopen(G, "w");
    fprintf(F, "set terminal png\n");
    fprintf(F, "set output \"trim-%08d.png\"\n",
            readID);
    fprintf(F, "plot \"trim-%08d.dat\" using 3:5 with linespoints, \"trim-%08d.dat\" using 10:12 with linespoints\n",
            readID,
            readID);
    AS_UTL_closeFile(F, G);


    system(S);
  }
#endif

  //fprintf(stderr, "BEST at %u position %u pt %u\n", best3score, trim3[best3pt], best3pt);

  //  Set trimming.  Be just a little aggressive, and get rid of an extra base or two, if possible.
  //  If not possible, the read is crap, and will be deleted anyway.

  fbgn = trim5[best5pt];
  fend = trim3[best3pt];

  if ((fbgn + 4 < fend) &&
      (2        < fend)) {
    fbgn += 2;
    fend -= 2;
  }

  //  Did we go outside the largestCovered region?  Just reset.  Ideally we'd assert, and crash.

  if (fbgn < lbgn)  fbgn = lbgn;
  if (lend < fend)  fend = lend;

  //  Lastly, did we just end up with a bogus trim?  There isn't any guard against picking the 5'
  //  trim point to the right of the 3' trim point.

#if 1
  if (fend < fbgn) {
    fprintf(stderr, "iid = %u\n", readID);
    fbgn = lbgn;
    fend = lend;
  }
#endif

  assert(fbgn <= fend);

  return(true);
}
