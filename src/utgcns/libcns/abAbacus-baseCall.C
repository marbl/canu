
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
 *    src/AS_CNS/BaseCall.C
 *    src/AS_CNS/BaseCall.c
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/BaseCall.C
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
 *    Brian P. Walenz from 2005-JUN-16 to 2013-AUG-01
 *      are Copyright 2005-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2005-SEP-29 to 2006-OCT-03
 *      are Copyright 2005-2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-FEB-27 to 2009-MAY-14
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUL-28
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "abAbacus.H"

#include <vector>

using namespace std;


void
abColumn::baseCallMajority(void) {
  uint32 bsSum[CNS_NUM_SYMBOLS] = {0};  //  Number of times we've seen this base
  uint32 qvSum[CNS_NUM_SYMBOLS] = {0};  //  Sum of their QVs

  for (uint32 ii=0; ii<_beadsLen; ii++) {
    char    base  = _beads[ii].base();
    uint8   qual  = _beads[ii].qual();
    uint32  bidx  = baseToIndex[base];

    if (bidx >= CNS_NUM_SYMBOLS)
      fprintf(stderr, "abColumn::baseCallMajority()--  For column %u, link %u, base '%c' (%d) is invalid.\n",
              position(), ii, _beads[ii].base(), _beads[ii].base());
    assert(bidx < CNS_NUM_SYMBOLS);

    bsSum[bidx] += 1;
    qvSum[bidx] += qual;
  }

  //  Find the best, and second best, ignore ties.

  uint32 bestIdx = 0;
  uint32 nextIdx = 1;

  for (uint32 i=1; i<CNS_NUM_SYMBOLS; i++) {
    if        (((bsSum[i] >  bsSum[bestIdx])) ||
               ((bsSum[i] >= bsSum[bestIdx]) && (qvSum[i] > qvSum[bestIdx]))) {
      nextIdx = bestIdx;
      bestIdx = i;
    } else if (((bsSum[i] >  bsSum[nextIdx])) ||
               ((bsSum[i] >= bsSum[nextIdx]) && (qvSum[i] > qvSum[nextIdx]))) {
      nextIdx = i;
    }
  }

  //  Original version set QV to zero.

  _call = indexToBase[bestIdx];
  _qual = 0;

  //  If the best is a gap, use the lowercase second best - the alignment will trest this specially

  if (_call == '-')
    _call = indexToBase[nextIdx] - 'A' + 'a';
}



uint32
baseToTauIndex(char base) {

  switch (base) {
  case '-':
    return(0);
    break;
  case 'A':
    return(1);
    break;
  case 'C':
    return(2);
    break;
  case 'G':
    return(3);
    break;
  case 'T':
    return(4);
    break;
  default:
    fprintf(stderr, "baseToTauIndex()-- invalid base '%c' %d\n", base, base);
    assert(0);
    break;
  }

  return(0);
}



void
abColumn::baseCallQuality(void) {
  char    consensusBase = '-';
  char    consensusQual = 0;

  //  The original versions classified reads as 'best allele', 'other allele' or 'guide allele'.
  //  Other allele was set if we were targetting a specific set of reads here (e.g., we had already
  //  clustered reads into haplotypes).  Guide allele was used if the read was not really a read
  //  (originally, these were reads that came from outside the project, but eventually, I think it
  //  came to also mean the read was a unitig surrogate).  All this was stripped out in early
  //  December 2015.  The pieces removed all mirrored what is done for the bReads.

  vector<uint16>  bReads;   uint32  bBaseCount[CNS_NUM_SYMBOLS] = {0};  uint32  bQVSum[CNS_NUM_SYMBOLS] = {0};  //  Best allele

  uint32 frag_cov = 0;

  // Scan a column of aligned bases.  Sort into three groups:
  //  - those corresponding to the reads of the best allele,
  //  - those corresponding to the reads of the other allele and
  //  - those corresponding to non-read fragments (aka guides)

  for (uint32 ii=0; ii<_beadsLen; ii++) {
    char    base = _beads[ii].base();
    uint8   qual = _beads[ii].qual();
    uint32  bidx = baseToIndex[base];  //  If we encode bases properly, this will go away.

    if (bidx >= CNS_NUM_SYMBOLS)
      fprintf(stderr, "abColumn::baseCallQuality()--  For column %u, link %u, base '%c' (%d) is invalid.\n",
              position(), ii, _beads[ii].base(), _beads[ii].base());
    assert(bidx < CNS_NUM_SYMBOLS);

    bBaseCount[bidx] += 1;   //  Could have saved to 'best', 'other' or 'guide' here.
    bQVSum[bidx]     += qual;

    bReads.push_back(ii);
  }

  double  cw[5]    = { 0.0, 0.0, 0.0, 0.0, 0.0 };      // "consensus weight" for a given base
  double  tau[5]   = { 1.0, 1.0, 1.0, 1.0, 1.0 };

  //  Compute tau based on real reads.

  for (uint32 cind=0; cind < bReads.size(); cind++) {
    char     base = _beads[ bReads[cind] ].base();
    uint8    qual = _beads[ bReads[cind] ].qual();

    if (qual == 0)
      qual += 5;

    tau[0] += (base == '-') ? PROB[qual] : EPROB[qual];
    tau[1] += (base == 'A') ? PROB[qual] : EPROB[qual];
    tau[2] += (base == 'C') ? PROB[qual] : EPROB[qual];
    tau[3] += (base == 'G') ? PROB[qual] : EPROB[qual];
    tau[4] += (base == 'T') ? PROB[qual] : EPROB[qual];

    //fprintf(stderr, "TAU7[%2d] %f %f %f %f %f qv %d\n", cind, tau[0], tau[1], tau[2], tau[3], tau[4], qual);

    consensusQual = qual;
  }

  //fprintf(stderr, "TAU8     %f %f %f %f %f\n", tau[0], tau[1], tau[2], tau[3], tau[4]);

  //  Occasionally we get a single read of coverage, and the base is an N, which we ignored above.
  //  This is probably of historical interest any more (it happened with 454 reads) but is left in
  //  because its a cheap and working and safe.

  if (bReads.size() == 0) {  //  + oReads.size() + gReads.size()
    _call = 'N';
    _qual = 0;
    return;
  }

  //  Do we need to scale before normalizing?  Anything out of bounds?  We'll try to scale the small
  //  values up to DBL_MIN without making the large values larger than DBL_MAX.  If we end up with
  //  some values still too small, oh well.  We have a winner (the large value) anyway!
  //
  //  Note, in log-space, min value (1e-309) is around -711.5, and max value (1e309) is around 711.5.

  //fprintf(stderr, "tau: %f %f %f %f %f\n", tau[0], tau[1], tau[2], tau[3], tau[4]);

  double  minTau =  DBL_MAX;
  double  maxTau = -DBL_MAX;

  minTau = MIN(minTau, tau[0]);
  minTau = MIN(minTau, tau[1]);
  minTau = MIN(minTau, tau[2]);
  minTau = MIN(minTau, tau[3]);
  minTau = MIN(minTau, tau[4]);

  maxTau = MAX(maxTau, tau[0]);
  maxTau = MAX(maxTau, tau[1]);
  maxTau = MAX(maxTau, tau[2]);
  maxTau = MAX(maxTau, tau[3]);
  maxTau = MAX(maxTau, tau[4]);

  //  Now that we know the min and max values, shift them as far positive as possible.  Ideally,
  //  this will just add an offset to bring the smallest value up to the minimum representable.

  double  minValue   = log(DBL_MIN) + DBL_EPSILON;
  double  maxValue   = log(DBL_MAX) - DBL_EPSILON;

  double  scaleValue = 0.0;

  if (minTau < minValue)
    scaleValue = minValue - minTau;

  if (maxTau + scaleValue > maxValue)
    scaleValue = maxValue - maxTau;

  tau[0] += scaleValue;
  tau[1] += scaleValue;
  tau[2] += scaleValue;
  tau[3] += scaleValue;
  tau[4] += scaleValue;

  //  It could however overflow the max (which is the value we care about), so any values still too
  //  low are thresholded.

  if (tau[0] < minValue)  tau[0] = minValue;
  if (tau[1] < minValue)  tau[1] = minValue;
  if (tau[2] < minValue)  tau[2] = minValue;
  if (tau[3] < minValue)  tau[3] = minValue;
  if (tau[4] < minValue)  tau[4] = minValue;

  //fprintf(stderr, "TAU9     %f %f %f %f %f value %f/%f tau %f/%f scale %f\n",
  //        tau[0], tau[1], tau[2], tau[3], tau[4], minValue, maxValue, minTau, maxTau, scaleValue);

  //  I give up.  I can't make the following asserts true when the scale value is reset to not
  //  exceed maxValue.  I'm off, somewhere, by 1e-13 (EPSILON is 2e-16, nowhere near).  My test case
  //  takes 18 minutes to get here, and I don't really want to slap in a bunch of logging.  So,
  //  thresholding it is.

  if (tau[0] > maxValue)  tau[0] = maxValue;
  if (tau[1] > maxValue)  tau[1] = maxValue;
  if (tau[2] > maxValue)  tau[2] = maxValue;
  if (tau[3] > maxValue)  tau[3] = maxValue;
  if (tau[4] > maxValue)  tau[4] = maxValue;

  assert(tau[0] <= maxValue);
  assert(tau[1] <= maxValue);
  assert(tau[2] <= maxValue);
  assert(tau[3] <= maxValue);
  assert(tau[4] <= maxValue);

  tau[0] = exp(tau[0]);
  tau[1] = exp(tau[1]);
  tau[2] = exp(tau[2]);
  tau[3] = exp(tau[3]);
  tau[4] = exp(tau[4]);

  //fprintf(stderr, "TAU10    %f %f %f %f %f\n", tau[0], tau[1], tau[2], tau[3], tau[4]);

  assert(tau[0] >= 0.0);
  assert(tau[1] >= 0.0);
  assert(tau[2] >= 0.0);
  assert(tau[3] >= 0.0);
  assert(tau[4] >= 0.0);

  cw[0] = tau[0] * 0.2;
  cw[1] = tau[1] * 0.2;
  cw[2] = tau[2] * 0.2;
  cw[3] = tau[3] * 0.2;
  cw[4] = tau[4] * 0.2;

  double  normalize  = cw[0] + cw[1] + cw[2] + cw[3] + cw[4];
  double  normalizeS = normalize;

  if (normalize > 0.0) {
    normalize = 1.0 / normalize;

    cw[0] *= normalize;
    cw[1] *= normalize;
    cw[2] *= normalize;
    cw[3] *= normalize;
    cw[4] *= normalize;
  }

  if ((cw[0] < 0.0) || (cw[1] < 0.0) || (cw[2] < 0.0) || (cw[3] < 0.0) || (cw[4] < 0.0))
    fprintf(stderr, "ERROR: cw[0-4] invalid: %f %f %f %f %f - tau %f %f %f %f %f - normalize %.60e %.60e\n",
            cw[0], cw[1], cw[2], cw[3], cw[4],
            tau[0], tau[1], tau[2], tau[3], tau[4],
            normalize, normalizeS);
  if (!(cw[0] >= 0.0) || !(cw[1] >= 0.0) || !(cw[2] >= 0.0) || !(cw[3] >= 0.0) || !(cw[4] >= 0.0))
    fprintf(stderr, "ERROR: cw[0-4] invalid: %f %f %f %f %f - tau %f %f %f %f %f - normalize %.60e %.60e\n",
            cw[0], cw[1], cw[2], cw[3], cw[4],
            tau[0], tau[1], tau[2], tau[3], tau[4],
            normalize, normalizeS);
  assert(cw[0] >= 0.0);
  assert(cw[1] >= 0.0);
  assert(cw[2] >= 0.0);
  assert(cw[3] >= 0.0);
  assert(cw[4] >= 0.0);

  double cwMax = DBL_MIN;
  uint32 cwIdx = 0;     //  Default is a gap.

  if (cw[0] > cwMax) { cwIdx = 0;  cwMax = cw[0]; }  //  '-'
  if (cw[1] > cwMax) { cwIdx = 1;  cwMax = cw[1]; }  //  'A'
  if (cw[2] > cwMax) { cwIdx = 2;  cwMax = cw[2]; }  //  'C'
  if (cw[3] > cwMax) { cwIdx = 3;  cwMax = cw[3]; }  //  'G'
  if (cw[4] > cwMax) { cwIdx = 4;  cwMax = cw[4]; }  //  'T'

  //  If cwMax = 0 then consensus is a gap.  Otherwise, we deterministically set it to A C G T.

#if 0
  fprintf(stderr, "prob('%c') %f %c\n", indexToBase[0], cw[0], (cw[0] == cwMax) ? '*' : ' ');
  fprintf(stderr, "prob('%c') %f %c\n", indexToBase[1], cw[1], (cw[1] == cwMax) ? '*' : ' ');
  fprintf(stderr, "prob('%c') %f %c\n", indexToBase[2], cw[2], (cw[2] == cwMax) ? '*' : ' ');
  fprintf(stderr, "prob('%c') %f %c\n", indexToBase[3], cw[3], (cw[3] == cwMax) ? '*' : ' ');
  fprintf(stderr, "prob('%c') %f %c\n", indexToBase[4], cw[4], (cw[4] == cwMax) ? '*' : ' ');
#endif

  consensusBase = indexToBase[cwIdx];

  //  Compute the QV.

  //  If cwMax is big, we've max'd out the QV and set it to the maximum.
  //
  if (cwMax >= 1.0 - DBL_EPSILON) {
    consensusQual = CNS_MAX_QV;
  }

  //  Otherwise compute the QV.  If there is more than one read, or we used the surrogate,
  //  we can compute from cwMax.  If only one read, use its qv.
  //
  else {
    int32  qual = consensusQual;

    if ((bReads.size() > 1) /*|| (used_surrogate == true)*/) {
      double dqv =  -10.0 * log10(1.0 - cwMax);

      qual = DBL_TO_INT(dqv);

      if (dqv - qual >= 0.50)
        qual++;
    }

    qual = MIN(CNS_MAX_QV, qual);
    qual = MAX(CNS_MIN_QV, qual);

    consensusQual = qual;
  }

  //  If no target allele, or this is the target allele, set the base.  Since there
  //  is (currently) always no target allele, we always set the base.

  //fprintf(stderr, "SET %u to %c qv %c\n", position, consensusBase, consensusQV);

  _call = consensusBase;
  _qual = consensusQual;
}



char
abColumn::baseCall(bool  highQuality) {

  if (highQuality)
    baseCallQuality();
  else
    baseCallMajority();

  return(_call);
}


