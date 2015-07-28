
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static char *rcsid = "$Id$";

#include "abAbacus.H"

#include <vector>

using namespace std;


void
abAbacus::baseCallMajority(abColID cid) {
  uint32 bsSum[CNS_NUM_SYMBOLS] = {0};  //  Number of times we've seen this base
  uint32 qvSum[CNS_NUM_SYMBOLS] = {0};  //  Sum of their QVs

  abColumn *column = getColumn(cid);
  abBead   *call   = getBead(column->callID());

  abColBeadIterator *cbi = createColBeadIterator(cid);

  for (abBeadID bid=cbi->next(); bid.isValid(); bid=cbi->next()) {
    abBead *bead = getBead(bid);
    char  bs     = getBase(bead->baseIdx());
    char  qv     = getQual(bead->baseIdx());

    bsSum[baseToIndex[bs]] += 1;
    qvSum[baseToIndex[bs]] += qv;
  }

  delete cbi;

  //  Find the best, ignore ties.

  uint32 bestIdx  = 0;

  for (uint32 i=0; i<CNS_NUM_SYMBOLS; i++)
    if (((bsSum[i] >  bsSum[bestIdx])) ||
        ((bsSum[i] >= bsSum[bestIdx]) && (qvSum[i] > qvSum[bestIdx])))
      bestIdx  = i;

  //  Original version set QV to zero.

  char  base = indexToBase[bestIdx];
  char  qv   = '0';

  setBase(call->baseIdx(), base);
  setQual(call->baseIdx(), qv);
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
abAbacus::baseCallQuality(abColID cid) {
  char    consensusBase = '-';
  char    consensusQV   = '0';

  vector<abBead *>  bReads;  uint32  bBaseCount[CNS_NUM_SYMBOLS] = {0};  uint32  bQVSum[CNS_NUM_SYMBOLS] = {0};  //  Best allele
  vector<abBead *>  oReads;  uint32  oBaseCount[CNS_NUM_SYMBOLS] = {0};  uint32  oQVSum[CNS_NUM_SYMBOLS] = {0};  //  Other allele
  vector<abBead *>  gReads;  uint32  gBaseCount[CNS_NUM_SYMBOLS] = {0};  uint32  gQVSum[CNS_NUM_SYMBOLS] = {0};  //  Guide allele

#if 0
  uint32 highest1_qv[CNS_NUM_SYMBOLS] = {0};
  uint32 highest2_qv[CNS_NUM_SYMBOLS] = {0};
#endif

  uint32 frag_cov = 0;

  bool   used_surrogate = false;

  abColumn *column  = getColumn(cid);
  abBead   *call    = getBead(column->callID());

  uint32  position  = getColumn(call->colIdx())->position();

  abColBeadIterator *cbi = createColBeadIterator(cid);

  // Scan a column of aligned bases.  Sort into three groups:
  //  - those corresponding to the reads of the best allele,
  //  - those corresponding to the reads of the other allele and
  //  - those corresponding to non-read fragments (aka guides)

  //fprintf(stderr, "POSITION %d - ", position);

  for (abBeadID bid=cbi->next(); bid.isValid(); bid=cbi->next()) {
    abBead *bead    = getBead(bid);
    char    base    = getBase(bead->baseIdx());
    int32   baseIdx = baseToIndex[base];
    int     qv      = getQual(bead->baseIdx()) - '0';

    //  Not a base call?  Skip it.
    if (base == 'N')
      continue;

    //  Not a read?  Save the guide base.
    if (getSequence(bead->seqIdx())->isRead() == false) {
      gBaseCount[baseIdx] += 1;
      gQVSum[baseIdx]     += qv;

      gReads.push_back(bead);

      //fprintf(stderr, " guide/%c/%d", base, qv);

      continue;
    }

    //  Otherwise, a real read, with a real base.

    frag_cov++;

#if 0
    //  Find the allele for this iid.  It searched a map of readIID to variant id.
    uint32  iid = getSequence(bead->seqIdx())->gkpIdent();
#endif

    //  If the allele is the target, or we're using all alleles, save to the 'best' list.
    //  Currently, this is the only case.

    if (1) {
      // Save the base in the majority allele.
      bBaseCount[baseIdx]++;
      bQVSum[baseIdx] += qv;
      bReads.push_back(bead);

      //fprintf(stderr, " base/%c/%d", base, qv);
    }

    else {
      //  Save the base in the 'other' allele.
      oBaseCount[baseIdx]++;
      oQVSum[baseIdx] += qv;
      oReads.push_back(bead);

      //fprintf(stderr, " other/%c/%d", base, qv);
    }

    //  Remember the two highest QVs

#if 0
    if (highest1_qv[baseIdx] < qv) {
      highest2_qv[baseIdx] = highest1_qv[baseIdx];
      highest1_qv[baseIdx] = qv;

    } else if (highest2_qv[baseIdx] < qv) {
      highest2_qv[baseIdx] = qv;
    }
#endif
  }

  //fprintf(stderr, "\n");

  delete cbi;

  double  cw[5]    = { 0.0, 0.0, 0.0, 0.0, 0.0 };      // "consensus weight" for a given base
  double  tau[5]   = { 1.0, 1.0, 1.0, 1.0, 1.0 };

  //  Compute tau based on guides

  //fprintf(stderr, "TAU1     %f %f %f %f %f\n", tau[0], tau[1], tau[2], tau[3], tau[4]);

  for (uint32 cind = 0; cind < gReads.size(); cind++) {
    abBead *gb   = gReads[cind];
    char    base = getBase(gb->baseIdx());
    uint32  qv   = getQual(gb->baseIdx()) - '0';

    used_surrogate = true;

    if (qv == 0)
      qv = 5;    /// HUH?!!

    tau[0] += (base == '-') ? PROB[qv] : EPROB[qv];
    tau[1] += (base == 'A') ? PROB[qv] : EPROB[qv];
    tau[2] += (base == 'C') ? PROB[qv] : EPROB[qv];
    tau[3] += (base == 'G') ? PROB[qv] : EPROB[qv];
    tau[4] += (base == 'T') ? PROB[qv] : EPROB[qv];

    //fprintf(stderr, "TAU2[%2d] %f %f %f %f %f qv %d\n", cind, tau[0], tau[1], tau[2], tau[3], tau[4], qv);

    consensusQV = qv;
  }

  //fprintf(stderr, "TAU3     %f %f %f %f %f\n", tau[0], tau[1], tau[2], tau[3], tau[4]);

  //  If other reads exist, reset.

  if (oReads.size() > 0)
    tau[0] = tau[1] = tau[2] = tau[3] = tau[4] = 1.0;

  //fprintf(stderr, "TAU4     %f %f %f %f %f\n", tau[0], tau[1], tau[2], tau[3], tau[4]);

  //  Compute tau based on others

  for (uint32 cind=0; cind < oReads.size(); cind++) {
    abBead  *gb   = oReads[cind];
    char     base = getBase(gb->baseIdx());
    int32    qv   = getQual(gb->baseIdx()) - '0';

    used_surrogate = false;

    if (qv == 0)
      qv = 5;

    tau[0] += (base == '-') ? PROB[qv] : EPROB[qv];
    tau[1] += (base == 'A') ? PROB[qv] : EPROB[qv];
    tau[2] += (base == 'C') ? PROB[qv] : EPROB[qv];
    tau[3] += (base == 'G') ? PROB[qv] : EPROB[qv];
    tau[4] += (base == 'T') ? PROB[qv] : EPROB[qv];

    //fprintf(stderr, "TAU5[%2d] %f %f %f %f %f qv %d\n", cind, tau[0], tau[1], tau[2], tau[3], tau[4], qv);

    consensusQV = qv;
  }

  //  If real reads exist, reset.

  if (bReads.size() > 0)
    tau[0] = tau[1] = tau[2] = tau[3] = tau[4] = 1.0;

  //fprintf(stderr, "TAU6     %f %f %f %f %f\n", tau[0], tau[1], tau[2], tau[3], tau[4]);

  //  Compute tau based on real reads.

  for (uint32 cind=0; cind < bReads.size(); cind++) {
    abBead  *gb   = bReads[cind];
    char     base = getBase(gb->baseIdx());
    int32    qv   = getQual(gb->baseIdx()) - '0';

    used_surrogate = false;

    if (qv == 0)
      qv += 5;

    tau[0] += (base == '-') ? PROB[qv] : EPROB[qv];
    tau[1] += (base == 'A') ? PROB[qv] : EPROB[qv];
    tau[2] += (base == 'C') ? PROB[qv] : EPROB[qv];
    tau[3] += (base == 'G') ? PROB[qv] : EPROB[qv];
    tau[4] += (base == 'T') ? PROB[qv] : EPROB[qv];

    //fprintf(stderr, "TAU7[%2d] %f %f %f %f %f qv %d\n", cind, tau[0], tau[1], tau[2], tau[3], tau[4], qv);

    consensusQV = qv;
  }

  //fprintf(stderr, "TAU8     %f %f %f %f %f\n", tau[0], tau[1], tau[2], tau[3], tau[4]);

  //  Occasionally we get a single read of coverage, and the base is an N, which we ignored above.

  if ((bReads.size() == 0) &&
      (oReads.size() == 0) &&
      (gReads.size() == 0)) {
    setBase(call->baseIdx(), 'N');
    setQual(call->baseIdx(), '0');

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
    consensusQV = CNS_MAX_QV + '0';
  }

  //  Otherwise compute the QV.  If there is more than one read, or we used the surrogate,
  //  we can compute from cwMax.  If only one read, use its qv.
  //
  else {
    int32  qv = consensusQV;

    if ((frag_cov > 1) || (used_surrogate == true)) {
      double dqv =  -10.0 * log10(1.0 - cwMax);

      qv = DBL_TO_INT(dqv);

      if (dqv - qv >= 0.50)
        qv++;
    }

    qv = MIN(CNS_MAX_QV, qv);
    qv = MAX(CNS_MIN_QV, qv);

    consensusQV = qv + '0';
  }


  //  If no target allele, or this is the target allele, set the base.  Since there
  //  is (currently) always no target allele, we always set the base.

  //fprintf(stderr, "SET %u to %c qv %c\n", position, consensusBase, consensusQV);

  setBase(call->baseIdx(), consensusBase);
  setQual(call->baseIdx(), consensusQV);
}



char
abAbacus::baseCall(abColID      cid,
                   bool         highQuality) {

  if (highQuality)
    baseCallQuality(cid);
  else
    baseCallMajority(cid);

  abColumn *column = getColumn(cid);
  abBead   *call   = getBead(column->callID());

  return(getBase(call->baseIdx()));
}


