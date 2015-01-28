
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





void
abAbacus::baseCallQuality(abColID cid) {
  char    consensusBase = '-';
  char    consensusQV   = '0';

  vector<abBead *>  bReads;  uint32  bBaseCount[CNS_NUM_SYMBOLS] = {0};  uint32  bQVSum[CNS_NUM_SYMBOLS] = {0};  //  Best allele
  vector<abBead *>  oReads;  uint32  oBaseCount[CNS_NUM_SYMBOLS] = {0};  uint32  oQVSum[CNS_NUM_SYMBOLS] = {0};  //  Other allele
  vector<abBead *>  gReads;  uint32  gBaseCount[CNS_NUM_SYMBOLS] = {0};  uint32  gQVSum[CNS_NUM_SYMBOLS] = {0};  //  Guide allele

  double  cw[5]    = { 0.0, 0.0, 0.0, 0.0, 0.0 };      // "consensus weight" for a given base
  double  tau[5]   = { 1.0, 1.0, 1.0, 1.0, 1.0 };

  uint32 highest1_qv[CNS_NUM_SYMBOLS] = {0};
  uint32 highest2_qv[CNS_NUM_SYMBOLS] = {0};

  uint32 frag_cov = 0;

  bool   used_surrogate = false;

  abColumn *column = getColumn(cid);
  abBead   *call   = getBead(column->callID());

  abColBeadIterator *cbi = createColBeadIterator(cid);

  // Scan a column of aligned bases.  Sort into three groups:
  //  - those corresponding to the reads of the best allele,
  //  - those corresponding to the reads of the other allele and
  //  - those corresponding to non-read fragments (aka guides)

  for (abBeadID bid=cbi->next(); bid.isValid(); bid=cbi->next()) {
    abBead *bead    = getBead(bid);
    char    base    = getBase(bead->baseIdx());
    int32   baseIdx = baseToIndex[base];
    int     qv      = getQual(bead->baseIdx()) - '0';

    if (base == 'N')
      continue;

    if (getSequence(bead->seqIdx())->isRead() == false) {
      gBaseCount[baseIdx] += 1;
      gQVSum[baseIdx]     += qv;

      gReads.push_back(bead);

      continue;
    }

    frag_cov++;

    uint32  iid     = getSequence(bead->seqIdx())->gkpIdent();

    //  Find the allele for this iid.  It searched a map of readIID to variant id.

    //  If the allele is the target, or we're using all alleles, save to the 'best' list.
    //  Currently, this is the only case.

    if (1) {
      bBaseCount[baseIdx]++;
      bQVSum[baseIdx] += qv;
      bReads.push_back(bead);
    }

    else {
      oBaseCount[baseIdx]++;
      oQVSum[baseIdx] += qv;
      oReads.push_back(bead);
    }

    //  Remember the two highest QVs

    if (highest1_qv[baseIdx] < qv) {
      highest2_qv[baseIdx] = highest1_qv[baseIdx];
      highest1_qv[baseIdx] = qv;

    } else if (highest2_qv[baseIdx] < qv) {
      highest2_qv[baseIdx] = qv;
    }
  }

  delete cbi;

  //  Compute tau based on guides

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

    consensusQV = qv;
  }

  //  If actual reads,reset.

  if (oReads.size() > 0)
    tau[0] = tau[1] = tau[2] = tau[3] = tau[4] = 1.0;

  //  Compute tau based on others

  for (uint32 cind=0; cind<oReads.size(); cind++) {
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

    consensusQV = qv;
  }

  //  If real reads, reset.

  if (bReads.size() > 0)
    tau[0] = tau[1] = tau[2] = tau[3] = tau[4] = 1.0;

  //  Compute tau based on real reads.

  for (uint32 cind=0; cind<bReads.size(); cind++) {
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

    consensusQV = qv;
  }

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

  double  scaleValue = 0.0;

  double  minValue   = log(DBL_MIN);
  double  maxValue   = log(DBL_MAX);

  double  minTau = DBL_MAX;
  double  maxTau = DBL_MIN;

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

  if (minTau < minValue)
    scaleValue = minValue - minTau;

  if (maxTau + scaleValue > maxValue)
    scaleValue = maxValue - maxTau;

  tau[0] = exp(tau[0] + scaleValue);
  tau[1] = exp(tau[1] + scaleValue);
  tau[2] = exp(tau[2] + scaleValue);
  tau[3] = exp(tau[3] + scaleValue);
  tau[4] = exp(tau[4] + scaleValue);

  double  normalize = 0.0;

  cw[0] = tau[0] * 0.2;    normalize += cw[0];
  cw[1] = tau[1] * 0.2;    normalize += cw[1];
  cw[2] = tau[2] * 0.2;    normalize += cw[2];
  cw[3] = tau[3] * 0.2;    normalize += cw[3];
  cw[4] = tau[4] * 0.2;    normalize += cw[4];

  if (normalize > 0.0) {
    normalize = 1.0 / normalize;

    cw[0] *= normalize;
    cw[1] *= normalize;
    cw[2] *= normalize;
    cw[3] *= normalize;
    cw[4] *= normalize;
  }

  double cwMax = 0.0;
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

  //  If there isn't a clear winner:
  //
  //    If there is more than one fragment, or we used a surrogate (zero fragments) compute the QV
  //    from cwMax.
  //
  //    Otherwise, there is only one fragment, and we should use that qv (saved in consensusQV).
  //

  if (cwMax < 1.0 - DBL_EPSILON) {
    int32  qv = consensusQV;

    if ((frag_cov != 1) || (used_surrogate == true)) {
      double dqv =  -10.0 * log10(1.0 - cwMax);

      qv = DBL_TO_INT(dqv);

      if (dqv - qv >= 0.50)
        qv++;
    }

    qv = MIN(CNS_MAX_QV, qv);
    qv = MAX(CNS_MIN_QV, qv);

    consensusQV = qv + '0';
  }

  else {
    consensusQV = CNS_MAX_QV + '0';
  }

  //  If no target allele, or this is the target allele, set the base.  Since there
  //  is (currently) always no target allele, we always set the base.

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


