
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

#if 0
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"
//#include "MicroHetREZ.H"
#include "AS_UTL_reverseComplement.H"
#endif

#include <vector>

using namespace std;


void
abAbacus::baseCallMajority(abColID cid) {
  uint32 bsSum[CNS_NP] = {0};  //  Number of times we've seen this base
  uint32 qvSum[CNS_NP] = {0};  //  Sum of their QVs

  abColumn *column = getColumn(cid);
  abBead   *call   = getBead(column->callID());

  abColBeadIterator *cbi = createColBeadIterator(cid);

  for (abBeadID bid=cbi->next(); bid.isValid(); bid=cbi->next()) {
    abBead *bead = getBead(bid);
    char  bs     = getBase(bead->baseIdx());
    char  qv     = getQual(bead->baseIdx());

    bsSum[RINDEX[bs]] += 1;
    qvSum[RINDEX[bs]] += qv;
  }

  //  Find the best, ignore ties.

  uint32 bestIdx  = 0;

  for (uint32 i=0; i<CNS_NALPHABET; i++)
    if (((bsSum[i] >  bsSum[bestIdx])) ||
        ((bsSum[i] >= bsSum[bestIdx]) && (qvSum[i] > qvSum[bestIdx])))
      bestIdx  = i;

  //  Original version set QV to zero.

  char  base = toupper(RALPHABET[bestIdx]);
  char  qv   = '0';

  setBase(call->baseIdx(), base);
  setQual(call->baseIdx(), qv);
}



void
abAbacus::baseCallQuality(abVarRegion   &vreg,
                          abColID        cid,
                          double        &var,
                          int32          target_allele,
                          bool           getScores,
                          int32          split_alleles,
                          int32          smooth_win) {

  char    consensusBase = '-';
  char    consensusQV   = '0';

  vector<abBead *>  bReads;  uint32  bBaseCount[CNS_NP] = {0};  uint32  bQVSum[CNS_NP] = {0};
  vector<abBead *>  oReads;  uint32  oBaseCount[CNS_NP] = {0};  uint32  oQVSum[CNS_NP] = {0};
  vector<abBead *>  gReads;  uint32  gBaseCount[CNS_NP] = {0};  uint32  gQVSum[CNS_NP] = {0};

  double  cw[5]    = { 0.0, 0.0, 0.0, 0.0, 0.0 };      // "consensus weight" for a given base
  double  tau[5]   = { 1.0, 1.0, 1.0, 1.0, 1.0 };

  uint32 highest1_qv[CNS_NP] = {0};
  uint32 highest2_qv[CNS_NP] = {0};

  uint32 frag_cov = 0;

  bool   used_surrogate = false;

  abColumn *column = getColumn(cid);
  abBead   *call   = getBead(column->callID());

  abColBeadIterator *cbi = createColBeadIterator(cid);


  // Scan a column of aligned bases (=beads).
  // Sort the beads into three groups:
  //      - those corresponding to the reads of the best allele,
  //      - those corresponding to the reads of the other allele and
  //      - those corresponding to non-read fragments (aka guides)

  for (abBeadID bid=cbi->next(); bid.isValid(); bid=cbi->next()) {
    abBead *bead    = getBead(bid);
    char    base    = getBase(bead->baseIdx());
    int32   baseIdx = RINDEX[base];
    int     qv      = getQual(bead->baseIdx()) - '0';

    if (base == 'N')
      continue;

    if (getSequence(bead->seqIdx())->isRead() == false) {
      //assert(type == AS_UNITIG);

      gBaseCount[baseIdx] += 1;
      gQVSum[baseIdx]     += qv;

      gReads.push_back(bead);

      continue;
    }

    frag_cov++;

    uint32  iid     = getSequence(bead->seqIdx())->gkpIdent();
    uint32  vregidx = 0;

    assert(vreg.nr >= 0);

    for (; vregidx < vreg.nr; vregidx++)
      if (iid == vreg.iids[vregidx])
        break;

    //  If there are vreg's, assert we found the iid in it
    assert((vreg.nr <= 0) || (vregidx < vreg.nr));

    // Will be used when detecting alleles

    if (getScores) {
      vreg.curr_bases.push_back(base);
      vreg.iids.push_back(iid);

      vreg.nb = vreg.curr_bases.size();
    }

    // Will be used when detecting variation

    if (((target_allele < 0)        ||   // use any allele
         (split_alleles == 0)       ||   // use any allele
         ((vreg.nr > 0)  &&
          (vreg.reads[vregidx].allele_id == target_allele)))) { // use the best allele
      bBaseCount[baseIdx]++;
      bQVSum[baseIdx] += qv;
      bReads.push_back(bead);
    } else {
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

  //  Compute tau based on guides
  //
  for (uint32 cind = 0; cind < gReads.size(); cind++) {
    abBead *gb   = gReads[cind];
    char    base = getBase(gb->baseIdx());
    uint32  qv   = getQual(gb->baseIdx()) - '0';

    used_surrogate = true;

    if (qv == 0)
      qv += 5;    /// HUH?!!

    tau[0] += (base == '-') ? PROB[qv] : EPROB[qv];
    tau[1] += (base == 'A') ? PROB[qv] : EPROB[qv];
    tau[2] += (base == 'C') ? PROB[qv] : EPROB[qv];
    tau[3] += (base == 'G') ? PROB[qv] : EPROB[qv];
    tau[4] += (base == 'T') ? PROB[qv] : EPROB[qv];

    consensusQV = qv;
  }

  //  If others, reset.
  //
  if (oReads.size() > 0)
    tau[0] = tau[1] = tau[2] = tau[3] = tau[4] = 1.0;

  //  Compute tau based on others
  //
  for (uint32 cind=0; cind<oReads.size(); cind++) {
    abBead  *gb   = oReads[cind];
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

  //  If real reads, reset.
  //
  if (bReads.size() > 0)
    tau[0] = tau[1] = tau[2] = tau[3] = tau[4] = 1.0;

  //  Compute tau based on real reads.
  //
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
  //
  if ((bReads.size() == 0) &&
      (oReads.size() == 0) &&
      (gReads.size() == 0)) {
    //fprintf(stderr, "No coverage for column=%d.  Assume it's an N in a single coverage area.\n", cid);

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
  fprintf(stderr, "prob('%c') %f %c\n", RALPHABET[0], cw[0], (cw[0] == cwMax) ? '*' : ' ');
  fprintf(stderr, "prob('%c') %f %c\n", RALPHABET[1], cw[1], (cw[1] == cwMax) ? '*' : ' ');
  fprintf(stderr, "prob('%c') %f %c\n", RALPHABET[2], cw[2], (cw[2] == cwMax) ? '*' : ' ');
  fprintf(stderr, "prob('%c') %f %c\n", RALPHABET[3], cw[3], (cw[3] == cwMax) ? '*' : ' ');
  fprintf(stderr, "prob('%c') %f %c\n", RALPHABET[4], cw[4], (cw[4] == cwMax) ? '*' : ' ');
#endif

  consensusBase = RALPHABET[cwIdx];

  //  If there isn't a clear winner
  if (cwMax < 1.0 - DBL_EPSILON) {
    int32  qv = consensusQV;

    //  If there is more than one fragment, or we used a surrogate (zero fragments) compute the QV
    //  from cwMax.  Otherwise, there is only one fragment, and we should use that qv (saved in
    //  consensusQV).
    //
    if ((frag_cov != 1) || (used_surrogate == true)) {
      double dqv =  -10.0 * log10(1.0 - cwMax);

      qv = DBL_TO_INT(dqv);

      if (dqv - qv >= 0.50)
        qv++;
    }

    qv = MIN(CNS_MAX_QV, qv);
    qv = MAX(CNS_MIN_QV, qv);

    consensusQV = qv + '0';
  } else {
    consensusQV = CNS_MAX_QV + '0';
  }

  setQual(call->baseIdx(), consensusQV);

  if ((target_allele  < 0) ||
      (target_allele == vreg.alleles[0].id)) {
    setBase(call->baseIdx(), consensusBase);
    setQual(call->baseIdx(), consensusQV);
  }

  // Detecting variation

  uint32   bReadCount = 0;
  uint32   ci         = 0;

  bReadCount += bBaseCount[0];  if (consensusBase == RALPHABET[0])  ci = 0;  //  '-' (RALPHABET)
  bReadCount += bBaseCount[1];  if (consensusBase == RALPHABET[1])  ci = 1;  //  'A'
  bReadCount += bBaseCount[2];  if (consensusBase == RALPHABET[2])  ci = 2;  //  'C'
  bReadCount += bBaseCount[3];  if (consensusBase == RALPHABET[3])  ci = 3;  //  'G'
  bReadCount += bBaseCount[4];  if (consensusBase == RALPHABET[4])  ci = 4;  //  'T'
  bReadCount += bBaseCount[5];  if (consensusBase == RALPHABET[5])  ci = 5;  //  'N'

  assert(RALPHABET[0] == '-');
  assert(RALPHABET[1] == 'A');
  assert(RALPHABET[2] == 'C');
  assert(RALPHABET[3] == 'G');
  assert(RALPHABET[4] == 'T');
  assert(RALPHABET[5] == 'N');

  uint32   sumQVall = 0;
  uint32   sumQVcns = 0;

  for (int32 bi=0; bi<5; bi++) {
    if ((bBaseCount[bi] < 2) ||
        (bBaseCount[ci] < 2))
      //  Not enough support for either the variation or the consensus.
      continue;

    if (consensusBase == RALPHABET[bi])
      //  This is the consensus base.
      continue;

    double aveQV = (double)bQVSum[bi] / bBaseCount[bi];
    uint32 sumQV = highest1_qv[bi] + highest2_qv[bi];

    bool   isGap = (consensusBase == '-') || (RALPHABET[bi] == '-');

    if (((isGap == false) && (aveQV >= MIN_AVE_QV_FOR_VARIATION)) ||
        ((isGap == true)  && (sumQV >= MIN_SUM_QVS_FOR_VARIATION))) {
      sumQVall  += bQVSum[bi];
      sumQVcns   = (RALPHABET[bi] == consensusBase) ? bQVSum[bi] : sumQVcns;
    }
  }

  if ((bReadCount > 1) &&
      (sumQVall > 0)) {
    double  ratio = (double)sumQVcns / sumQVall;

    var = ((smooth_win > 0) && (consensusBase == '-')) ? (ratio - 1.0) : (1.0 - ratio);
  } else {
    var = ((smooth_win > 0) && (consensusBase == '-')) ? -2.0 : 0.0;
  }
}



char
abAbacus::baseCall(abVarRegion &vreg,
                   abColID      cid,
                   bool         highQuality,
                   double      &var,
                   int32        target_allele,
                   bool         getScores,
                   int32        split_alleles,
                   int32        smooth_win) {

  abColumn *column = getColumn(cid);
  abBead   *call   = getBead(column->callID());

  //  NOTE: negative target_allele means the the alleles will be used
  //  Hardcoded in baseCallQuality
  assert(target_allele == -1);

  var     = 0.0;
  vreg.nb = 0;

  if (highQuality)
    baseCallQuality(vreg,
                    cid,
                    var,
                    target_allele,
                    getScores,
                    split_alleles,
                    smooth_win);
  else
    baseCallMajority(cid);

  return(getBase(call->baseIdx()));
}
