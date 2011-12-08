
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

static char *rcsid = "$Id: BaseCall.c,v 1.5 2011-12-08 00:24:36 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"

VA_DEF(int16)

void
BaseCall(int32 cid,
         int32 quality,
         double *var,
         VarRegion  *vreg,
         int32 target_allele,
         char *cons_base,
         int32 get_scores,
         CNS_Options *opp) {

  // NOTE: negative target_allele means the the alleles will be used

  Column *column=GetColumn(columnStore,cid);
  Bead *call = GetBead(beadStore, column->call);
  Bead *bead;
  int32 best_read_base_count[CNS_NP]  = {0};
  int32 other_read_base_count[CNS_NP] = {0};
  int32 guide_base_count[CNS_NP]      = {0};

  char bases[CNS_NALPHABET] = {'-', 'A', 'C', 'G', 'T', 'N'};
  int32 best_read_qv_count[CNS_NP] = {0};
  int32 other_read_qv_count[CNS_NP] = {0};
  int32 highest_qv[CNS_NP] = {0};
  int32 highest2_qv[CNS_NP] = {0};

  int32 b_read_depth=0, o_read_depth=0, guide_depth=0;
  int32 score=0;
  int32 bi, bi_cons=0;
  beadIdx bid;
  int32 iid = 0;
  char cqv, cbase;
  int32 qv = 0;
  static  double cw[CNS_NP];      // "consensus weight" for a given base
  static  double tau[CNS_NP];
  int32            tauValid = 0;
  FragType type;
  UnitigType utype;
  ColumnBeadIterator ci;
  int32 used_surrogate=0;
  int32 sum_qv_cbase=0, sum_qv_all=0;

  vreg->nb = 0;

  //  Make sure that we have valid options here, we then reset the
  //  pointer to the freshly copied options, so that we can always
  //  assume opp is a valid pointer
  //
  CNS_Options  opp_private;
  if (opp == NULL) {
    opp_private.split_alleles   = CNS_OPTIONS_SPLIT_ALLELES_DEFAULT;
    opp_private.smooth_win      = CNS_OPTIONS_MIN_ANCHOR_DEFAULT;
    opp = &opp_private;
  }

  CreateColumnBeadIterator(cid, &ci);

  *var = 0.;

  if (quality > 0) {
    static int32 guides_alloc=0;
    static VarArrayBead  *guides;
    static VarArrayBead  *b_reads;
    static VarArrayBead  *o_reads;
    static VarArrayint16 *tied;
    uint32 bmask;
    int32    num_b_reads, num_o_reads, num_guides;
    Bead  *gb;
    int32    cind;
    double tmpqv;
    int16  bi;
    int32    b_read_count = 0;
    int32    frag_cov=0;
    int16  max_ind=0;
    double max_cw=0.0;   // max of "consensus weights" of all bases
    double normalize=0.;
    int32    nr=0, max_nr=128;

    if (!guides_alloc) {
      guides = CreateVA_Bead(16);
      b_reads  = CreateVA_Bead(16);
      o_reads  = CreateVA_Bead(16);
      tied   = CreateVA_int16(32);
      guides_alloc = 1;
    } else {
      ResetBead(guides);
      ResetBead(b_reads);
      ResetBead(o_reads);
      Resetint16(tied);
    }
    for (bi=0;bi<CNS_NP;bi++) {
      tau[bi] = 1.0;
    }

    // Scan a column of aligned bases (=beads).
    // Sort the beads into three groups:
    //      - those corresponding to the reads of the best allele,
    //      - those corresponding to the reads of the other allele and
    //      - those corresponding to non-read fragments (aka guides)
    while ( (bid = NextColumnBead(&ci)).isValid()) {
      bead =  GetBead(beadStore,bid);
      cbase = *Getchar(sequenceStore,bead->soffset);    // current base
      qv = (int) ( *Getchar(qualityStore,bead->soffset)-'0');
      if ( cbase == 'N' ) {
        // skip 'N' base calls
        // fprintf(stderr,
        //    "encountered 'n' base in fragment data at column cid=%d\n",
        //    cid);
        continue;
      }
      bmask = AMASK[RINDEX[cbase]];
      type  = GetFragment(fragmentStore,bead->frag_index)->type;
      iid   = GetFragment(fragmentStore,bead->frag_index)->iid;

      if ((type == AS_READ)   ||
          (type == AS_EXTR)   ||
          (type == AS_TRNR)) {
        int32 k;

        for (k=0; k<vreg->nr; k++)
          if (iid == vreg->iids[k])
            break;

        assert((vreg->nr <= 0) || (k < vreg->nr));

        // Will be used when detecting alleles
        if (target_allele < 0 && get_scores) {
          vreg->curr_bases[vreg->nb] = cbase;
          vreg->iids[vreg->nb]  = iid;
          vreg->nb++;
          if (vreg->nb == vreg->max_nr) {
            vreg->max_nr += INITIAL_NR;
            vreg->curr_bases = (char *)safe_realloc(vreg->curr_bases, vreg->max_nr*sizeof(char));
            vreg->iids = (int32 *)safe_realloc(vreg->iids, vreg->max_nr*sizeof(int32));
          }
        }

        // Will be used when detecting variation
        if (((target_allele < 0)  ||   // use any allele
             !opp->split_alleles  ||   // use any allele
             (vreg->nr >  0  &&
              vreg->reads[k].allele_id == target_allele))) { // use the best allele
          best_read_base_count[RINDEX[cbase]]++;
          best_read_qv_count[RINDEX[cbase]] += qv;
          AppendBead(b_reads, bead);
        } else {
          other_read_base_count[RINDEX[cbase]]++;
          other_read_qv_count[RINDEX[cbase]] += qv;
          AppendBead(o_reads, bead);
        }

        if (highest_qv[RINDEX[cbase]] < qv) {
          highest2_qv[RINDEX[cbase]] = highest_qv[RINDEX[cbase]];
          highest_qv[RINDEX[cbase]] = qv;
        } else if (highest_qv[RINDEX[cbase]] >= qv &&
                   highest2_qv[RINDEX[cbase]] < qv) {
          highest2_qv[RINDEX[cbase]] = qv;
        }
      } else {
        guide_base_count[RINDEX[cbase]]++;
        AppendBead(guides, bead);
      }

      if ( type != AS_UNITIG ) {
        frag_cov++;
      }
    }

    b_read_depth = GetNumBeads(b_reads);
    o_read_depth = GetNumBeads(o_reads);
    guide_depth  = GetNumBeads(guides);

    // For each base, calculate tau
    // It will be used to calculate cw
    if (b_read_depth > 0) {
      for (cind = 0; cind < b_read_depth; cind++) {
        gb = GetBead(b_reads, cind);
        cbase = *Getchar(sequenceStore,gb->soffset);
        qv = (int) ( *Getchar(qualityStore,gb->soffset)-'0');
        if ( qv == 0 )
          qv += 5;
        bmask = AMASK[RINDEX[cbase]];
        for (bi=0;bi<CNS_NP;bi++) {
          if ( (bmask>>bi) & 1 ) {
            tau[bi]*= PROB[qv];
          } else {
            tau[bi]*= (double) TAU_MISMATCH * EPROB[qv];
          }
        }
        tauValid = 1;
      }
    } else {
      for (cind = 0; cind < o_read_depth; cind++) {
        gb = GetBead(o_reads, cind);
        cbase = *Getchar(sequenceStore,gb->soffset);
        qv = (int) ( *Getchar(qualityStore,gb->soffset)-'0');
        if ( qv == 0 )
          qv += 5;
        bmask = AMASK[RINDEX[cbase]];
        for (bi=0;bi<CNS_NP;bi++) {
          if ( (bmask>>bi) & 1 ) {
            tau[bi]*= PROB[qv];
          } else {
            tau[bi]*= (double) TAU_MISMATCH * EPROB[qv];
          }
        }
        tauValid = 1;
      }
    }

    // If there are no reads, assume it is because a surrogate is there.

    if (b_read_depth == 0 && o_read_depth == 0) {
      for (cind = 0; cind < guide_depth; cind++) {
        gb = GetBead(guides,cind);

        type  = GetFragment(fragmentStore,gb->frag_index)->type;
        utype = GetFragment(fragmentStore,gb->frag_index)->utype;

#if 0
        //  If your assembly absolutely needs to get done, you can
        //  disable this block, and we'll use whatever unitig is
        //  here to get the consensus.
        //
        //  Originally, we only wanted to use the unitig consensus
        //  for places where a surrogate could be.  And so, for any
        //  unitig that isn't a surrogate-type, we skip.
        //
        if ((type == AS_UNITIG) &&
            (utype != AS_STONE_UNITIG) &&
            (utype != AS_PEBBLE_UNITIG) &&
            (utype != AS_OTHER_UNITIG))
          //  Not a surrogate unitig!
          continue;
#else
        //fprintf(stderr, "WARNING: no coverage for column=%d (b_read_depth=%d o_read_depth=%d guide_depth=%d)\n",
        //        cid, b_read_depth, o_read_depth, guide_depth);
#endif

        used_surrogate=1;

        cbase = *Getchar(sequenceStore,gb->soffset);
        qv = (int) ( *Getchar(qualityStore, gb->soffset)-'0');
        if ( qv == 0 )
          qv += 5;
        bmask = AMASK[RINDEX[cbase]];
        for (bi=0; bi<CNS_NP; bi++) {
          if ( (bmask>>bi) & 1 ) {
            tau[bi] *= PROB[qv];
          } else {
            tau[bi] *= (double) TAU_MISMATCH * EPROB[qv];
          }
        }
        tauValid = 1;
      }
    }

    //  BPW has seen a low-coverage 454 assembly (gs20) where some
    //  unitigs had one read of coverage, but the base in the read
    //  was an N (use the "-v 4" option to consensus to see the
    //  multialign).  The only recognizable pattern was all read
    //  depths were zero.  The N showed up as base_count[5] == 1,
    //  all others were zero (others being gap, a, c, g, t).
    //
    //  Warn, and continue.
    //
    if ((tauValid == 0) && (b_read_depth == 0) && (o_read_depth == 0) && (guide_depth == 0)) {
      //fprintf(stderr, "No coverage (tau not valid) for column=%d.  Assume it's an N in a single coverage area.\n",
      //        cid);
      tauValid = 1;
    }

    //  If you hit this assert (and are working on a contig), you
    //  can enable the above ifdef block to use whatever unitig is
    //  here.  It's not enabled by default because it should never
    //  happen.
    //
    assert(tauValid);


    max_ind = 0;      // consensus is gap
    max_cw  = 0.0;

    //  This is gross.

    for (bi=0; bi<5; bi++) {
      cw[bi]     = tau[bi] * 0.2;
      normalize += cw[bi];
    }

    if (normalize)
      normalize = 1./normalize;

    // Calculate max_ind as {i | cw[i] -> max_cw}
    // Store all other indexes { i | cw[i] == max_cw } in VA Array tied
    for (bi=0; bi<5; bi++) {
      cw[bi] *= normalize;

      if (cw[bi] > max_cw + ZERO_PLUS) {
        max_ind = bi;
        max_cw = cw[bi];
        Resetint16(tied);
      } else if (DBL_EQ_DBL(cw[bi], max_cw)) {
        Appendint16(tied,&bi);
      }
    }

    // If max_cw == 0, then consensus base call will be a gap
    // (max_ind==0)
    //
    // Otherwise, it will be selected RANDOMLY (!!!) from all
    // {i|cw[i]==max_cw}

    if (DBL_EQ_DBL(max_cw, (double)0.0)) {
      max_ind = 0;      // consensus is gap
    } else {
      if (GetNumint16s(tied)> 0) {
        Appendint16(tied, &max_ind);
        max_ind = *Getint16(tied,1);
        max_cw = cw[max_ind];
      }
    }


#if 0
    fprintf(stdout,"calculated probabilities:\n");
    for (bi=0;bi<CNS_NP;bi++)
      fprintf(stdout,"%c = %16.8f %c\n", RALPHABET[bi],cw[bi], (bi == max_ind) ? '*' : ' ');
#endif

    // Set the consensus base quality value
    cbase = RALPHABET[max_ind];      // consensus base
    if (DBL_EQ_DBL(max_cw, (double)1.0)) {
      cqv = CNS_MAX_QV+'0';
      Setchar(qualityStore, call->soffset, &cqv);
    } else {
      if ( frag_cov != 1 || used_surrogate) {
        tmpqv =  -10.0 * log10(1.0-max_cw);
        qv = DBL_TO_INT(tmpqv);
        if ((tmpqv - qv)>=.50)
          qv++;
      }

      if      (qv > CNS_MAX_QV)
        cqv = '0' + CNS_MAX_QV;
      else if (qv < CNS_MIN_QV)
        cqv = '0' + CNS_MIN_QV;
      else
        cqv = '0' + qv;
    }


    //  BPW wants to ensure that QV is terrible if there is no read
    //  coverage -- see comments above.  We also force the basecall
    //  to N, probably causing future problems.
    //
    if ((b_read_depth == 0) && (o_read_depth == 0) && (guide_depth == 0)) {
      cbase = 'N';
      cqv   = '0';
    }

    *cons_base = cbase;

    if (target_allele <  0 || target_allele == vreg->alleles[0].id) {
      Setchar(sequenceStore, call->soffset, &cbase);
      Setchar(qualityStore, call->soffset, &cqv);
    }

    // Detecting variation
    for (bi=0; bi<CNS_NALPHABET-1; bi++) {
      b_read_count += best_read_base_count[bi];
      if (*cons_base == bases[bi])
        bi_cons = bi;
    }

    for (bi=0; bi<CNS_NALPHABET-1; bi++) {
      // NALAPHBET-1 to exclude "n" base call
      bmask = AMASK[bi];  // mask for indicated base
      if ( ! ((bmask>>max_ind) & 1) ) {
        // penalize only if base in not represented in call
        score += best_read_base_count[bi] + other_read_base_count[bi]
          + guide_base_count[bi];
      }
      // To be considered, varied base should be confirmed by another base
      // and either overall quality should be >= 60
      // (Granger's suggestion - GD)
      //
      if (*cons_base != '-' && bases[bi] != '-' && *cons_base != bases[bi] &&
          best_read_base_count[bi     ] > 1 &&   // variation represented by >=2 bases
          best_read_base_count[bi_cons] > 1 &&   // consensus represented by >=2 bases
          (float)best_read_qv_count[bi]/(float)best_read_base_count[bi] >= MIN_AVE_QV_FOR_VARIATION) {
        sum_qv_all += best_read_qv_count[bi];
        if (RALPHABET[bi] == cbase)
          sum_qv_cbase = best_read_qv_count[bi];
      }
      else if ((*cons_base == '-' || bases[bi] == '-') && *cons_base != bases[bi] &&
               best_read_base_count[bi     ] > 1 && // variation represented by >=2 bases
               best_read_base_count[bi_cons] > 1 && // consensus represented by >=2 bases
               highest_qv[bi] + highest2_qv[bi] >= MIN_SUM_QVS_FOR_VARIATION) {
        sum_qv_all += best_read_qv_count[bi];
        if (RALPHABET[bi] == cbase)
          sum_qv_cbase = best_read_qv_count[bi];
      }
    }

    if ((b_read_count <= 1 ) || (sum_qv_all == 0)) {
      *var = 0.;
      if (opp->smooth_win > 0 && cbase == '-')
        *var = -2;
    } else {
      *var = 1. - (double)sum_qv_cbase / (double)sum_qv_all;
      if (opp->smooth_win > 0 && cbase == '-') {
        *var = - (*var);
      }
    }

  } else if (quality == 0 ) {
    int32 max_count=0,max_index=-1,tie_count=0;
    int32 tie_breaker, max_tie, i;

    CreateColumnBeadIterator(cid, &ci);

    while ( (bid = NextColumnBead(&ci)).isValid() ) {
      bead = GetBead(beadStore,bid);
      cbase = *Getchar(sequenceStore,bead->soffset);
      qv = (int) ( *Getchar(qualityStore, bead->soffset)-'0');
      type = GetFragment(fragmentStore,bead->frag_index)->type;
      if (type  != AS_READ &&
          type  != AS_EXTR &&
          type  != AS_TRNR ) {
        guide_base_count[RINDEX[cbase]]++;
      } else {
        best_read_base_count[RINDEX[cbase]]++;
      }
    }
    for (i=0; i<CNS_NALPHABET; i++) {
      if (best_read_base_count[i]+guide_base_count[i] > max_count) {
        max_count = best_read_base_count[i] + guide_base_count[i];
        max_index = i;
      }
    }
    if ( best_read_base_count[max_index] + guide_base_count[max_index] >
         (b_read_depth                 + guide_depth)/2 ) {
      tie_count = 0;
    } else {
      for (i=0;i<CNS_NALPHABET;i++) {
        if (best_read_base_count[i]+guide_base_count[i] == max_count) {
          max_index = i;
          tie_count++;
        }
      }
    }

    max_tie=-1;

    if ( tie_count > 1 ) {
      for (i=1;i<CNS_NALPHABET;i++) {
        //  i starts at 1 to prevent ties from being broken with '-'
        if ( best_read_base_count[i]+guide_base_count[i] == max_count ) {
          tie_breaker = random();
          if (tie_breaker > max_tie) {
            max_tie = tie_breaker;
            max_index = i;
          }
        }
      }
    }

    cbase=toupper(RALPHABET[max_index]);
    Setchar(sequenceStore, call->soffset, &cbase);
    cqv = 0 + '0';
    Setchar(qualityStore, call->soffset, &cqv);

    for (bi=0;bi<CNS_NALPHABET;bi++)
      if (bi != RINDEX[cbase])
        score += best_read_base_count[bi]+guide_base_count[bi];

  } else if (quality == -1 ) {
    // here, just promote the aligned fragment's seq and quality to the basecall
    char bqv;
    bid = NextColumnBead(&ci);
    bead =  GetBead(beadStore,bid);
    cbase = *Getchar(sequenceStore, bead->soffset);
    bqv  = *Getchar(qualityStore,bead->soffset);
    Setchar(sequenceStore, call->soffset, &cbase);
    Setchar(qualityStore, call->soffset, &bqv);
  }
}
