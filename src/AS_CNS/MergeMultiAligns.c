
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

static char *rcsid = "$Id: MergeMultiAligns.c,v 1.4 2009-07-30 10:42:56 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"

#undef DEBUG_MERGEMULTIALIGNS

static
MultiAlignT *
MergeMultiAligns(tSequenceDB *sequenceDBp,
                 gkStore *frag_store,
                 VA_TYPE(IntMultiPos) *positions,
                 int quality,
                 int verbose,
                 CNS_Options *opp) {
  // frag_store needed? no

  // C----------------------------C
  // u-------u     u---------u
  //        u-------u       u-----u
  //                             C----------------------------C
  //                       +     u----------------------------u
  MultiAlignT *cma;
  MANode *ma;
  int num_contigs;
  int32 num_columns=0;
  int complement;
  int32 fid,i,align_to;
  IntMultiPos *cpositions;
  SeqInterval *offsets;
  static VA_TYPE(int32) *trace=NULL;

  //  We need to reset the global sequenceDB pointer -- if we call
  //  this from anything but consensus, the global pointer isn't set.
  //
  sequenceDB = sequenceDBp;

  num_contigs = GetNumIntMultiPoss(positions);
  cpositions = GetIntMultiPos(positions,0);
  allow_neg_hang=0;
  USE_SDB=1;

  offsets = (SeqInterval *) safe_calloc(num_contigs,sizeof(SeqInterval));
  for (i=0;i<num_contigs;i++) {
    num_columns = ( cpositions[i].position.bgn>num_columns)? cpositions[i].position.bgn : num_columns;
    num_columns = ( cpositions[i].position.end>num_columns)? cpositions[i].position.end : num_columns;
  }

  gkpStore = frag_store;
  ResetStores(num_contigs,num_columns);

  if (num_contigs == 1) {
    cma = loadMultiAlignTFromSequenceDB(sequenceDBp, cpositions[0].ident, FALSE);
    safe_free(offsets);
    return cma;
  }

  for (i=0;i<num_contigs;i++) {
    complement = (cpositions[i].position.bgn<cpositions[i].position.end) ? 0 : 1;
    fid = AppendFragToLocalStore(cpositions[i].type,
                                 cpositions[i].ident,
                                 complement,
                                 0,
                                 AS_OTHER_UNITIG, NULL);
    offsets[fid].bgn = complement?cpositions[i].position.end:cpositions[i].position.bgn;
    offsets[fid].end = complement?cpositions[i].position.bgn:cpositions[i].position.end;

    if (VERBOSE_MULTIALIGN_OUTPUT)
      fprintf(stderr,"MergeMultiAligns()-- id=%10d %s %c %12d %12d\n",
              cpositions[i].ident,
              (complement) ? "<----" : "---->",
              cpositions[i].type,
              offsets[fid].bgn,
              offsets[fid].end);
  }

  ma = CreateMANode(cpositions[0].ident);

  if (trace == NULL)
    trace = CreateVA_int32(AS_READ_MAX_LEN);
  else
    ResetVA_int32(trace);

  SeedMAWithFragment(ma->lid, GetFragment(fragmentStore,0)->lid,0,opp);

  // Now, loop on remaining fragments, aligning to:
  //    a)  containing frag (if contained)
  // or b)  previously aligned frag
  for (i=1;i<num_contigs;i++) {
    int ahang,bhang,ovl;
    OverlapType otype;
    int olap_success=0;
    int try_contained=0;
    Fragment *afrag = NULL;
    Fragment *bfrag = GetFragment(fragmentStore,i);

    align_to = i-1;

    while (!olap_success) {
      //  Skip contained stuff.
      while ((align_to > 0) &&
             ((GetFragment(fragmentStore, align_to)->is_contained) ||
              (GetFragment(fragmentStore, align_to)->container_iid > 0)))
        align_to--;

      if (align_to < 0) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MergeMultiAligns()-- unable to find uncontained contig upstream from current contig %d.  Fail.\n",
                  bfrag->iid);
        break;
      }

      afrag = GetFragment(fragmentStore, align_to);

      //  ECR violates the usual "positive ahang" rule.
      //assert(offsets[afrag->lid].bgn <= offsets[bfrag->lid].bgn);

      //  This code copied from MultiAlignUnitig.  Placement of some
      //  contained (and others?)  sequences seems screwed up, so
      //  we're using the length instead.

      ahang = offsets[bfrag->lid].bgn - offsets[afrag->lid].bgn;
      bhang = offsets[bfrag->lid].end - offsets[afrag->lid].end;

      if (offsets[afrag->lid].bgn < offsets[bfrag->lid].bgn)
        if (offsets[afrag->lid].end < offsets[bfrag->lid].end)
          ovl = offsets[afrag->lid].end - offsets[bfrag->lid].bgn;
        else
          //ovl = offsets[bfrag->lid].end - offsets[bfrag->lid].bgn;
          ovl = bfrag->length;
      else
        if (offsets[afrag->lid].end < offsets[bfrag->lid].end)
          //ovl = offsets[afrag->lid].end - offsets[afrag->lid].bgn;
          ovl = afrag->length;
        else
          ovl = offsets[bfrag->lid].end - offsets[afrag->lid].bgn;

      //  End of copy

      if (VERBOSE_MULTIALIGN_OUTPUT)
        fprintf(stderr, "MergeMultiAligns()-- ovl=%d from a %d-%d b %d-%d hangs %d %d\n",
                ovl,
                offsets[afrag->lid].bgn, offsets[afrag->lid].end,
                offsets[bfrag->lid].bgn, offsets[bfrag->lid].end,
                ahang, bhang);

      if (ovl <= 0) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MergeMultiAligns()-- uncontained contig upstream is found, but positions indicate no overlap between contigs %d and %d.  Fail.", afrag->iid, bfrag->iid);

        DeleteMANode(ma->lid);
        safe_free(offsets);
        return NULL;
      }

      olap_success = GetAlignmentTrace(afrag->lid, NULL, bfrag->lid, &ahang, &bhang, ovl, trace, &otype, DP_Compare, DONT_SHOW_OLAP, 0, AS_MERGE, AS_CGW_ERROR_RATE);

      if (!olap_success)
        olap_success = GetAlignmentTrace(afrag->lid, NULL, bfrag->lid, &ahang, &bhang, ovl, trace, &otype, Local_Overlap_AS_forCNS, DONT_SHOW_OLAP, 0, AS_MERGE, AS_CGW_ERROR_RATE);

      if (!olap_success) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MergeMultiAligns()-- positions of contigs %d and %d overlap, but GetAlignmentTrace() fails.\n",
                  afrag->iid, bfrag->iid);
        break;
      }
    }

    if (!olap_success) {
      if (VERBOSE_MULTIALIGN_OUTPUT)
        fprintf(stderr, "MergeMultiAligns()-- failed to find overlap between contigs %d and %d.  Fail.\n",
                afrag->iid, bfrag->iid);
      DeleteMANode(ma->lid);
      safe_free(offsets);
      return NULL;
    }

    if (otype == AS_CONTAINMENT) {
      bfrag->is_contained = 1;
      if (bfrag->container_iid == 0)
        bfrag->container_iid = 1;  //  Not sure why 1 and not afrag->iid
    }

    ApplyAlignment(afrag->lid, 0, NULL, bfrag->lid, ahang, Getint32(trace,0));
  } /* loop through all contigs */

  {
    IntMultiVar *vl;
    int32 nv;
    RefreshMANode(ma->lid, 0, opp, &nv, &vl, 0, 0);
  }

  // Now, want to generate a new MultiAlignT which merges the u_list and f_list of the contigs
  // merge the f_lists and u_lists by cloning and concating (or constructing dummy, when dealing with single read

  int ifrag;
  int iunitig;
  IntMultiPos *imp;
  IntUnitigPos *iup;

  cma = CreateMultiAlignT();
  cma->consensus = CreateVA_char(GetMANodeLength(ma->lid)+1);
  cma->quality   = CreateVA_char(GetMANodeLength(ma->lid)+1);

  GetMANodeConsensus(ma->lid, cma->consensus, cma->quality);

  // no deltas required at this stage
  cma->fdelta = CreateVA_int32(0);
  cma->udelta = CreateVA_int32(0);

  if ((cpositions[0].type == AS_UNITIG) || (cpositions[0].type == AS_CONTIG)) {
    MultiAlignT *ma = loadMultiAlignTFromSequenceDB(sequenceDBp, cpositions[0].ident, cpositions[0].type == AS_UNITIG);
    cma->f_list = Clone_VA(ma->f_list);
    cma->v_list = Clone_VA(ma->v_list);
    cma->u_list = Clone_VA(ma->u_list);
  } else {
    assert(cpositions[0].type == AS_READ);
    cma->f_list = CreateVA_IntMultiPos(0);
    cma->v_list = CreateVA_IntMultiVar(0);
    cma->u_list = CreateVA_IntUnitigPos(0);
    AppendVA_IntMultiPos(cma->f_list,cpositions+0);
  }

  for (i=1;i<num_contigs;i++) {
    if ((cpositions[i].type == AS_UNITIG) || (cpositions[i].type == AS_CONTIG)) {
      MultiAlignT *ma = loadMultiAlignTFromSequenceDB(sequenceDBp, cpositions[i].ident, cpositions[i].type == AS_UNITIG);
      ConcatVA_IntMultiPos(cma->f_list,ma->f_list);
      ConcatVA_IntMultiPos(cma->v_list,ma->v_list);
      ConcatVA_IntUnitigPos(cma->u_list,ma->u_list);
    } else {
      assert(cpositions[i].type == AS_READ);
      AppendVA_IntMultiPos(cma->f_list,cpositions+i);
    }
  }

  ifrag=0;
  iunitig=0;
  for (i=0;i<num_contigs;i++) {
    Fragment *cfrag=GetFragment(fragmentStore,i);  /* contig pseudo-frag */

    if ((cfrag->type == AS_UNITIG) || (cfrag->type == AS_CONTIG)) {
      CNS_AlignedContigElement *components=GetCNS_AlignedContigElement(fragment_positions,cfrag->components);
      CNS_AlignedContigElement *compci;

      int ci=0;
      int32 len,bgn,end,left,right,tmp;


      //  cfrag->length seems to be including gaps (i.e., out of
      //  consensus), but we need to use the length it is placed at to
      //  reverse-complement positions.  So, use the difference
      //  between the largest and smallest coordiante as the length.
      
      len = cfrag->length;
      bgn = len;
      end = 0;

      for (ci =0; ci<cfrag->n_components; ci++) {
        compci = &components[ci];

        if (compci->position.bgn < bgn)  bgn = compci->position.bgn;
        if (compci->position.end < bgn)  bgn = compci->position.end;
        if (compci->position.bgn > end)  end = compci->position.bgn;
        if (compci->position.end > end)  end = compci->position.end;
      }

      len = end - bgn;

      // make adjustments to positions

      ci = 0;
      while (ci < cfrag->n_components) {
        compci = &components[ci];

#ifdef DEBUG_MERGEMULTIALIGNS
        if (compci->frg_or_utg==CNS_ELEMENT_IS_UNITIG)
          fprintf(stderr, "compci complement=%d length=%d bgn=%d end=%d\n", cfrag->complement, cfrag->length, compci->position.bgn, compci->position.end);
#endif

        if ( cfrag->complement ) {
          bgn = len - compci->position.bgn;
          end = len - compci->position.end;
        } else {
          bgn = compci->position.bgn;
          end = compci->position.end;
        }

        left  = (bgn < end) ? bgn : end;
        right = (bgn < end) ? end : bgn;

        if (left < 0)
          left = 0;
        if (right > len)
          right = len;

#ifdef DEBUG_MERGEMULTIALIGNS
        if (compci->frg_or_utg==CNS_ELEMENT_IS_UNITIG)
          fprintf(stderr, "left=%d right=%d bgn=%d end=%d\n", left, right, bgn, end);
#endif

        left  = GetColumn(columnStore, GetBead(beadStore,cfrag->firstbead + left)   ->column_index)->ma_index;
        right = GetColumn(columnStore, GetBead(beadStore,cfrag->firstbead + right-1)->column_index)->ma_index + 1;

        tmp = bgn;
        bgn = (bgn < end) ? left  : right;
        end = (tmp < end) ? right : left;

#ifdef DEBUG_MERGEMULTIALIGNS
        if (compci->frg_or_utg==CNS_ELEMENT_IS_UNITIG)
          fprintf(stderr, "left=%d right=%d bgn=%d end=%d\n", left, right, bgn, end);
#endif

        if (compci->frg_or_utg==CNS_ELEMENT_IS_UNITIG) {
          iup = GetIntUnitigPos(cma->u_list,iunitig);
          iup->position.bgn = bgn;
          iup->position.end = end;
          iup->delta_length = 0;
          iup->delta = NULL;

#ifdef DEBUG_MERGEMULTIALIGNS
          fprintf(stderr, "Placing IUP  "F_CID" at "F_S32","F_S32" based on positions "F_S32","F_S32" (compl %d length %d within input parent)\n",
                  iup->ident, bgn, end, compci->position.bgn, compci->position.end, cfrag->complement, len);
#endif

          ci++;
          iunitig++;
        } else {
          imp = GetIntMultiPos(cma->f_list,ifrag);
          imp->ident = compci->idx.fragment.frgIdent;
          imp->sourceInt = compci->idx.fragment.frgSource;
          imp->position.bgn = bgn;
          imp->position.end = end;
          imp->delta_length = 0;
          imp->delta = NULL;

#ifdef DEBUG_MERGEMULTIALIGNS
          //  Generally not interesting.
          //fprintf(stderr, "Placing IMP1 "F_CID" at "F_S32","F_S32" based on positions "F_S32","F_S32" (compl %d length %d within input parent)\n",
          //        imp->ident, bgn, end, compci->position.bgn, compci->position.end, cfrag->complement, len);
#endif

          ci++;
          ifrag++;
        }
      }
    } else {

      int32 bgn,end;

      assert(cfrag->type == AS_READ);

      //  cfrag is a fragment, so length should be ungapped

      bgn = GetBead(beadStore,cfrag->firstbead)->column_index;
      end = GetBead(beadStore,cfrag->firstbead + cfrag->length - 1 )->column_index + 1;

      if(cfrag->complement){
        int32 tmp = bgn;
        bgn = end;
        end = tmp;
      }

      imp = GetIntMultiPos(cma->f_list,ifrag);
      imp->position.bgn = bgn;
      imp->position.end = end;

#ifdef DEBUG_MERGEMULTIALIGNS
      //  Generally not interesting.
      //fprintf(stderr, "Placing IMP2 "F_CID" at "F_S32","F_S32" based on positions "F_S32","F_S32" (compl %d length %d within input parent)\n",
      //        imp->ident, bgn,end, offsets[i].bgn, offsets[i].end, cfrag->complement, cfrag->length);
#endif

      ifrag++;
    }
  }

  DeleteMANode(ma->lid);
  safe_free(offsets);

  return cma;
}


// MergeMultiAlignsFast_new is the original CGW/CNS interface for
// contigging and is now a wrapper around the more modern
// MergeMultiAligns which allows "contained" relationships among the
// input contigs

MultiAlignT *
MergeMultiAlignsFast_new(tSequenceDB *sequenceDBp,
                         gkStore *frag_store,
                         VA_TYPE(IntElementPos) *positions,
                         int quality,
                         int verbose,
                         CNS_Options *opp) {

  static VA_TYPE(IntMultiPos) *mpositions=NULL;
  static IntMultiPos mpos;

  IntElementPos *epos = GetIntElementPos(positions,0);
  int i;

#ifdef DEBUG_MERGEMULTIALIGNS
  VERBOSE_MULTIALIGN_OUTPUT = 1;
#endif

  if (mpositions == NULL )
    mpositions = CreateVA_IntMultiPos(32);

  ResetVA_IntMultiPos(mpositions);

  mpos.contained    = 0;
  mpos.delta_length = 0;
  mpos.delta        = NULL;

  for (i=0; i<GetNumIntElementPoss(positions); i++, epos++) {
    mpos.type     = epos->type;
    mpos.ident    = epos->ident;
    mpos.position = epos->position;

    AppendVA_IntMultiPos(mpositions,&mpos);
  }

  allow_neg_hang = 0;

  return(MergeMultiAligns(sequenceDBp, frag_store, mpositions, quality, verbose, opp));
}
