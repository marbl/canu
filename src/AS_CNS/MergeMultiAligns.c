
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

static char *rcsid = "$Id: MergeMultiAligns.c,v 1.12 2012-09-13 09:56:07 brianwalenz Exp $";

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

// C----------------------------C
// u-------u     u---------u
//        u-------u       u-----u
//                             C----------------------------C
//                       +     u----------------------------u

MultiAlignT *
MergeMultiAlignsFast_new(VA_TYPE(IntElementPos) *positions, CNS_Options *opp) {
  int32          num_contigs = GetNumIntElementPoss(positions);
  IntElementPos *cpositions  = GetIntElementPos(positions, 0);

#ifdef DEBUG_MERGEMULTIALIGNS
  VERBOSE_MULTIALIGN_OUTPUT = 1;
#endif

  if (num_contigs == 1) {
    MultiAlignT *cma = tigStore->loadMultiAlign(cpositions[0].ident, FALSE);

    if (cma == NULL)
      fprintf(stderr, "MergeMultiAlignsFast_new() returns null cma from existing contig %d?\n", cpositions[0].ident);

    return(cma);
  }

  allow_neg_hang = 0;

  int32        num_bases   = 0;
  int32        num_columns = 0;

  for (int32 i=0; i<num_contigs; i++) {
    int32 flen   = (cpositions[i].position.bgn < cpositions[i].position.end) ? (cpositions[i].position.end < cpositions[i].position.bgn) : (cpositions[i].position.bgn - cpositions[i].position.end);
    num_bases   += flen + 2 * AS_CNS_ERROR_RATE * flen;

    num_columns = (cpositions[i].position.bgn > num_columns) ? cpositions[i].position.bgn : num_columns;
    num_columns = (cpositions[i].position.end > num_columns) ? cpositions[i].position.end : num_columns;
  }

  ResetStores(num_bases, num_contigs, num_columns);

  SeqInterval *offsets     = (SeqInterval *)safe_calloc(num_contigs,sizeof(SeqInterval));

  for (int32 i=0; i<num_contigs; i++) {
    int32 complement = (cpositions[i].position.bgn<cpositions[i].position.end) ? 0 : 1;
    int32 fid        = AppendFragToLocalStore(cpositions[i].type,
                                              cpositions[i].ident,
                                              complement,
                                              0,
                                              AS_OTHER_UNITIG);

    assert(cpositions[i].type == AS_CONTIG);

    offsets[fid].bgn = complement ? cpositions[i].position.end : cpositions[i].position.bgn;
    offsets[fid].end = complement ? cpositions[i].position.bgn : cpositions[i].position.end;

    if (VERBOSE_MULTIALIGN_OUTPUT)
      fprintf(stderr,"MergeMultiAlignsFast_new()-- id=%10d %s %c %12d %12d\n",
              cpositions[i].ident,
              (complement) ? "<----" : "---->",
              cpositions[i].type,
              offsets[fid].bgn,
              offsets[fid].end);
  }

  MANode *manode = CreateMANode(cpositions[0].ident);

  SeedMAWithFragment(manode->lid, GetFragment(fragmentStore,0)->lid,opp);

  // Now, loop on remaining fragments, aligning to:
  //    a)  containing frag (if contained)
  // or b)  previously aligned frag

  VA_TYPE(int32) *trace = CreateVA_int32(AS_READ_MAX_NORMAL_LEN);

  for (int32 i=1; i<num_contigs; i++) {
    int32          ahang,bhang,ovl;
    OverlapType    otype;
    int32          olap_success=0;
    int32          try_contained=0;
    Fragment      *afrag = NULL;
    Fragment      *bfrag = GetFragment(fragmentStore,i);

    int32          align_to = i-1;

    while (!olap_success) {
      //  Skip contained stuff.
      while ((align_to > 0) &&
             ((GetFragment(fragmentStore, align_to)->is_contained) ||
              (GetFragment(fragmentStore, align_to)->container_iid > 0)))
        align_to--;

      if (align_to < 0) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MergeMultiAlignsFast_new()-- unable to find uncontained contig upstream from current contig %d.  Fail.\n",
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
        fprintf(stderr, "MergeMultiAlignsFast_new()-- ovl=%d from a %d-%d b %d-%d hangs %d %d\n",
                ovl,
                offsets[afrag->lid].bgn, offsets[afrag->lid].end,
                offsets[bfrag->lid].bgn, offsets[bfrag->lid].end,
                ahang, bhang);

      if (ovl <= 0) {
        fprintf(stderr, "MergeMultiAlignsFast_new()-- uncontained contig upstream is found, but positions indicate no overlap between contigs %d and %d.\n",
                afrag->iid, bfrag->iid);

        DeleteMANode(manode->lid);
        safe_free(offsets);
        return(NULL);
      }

      //  Use driver here?  Probably not; that will increase the erate and try other tricks we don't
      //  really want to have done *yet).

      olap_success = GetAlignmentTrace(afrag->lid, NULL, bfrag->lid, &ahang, &bhang, ovl, trace, &otype, DP_Compare, DONT_SHOW_OLAP, 0, GETALIGNTRACE_MERGE, AS_CGW_ERROR_RATE);

      if (!olap_success)
        olap_success = GetAlignmentTrace(afrag->lid, NULL, bfrag->lid, &ahang, &bhang, ovl, trace, &otype, Local_Overlap_AS_forCNS, DONT_SHOW_OLAP, 0, GETALIGNTRACE_MERGE, AS_CGW_ERROR_RATE);

      if (!olap_success) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MergeMultiAlignsFast_new()-- positions of contigs %d and %d overlap, but GetAlignmentTrace() fails.\n",
                  afrag->iid, bfrag->iid);
        break;
      }
    }

    if (!olap_success) {
      fprintf(stderr, "MergeMultiAlignsFast_new()-- failed to find overlap between contigs %d and %d.\n",
              afrag->iid, bfrag->iid);
      DeleteMANode(manode->lid);
      safe_free(offsets);
      return(NULL);
    }

    if (otype == AS_CONTAINMENT) {
      bfrag->is_contained = 1;
      if (bfrag->container_iid == 0)
        bfrag->container_iid = 1;  //  Not sure why 1 and not afrag->iid
    }

    ApplyAlignment(afrag->lid, 0, NULL, bfrag->lid, ahang, Getint32(trace,0));
  } /* loop through all contigs */

  Delete_VA(trace);

  {
    IntMultiVar *vl;
    int32 nv;
    RefreshMANode(manode->lid, 0, opp, &nv, &vl, 0, 0);
  }

  // Now, want to generate a new MultiAlignT which merges the u_list and f_list of the contigs
  // merge the f_lists and u_lists by cloning and concating (or constructing dummy, when dealing with single read

  MultiAlignT *cma = CreateEmptyMultiAlignT();

  GetMANodeConsensus(manode->lid, cma->consensus, cma->quality);

  for (int32 i=0; i<num_contigs; i++) {
    MultiAlignT *ma = tigStore->loadMultiAlign(cpositions[i].ident, FALSE);
    ConcatVA_IntMultiPos(cma->f_list,  ma->f_list);
    ConcatVA_IntUnitigPos(cma->u_list, ma->u_list);
    ConcatVA_IntMultiVar(cma->v_list,  ma->v_list);
  }

  //  These are indices into the now merged cma->f_list and cma->u_list
  int32 ifrag   = 0;
  int32 iunitig = 0;

  for (int32 i=0; i<num_contigs; i++) {
    Fragment                 *cfrag      = GetFragment(fragmentStore,i);
    CNS_AlignedContigElement *components = GetCNS_AlignedContigElement(fragment_positions, cfrag->components);
    CNS_AlignedContigElement *compci;

    assert(cfrag->type == AS_CONTIG);

    // make adjustments to positions

    for (int32 ci=0; ci < cfrag->n_components; ci++) {
      compci = &components[ci];

#ifdef DEBUG_MERGEMULTIALIGNS
      if (compci->frg_or_utg==CNS_ELEMENT_IS_UNITIG)
        fprintf(stderr, "compci complement=%d length=%d bgn=%d end=%d\n", cfrag->complement, cfrag->length, compci->position.bgn, compci->position.end);
#endif

      int32 bgn = compci->position.bgn;
      int32 end = compci->position.end;

      if (cfrag->complement) {
        bgn = cfrag->length - compci->position.bgn;
        end = cfrag->length - compci->position.end;
      }

      int32 left  = (bgn < end) ? bgn : end;
      int32 right = (bgn < end) ? end : bgn;

      if (left < 0)
        left = 0;
      if (right > cfrag->length)
        right = cfrag->length;

#ifdef DEBUG_MERGEMULTIALIGNS
      if (compci->frg_or_utg==CNS_ELEMENT_IS_UNITIG)
        fprintf(stderr, "left=%d right=%d bgn=%d end=%d\n", left, right, bgn, end);
#endif

      left  = GetColumn(columnStore, GetBead(beadStore,cfrag->firstbead.get() + left)   ->column_index)->ma_index;
      right = GetColumn(columnStore, GetBead(beadStore,cfrag->firstbead.get() + right-1)->column_index)->ma_index +1;

      int32 tmp = bgn;

      bgn = (bgn < end) ? left  : right;
      end = (tmp < end) ? right : left;

#ifdef DEBUG_MERGEMULTIALIGNS
      if (compci->frg_or_utg==CNS_ELEMENT_IS_UNITIG)
        fprintf(stderr, "left=%d right=%d bgn=%d end=%d\n", left, right, bgn, end);
#endif

      if (compci->frg_or_utg==CNS_ELEMENT_IS_UNITIG) {
        IntUnitigPos *iup = GetIntUnitigPos(cma->u_list, iunitig);

#ifdef DEBUG_MERGEMULTIALIGNS
        fprintf(stderr, "Placing IUP  "F_CID" from "F_S32","F_S32" to "F_S32","F_S32" based on positions "F_S32","F_S32" (compl %d length %d within input parent)\n",
                iup->ident,
                iup->position.bgn, iup->position.end,
                bgn, end,
                compci->position.bgn, compci->position.end, cfrag->complement, cfrag->length);
#endif

        iup->position.bgn = bgn;
        iup->position.end = end;
        iup->delta_length = 0;
        iup->delta = NULL;

        iunitig++;
      } else {
        IntMultiPos *imp = GetIntMultiPos(cma->f_list, ifrag);
        assert(imp->ident == compci->idx.fragment.frgIdent);

        imp->ident = compci->idx.fragment.frgIdent;
        imp->position.bgn = bgn;
        imp->position.end = end;
        imp->delta_length = 0;
        imp->delta = NULL;

#ifdef DEBUG_MERGEMULTIALIGNS
        //  Generally not interesting.
        //fprintf(stderr, "Placing IMP1 "F_CID" at "F_S32","F_S32" based on positions "F_S32","F_S32" (compl %d length %d within input parent)\n",
        //        imp->ident, bgn, end, compci->position.bgn, compci->position.end, cfrag->complement, len);
#endif

        ifrag++;
      }
    }  //  over all components in the contig
  }  //  over all contigs

  DeleteMANode(manode->lid);
  safe_free(offsets);

  if (cma == NULL)
    fprintf(stderr, "MergeMultiAlignsFast_new() returns null cma?\n");

  return(cma);
}
