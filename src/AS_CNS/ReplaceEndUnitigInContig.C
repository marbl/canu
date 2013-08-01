
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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"
#include "MicroHetREZ.H"
#include "AS_UTL_reverseComplement.H"

#undef DEBUG_POSITIONS

static
void
PrintIMPInfo(FILE *print, int32 nfrags, IntMultiPos *imps) {
  uint32 bgn,end;
  for (int32 i=0;i<nfrags;i++) {
    bgn=imps->position.bgn;
    end=imps->position.end;
    if ( bgn < end )
      fprintf(print,"%12d F %c %10d, %10d -->\n",imps->ident,imps->type,bgn,end);
    else
      fprintf(print,"%12d F %c %10d, %10d <--\n",imps->ident,imps->type,end,bgn);
    imps++;
  }
}


static
void
PrintIUPInfo(FILE *print, int32 nfrags, IntUnitigPos *iups) {
  uint32 bgn,end;
  for (int32 i=0;i<nfrags;i++) {
    bgn=iups->position.bgn;
    end=iups->position.end;
    if ( bgn < end )
      fprintf(print,"%12d U %c %10d, %10d -->\n",iups->ident,iups->type,bgn,end);
    else
      fprintf(print,"%12d U %c %10d, %10d <--\n",iups->ident,iups->type,end,bgn);
    iups++;
  }
}


MultiAlignT *
ReplaceEndUnitigInContig(uint32 contig_iid,
                         uint32 unitig_iid,
                         int32 extendingLeft,
                         CNS_Options *opp) {
  int32 cid,tid; // local id of contig (cid), and unitig(tid)
  int32 aid,bid;
  int32 i,num_unitigs;
  MultiAlignT *oma;
  MultiAlignT *cma;
  IntUnitigPos *u_list;
  IntMultiPos *f_list;
  IntMultiVar  *v_list;
  int32 append_left=0;
  int32 num_frags=0;
  int32 complement=0;
  MANode *ma;
  Fragment *cfrag;
  Fragment *tfrag = NULL;
  static VA_TYPE(int32) *trace=NULL;

  oma =  tigStore->loadMultiAlign(contig_iid, FALSE);

  ResetStores(2 * GetNumchars(oma->consensus),
              2,
              2 * GetNumchars(oma->consensus));

  num_unitigs = GetNumIntUnitigPoss(oma->u_list);
  num_frags   = GetNumIntMultiPoss(oma->f_list);

  u_list = GetIntUnitigPos(oma->u_list,0);
  f_list = GetIntMultiPos(oma->f_list,0);
  v_list = GetIntMultiVar(oma->v_list,0);

  //PrintIMPInfo(stderr, num_frags,   f_list);
  //PrintIUPInfo(stderr, num_unitigs, u_list);

  // capture the consensus sequence of the original contig and put into local "fragment" format
  cid = AppendFragToLocalStore(AS_CONTIG,
                               contig_iid,
                               0,
                               0,
                               AS_OTHER_UNITIG);

  fprintf(stderr,"ReplaceEndUnitigInContig()-- contig %d unitig %d isLeft(%d)\n",
          contig_iid,unitig_iid,extendingLeft);

  //  The only real value-added from ReplaceUnitigInContig is a new consensus sequence for the contig
  //  some adjustments to positions go along with this, but the real compute is an alignment
  //  between the old contig consensus and the updated unitig
  //
  //  first we want to determine whether unitig is on left or right of contig,
  //  so that alignment can be done with a positive ahang
  //  if u is at left, i.e.:
  //
  //  C---------------C
  //  u------u
  //  then initialize new alignment with unitig, and add contig, else
  //
  //  if u is at right, i.e.:
  //
  //  C---------------C
  //           u------u
  //  then initialize new alignment with contig, and add unitig, else

  ma = CreateMANode(0);

  if ( trace == NULL )
    trace = CreateVA_int32(AS_READ_MAX_NORMAL_LEN);
  ResetVA_int32(trace);

  {
    int32 ahang,bhang,pos_offset=0;
    int32 tigs_adjusted_pos=0;
    OverlapType otype;
    int32 olap_success=0;
    cfrag=GetFragment(fragmentStore,cid);
    for(i=0;i<num_unitigs;i++) {
      uint32 id=u_list[i].ident;
      if ( id == unitig_iid ) {
        int32 bgn=u_list[i].position.bgn;
        int32 end=u_list[i].position.end;
        int32 complement_tmp=(bgn<end)?0:1;
        int32 left=(complement_tmp)?end:bgn;
        int32 right=(complement_tmp)?bgn:end;
        complement=complement_tmp;
        tid = AppendFragToLocalStore(AS_UNITIG,
                                     id,
                                     complement,
                                     0,
                                     AS_OTHER_UNITIG);
        tfrag=GetFragment(fragmentStore,tid);

        if ( extendingLeft ) {
          // need to set aid to unitig to preserve positive ahang -- and we now should always
          // have a bhang of zero.
          append_left=1;
          aid=tid;
          bid=cid;
          // and ahang estimate is the diff in size between
          // new unitig (GetFragment(fragmentStore,tid)->length) and old unitig (right-left)
          ahang = GetFragment(fragmentStore,tid)->length - (right-left);
          bhang = 0;
        }  else {
          //  --------
          //       ---+++
          //  We extended the unitig by "+++".  The ahang is just the
          //  start position of the original placement, and the bhang
          //  is the amount extended (as above).
          aid=cid;
          bid=tid;
          ahang = left;
          bhang = GetFragment(fragmentStore,tid)->length - (right-left);
        }
        SeedMAWithFragment(ma->lid, aid, opp);

        //  The expected length of this alignment is always the length of the original unitig.
        int32 ovl = right - left;

        olap_success = GetAlignmentTrace(aid, 0, bid, &ahang, &bhang, ovl, trace, &otype, DP_Compare, DONT_SHOW_OLAP, 0, GETALIGNTRACE_MERGE, AS_CGW_ERROR_RATE);

        if (!olap_success)
          olap_success = GetAlignmentTrace(aid, 0, bid, &ahang, &bhang, ovl, trace, &otype, Local_Overlap_AS_forCNS, DONT_SHOW_OLAP, 0, GETALIGNTRACE_MERGE, AS_CGW_ERROR_RATE);

        //  If the alignment fails -- usually because the ahang is
        //  negative -- return an empty alignment.  This causes
        //  extendClearRanges (the sole user of this function) to
        //  gracefully handle the failure.
        //
        if (olap_success == 0) {
          return(NULL);
          assert(olap_success);
        }

        ApplyAlignment(aid, 0, NULL, bid, ahang, Getint32(trace,0));
        RefreshMANode(ma->lid, 0, opp, NULL, NULL, 0, 0);

        //PrintAlignment(stderr,ma->lid,0,-1);

        break;
      }
    }
  }

  // Now, want to generate a new MultiAlignT which is an appropriate adjustment of original
  cma = CreateMultiAlignT();

  cma->maID = -1;
  cma->data = oma->data;

  cma->consensus = CreateVA_char(GetMANodeLength(ma->lid)+1);
  cma->quality   = CreateVA_char(GetMANodeLength(ma->lid)+1);

  GetMANodeConsensus(ma->lid, cma->consensus, cma->quality);

  // no deltas required at this stage
  // merge the f_lists and u_lists by cloning and concating

  cma->f_list = Clone_VA(oma->f_list);
  cma->u_list = Clone_VA(oma->u_list);
  cma->v_list = Clone_VA(oma->v_list);

  cma->fdelta = CreateVA_int32(0);
  cma->udelta = CreateVA_int32(0);

  {
    CNS_AlignedContigElement *components;
    CNS_AlignedContigElement *tcomponents;
    CNS_AlignedContigElement *contig_component;
    CNS_AlignedContigElement *aligned_component;
    int32 ifrag=0;
    int32 iunitig=0;
    IntMultiPos *imp;
    IntUnitigPos *iup;
    Fragment *frag;
    int32 ci=0;
    int32 tc=0; //unitig component index
    int32 bgn,end,left,right,tmp;
    int32 range_bgn=0,range_end=0,new_tig=0;
    components=GetCNS_AlignedContigElement(fragment_positions,cfrag->components);
    tcomponents=GetCNS_AlignedContigElement(fragment_positions,tfrag->components);
    // make adjustments to positions
    if ( append_left) {
      // fragments within unitig are 0 to tfrag->n_components
      // and cfrag->n_components-num_unitigs
      range_bgn = 0;
      range_end = tfrag->n_components-1;
      new_tig=cfrag->n_components-num_unitigs;
    } else {  // changed unitig on right
      // fragments within unitig are (num_frags-tfrag->n_components) to num_frags
      // and cfrag->n_components-1;
      range_bgn = (num_frags-(tfrag->n_components-1));
      range_end = num_frags;
      new_tig=cfrag->n_components-1;
    }
    while (ci < cfrag->n_components) {
      contig_component = &components[ci];
      if ( contig_component->frg_or_utg == CNS_ELEMENT_IS_FRAGMENT && contig_component->idx.fragment.frgInUnitig == unitig_iid ) {
        aligned_component = &tcomponents[tc++];
        if ( complement ) {
          bgn = tfrag->length-aligned_component->position.bgn;
          end = tfrag->length-aligned_component->position.end;
        } else {
          bgn = aligned_component->position.bgn;
          end = aligned_component->position.end;
        }
        frag = tfrag;
#ifdef DEBUG_POSITIONS
        fprintf(stderr,"compci->idx %12d bgn: %10d end: %10d\n",ci,bgn,end);
#endif
      } else if ( ci == new_tig ) {
        aligned_component =  &tcomponents[tc++];
        if ( complement ) {
          bgn = tfrag->length-aligned_component->position.bgn;
          end = tfrag->length-aligned_component->position.end;
        } else {
          bgn = aligned_component->position.bgn;
          end = aligned_component->position.end;
        }
        frag = tfrag;
#ifdef DEBUG_POSITIONS
        fprintf(stderr,"compci->idx %12d bgn: %10d end: %10d\n",ci,bgn,end);
#endif
      } else {
        aligned_component =  contig_component;
        bgn = aligned_component->position.bgn;
        end = aligned_component->position.end;
        frag = cfrag;
#ifdef DEBUG_POSITIONS
        fprintf(stderr,"compci->idx %12d bgn: %10d end: %10d\n",ci,bgn,end);
#endif
      }
      left = (bgn<end)?bgn:end;
      right = (bgn<end)?end:bgn;
      //if ( ci == new_tig ) {
      //    left = 0;
      //    right = frag->length;
      //}
      left = GetColumn(columnStore, GetBead(beadStore,frag->firstbead.get() + left   )->column_index)->ma_index;
      right= GetColumn(columnStore, GetBead(beadStore,frag->firstbead.get() + right-1)->column_index)->ma_index + 1;
      tmp = bgn;
      bgn = (bgn<end)?left:right;
      end = (tmp<end)?right:left;
      if (aligned_component->frg_or_utg==CNS_ELEMENT_IS_UNITIG) {
        iup = GetIntUnitigPos(cma->u_list,iunitig);
        iup->position.bgn = bgn;
        iup->position.end = end;
        iup->delta_length = 0;
        iup->delta = NULL;
#ifdef DEBUG_POSITIONS
        fprintf(stderr," element %d at %d,%d\n",
                ci,bgn,end);
#endif
        ci++;iunitig++;
      } else {
        imp = GetIntMultiPos(cma->f_list,ifrag);
        imp->ident = aligned_component->idx.fragment.frgIdent;
        imp->contained = aligned_component->idx.fragment.frgContained;
        imp->position.bgn = bgn;
        imp->position.end = end;
#ifdef DEBUG_POSITIONS
        fprintf(stderr," element %d at %d,%d\n", ci,bgn,end);
#endif
        imp->delta_length = 0;
        imp->delta = NULL;
        ci++;ifrag++;
      }
    }
  }
  DeleteMANode(ma->lid);
  return cma;
}
