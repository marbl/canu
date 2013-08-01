
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


#define MAX_SURROGATE_FUDGE_FACTOR 5000


//  Probably should be listed with FragmentMap, but it's only used here.
HashTable_AS          *fragmentToIMP = NULL;


static
void
PlaceFragments(int32 fid,
               IntUnitigPos *aiup,
               CNS_Options  *opp) {

  Fragment                 *afrag = GetFragment(fragmentStore,fid);
  Fragment                 *bfrag = NULL;

  CNS_AlignedContigElement *belem = GetCNS_AlignedContigElement(fragment_positions, afrag->components);

  if (afrag->n_components == 0)
    return;

  VA_TYPE(int32) *trace = CreateVA_int32(AS_READ_MAX_NORMAL_LEN+1);

  for (; belem->frg_or_utg == CNS_ELEMENT_IS_FRAGMENT; belem++) {

    if (FALSE == ExistsInHashTable_AS(fragmentMap, belem->idx.fragment.frgIdent, 0))
      //  Fragment is not in the contigs f_list.  It is an unplaced read from a surrogate.
      continue;

    // if it exists in the fragmentMap it should exist in this map as well since it was added at the same time
    // look up where this fragment is placed within the entire contig, see if that matches where we're about to place it
    // this is necessary for surrogates that are multiply placed in a single contig for example:
    // contig: --------------*****--------*****--------> where ***** represents a surrogate
    //                        ---->                      readA
    // when placing readA within the surrogate unitig, we see if readA belongs in surrogate instance A or B
    // by computing the position of the unitig within the contig and adding ahang to it
    // if this computed position matches the position that the IMP record retrieved below tells us, proceed, otherwise skip placement
    IntMultiPos *bimp = (IntMultiPos *)LookupValueInHashTable_AS(fragmentToIMP, belem->idx.fragment.frgIdent, 0);
    int32  bbgn = (bimp->position.bgn < bimp->position.end ? bimp->position.bgn : bimp->position.end);
    int32  abgn = (aiup->position.bgn < aiup->position.end ? aiup->position.bgn : aiup->position.end);

    int32 fcomplement = afrag->complement;
    int32 bcomplement = (belem->position.bgn < belem->position.end) ? 0 : 1;

    int32           ahang = 0;
    int32           bhang = 0;
    int32           ovl   = 0;
    OverlapType   otype;

    //  all of fid's component frags will be aligned to it (not to
    //  each other)
    //
    //            fcomplement==0                                fcomplement==1
    //
    //        A)       fid                                  C)     fid
    //          ------------------>                            <----------------
    //          --->                                                        <---
    //           bid (bcomplement==0)                                       bid
    //
    //        B)       fid                                  D)     fid
    //          ------------------>                            <----------------
    //          <---                                                        --->
    //           bid (bcomplement==1)                                       bid
    //

    //  The afrag is a unitig, so b is always contained (length of
    //  overlap is length of belem).

    if        (fcomplement && bcomplement) {
      ahang = afrag->length - belem->position.bgn; /* Case D */
      bhang = belem->position.end - afrag->length;
      ovl   = belem->position.bgn - belem->position.end;
    } else if (fcomplement && !bcomplement) {
      ahang = afrag->length - belem->position.end; /* Case C */
      bhang = belem->position.bgn - afrag->length;
      ovl   = belem->position.end - belem->position.bgn;
    } else if (!fcomplement && bcomplement) {
      ahang = belem->position.end;                 /* Case B */
      bhang = belem->position.bgn - afrag->length;
      ovl   = belem->position.bgn - belem->position.end;
    } else {
      ahang = belem->position.bgn;                 /* Case A */
      bhang = belem->position.end - afrag->length;
      ovl   = belem->position.end - belem->position.bgn;
    }

    assert(ahang >= 0);
    assert(bhang <= 0);
    assert(ovl   >  0);

    if (aiup->num_instances > 1 && abs(ahang + abgn - bbgn) > MAX_SURROGATE_FUDGE_FACTOR) { 
      if (VERBOSE_MULTIALIGN_OUTPUT)
        fprintf(stderr, "Not placing fragment %d into unitig %d because the positions (%d, %d) do not match (%d, %d)\n",
                belem->idx.fragment.frgIdent, afrag->iid,
                bimp->position.bgn, bimp->position.end,
                ahang + GetColumn(columnStore,(GetBead(beadStore,afrag->firstbead.get()                ))->column_index)->ma_index,
                bhang + GetColumn(columnStore,(GetBead(beadStore,afrag->firstbead.get()+afrag->length-1))->column_index)->ma_index+1);
      continue;
    }

    int32 blid = AppendFragToLocalStore(belem->idx.fragment.frgType,
                                      belem->idx.fragment.frgIdent,
                                      (bcomplement != fcomplement),
                                      belem->idx.fragment.frgContained,
                                      AS_OTHER_UNITIG);

    afrag = GetFragment(fragmentStore, fid);  // AppendFragToLocalStore can change the pointer on us.
    bfrag = GetFragment(fragmentStore, blid);

    if (!GetAlignmentTraceDriver(afrag, NULL, bfrag, &ahang, &bhang, ovl, trace, &otype, GETALIGNTRACE_CONTIGF, 0)) {

      //if (!GetAlignmentTrace(afrag->lid, 0, blid, &ahang, &bhang, ovl, trace, &otype, DP_Compare,              DONT_SHOW_OLAP, 0, AS_CONSENSUS, AS_CNS_ERROR_RATE) &&
      //  !GetAlignmentTrace(afrag->lid, 0, blid, &ahang, &bhang, ovl, trace, &otype, Local_Overlap_AS_forCNS, DONT_SHOW_OLAP, 0, AS_CONSENSUS, AS_CNS_ERROR_RATE)) {

      Bead   *afirst = GetBead(beadStore, afrag->firstbead.get() + ahang);
      Column *col    = GetColumn(columnStore, afirst->column_index);
      MANode *manode = GetMANode(manodeStore, col->ma_id);

      RefreshMANode(manode->lid, 0, opp, NULL, NULL, 0, 0);  //  BPW not sure why we need this

      fprintf(stderr, "Could (really) not find overlap between %d (%c) and %d (%c) estimated ahang: %d (ejecting frag %d from contig)\n",
              afrag->iid, afrag->type, belem->idx.fragment.frgIdent, belem->idx.fragment.frgType, ahang, belem->idx.fragment.frgIdent);

      GetFragment(fragmentStore,blid)->deleted = 1;
    } else {
      ApplyAlignment(afrag->lid, 0, NULL, blid, ahang, Getint32(trace,0));
    }
  }  //  over all fragments

  Delete_VA(trace);
}


bool
MultiAlignContig(MultiAlignT  *ma,
                 gkStore      *UNUSED,
                 CNS_Options  *opp) {
  int32        num_bases     = 0;
  int32        num_unitigs   = GetNumIntUnitigPoss(ma->u_list);
  int32        num_frags     = GetNumIntMultiPoss(ma->f_list);
  int32        num_columns   = 0;

  IntMultiPos  *flist    = GetIntMultiPos(ma->f_list, 0);
  IntUnitigPos *ulist    = GetIntUnitigPos(ma->u_list, 0);
  IntMultiVar  *vlist    = GetIntMultiVar(ma->v_list, 0);

  SeqInterval  *offsets       = (SeqInterval *) safe_calloc(num_unitigs,sizeof(SeqInterval));

  for (int32 i=0;i<num_unitigs;i++) {
    int32 flen   = (ulist[i].position.bgn < ulist[i].position.end) ? (ulist[i].position.end < ulist[i].position.bgn) : (ulist[i].position.bgn - ulist[i].position.end);
    num_bases   += flen + 2 * AS_CNS_ERROR_RATE * flen;

    num_columns = (ulist[i].position.bgn > num_columns) ? ulist[i].position.bgn : num_columns;
    num_columns = (ulist[i].position.end > num_columns) ? ulist[i].position.end : num_columns;

    //fprintf(stderr, "CTG %d UTG %d %d-%d\n",
    //        ma->maID, ulist[i].ident, ulist[i].position.bgn, ulist[i].position.end);
  }

  for (int32 i=0;i<num_frags;i++) {
    int32 flen   = (flist[i].position.bgn < flist[i].position.end) ? (flist[i].position.end < flist[i].position.bgn) : (flist[i].position.bgn - flist[i].position.end);
    num_bases   += flen + 2 * AS_CNS_ERROR_RATE * flen;
  }

  ResetStores(num_bases, num_unitigs, num_columns);

  fragmentMap   = CreateScalarHashTable_AS();
  fragmentToIMP = CreateScalarHashTable_AS();

  for (int32 i=0; i<num_frags; i++) {

    //  Add all fragments in the contigs f_list to the fragmentMap.  This tells us if a fragment is
    //  not placed in a surrogate (because they aren't in the contigs f_list, but will appear in a
    //  surrogate unitigs f_list).
    //
    if (HASH_SUCCESS != InsertInHashTable_AS(fragmentMap, flist[i].ident, 0, 1, 0)) {
      fprintf(stderr, "MultiAlignContig()-- Contig %d FAILED.  Fragment %d is a duplicate.\n",
              ma->maID, flist[i].ident);
      return(false);
    }

    // SK store IID to IMP message mapping
    InsertInHashTable_AS(fragmentToIMP, flist[i].ident, 0, (uint64)&flist[i], 0);
  }

  for (int32 i=0;i<num_unitigs;i++) {
    uint32 complement = (ulist[i].position.bgn<ulist[i].position.end)?0:1;
    uint32 fid = AppendFragToLocalStore(AS_UNITIG,
                                 ulist[i].ident,
                                 complement,
                                 0,
                                 ulist[i].type);
    offsets[fid].bgn = complement?ulist[i].position.end:ulist[i].position.bgn;
    offsets[fid].end = complement?ulist[i].position.bgn:ulist[i].position.end;
  }

  MANode *manode = CreateMANode(ma->maID);

  // Seed multiAlignment with 1st fragment of 1st unitig

  SeedMAWithFragment(manode->lid, GetFragment(fragmentStore,0)->lid, opp);
  PlaceFragments(GetFragment(fragmentStore,0)->lid, ulist + GetFragment(fragmentStore,0)->lid, opp);

  // Now, loop on remaining fragments, aligning to:
  //    a)  containing frag (if contained)
  // or b)  previously aligned frag

  VA_TYPE(int32) *trace = CreateVA_int32(AS_READ_MAX_NORMAL_LEN+1);

  for (int32 i=1;i<num_unitigs;i++) {
    Fragment *afrag = NULL;
    Fragment *bfrag = GetFragment(fragmentStore,i);

    int32    ahang  = 0;
    int32    bhang  = 0;
    int32    ovl    = 0;
    int32    alid   = 0;
    int32    blid   = bfrag->lid;

    OverlapType otype;

    int32 olap_success  = 0;
    int32 try_contained = 0;
    int32 align_to      = i - 1;

    Fragment *afrag_first = NULL;
    int32       ahang_first = 0;
    int32       bhang_first = 0;

    while (!olap_success) {
    nextFrag:

      if (try_contained == 0)
        //  Skip contained stuff.
        while ((align_to > 0) &&
               ((GetFragment(fragmentStore, align_to)->is_contained) ||
                (GetFragment(fragmentStore, align_to)->container_iid > 0)))
          align_to--;

      if (align_to < 0) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignContig: hit the beginning of unitig list: no unitig upstream overlaps with current unitig %d\n", bfrag->iid);

        if (try_contained == 0) {
          if (VERBOSE_MULTIALIGN_OUTPUT)
            fprintf(stderr, "MultiAlignContig: trying contained afrags for bfrag %d\n", bfrag->iid);
          try_contained = 1;
          align_to      = i-1;
          goto nextFrag;
        }

        break;
      }

      afrag = GetFragment(fragmentStore, align_to);
      alid  = afrag->lid;

      ahang = offsets[blid].bgn - offsets[alid].bgn;
      bhang = offsets[blid].end - offsets[alid].end;

      if (afrag_first == NULL) {
        afrag_first = afrag;
        ahang_first = ahang;
        bhang_first = bhang;
      }

      //  This code copied from MultiAlignUnitig.

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

      if (ovl <= 0) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignContig: positions of afrag %d and bfrag %d do not overlap.  Proceed to the next upstream afrag\n", afrag->iid, bfrag->iid);
        align_to--;
        goto nextFrag;
      }

      olap_success = GetAlignmentTraceDriver(afrag, NULL,
                                             bfrag,
                                             &ahang, &bhang, ovl,
                                             trace,
                                             &otype,
                                             GETALIGNTRACE_CONTIGU,
                                             (blid + 1 < num_unitigs) ? (offsets[blid + 1].bgn - offsets[blid].bgn) : 800);

      //  Nope, fail.
      if (!olap_success) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignContig: Positions of afrag %d (%c) and bfrag %d (%c) overlap, but GetAlignmentTrace returns no overlap success.\n",
                  afrag->iid, afrag->type, bfrag->iid, bfrag->type);

        align_to--;

        if ((align_to < 0) && (!try_contained)) {
          if (VERBOSE_MULTIALIGN_OUTPUT)
            fprintf(stderr, "MultiAlignContig: Try contained afrags for bfrag %d\n", bfrag->iid);
          try_contained = 1;
          align_to = i-1;
        }
      }
    }  //  while !olap_success


    if ((!olap_success) && (FORCE_UNITIG_ABUT == 0)) {
      fprintf(stderr,"MultiAlignContig: Could (really) not find overlap between %d (%c) and %d (%c), estimated ahang %d\n",
              afrag->iid,afrag->type,bfrag->iid,bfrag->type, ahang);
      fprintf(stderr,"MultiAlignContig: You can (possibly) force these to abut with '-D forceunitigabut', but that code is buggy at best.\n");
      goto returnFailure;
    }

#if 1
    if ((!olap_success) && (FORCE_UNITIG_ABUT == 1)) {
      if (afrag_first) {
        afrag = afrag_first;
        ahang = ahang_first;
        bhang = bhang_first;
      } else {
        //  Dang, we're really screwed.  Nobody overlapped with us.
        //  Cross our fingers and find the closest end point.
        //
        int32   maxOvl = -offsets[blid].bgn;

        //if (VERBOSE_MULTIALIGN_OUTPUT)
        //  fprintf(stderr, "MultiAlignContig:  YIKES!  Your unitig doesn't overlap with anything!  Picking the closest thing!\n");

        align_to = i-1;

        while (align_to >= 0) {
          if ((try_contained == 0) &&
              ((GetFragment(fragmentStore, align_to)->is_contained) ||
               (GetFragment(fragmentStore, align_to)->container_iid > 0))) {
            //  NOP!  Found a contained frag, and we want to skip it.
          } else if (maxOvl < offsets[alid].end - offsets[blid].bgn) {
            afrag  = GetFragment(fragmentStore, align_to);
            alid   = afrag->lid;
            ahang  = offsets[blid].bgn - offsets[alid].bgn;
            maxOvl = offsets[alid].end - offsets[blid].bgn;

            //fprintf(stderr, "MultiAlignContig:  RESET align_to=%d alid=%d maxOvl=%d ahang=%d\n", align_to, alid, maxOvl, ahang);
          }

          align_to--;
        }  //  while align_to >= 0
      }

      fprintf(stderr, "MultiAlignContig:  Forcing abut between afrag %d (%c) and bfrag %d (%c) in contig %d.\n",
              afrag->iid, afrag->type, bfrag->iid, bfrag->type, ma->maID);

      //  Force a 1bp overlap.  We'd like to strictly abut, but ApplyAlignment() requires that there
      //  be an overlap, and removing checks for that seem like a bad idea.
      //
      ahang = afrag->length - 1;

      otype = AS_DOVETAIL;

      int32 zero = 0;

      ResetVA_int32(trace);
      AppendVA_int32(trace, &zero);

      assert(*Getint32(trace,0) == 0);
      assert(GetNumint32s(trace) == 1);
    }
#endif

    //  Unitig is placed, or we just forced it to be placed.

    if (otype == AS_CONTAINMENT) {
      bfrag->is_contained = 1;
      if (bfrag->container_iid == 0)
        bfrag->container_iid = 1;  //  Not sure why 1 and not afrag->iid
    }

    ApplyAlignment(afrag->lid, 0, NULL, bfrag->lid, ahang, Getint32(trace,0));
    PlaceFragments(bfrag->lid, ulist + bfrag->lid, opp);
  }  //  over all unitigs

  // Now, must find fragments in regions of overlapping unitigs, and adjust
  // their alignments as needed
  RefreshMANode(manode->lid, 0, opp, NULL, NULL, 0, 0);

  //fprintf(stderr,"MultiAlignContig: Initial pairwise induced alignment\n");
  //PrintAlignment(stderr,manode->lid,0,-1);

  AbacusRefine(manode,0,-1,CNS_SMOOTH, opp);
  MergeRefine(manode->lid, NULL, 0, opp, 1);
  AbacusRefine(manode,0,-1,CNS_POLYX, opp);

  //fprintf(stderr,"MultiAlignContig: POLYX refined alignment\n");
  //PrintAlignment(stderr,manode->lid,0,-1);

  {
    IntMultiVar  *vl = NULL;
    int32         nv = 0;

    RefreshMANode(manode->lid, 0, opp, &nv, &vl, 0, 0);
    AbacusRefine(manode,0,-1,CNS_INDEL, opp);
    MergeRefine(manode->lid, ma->v_list, 0, opp, 2);
  }

  //fprintf(stderr,"MultiAlignContig: Final refined alignment\n");
  //PrintAlignment(stderr,manode->lid,0,-1);

  //if (num_frags == 0)
  //  PrintAlignment(stderr,manode->lid,0,-1);

  GetMANodeConsensus(manode->lid, ma->consensus, ma->quality);
  GetMANodePositions(manode->lid, ma);

  DeleteMANode(manode->lid);

  safe_free(offsets);
  Delete_VA(trace);
  DeleteHashTable_AS(fragmentMap);  
  fragmentMap = NULL;
  DeleteHashTable_AS(fragmentToIMP);
  fragmentToIMP = NULL;
  return(true);

 returnFailure:
  safe_free(offsets);
  Delete_VA(trace);
  DeleteHashTable_AS(fragmentMap);
  fragmentMap = NULL;
  DeleteHashTable_AS(fragmentToIMP);
  fragmentToIMP = NULL;
  return(false);
}
