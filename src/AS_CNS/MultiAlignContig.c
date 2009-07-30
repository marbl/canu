
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

static char *rcsid = "$Id: MultiAlignContig.c,v 1.4 2009-07-30 10:42:56 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"


#define MAX_SURROGATE_FUDGE_FACTOR 5000


//  Probably should be listed with FragmentMap, but it's only used here.
HashTable_AS          *fragmentToIMP = NULL;


static
void
PlaceFragments(int32 fid,
               IntUnitigPos *aiup,
               CNS_Options  *opp) {

  Fragment                 *afrag = GetFragment(fragmentStore,fid);
  CNS_AlignedContigElement *bfrag = GetCNS_AlignedContigElement(fragment_positions, afrag->components);

  if (afrag->n_components == 0)
    return;

  VA_TYPE(int32) *trace = CreateVA_int32(AS_READ_MAX_LEN+1);

  for (; bfrag->frg_or_utg == CNS_ELEMENT_IS_FRAGMENT; bfrag++) {
    if (!ExistsInHashTable_AS(fragmentMap, bfrag->idx.fragment.frgIdent, 0))
      continue;

    // if it exists in the fragmentMap it should exist in this map as well since it was added at the same time
    // look up where this fragment is placed within the entire contig, see if that matches where we're about to place it
    // this is necessary for surrogates that are multiply placed in a single contig for example:
    // contig: --------------*****--------*****--------> where ***** represents a surrogate
    //                        ---->                      readA
    // when placing readA within the surrogate unitig, we see if readA belongs in surrogate instance A or B
    // by computing the position of the unitig within the contig and adding ahang to it
    // if this computed position matches the position that the IMP record retrieved below tells us, proceed, otherwise skip placement
    IntMultiPos *bimp = (IntMultiPos *)LookupValueInHashTable_AS(fragmentToIMP, bfrag->idx.fragment.frgIdent, 0);
    int32  bbgn = (bimp->position.bgn < bimp->position.end ? bimp->position.bgn : bimp->position.end);
    int32  abgn = (aiup->position.bgn < aiup->position.end ? aiup->position.bgn : aiup->position.end);

    int fcomplement = afrag->complement;
    int bcomplement = (bfrag->position.bgn < bfrag->position.end) ? 0 : 1;

    int           ahang = 0;
    int           bhang = 0;
    int           ovl   = 0;
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
    //  overlap is length of bfrag).

    if        (fcomplement && bcomplement) {
      ahang = afrag->length - bfrag->position.bgn; /* Case D */
      bhang = bfrag->position.end - afrag->length;
      ovl   = bfrag->position.bgn - bfrag->position.end;
    } else if (fcomplement && !bcomplement) {
      ahang = afrag->length - bfrag->position.end; /* Case C */
      bhang = bfrag->position.bgn - afrag->length;
      ovl   = bfrag->position.end - bfrag->position.bgn;
    } else if (!fcomplement && bcomplement) {
      ahang = bfrag->position.end;                 /* Case B */
      bhang = bfrag->position.bgn - afrag->length;
      ovl   = bfrag->position.bgn - bfrag->position.end;
    } else {
      ahang = bfrag->position.bgn;                 /* Case A */
      bhang = bfrag->position.end - afrag->length;
      ovl   = bfrag->position.end - bfrag->position.bgn;
    }

    assert(ahang >= 0);
    assert(bhang <= 0);
    assert(ovl   >  0);

    if (aiup->num_instances > 1 && abs(ahang + abgn - bbgn) > MAX_SURROGATE_FUDGE_FACTOR) { 
      if (VERBOSE_MULTIALIGN_OUTPUT)
        fprintf(stderr, "Not placing fragment %d into unitig %d because the positions (%d, %d) do not match (%d, %d)\n",
                bfrag->idx.fragment.frgIdent, afrag->iid,
                bimp->position.bgn, bimp->position.end,
                ahang + GetColumn(columnStore,(GetBead(beadStore,afrag->firstbead))->column_index)->ma_index,
                bhang + GetColumn(columnStore,(GetBead(beadStore,afrag->firstbead+afrag->length-1))->column_index)->ma_index+1);
      continue;
    }

    int blid = AppendFragToLocalStore(bfrag->idx.fragment.frgType,
                                      bfrag->idx.fragment.frgIdent,
                                      (bcomplement != fcomplement),
                                      bfrag->idx.fragment.frgContained,
                                      AS_OTHER_UNITIG, NULL);

    afrag = GetFragment(fragmentStore, fid);  // AppendFragToLocalStore can change the pointer on us.

    if (!GetAlignmentTrace(afrag->lid, 0, blid, &ahang, &bhang, ovl, trace, &otype, DP_Compare,              DONT_SHOW_OLAP, 0, AS_CONSENSUS, AS_CNS_ERROR_RATE) &&
        !GetAlignmentTrace(afrag->lid, 0, blid, &ahang, &bhang, ovl, trace, &otype, Local_Overlap_AS_forCNS, DONT_SHOW_OLAP, 0, AS_CONSENSUS, AS_CNS_ERROR_RATE)) {

      Bead   *afirst = GetBead(beadStore, afrag->firstbead + ahang);
      Column *col    = GetColumn(columnStore, afirst->column_index);
      MANode *manode = GetMANode(manodeStore, col->ma_id);

      RefreshMANode(manode->lid, 0, opp, NULL, NULL, 0, 0);  //  BPW not sure why we need this

      fprintf(stderr, "Could (really) not find overlap between %d (%c) and %d (%c) estimated ahang: %d (ejecting frag %d from contig)\n",
              afrag->iid, afrag->type, bfrag->idx.fragment.frgIdent, bfrag->idx.fragment.frgType, ahang, bfrag->idx.fragment.frgIdent);

      GetFragment(fragmentStore,blid)->deleted = 1;
    } else {
      ApplyAlignment(afrag->lid, 0, NULL, blid, ahang, Getint32(trace,0));
    }
  }  //  over all fragments

  Delete_VA(trace);
}


int
MultiAlignContig(IntConConMesg *contig,
                 VA_TYPE(char) *sequence,
                 VA_TYPE(char) *quality,
                 VA_TYPE(int32) *deltas,
                 CNS_PrintKey printwhat,
                 CNS_Options *opp) {

  MANode       *ma=NULL;
  int           num_unitigs=0;
  int           num_frags=0;
  int           num_columns=0;
  int           complement=0;
  int           forced_contig=0;
  int           fid=0;
  int           i=0;
  int           align_to=0;
  IntUnitigPos *upositions=NULL;
  SeqInterval  *offsets=NULL;

  VA_TYPE(int32) *trace = NULL;

  IntMultiVar  *vl = NULL;
  int32         nv = 0;

  if (contig == NULL)
    return(FALSE);

  num_unitigs            = contig->num_unitigs;
  num_frags              = contig->num_pieces;
  upositions             = contig->unitigs;

  offsets = (SeqInterval *) safe_calloc(num_unitigs,sizeof(SeqInterval));
  for (i=0;i<num_unitigs;i++) {
    num_columns = (upositions[i].position.bgn>num_columns)? upositions[i].position.bgn : num_columns;
    num_columns = (upositions[i].position.end>num_columns)? upositions[i].position.end : num_columns;
  }

  ResetStores(num_unitigs, num_columns);

  fragmentMap   = CreateScalarHashTable_AS();
  fragmentToIMP = CreateScalarHashTable_AS();
  for (i=0;i<num_frags;i++) {
    if (ExistsInHashTable_AS (fragmentMap, contig->pieces[i].ident, 0)) {
      fprintf(stderr, "MultiAlignContig: Failure to insert ident %d in fragment hashtable, already present\n",contig->pieces[i].ident);
      //  If in a rush, just skip this frag, but you're better off fixing the input.
      assert(0);
    }

    //  The '1' value stored here is used by GetMANodePositions().
    InsertInHashTable_AS(fragmentMap, contig->pieces[i].ident, 0, 1, 0);
    
    // SK store IID to IMP message mapping
    InsertInHashTable_AS(fragmentToIMP, contig->pieces[i].ident, 0, (uint64)&contig->pieces[i], 0);
  }

  for (i=0;i<num_unitigs;i++) {
    complement = (upositions[i].position.bgn<upositions[i].position.end)?0:1;
    fid = AppendFragToLocalStore(AS_UNITIG,
                                 upositions[i].ident,
                                 complement,
                                 0,
                                 upositions[i].type, unitigStore);
    offsets[fid].bgn = complement?upositions[i].position.end:upositions[i].position.bgn;
    offsets[fid].end = complement?upositions[i].position.bgn:upositions[i].position.end;
  }

  if (DUMP_UNITIGS_IN_MULTIALIGNCONTIG > 0) {
    for (i=0; i<num_unitigs; i++) {
      Fragment *f = GetFragment(fragmentStore,i);
      char     *s = Getchar(sequenceStore,f->sequence);
      fprintf(stderr, ">unitig-%d\n%s\n", f->iid, s);
    }
  }

  ma = CreateMANode(contig->iaccession);

  trace = CreateVA_int32(AS_READ_MAX_LEN+1);

  // Seed multiAlignment with 1st fragment of 1st unitig

  SeedMAWithFragment(ma->lid, GetFragment(fragmentStore,0)->lid, 0, opp);
  PlaceFragments(GetFragment(fragmentStore,0)->lid, &contig->unitigs[GetFragment(fragmentStore,0)->lid], opp);

  // Now, loop on remaining fragments, aligning to:
  //    a)  containing frag (if contained)
  // or b)  previously aligned frag

  for (i=1;i<num_unitigs;i++) {
    Fragment *afrag = NULL;
    Fragment *bfrag = GetFragment(fragmentStore,i);

    int    ahang  = 0;
    int    bhang  = 0;
    int    ovl    = 0;
    int    alid   = 0;
    int    blid   = bfrag->lid;

    OverlapType otype;

    int olap_success  = 0;
    int try_contained = 0;
    int align_to      = i - 1;

    Fragment *afrag_first = NULL;
    int       ahang_first = 0;
    int       bhang_first = 0;

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
                                             'c',
                                             (blid + 1 < num_unitigs) ? (offsets[blid + 1].bgn - offsets[blid].bgn) : 800);

      //  Nope, fail.
      if (!olap_success) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignContig: Positions of afrag %d (%c) and bfrag %d (%c) overlap, but GetAlignmentTrace returns no overlap success.\n",
                  afrag->iid,afrag->type,bfrag->iid,bfrag->type, ahang);

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
      forced_contig = 1;

      if (afrag_first) {
        afrag = afrag_first;
        ahang = ahang_first;
        bhang = bhang_first;
      } else {
        //  Dang, we're really screwed.  Nobody overlapped with us.
        //  Cross our fingers and find the closest end point.
        //
        int   maxOvl = -offsets[blid].bgn;

        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignContig:  YIKES!  Your unitig doesn't overlap with anything!  Picking the closest thing!\n");

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

            fprintf(stderr, "MultiAlignContig:  RESET align_to=%d alid=%d maxOvl=%d ahang=%d\n", align_to, alid, maxOvl, ahang);
          }

          align_to--;
        }  //  while align_to >= 0
      }

      if (VERBOSE_MULTIALIGN_OUTPUT)
        fprintf(stderr, "MultiAlignContig:  Forcing abut between afrag %d (%c) and bfrag %d (%c).\n",
                afrag->iid, afrag->type, bfrag->iid, bfrag->type);

      //  If our ahang is too big, force a 20bp overlap.
      //
      if (ahang + 20 > afrag->length) {
        fprintf(stderr, "MultiAlignContig: RESET ahang from %d to %d\n", ahang, afrag->length-20);
        ahang = afrag->length - 20;
      }

      //  And now, after all this work, we probably will just die in
      //  ApplyAlignment() ten lines below.  Too bad.....

      otype = AS_DOVETAIL;
    }
#endif

    //  Unitig is placed, or we just forced it to be placed.

    if (otype == AS_CONTAINMENT) {
      bfrag->is_contained = 1;
      if (bfrag->container_iid == 0)
        bfrag->container_iid = 1;  //  Not sure why 1 and not afrag->iid
    }

    ApplyAlignment(afrag->lid, 0, NULL, bfrag->lid, ahang, Getint32(trace,0));
    PlaceFragments(bfrag->lid, &contig->unitigs[bfrag->lid], opp);
  }  //  over all unitigs

  // Now, must find fragments in regions of overlapping unitigs, and adjust
  // their alignments as needed
  RefreshMANode(ma->lid, 0, opp, NULL, NULL, 0, 0);

  if (printwhat == CNS_VERBOSE) {
    fprintf(stderr,"MultiAlignContig: Initial pairwise induced alignment\n");
    PrintAlignment(stderr,ma->lid,0,-1,printwhat);
  }

  AbacusRefine(ma,0,-1,CNS_SMOOTH, opp);
  MergeRefine(ma->lid, NULL, NULL, 0, opp, 1);
  AbacusRefine(ma,0,-1,CNS_POLYX, opp);

  if (printwhat == CNS_VERBOSE) {
    fprintf(stderr,"MultiAlignContig: POLYX refined alignment\n");
    PrintAlignment(stderr,ma->lid,0,-1,printwhat);
  }

  RefreshMANode(ma->lid, 0, opp, &nv, &vl, 0, 0);
  AbacusRefine(ma,0,-1,CNS_INDEL, opp);
  MergeRefine(ma->lid, &(contig->v_list), &(contig->num_vars), 0, opp, 2);

  if ((printwhat == CNS_VERBOSE) ||
      (printwhat == CNS_VIEW_UNITIG)) {
    fprintf(stderr,"MultiAlignContig: Final refined alignment\n");
    PrintAlignment(stderr,ma->lid,0,-1,printwhat);
  }

  if (num_frags == 0)
    PrintAlignment(stderr,ma->lid,0,-1,printwhat);

  GetMANodeConsensus(ma->lid,sequence,quality);
  contig->consensus  = Getchar(sequence,0);
  contig->quality    = Getchar(quality,0);
  contig->num_pieces = GetMANodePositions(ma->lid, num_frags, contig->pieces, num_unitigs, contig->unitigs, deltas);
  contig->length     = GetNumchars(sequence)-1;
  contig->forced     = forced_contig;

  DeleteMANode(ma->lid);

  safe_free(offsets);
  Delete_VA(trace);
  DeleteHashTable_AS(fragmentMap);  
  fragmentMap = NULL;
  DeleteHashTable_AS(fragmentToIMP);
  fragmentToIMP = NULL;
  return(TRUE);

 returnFailure:
  safe_free(offsets);
  Delete_VA(trace);
  DeleteHashTable_AS(fragmentMap);
  fragmentMap = NULL;
  DeleteHashTable_AS(fragmentToIMP);
  fragmentToIMP = NULL;
  return(FALSE);
}
