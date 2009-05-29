
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

static char *rcsid = "$Id: MultiAlignUnitig.c,v 1.1 2009-05-29 17:29:19 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"



static
int
MANode2Array(MANode *ma, int *depth, char ***array, int ***id_array,
                 int show_cel_status) {
  char **multia;
  int **ia;
  int length = GetNumColumns(ma->columns);
  // find max column depth.
  int max_depth=0;
  int col_depth;
  int column_index;
  Column *col;
  char laneformat[40];
  int num_frags=GetNumFragments(fragmentStore);
  Fragment *frag;
  int fid;
  int *rowptr,*row_assign;
  int ir,fbgn,fend;
  int i;
  *depth =  0;
  for (column_index = ma->first;column_index != -1;  ) {
    col = GetColumn(columnStore, column_index);
    if ( col != NULL ) {
      col_depth = GetDepth(col);
      max_depth = (col_depth > max_depth)?col_depth:max_depth;
    }
    column_index = col->next;
  }
  *depth = 2*max_depth; // rough estimate. first pack rows, then adjust to actual consumed rows
  rowptr = (int *)safe_malloc((*depth)*sizeof(int));
  row_assign = (int *)safe_malloc(num_frags*sizeof(int));
  for (ir=0;ir<*depth;ir++) rowptr[ir] = 0;
  for (ir=0;ir<num_frags;ir++) row_assign[ir] = -1;
  frag = GetFragment(fragmentStore,0);
  // setup the packing
  for ( fid=0;fid<num_frags;fid++ ) {
    if ( frag->type != AS_UNITIG ) {
      fbgn = GetColumn(columnStore,(GetBead(beadStore,frag->firstbead))->column_index)->ma_index;
      fend = GetColumn(columnStore,
                       (GetBead(beadStore,frag->firstbead+frag->length-1))->column_index)->ma_index+1;
      for (ir=0;ir<*depth;ir++) {
        if (fbgn <  rowptr[ir] ) continue;
        rowptr[ir] = fend;
        row_assign[fid] = ir;
        break;
      }
      if (row_assign[fid] <= -1)
        {
          *depth += max_depth;
          rowptr = (int *)safe_realloc(rowptr, (*depth)*sizeof(int));
          fid--;
          continue;
        }
    }
    frag++;
  }
  // now, find out actual depth
  max_depth = 0;
  for (ir=0;ir<*depth;ir++) {
    if (rowptr[ir] == 0 ) {
      max_depth = ir+1;
      break;
    }
  }
  if ( max_depth == 0 ) max_depth = ir;
  *depth = max_depth;
  multia = (char **)safe_malloc(2*(*depth)*sizeof(char *));
  ia = (int **)safe_malloc((*depth)*sizeof(int *));
  sprintf(laneformat,"%%%ds",length);
  {int j;
    for (i=0;i<(*depth);i++) {
      ia[i] = (int *) safe_malloc( length*sizeof(int));
      for (j=0;j<length;j++) ia[i][j] = 0;
    }
  }
  for (i=0;i<2*(*depth);i++) {
    multia[i] = (char *) safe_malloc((length+1)*sizeof(char));
    sprintf(multia[i],laneformat," ");
    *(multia[i]+length) = '\0';
  }
  {
    Bead *fb;
    FragmentBeadIterator fi;
    int bid;
    char bc,bq;
    Column *bcolumn;
    int ma_index;

    frag = GetFragment(fragmentStore,0);
    for ( fid=0;fid<num_frags;fid++ ) {
      if ( frag->type != AS_UNITIG ) {
        ir = row_assign[fid];
        fb = GetBead(beadStore,frag->firstbead);
        bcolumn =  GetColumn(columnStore,fb->column_index);

        CreateFragmentBeadIterator(fid,&fi);

        while ( (bid = NextFragmentBead(&fi)) != -1 ) {
          fb = GetBead(beadStore,bid);
          bc = *Getchar(sequenceStore,fb->soffset);
          bq = *Getchar(qualityStore,fb->soffset);
          bcolumn =  GetColumn(columnStore,fb->column_index);
          ma_index = bcolumn->ma_index;
          // find the first open row here, and put in the sequence/quality/ident
          multia[2*ir][ma_index] = bc;
          multia[2*ir+1][ma_index] = bq;
          ia[ir][ma_index] = frag->iid;
        }
      }
      frag++;
    }
  }
  *array = multia;
  *id_array = ia;
  safe_free(rowptr);
  safe_free(row_assign);
  return 1;
}



int
MultiAlignUnitig(IntUnitigMesg   *unitig,
                 GateKeeperStore *fragStore,
                 VA_TYPE(char)   *sequence,
                 VA_TYPE(char)   *quality,
                 VA_TYPE(int32)  *deltas,
                 CNS_PrintKey     printwhat,
                 CNS_Options     *opp) {

  int32 i;

  VA_TYPE(int32) *trace   = NULL;
  MANode         *ma      = NULL;
  SeqInterval    *offsets = NULL;

  gkpStore          = fragStore;

  {
    int32 num_columns = 0;

    for (i=0;i<unitig->num_frags;i++) {
      num_columns = (unitig->f_list[i].position.bgn > num_columns) ? unitig->f_list[i].position.bgn : num_columns;
      num_columns = (unitig->f_list[i].position.end > num_columns) ? unitig->f_list[i].position.end : num_columns;
    }

    ResetStores(unitig->num_frags, num_columns);
  }

  fragmentMap = CreateScalarHashTable_AS();

  //  Magic initialization prevents us calling CreateMANode() until now.

  trace   = CreateVA_int32(2 * AS_READ_MAX_LEN);
  ma      = CreateMANode(unitig->iaccession);
  offsets = (SeqInterval *)safe_calloc(unitig->num_frags, sizeof(SeqInterval));

  assert(ma->lid == 0);

  for (i=0; i<unitig->num_frags; i++) {
    int complement = (unitig->f_list[i].position.bgn < unitig->f_list[i].position.end) ? 0 : 1;
    int fid;

    if (unitig->f_list[i].type != AS_READ) {
      fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Fragment %d is not a read.\n",
              unitig->iaccession, unitig->f_list[i].ident);
      goto MAUnitigFailure;
    }

    if (HASH_SUCCESS != InsertInHashTable_AS(fragmentMap,unitig->f_list[i].ident, 0, 1, 0)) {
      fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Fragment %d is a duplicate.\n",
              unitig->iaccession, unitig->f_list[i].ident);
      goto MAUnitigFailure;
    }

    // This guy allocates and initializes the beads for each
    // fragment.  Beads are not fully inserted in the abacus here.

    fid = AppendFragToLocalStore(unitig->f_list[i].type,
                                 unitig->f_list[i].ident,
                                 complement,
                                 unitig->f_list[i].contained,
                                 AS_OTHER_UNITIG, NULL);

    offsets[fid].bgn = complement ? unitig->f_list[i].position.end : unitig->f_list[i].position.bgn;
    offsets[fid].end = complement ? unitig->f_list[i].position.bgn : unitig->f_list[i].position.end;

    //  If this is violated, then align_to is not a valid index into unitig->f_list.
    assert(fid == i);

    if (VERBOSE_MULTIALIGN_OUTPUT)
      fprintf(stderr,"MultiAlignUnitig()-- Added fragment %d (%d-%d) in unitig %d to store at local id %d.\n",
              unitig->f_list[i].ident, unitig->f_list[i].position.bgn, unitig->f_list[i].position.end, unitig->iaccession, fid);
  }

  SeedMAWithFragment(ma->lid, GetFragment(fragmentStore,0)->lid,0, opp);

  // Now, loop on remaining fragments, aligning to:
  //    a)  containing frag (if contained)
  // or b)  previously aligned frag

  for (i=1; i<unitig->num_frags; i++) {
    int         align_to     = i-1;
    Fragment   *afrag        = NULL;                           //  Last frag
    Fragment   *bfrag        = GetFragment(fragmentStore, i);  //  This frag

    int         ahang        = 0;
    int         bhang        = 0;

    int         olap_success = 0;
    OverlapType otype        = 0;

    int         ovl          = 0;

    if (VERBOSE_MULTIALIGN_OUTPUT) {
      fprintf(stderr, "\n");
      fprintf(stderr, "MultiAlignUnitig()-- processing fragment %d iid=%d pos=%d,%d parent=%d,%d,%d contained=%d\n",
              i,
              unitig->f_list[i].ident,
              unitig->f_list[i].position.bgn,
              unitig->f_list[i].position.end,
              unitig->f_list[i].parent,
              unitig->f_list[i].ahang,
              unitig->f_list[i].bhang,
              unitig->f_list[i].contained);
    }

    //  Scan for the thickest overlap to an existing fragment.  If our
    //  parent is NOT this fragment, we may be in for alignment
    //  troubles.
    //


    //  If we have a parent, assume the hangs are correct and just
    //  align to it.  If that works, we're done.
    //
    //  But first, we need to check that the parent is in fact our
    //  thickest overlap.  If not, we'll introduce too many gaps in
    //  the consensus to use this method.
    //
    if (unitig->f_list[i].parent > 0) {
      int  numDovetail = 0;

      //  Search for the parent fragment -- it MUST be the first
      //  non-contained we see, otherwise, there is a thicker overlap
      //  we should be using.
      //
      for (align_to = i-1; align_to >= 0; align_to--) {
        afrag = GetFragment(fragmentStore, align_to);

        //  Keep track of the number of non-contained fragments we've
        //  encountered.
        if (unitig->f_list[align_to].contained == 0)
          numDovetail++;

        //  Did we find the parent?
        if (unitig->f_list[i].parent == afrag->iid) {
          ahang = unitig->f_list[i].ahang;
          bhang = unitig->f_list[i].bhang;

          //  Here, we trust the hangs more than the placement, but we
          //  cannot compute the overlap length from just the hangs (OK,
          //  we could, using the fragment length, I guess).

          if (ahang > 0)
            if (bhang > 0)
              //  normal dovetail
              ovl = offsets[afrag->lid].end - offsets[bfrag->lid].bgn;
            else
              //  b is contained in a
              //ovl = offsets[bfrag->lid].end - offsets[bfrag->lid].bgn;
              ovl = bfrag->length;
          else
            //  We don't expect these to occur.
            if (bhang > 0)
              //  a is contained in b
              //ovl = offsets[afrag->lid].end - offsets[afrag->lid].bgn;
              ovl = afrag->length;
            else
              //  anti normal dovetail
              ovl = offsets[bfrag->lid].end - offsets[afrag->lid].bgn;

          //  If this is the first non-contained fragment we've
          //  encountered, OR the fragment is contained, try the
          //  alignment.
          //
          if ((numDovetail == 1) || (unitig->f_list[i].contained == afrag->iid)) {
            if (VERBOSE_MULTIALIGN_OUTPUT)
              fprintf(stderr, "MultiAlignUnitig()-- (par) ovl=%d (afrag: lid=%d %d-%d  bfrag: lid=%d %d-%d  hangs: %d %d\n",
                      ovl,
                      afrag->lid, offsets[afrag->lid].bgn, offsets[afrag->lid].end,
                      bfrag->lid, offsets[bfrag->lid].bgn, offsets[bfrag->lid].end,
                      ahang, bhang);

            olap_success = GetAlignmentTraceDriver(afrag, 0,
                                                   bfrag,
                                                   &ahang, &bhang, ovl,
                                                   trace,
                                                   &otype,
                                                   'u',
                                                   0);
          } else {
            if (VERBOSE_MULTIALIGN_OUTPUT)
              fprintf(stderr, "MultiAlignUnitig()-- (par) ovl=%d (afrag: lid=%d %d-%d  bfrag: lid=%d %d-%d  hangs: %d %d  NOT FIRST DOVETAIL, skip\n",
                      ovl,
                      afrag->lid, offsets[afrag->lid].bgn, offsets[afrag->lid].end,
                      bfrag->lid, offsets[bfrag->lid].bgn, offsets[bfrag->lid].end,
                      ahang, bhang);
          }

          //  We're done with the loop now.  Continue; if we found an
          //  overlap, we skip the while loop below and process this
          //  overlap.

          break;
        }  //  Found the parent
      }  //  Over all potential parent fragments
    }  //  No parent


    //  Start over.
    align_to = i-1;

    while (! olap_success) {
      int         allow_neg_hang_once = 0;
      int         ahangvalid   = 0;

      if (align_to < 0) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignUnitig()-- hit the beginning of fragment list: no fragment upstream overlaps with current fragment %d\n",
                  bfrag->iid);
        break;
      }


      //  Grab the candidate fragment.
      //
      afrag = GetFragment(fragmentStore, align_to);


      //  If there is an overlap store, and our parent isn't this fragment, query
      //  the store to find the overlap to get the correct hangs.
      //
      if ((ahangvalid == 0) &&
          (unitig->f_list[i].parent != afrag->iid) &&
          (ovlStore)) {
        OVSoverlap  ovsOverlap;
        int         foundOverlap = 0;

        fprintf(stderr, "get afrag=%d from store.  bfrag=%d\n", afrag->iid, bfrag->iid);

        AS_OVS_setRangeOverlapStore(ovlStore, afrag->iid, afrag->iid);

        while ((ahangvalid == 0) &&
               (AS_OVS_readOverlapFromStore(ovlStore, &ovsOverlap, AS_OVS_TYPE_OVL))) {
          if (ovsOverlap.b_iid != bfrag->iid)
            continue;

          if (unitig->f_list[align_to].position.bgn < unitig->f_list[align_to].position.end) {
            ahang = ovsOverlap.dat.ovl.a_hang;
            bhang = ovsOverlap.dat.ovl.b_hang;
          } else {
            ahang = -ovsOverlap.dat.ovl.b_hang;
            bhang = -ovsOverlap.dat.ovl.a_hang;
          }

          if (ahang > 0)
            if (bhang > 0)
              //  normal dovetail
              ovl = offsets[afrag->lid].end - offsets[bfrag->lid].bgn;
            else
              //  b is contained in a
              //ovl = offsets[bfrag->lid].end - offsets[bfrag->lid].bgn;
              ovl = bfrag->length;
          else
            //  We don't expect these to occur.
            if (bhang > 0)
              //  a is contained in b
              //ovl = offsets[afrag->lid].end - offsets[afrag->lid].bgn;
              ovl = afrag->length;
            else
              //  anti normal dovetail
              ovl = offsets[bfrag->lid].end - offsets[afrag->lid].bgn;

          if (VERBOSE_MULTIALIGN_OUTPUT)
            fprintf(stderr, "MultiAlignUnitig()-- (ovs) ovl=%d (afrag: lid=%d %d-%d  bfrag: lid=%d %d-%d  hangs: %d %d\n",
                    ovl,
                    afrag->lid, offsets[afrag->lid].bgn, offsets[afrag->lid].end,
                    bfrag->lid, offsets[bfrag->lid].bgn, offsets[bfrag->lid].end,
                    ahang, bhang);

          //  SPECIAL CASE;  This is a REAL OVERLAP, so allow the negative hang.
          if (ahang < 0)
            allow_neg_hang_once = TRUE;

          ahangvalid = 3;
        }
      }

      //  If the parent is defined, then unitigger/bog should have
      //  filled in the correct ahang and bhang.
      //
      if ((ahangvalid == 0) &&
          (unitig->f_list[i].parent > 0) &&
          (unitig->f_list[i].parent == afrag->iid)) {
        ahang = unitig->f_list[i].ahang;
        bhang = unitig->f_list[i].bhang;

        //  Here, we trust the hangs more than the placement, but we
        //  cannot compute the overlap length from just the hangs (OK,
        //  we could, using the fragment length, I guess).

        if (ahang > 0)
          if (bhang > 0)
            //  normal dovetail
            ovl = offsets[afrag->lid].end - offsets[bfrag->lid].bgn;
          else
            //  b is contained in a
            //ovl = offsets[bfrag->lid].end - offsets[bfrag->lid].bgn;
            ovl = bfrag->length;
        else
          //  We don't expect these to occur.
          if (bhang > 0)
            //  a is contained in b
            //ovl = offsets[afrag->lid].end - offsets[afrag->lid].bgn;
            ovl = afrag->length;
          else
            //  anti normal dovetail
            ovl = offsets[bfrag->lid].end - offsets[afrag->lid].bgn;

        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignUnitig()-- (new) ovl=%d (afrag: lid=%d %d-%d  bfrag: lid=%d %d-%d  hangs: %d %d\n",
                  ovl,
                  afrag->lid, offsets[afrag->lid].bgn, offsets[afrag->lid].end,
                  bfrag->lid, offsets[bfrag->lid].bgn, offsets[bfrag->lid].end,
                  ahang, bhang);

        ahangvalid = 2;
      }

      //  Otherwise, fallback to the layout to guesstimate the hangs.
      //
      if (ahangvalid == 0) {
        //  This code copied into MergeMultiAligns & MultiAlignContig
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
          fprintf(stderr, "MultiAlignUnitig()-- (old) ovl=%d (afrag: lid=%d %d-%d  bfrag: lid=%d %d-%d  hangs: %d %d\n",
                  ovl,
                  afrag->lid, offsets[afrag->lid].bgn, offsets[afrag->lid].end,
                  bfrag->lid, offsets[bfrag->lid].bgn, offsets[bfrag->lid].end,
                  ahang, bhang);

        ahangvalid = 1;
      }


      //  If there is no overlap in the layout, fail.
      //
      if (ovl <= 0) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignUnitig()-- positions of afrag %d and bfrag %d do not overlap; proceed to the next upstream afrag\n",
                  afrag->iid, bfrag->iid);
        align_to--;
        continue;
      }


      // Make sure ahang is above the cutoff value.  If it's not, may
      // need to sort fragments begfore processing.
      //
      if ((ahang < CNS_NEG_AHANG_CUTOFF) && (allow_neg_hang == FALSE) && (allow_neg_hang_once == FALSE)) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignUnitig()-- too negative ahang is detected for afrag %d and bfrag %d; proceed to the next upstraem afrag\n",
                  afrag->iid,  bfrag->iid);
        align_to--;
        continue;
      }


      //  If the bfrag is marked as contained, require that it be
      //  contained in the overlap.  This test is after the negative
      //  ahang test, and we allow whatever ahang passed.  We only
      //  need to check that the bhang indicates containment (and that
      //  we're not to the expected container).
      //
      if ((bfrag->container_iid > 0) &&
          (bhang > 0) &&
          (bfrag->container_iid != afrag->iid)){
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignUnitig()-- afrag %d is not container of contained bfrag %d; proceed to the next upstraem afrag\n",
                  afrag->iid, bfrag->iid);
        align_to--;
        continue;
      }


      //  If the bfrag is not contained, prevent alignments to
      //  contained fragments, unless that is our parent.
      //
      if ((bfrag->container_iid == 0) &&
          (afrag->container_iid  > 0) &&
          ((allow_contained_parent == 0) || (unitig->f_list[i].parent != afrag->iid))) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignUnitig()-- afrag %d is contained, and not parent of uncontained bfrag %d; proceed to the next upstraem afrag\n",
                  afrag->iid, bfrag->iid);
        align_to--;
        continue;
      }


      if (VERBOSE_MULTIALIGN_OUTPUT)
        fprintf(stderr, "MultiAlignUnitig()-- Aligning frag #%d (iid %d, range %d,%d) to afrag iid %d range %d,%d -- ovl=%d ahang=%d\n",
                unitig->f_list[i].ident,
                bfrag->iid,
                offsets[bfrag->lid].bgn,
                offsets[bfrag->lid].end,
                afrag->iid,
                offsets[afrag->lid].bgn,
                offsets[afrag->lid].end,
                ovl,
                ahang);


      olap_success = GetAlignmentTraceDriver(afrag, 0,
                                             bfrag,
                                             &ahang, &bhang, ovl,
                                             trace,
                                             &otype,
                                             'u',
                                             0);

      if (!olap_success) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr, "MultiAlignUnitig()-- positions of bfrag %d (%c) and afrag %d (%c) overlap, but no overlap found; estimated ahang: %d%s\n",
                  bfrag->iid, bfrag->type,
                  afrag->iid, afrag->type,
                  ahang,
                  (bfrag->container_iid)?" (reported contained)":"");
        align_to--;
      }
    }   //  while (!olap_success)


    if (!olap_success) {
      fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Overlap not found between %sbfrag %d (%c) and afrag %d (%c)%s.\n",
              unitig->iaccession,
              (bfrag->container_iid == 0) ? "" : "contained ",
              bfrag->iid, bfrag->type,
              afrag->iid, afrag->type,
              (bfrag->container_iid == 0) ? "" : ((bfrag->container_iid == afrag->iid) ? ", the container" : ", not the container"));

      goto MAUnitigFailure;
    }


    // Without this marking, the multialignment is likely to have
    // pieces which are not properly aligned, and which will appear as
    // block indels (large gap-blocks) which will foil future overlaps
    // involving the consensus sequence of this "reformed" unitig
    //
    if (otype == AS_CONTAINMENT) {
      bfrag->is_contained = 1;
      if (bfrag->container_iid == 0)
        bfrag->container_iid = 1;  //  Not sure why 1 and not afrag->iid
    }

    if (VERBOSE_MULTIALIGN_OUTPUT)
      fprintf(stderr, "MultiAlignUnitig()--  bfrag %d (%c) is %s in afrag %d (%c).\n",
              bfrag->iid, bfrag->type,
              (otype == AS_CONTAINMENT) ? "contained" : "not contained",
              afrag->iid, afrag->type);


    //  Update parent and hangs to reflect the overlap that succeeded.
    //
    if ((unitig->f_list[i].parent != afrag->iid) ||
        (unitig->f_list[i].ahang  != ahang) ||
        (unitig->f_list[i].bhang  != bhang) ||
        ((otype == AS_CONTAINMENT) && (unitig->f_list[i].contained != afrag->iid))) {
      if (unitig->f_list[i].parent)
        fprintf(stderr, "MultiAlignUnitig()-- update parent from id=%d %d,%d to id=%d %d,%d\n",
                unitig->f_list[i].parent, unitig->f_list[i].ahang, unitig->f_list[i].bhang,
                afrag->iid, ahang, bhang);

      unitig->f_list[i].parent = afrag->iid;
      unitig->f_list[i].ahang  = ahang;
      unitig->f_list[i].bhang  = bhang;

      if (otype == AS_CONTAINMENT)
        unitig->f_list[i].contained = afrag->iid;
    }

    ApplyAlignment(afrag->lid, 0, bfrag->lid, ahang, Getint32(trace, 0));
  }  //  over all frags

  RefreshMANode(ma->lid, 0, opp, NULL, NULL, 0, 0);

  AbacusRefine(ma,0,-1,CNS_SMOOTH, opp);
  MergeRefine(ma->lid, NULL, NULL, 1, opp, 1);

  AbacusRefine(ma,0,-1,CNS_POLYX, opp);
  MergeRefine(ma->lid, NULL, NULL, 1, opp, 1);

  AbacusRefine(ma,0,-1,CNS_INDEL, opp);
  MergeRefine(ma->lid, NULL, NULL, 1, opp, 1);

  GetMANodeConsensus(ma->lid,sequence,quality);
  GetMANodePositions(ma->lid, unitig->num_frags, unitig->f_list, 0, NULL, deltas);

  unitig->consensus = Getchar(sequence,0);
  unitig->quality   = Getchar(quality,0);
  unitig->length    = GetNumchars(sequence)-1;


  if ((printwhat == CNS_VERBOSE) ||
      (printwhat == CNS_VIEW_UNITIG))
    PrintAlignment(stderr,ma->lid,0,-1,printwhat);


  //  Fix up our containments.  Sometimes consensus will move a
  //  fragment out of the containment relationship, and this causes
  //  headaches later on.
  //
  {
    int   frag = 0;
    int   cntr = 0;

    for (frag=1; frag < unitig->num_frags; frag++) {
      if (unitig->f_list[frag].contained) {
        for (cntr=frag-1; cntr>=0; cntr--) {
          if (unitig->f_list[frag].contained == unitig->f_list[cntr].ident) {
            int cbeg = MIN(unitig->f_list[cntr].position.bgn, unitig->f_list[cntr].position.end);
            int fbeg = MIN(unitig->f_list[frag].position.bgn, unitig->f_list[frag].position.end);

            int cend = MAX(unitig->f_list[cntr].position.bgn, unitig->f_list[cntr].position.end);
            int fend = MAX(unitig->f_list[frag].position.bgn, unitig->f_list[frag].position.end);

            if ((fbeg < cbeg) || (cend > cend)) {
              //fprintf(stderr, "WARNING: FRAG %d (%d,%d) is NOT contained in %d (%d,%d)\n",
              //        unitig->f_list[frag].ident, unitig->f_list[frag].position.bgn, unitig->f_list[frag].position.end,
              //        unitig->f_list[cntr].ident, unitig->f_list[cntr].position.bgn, unitig->f_list[cntr].position.end);

              unitig->f_list[frag].contained = 0;
            }
          }
        }
      }
    }
  }



  //  While we have fragments in memory, compute the microhet
  //  probability.  Ideally, this would be done in CGW when loading
  //  unitigs (the only place the probability is used) but the code
  //  wants to load sequence and quality for every fragment, and
  //  that's too expensive.
  {
    int    depth  = 0;
    char **multia = NULL;
    int  **id_array = NULL;

    MANode2Array(ma, &depth, &multia, &id_array,0);

    unitig->microhet_prob = AS_REZ_MP_MicroHet_prob(multia, id_array, gkpStore, unitig->length, depth);

    for (i=0;i<depth;i++) {
      safe_free(multia[2*i]);
      safe_free(multia[2*i+1]);
      safe_free(id_array[i]);
    }
    safe_free(multia);
    safe_free(id_array);
  }

 MAUnitigSuccess:
  DeleteVA_int32(trace);
  DeleteHashTable_AS(fragmentMap);  fragmentMap = NULL;
  DeleteMANode(ma->lid);

  safe_free(offsets);

  return(TRUE);

 MAUnitigFailure:
  DeleteVA_int32(trace);
  DeleteHashTable_AS(fragmentMap);  fragmentMap = NULL;
  DeleteMANode(ma->lid);

  safe_free(offsets);

  return(FALSE);
}


