
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

static char *rcsid = "$Id: MultiAlignUnitig.c,v 1.5 2009-06-22 12:40:58 brianwalenz Exp $";

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
                 gkStore         *fragStore,
                 VA_TYPE(char)   *sequence,
                 VA_TYPE(char)   *quality,
                 VA_TYPE(int32)  *deltas,
                 CNS_PrintKey     printwhat,
                 CNS_Options     *opp) {

  int32 i;

  VA_TYPE(int32) *trace   = NULL;
  MANode         *ma      = NULL;
  SeqInterval    *offsets = NULL;
  SeqInterval    *placed  = NULL;

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
  placed  = (SeqInterval *)safe_calloc(unitig->num_frags, sizeof(SeqInterval));

  assert(ma->lid == 0);

  uint32    frankensteinLen = 0;
  uint32    frankensteinMax = 1024 * 1024;
  char     *frankenstein    = (char *)safe_malloc(sizeof(char) * frankensteinMax);
  int32    *frankensteinBof = (int32 *)safe_malloc(sizeof(int32) * frankensteinMax);

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

    // This guy allocates and initializes the beads for each fragment.  Beads are not fully inserted
    // in the abacus here.

    fid = AppendFragToLocalStore(unitig->f_list[i].type,
                                 unitig->f_list[i].ident,
                                 complement,
                                 unitig->f_list[i].contained,
                                 AS_OTHER_UNITIG, NULL);

    offsets[fid].bgn = complement ? unitig->f_list[i].position.end : unitig->f_list[i].position.bgn;
    offsets[fid].end = complement ? unitig->f_list[i].position.bgn : unitig->f_list[i].position.end;

    placed[fid].bgn  = 0;
    placed[fid].end  = 0;

    //  If this is violated, then the implicit map from offsets[] and placed[] to unitig->f_list is
    //  incorrect.
    assert(fid == i);

    if (VERBOSE_MULTIALIGN_OUTPUT)
      fprintf(stderr,"MultiAlignUnitig()-- Added fragment %d (%d-%d) in unitig %d to store at local id %d.\n",
              unitig->f_list[i].ident, unitig->f_list[i].position.bgn, unitig->f_list[i].position.end, unitig->iaccession, fid);
  }

  SeedMAWithFragment(ma->lid, GetFragment(fragmentStore,0)->lid,0, opp);

  //  Save columns
  {
    int32   bidx = GetFragment(fragmentStore, 0)->firstbead;
    Bead   *bead = GetBead(beadStore, bidx);

    while (bead) {
      frankenstein   [frankensteinLen] = *Getchar(sequenceStore, bead->soffset);
      frankensteinBof[frankensteinLen] = bead->boffset;

      frankensteinLen++;

      bead = (bead->next == -1) ? NULL : GetBead(beadStore, bead->next);
    }

    frankenstein   [frankensteinLen] = 0;
    frankensteinBof[frankensteinLen] = -1;

    placed[0].bgn = 0;
    placed[0].end = frankensteinLen;
  }



  for (i=1; i<unitig->num_frags; i++) {
    int         piid         = 0;
    int         pbeg         = 0;
    int         pend         = 0;

    Fragment   *afrag        = NULL;
    Fragment   *bfrag        = GetFragment(fragmentStore, i);  //  This frag

    int         ahang        = 0;  //  Hangs relative to frankenstein
    int         bhang        = 0;

    int         beg          = 0;  //  Positions in frankenstein
    int         end          = 0;

    int         olap_success = 0;
    OverlapType otype;

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

#define SHOW_PLACEMENT_BEFORE
#ifdef SHOW_PLACEMENT_BEFORE
      for (uint32 x=0; x<=i; x++)
        fprintf(stderr, "MultiAlignUnitig()-- %3d  f_list %6d,%6d  offsets %6d,%6d  placed %6d,%6d\n",
                x,
                unitig->f_list[x].position.bgn, unitig->f_list[x].position.end,
                offsets[x].bgn, offsets[x].end,
                placed[x].bgn, placed[x].end);
#endif
    }

    //  Place the fragment in the frankenstein using the parent and hangs.  If no parent supplied,
    //  fallback to the positions.
    //
    //   ---------------------------------
    //              |-------------|  piid == afrag
    //                 |----|        i
    //

    if (unitig->f_list[i].parent > 0) {
      for (piid = i-1; piid >= 0; piid--) {
        afrag = GetFragment(fragmentStore, piid);

        if (unitig->f_list[i].parent == afrag->iid) {
          beg = placed[piid].bgn + unitig->f_list[i].ahang;
          end = placed[piid].end + unitig->f_list[i].bhang;
          ovl = end - beg;

          ahang = beg;
          bhang = end - frankensteinLen;

          //fprintf(stderr, "PLACE(1)-- beg,end %d,%d  hangs %d,%d  fLen %d\n",
          //        beg, end, ahang, bhang, frankensteinLen);

          //  HACK.  If the positions don't agree, move along.
          int bdiff = beg - offsets[i].bgn;
          int ediff = end - offsets[i].end;

          bdiff = (bdiff < 0) ? -bdiff : bdiff;
          ediff = (ediff < 0) ? -ediff : ediff;

          if ((bdiff < 300) && (ediff < 300))
            goto doAlign;
          else
            fprintf(stderr, "PLACE(1)-- Change too big; expected %d,%d got %d,%d\n",
                    offsets[i].bgn, offsets[i].end,
                    beg, end);
        }
      }
    }


    if (unitig->f_list[i].contained != 0) {
      for (piid = i-1; piid >= 0; piid--) {
        afrag = GetFragment(fragmentStore, piid);

        if (unitig->f_list[i].contained == afrag->iid) {
          beg = placed[piid].bgn + offsets[i].bgn - offsets[afrag->lid].bgn;
          end = placed[piid].end + offsets[i].end - offsets[afrag->lid].end;
          ovl = end - beg;

          ahang = beg;
          bhang = end - frankensteinLen;

          //fprintf(stderr, "PLACE(2)-- beg,end %d,%d  hangs %d,%d  fLen %d\n",
          //        beg, end, ahang, bhang, frankensteinLen);

          goto doAlign;
        }
      }
    }


    if (unitig->f_list[i].contained == 0) {
      int   thickest    = -1;
      int   thickestLen = 0;

      //  Find the thickest piid overlap
      for (piid = i-1; piid >= 0; piid--) {
        afrag = GetFragment(fragmentStore, piid);

        if ((offsets[i].bgn < offsets[afrag->lid].end) ||
            (offsets[i].end > offsets[afrag->lid].bgn)) {
          beg = placed[piid].bgn + offsets[i].bgn - offsets[afrag->lid].bgn;
          end = placed[piid].end + offsets[i].end - offsets[afrag->lid].end;
          ovl = end - beg;

          if (thickestLen < ovl) {
            thickest    = piid;
            thickestLen = ovl;
          }
        }

        assert(thickest >= 0);

        beg = placed[thickest].bgn + offsets[i].bgn - offsets[afrag->lid].bgn;
        end = placed[thickest].end + offsets[i].end - offsets[afrag->lid].end;
        ovl = end - beg;
        
        ahang = beg;
        bhang = end - frankensteinLen;

        //fprintf(stderr, "PLACE(3)-- beg,end %d,%d  hangs %d,%d  fLen %d\n",
        //        beg, end, ahang, bhang, frankensteinLen);

        goto doAlign;
      }
    }

    assert(0);
    
  doAlign:

#define AHANGOFFSET
#ifdef AHANGOFFSET
    int32  ahang_offset = MAX(0, ahang - 100);
    int32  bhang_offset = 0;
    int32  bhang_posn   = frankensteinLen;  //  unless reset in the test below, use only for diagnostic
    char   bhang_base   = 0;

    if (bhang + 100 < 0) {
      bhang_offset = bhang + 100;
      bhang_posn   = frankensteinLen + bhang_offset;
      bhang_base   = frankenstein[bhang_posn];
      frankenstein[bhang_posn] = 0;
    }

    ahang -= ahang_offset;
    bhang -= bhang_offset;
#else
    int32 ahang_offset = 0;
#endif

    if (VERBOSE_MULTIALIGN_OUTPUT)
      fprintf(stderr, "MultiAlignUnitig()-- Aligning bfrag iid %d utgpos %d,%d to afrag iid %d utgpos %d,%d -- frankenstein=%d,%d ovl=%d hang=%d,%d\n",
              bfrag->iid,
              offsets[bfrag->lid].bgn,
              offsets[bfrag->lid].end,
              afrag->iid,
              offsets[afrag->lid].bgn,
              offsets[afrag->lid].end,
              ahang_offset, bhang_posn,
              ovl,
              ahang, bhang);


    olap_success = GetAlignmentTraceDriver(0, frankenstein + ahang_offset,
                                           bfrag,
                                           &ahang, &bhang, ovl,
                                           trace,
                                           &otype,
                                           'u',
                                           0);
    //fprintf(stderr, "GetAlignmentTrace returned success=%d hang=%d,%d\n",
    //        olap_success, ahang, bhang);

    //  Try again allowing negative hangs.
    if (!olap_success) {
      allow_neg_hang = 1;
      olap_success = GetAlignmentTraceDriver(0, frankenstein + ahang_offset,
                                             bfrag,
                                             &ahang, &bhang, ovl,
                                             trace,
                                             &otype,
                                             'u',
                                             0);
      //fprintf(stderr, "GetAlignmentTrace returned success=%d hang=%d,%d\n",
      //        olap_success, ahang, bhang);
      allow_neg_hang = 0;
    }


#ifdef AHANGOFFSET
    //  If we've set the base, we've trimmed the end of frankenstein
    //  to limit the alignment.  Restore it to what it was before.
    //
    if (bhang_base)
      frankenstein[bhang_posn] = bhang_base;

    //  Update the trace to account for the bases in A we ignored in GetAlignmentTrace()....but
    //  aren't going to ignore when we ApplyAlignment().
    //
    for (int32 *t = Getint32(trace, 0); (t != NULL) && (*t != 0); t++)
      if (*t < 0)
        *t -= ahang_offset;

    //  Make sure we didn't bump against what we trimmed off
    assert((ahang >= 0) || (ahang_offset == 0));
    assert((bhang <= 0) || (bhang_offset == 0));

    ahang += ahang_offset;
    bhang += bhang_offset;
#endif

    if (!olap_success) {
      fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Overlap not found between %sbfrag %d (%c) and afrag %d (%c)%s.\n",
              unitig->iaccession,
              (bfrag->container_iid == 0) ? "" : "contained ",
              bfrag->iid, bfrag->type,
              afrag->iid, afrag->type,
              (bfrag->container_iid == 0) ? "" : ((bfrag->container_iid == afrag->iid) ? ", the container" : ", not the container"));

      //fprintf(stderr, ">frankenstein\n%s\n", frankenstein);
      //assert(0);

      goto MAUnitigFailure;
    }


    //  Add the alignment to abacus
    //
    ApplyAlignment(-1,
                   frankensteinLen, frankensteinBof,
                   i,
                   ahang, Getint32(trace, 0));


    //  Update the placement of this read in the frankenstein
    //
    placed[i].bgn = ahang;
    placed[i].end = frankensteinLen + bhang;
    //fprintf(stderr, "PLACE(4)-- set %d to %d,%d\n", i, placed[i].bgn, placed[i].end);


    //  Update parent and hangs to reflect the overlap that succeeded.
    //
    if ((unitig->f_list[i].parent != afrag->iid) ||
        (unitig->f_list[i].ahang  != placed[i].bgn - placed[piid].bgn) ||
        (unitig->f_list[i].bhang  != placed[i].end - placed[piid].end) ||
        ((otype == AS_CONTAINMENT) && (unitig->f_list[i].contained != afrag->iid))) {
      if ((unitig->f_list[i].parent) &&
          (VERBOSE_MULTIALIGN_OUTPUT))
        fprintf(stderr, "MultiAlignUnitig()-- update parent from id=%d %d,%d to id=%d %d,%d\n",
                unitig->f_list[i].parent, unitig->f_list[i].ahang, unitig->f_list[i].bhang,
                afrag->iid,
                placed[i].bgn - placed[piid].bgn,
                placed[i].end - placed[piid].end);

      unitig->f_list[i].parent = afrag->iid;
      unitig->f_list[i].ahang  = placed[i].bgn - placed[piid].bgn;
      unitig->f_list[i].bhang  = placed[i].end - placed[piid].end;

      if (otype == AS_CONTAINMENT)
        unitig->f_list[i].contained = afrag->iid;
    }


    //
    //  Extend the frankenstein.  Son of Frankenstein!
    //
    //  We know the last bead in the current frankenstein.  We use that to find the first bead in
    //  the new sequence, then march along the new sequence copying column IDs and bases to
    //  frankenstein.
    //
    //  Details (for bhangs):
    //
    //  Grab the last bead of the current frankenstein.  That bead should be the first thing added
    //  to a column, and so should be on the bottom of the column (that's the assert -- if not, grab
    //  the column, and find the last bead).  Then, search up the column for the first bead from the
    //  current fragment.
    //
    //  This bead is the last thing in the current frankenstein.  Move one spot to the right, now
    //  we're at the first thing we need to add to frankenstein.  Walk along the beads, adding to
    //  frankenstein, until there are no more beads.
    //
    //  Details (for ahangs):
    //
    //  Very similar to the bhang case but complicated in that we push bases onto the start of
    //  frankenstein, and so we must also update the fragment position mapping array.

    if (bhang > 0) {
      int32   bidx = frankensteinBof[frankensteinLen-1];
      Bead   *bead = GetBead(beadStore, bidx);

      assert(bead->down == -1);  //  Just searching up will break if this triggers

      while ((bead) && (bead->frag_index != i))
        bead = (bead->up == -1) ? NULL : GetBead(beadStore, bead->up);

      assert((bead) && (bead->frag_index == i));  //  Never found the correct fragment?!?

      for (bead = (bead->next == -1) ? NULL : GetBead(beadStore, bead->next);
           bead;
           bead = (bead->next == -1) ? NULL : GetBead(beadStore, bead->next)) {
        frankenstein   [frankensteinLen] = *Getchar(sequenceStore, bead->soffset);
        frankensteinBof[frankensteinLen] = bead->boffset;

#undef SHOW_FRANKENSTEIN_ADD
#ifdef SHOW_FRANKENSTEIN_ADD
        fprintf(stderr, "frankenstein[%4d]-- Added bead %d,%c frag=%d\n",
                frankensteinLen,
                frankensteinBof[frankensteinLen],
                frankenstein   [frankensteinLen],
                bead->frag_index);
#endif

        frankensteinLen++;

        if (frankensteinLen > frankensteinMax) {
          //  Just being lazy; need to reallocate this.
          assert(frankensteinLen < frankensteinMax);
        }
      }

      frankenstein   [frankensteinLen] = 0;
      frankensteinBof[frankensteinLen] = -1;
    }  //  End of extending to the right.


    if (ahang < 0) {
      if (frankensteinLen + -ahang > frankensteinMax) {
        //  Just being lazy; need to reallocate this.
        assert(frankensteinLen + -ahang < frankensteinMax);
      }

      for (int32 x=frankensteinLen; x>=0; x--) {
        frankenstein   [x + -ahang] = frankenstein   [x];
        frankensteinBof[x + -ahang] = frankensteinBof[x];
      }
      frankensteinLen += -ahang;

      for (int32 x=0; x<-ahang; x++) {
        frankenstein   [x] = 0;
        frankensteinBof[x] = -1;
      }

      for (int32 x=0; x<=i; x++) {
        placed[x].bgn  += -ahang;
        placed[x].end  += -ahang;
      }

      int32   bidx = frankensteinBof[-ahang];
      Bead   *bead = GetBead(beadStore, bidx);

      assert(bead->down == -1);  //  Just searching up will break if this triggers
      assert(bead->prev == -1);  //  Should be the first bead in the frankenstein

      while ((bead) && (bead->frag_index != i))
        bead = (bead->up == -1) ? NULL : GetBead(beadStore, bead->up);

      assert((bead) && (bead->frag_index == i));  //  Never found the correct fragment?!?

      while (bead->prev != -1) {
        //fprintf(stderr, "prev bead: boffset %d prev %d\n", bead->boffset, bead->prev);
        bead = GetBead(beadStore, bead->prev);
      }

      assert((bead) && (bead->frag_index == i));  //  Never found the correct fragment?!?

      for (int x=0; x<-ahang; x++) {
        frankenstein   [x] = *Getchar(sequenceStore, bead->soffset);
        frankensteinBof[x] = bead->boffset;

        bead = GetBead(beadStore, bead->next);
      }
    }  //  End of extending to the left.
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
  safe_free(placed);
  safe_free(frankenstein);
  safe_free(frankensteinBof);

  return(TRUE);

 MAUnitigFailure:
  DeleteVA_int32(trace);
  DeleteHashTable_AS(fragmentMap);  fragmentMap = NULL;
  DeleteMANode(ma->lid);

  safe_free(offsets);
  safe_free(placed);
  safe_free(frankenstein);
  safe_free(frankensteinBof);

  return(FALSE);
}


