
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

static char *rcsid = "$Id: MultiAlignUnitig.c,v 1.11 2009-07-07 17:29:04 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"

#undef SHOW_PLACEMENT_BEFORE
#undef SHOW_PLACEMENT

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


class unitigConsensus {
public:
  unitigConsensus(IntUnitigMesg *unitig_, CNS_Options *opp_) {
    unitig = unitig_;
    opp    = opp_;
    trace  = NULL;
    ma     = NULL;
    offsets = NULL;
    placed  = NULL;
    ovl     = 0;
    ahang   = 0;
    bhang   = 0;
    tiid    = 0;
    piid    = 0;
    frankensteinLen = 0;
    frankensteinMax = 0;
    frankenstein    = NULL;
    frankensteinBof = NULL;
  };

  ~unitigConsensus() {
    DeleteVA_int32(trace);
    DeleteHashTable_AS(fragmentMap);  fragmentMap = NULL;
    DeleteMANode(ma->lid);

    safe_free(offsets);
    safe_free(placed);
    safe_free(frankenstein);
    safe_free(frankensteinBof);
  };

  int  initialize(void); 

  void reportStartingWork(void);
  void reportFailure(void);

  int  moreFragments(void)  { tiid++;  return (tiid < unitig->num_frags); };

  int  computePositionFromParent(void);
  int  computePositionFromContainer(void);
  int  computePositionFromLayout(void);
  int  computePositionFromAlignment(void);
  int  computePositionFromConsensus(void);

  int  alignFragment(void);
  void applyAlignment(void);

  void rebuildFrankensteinFromFragment(void);

  void generateConsensus(VA_TYPE(char)   *sequence,
                         VA_TYPE(char)   *quality,
                         VA_TYPE(int32)  *deltas,
                         CNS_PrintKey     printwhat);

private:
  IntUnitigMesg  *unitig;
  CNS_Options    *opp;

  VA_TYPE(int32) *trace;
  MANode         *ma;
  SeqInterval    *offsets;
  SeqInterval    *placed;

  int32           ovl;    //  Expected overlap in bases to the frankenstein
  int32           ahang;  //  Expected hangs to the frankenstein
  int32           bhang;

  int32           tiid;   //  This frag IID

  int32           piid;   //  Parent frag IID

  uint32          frankensteinLen;
  uint32          frankensteinMax;
  char           *frankenstein;
  int32          *frankensteinBof;
};


void
unitigConsensus::reportStartingWork(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "MultiAlignUnitig()-- processing fragment %d iid=%d pos=%d,%d parent=%d,%d,%d contained=%d\n",
          tiid,
          unitig->f_list[tiid].ident,
          unitig->f_list[tiid].position.bgn,
          unitig->f_list[tiid].position.end,
          unitig->f_list[tiid].parent,
          unitig->f_list[tiid].ahang,
          unitig->f_list[tiid].bhang,
          unitig->f_list[tiid].contained);

#ifdef SHOW_PLACEMENT_BEFORE
  for (uint32 x=0; x<=tiid; x++)
    fprintf(stderr, "MultiAlignUnitig()-- %3d  f_list %6d,%6d  offsets %6d,%6d  placed %6d,%6d\n",
            x,
            unitig->f_list[x].position.bgn, unitig->f_list[x].position.end,
            offsets[x].bgn, offsets[x].end,
            placed[x].bgn, placed[x].end);
#endif
}


void
unitigConsensus::reportFailure(void) {
  fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Could not align fragment %d.\n",
          unitig->iaccession, unitig->f_list[tiid].ident);
  //fprintf(stderr, ">frankenstein\n%s\n", frankenstein);
}

int
unitigConsensus::initialize(void) {

  int32 num_columns = 0;

  for (int32 i=0; i<unitig->num_frags; i++) {
    num_columns = (unitig->f_list[i].position.bgn > num_columns) ? unitig->f_list[i].position.bgn : num_columns;
    num_columns = (unitig->f_list[i].position.end > num_columns) ? unitig->f_list[i].position.end : num_columns;
  }

  ResetStores(unitig->num_frags, num_columns);

  //  Magic initialization (in ResetStores()) prevents us calling CreateMANode() until now.

  trace   = CreateVA_int32(2 * AS_READ_MAX_LEN);
  ma      = CreateMANode(unitig->iaccession);
  offsets = (SeqInterval *)safe_calloc(unitig->num_frags, sizeof(SeqInterval));
  placed  = (SeqInterval *)safe_calloc(unitig->num_frags, sizeof(SeqInterval));

  assert(ma->lid == 0);

  frankensteinLen = 0;
  frankensteinMax = 1024 * 1024;
  frankenstein    = (char *)safe_malloc(sizeof(char) * frankensteinMax);
  frankensteinBof = (int32 *)safe_malloc(sizeof(int32) * frankensteinMax);

  for (int32 i=0; i<unitig->num_frags; i++) {
    int32 complement = (unitig->f_list[i].position.bgn < unitig->f_list[i].position.end) ? 0 : 1;
    int32 fid;

    if (unitig->f_list[i].type != AS_READ) {
      fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Fragment %d is not a read.\n",
              unitig->iaccession, unitig->f_list[i].ident);
      return(false);
    }

    if (HASH_SUCCESS != InsertInHashTable_AS(fragmentMap,unitig->f_list[i].ident, 0, 1, 0)) {
      fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Fragment %d is a duplicate.\n",
              unitig->iaccession, unitig->f_list[i].ident);
      return(false);
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

  return(true);
}



//  Place the fragment in the frankenstein using the parent and hangs.  If no parent supplied,
//  fallback to the positions.
//
//   ---------------------------------
//              |-------------|  piid == afrag
//                 |----|        i
//
int
unitigConsensus::computePositionFromParent(void) {

  if (unitig->f_list[tiid].parent == 0)
    return(false);

  for (piid = tiid-1; piid >= 0; piid--) {
    Fragment *afrag = GetFragment(fragmentStore, piid);

    if (unitig->f_list[tiid].parent == afrag->iid) {
      int32 beg = placed[piid].bgn + unitig->f_list[tiid].ahang;
      int32 end = placed[piid].end + unitig->f_list[tiid].bhang;

      ovl   = end - beg;
      ahang = beg;
      bhang = end - frankensteinLen;

#ifdef SHOW_PLACEMENT
      fprintf(stderr, "PLACE(1)-- beg,end %d,%d  hangs %d,%d  fLen %d\n",
              beg, end, ahang, bhang, frankensteinLen);
#endif

      //  HACK.  If the positions don't agree, move along.  BOG sometimes supplies the wrong
      //  parent for a read.
      int32 bdiff = beg - offsets[tiid].bgn;
      int32 ediff = end - offsets[tiid].end;

      bdiff = (bdiff < 0) ? -bdiff : bdiff;
      ediff = (ediff < 0) ? -ediff : ediff;

      if ((bdiff < 300) && (ediff < 300)) {
        return(true);
      } else {
        fprintf(stderr, "PLACE(1)-- Change too big; expected %d,%d got %d,%d\n",
                offsets[tiid].bgn, offsets[tiid].end,
                beg, end);
        return(false);
      }
    }
  }
  return(false);
}



int
unitigConsensus::computePositionFromContainer(void) {

  if (unitig->f_list[tiid].contained == 0)
    return(false);

  for (piid = tiid-1; piid >= 0; piid--) {
    Fragment *afrag = GetFragment(fragmentStore, piid);

    if (unitig->f_list[tiid].contained == afrag->iid) {
      int32 beg = placed[piid].bgn + offsets[tiid].bgn - offsets[afrag->lid].bgn;
      int32 end = placed[piid].end + offsets[tiid].end - offsets[afrag->lid].end;

      ovl   = end - beg;
      ahang = beg;
      bhang = end - frankensteinLen;

#ifdef SHOW_PLACEMENT
      fprintf(stderr, "PLACE(2)-- beg,end %d,%d  hangs %d,%d  fLen %d\n",
              beg, end, ahang, bhang, frankensteinLen);
#endif

      return(true);
    }
  }

  return(false);
}


int
unitigConsensus::computePositionFromLayout(void) {
  int32   thickestLen = 0;

  //  Find the thickest qiid overlap
  for (uint32 qiid = tiid-1; qiid >= 0; qiid--) {
    if ((offsets[tiid].bgn < offsets[qiid].end) ||
        (offsets[tiid].end > offsets[qiid].bgn)) {
      int32 beg = placed[qiid].bgn + offsets[tiid].bgn - offsets[qiid].bgn;
      int32 end = placed[qiid].end + offsets[tiid].end - offsets[qiid].end;

      int32 ooo = MIN(end, frankensteinLen) - beg;

      if (thickestLen < ooo) {
        thickestLen = ooo;

        ovl   = MIN(end, frankensteinLen) - beg;
        ahang = beg;
        bhang = end - frankensteinLen;

        piid  = qiid;
      }
    }

    //  Oh crap, no overlap in unitig coords?  Something broke somewhere else.
    assert(thickestLen > 0);

#ifdef SHOW_PLACEMENT
    fprintf(stderr, "PLACE(3)-- beg,end %d,%d  hangs %d,%d  fLen %d\n",
            ahang, bhang + frankensteinLen, ahang, bhang, frankensteinLen);
#endif

    return(true);
  }

  return(false);
}



int
unitigConsensus::computePositionFromAlignment(void) {

  //  Occasionally we get a fragment that just refuses to go in the correct spot.  Search for the
  //  correct placement in all of frankenstein, update ahang,bhang and retry.

  Overlap  *O           = NULL;
  double    thresh      = 1e-6;
  int32     minlen      = AS_OVERLAP_MIN_LEN;

  char     *fragment    = Getchar(sequenceStore, GetFragment(fragmentStore, tiid)->sequence);
  int32     fragmentLen = strlen(fragment);

  if (VERBOSE_MULTIALIGN_OUTPUT)
    fprintf(stderr, "MultiAlignUnitig()--  Attempting to find placement from frankenstein using DP_Compare.\n");
  O = DP_Compare(frankenstein,
                 fragment,
                 -fragmentLen, frankensteinLen,  //  ahang bounds
                 frankensteinLen, fragmentLen,   //  length of fragments
                 0,
                 AS_CNS_ERROR_RATE, thresh, minlen,
                 AS_FIND_ALIGN);

  if (O == NULL) {
    if (VERBOSE_MULTIALIGN_OUTPUT)
      fprintf(stderr, "MultiAlignUnitig()--  Attempting to find placement from frankenstein using Local_Overlap.\n");
    O = Local_Overlap_AS_forCNS(frankenstein,
                                fragment,
                                -fragmentLen, frankensteinLen,  //  ahang bounds
                                frankensteinLen, fragmentLen,   //  length of fragments
                                0,
                                AS_CNS_ERROR_RATE, thresh, minlen,
                                AS_FIND_ALIGN);
  }

  if (O == NULL)
    return(false);

  //  From the overlap and existing placements, find the thickest overlap, to set the piid and
  //  hangs, then reset the original placement based on that parents original placement.

  placed[tiid].bgn = O->begpos;
  placed[tiid].end = O->endpos + frankensteinLen;

  int32   thickestLen = 0;

  for (uint32 qiid = tiid-1; qiid >= 0; qiid--) {
    if ((placed[tiid].bgn < placed[qiid].end) ||
        (placed[tiid].end > placed[qiid].bgn)) {
      int32 beg = placed[tiid].bgn;
      int32 end = placed[tiid].end;

      int32 ooo = MIN(end, frankensteinLen) - beg;

      if (thickestLen < ooo) {
        thickestLen = ooo;

        ovl   = MIN(end, frankensteinLen) - beg;
        ahang = beg;
        bhang = end - frankensteinLen;

        piid  = qiid;
      }
    }

    //  Oh crap, no overlap in unitig coords?  Something broke somewhere else.
    assert(thickestLen > 0);

    offsets[tiid].bgn = offsets[piid].bgn + placed[tiid].bgn - placed[piid].bgn;
    offsets[tiid].end = offsets[piid].end + placed[tiid].end - placed[piid].end;

#ifdef SHOW_PLACEMENT
    fprintf(stderr, "PLACE(5)-- beg,end %d,%d  hangs %d,%d  fLen %d\n",
            ahang, bhang + frankensteinLen, ahang, bhang, frankensteinLen);
#endif

    return(true);
  }

  return(false);
}


int
unitigConsensus::computePositionFromConsensus(void) {

  //  Run abacus to rebuild an intermediate consensus sequence.  VERY expensive, and doesn't
  //  update the placed[] array...but the changes shouldn't be huge.

  RefreshMANode(ma->lid, 0, opp, NULL, NULL, 0, 0);

  //  Are all three needed??

  AbacusRefine(ma,0,-1,CNS_SMOOTH, opp);
  MergeRefine(ma->lid, NULL, NULL, 1, opp, 1);

  AbacusRefine(ma,0,-1,CNS_POLYX, opp);
  MergeRefine(ma->lid, NULL, NULL, 1, opp, 1);

  AbacusRefine(ma,0,-1,CNS_INDEL, opp);
  MergeRefine(ma->lid, NULL, NULL, 1, opp, 1);

  //  Extract the consensus sequence

  char          *newFrank    = (char *)safe_calloc(frankensteinMax, sizeof(char));
  int32         *newFrankBof = (int32 *)safe_calloc(frankensteinMax, sizeof(int32));
  int32          newFrankLen = 0;
  SeqInterval   *newPlaced   = (SeqInterval *)safe_calloc(unitig->num_frags, sizeof(SeqInterval));

  for (int32 i=0; i<unitig->num_frags; i++) {
    newPlaced[i].bgn = frankensteinMax;
    newPlaced[i].end = 0;
  }

  char          *oldFrank    = frankenstein;
  int32         *oldFrankBof = frankensteinBof;
  int32          oldFrankLen = frankensteinLen;
  SeqInterval   *oldPlaced   = placed;

  ConsensusBeadIterator  bi;
  int32                  bid;

  CreateConsensusBeadIterator(ma->lid, &bi);

  while ((bid = NextConsensusBead(&bi)) != -1) {
    Bead *bead = GetBead(beadStore, bid);
    char  cnsc = *Getchar(sequenceStore, bead->soffset);

    if (cnsc == '-')
      continue;

    while (bead->down != -1) {
      if (bead->frag_index >= 0) {
        if (newPlaced[bead->frag_index].bgn > newFrankLen)
          newPlaced[bead->frag_index].bgn = newFrankLen;
        if (newPlaced[bead->frag_index].end < newFrankLen)
          newPlaced[bead->frag_index].end = newFrankLen;
      }

      bead = GetBead(beadStore, bead->down);
    }

    newFrank   [newFrankLen] = cnsc;
    newFrankBof[newFrankLen] = bead->boffset;
    newFrankLen++;
  }

  newFrank   [newFrankLen] = 0;
  newFrankBof[newFrankLen] = -1;

  //  Call the usual position-from-frankenstein function, but swap in the new one.  This might be
  //  slightly screwing up our placement.  Hopefully that error won't propagate.

  frankenstein    = newFrank;
  frankensteinBof = newFrankBof;
  frankensteinLen = newFrankLen;
  placed          = newPlaced;

  int    success = computePositionFromAlignment();

  if (success) {
    safe_free(oldFrank);
    safe_free(oldFrankBof);
    safe_free(oldPlaced);
  } else {
    frankenstein    = oldFrank;
    frankensteinBof = oldFrankBof;
    frankensteinLen = oldFrankLen;
    placed          = oldPlaced;

    safe_free(newFrank);
    safe_free(newFrankBof);
    safe_free(newPlaced);
  }

  return(success);
}




int
unitigConsensus::alignFragment(void) {

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
    fprintf(stderr, "MultiAlignUnitig()-- Aligning fragment %d utgpos %d,%d to frankenstein=%d,%d ovl=%d hang=%d,%d\n",
            tiid,
            offsets[tiid].bgn,
            offsets[tiid].end,
            ahang_offset, bhang_posn,
            ovl,
            ahang, bhang);

  Fragment   *bfrag   = GetFragment(fragmentStore, tiid);

  int32       success = false;
  OverlapType otype;

  if (success == false) {
    success = GetAlignmentTraceDriver(NULL, frankenstein + ahang_offset, bfrag, &ahang, &bhang, ovl, trace, &otype, 'u', 0);
  }

  //  Try again allowing negative hangs.
  //
  if (success == false) {
    allow_neg_hang = 1;

    success = GetAlignmentTraceDriver(NULL, frankenstein + ahang_offset, bfrag, &ahang, &bhang, ovl, trace, &otype, 'u', 0);

    allow_neg_hang = 0;
  }

  //  Occasionally we get a short contained fragment that has just enough errors to frankenstein,
  //  and is just short enough, to have a high error alignment.  Allow a high error alignment.
  //
  if (success == false) {
    double  origErate = AS_CNS_ERROR_RATE;

    AS_CNS_ERROR_RATE = MIN(AS_MAX_ERROR_RATE, 2 * AS_CNS_ERROR_RATE);

    success = GetAlignmentTraceDriver(NULL, frankenstein + ahang_offset, bfrag, &ahang, &bhang, ovl, trace, &otype, 'u', 0);

    AS_CNS_ERROR_RATE = origErate;
  }


#ifdef AHANGOFFSET
  //  If we've set the base, we've trimmed the end of frankenstein to limit the alignment.  Restore
  //  it to what it was before.
  //
  if (bhang_base)
    frankenstein[bhang_posn] = bhang_base;

  //  Update the trace to account for the bases in A we ignored in GetAlignmentTrace()....but aren't
  //  going to ignore when we ApplyAlignment().
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

  return(success);
}




void
unitigConsensus::applyAlignment(void) {

  //  Add the alignment to abacus
  //
  ApplyAlignment(-1,
                 frankensteinLen, frankensteinBof,
                 tiid,
                 ahang, Getint32(trace, 0));


  //  Update the placement of this read in the frankenstein
  //
  placed[tiid].bgn = ahang;
  placed[tiid].end = frankensteinLen + bhang;

#ifdef SHOW_PLACEMENT
  fprintf(stderr, "PLACE(4)-- set %d to %d,%d\n", tiid, placed[tiid].bgn, placed[tiid].end);
#endif

  //  Update parent and hangs to reflect the overlap that succeeded.
  //
  //  Containment is obvious for the bhang; if the ahang is negative, we
  //  are not contained.
  //
  //  We should probably reest everything for negative ahangs...
  //
  unitig->f_list[tiid].parent    = unitig->f_list[piid].ident;
  unitig->f_list[tiid].ahang     = placed[tiid].bgn - placed[piid].bgn;
  unitig->f_list[tiid].bhang     = placed[tiid].end - placed[piid].end;
  unitig->f_list[tiid].contained = (bhang > 0) ? 0 : unitig->f_list[piid].ident;
  unitig->f_list[tiid].contained = (ahang < 0) ? 0 : unitig->f_list[tiid].contained;


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

    while ((bead) && (bead->frag_index != tiid))
      bead = (bead->up == -1) ? NULL : GetBead(beadStore, bead->up);

    assert((bead) && (bead->frag_index == tiid));  //  Never found the correct fragment?!?

    for (bead = (bead->next == -1) ? NULL : GetBead(beadStore, bead->next);
         bead;
         bead = (bead->next == -1) ? NULL : GetBead(beadStore, bead->next)) {
      char ch = *Getchar(sequenceStore, bead->soffset);
      if (ch != '-') {
        frankenstein   [frankensteinLen] = ch;
        frankensteinBof[frankensteinLen] = bead->boffset;
        frankensteinLen++;
      } else {
        //  Here, ch should never be a gap, since we're just tacking on completely new
        //  sequence.
      }

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

    //  Make space for the new stuff
    for (int32 x=frankensteinLen; x>=0; x--) {
      frankenstein   [x + -ahang] = frankenstein   [x];
      frankensteinBof[x + -ahang] = frankensteinBof[x];
    }
    frankensteinLen += -ahang;

    //  Zero out the new stuff, temporarily.
    for (int32 x=0; x<-ahang; x++) {
      frankenstein   [x] = 0;
      frankensteinBof[x] = -1;
    }

    //  Adjust positions.
    for (int32 x=0; x<=tiid; x++) {
      placed[x].bgn  += -ahang;
      placed[x].end  += -ahang;
    }

    int32   bidx = frankensteinBof[-ahang];
    Bead   *bead = GetBead(beadStore, bidx);

    assert(bead->down == -1);  //  Just searching up will break if this triggers
    assert(bead->prev == -1);  //  Should be the first bead in the frankenstein

    while ((bead) && (bead->frag_index != tiid))
      bead = (bead->up == -1) ? NULL : GetBead(beadStore, bead->up);

    assert((bead) && (bead->frag_index == tiid));  //  Never found the correct fragment?!?

    while (bead->prev != -1) {
      //fprintf(stderr, "prev bead: boffset %d prev %d\n", bead->boffset, bead->prev);
      bead = GetBead(beadStore, bead->prev);
    }

    assert((bead) && (bead->frag_index == tiid));  //  Never found the correct fragment?!?

    //  Append the new stuff, and VERY IMPORTANT, steal the former first position from the original
    //  (that's the <= test).  If this isn't done, we eventually hit the assert at the end of
    //  findBeadInColumn, since frankenstein doesn't have a "dovetail overlap":
    //
    //  FRANK   BAAAAAA...
    //  A        ------...
    //  B       -------...
    //
    //  We need to move from the first base in B, up the abacus to fragment A, but can't do that
    //  since fragment A has nothing in abacus.  To prevent this, we switch the second column from A
    //  to B.
    //
    for (int32 x=0; x<=-ahang; x++) {
      frankenstein   [x] = *Getchar(sequenceStore, bead->soffset);
      frankensteinBof[x] = bead->boffset;

      bead = GetBead(beadStore, bead->next);
    }
  }  //  End of extending to the left.
}


//  Attempt to reduce the amount of fragmentation and noise near the 5' end of frankenstein by
//  rebuilding using all of the last fragment if that last was not contained.
//
//  Rebuild:            No rebuild:         No rebuild:
//
//  ----------          ----------------      --------
//    -----------          ---------        -------
//
//  This seems to break ApplyAlignment.  Perhaps we violate some unknown constraint
//  about gap placement?  Additionally, it doesn't take into account the orientation
//  of fragments, and without doing that, this seems to result in a net loss in
//  the quality of frankenstein.
//
//  If -----> is a fragment from 5' to 3', then the '>' end is of lower quality.
//  Assuming thick overlaps:
//
//                  append                   replace
//  ------->        adds low quality         improves the overlapping sequence
//    -------->     bases
//
//  ------->        adds high quality        changes cancel out
//    <--------     bases
//
//  <-------        adds low q               changes cancel out
//    -------->
//
//  <-------        adds high q              replaces high q overlap with low q
//    <--------
//
//
void
unitigConsensus::rebuildFrankensteinFromFragment(void) {

  if ((ahang <= 0) ||
      (bhang <= 0))
    return;

  Bead  *bead = GetBead(beadStore, frankensteinBof[frankensteinLen-1]);

  //  If the last bead is in this fragment, then we used at least one base of read i in the
  //  frankenstein.  Back up until we find the first column used in this alignment, search up
  //  for the Bof bead (which might be a gap) and rebuild.
  //
  //  Once we find the start of this read, we chop off all columns equal to or larger than
  //  the column we are in, and reform using this read.

  assert(bead->frag_index == tiid);

  //  Potentially, we can skip all this searching if we assume frankenstein[ahang] is a bead in
  //  the first column with this read.

  //  Move to the bottom of the column....
  while ((bead) && (bead->down != -1))
    bead = GetBead(beadStore, bead->down);

  //  Then search up for the fragment we just aligned.
  while ((bead) && (bead->frag_index != tiid))
    bead = (bead->up == -1) ? NULL : GetBead(beadStore, bead->up);

  assert(bead);

  //  Move to the start of the fragment.
  while ((bead) && (bead->prev != -1))
    bead = GetBead(beadStore, bead->prev);

  assert(bead);
  assert(bead->prev == -1);
  assert(bead->frag_index == tiid);

  //  Trim back frankenstein until just before the start of the fragment.
  int32 fidx = frankensteinLen - 1;
  while ((fidx >= 0) && (GetBead(beadStore, frankensteinBof[fidx])->column_index >= bead->column_index))
    fidx--;

  frankensteinLen = fidx + 1;

  //  And append bases for this fragment.
  for (;
       bead;
       bead = (bead->next == -1) ? NULL : GetBead(beadStore, bead->next)) {
    char ch = *Getchar(sequenceStore, bead->soffset);

    if (ch != '-') {
      assert(bead->frag_index == tiid);

      frankenstein   [frankensteinLen] = ch;
      frankensteinBof[frankensteinLen] = bead->boffset;
      frankensteinLen++;
    } else {
      //  Here, ch CAN be a gap, since we're adding in sequence from a multialignment.
    }

    if (frankensteinLen > frankensteinMax) {
      //  Just being lazy; need to reallocate this.
      assert(frankensteinLen < frankensteinMax);
    }
  }

  frankenstein   [frankensteinLen] = 0;
  frankensteinBof[frankensteinLen] = -1;
}




void
unitigConsensus::generateConsensus(VA_TYPE(char)   *sequence,
                 VA_TYPE(char)   *quality,
                 VA_TYPE(int32)  *deltas,
                 CNS_PrintKey     printwhat) {

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


  //  While we have fragments in memory, compute the microhet probability.  Ideally, this would be
  //  done in CGW when loading unitigs (the only place the probability is used) but the code wants
  //  to load sequence and quality for every fragment, and that's too expensive.
  {
    int    depth  = 0;
    char **multia = NULL;
    int  **id_array = NULL;

    MANode2Array(ma, &depth, &multia, &id_array,0);

    unitig->microhet_prob = AS_REZ_MP_MicroHet_prob(multia, id_array, gkpStore, unitig->length, depth);

    for (int32 i=0;i<depth;i++) {
      safe_free(multia[2*i]);
      safe_free(multia[2*i+1]);
      safe_free(id_array[i]);
    }
    safe_free(multia);
    safe_free(id_array);
  }
}


int
MultiAlignUnitig(IntUnitigMesg   *unitig,
                 gkStore         *fragStore,
                 VA_TYPE(char)   *sequence,
                 VA_TYPE(char)   *quality,
                 VA_TYPE(int32)  *deltas,
                 CNS_PrintKey     printwhat,
                 CNS_Options     *opp) {

  gkpStore    = fragStore;
  fragmentMap = CreateScalarHashTable_AS();

  unitigConsensus uc(unitig, opp);

  if (uc.initialize() == FALSE)
    return(FALSE);

  while (uc.moreFragments()) {
    if (VERBOSE_MULTIALIGN_OUTPUT)
      uc.reportStartingWork();

    int success = false;

    if ((success == false) && (uc.computePositionFromParent()))
      success = uc.alignFragment();

    if ((success == false) && (uc.computePositionFromContainer()))
      success = uc.alignFragment();

    if ((success == false) && (uc.computePositionFromLayout()))
      success = uc.alignFragment();

    if ((success == false) && (uc.computePositionFromAlignment()))
      success = uc.alignFragment();

    if ((success == false) && (uc.computePositionFromConsensus()))
      success = uc.alignFragment();

    if (success == false) {
      uc.reportFailure();
      return(FALSE);
    }

    uc.applyAlignment();
  }

  uc.generateConsensus(sequence, quality, deltas, printwhat);

  return(TRUE);
}
