
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

static char *rcsid = "$Id: MultiAlignUnitig.c,v 1.44 2011-11-23 07:25:38 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"

#define SHOW_ALGORITHM         2
#define SHOW_PLACEMENT_BEFORE  3
#define SHOW_PLACEMENT         3


static
int
MANode2Array(MANode *ma, int32 *depth, char ***array, int32 ***id_array,
             int32 show_cel_status) {
  char **multia;
  int32 **ia;
  int32 length = GetNumColumns(ma->columnList);
  // find max column depth.
  int32 max_depth=0;
  int32 col_depth;
  int32 column_index;
  Column *col;
  char laneformat[40];
  int32 num_frags=GetNumFragments(fragmentStore);
  Fragment *frag;
  int32 fid;
  int32 *rowptr,*row_assign;
  int32 ir,fbgn,fend;
  int32 i;
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
  rowptr = (int32 *)safe_malloc((*depth)*sizeof(int));
  row_assign = (int32 *)safe_malloc(num_frags*sizeof(int));
  for (ir=0;ir<*depth;ir++) rowptr[ir] = 0;
  for (ir=0;ir<num_frags;ir++) row_assign[ir] = -1;
  frag = GetFragment(fragmentStore,0);
  // setup the packing
  for ( fid=0;fid<num_frags;fid++ ) {
    if ( frag->type != AS_UNITIG ) {
      fbgn = GetColumn(columnStore, (GetBead(beadStore,frag->firstbead                     ))->column_index)->ma_index;
      fend = GetColumn(columnStore, (GetBead(beadStore,frag->firstbead.get()+frag->length-1))->column_index)->ma_index+1;
      for (ir=0;ir<*depth;ir++) {
        if (fbgn <  rowptr[ir] ) continue;
        rowptr[ir] = fend;
        row_assign[fid] = ir;
        break;
      }
      if (row_assign[fid] <= -1)
        {
          *depth += max_depth;
          rowptr = (int32 *)safe_realloc(rowptr, (*depth)*sizeof(int));
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
  ia = (int32 **)safe_malloc((*depth)*sizeof(int32 *));
  sprintf(laneformat,"%%%ds",length);
  {int32 j;
    for (i=0;i<(*depth);i++) {
      ia[i] = (int32 *) safe_malloc( length*sizeof(int));
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
    beadIdx bid;
    char bc,bq;
    Column *bcolumn;
    int32 ma_index;

    frag = GetFragment(fragmentStore,0);
    for ( fid=0;fid<num_frags;fid++ ) {
      if ( frag->type != AS_UNITIG ) {
        ir = row_assign[fid];
        fb = GetBead(beadStore,frag->firstbead);
        bcolumn =  GetColumn(columnStore,fb->column_index);

        CreateFragmentBeadIterator(fid,&fi);

        while ( (bid = NextFragmentBead(&fi)) .isValid() ) {
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
  unitigConsensus(MultiAlignT *ma_, CNS_Options *opp_) {
    ma       = ma_;
    fraglist = GetVA_IntMultiPos(ma->f_list, 0);
    numfrags = GetNumIntMultiPoss(ma->f_list);
    opp      = opp_;
    trace    = NULL;
    manode   = NULL;
    utgpos   = NULL;
    cnspos   = NULL;
    tiid     = 0;
    piid     = -1;

    frankensteinLen = 0;
    frankensteinMax = 0;
    frankenstein    = NULL;
    frankensteinBof = NULL;

    fragmentMap = CreateScalarHashTable_AS();
  };

  ~unitigConsensus() {
    DeleteVA_int32(trace);
    DeleteHashTable_AS(fragmentMap);  fragmentMap = NULL;
    if (manode)
      DeleteMANode(manode->lid);

    safe_free(utgpos);
    safe_free(cnspos);
    safe_free(frankenstein);
    safe_free(frankensteinBof);
  };

  int32  initialize(void); 

  void   reportStartingWork(void);
  void   reportFailure(void);

  int32  moreFragments(void)  { tiid++;  return (tiid < numfrags); };

  int32  computePositionFromParent(bool doContained);
  int32  computePositionFromLayout(void);
  int32  computePositionFromAlignment(void);

  void   rebuild(bool recomputeFullConsensus);

  bool   rejectAlignment(bool allowBhang, ALNoverlap *O);
  int32  alignFragmentToFragments(void);

  int32  alignFragment(void);
  void   applyAlignment(int32 frag_aiid=-1, int32 frag_ahang=0, int32 *frag_trace=NULL);

  void   generateConsensus(CNS_PrintKey     printwhat);

  bool   fixFailures(void);

private:
  MultiAlignT    *ma;
  IntMultiPos    *fraglist;
  int32           numfrags;

  CNS_Options    *opp;

  VA_TYPE(int32) *trace;
  int32           traceBgn;

  MANode         *manode;
  SeqInterval    *utgpos;   //  Original unitigger location, DO NOT MODIFY
  SeqInterval    *cnspos;   //  Actual location in frankenstein.

  //int32           ovl;    //  Expected overlap in bases to the frankenstein
  //int32           ahang;  //  Expected hangs to the frankenstein
  //int32           bhang;

  int32           tiid;   //  This frag IID
  int32           piid;   //  Parent frag IID - if -1, not valid

  int32           frankensteinLen;
  int32           frankensteinMax;
  char           *frankenstein;
  beadIdx         *frankensteinBof;
};


void
unitigConsensus::reportStartingWork(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "MultiAlignUnitig()-- processing fragment mid %d pos %d,%d parent %d,%d,%d contained %d\n",
          fraglist[tiid].ident,
          fraglist[tiid].position.bgn,
          fraglist[tiid].position.end,
          fraglist[tiid].parent,
          fraglist[tiid].ahang,
          fraglist[tiid].bhang,
          fraglist[tiid].contained);

  if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_PLACEMENT_BEFORE)
    for (int32 x=0; x<=tiid; x++)
      fprintf(stderr, "MultiAlignUnitig()-- mid %3d  f_list %6d,%6d  utgpos %6d,%6d  cnspos %6d,%6d\n",
              fraglist[x].ident,
              fraglist[x].position.bgn, fraglist[x].position.end,
              utgpos[x].bgn, utgpos[x].end,
              cnspos[x].bgn, cnspos[x].end);
}


void
unitigConsensus::reportFailure(void) {
  fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Could not align fragment %d.\n",
          ma->maID, fraglist[tiid].ident);
  //fprintf(stderr, ">frankenstein\n%s\n", frankenstein);
}

int
unitigConsensus::initialize(void) {

  int32 num_columns = 0;
  int32 num_bases   = 0;

  if (numfrags == 0)
    return(false);

  for (int32 i=0; i<numfrags; i++) {
    int32 flen   = (fraglist[i].position.bgn < fraglist[i].position.end) ? (fraglist[i].position.end < fraglist[i].position.bgn) : (fraglist[i].position.bgn - fraglist[i].position.end);
    num_bases   += (int32)ceil(flen + 2 * AS_CNS_ERROR_RATE * flen);

    num_columns  = (fraglist[i].position.bgn > num_columns) ? fraglist[i].position.bgn : num_columns;
    num_columns  = (fraglist[i].position.end > num_columns) ? fraglist[i].position.end : num_columns;
  }

  ResetStores(num_bases, numfrags, num_columns);

  //  Magic initialization (in ResetStores()) prevents us calling CreateMANode() until now.

  trace    = CreateVA_int32(2 * AS_READ_MAX_NORMAL_LEN);
  traceBgn = 0;

  manode   = CreateMANode(ma->maID);
  utgpos   = (SeqInterval *)safe_calloc(numfrags, sizeof(SeqInterval));
  cnspos   = (SeqInterval *)safe_calloc(numfrags, sizeof(SeqInterval));

  assert(manode->lid == 0);

  frankensteinLen = 0;
  frankensteinMax = MAX(1024 * 1024, 2 * num_columns);
  frankenstein    = (char    *)safe_malloc(sizeof(char)    * frankensteinMax);
  frankensteinBof = (beadIdx *)safe_malloc(sizeof(beadIdx) * frankensteinMax);

  for (int32 i=0; i<numfrags; i++) {
    int32 complement = (fraglist[i].position.bgn < fraglist[i].position.end) ? 0 : 1;
    int32 fid;

    if (fraglist[i].type != AS_READ) {
      fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Fragment %d is not a read.\n",
              ma->maID, fraglist[i].ident);
      return(false);
    }

    if (HASH_SUCCESS != InsertInHashTable_AS(fragmentMap,fraglist[i].ident, 0, 1, 0)) {
      fprintf(stderr, "MultiAlignUnitig()-- Unitig %d FAILED.  Fragment %d is a duplicate.\n",
              ma->maID, fraglist[i].ident);
      return(false);
    }

    // This guy allocates and initializes the beads for each fragment.  Beads are not fully inserted
    // in the abacus here.

    fid = AppendFragToLocalStore(fraglist[i].type,
                                 fraglist[i].ident,
                                 complement,
                                 fraglist[i].contained,
                                 AS_OTHER_UNITIG);

    utgpos[fid].bgn = complement ? fraglist[i].position.end : fraglist[i].position.bgn;
    utgpos[fid].end = complement ? fraglist[i].position.bgn : fraglist[i].position.end;

    cnspos[fid].bgn  = 0;
    cnspos[fid].end  = 0;

    //  If this is violated, then the implicit map from utgpos[] and cnspos[] to fraglist is
    //  incorrect.
    assert(fid == i);

    //if (VERBOSE_MULTIALIGN_OUTPUT)
    //  fprintf(stderr,"MultiAlignUnitig()-- Added fragment mid %d pos %d,%d in unitig %d to store at local id %d.\n",
    //          fraglist[i].ident, fraglist[i].position.bgn, fraglist[i].position.end, ma->maID, fid);
  }

  SeedMAWithFragment(manode->lid, GetFragment(fragmentStore,0)->lid,0, opp);

  //  Save columns
  {
    beadIdx bidx = GetFragment(fragmentStore, 0)->firstbead;
    Bead   *bead = GetBead(beadStore, bidx);

    while (bead) {
      frankenstein   [frankensteinLen] = *Getchar(sequenceStore, bead->soffset);
      frankensteinBof[frankensteinLen] = bead->boffset;

      frankensteinLen++;

      bead = (bead->next.isInvalid()) ? NULL : GetBead(beadStore, bead->next);
    }

    frankenstein   [frankensteinLen] = 0;
    frankensteinBof[frankensteinLen] = beadIdx();

    cnspos[0].bgn = 0;
    cnspos[0].end = frankensteinLen;
  }

  return(true);
}



int
unitigConsensus::computePositionFromParent(bool doContained) {

  assert(piid == -1);

  AS_IID parent = (doContained == false) ? fraglist[tiid].parent : fraglist[tiid].contained;

  if (parent == 0)
    goto computePositionFromParentFail;

  if ((doContained == true) && (fraglist[tiid].parent == fraglist[tiid].contained))
    //  Already tried the parent, no need to try it again.
    goto computePositionFromParentFail;

  for (piid = tiid-1; piid >= 0; piid--) {
    Fragment *afrag = GetFragment(fragmentStore, piid);

    if (parent != afrag->iid)
      //  Not the parent.
      continue;

    if ((cnspos[piid].bgn == 0) &&
        (cnspos[piid].end == 0))
      //  Is the parent, but that isn't placed.
      goto computePositionFromParentFail;

    if ((utgpos[piid].end < utgpos[tiid].bgn) ||
        (utgpos[tiid].end < utgpos[piid].bgn)) {
      //  Is the parent, parent is placed, but the parent doesn't agree with the placement.
      if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_PLACEMENT)
        fprintf(stderr, "computePositionFromParent()-- parent %d at utg %d,%d doesn't agree with my utg %d,%d\n",
                parent,
                utgpos[piid].bgn, utgpos[piid].end,
                utgpos[tiid].bgn, utgpos[tiid].end);
      goto computePositionFromParentFail;
    }

    cnspos[tiid].bgn = cnspos[piid].bgn + fraglist[tiid].ahang;
    cnspos[tiid].end = cnspos[piid].end + fraglist[tiid].bhang;

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_PLACEMENT)
      fprintf(stderr, "computePositionFromParent()-- parent %d at %d,%d --> beg,end %d,%d (fLen %d)\n",
              parent,
              cnspos[piid].bgn, cnspos[piid].end,
              cnspos[tiid].bgn, cnspos[tiid].end,
              frankensteinLen);
    return(true);
  }

 computePositionFromParentFail:
  cnspos[tiid].bgn = 0;
  cnspos[tiid].end = 0;

  piid = -1;

  if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
    fprintf(stderr, "computePositionFromParent()-- Returns fail.\n");
  return(false);
}



int
unitigConsensus::computePositionFromLayout(void) {
  int32   thickestLen = 0;

  assert(piid == -1);

  //  Find the thickest qiid overlap to any cnspos fragment
  for (int32 qiid = tiid-1; qiid >= 0; qiid--) {
    if ((utgpos[tiid].bgn < utgpos[qiid].end) &&
        (utgpos[tiid].end > utgpos[qiid].bgn) &&
        ((cnspos[qiid].bgn != 0) ||
         (cnspos[qiid].end != 0))) {
      cnspos[tiid].bgn = cnspos[qiid].bgn + utgpos[tiid].bgn - utgpos[qiid].bgn;
      cnspos[tiid].end = cnspos[qiid].end + utgpos[tiid].end - utgpos[qiid].end;

      int32 ooo = MIN(cnspos[tiid].end, frankensteinLen) - cnspos[tiid].bgn;

      if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_PLACEMENT)
        fprintf(stderr, "computePositionFromLayout()-- layout %d at utg %d,%d cns %d,%d --> utg %d,%d cns %d,%d -- overlap %d\n",
                fraglist[qiid].ident,
                utgpos[qiid].bgn, utgpos[qiid].end, cnspos[qiid].bgn, cnspos[qiid].end,
                utgpos[tiid].bgn, utgpos[tiid].end, cnspos[tiid].bgn, cnspos[tiid].end,
                ooo);

      //  Occasionally we see an overlap in the original placement (utgpos overlap) by after
      //  adjusting our fragment to the frankenstein position, we no longer have an overlap.  This
      //  seems to be caused by a bad original placement.
      //
      //  Example:
      //  utgpos[a] = 13480,14239    cnspos[a] = 13622,14279
      //  utgpos[b] = 14180,15062
      //
      //  Our placement is 200bp different at the start, but close at the end.  When we compute the
      //  new start placement, it starts after the end of the A read -- the utgpos say the B read
      //  starts 700bp after the A read, which is position 13622 + 700 = 14322....50bp after A ends.

      if ((cnspos[tiid].bgn < frankensteinLen) &&
          (thickestLen < ooo)) {
        thickestLen = ooo;

        int32 ovl   = ooo;
        int32 ahang = cnspos[tiid].bgn;
        int32 bhang = cnspos[tiid].end - frankensteinLen;

        piid  = qiid;
      }
    }
  }

  //  If we have a VALID thickest placement, use that (recompute the placement that is likely
  //  overwritten -- ahang, bhang and piid are still correct).

  if (thickestLen >= AS_OVERLAP_MIN_LEN) {
    assert(piid != -1);

    cnspos[tiid].bgn = cnspos[piid].bgn + utgpos[tiid].bgn - utgpos[piid].bgn;
    cnspos[tiid].end = cnspos[piid].end + utgpos[tiid].end - utgpos[piid].end;

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_PLACEMENT)
      fprintf(stderr, "computePositionFromLayout()-- layout %d at %d,%d --> beg,end %d,%d (fLen %d)\n",
              fraglist[piid].ident,
              cnspos[piid].bgn, cnspos[piid].end,
              cnspos[tiid].bgn, cnspos[tiid].end,
              frankensteinLen);

    return(true);
  }

  cnspos[tiid].bgn = 0;
  cnspos[tiid].end = 0;

  piid = -1;

  if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
    fprintf(stderr, "computePositionFromLayout()-- Returns fail (not found or not thick enough).\n");
  return(false);
}



//  Occasionally we get a fragment that just refuses to go in the correct spot.  Search for the
//  correct placement in all of frankenstein, update ahang,bhang and retry.
//
//  We don't expect to have big negative ahangs, and so we don't allow them.  To unlimit this, use
//  "-fragmentLen" instead of the arbitrary cutoff below.
int
unitigConsensus::computePositionFromAlignment(void) {

  assert(piid == -1);

  ALNoverlap  *O           = NULL;
  double       thresh      = 1e-3;
  int32        minlen      = AS_OVERLAP_MIN_LEN;
  int32        ahanglimit  = -10;

  char        *fragment    = Getchar(sequenceStore, GetFragment(fragmentStore, tiid)->sequence);
  int32        fragmentLen = strlen(fragment);

  O = DP_Compare(frankenstein,
                 fragment,
                 ahanglimit, frankensteinLen,  //  ahang bounds
                 frankensteinLen, fragmentLen,   //  length of fragments
                 0,
                 AS_CNS_ERROR_RATE, thresh, minlen,
                 AS_FIND_ALIGN);

  if (O == NULL)
    O = Local_Overlap_AS_forCNS(frankenstein,
                                fragment,
                                ahanglimit, frankensteinLen,  //  ahang bounds
                                frankensteinLen, fragmentLen,   //  length of fragments
                                0,
                                AS_CNS_ERROR_RATE, thresh, minlen,
                                AS_FIND_ALIGN);

  if (O == NULL) {
    cnspos[tiid].bgn = 0;
    cnspos[tiid].end = 0;

    piid = -1;

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "computePositionFromAlignment()-- Returns fail.\n");
    return(false);
  }

  //  From the overlap and existing placements, find the thickest overlap, to set the piid and
  //  hangs, then reset the original placement based on that parents original placement.
  //
  //  To work with fixFailures(), we need to scan the entire fragment list.  This isn't so
  //  bad, really, since before we were scanning (on average) half of it.
  //
  cnspos[tiid].bgn = O->begpos;
  cnspos[tiid].end = O->endpos + frankensteinLen;
  //fprintf(stderr, "cnspos[%3d] mid %d %d,%d\n", tiid, fraglist[tiid].ident, cnspos[tiid].bgn, cnspos[tiid].end);

  int32   thickestLen = 0;

  for (int32 qiid = numfrags-1; qiid >= 0; qiid--) {
    if ((tiid != qiid) &&
        (cnspos[tiid].bgn < cnspos[qiid].end) &&
        (cnspos[tiid].end > cnspos[qiid].bgn)) {
      int32 ooo = (MIN(cnspos[tiid].end, cnspos[qiid].end) -
                   MAX(cnspos[tiid].bgn, cnspos[qiid].bgn));

      if (thickestLen < ooo) {
        thickestLen = ooo;

        int32 ovl   = ooo;
        int32 ahang = cnspos[tiid].bgn;
        int32 bhang = cnspos[tiid].end - frankensteinLen;

        piid  = qiid;
      }
    }
  }

  if (thickestLen > 0) {
    assert(piid != -1);

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_PLACEMENT)
      fprintf(stderr, "computePositionFromAlignment()-- layout %d at %d,%d --> beg,end %d,%d (fLen %d)\n",
              fraglist[piid].ident,
              cnspos[piid].bgn, cnspos[piid].end,
              cnspos[tiid].bgn, cnspos[tiid].end,
              frankensteinLen);

    return(true);
  }

  cnspos[tiid].bgn = 0;
  cnspos[tiid].end = 0;

  piid = -1;

  if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
    fprintf(stderr, "computePositionFromAlignment()-- Returns fail.\n");
  return(false);
}


void
unitigConsensus::rebuild(bool recomputeFullConsensus) {

  if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
    fprintf(stderr, "rebuild()--\n");

  if (recomputeFullConsensus == false) {
    //  For each column, vote for the most common base.
    int32 cid   = manode->first;
    int32 index = 0;
    
    Resetint32(manode->columnList);

    while (cid  > -1) {
      Column *column = GetColumn(columnStore, cid);
      int32   nA = GetColumnBaseCount(column, 'A');
      int32   nC = GetColumnBaseCount(column, 'C');
      int32   nG = GetColumnBaseCount(column, 'G');
      int32   nT = GetColumnBaseCount(column, 'T');
      int32   n_ = GetColumnBaseCount(column, '-');
      int32   nn = 0;

      Bead   *bead = GetBead(beadStore, column->call);
      char    call = '-';

      if (nA > nn) { nn = nA;  call = 'A'; }
      if (nC > nn) { nn = nC;  call = 'C'; }
      if (nG > nn) { nn = nG;  call = 'G'; }
      if (nT > nn) { nn = nT;  call = 'T'; }

      Setchar(sequenceStore, bead->soffset, &call);

      //  This is extracted from RefreshMANode()
      column->ma_index = index++;
      AppendVA_int32(manode->columnList, &cid);

      cid = column->next;
    }

  } else {
    //  Run abacus to rebuild an intermediate consensus sequence.  VERY expensive, and doesn't
    //  update the cnspos[] array...but the changes shouldn't be huge.
    RefreshMANode(manode->lid, 0, opp, NULL, NULL, 0, 0);

    //  Are all three needed??

    AbacusRefine(manode,0,-1,CNS_SMOOTH, opp);
    MergeRefine(manode->lid, NULL, 1, opp, 1);

    AbacusRefine(manode,0,-1,CNS_POLYX, opp);
    MergeRefine(manode->lid, NULL, 1, opp, 1);

    AbacusRefine(manode,0,-1,CNS_INDEL, opp);
    MergeRefine(manode->lid, NULL, 1, opp, 1);
  }

  //  Extract the consensus sequence.  Frankenstein is ponting to a consensus bead, not a fragment bead.

  ConsensusBeadIterator  bi;
  beadIdx                bid;

  int32                  gapToUngapLen = 0;

  CreateConsensusBeadIterator(manode->lid, &bi);

  frankensteinLen = 0;

  while ((bid = NextConsensusBead(&bi)) .isValid()) {
    Bead *bead = GetBead(beadStore, bid);
    char  cnsc = *Getchar(sequenceStore, bead->soffset);

    gapToUngapLen++;

    if (cnsc == '-')
      continue;

    while (frankensteinLen >= frankensteinMax) {
      frankensteinMax *= 2;
      frankenstein     = (char    *)safe_realloc(frankenstein,    sizeof(char)    * frankensteinMax);
      frankensteinBof  = (beadIdx *)safe_realloc(frankensteinBof, sizeof(beadIdx) * frankensteinMax);
    }
    assert(frankensteinLen < frankensteinMax);

    frankenstein   [frankensteinLen] = cnsc;
    frankensteinBof[frankensteinLen] = bead->boffset;
    frankensteinLen++;
  }

  frankenstein   [frankensteinLen] = 0;
  frankensteinBof[frankensteinLen] = beadIdx();

  //fprintf(stderr, "AFTER REBUILD %s\n", frankenstein);

  //  Update the positions.  This is the same method as used by GetMANodePositions to update the
  //  unitig f_list at the end, except we need to translate the ma_index from gapped to ungapped
  //  coordinates.

  int32  *gapToUngap = new int32 [gapToUngapLen];

  for (int32 i=0; i<gapToUngapLen; i++)
    gapToUngap[i] = -1;

  for (int32 i=0; i<frankensteinLen; i++) {
    Bead   *bead = GetBead(beadStore, frankensteinBof[i]);
    Column *col  = GetColumn(columnStore, bead->column_index);

    assert(col->ma_index >= 0);
    assert(col->ma_index < gapToUngapLen);

    gapToUngap[col->ma_index] = i;
  }

  for (int32 lg=0, i=0; i<gapToUngapLen; i++) {
    if (gapToUngap[i] == -1)
      gapToUngap[i] = lg;
    else
      lg = gapToUngap[i];
  }

  //  Update the position of each fragment in the consensus sequence.
  //  with, skip it.

  for (int32 i=0; i<=tiid; i++) {
    if ((cnspos[i].bgn == 0) &&
        (cnspos[i].end == 0))
      //  Uh oh, not placed originally.
      continue;

    Fragment *frg  = GetFragment(fragmentStore, i);
    Bead     *frst = GetBead(beadStore, frg->firstbead);
    Bead     *last = GetBead(beadStore, frg->firstbead.get() + frg->length - 1);

    int32     frstIdx = GetColumn(columnStore, frst->column_index)->ma_index;
    int32     lastIdx = GetColumn(columnStore, last->column_index)->ma_index;

    assert(frstIdx >= 0);
    assert(lastIdx >= 0);

    assert(frstIdx < gapToUngapLen);
    assert(lastIdx < gapToUngapLen);

    cnspos[i].bgn = gapToUngap[frstIdx];
    cnspos[i].end = gapToUngap[lastIdx] + 1;
  }

  //  Finally, update the parent/hang of the fragment we just placed.

  if (piid >= 0) {
    fraglist[tiid].parent    = fraglist[piid].ident;
    fraglist[tiid].ahang     = cnspos[tiid].bgn - cnspos[piid].bgn;
    fraglist[tiid].bhang     = cnspos[tiid].end - cnspos[piid].end;
    fraglist[tiid].contained = (fraglist[tiid].bhang > 0) ? 0 : fraglist[piid].ident;
    fraglist[tiid].contained = (fraglist[tiid].ahang < 0) ? 0 : fraglist[tiid].contained;
  }

  piid = -1;

  delete [] gapToUngap;
}


int
unitigConsensus::alignFragment(void) {
  double        origErate    = AS_CNS_ERROR_RATE;
  int32         bgnExtra     = 0;
  int32         endExtra     = 0;

  ALNoverlap  *O           = NULL;
  double       thresh      = 1e-3;
  int32        minlen      = AS_OVERLAP_MIN_LEN;

  assert((cnspos[tiid].bgn != 0) || (cnspos[tiid].end != 0));
  assert(piid != -1);

  //  Compute how much extra consensus sequence to align to.  If we allow too much, we risk placing
  //  the fragment in an incorrect location -- worst case is that we find a longer higher scoring
  //  alignment, but lower identity that is rejected.  If we allow too little, we leave unaligned
  //  bits on the end of the fragment -- worst case, we have 1 extra base which might add an
  //  unnecessary gap.
  //
  //  The default used to be 100bp, which made sense for Sanger reads (maybe) but makes no sense for
  //  short reads of length say 100bp.
  //
  //  We are placing this read relative to some other read:
  //
  //      anchoring fragment   -----------------------------------------    ediff
  //      new fragment           bdiff   ------------------------------------------------
  //
  //  So we should allow AS_CNS_ERROR_RATE indel in those relative positionings.
  //

  bgnExtra = (int32)ceil(AS_CNS_ERROR_RATE * (cnspos[tiid].bgn - cnspos[piid].bgn));
  endExtra = (int32)ceil(AS_CNS_ERROR_RATE * (cnspos[tiid].end - cnspos[piid].end));

  if (bgnExtra < 0)  bgnExtra = -bgnExtra;
  if (endExtra < 0)  endExtra = -endExtra;

  if (bgnExtra < 10) bgnExtra = 10;
  if (endExtra < 10) endExtra = 10;

 alignFragmentAgain:
  if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
    fprintf(stderr, "alignFragment()-- Allow bgnExtra=%d and endExtra=%d\n",
            bgnExtra, endExtra);

  int32 frankBgn     = MAX(0, cnspos[tiid].bgn - bgnExtra);   //  Start position in frankenstein
  int32 frankEnd     = frankensteinLen;                       //  Truncation of frankenstein
  char  frankEndBase = 0;                                     //  Saved base from frankenstein

  bool  allowBhang   = true;
  bool  tryAgain     = false;

  //  If the bhang + extra are negative then we need to truncate frankenstein so we don't align
  //  to the wrong spot.
  //
  if (cnspos[tiid].end + endExtra < frankEnd) {
    frankEnd     = cnspos[tiid].end + endExtra;
    frankEndBase = frankenstein[frankEnd];
    frankenstein[frankEnd] = 0;
    allowBhang   = false;
  }

  char      *aseq  = frankenstein + frankBgn;
  int32      alen  = frankEnd - frankBgn;

  char      *bseq  = Getchar(sequenceStore, GetFragment(fragmentStore, tiid)->sequence);  //  The fragment
  Fragment  *bfrag = GetFragment(fragmentStore, tiid);
  int32      blen  = bfrag->length;

  //  Now just fish for an alignment that is decent.  This is mostly straight from alignFragmentToFragment.

#if 0
  if (O == NULL) {
    O = DP_Compare(aseq,
                   bseq,
                   0, alen,            //  ahang bounds
                   alen, blen,         //  ahang, bhang exclusion are unused here
                   0,
                   AS_CNS_ERROR_RATE + 0.02, thresh, minlen,
                   AS_FIND_ALIGN);
    if ((O) && (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)) {
      fprintf(stderr, "DP_Compare found:\n");
      PrintALNoverlap(stderr, aseq, bseq, O);
    }
  }
  if (rejectAlignment(allowBhang, O))
    O = NULL;
#endif

#if 0
  if (O == NULL) {
    O = Local_Overlap_AS_forCNS(aseq,
                                bseq,
                                0, alen,            //  ahang bounds
                                alen, blen,         //  ahang, bhang exclusion are unused here
                                0,
                                AS_CNS_ERROR_RATE + 0.02, thresh, minlen,
                                AS_FIND_ALIGN);
    if ((O) && (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)) {
      fprintf(stderr, "Local_Overlap found:\n");
      PrintALNoverlap(stderr, aseq, bseq, O);
    }
  }
  if (rejectAlignment(allowBhang, O))
    O = NULL;
#endif

  if (O == NULL) {
    O = Optimal_Overlap_AS_forCNS(aseq,
                                  bseq,
                                  0, alen,            //  ahang bounds are unused here
                                  alen, blen,         //  ahang, bhang exclusion
                                  0,
                                  AS_CNS_ERROR_RATE + 0.02, thresh, minlen,
                                  AS_FIND_ALIGN);
    if ((O) && (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)) {
      fprintf(stderr, "Optimal_Overlap found:\n");
      PrintALNoverlap(stderr, aseq, bseq, O);
    }
  }
  if ((O) && (O->begpos < 0)) {
    bgnExtra += 2 * -O->begpos;
    tryAgain = true;
  }
  if ((O) && (O->endpos > 0) && (allowBhang == false)) {
    endExtra += 2 * O->endpos;
    tryAgain = true;
  }
  if (rejectAlignment(allowBhang, O))
    O = NULL;

  //  Restore the base we might have removed from frankenstein.
  if (frankEndBase)
    frankenstein[frankEnd] = frankEndBase;

  if (O) {
    Resetint32(trace);

    traceBgn = frankBgn + O->begpos;

    for (int32 *t = O->trace; (t != NULL) && (*t != 0); t++) {
      if (*t < 0)
        *t -= frankBgn;
      AppendVA_int32(trace, t);
    }

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "alignFragment()-- Alignment succeeded.\n");

    return(true);
  }

  if (tryAgain)
    goto alignFragmentAgain;

  //  No alignment.  Dang.
  cnspos[tiid].bgn = 0;
  cnspos[tiid].end = 0;

  piid = -1;

  if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
    fprintf(stderr, "alignFragment()-- No alignment found.\n");

  return(false);
}


bool
unitigConsensus::rejectAlignment(bool allowBhang, ALNoverlap *O) {

  if (O == NULL) {
    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "alignFragmentToFragments()-- No alignment found.\n");
    return(true);
  }

  //  Negative ahang?  Nope, don't want it.
  if (O->begpos < 0) {
    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "alignFragmentToFragments()-- No alignment found -- begpos = %d (negative ahang).\n", O->begpos);
    return(true);
  }

  //  Positive bhang and not the last fragment?  Nope, don't want it.
  if ((O->endpos > 0) && (allowBhang == false)) {
    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "alignFragmentToFragments()-- No alignment found -- endpos = %d (positive bhang not allowed).\n", O->endpos);
    return(true);
  }

  //  Too noisy?  Nope, don't want it.
  if (((double)O->diffs / (double)O->length) > AS_CNS_ERROR_RATE) {
    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "alignFragmentToFragments()-- No alignment found -- erate %f > max allowed %f.\n",
              (double)O->diffs / (double)O->length, AS_CNS_ERROR_RATE);
    return(true);
  }

  //  Too short?  Nope, don't want it.
  if (O->length < AS_OVERLAP_MIN_LEN) {
    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "alignFragmentToFragments()-- No alignment found -- too short %d < min allowed %d.\n",
              O->length, AS_OVERLAP_MIN_LEN);
    return(true);
  }

  return(false);
}




int
unitigConsensus::alignFragmentToFragments(void) {
  char     *fragment    = Getchar(sequenceStore, GetFragment(fragmentStore, tiid)->sequence);
  int32     fragmentLen = strlen(fragment);

  for (int32 qiid = tiid-1; qiid >= 0; qiid--) {
    ALNoverlap  *O           = NULL;
    double       thresh      = 1e-3;
    int32        minlen      = AS_OVERLAP_MIN_LEN;

    //  If the current fragment is not contained and the target fragment doesn't extend to the end
    //  of frankenstein, don't even bother aligning.  Any alignment we'd get would have to be
    //  contained else we'll insert huge gaps into the multialign, but the fragment isn't marked as
    //  contained.
    //
    if ((fraglist[tiid].contained == 0) &&
        (cnspos[qiid].end != frankensteinLen))
      continue;

    //  If the target fragment was skipped (failOnFirstFailure), we can't use it for alignment into
    //  the unitig.  Skip it.
    //
    if ((cnspos[qiid].bgn == 0) && (cnspos[qiid].end == 0))
      continue;

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "alignFragmentToFragments()--  Testing vs %d\n", fraglist[qiid].ident);

    char      *aseq = Getchar(sequenceStore, GetFragment(fragmentStore, qiid)->sequence);
    char      *bseq = Getchar(sequenceStore, GetFragment(fragmentStore, tiid)->sequence);

    int32      alen = GetFragment(fragmentStore, qiid)->length;
    int32      blen = GetFragment(fragmentStore, tiid)->length;

    bool       allowBhang = (cnspos[qiid].end == frankensteinLen);

    //  Lots of useless diagnostics here.

    //if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM) {
    //  fprintf(stderr, "A idx=%d id=%d len=%d %s\n", qiid, fraglist[qiid].ident, alen, aseq);
    //  fprintf(stderr, "B idx=%d id=%d len=%d %s\n", tiid, fraglist[tiid].ident, blen, bseq);
    //}

    //  Go fishing for an alignment.

    if (O == NULL) {
      O = DP_Compare(aseq,
                     bseq,
                     0, alen,            //  ahang bounds
                     alen, blen,         //  ahang, bhang exclusion are unused here
                     0,
                     AS_CNS_ERROR_RATE + 0.02, thresh, minlen,
                     AS_FIND_ALIGN);
      if ((O) && (VERBOSE_MULTIALIGN_OUTPUT)) {
        fprintf(stderr, "DP_Compare found:\n");
        PrintALNoverlap(stderr, aseq, bseq, O);
      }
    }
    if (rejectAlignment(allowBhang, O))
      O = NULL;

    if (O == NULL) {
      O = Local_Overlap_AS_forCNS(aseq,
                                  bseq,
                                  0, alen,            //  ahang bounds
                                  alen, blen,         //  ahang, bhang exclusion are unused here
                                  0,
                                  AS_CNS_ERROR_RATE + 0.02, thresh, minlen,
                                  AS_FIND_ALIGN);
      if ((O) && (VERBOSE_MULTIALIGN_OUTPUT)) {
        fprintf(stderr, "Local_Overlap found:\n");
        PrintALNoverlap(stderr, aseq, bseq, O);
      }
    }
    if (rejectAlignment(allowBhang, O))
      O = NULL;

    if (O == NULL) {
      O = Optimal_Overlap_AS_forCNS(aseq,
                                    bseq,
                                    0, alen,            //  ahang bounds are unused here
                                    alen, blen,         //  ahang, bhang exclusion
                                    0,
                                    AS_CNS_ERROR_RATE + 0.02, thresh, minlen,
                                    AS_FIND_ALIGN);
      if ((O) && (VERBOSE_MULTIALIGN_OUTPUT)) {
        fprintf(stderr, "Optimal_Overlap found:\n");
        PrintALNoverlap(stderr, aseq, bseq, O);
      }
    }
    if (rejectAlignment(allowBhang, O))
      continue;


    //  Make up plausible guesses for where this fragment was located.  This is needed by rebuild(),
    //  which expects all placed fragments to have at least a guess.

    cnspos[tiid].bgn = cnspos[qiid].bgn + O->begpos;
    cnspos[tiid].end = cnspos[qiid].end + O->endpos;

    piid = qiid;

    //  Add the alignment to abacus.

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_PLACEMENT)
      fprintf(stderr, "alignFragmentToFragments()-- parent %d at %d,%d --> beg,end %d,%d (fLen %d)\n",
              fraglist[qiid].ident,
              cnspos[qiid].bgn, cnspos[qiid].end,
              cnspos[tiid].bgn, cnspos[tiid].end,
              frankensteinLen);

    applyAlignment(qiid, O->begpos, O->trace);
    rebuild(false);
    return(true);
  }

  cnspos[tiid].bgn = 0;
  cnspos[tiid].end = 0;

  piid = -1;

  if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
    fprintf(stderr, "alignFragmentToFragments()-- Found no alignments.\n");
  return(false);
}




void
unitigConsensus::applyAlignment(int32 frag_aiid, int32 frag_ahang, int32 *frag_trace) {

  
  if (frag_aiid >= 0) {
    //  Aligned to a fragent
    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "applyAlignment()-- frag_aiid=%d frag_ahang=%d\n", frag_aiid, frag_ahang);
    ApplyAlignment(frag_aiid,
                   0, NULL,
                   tiid,
                   frag_ahang, frag_trace);

  } else {
    //  Aligned to frankenstein
    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "applyAlignment()-- frankenstein\n");
    ApplyAlignment(-1,
                   frankensteinLen, frankensteinBof,
                   tiid,
                   traceBgn, Getint32(trace, 0));
  }
}







void
unitigConsensus::generateConsensus(CNS_PrintKey     printwhat) {

  RefreshMANode(manode->lid, 0, opp, NULL, NULL, 0, 0);

  AbacusRefine(manode,0,-1,CNS_SMOOTH, opp);
  MergeRefine(manode->lid, NULL, 1, opp, 1);

  AbacusRefine(manode,0,-1,CNS_POLYX, opp);
  MergeRefine(manode->lid, NULL, 1, opp, 1);

  AbacusRefine(manode,0,-1,CNS_INDEL, opp);
  MergeRefine(manode->lid, NULL, 1, opp, 1);

  GetMANodeConsensus(manode->lid, ma->consensus, ma->quality);
  GetMANodePositions(manode->lid, ma);

  //  Although we generally don't care about delta values during assembly, we need them for the
  //  output, and this is the only time we compute them.  So, we've gotta hang on to them.
  //
  //for (int32 i=0; i<numfrags; i++) {
  //  fraglist[i].delta_length = 0;
  //  fraglist[i].delta        = NULL;
  //}

  //  Update or create the unitig in the MultiAlignT.

  if (GetNumIntUnitigPoss(ma->u_list) == 0) {
    IntUnitigPos  iup;

    iup.type           = AS_OTHER_UNITIG;
    iup.ident          = ma->maID;
    iup.position.bgn   = 0;
    iup.position.end   = GetMultiAlignLength(ma);
    iup.num_instances  = 0;
    iup.delta_length   = 0;
    iup.delta          = NULL;

    AppendIntUnitigPos(ma->u_list, &iup);
  } else {
    IntUnitigPos  *iup = GetIntUnitigPos(ma->u_list, 0);

    iup->position.bgn = 0;
    iup->position.end = GetMultiAlignLength(ma);
  }


  if ((printwhat == CNS_VERBOSE) ||
      (printwhat == CNS_VIEW_UNITIG))
    PrintAlignment(stderr,manode->lid,0,-1,printwhat);

  //  While we have fragments in memory, compute the microhet probability.  Ideally, this would be
  //  done in CGW when loading unitigs (the only place the probability is used) but the code wants
  //  to load sequence and quality for every fragment, and that's too expensive.
  {
    int32    depth  = 0;
    char **multia = NULL;
    int32  **id_array = NULL;

    MANode2Array(manode, &depth, &multia, &id_array,0);

    ma->data.unitig_microhet_prob = AS_REZ_MP_MicroHet_prob(multia, id_array, gkpStore, frankensteinLen, depth);

    for (int32 i=0;i<depth;i++) {
      safe_free(multia[2*i]);
      safe_free(multia[2*i+1]);
      safe_free(id_array[i]);
    }
    safe_free(multia);
    safe_free(id_array);
  }
}


bool
unitigConsensus::fixFailures(void) {
  int32  numFails = 0;

  for (tiid=0; tiid<numfrags; tiid++) {
    if ((cnspos[tiid].bgn == 0) && (cnspos[tiid].end == 0))
      numFails++;
  }

  for (tiid=0; tiid<numfrags; tiid++) {
    if ((cnspos[tiid].bgn != 0) || (cnspos[tiid].end != 0))
      //  Fragment already cnspos.
      continue;

    if (VERBOSE_MULTIALIGN_OUTPUT)
      reportStartingWork();

    if (computePositionFromParent(false) && alignFragment())  goto applyAlignmentAgain;
    if (computePositionFromParent(true)  && alignFragment())  goto applyAlignmentAgain;
    if (computePositionFromLayout()      && alignFragment())  goto applyAlignmentAgain;
    if (computePositionFromAlignment()   && alignFragment())  goto applyAlignmentAgain;

    //  Dang, still failed to align the fragment.

    return(false);

  applyAlignmentAgain:
    fprintf(stderr, "fixFailures()-- FIXING fragment %d from utgpos of %d,%d to cnspos of %d,%d\n",
            fraglist[tiid].ident,
            utgpos[tiid].bgn, utgpos[tiid].end,
            cnspos[tiid].bgn, cnspos[tiid].end);
  }

  //  Everything cnspos.  Update the ORIGINAL unitig placement.  We'll rerun it outside this
  //  function.

  for (int32 i=0; i<numfrags; i++) {
    if (fraglist[i].position.bgn < fraglist[i].position.end) {
      fraglist[i].position.bgn = cnspos[i].bgn;
      fraglist[i].position.end = cnspos[i].end;
    } else {
      fraglist[i].position.bgn = cnspos[i].end;
      fraglist[i].position.end = cnspos[i].bgn;
    }
  }

  //  Do a quick insertion sort to reorder the ORIGINAL unitig.

  for (int32 i=1; i<numfrags; i++) {
    IntMultiPos  Z;
    IntMultiPos  T = fraglist[i];
    int32        j = i-1;
    bool         C = true;

    memset(&Z, 0, sizeof(IntMultiPos));

    while (C) {
      if (MIN(fraglist[j].position.bgn, fraglist[j].position.end) > MIN(T.position.bgn, T.position.end)) {
        fraglist[j+1] = fraglist[j];
        fraglist[j]   = Z;
        j--;
      } else {
        C = false;
      }
      if (j < 0)
        C = false;
    }

    fraglist[j+1] = T;
  }

  //  Report what we have now.

#if 1
  for (int32 i=0; i<numfrags; i++)
    fprintf(stderr, "fixFailures()-- %d ident %d contained %d parent %d %d,%d position %d,%d\n",
            i,
            fraglist[i].ident,
            fraglist[i].contained,
            fraglist[i].parent,
            fraglist[i].ahang,
            fraglist[i].bhang,
            fraglist[i].position.bgn,
            fraglist[i].position.end);
#endif

  //  And return success!

  return(true);
}


int
MultiAlignUnitig(MultiAlignT     *ma,
                 gkStore         *fragStore,
                 CNS_PrintKey     printwhat,
                 CNS_Options     *opp) {
  double             origErate          = AS_CNS_ERROR_RATE;
  bool               failOnFirstFailure = false;
  bool               failuresToFix      = false;
  unitigConsensus   *uc                 = NULL;

 tryAgain:
  origErate          = AS_CNS_ERROR_RATE;
  uc                 = new unitigConsensus(ma, opp);

  if (uc->initialize() == FALSE)
    goto returnFailure;

  while (uc->moreFragments()) {
    if (VERBOSE_MULTIALIGN_OUTPUT)
      uc->reportStartingWork();

    if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    if (uc->alignFragmentToFragments())
      continue;

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "MultiAlignUnitig()-- recompute full consensus\n");
    uc->rebuild(true);

    if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    if (uc->alignFragmentToFragments())
      continue;

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "MultiAlignUnitig()-- increase allowed error rate from %f to %f\n", AS_CNS_ERROR_RATE, MIN(AS_MAX_ERROR_RATE, 1.3333 * AS_CNS_ERROR_RATE));

    AS_CNS_ERROR_RATE = MIN(AS_MAX_ERROR_RATE, 1.3333 * AS_CNS_ERROR_RATE);

    if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    if (uc->alignFragmentToFragments())
      continue;

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "MultiAlignUnitig()-- increase allowed error rate from %f to %f\n", AS_CNS_ERROR_RATE, MIN(AS_MAX_ERROR_RATE, 2.0 * AS_CNS_ERROR_RATE));

    AS_CNS_ERROR_RATE = MIN(AS_MAX_ERROR_RATE, 2.0 * AS_CNS_ERROR_RATE);

    if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    if (uc->alignFragmentToFragments())
      continue;

    //  Failed to align the fragment.  Dang.  Either fail immediately, or keep going expecting to
    //  fix things later.

    uc->reportFailure();

    if (failOnFirstFailure)
      goto returnFailure;

    failuresToFix = true;
    continue;

  applyAlignment:
    AS_CNS_ERROR_RATE = origErate;

    uc->applyAlignment();
    uc->rebuild(false);
  }

#if 1
  if (failuresToFix) {
    fprintf(stderr, "MultiAlignUnitig()-- WARNING!  Attempting to resolve alignment failures for unitig %d.\n", ma->maID);

    if (uc->fixFailures() == false) {
      fprintf(stderr, "MultiAlignUnitig()-- ERROR:  Failed to resolve alignment failures for unitig %d.\n", ma->maID);
      goto returnFailure;
    }

    fprintf(stderr, "MultiAlignUnitig()-- WARNING!  Failures resolved for unitig %d.  Attempt consensus again.\n", ma->maID);

    //  Fixed the placements.  Try again.

    delete uc;

    failOnFirstFailure = true;
    failuresToFix      = false;

    goto tryAgain;
  }
#endif

  uc->generateConsensus(printwhat);

  delete uc;
  return(TRUE);

 returnFailure:
  delete uc;
  return(FALSE);
}
