
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

static char *rcsid = "$Id: MultiAlignment_CNS.c,v 1.256 2011-01-03 03:07:16 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"


//
// Persistent store of the fragment data (produced upstream)
//
gkStore               *gkpStore      = NULL;
OverlapStore          *ovlStore      = NULL;
MultiAlignStore       *tigStore      = NULL;

HashTable_AS          *fragmentMap   = NULL;


//
// Stores for the sequence/quality/alignment information
// (reset after each multialignment)
//
VA_TYPE(char) *sequenceStore = NULL;
VA_TYPE(char) *qualityStore  = NULL;
VA_TYPE(Bead) *beadStore     = NULL;

//
// Local stores for
//      fragment information:
//                indices into sequence/quality stores
//                index into "bead" store for alignment information
//      column information:
//                basecall, profile, index in multialignment
//                indexed pointers to next and previous columns
//      multialignment information:
//                first and last column, profile, index in multialignment
//                VA of all component columns
//
// (All are reset after each multialignment)
//
VA_TYPE(Fragment) *fragmentStore = NULL;
VA_TYPE(Column)   *columnStore   = NULL;
VA_TYPE(MANode)   *manodeStore   = NULL;

int32 thisIsConsensus = 0;

//
// Convenience arrays for misc. fragment information
// (All are reset after each multialignment)
//
VA_TYPE(int32) *fragment_indices  = NULL;
VA_TYPE(int32) *abacus_indices    = NULL;

VA_TYPE(CNS_AlignedContigElement) *fragment_positions = NULL;

int64 gaps_in_alignment = 0;

int32 allow_neg_hang         = 0;


// Variables used to compute general statistics

int32 NumColumnsInUnitigs = 0;
int32 NumRunsOfGapsInUnitigReads = 0;
int32 NumGapsInUnitigs = 0;
int32 NumColumnsInContigs = 0;
int32 NumRunsOfGapsInContigReads = 0;
int32 NumGapsInContigs = 0;
int32 NumAAMismatches = 0; // mismatches b/w consensi of two different alleles
int32 NumVARRecords = 0;
int32 NumVARStringsWithFlankingGaps = 0;
int32 NumUnitigRetrySuccess = 0;
int32 contig_id = 0;

//
//  Tables to facilitate SNP Basecalling
//

double EPROB[CNS_MAX_QV-CNS_MIN_QV+1] = {0};  // prob of error for each quality value
double PROB[CNS_MAX_QV-CNS_MIN_QV+1]  = {0};  // prob of correct call for each quality value (1-eprob)
int32    RINDEX[RINDEXMAX]              = {0};

char ALPHABET[6] = {'-','a','c','g','t','n'};

char RALPHABET[CNS_NP] = {'-','A','C','G','T','N',
                          'a','c','g','t',   // -A, -C, -G, -T
                          'M','R','W',   //     AC, AG, AT
                          'S','Y',   //         CG, CT
                          'K',   //             GT
                          'm','r','w',   //    -AC,-AG,-AT
                          's','y',   //        -CG,-CT
                          'k',   //            -GT
                          'V','H','D','B',   //ACG,ACT,AGT,CGT
                          'v','h','d','b',   //-ACG,-ACT,-AGT,-CGT
                          'X','x'};// ACGT, -ACGT, ??

char RALPHABETC[CNS_NP] = {'-','T','G','C','A','N',
                           't','g','c','a',   // -A, -C, -G, -T
                           'K','Y','W',   //     AC, AG, AT
                           'S','R',   //         CG, CT
                           'M',   //             GT
                           'k','y','w',   //    -AC,-AG,-AT
                           's','r',   //        -CG,-CT
                           'm',   //            -GT
                           'B','D','H','V',   //ACG,ACT,AGT,CGT
                           'b','d','h','v',   //-ACG,-ACT,-AGT,-CGT
                           'X','x'};// ACGT,-ACGT, ??

double TAU_MISMATCH = (double)1./(5. - 1.);

uint32 AMASK[5] = {013607700741, // -
                   015670707042, // a
                   016733131104, // c
                   017355252210, // g
                   017566464420};



//  Define this to dump multifasta to stderr of the unitigs
//  we try to align in MultiAlignContig().  Was useful for
//  debugging bad layout.
//
int32 DUMP_UNITIGS_IN_MULTIALIGNCONTIG = 0;

// Be noisy when doing multi alignments - this used to be a #ifdef,
// which made it difficult to switch on in the middle of a debug.
//
int32 VERBOSE_MULTIALIGN_OUTPUT = 0;

//  If non-zero, we'll force-abut unitigs that don't align together.
//  Typically, these are caused by microscopic overlaps between
//  unitigs -- certainly less than 20bp long.
//
int32 FORCE_UNITIG_ABUT = 0;


//  This is called in ResetStores -- which is called before any
//  consensus work is done.
static
void
InitializeAlphTable(void) {
  int32 i;

  if (RINDEX[0] == 31)
    return;

  for (i=0; i<RINDEXMAX; i++)
    RINDEX[i] = 31;

  for(i=0; i<CNS_NP; i++)
    RINDEX[(int)RALPHABET[i]] = i;

  int32 qv=CNS_MIN_QV;

  for (i=0; i<CNS_MAX_QV-CNS_MIN_QV+1; i++) {
    EPROB[i]= pow(10, -qv/10.);
    PROB[i] = (1.0 - EPROB[i]);
    qv++;
  }
}



////////////////////////////////////////
//
//  Base Count data structure
//
//external
int
IncBaseCount(BaseCount *b, char c) {
  int32 i= RINDEX[c];
  if (c == 'N' || c == 'n' ) i=5;
  b->depth++;
  if( i<0 || i>5 ){
    fprintf(stderr, "IncBaseCount i out of range (possibly non ACGTN letter?)");    
    assert(0);
  }
  return b->count[i]++;
}

//external
int
DecBaseCount(BaseCount *b, char c) {
  int32 i= RINDEX[c];
  if (c == 'N' || c == 'n' ) i=5;
  b->depth--;
  if( i<0 || i>5 ){
    fprintf(stderr, "DecBaseCount i out of range");
    assert(0);
  }
  return b->count[i]--;
}

//external
int
GetBaseCount(BaseCount *b, char c) {
  int32 i= RINDEX[c];
  if (c == 'N' || c == 'n' ) i=5;
  return b->count[i];
}

//external
int
GetColumnBaseCount(Column *b, char c) {
  return GetBaseCount(&b->base_count,c);
}

//external
int
GetDepth(Column *c) {
  return c->base_count.depth;
}

//external
void
ResetBaseCount(BaseCount *b) {
  memset(b,'\0',sizeof(BaseCount));
}

//external
void
ShowBaseCount(BaseCount *b) {
  int32 i;
  fprintf(stderr,"%d total\n",b->depth);
  for (i=0;i<CNS_NALPHABET;i++) {
    fprintf(stderr,"%c\t",ALPHABET[i]);
  }
  fprintf(stderr,"\n");
  for (i=0;i<CNS_NALPHABET;i++) {
    fprintf(stderr,"%d\t",b->count[i]);
  }
  fprintf(stderr,"\n");
}

//external
void
ShowBaseCountPlain(FILE *out,BaseCount *b) {
  int32 i;
  fprintf(out,"%d\t",b->depth);
  for (i=0;i<CNS_NALPHABET;i++) {
    fprintf(out,"%d\t",b->count[i]);
  }
}

//external
char
GetMaxBaseCount(BaseCount *b,int32 start_index) {  // start at 1 to disallow gap
  int32 max_index = start_index,i;
  int32 tied = 0,tie_breaker,max_tie=0;
  for (i=start_index;i<CNS_NALPHABET-1;i++) {
    if (b->count[i] > b->count[max_index] ) {
      max_index = i;
      tied = 0;
    } else if ( b->count[i] == b->count[max_index]) {
      tied++;
    }
  }
  if ( tied > 1 ) {
    for (i=1;i<CNS_NALPHABET-1;i++) { /* i starts at 1 to prevent ties */
      /* from being broken with '-'    */
      if ( b->count[i] == b->count[max_index] ) {
        /* Break unresolved ties with random numbers: */
        tie_breaker = random();
        if (tie_breaker > max_tie) {
          max_tie = tie_breaker;
          max_index = i;
        }
      }
    }
  }
  return toupper(ALPHABET[max_index]);
}

//external
void
CheckColumnBaseCount(Column *c) {
  int32 counts[256] = {0};

  if (c->next == -1)
    return;

  Bead *cbead = GetBead(beadStore,c->call);

  while (cbead->down.isValid()) {
    cbead = GetBead(beadStore,cbead->down);
    counts[*Getchar(sequenceStore,cbead->soffset)]++;
  }

  if (counts['A'] != GetColumnBaseCount(c, 'A'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d A %d != %d\n", c->lid, counts['A'], GetColumnBaseCount(c, 'A'));
  if (counts['C'] != GetColumnBaseCount(c, 'C'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d C %d != %d\n", c->lid, counts['C'], GetColumnBaseCount(c, 'C'));
  if (counts['G'] != GetColumnBaseCount(c, 'G'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d G %d != %d\n", c->lid, counts['G'], GetColumnBaseCount(c, 'G'));
  if (counts['T'] != GetColumnBaseCount(c, 'T'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d T %d != %d\n", c->lid, counts['T'], GetColumnBaseCount(c, 'T'));
  if (counts['-'] != GetColumnBaseCount(c, '-'))
    fprintf(stderr, "CheckColumnBaseCount()-- cid=%d - %d != %d\n", c->lid, counts['-'], GetColumnBaseCount(c, '-'));
}



////////////////////////////////////////
//
//  Basic MANode operations
//
//external
MANode *
CreateMANode(int32 iid) {
  MANode ma;
  ma.lid     = GetNumMANodes(manodeStore);
  ma.iid     = iid;
  ma.first   = -1;
  ma.last    = -1;
  ma.columns = CreateVA_int32(GetAllocatedColumns(columnStore));

  AppendVA_MANode(manodeStore, &ma);

  return(GetMANode(manodeStore, ma.lid));
}

//external
void
DeleteMANode(int32 iid) {
  // Columns are in the columnStore, not under our control
  DeleteVA_int32(GetMANode(manodeStore,iid)->columns);
}

//external
int32
GetMANodeLength(int32 mid) {
  MANode *ma = GetMANode(manodeStore,mid);
  if ((ma) == NULL) return -1;
  return GetNumint32s(ma->columns);
}


//external
void
SeedMAWithFragment(int32 mid,
                   int32 fid,
                   int32 quality,
                   CNS_Options *opp) {
  Fragment *fragment = GetFragment(fragmentStore,fid);
  assert(fragment != NULL);

  MANode *ma = GetMANode(manodeStore, mid);
  assert(ma != NULL);

  FragmentBeadIterator fi;

  CreateFragmentBeadIterator(fid, &fi);

  // bid is the offset of the Bead seeding the column
  beadIdx bid = NextFragmentBead(&fi);

  Bead   *bead   = GetBead(beadStore,bid);
  Column *column = CreateColumn(bid);

  column->ma_id    =  mid;
  column->ma_index =  0;

  AddColumnToMANode(mid,*column);

  int32 cid = column->lid;

  while ( (bid = NextFragmentBead(&fi)) .isValid() )
    cid = ColumnAppend(cid, bid);

  fragment->manode=mid;

  if (quality > 0)
    RefreshMANode(mid, quality, opp, NULL, NULL, 1, 0);
  else
    RefreshMANode(mid, quality, opp, NULL, NULL, 0, 0);
}

//external
int
GetMANodeConsensus(int32 mid, VA_TYPE(char) *sequence, VA_TYPE(char) *quality) {
  ConsensusBeadIterator bi;
  Bead   *bead;
  beadIdx  bid;
  int32   length=GetMANodeLength(mid);
  int32   i=0;

  ResetVA_char(sequence);
  EnableRangeVA_char(sequence,length+1);

  ResetVA_char(quality);
  EnableRangeVA_char(quality,length+1);

  CreateConsensusBeadIterator(mid, &bi);

  while ( (bid = NextConsensusBead(&bi)) .isValid() ) {
    bead = GetBead(beadStore,bid);
    SetVA_char(sequence,i,Getchar(sequenceStore,bead->soffset));
    SetVA_char(quality,i,Getchar(qualityStore,bead->soffset));
    i++;
  }
  return length;
}


//  Used in GetMANodePositions
static
int32
GetFragmentDeltas(int32 fid, VA_TYPE(int32) *deltas, int32 length) {
  beadIdx               bid;
  FragmentBeadIterator fi;
  int32                index = 0;
  int32                added = 0;

  CreateFragmentBeadIterator(fid, &fi);

  //  index < length eliminates any endgaps from the delta list KAR, 09/19/02

  while (((bid = NextFragmentBead(&fi)) .isValid()) && (index < length)) {
    if (*Getchar(sequenceStore, GetBead(beadStore,bid)->soffset) == '-') {
      Appendint32(deltas, &index);
      added++;
    } else {
      index++;
    }
  }

  return(added);
}

//external
int
GetMANodePositions(int32        mid,
                   MultiAlignT *ma) {

  MANode *manode = GetMANode(manodeStore, mid);

  int32  n_frags   = 0;
  int32  n_unitigs = 0;

  ResetVA_int32(ma->fdelta);
  ResetVA_int32(ma->udelta);

  for (uint32 i=0; i<GetNumFragments(fragmentStore); i++) {
    Fragment *fragment = GetFragment(fragmentStore, i);

    //fprintf(stderr, "GetMANodePositions()--  frag %d ident %d deleted %d\n",
    //        i, fragment->iid, fragment->deleted);

    if (fragment->deleted)
      continue;

    assert(fragment->manode == mid);

    int32 bgn = GetColumn(columnStore, (GetBead(beadStore,fragment->firstbead.get()                       ))->column_index)->ma_index;
    int32 end = GetColumn(columnStore, (GetBead(beadStore,fragment->firstbead.get() + fragment->length - 1))->column_index)->ma_index + 1;

    if (fragment->type == AS_READ) {
      if (FALSE == ExistsInHashTable_AS (fragmentMap, fragment->iid, 0))
        //  Fragment is not in the contig f_list; is in a surrogate.
        continue;

      if (1 != LookupValueInHashTable_AS (fragmentMap, fragment->iid, 0))
        //  Attempting to place a surrogate fragment more than once.
        continue;

      //  Indicate we've placed the fragment.
      ReplaceInHashTable_AS(fragmentMap, fragment->iid, 0, 2, 0);

      IntMultiPos *imp = GetIntMultiPos(ma->f_list, n_frags++);

      imp->ident        = fragment->iid;
      imp->type         = fragment->type;
      imp->position.bgn = (fragment->complement) ? end : bgn;
      imp->position.end = (fragment->complement) ? bgn : end;
      imp->delta        = NULL;
      imp->delta_length = GetFragmentDeltas(i, ma->fdelta, fragment->length);
    }

    if (fragment->type  == AS_UNITIG) {
      IntUnitigPos  *iup = GetIntUnitigPos(ma->u_list, n_unitigs++);

      assert(iup->ident == fragment->iid);

      iup->position.bgn = (fragment->complement) ? end : bgn;
      iup->position.end = (fragment->complement) ? bgn : end;
      iup->delta        = NULL;
      iup->delta_length = GetFragmentDeltas(i, ma->udelta, fragment->length);
    }
  }

  //  Because contig consensus might have ejected fragments that don't align, the new list can be
  //  shorter than the original list.
  //
  ResetToRangeVA_IntMultiPos(ma->f_list, n_frags);

  //  Set delta pointers into the VA.

  int32  fdeltapos = 0;
  int32  udeltapos = 0;

  n_frags   = 0;
  n_unitigs = 0;

  for (uint32 i=0; i<GetNumFragments(fragmentStore); i++) {
    Fragment *fragment = GetFragment(fragmentStore, i);

    if (fragment->deleted)
      continue;

    if (fragment->type == AS_READ) {
      if (FALSE == ExistsInHashTable_AS(fragmentMap, fragment->iid, 0))
        continue;

      // all of the contig's fragments should've had their value set to 2 in previous block

      assert(2 == LookupValueInHashTable_AS(fragmentMap, fragment->iid, 0));
      DeleteFromHashTable_AS(fragmentMap, fragment->iid, 0);

      IntMultiPos *imp = GetIntMultiPos(ma->f_list, n_frags++);

      imp->delta = (imp->delta_length == 0) ? NULL : Getint32(ma->fdelta, fdeltapos);

      fdeltapos += imp->delta_length;
    }

    if (fragment->type == AS_UNITIG) {
      IntUnitigPos  *iup = GetIntUnitigPos(ma->u_list, n_unitigs++);

      iup->delta = (iup->delta_length == 0) ? NULL : Getint32(ma->udelta, udeltapos);

      udeltapos += iup->delta_length;
    }
  }

  return n_frags;
}



////////////////////////////////////////
//
//  Iterators over Abacus structures
//
//external
void
CreateColumnBeadIterator(int32 cid,ColumnBeadIterator *bi) {
  Column *column = GetColumn(columnStore,cid);
  assert(column != NULL);
  bi->column = *column;
  bi->bead = bi->column.call;
}

//external
beadIdx
NextColumnBead(ColumnBeadIterator *bi) {
  beadIdx nid;

  if (bi->bead.isValid()) {
    Bead *bead = GetBead(beadStore, bi->bead);
    nid = bead->down;
    bi->bead = nid;
  }
  return nid;
}



//external
void
NullifyFragmentBeadIterator(FragmentBeadIterator *bi) {
  bi->fragment = *GetFragment(fragmentStore,0);
  bi->bead     = beadIdx();
  bi->isNull   = true;
}

//external
int32
IsNULLIterator(FragmentBeadIterator *bi) {
  return(bi->isNull);
}

//external
void
CreateFragmentBeadIterator(int32 fid,FragmentBeadIterator *bi) {
  Fragment *fragment = GetFragment(fragmentStore,fid);
  assert(fragment != NULL);
  bi->fragment = *fragment;
  bi->bead     = bi->fragment.firstbead;
  bi->isNull   = false;
}

//external
beadIdx
NextFragmentBead(FragmentBeadIterator *bi) {
  beadIdx nid;
  assert(bi->isNull == false);
  if (bi->bead.isValid()) {
    Bead *bead = GetBead(beadStore, bi->bead);
    nid = bead->boffset;
    bi->bead = bead->next;
  }
  return nid;
}



//external
void
CreateConsensusBeadIterator(int32 mid,ConsensusBeadIterator *bi) {
  Column *first = GetColumn(columnStore,(GetMANode(manodeStore,mid))->first);
  assert(first != NULL);
  bi->manode_id = mid;
  bi->bead = first->call;
}

//external
beadIdx
NextConsensusBead(ConsensusBeadIterator *bi) {
  beadIdx nid;
  if (bi->bead.isValid()) {
    Bead *bead = GetBead(beadStore, bi->bead);
    nid = bead->boffset;
    bi->bead = bead->next;
  }
  return nid;
}




////////////////////////////////////////
//
//  Bead Manipulation
//
//external
void
ClearBead(beadIdx bid) {
  Bead *b = GetBead(beadStore,bid);

  b->boffset      = beadIdx();
  b->soffset      = seqIdx();
  b->foffset      = -1;
  b->prev         = beadIdx();
  b->next         = beadIdx();
  b->up           = beadIdx();
  b->down         = beadIdx();
  b->frag_index   = -1;
  b->column_index = -1;
}


//external
void
AlignBeadToColumn(int32 cid, beadIdx bid, char *label) {
  Column *column=GetColumn(columnStore,cid);
  Bead *call = GetBead(beadStore,column->call);
  Bead *first = GetBead(beadStore,call->down);
  Bead *align = GetBead(beadStore,bid);

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "AlignBeadToColumn()-- %s frag=%d bead=%d,%c moving from column=%d to column=%d\n", label, align->frag_index, bid, *Getchar(sequenceStore,align->soffset), align->column_index, cid);
#endif

  align->down         = first->boffset;
  align->up           = call->boffset;
  call->down          = align->boffset;
  first->up           = align->boffset;
  align->column_index = cid;

  IncBaseCount(&column->base_count,*Getchar(sequenceStore,align->soffset));
}


// remove bid from it's column, returning the next bead up in the column
//
//external
beadIdx
UnAlignBeadFromColumn(beadIdx bid) {
  Bead *bead = GetBead(beadStore,bid);

  if (bead->column_index == -1)
    return beadIdx();

  Column *column = GetColumn(columnStore,bead->column_index);
  Bead   *upbead = GetBead(beadStore,bead->up);
  char    bchar  = *Getchar(sequenceStore,bead->soffset);

  upbead->down = bead->down;

  if (bead->down.isValid() )
    GetBead(beadStore, bead->down)->up = upbead->boffset;

  DecBaseCount(&column->base_count,bchar);

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "UnAlignBeadFromColumn()-- frag=%d bead=%d leaving column=%d\n", bead->frag_index, bead->boffset, bead->column_index);
#endif

  bead->up   = beadIdx();
  bead->down = beadIdx();
  bead->column_index = -1;

  return upbead->boffset;
}


//external
beadIdx
UnAlignTrailingGapBeads(beadIdx bid) {
  // remove bid from it's column, returning the prev or next bead in the fragment
  Bead *bead = GetBead(beadStore,bid);
  Bead *upbead,*prevbead,*nextbead;
  beadIdx anchor;
  Column *column;
  char bchar;

  // find direction to remove
  anchor = bead->prev;
  while ( bead->next.isValid() && *Getchar(sequenceStore,(GetBead(beadStore,bead->next))->soffset) == '-' ) {
    bead = GetBead(beadStore,bead->next);
  }
  if (bead->next.isValid() ) {
    anchor = bead->next;
    while (bead->prev.isValid() && *Getchar(sequenceStore,(GetBead(beadStore,bead->prev))->soffset) == '-' ) {
      bead = GetBead(beadStore,bead->prev);
    }
  }
  while ( bead->boffset != anchor) {
    column = GetColumn(columnStore,bead->column_index);
    upbead = GetBead(beadStore,bead->up);
    bchar = *Getchar(sequenceStore,bead->soffset);
    if( bchar != '-'){
      fprintf(stderr, "UnAlignTrailingGapBead bchar is not a gap");
      assert(0);
    }
    upbead->down = bead->down;
    if (bead->down.isValid() ) {
      GetBead(beadStore, bead->down)->up = upbead->boffset;
    }
    DecBaseCount(&column->base_count,bchar);

#ifdef DEBUG_ABACUS_ALIGN
    fprintf(stderr, "UnAlignTrailingGapBeads()-- frag=%d bead=%d leaving column=%d\n", bead->frag_index, bead->boffset, bead->column_index);
#endif

    bead->up   = beadIdx();
    bead->down = beadIdx();
    bead->column_index = -1;
    if ( bead->next.isInvalid() ) {
      prevbead = GetBead(beadStore,bead->prev);
      prevbead->next = beadIdx();
      bead->prev = beadIdx();
      bead = GetBead(beadStore,prevbead->boffset);
    } else {
      nextbead = GetBead(beadStore,bead->next);
      nextbead->prev = beadIdx();
      bead->next = beadIdx();
      bead = GetBead(beadStore,nextbead->boffset);
    }
  }
  return anchor;
}


//external
void
LateralExchangeBead(beadIdx lid, beadIdx rid) {
  Bead rtmp; // this is just some tmp space for the swap

  //  This function swaps the contents of two beads, ensuring that
  //  there are only gaps between them.
  //
  //  HORRIBLY complicated because ApplyAbacus() and MergeCompatible()
  //  hold on to pointers to beads.  It would have been much simpler
  //  to just swap the soffset and foffset, leaving EVERYTHING ELSE
  //  exactly the same.

  Bead *leftbead = GetBead(beadStore,lid);
  Bead *rightbead = GetBead(beadStore,rid);

  Column *leftcolumn = GetColumn(columnStore,leftbead->column_index);
  Column *rightcolumn = GetColumn(columnStore,rightbead->column_index);

  char leftchar = *Getchar(sequenceStore,leftbead->soffset);
  char rightchar = *Getchar(sequenceStore,rightbead->soffset);

  // now, verify that left and right are either
  // a) neighbors, or b) have only '-'s intervening

  {
    Bead *ibead   = leftbead;
    int32   failure = 0;
    int32   limit   = 20;

    while (ibead->next.isValid()) {
      ibead = GetBead(beadStore,ibead->next);

      if (ibead->boffset == rid)
        break;

      if( *Getchar(sequenceStore,ibead->soffset) != '-')
        failure++;
    }

    if (failure) {
      ibead = leftbead;

      while (ibead->next.isValid()) {
        ibead = GetBead(beadStore,ibead->next);

        fprintf(stderr, "bead %c boffset="F_U64" prev="F_U64" next="F_U64" up="F_U64" down="F_U64" fragindex=%d colulmnindex=%d\n",
                *Getchar(sequenceStore,ibead->soffset),
                (uint64)ibead->boffset.get(),
                (uint64)ibead->prev.get(),
                (uint64)ibead->next.get(),
                (uint64)ibead->up.get(),
                (uint64)ibead->down.get(),
                ibead->frag_index,
                ibead->column_index);

        if (ibead->boffset == rid)
          break;

        if (limit-- == 0)
          break;
      }

      fprintf(stderr, "LateralExchangeBead can't exchange bead "F_U64" with "F_U64"; bases in between!\n",
              (uint64)lid.get(),
              (uint64)rid.get());
      assert(failure == 0);
    }
  }

  rtmp = *rightbead;

  rightbead->up = leftbead->up;
  rightbead->down = leftbead->down;
  rightbead->prev = leftbead->prev;
  rightbead->next = leftbead->next;
  if ( rightbead->up.isValid() ) (GetBead(beadStore,rightbead->up))->down = rid;
  if ( rightbead->down.isValid())  (GetBead(beadStore,rightbead->down))->up = rid;
  if ( rightbead->prev.isValid())  (GetBead(beadStore,rightbead->prev))->next = rid;

  leftbead->up = rtmp.up;
  leftbead->down = rtmp.down;
  leftbead->next = rtmp.next;
  leftbead->prev = rtmp.prev;
  if ( leftbead->up.isValid() ) (GetBead(beadStore,leftbead->up))->down = lid;
  if ( leftbead->down.isValid())  (GetBead(beadStore,leftbead->down))->up = lid;
  if ( leftbead->next.isValid())  (GetBead(beadStore,leftbead->next))->prev = lid;

  // now, handle separately cases of a) left and right are adjacent, and b) gaps intervene
  if ( rtmp.prev == lid) {
    rightbead->next = lid;
    leftbead->prev = rid;
  } else {
    if ( rightbead->next.isValid())  (GetBead(beadStore,rightbead->next))->prev = rid;
    if ( leftbead->prev.isValid())  (GetBead(beadStore,leftbead->prev))->next = lid;
  }

  rightbead->column_index = leftbead->column_index;
  leftbead->column_index = rtmp.column_index;

  // change basecounts for affected columns
  DecBaseCount(&leftcolumn->base_count,leftchar);
  IncBaseCount(&leftcolumn->base_count,rightchar);
  DecBaseCount(&rightcolumn->base_count,rightchar);
  IncBaseCount(&rightcolumn->base_count,leftchar);
}



//external
beadIdx
AppendGapBead(beadIdx bid) {
  // The gap will appear immediately following bid
  Bead *prev = GetBead(beadStore,bid);
  Bead bead;
  char base='-';
  char qv;

  bead.boffset.set(GetNumBeads(beadStore));
  bead.soffset.set(GetNumchars(sequenceStore));
  bead.foffset = prev->foffset+1;
  bead.up = beadIdx();
  bead.down = beadIdx();
  bead.frag_index = prev->frag_index;
  bead.column_index = -1;
  bead.next = prev->next;
  bead.prev = prev->boffset;
  prev->next = bead.boffset;
  qv = *Getchar(qualityStore,prev->soffset);
  if (bead.next.isValid() ) {
    Bead *next = GetBead(beadStore,bead.next);
    char nqv = *Getchar(qualityStore,next->soffset);
    next->prev = bead.boffset;
    if (nqv < qv ) qv = nqv;
    if ( qv == '0'  ) {
      qv = '0' + 5;
    }
  }
  AppendVA_char(sequenceStore,&base);
  AppendVA_char(qualityStore,&qv);
  AppendVA_Bead(beadStore,&bead);
  gaps_in_alignment++;
  return bead.boffset;
}

//external
beadIdx
PrependGapBead(beadIdx bid) {
  // The gap will appear immediately before bid
  Bead *next = GetBead(beadStore,bid);
  Bead bead;
  char base='-';
  char qv;

  assert(next->frag_index >= 0);

  bead.boffset.set(GetNumBeads(beadStore));
  bead.soffset.set(GetNumchars(sequenceStore));
  bead.foffset = next->foffset;
  bead.up = beadIdx();
  bead.down = beadIdx();
  bead.frag_index = next->frag_index;
  bead.column_index = -1;
  bead.next = bid;
  bead.prev = next->prev;
  next->prev = bead.boffset;
  qv = *Getchar(qualityStore,next->soffset);
  if (bead.prev.isValid() ) {
    Bead *prev = GetBead(beadStore,bead.prev);
    char nqv = *Getchar(qualityStore,prev->soffset);
    prev->next = bead.boffset;
    if (nqv < qv ) qv = nqv;
    if ( qv == '0'  ) {
      qv = '0' + 5;
    }
  }
  AppendVA_char(sequenceStore,&base);
  AppendVA_char(qualityStore,&qv);
  AppendVA_Bead(beadStore,&bead);
  gaps_in_alignment++;
  return bead.boffset;
}


////////////////////////////////////////
//
//  Columns
//
//external
Column *
CreateColumn(beadIdx bid) {
  Column column;
  Bead call;
  Bead *head;

  // create a new column, seeded with the bead bid

  column.lid = GetNumColumns(columnStore);
  column.prev = -1;
  column.next = -1;
  column.call.set(GetNumBeads(beadStore));
  column.ma_index = -1;
  ResetBaseCount(&column.base_count);
  call.boffset = column.call;
  call.foffset = 0;
  call.soffset.set(GetNumchars(sequenceStore));
  call.down = bid;
  call.up = beadIdx();
  call.prev = beadIdx();
  call.next = beadIdx();
  call.frag_index = -1;
  call.column_index = column.lid;
  AppendVA_Bead(beadStore,&call);
  AppendVA_char(sequenceStore,"n");
  AppendVA_char(qualityStore,"0");
  head = GetBead(beadStore,bid);
  head->up = call.boffset;
  head->column_index = column.lid;
  IncBaseCount(&column.base_count,*Getchar(sequenceStore,head->soffset));
  AppendVA_Column(columnStore, &column);
#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "CreateColumn()-- Added consensus call bead="F_BEADIDX" to column="F_BEADIDX" for existing bead="F_BEADIDX"\n", call.boffset, column.lid, head->boffset);
#endif
  return GetColumn(columnStore, column.lid);
}


//external
void
AddColumnToMANode(int32 ma, Column column) {
  MANode *manode = GetMANode(manodeStore,ma);

  Appendint32(manode->columns,&column.lid);

  if (column.next == -1 )
    manode->last = column.lid;

  if (column.prev == -1 )
    manode->first = column.lid;
}



//external
int32
ColumnAppend(int32 cid, beadIdx bid) {
  // bid is the offset of the Bead seeding the column

  ColumnBeadIterator ci;
  beadIdx nid;

  Bead *bead = GetBead(beadStore,bid);
  assert(bead != NULL);

  Column *column = CreateColumn(bid);
  assert(column != NULL);

#ifdef DEBUG_ABACUS_ALIGN
  fprintf(stderr, "ColumnAppend()-- adding column "F_BEADIDX" for bid="F_BEADIDX" after column cid=%d\n", column->lid, bid, cid);
#endif

  Bead   *call     = GetBead(beadStore,column->call);
  Column *prev     = GetColumn(columnStore,cid);
  Bead   *prevcall = GetBead(beadStore,prev->call);

  column->next = prev->next;
  column->prev = cid;
  call->next = prevcall->next;
  call->prev = prevcall->boffset;
  prev->next = column->lid;
  prevcall->next = call->boffset;

  if (column->next != -1)
    GetColumn(columnStore,column->next)->prev = column->lid;

  if (call->next.isValid())
    GetBead(beadStore,call->next)->prev = call->boffset;

  CreateColumnBeadIterator(cid, &ci);

  while ( (nid = NextColumnBead(&ci)).isValid() ) {
    bead = GetBead(beadStore,nid);
    if ((bead->next.isValid()) &&
        (bead->next != bid))
      AlignBeadToColumn(column->lid, AppendGapBead(nid), "ColumnAppend()");
  }
  column->ma_id =  prev->ma_id;
  column->ma_index =  prev->ma_index + 1;
  AddColumnToMANode(column->ma_id,*column);
  return column->lid;
}


//external -- unused, but looks handy
#if 0
void
ShowColumn(beadIdx cid) {
  Column *column = GetColumn(columnStore,cid);
  Bead *call;
  Bead *bead;
  FragType type;
  UnitigType utype;
  ColumnBeadIterator ci;
  beadIdx bid;

  CreateColumnBeadIterator(cid,&ci);

  call = GetBead(beadStore,column->call);
  fprintf(stderr,"\nstore_index: %-20d ( prev: "F_BEADIDX" next: "F_BEADIDX")\n",column->lid,column->prev,column->next);
  fprintf(stderr,"ma_index:    %-20d\n",column->ma_index);
  fprintf(stderr,"------------------\n");
  fprintf(stderr,"composition:\n");
  while ( (bid = NextColumnBead(&ci)).isValid() ) {
    bead = GetBead(beadStore,bid);
    type = GetFragment(fragmentStore,bead->frag_index)->type;
    utype = GetFragment(fragmentStore,bead->frag_index)->utype;
    fprintf(stderr,"             %c /%c (%10d) <-- %d iid:%d cid:%d UDLR:%d %d %d %d type:%c utype:%c\n",
            *Getchar(sequenceStore,bead->soffset),
            *Getchar(qualityStore,bead->soffset),
            bid,
            bead->frag_index,
            GetFragment(fragmentStore,bead->frag_index)->iid,
            bead->column_index,
            bead->up,
            bead->down,
            bead->prev,
            bead->next,
            type,
            utype);

    assert(bead->column_index == cid);
  }
  fprintf(stderr,"------------------\n");
  fprintf(stderr,"call:        %c /%c\n",toupper(*Getchar(sequenceStore,call->soffset)),*Getchar(qualityStore,call->soffset));
}
#endif


////////////////////////////////////////
//
//  Data Management
//
void
ResetStores(int32 num_bases, int32 num_frags, int32 num_columns) {

  if (fragmentStore == NULL) {
    InitializeAlphTable();

    fragmentStore      = CreateVA_Fragment(num_frags);
    fragment_indices   = CreateVA_int32(num_frags);
    abacus_indices     = CreateVA_int32(50000);
    fragment_positions = CreateVA_CNS_AlignedContigElement(2 * num_frags);
    sequenceStore      = CreateVA_char(num_bases);
    qualityStore       = CreateVA_char(num_bases);
    columnStore        = CreateVA_Column(num_columns);
    beadStore          = CreateVA_Bead(num_bases + num_columns);
    manodeStore        = CreateVA_MANode(1);
  }

  ResetVA_Fragment(fragmentStore);
  MakeRoom_VA(fragmentStore, num_frags);

  ResetVA_int32(fragment_indices);
  MakeRoom_VA(fragment_indices, num_frags);

  ResetVA_int32(abacus_indices);

  ResetVA_CNS_AlignedContigElement(fragment_positions);
  MakeRoom_VA(fragment_positions, 2 * num_frags);

  ResetVA_char(sequenceStore);
  MakeRoom_VA(sequenceStore, num_bases);

  ResetVA_char(qualityStore);
  MakeRoom_VA(qualityStore, num_bases);

  ResetVA_Column(columnStore);
  MakeRoom_VA(columnStore, num_columns);

  ResetVA_Bead(beadStore);
  MakeRoom_VA(beadStore, num_bases + num_columns);

  ResetVA_MANode(manodeStore);

  gaps_in_alignment = 0;
}


//  Used in AppendFragToLocalStore
static
int
SetUngappedFragmentPositions(FragType type,int32 n_frags, MultiAlignT *uma) {

  int32 num_frags   = GetNumIntMultiPoss(uma->f_list);
  int32 num_unitigs = GetNumIntUnitigPoss(uma->u_list);

  HashTable_AS *unitigFrags = CreateScalarHashTable_AS();

  int32 num_columns   = GetMultiAlignLength(uma);
  int32 ungapped_pos  = 0;

  int32 *gapped_positions = (int32 *)safe_malloc(sizeof(int32) * (num_columns + 1));
  char  *consensus        = Getchar(uma->consensus,0);

  for (int32 i=0; i<num_columns+1; i++) {
    gapped_positions[i] = ungapped_pos;

    if (consensus[i] != '-')
      ungapped_pos++;
  }

  //  Remember the first fragment we add.
  int32 first_frag = GetNumCNS_AlignedContigElements(fragment_positions);

  for (int32 ifrag=0; ifrag<num_frags; ifrag++) {
    CNS_AlignedContigElement epos;
    IntMultiPos *frag = GetIntMultiPos(uma->f_list, ifrag);

    if (ExistsInHashTable_AS(unitigFrags, frag->ident, 0)) {
      fprintf(stderr,"SetUngappedFragmentPositions()-- ident %d already in hashtable\n", frag->ident);
      assert(0);
    }
    if (HASH_SUCCESS != InsertInHashTable_AS(unitigFrags, frag->ident, 0, 1, 0)) {
      fprintf(stderr,"SetUngappedFragmentPositions()-- Failure to insert ident %d in hashtable\n", frag->ident);
      assert(0);
    }

    assert(frag->position.bgn >= 0);
    assert(frag->position.bgn < num_columns + 1);
    assert(frag->position.end >= 0);
    assert(frag->position.end < num_columns + 1);

    epos.frg_or_utg                  = CNS_ELEMENT_IS_FRAGMENT;
    epos.idx.fragment.frgIdent       = frag->ident;
    epos.idx.fragment.frgType        = frag->type;
    epos.idx.fragment.frgContained   = frag->contained;
    epos.idx.fragment.frgInUnitig    = (type == AS_CONTIG) ? -1 : uma->maID;
    epos.position.bgn                = gapped_positions[frag->position.bgn];
    epos.position.end                = gapped_positions[frag->position.end];

    //fprintf(stderr, "SetUngappedFragmentPositions()-- FRG id=%d type=%c pos=%d,%d (orig pos=%d,%d)\n",
    //        frag->ident, frag->type, epos.position.bgn, epos.position.end, frag->position.bgn, frag->position.end);

    //  Adjust the ungapped position if we fall within a gap
    //
    if (epos.position.bgn == epos.position.end) {
      fprintf(stderr,"SetUngappedFragmentPositions()-- Encountered bgn==end=="F_S32" in ungapped coords within SetUngappedFragmentPositions for "F_CID "(gapped coords "F_S32","F_S32")\n",
              epos.position.bgn,frag->ident,frag->position.bgn,frag->position.end);
      assert(frag->position.bgn != frag->position.end);

      if (frag->position.bgn < frag->position.end) {
        if (epos.position.bgn > 0)
          epos.position.bgn--;
        else
          epos.position.end++;
      } else {
        if (epos.position.end > 0)
          epos.position.end--;
        else
          epos.position.bgn++;
      }
      fprintf(stderr,"SetUngappedFragmentPositions()--   Reset to "F_S32","F_S32"\n",
              epos.position.bgn,
              epos.position.end);
    }

    AppendVA_CNS_AlignedContigElement(fragment_positions, &epos);
  }


  for (int32 ifrag=0; ifrag < num_unitigs; ifrag++){
    CNS_AlignedContigElement epos;
    IntUnitigPos *unitig = GetIntUnitigPos(uma->u_list, ifrag);

    epos.frg_or_utg           = CNS_ELEMENT_IS_UNITIG;
    epos.idx.unitig.utgIdent  = unitig->ident;
    epos.idx.unitig.utgType   = unitig->type;
    epos.position.bgn         = gapped_positions[unitig->position.bgn];
    epos.position.end         = gapped_positions[unitig->position.end];

    //fprintf(stderr, "SetUngappedFragmentPositions()-- UTG id=%d type=%c pos=%d,%d (orig pos=%d,%d)\n",
    //        unitig->ident, unitig->type, epos.position.bgn, epos.position.end, unitig->position.bgn, unitig->position.end);

    AppendVA_CNS_AlignedContigElement(fragment_positions,&epos);
  }

  //  This is used only by ReplaceEndUnitigInContig().  Mark fragments in the "anchoring" contig
  //  that belong to this unitig.
  //
  if (type != AS_CONTIG) {
    Fragment *anchor = GetFragment(fragmentStore,0);

    if ((anchor != NULL) &&
        (anchor->type == AS_CONTIG)) {
      CNS_AlignedContigElement *af = GetCNS_AlignedContigElement(fragment_positions, anchor->components);

      for (int32 ifrag=0; ifrag < anchor->n_components; ifrag++, af++) {
        if ((af->frg_or_utg == CNS_ELEMENT_IS_FRAGMENT) &&
            (ExistsInHashTable_AS(unitigFrags, af->idx.fragment.frgIdent, 0)))
          af->idx.fragment.frgInUnitig = uma->maID;
      }
    }
  }

  DeleteHashTable_AS(unitigFrags);
  safe_free(gapped_positions);

  return first_frag;
}



int32
AppendFragToLocalStore(FragType          type,
                       int32               iid,
                       int32               complement,
                       int32               contained,
                       UnitigType        utype) {

  char seqbuffer[AS_READ_MAX_NORMAL_LEN+1];
  char qltbuffer[AS_READ_MAX_NORMAL_LEN+1];
  char *sequence = NULL,*quality = NULL;
  static VA_TYPE(char) *ungappedSequence = NULL;
  static VA_TYPE(char) *ungappedQuality  = NULL;
  Fragment fragment;
  uint32 clr_bgn, clr_end;
  static gkFragment fsread;  //  static for performance only
  MultiAlignT *uma = NULL;

  if (ungappedSequence == NULL) {
    ungappedSequence = CreateVA_char(0);
    ungappedQuality  = CreateVA_char(0);
  }

  switch (type) {
    case AS_READ:
    case AS_EXTR:
    case AS_TRNR:
      gkpStore->gkStore_getFragment(iid,&fsread,GKFRAGMENT_QLT);

      fsread.gkFragment_getClearRegion(clr_bgn, clr_end);

      strcpy(seqbuffer, fsread.gkFragment_getSequence());
      strcpy(qltbuffer, fsread.gkFragment_getQuality());

#ifdef PRINTUIDS
      fragment.uid = gkFragment_getReadUID(&fsread);
#endif
      fragment.type = AS_READ;
      fragment.source = NULL;
      seqbuffer[clr_end] = '\0';
      qltbuffer[clr_end] = '\0';
      sequence = &seqbuffer[clr_bgn];
      quality = &qltbuffer[clr_bgn];
      fragment.length = (int32) (clr_end - clr_bgn);
      fragment.n_components = 0;  // no component frags or unitigs
      fragment.components = -1;

      //fprintf(stderr, "AppendFragToLocalStore()-- FRG %d len=%d clr=%d,%d\n", iid, fragment.length, clr_bgn, clr_end);
      break;
    case AS_UNITIG:
    case AS_CONTIG:
      if (tigStore)
        uma = tigStore->loadMultiAlign(iid, type == AS_UNITIG);
      if (uma == NULL)
        fprintf(stderr,"Lookup failure in CNS: MultiAlign for unitig %d could not be found.\n",iid);
      assert(uma != NULL);

      //  Contigs used to be added gapped, unitigs as ungapped.
      //  This caused no end of trouble in MergeMultiAligns and
      //  ReplaceEndUnitigInContig.

      ResetVA_char(ungappedSequence);
      ResetVA_char(ungappedQuality);

      GetMultiAlignUngappedConsensus(uma, ungappedSequence, ungappedQuality);

      sequence = Getchar(ungappedSequence,0);
      quality = Getchar(ungappedQuality,0);

      fragment.length = GetMultiAlignUngappedLength(uma);

      fragment.utype = (type == AS_UNITIG) ? utype : AS_OTHER_UNITIG;

      fragment.n_components = GetNumIntMultiPoss(uma->f_list) + GetNumIntUnitigPoss(uma->u_list);
      fragment.components   = SetUngappedFragmentPositions(type, fragment.n_components, uma);

      //fprintf(stderr, "AppendFragToLocalStore()-- TIG %d len=%d\n", iid, fragment.length);
      break;
    default:
      fprintf(stderr, "AppendFragToLocalStore invalid FragType");
      assert(0);
  }

  if (complement)
    reverseComplement(sequence, quality, strlen(sequence));

  fragment.lid = GetNumFragments(fragmentStore);
  fragment.iid = iid;
  fragment.type = type;
  fragment.complement = complement;
  fragment.container_iid = contained;
  fragment.is_contained = (contained > 0) ? 1 : 0;
  fragment.deleted = 0;
  fragment.manode = -1;

  fragment.sequence.set(GetNumchars(sequenceStore));
  fragment.firstbead.set(GetNumBeads(beadStore));

  AppendRangechar(sequenceStore, fragment.length + 1, sequence);
  AppendRangechar(qualityStore,  fragment.length + 1, quality);

  {
    Bead bead;

    beadIdx boffset = fragment.firstbead;
    seqIdx  soffset = fragment.sequence;
    int32  foffset;

    bead.boffset      = beadIdx();
    bead.soffset      = seqIdx();
    bead.foffset      = -1;
    bead.prev         = beadIdx();
    bead.next         = beadIdx();
    bead.up           = beadIdx();
    bead.down         = beadIdx();
    bead.frag_index   = fragment.lid;
    bead.column_index = -1;

    for (foffset = 0; foffset < fragment.length; foffset++ ) {
      bead.foffset = foffset;
      bead.boffset.set(boffset.get() + foffset);
      bead.soffset.set(soffset.get() + foffset);

      bead.next.set(bead.boffset.get() + 1);
      bead.prev.set(bead.boffset.get() - 1);

      if (foffset == fragment.length - 1)
        bead.next = beadIdx();

      if (foffset == 0)
        bead.prev = beadIdx();

      SetBead(beadStore, bead.boffset.get(), &bead);
    }
  }

  AppendVA_Fragment(fragmentStore,&fragment);

  return fragment.lid;
}







////////////////////////////////////////
//
//  Variation Functions
//
//external
void
AllocateDistMatrix(VarRegion  *vreg, int32 init) {

  vreg->dist_matrix = (int32 **)safe_calloc(vreg->nr, sizeof(int32 *));
  for (int32 j=0; j<vreg->nr; j++) {
      vreg->dist_matrix[j] = (int32 *)safe_calloc(vreg->nr, sizeof(int));
      for (int32 k=0; k<vreg->nr; k++)
        vreg->dist_matrix[j][k] = init;
    }
}



#ifdef DEBUG_VAR_RECORDS
//external
void
OutputDistMatrix(FILE *fout, VarRegion  *vreg) {

  fprintf(fout, "Distance matrix=\n");
  for (int32 j=0; j<vreg->nr; j++) {
      for (int32 k=0; k<vreg->nr; k++)
        fprintf(fout, " %d", vreg->dist_matrix[j][k]);
      fprintf(fout, "\n");
    }
}
#endif



//  Used in PopulateDistMatrix
static
int
GetDistanceBetweenReads(char *read1, char *read2, int32 len) {
  int32 i, j, k, uglen1=0, uglen2=0, uglen;
  int32 dist, gapped_dist = 0, ungapped_dist = 0;
  char *ugread1 = (char*)safe_malloc(len*sizeof(char));
  char *ugread2 = (char*)safe_malloc(len*sizeof(char));

  // Compute gapped distance
  for (k=0; k<len; k++) {
    if (read1[k] != read2[k])
      gapped_dist++;

    if (read1[k] != '-') {
        ugread1[uglen1] = read1[k];
        uglen1++;
      }
    if (read2[k] != '-') {
        ugread2[uglen2] = read2[k];
        uglen2++;
      }
  }

  uglen = (uglen1<uglen2) ? uglen2:uglen1;
  for (k=0; k<uglen; k++) {
      // Compute ungapped distance
      if (k<uglen1 && k<uglen2 && ugread1[k] != ugread2[k])
        ungapped_dist++;
      else if (k <uglen1 && k>=uglen2)
        ungapped_dist++;
      else if (k>=uglen1 && k <uglen2)
        ungapped_dist++;
    }
  dist = (gapped_dist < ungapped_dist) ? gapped_dist : ungapped_dist;
  safe_free(ugread1);
  safe_free(ugread2);
  return dist;
}


//external
void
PopulateDistMatrix(Read *reads, int32 len, VarRegion  *vreg) {
  int32 i, j;

  // Update the matrix
  for (i=0; i<vreg->nr; i++) {
    for (j=i; j<vreg->nr; j++) {
      vreg->dist_matrix[i][j] = GetDistanceBetweenReads(reads[i].bases, reads[j].bases, len);
      vreg->dist_matrix[j][i] = vreg->dist_matrix[i][j];
    }
  }
}



//external
void
OutputReads(FILE *fout, Read *reads, int32 nr, int32 width) {
  int32 i, j;
  fprintf(fout, "\nReads =\n");

  for (i=0; i<nr; i++) {
    fprintf(fout, "%d   ", reads[i].allele_id);
    for (j=0; j<width; j++)
      fprintf(fout, "%c", reads[i].bases[j]);
    fprintf(fout, "\n");
  }
  fprintf(fout, "\n\n");
}

//external
void
OutputAlleles(FILE *fout, VarRegion *vreg) {
  int32 i, j;
  fprintf(fout,   "Outputting alleles:\n");
  fprintf(fout,   "nr= %d na= %d nca= %d\n", vreg->nr, vreg->na, vreg->nca);
  fprintf(fout,   "Num_reads= ");
  for (i=0; i<vreg->na; i++)
      fprintf(fout,   "%d ", vreg->alleles[i].num_reads);
  fprintf(fout,   "\n");
  fprintf(fout,   "Weights= ");
  for (i=0; i<vreg->na; i++)
      fprintf(fout,   "%d ", vreg->alleles[i].weight);
  fprintf(fout,   "\n");
  fprintf(fout,   "Reads= \n");
  for (i=0; i<vreg->na; i++) {
      fprintf(fout,   "   Allele order= %d, id= %d:\n", i, vreg->alleles[i].id);
      for (j=0; j<vreg->alleles[i].num_reads; j++) {
          int32 k, read_id = vreg->alleles[i].read_ids[j];
          int32 len = vreg->end-vreg->beg;
          fprintf(fout,   "    %d   ", read_id);
          for (k=0; k<len; k++)
            fprintf(fout,   "%c", vreg->reads[read_id].bases[k]);
          fprintf(fout,   "   %d\n", vreg->alleles[i].read_iids[j]);
        }
    }
}


// Allocate memrory for reads
//external
void
AllocateMemoryForReads(Read **reads, int32 nr, int32 len,
                       int32 default_qv) {
  int32 i, j;

  assert(nr > 0);

  *reads = (Read *)safe_malloc(nr*sizeof(Read));
  for (i=0; i<nr; i++) {
      (*reads)[i].allele_id = -1;
      (*reads)[i].ave_qv = 0.;
      (*reads)[i].bases = (char *)safe_malloc(len*sizeof(char));
      (*reads)[i].qvs   = (int32  *)safe_malloc(len*sizeof(int32 ));
      for(j=0; j<len; j++) {
          (*reads)[i].bases[j] = '-';
          (*reads)[i].qvs[j] = default_qv;
        }
    }
}

// Allocate memrory for alleles
//external
void
AllocateMemoryForAlleles(Allele **alleles, int32 nr, int32 *na) {
  int32 j;

  assert(nr > 0);

  *na = 0;
  *alleles = (Allele *)safe_calloc(nr, sizeof(Allele));
  for (j=0; j<nr; j++) {
      (*alleles)[j].id = -1;
      (*alleles)[j].weight = 0;
      (*alleles)[j].read_ids = (int32 *)safe_calloc(nr, sizeof(int));
      (*alleles)[j].read_iids = (int32 *)safe_calloc(nr, sizeof(int32));
    }
}


// Reverse sort confirmed alleles by ungapped length
//external
 void
SortAllelesByLength(Allele *alleles, int32 num_alleles, Read *reads) {
  int32 i, j, best_id;
  Allele temp;

  for (i=0; i<num_alleles; i++) {
      int32 best_uglen = alleles[i].uglen;
      best_id = -1;
      for (j=i+1; j<num_alleles; j++) {
          if (best_uglen  < alleles[j].uglen ) {
              best_uglen  = alleles[j].uglen ;
              best_id = j;
            }
        }
      if (best_id >= 0) {
          temp       = alleles[i];
          alleles[i] = alleles[best_id];
          alleles[best_id] = temp;
        }
    }
  // Update allele_id of reads
  for (i=0; i<num_alleles; i++) {
      for (j=0; j<alleles[i].num_reads; j++) {
          int32 read_id = alleles[i].read_ids[j];
          reads[read_id].allele_id = i;
        }
    }
}


// Reverse sort by weight
//external
 void
SortAllelesByWeight(Allele *alleles, int32 num_alleles, Read *reads) {
  int32 i, j, best_id;
  Allele temp;

  for (i=0; i<num_alleles; i++) {
      int32 best_weight = alleles[i].weight;
      best_id = -1;
      for (j=i+1; j<num_alleles; j++) {
          if (best_weight < alleles[j].weight) {
              best_weight = alleles[j].weight;
              best_id = j;
            }
        }
      if (best_id >= 0) {
          temp       = alleles[i];
          alleles[i] = alleles[best_id];
          alleles[best_id] = temp;
        }
    }
  // Update allele_id of reads
  for (i=0; i<num_alleles; i++) {
      for (j=0; j<alleles[i].num_reads; j++) {
          int32 read_id = alleles[i].read_ids[j];
          reads[read_id].allele_id = i;
        }
    }
}

// Sort confirmed alleles according to their mapping
// between two "phased" VAR records
//external
void
SortAllelesByMapping(Allele *alleles, int32 nca, Read *reads, int32 *allele_map) {
  int32 i, j, k;
  Allele temp;

  for (i=0; i<nca; i++) {
      // j is id of the allele that should be at i-th place
      for (j=0; j<nca; j++)
        if (allele_map[j] == i) break;

      for (k=i; k<nca; k++) {
          if (alleles[k].id == j) {
              temp       = alleles[i];
              alleles[i] = alleles[k];
              alleles[k] = temp;
              break;
            }
        }
    }

  // Update allele_ids
  for (i=0; i<nca; i++) {
      alleles[i].id = i;
      for (j=0; j<alleles[i].num_reads; j++) {
          int32 read_id = alleles[i].read_ids[j];
          reads[read_id].allele_id = i;
        }
    }
}



/*******************************************************************************
 * Function: ClusterReads
 * Purpose:  detect allele and split reads between the alleles
 *******************************************************************************
 */
//external
void
ClusterReads(Read *reads, int32 nr, Allele *alleles, int32 *na, int32 *nca, int32 **dist_matrix) {
  int32 aid, anr, row, col;

  *na = 0;

  // Iniytialize alleles

  // Process zero elements first
  for (row=0; row<nr; row++) {
      for (col=row+1; col<nr; col++) {
          if (dist_matrix[row][col]!=0)
            continue;

          if (reads[row].allele_id < 0 &&
              reads[col].allele_id < 0) {
              // New allele
              reads[row].allele_id = *na;
              reads[col].allele_id = *na;
              alleles[*na].weight =
                ROUND(reads[row].ave_qv) + ROUND(reads[col].ave_qv);
              alleles[*na].uglen  = reads[row].uglen;
              alleles[*na].read_ids[0] = row;
              alleles[*na].read_ids[1] = col;
              alleles[*na].read_iids[0] = reads[row].iid;
              alleles[*na].read_iids[1] = reads[col].iid;
              alleles[*na].num_reads = 2;
              alleles[*na].id = *na;
              (*na)++;
            } else if (reads[row].allele_id < 0 &&
                   reads[col].allele_id >=0) {
              // Already existing allele
              aid = reads[col].allele_id;
              reads[row].allele_id = aid;
              alleles[aid].weight += ROUND(reads[row].ave_qv);
              anr = alleles[aid].num_reads;
              alleles[aid].read_ids[anr] = row;
              alleles[aid].read_iids[anr] = reads[row].iid;
              alleles[aid].num_reads++;
            } else if (reads[row].allele_id >=0 &&
                   reads[col].allele_id < 0) {
              // Already existing allele
              aid = reads[row].allele_id;
              reads[col].allele_id = aid;
              alleles[aid].weight += ROUND(reads[col].ave_qv);
              anr = alleles[aid].num_reads;
              alleles[aid].read_ids[anr] = col;
              alleles[aid].read_iids[anr] = reads[col].iid;
              alleles[aid].num_reads++;
            }
        }
    }

  *nca = *na;

  //Now process the remaining reads; assign each to its "own" allele
  for (row=0; row<nr; row++) {
      if (reads[row].allele_id < 0) {
          // New allele
          reads[row].allele_id      = *na;
          alleles[*na].weight       = ROUND(reads[row].ave_qv);
          alleles[*na].uglen        = reads[row].uglen;
          alleles[*na].read_ids[0]  = row;
          alleles[*na].read_iids[0] = reads[row].iid;
          alleles[*na].num_reads    = 1;
          alleles[*na].id           = *na;
          (*na)++;
        }
    }
}
