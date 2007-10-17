
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
/*********************************************************************
   Module:       MultiAlignment_CNS.c
   Description:  multialignment and associated functions
   Assumptions:  
*********************************************************************/

static char CM_ID[] = "$Id: MultiAlignment_CNS.c,v 1.170 2007-10-17 22:34:58 gdenisov Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MicroHetREZ.h"


/* Controls for the DP_Compare and Realignment schemes */
//#define USE_AFFINE_OVERLAP
//#define USE_LOCAL_OVERLAP
#ifdef USE_AFFINE_OVERLAP
//#define COMPARE_FUNC Affine_Overlap_AS_forCNS
#define CMPFNC "Affine_Overlap_AS_forCNS"
#include "AS_ALN_forcns.h"
#define HANDLE_OVERLAP_INITIAL_GAPS
#elif defined(USE_LOCAL_OVERLAP)
// #define COMPARE_FUNC Local_Overlap_AS_forCNS
#define CMPFNC "Local_Overlap_AS_forCNS"
#include "AS_ALN_forcns.h"
#define HANDLE_OVERLAP_INITIAL_GAPS
#else
//#define COMPARE_FUNC DP_Compare
#define CMPFNC "DP_Compare(nonaffine)"
#endif

#define ALT_QV_THRESH                      30
#define IDENT_NAMESPACE                     1
#define DONT_SHOW_OLAP                      0
#define MIN_AVE_QV_FOR_VARIATION           21
#define MIN_SUM_QVS_FOR_VARIATION          60
#define QV_FOR_MULTI_GAP                   14
#define SHOW_OLAP                           1
#undef  ALIGN_TO_CONSENSUS
#define PRINTUIDS

#define CNS_DP_THRESH                       1e-6
#define CNS_DP_MINLEN                      30

#define CNS_NEG_AHANG_CUTOFF               -5
#define CNS_MAX_ALIGN_SLIP                 20
#define INITIAL_NR                        100
#define MAX_ALLOWED_MA_DEPTH               40
#define MAX_EXTEND_LENGTH                2048
#define MAX_WINDOW_FOR_ABACUS_REFINE      100
#define SHOW_ABACUS                         0
#define SHOW_ALLELES                        0
#define SHOW_READS                          0
#define SHOW_CONFIRMED_READS                0
#define STABWIDTH                           6
#define DEBUG_ABACUS                        0
#define DEBUG_CONSENSUS_CALL                0
#define DEBUG_VAR_RECORDS                   0
#define OUTPUT_DISTANCE_MATRIX              0

// Parameters used by Abacus processing code
#define MSTRING_SIZE                        3
#define MAX_SIZE_OF_ADJUSTED_REGION         5



//====================================================================
// Persistent store of the fragment data (produced upstream)
//
GateKeeperStore       *gkpStore    = NULL;
tSequenceDB           *sequenceDB  = NULL;
HashTable_AS          *fragmentMap = NULL;
MultiAlignStoreT      *unitigStore = NULL;


//====================================================================
// Stores for the sequence/quality/alignment information
// (reset after each multialignment)
//
VA_TYPE(char) *sequenceStore = NULL;
VA_TYPE(char) *qualityStore  = NULL;
VA_TYPE(Bead) *beadStore     = NULL;

//====================================================================
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

ALIGN_TYPE ALIGNMENT_CONTEXT = AS_NOALIGNTYPE;

int USE_SDB = 0;


//====================================================================
// Convenience arrays for misc. fragment information
// (All are reset after each multialignment)
//
VA_DEF(PtrT)

VA_TYPE(int32) *fragment_indices  = NULL;
VA_TYPE(int32) *abacus_indices    = NULL;
VA_TYPE(PtrT) *fragment_positions = NULL;

int64 gaps_in_alignment = 0;

int allow_forced_frags   = 0;
int allow_neg_hang       = 0;


extern int MaxBegGap;       // [ init value is 200; this could be set to the amount you extend the clear
                            // range of seq b, plus 10 for good measure]
extern int MaxEndGap;       // [ init value is 200; this could be set to the amount you extend the
                            // clear range of seq a, plus 10 for good measure]

VA_DEF(int16);

// Variables used to phase VAR records
static int32 vreg_id  = 0; // id of a VAR record

// Variables used to compute general statistics

int NumColumnsInUnitigs = 0;
int NumRunsOfGapsInUnitigReads = 0;
int NumGapsInUnitigs = 0;
int NumColumnsInContigs = 0;
int NumRunsOfGapsInContigReads = 0;
int NumGapsInContigs = 0;
int NumAAMismatches = 0; // mismatches b/w consensi of two different alleles
int NumVARRecords = 0;
int NumVARStringsWithFlankingGaps = 0;
int NumUnitigRetrySuccess = 0;
int contig_id = 0;

//*********************************************************************************
//  Tables to facilitate SNP Basecalling
//*********************************************************************************

double EPROB[CNS_MAX_QV-CNS_MIN_QV+1] = {0};  // prob of error for each quality value
double PROB[CNS_MAX_QV-CNS_MIN_QV+1]  = {0};  // prob of correct call for each quality value (1-eprob)
int    RINDEX[128] = {0};

#ifdef AS_ENABLE_SOURCE
char SRCBUFFER[2048] = {0};
#endif

// Utility variable to control width of "pages" of PrintAlignment output
int ALNPAGEWIDTH=100;

//  Define this to dump multifasta to stderr of the unitigs
//  we try to align in MultiAlignContig().  Was useful for
//  debugging bad layout.
//
int DUMP_UNITIGS_IN_MULTIALIGNCONTIG = 0;

// Be noisy when doing multi alignments - this used to be a #ifdef,
// which made it difficult to switch on in the middle of a debug.
//
int VERBOSE_MULTIALIGN_OUTPUT = 0;

//  If non-zero, we'll force-abut unitigs that don't align together.
//  Typically, these are caused by microscopic overlaps between
//  unitigs -- certainly less than 20bp long.
//
int FORCE_UNITIG_ABUT = 0;

//  Run-time selectable clear range; for unitigs, we use
//  AS_READ_CLEAR_OBT (the same as unitigger); for contigs, we use
//  AS_READ_CLEAR_LATEST.
//
int clear_range_to_use = AS_READ_CLEAR_LATEST;


int isRead(FragType type){
  switch(type){
    case  AS_READ   :
    case  AS_EXTR   :
    case  AS_TRNR   :
      return 1;
    default:
      return 0;
  }
}

int isChunk(FragType type){
  switch(type){
    case AS_UNITIG  :
    case AS_CONTIG  :
      return 1;
    default:
      return 0;
  }
}

int InitializeAlphTable(void) {
  int count=sizeof(RINDEX)/sizeof(int);
  int i;
  if (!RALPH_INIT) {
    for(i=0;i<count;i++) { 
      RINDEX[i] = 31;
    }
    count=sizeof(RALPHABET)/sizeof(char);
    for(i=0;i<count;i++) { 
      RINDEX[(int)RALPHABET[i]] = i;
    }

    { int qv=CNS_MIN_QV;
      for (i=0;i<CNS_MAX_QV-CNS_MIN_QV+1;i++) {
        EPROB[i]= pow(10,-qv/10.);
        PROB[i] = (1.0 - EPROB[i]);
        qv++;
      }
    }
    // show bits of mask for testing
    //   for (i=0;i<5;i++) {
    //     int j;
    //     for (j=0;j<30;j++) {
    //       fprintf(stderr,"%1d",((AMASK[i]>>j) & 1)?1:0);
    //     }
    //       fprintf(stderr,"\n");
    //   }
  }
  return 1;
}

//*********************************************************************************
// Basic conversions, etc. for sequence manipulation
//*********************************************************************************

int RBaseToInt(char c) {  
  // Translate characters representing base reads to array index
  return RINDEX[(int)c];
}

static char RIntToBase(int i) { return RALPHABET[i]; }

static char RBaseComplement(char c) {  
  /* Translate characters representing base reads to array index */
  return RALPHABETC[RINDEX[(int)c]];
}

int BaseToInt(char c) {  
  // Translate characters representing base reads to array index
  return RINDEX[(int)c];
}

static char IntToBase(int i) { return RALPHABET[i]; }
//static char IntToBase(int i) { return ALPHABET[i]; }

char BaseComplement(char c) {  
  return RALPHABETC[RINDEX[(int)c]];
  /* Translate characters representing base reads to array index */
}

void SequenceComplement(char *sequence, char *quality) {
  char *s, *t;
  int len,c;
  InitializeAlphTable();
  len = strlen(sequence);
  s = sequence;
  t = sequence + (len-1);
  while (s < t) { 
    c = *s;
    *s++ = BaseComplement(*t);
    *t-- = BaseComplement(c);
  }
  if (s == t) {
    *s = BaseComplement(*s);
  }
  if (quality != NULL) {
    s = quality;
    t = quality + (len-1);
    while (s < t) { 
      c = *s;
      *s++ = *t;
      *t-- = c;
    }
  }
}

//*********************************************************************************
// Manipulation of the BaseCount struct which reflects a char profile
//*********************************************************************************

int IncBaseCount(BaseCount *b, char c) {
  int i= BaseToInt(c);
  if (c == 'N' || c == 'n' ) i=5;
  b->depth++;
  if( i<0 || i>5 ){
    fprintf(stderr, "IncBaseCount i out of range (possibly non ACGTN letter?)");
    assert(0);
  }
  return b->count[i]++;
}

int DecBaseCount(BaseCount *b, char c) {
  int i= BaseToInt(c);
  if (c == 'N' || c == 'n' ) i=5;
  b->depth--;
  if( i<0 || i>5 ){
    fprintf(stderr, "DecBaseCount i out of range");
    assert(0);
  }
  return b->count[i]--;
}

int GetBaseCount(BaseCount *b, char c) {
  int i= BaseToInt(c);
  if (c == 'N' || c == 'n' ) i=5;
  return b->count[i];
}

int GetColumnBaseCount(Column *b, char c) {
  return GetBaseCount(&b->base_count,c);
}

int GetDepth(Column *c) {
  return c->base_count.depth;
}

void ResetBaseCount(BaseCount *b) {
  memset(b,'\0',sizeof(BaseCount)); 
}

void ShowBaseCount(BaseCount *b) {
  int i;
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

void ShowBaseCountPlain(FILE *out,BaseCount *b) {
  int i;
  fprintf(out,"%d\t",b->depth);
  for (i=0;i<CNS_NALPHABET;i++) {
    fprintf(out,"%d\t",b->count[i]);
  }
}

char GetConfMM(BaseCount *b,int mask) {  // mask out the consensus base
  int i;
  for (i=0;i<CNS_NALPHABET-1;i++) {
    if ( i==mask ) {
      continue; 
    }
    if ( b->count[i] >= 2 ) {
      return toupper(ALPHABET[i]);
    }
  }
  return  toupper(ALPHABET[mask]); // return the consensus base if there's no confirmed mismatch
}

char GetMaxBaseCount(BaseCount *b,int start_index) {  // start at 1 to disallow gap
  int max_index = start_index,i;
  int tied = 0,tie_breaker,max_tie=0;
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

  
//*********************************************************************************
// MANode (MultiAlignmentNode) creation
//*********************************************************************************

MANode * CreateMANode(int32 iid){
  MANode ma;
  ma.lid = GetNumMANodes(manodeStore);
  ma.iid = iid;
  ma.first = -1;
  ma.last = -1;
  ma.columns = CreateVA_int32(GetAllocatedColumns(columnStore));
  AppendVA_MANode(manodeStore,&ma);
  return GetMANode(manodeStore,ma.lid);
}

void DeleteMANode(int32 iid){
  MANode *ma=GetMANode(manodeStore,iid);
  // Columns are in the columnStore, which is automatically refreshed
  DeleteVA_int32(ma->columns);
}

int32 GetMANodeLength(int32 mid) {
  MANode *ma = GetMANode(manodeStore,mid);
  if ((ma) == NULL) return -1;
  return GetNumint32s(ma->columns);
}

//*********************************************************************************
//  Reset all the local stores
//*********************************************************************************

void ResetStores(int32 num_frags, int32 num_columns) {

  if ( fragmentStore == NULL ) {
    fragmentStore = CreateVA_Fragment(num_frags);
  } else {
    ResetVA_Fragment(fragmentStore);
    MakeRoom_VA(fragmentStore,num_frags,FALSE);
  }
   
  if ( fragment_indices == NULL ) {
    fragment_indices = CreateVA_int32(num_frags);
    abacus_indices = CreateVA_int32(50000);
  } else {
    ResetVA_int32(fragment_indices);
    MakeRoom_VA(fragment_indices,num_frags,FALSE);
    ResetVA_int32(abacus_indices);
  }
   
  if ( fragment_positions == NULL ) {
    fragment_positions = CreateVA_CNS_AlignedContigElement(2*num_frags);
  } else {
    ResetVA_CNS_AlignedContigElement(fragment_positions);
    MakeRoom_VA(fragment_positions,2*num_frags,FALSE);
  }
   
  if ( sequenceStore == NULL ) {
    sequenceStore = CreateVA_char(2048*num_frags);
    qualityStore = CreateVA_char(2048*num_frags);
  } else {
    ResetVA_char(sequenceStore);
    MakeRoom_VA(sequenceStore,2048*num_frags,FALSE);
    ResetVA_char(qualityStore);
    MakeRoom_VA(qualityStore,2048*num_frags,FALSE);
  }
   
  if ( columnStore == NULL ) {
    columnStore = CreateVA_Column(num_columns);
  } else {
    ResetVA_Column(columnStore);
    MakeRoom_VA(columnStore,num_columns,FALSE);
  }
   
  if ( beadStore == NULL ) {
    beadStore = CreateVA_Bead(2048*num_frags+num_columns);
  } else {
    ResetVA_Bead(beadStore);
    MakeRoom_VA(beadStore,2048*num_frags+num_columns,FALSE);
  }
   
  if ( manodeStore == NULL ) {
    manodeStore = CreateVA_MANode(1);
  } else {
    ResetVA_MANode(manodeStore);
  }
  gaps_in_alignment=0;
}

//*********************************************************************************
// Iterator for Column bases (called "beads")
//*********************************************************************************
int CreateColumnBeadIterator(int32 cid,ColumnBeadIterator *bi) {
  Column *column = GetColumn(columnStore,cid);
  if ( column == NULL ) { return 0;}
  bi->column = *column;
  bi->bead = bi->column.call;
  return 1;
}

int32 NextColumnBead(ColumnBeadIterator *bi) {
  int32 nid;
  Bead *bead;
  if (bi->bead == -1 ) {
    return -1;
  }
  bead = GetBead(beadStore, bi->bead);
  nid = bead->down;
  bi->bead = nid;
  return nid;
}

//*********************************************************************************
// Iterator for Fragment bases (called "beads")
//*********************************************************************************
int NullifyFragmentBeadIterator(FragmentBeadIterator *bi) {
  bi->fragment = *GetFragment(fragmentStore,0);
  bi->bead = -2;
  return 1;
}

int IsNULLIterator(FragmentBeadIterator *bi) {
  return ( bi->bead == -2 );
}

int CreateFragmentBeadIterator(int32 fid,FragmentBeadIterator *bi) {
  Fragment *fragment = GetFragment(fragmentStore,fid);
  if ( fragment == NULL ) { return 0;}
  bi->fragment = *fragment;
  bi->bead = bi->fragment.beads;
  return 1;
}

int32 NextFragmentBead(FragmentBeadIterator *bi) {
  int32 nid;
  Bead *bead;
  if (bi->bead == -1 ) {
    return -1;
  }
  bead = GetBead(beadStore, bi->bead);
  nid = bead->boffset;
  bi->bead = bead->next;
  return nid;
}

//*********************************************************************************
// Iterator for Consensus bases (called "beads")
//*********************************************************************************
int CreateConsensusBeadIterator(int32 mid,ConsensusBeadIterator *bi) {
  Column *first = GetColumn(columnStore,(GetMANode(manodeStore,mid))->first);
  bi->manode_id = mid;
  bi->bead = first->call;
  return 1;
}

int32 NextConsensusBead(ConsensusBeadIterator *bi) {
  int32 nid;
  Bead *bead;
  if (bi->bead == -1 ) {
    return -1;
  }
  bead = GetBead(beadStore, bi->bead);
  nid = bead->boffset;
  bi->bead = bead->next;
  return nid;
}

//*********************************************************************************
// Iterator for columns 
//*********************************************************************************
int CreateColumnIterator(int32 cid,ColumnIterator *ci) {
  GetColumn(columnStore,cid);
  ci->cid = cid;
  return 1;
}

int32 NextColumn(ColumnIterator *ci) {
  int32 nid;
  Column *column;
  if (ci->cid == -1 ) {
    return -1;
  }
  column = GetColumn(columnStore, ci->cid);
  nid = ci->cid;
  ci->cid = column->next;
  return nid;
}

//*********************************************************************************
// Insert a "gap bead" in a bead sequence (return the new bead's id)
//    int32 AppendGapBead(int32 bid);
//    int32 PrependGapBead(int32 bid);
//*********************************************************************************

int32 AppendGapBead(int32 bid) {
  // The gap will appear immediately following bid
  Bead *prev = GetBead(beadStore,bid);
  Bead bead;
  char base='-';
  char qv;

  if (prev == NULL ) {
    fprintf(stderr, "AppendGapBead prev==NULL");
    assert(0);
  }
  bead.boffset = GetNumBeads(beadStore);
  bead.soffset = GetNumchars(sequenceStore);
  bead.foffset = prev->foffset+1;
  bead.up = -1;
  bead.down = -1;
  bead.frag_index = prev->frag_index;
  bead.column_index = -1;
  bead.next = prev->next;
  bead.prev = prev->boffset;
  prev->next = bead.boffset;
  qv = *Getchar(qualityStore,prev->soffset);
  if (bead.next != -1 ) {
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

int32 PrependGapBead(int32 bid) {
  // The gap will appear immediately before bid
  Bead *next = GetBead(beadStore,bid);
  Bead bead;
  char base='-';
  char qv;

  assert(next != NULL);
  bead.boffset = GetNumBeads(beadStore);
  bead.soffset = GetNumchars(sequenceStore);
  bead.foffset = next->foffset;
  bead.up = -1;
  bead.down = -1;
  bead.frag_index = next->frag_index;
  bead.column_index = -1;
  bead.next = bid;
  bead.prev = next->prev;
  next->prev = bead.boffset;
  qv = *Getchar(qualityStore,next->soffset);
  if (bead.prev != -1 ) {
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

//  The corresponding SetGappedFragmentPositions exists in version
//  1.64, from mid-September 2007.

static
int
SetUngappedFragmentPositions(FragType type,int32 n_frags, MultiAlignT *uma) {
  char *consensus                  = Getchar(uma->consensus,0);
  VA_TYPE(int32) *gapped_positions = NULL;

  int num_frags;
  int num_unitigs;
  int32 ifrag,first_frag,last_frag;
  IntMultiPos *frag;
  IntUnitigPos *unitig;
  CNS_AlignedContigElement epos;
  HashTable_AS *unitigFrags;
  int hash_rc;

  num_frags   = GetNumIntMultiPoss(uma->f_list);
  num_unitigs = GetNumIntUnitigPoss(uma->u_list);

  unitigFrags = CreateScalarHashTable_AS(2 * (num_frags + num_unitigs));

  //  Earlier versions of this routine were extremely paranoid, and
  //  checked that a fragment/unitig ended at the end of the conitg --
  //  the thinking being that then the contig is covered by fragments.
  //  At least, that's the best I can guess from the bizarre code.
  //  This test was occasionally false, if we are called from
  //  ReplaceEndUnitigInContig.

  {
    int num_columns   = GetMultiAlignLength(uma);
    int ungapped_pos  = 0;
    int i;

    gapped_positions = CreateVA_int32(num_columns+1);

    for (i=0; i<num_columns+1;i++) {
      SetVA_int32(gapped_positions, i, &ungapped_pos);

      if (consensus[i] != '-')
        ungapped_pos++; 
    }
  }

  first_frag=GetNumCNS_AlignedContigElements(fragment_positions);

  for (ifrag=0;ifrag<num_frags;ifrag++){
    frag = GetIntMultiPos(uma->f_list,ifrag);
    epos.frg_or_utg = CNS_ELEMENT_IS_FRAGMENT;
    epos.idx.fragment.frgIdent = frag->ident;
    if (ExistsInHashTable_AS(unitigFrags, frag->ident, 0)) {
      fprintf(stderr,"SetUngappedFragmentPositions()-- ident %d already in hashtable\n",frag->ident);
      assert(0);
    }
    if (HASH_SUCCESS != InsertInHashTable_AS(unitigFrags, frag->ident, 0, 1, 0)) {
      fprintf(stderr,"SetUngappedFragmentPositions()-- Failure to insert ident %d in hashtable\n",frag->ident); 
      assert(0);
    }
    epos.idx.fragment.frgType = frag->type;
    epos.idx.fragment.frgContained = frag->contained;
    epos.idx.fragment.frgInUnitig = (type == AS_CONTIG)?-1:uma->maID;
    epos.idx.fragment.frgSource = frag->sourceInt;
    epos.position.bgn = *Getint32(gapped_positions,frag->position.bgn);
    epos.position.end = *Getint32(gapped_positions,frag->position.end);
    if(epos.position.bgn==epos.position.end){
      fprintf(stderr,"SetUngappedFragmentPositions()-- Encountered bgn==end==" F_COORD " in ungapped coords within SetUngappedFragmentPositions for " F_CID "(gapped coords " F_COORD "," F_COORD ")\n",
              epos.position.bgn,frag->ident,frag->position.bgn,frag->position.end);
      assert(frag->position.bgn!=frag->position.end);
      if(frag->position.bgn<frag->position.end){
        if(epos.position.bgn>0)
          epos.position.bgn--;
        else
          epos.position.end++;
      } else {
        if(epos.position.end>0)
          epos.position.end--;
        else
          epos.position.bgn++;
      }	 
      fprintf(stderr,"SetUngappedFragmentPositions()--   Reset to " F_COORD "," F_COORD "\n",
              epos.position.bgn,
              epos.position.end);
    }
    AppendVA_CNS_AlignedContigElement (fragment_positions,&epos);
  }
  last_frag = GetNumCNS_AlignedContigElements(fragment_positions)-1;
  for (ifrag=0;ifrag<num_unitigs;ifrag++){
    unitig = GetIntUnitigPos(uma->u_list,ifrag);
    epos.frg_or_utg = CNS_ELEMENT_IS_UNITIG;
    epos.idx.unitig.utgIdent = unitig->ident;
    epos.idx.unitig.utgType = unitig->type;
    epos.idx.unitig.utgFirst = first_frag;
    epos.idx.unitig.utgLast = last_frag;
    //epos.contained = 0;
    //epos.source = NULL;
    epos.position.bgn = *Getint32(gapped_positions,unitig->position.bgn);
    epos.position.end = *Getint32(gapped_positions,unitig->position.end);
    AppendVA_CNS_AlignedContigElement(fragment_positions,&epos);
  }
  if (type != AS_CONTIG) { 
    Fragment *anchor = GetFragment(fragmentStore,0);
    CNS_AlignedContigElement *anchor_frag;

    if ( anchor != NULL && anchor->type == AS_CONTIG ) {
      // mark fragments in "anchoring" contig that belong to this unitig
      uint32 first_id,last_id;
      int in_unitig_frags=0;
      first_id = GetCNS_AlignedContigElement(fragment_positions,first_frag)->idx.fragment.frgIdent;
      last_id = GetCNS_AlignedContigElement(fragment_positions,last_frag)->idx.fragment.frgIdent;
      anchor_frag=GetCNS_AlignedContigElement(fragment_positions,anchor->components);
      for (ifrag=0;ifrag<anchor->n_components;ifrag++,anchor_frag++) { 
        if ( anchor_frag->frg_or_utg == CNS_ELEMENT_IS_FRAGMENT ) {
          if (ExistsInHashTable_AS(unitigFrags, anchor_frag->idx.fragment.frgIdent, 0)) {
            anchor_frag->idx.fragment.frgInUnitig=uma->maID;
            in_unitig_frags++;
          }
        }
      }
      fprintf(stderr,"SetUngappedFragmentPositions()-- Marked %d fragments as belonging to unitig %d\n",in_unitig_frags,uma->maID);
    }
  }
  DeleteHashTable_AS(unitigFrags);
  DeleteVA_int32(gapped_positions);
  return first_frag;
}


static
void PrintIMPInfo(FILE *print, int32 nfrags, IntMultiPos *imps) {
  int i;
  uint32 bgn,end;
  for (i=0;i<nfrags;i++) {
    bgn=imps->position.bgn; 
    end=imps->position.end; 
    if ( bgn < end ) {
      fprintf(print,"%12d F %c %10d, %10d -->\n",imps->ident,imps->type,bgn,end);
    } else {
      fprintf(print,"%12d F %c %10d, %10d <--\n",imps->ident,imps->type,end,bgn);
    }
    imps++;
  }
}

static
void PrintIUPInfo(FILE *print, int32 nfrags, IntUnitigPos *iups) {
  int i;
  uint32 bgn,end;
  for (i=0;i<nfrags;i++) {
    bgn=iups->position.bgn; 
    end=iups->position.end; 
    if ( bgn < end ) {
      fprintf(print,"%12d U %c %10d, %10d -->\n",iups->ident,iups->type,bgn,end);
    } else {
      fprintf(print,"%12d U %c %10d, %10d <--\n",iups->ident,iups->type,end,bgn);
    }
    iups++;
  }
}


//*********************************************************************************
// Add a fragment to the basic local store for fragment data
//*********************************************************************************

static
int32 AppendFragToLocalStore(FragType          type,
                             int               iid,
                             int               complement,
                             int               contained,
                             UnitigType        utype,
                             MultiAlignStoreT *multialignStore) {

  char seqbuffer[AS_READ_MAX_LEN+1];
  char qltbuffer[AS_READ_MAX_LEN+1];
  char *sequence = NULL,*quality = NULL;
  static VA_TYPE(char) *ungappedSequence=NULL,*ungappedQuality=NULL;
  Fragment fragment;
  uint clr_bgn, clr_end;
  static fragRecord fsread;  //  static for performance only

  if (ungappedSequence== NULL ) {
    ungappedSequence = CreateVA_char(0);
    ungappedQuality = CreateVA_char(0);
  } else {
    ResetVA_char(ungappedSequence);
    ResetVA_char(ungappedQuality);
  }
  switch (type) {
    case AS_READ:
    case AS_EXTR:
    case AS_TRNR:
      getFrag(gkpStore,iid,&fsread,FRAG_S_QLT);

      clr_bgn = getFragRecordClearRegionBegin(&fsread, clear_range_to_use);
      clr_end = getFragRecordClearRegionEnd  (&fsread, clear_range_to_use);

      strcpy(seqbuffer, getFragRecordSequence(&fsread));
      strcpy(qltbuffer, getFragRecordQuality(&fsread));

      fragment.uid = getFragRecordUID(&fsread);
      //getReadType_ReadStruct(&fsread, &fragment.type);
      fragment.type = AS_READ;
      fragment.source = NULL;
      seqbuffer[clr_end] = '\0';
      qltbuffer[clr_end] = '\0';
      sequence = &seqbuffer[clr_bgn];
      quality = &qltbuffer[clr_bgn];
      fragment.length = (int32) (clr_end - clr_bgn);
      fragment.n_components = 0;  // no component frags or unitigs
      fragment.components = -1;
      break;
    case AS_UNITIG:
    case AS_CONTIG:
      {
        MultiAlignT *uma = NULL;
        if (USE_SDB) {
          assert(sequenceDB != NULL);
          uma = loadMultiAlignTFromSequenceDB(sequenceDB, iid, type == AS_UNITIG);
        } else { 
          uma = GetMultiAlignInStore(multialignStore,iid); 
        }
        if (uma == NULL)
          fprintf(stderr,"Lookup failure in CNS: Unitig %d could not be found in multiAlignStore / unitigStore.\n",iid);
        assert(uma != NULL);

        if (type == AS_CONTIG  && ALIGNMENT_CONTEXT != AS_MERGE) {
          sequence = Getchar(uma->consensus,0);
          quality = Getchar(uma->quality,0);
          fragment.length = GetMultiAlignLength(uma);
        } else {
          GetMultiAlignUngappedConsensus(uma, ungappedSequence, ungappedQuality);
          sequence = Getchar(ungappedSequence,0);
          quality = Getchar(ungappedQuality,0);
          fragment.length = GetMultiAlignUngappedLength(uma);
        }

        if (type == AS_UNITIG) { 
          fragment.utype = utype;
        } else {
          fragment.utype = AS_OTHER_UNITIG; 
        }
  
        //if (type == AS_UNITIG) {
        //  fprintf(stderr,"Unitig fragments for unitig %d:\n",iid);
        //  PrintIMPInfo(stderr,GetNumIntMultiPoss(uma->f_list), GetIntMultiPos(uma->f_list,0) );
        //}

        fragment.n_components = GetNumIntMultiPoss(uma->f_list) + GetNumIntUnitigPoss(uma->u_list);
        fragment.components   = SetUngappedFragmentPositions(type, fragment.n_components, uma);
        break;
      }
    default:
      {
        fprintf(stderr, "AppendFragToLocalStore invalid FragType");
        assert(0);
      }
  }
  if (complement) {
    SequenceComplement(sequence, quality);
  }
  fragment.lid = GetNumFragments(fragmentStore);
  fragment.iid = iid;
  fragment.type = type;
  fragment.complement = complement;
  fragment.contained = contained; 
  fragment.deleted = 0; 
  fragment.manode = -1; 
  fragment.sequence = GetNumchars(sequenceStore);
  fragment.quality = GetNumchars(qualityStore);
  fragment.beads = GetNumBeads(beadStore);
  AppendRangechar(sequenceStore, fragment.length + 1, sequence);
  AppendRangechar(qualityStore, fragment.length + 1, quality);

  {
    Bead bead;
    int32 boffset = fragment.beads;
    int32 soffset = fragment.sequence;
    int32 foffset;

    bead.up = -1;
    bead.down = -1;
    bead.frag_index = fragment.lid;
    bead.column_index = -1;

    for (foffset = 0; foffset < fragment.length; foffset++ ) {
      bead.foffset = foffset;
      bead.boffset = boffset+foffset;
      bead.soffset = soffset+foffset;
      bead.next = bead.boffset +1;
      bead.prev = bead.boffset -1;
      if ( foffset == fragment.length-1) bead.next = -1;
      if ( foffset == 0 ) bead.prev = -1;
      SetVA_Bead(beadStore,boffset+foffset,&bead);
    }
  }
  AppendVA_Fragment(fragmentStore,&fragment);
  //DeleteVA_char(ungappedSequence);
  //DeleteVA_char(ungappedQuality);
  return fragment.lid;
}

int32 AppendArtificialFragToLocalStore(FragType type, int32 iid, int complement,int32 contained,
                                       UnitigType utype, char *seq, char *qlt, int len);

//*********************************************************************************
// Basic manipulation of Bead data 
//*********************************************************************************

int32 AlignBead(int32 cid, int32 bid) {
  Column *column=GetColumn(columnStore,cid);
  Bead *call, *first, *align;
  assert(column != NULL );
  call = GetBead(beadStore,column->call);
  assert(call != NULL);
  first = GetBead(beadStore,call->down);
  assert(first != NULL);
  align = GetBead(beadStore,bid);
  assert(align != NULL);
  align->down = first->boffset;
  align->up = call->boffset;
  call->down = align->boffset;
  first->up = align->boffset;
  align->column_index = cid;
  IncBaseCount(&column->base_count,*Getchar(sequenceStore,align->soffset));
  return bid;
}

int32 UnAlignBead(int32 bid) {
  // remove bid from it's column, returning the next bead up in the column
  Bead *bead = GetBead(beadStore,bid);
  Bead *upbead;
  Column *column;
  char bchar;
  assert(bead != NULL);
  if (bead->column_index == -1 ) return -1;
  column = GetColumn(columnStore,bead->column_index);
  upbead = GetBead(beadStore,bead->up);
  bchar = *Getchar(sequenceStore,bead->soffset);
  upbead->down = bead->down;
  if (bead->down != -1 ) {
    GetBead(beadStore, bead->down)->up = upbead->boffset;
  }
  DecBaseCount(&column->base_count,bchar);
  bead->up = -1;
  bead->down = -1;
  bead->column_index = -1;
  return upbead->boffset;
}

int32 RemoveBeadFromFragment(int32 bid) {
  // remove bid from it's fragment, returning the next bead in the fragment
  Bead *bead = GetBead(beadStore,bid);
  Bead *nextbead;
  Bead *prevbead;
  assert(bead != NULL);
  if ( bead->next > -1) {
    nextbead = GetBead(beadStore,bead->next);
    nextbead->prev = bead->prev;
  } 
  if ( bead->prev > -1) {
    prevbead = GetBead(beadStore,bead->prev);
    prevbead->next = bead->next;
  } 
  return bead->next;
}

int32 UnAlignFragment(int32 fid) {
  Fragment *frag=GetFragment(fragmentStore,fid);
  Bead *bead; 
  int32 next_bid;
  assert(frag != NULL);
  bead = GetBead(beadStore,frag->beads);    
  assert(bead != NULL);
  next_bid = bead->next;  
  while (next_bid > 0 ) {
    UnAlignBead(bead->boffset);
    if ( *Getchar(sequenceStore,bead->soffset) == '-' ) {
      // remove the gap bead from the fragment 
      RemoveBeadFromFragment(bead->boffset);
    }
    bead = GetBead(beadStore,next_bid);
    next_bid = bead->next;
  }
  UnAlignBead(bead->boffset);
  frag->deleted = 1;
  frag->manode = -1;
  return 1;
}

int32 UnAlignTrailingGapBeads(int32 bid) {
  // remove bid from it's column, returning the prev or next bead in the fragment
  Bead *bead = GetBead(beadStore,bid);
  Bead *upbead,*prevbead,*nextbead;
  int32 anchor;
  Column *column;
  char bchar;
  assert(bead != NULL);
  // find direction to remove
  anchor = bead->prev;
  while ( bead->next != -1 && *Getchar(sequenceStore,(GetBead(beadStore,bead->next))->soffset) == '-' ) {
    bead = GetBead(beadStore,bead->next);
  }
  if (bead->next != -1 ) {
    anchor = bead->next;
    while (bead->prev != -1 && *Getchar(sequenceStore,(GetBead(beadStore,bead->prev))->soffset) == '-' ) {
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
    if (bead->down != -1 ) {
      GetBead(beadStore, bead->down)->up = upbead->boffset;
    }
    DecBaseCount(&column->base_count,bchar);
    bead->up = -1;
    bead->down = -1;
    bead->column_index = -1;
    if ( bead->next == -1 ) {
      prevbead = GetBead(beadStore,bead->prev);
      prevbead->next = -1;
      bead->prev = -1;
      bead = GetBead(beadStore,prevbead->boffset);
    } else {
      nextbead = GetBead(beadStore,bead->next);
      nextbead->prev = -1;
      bead->next = -1;
      bead = GetBead(beadStore,nextbead->boffset);
    }
  }
  return anchor;
}

int32 LateralExchangeBead(int32 lid, int32 rid) {
  // returned value is new leftmost bead id
  Bead *leftbead, *rightbead, *ibead;
  Column *leftcolumn, *rightcolumn;
  Bead rtmp; // this is just some tmp space for the swap
  char leftchar, rightchar;
  leftbead = GetBead(beadStore,lid);
  rightbead = GetBead(beadStore,rid);
  assert(leftbead != NULL);
  assert(rightbead != NULL);
  leftcolumn = GetColumn(columnStore,leftbead->column_index);
  rightcolumn = GetColumn(columnStore,rightbead->column_index);
  assert(leftcolumn != NULL);
  assert(rightcolumn != NULL);
  leftchar = *Getchar(sequenceStore,leftbead->soffset);
  rightchar = *Getchar(sequenceStore,rightbead->soffset);

  // now, verify that left and right are either
  // a) neighbors, or b) have only '-'s intervening
  ibead = leftbead;
  while ( ibead->next > -1) {
    ibead = GetBead(beadStore,ibead->next);
    if (ibead->boffset == rid ) break;
  
    if( *Getchar(sequenceStore,ibead->soffset) != '-') {
      fprintf(stderr, "LateralExchangeBead exchangebead!='-'");
      assert(0);
    }
  }
  rtmp = *rightbead;
  rightbead->up = leftbead->up;
  rightbead->down = leftbead->down;
  rightbead->prev = leftbead->prev;
  rightbead->next = leftbead->next;
  if ( rightbead->up != -1 ) (GetBead(beadStore,rightbead->up))->down = rid;
  if ( rightbead->down != -1)  (GetBead(beadStore,rightbead->down))->up = rid;
  if ( rightbead->prev != -1)  (GetBead(beadStore,rightbead->prev))->next = rid;
  leftbead->up = rtmp.up;
  leftbead->down = rtmp.down;
  leftbead->next = rtmp.next;
  leftbead->prev = rtmp.prev;
  if ( leftbead->up != -1 ) (GetBead(beadStore,leftbead->up))->down = lid;
  if ( leftbead->down != -1)  (GetBead(beadStore,leftbead->down))->up = lid;
  if ( leftbead->next != -1)  (GetBead(beadStore,leftbead->next))->prev = lid;
  // now, handle separately cases of a) left and right are adjacent, and b) gaps intervene
  if ( rtmp.prev == lid) {
    rightbead->next = lid;
    leftbead->prev = rid;
  } else {
    if ( rightbead->next != -1)  (GetBead(beadStore,rightbead->next))->prev = rid;
    if ( leftbead->prev != -1)  (GetBead(beadStore,leftbead->prev))->next = lid;
  }

  rightbead->column_index = leftbead->column_index;
  leftbead->column_index = rtmp.column_index;
  // change basecounts for affected columns
  DecBaseCount(&leftcolumn->base_count,leftchar);
  IncBaseCount(&leftcolumn->base_count,rightchar);
  DecBaseCount(&rightcolumn->base_count,rightchar);
  IncBaseCount(&rightcolumn->base_count,leftchar);
  return rid;
}

int32 LeftEndShiftBead(int32 bid, int32 eid) {
  //  Relationship must be one of:
  //
  //  a) end gap moving left:
  // 
  //     X > A > B > C > ... > -   becomes  X - A B C ... 
  //         ^________________/
  //
  //  b) non-gap moving left across only gap characters
  //    (more efficient special case, since first gap and last
  //     character can just be exchanged)
  // 
  //     X > - > - > - > ... > A   becomes  X A - - - ...
  //         ^________________/

  Bead *shift = GetBead(beadStore,eid);
  int32 aid = (GetBead(beadStore,bid))->prev;
  assert(shift != NULL);
  if ( *Getchar(sequenceStore,shift->soffset) != '-' ) {
    // assume first and internal characters are gaps
    LateralExchangeBead(bid, eid);
    return eid;
  }   else {
    while ( shift->prev != aid ) {
      LateralExchangeBead(shift->prev, shift->boffset); 
    }
    return shift->boffset;
  }
}

int32 RightEndShiftBead(int32 bid, int32 eid) {
  //  Relationship must be one of:
  //
  //  a) end gap moving left:
  // 
  //      - > A > B > ... > C > X  becomes  A B ... C - X
  //      \_________________^
  //
  //  b) non-gap moving left across only gap characters
  //    (more efficient special case, since first gap and last
  //     character can just be exchanged)
  // 
  //      A > - > - > ... > - > X  becomes  - - - ... A X
  //       \________________^

  Bead *shift = GetBead(beadStore,bid);
  int32 aid = (GetBead(beadStore,eid))->next;
  int32 rid; 
  assert(shift != NULL);
  if ( *Getchar(sequenceStore,shift->soffset) != '-' ) {
    // assume last and internal characters are gaps
    LateralExchangeBead(bid, eid);
    return eid;
  }   else {
    rid = shift->next;
    while ( shift->next != aid ) {
      LateralExchangeBead(shift->boffset, shift->next); 
    }
    return rid;
  }
}
   
//*********************************************************************************
// Basic manipulation of Column data 
//*********************************************************************************

Column * CreateColumn(int32 bid){
  // create a new column, seeded with the bead bid
  Column column;
  Bead call;
  Bead *head;
  
  column.lid = GetNumColumns(columnStore);
  column.prev = -1;
  column.next = -1;
  column.call = GetNumBeads(beadStore);
  column.ma_index = -1;
  ResetBaseCount(&column.base_count);
  call.boffset = column.call;
  call.foffset = 0;
  call.soffset = GetNumchars(sequenceStore);
  call.down = bid;
  call.up = -1;
  call.prev = -1;
  call.next = -1;
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
  return GetColumn(columnStore, column.lid);
}

int AddColumnToMANode(int32 ma, Column column){
  MANode *manode = GetMANode(manodeStore,ma);
  Appendint32(manode->columns,&column.lid);
  if (column.next == -1 ) {
    manode->last = column.lid;
  }
  if (column.prev == -1 ) {
    manode->first = column.lid;
  }
  return 1;
}

int32 ColumnAppend(int32 cid, int32 bid) {
  // bid is the offset of the Bead seeding the column
  Column *column;
  Column *prev;
  Column *next;
  ColumnBeadIterator ci;
  int32 nid;
  Bead *bead = GetBead(beadStore,bid);
  Bead *call,*prevcall,*nextcall;
  // make sure this bead exists before continuing
  assert(bead != NULL);
  column = CreateColumn(bid);
  // make sure this column exists before continuing
  assert(column != NULL);
  call = GetBead(beadStore,column->call);
  prev = GetColumn(columnStore,cid);
  prevcall = GetBead(beadStore,prev->call);
  column->next = prev->next;
  column->prev = cid;
  call->next = prevcall->next;
  call->prev = prevcall->boffset;
  prev->next = column->lid;
  prevcall->next = call->boffset; 
  if ( column->next != -1 ) {
    next = GetColumn(columnStore,column->next);
    next->prev = column->lid;
  }
  if ( call->next != -1 ) {
    nextcall = GetBead(beadStore,call->next);
    nextcall->prev = call->boffset;
  }
  if(! CreateColumnBeadIterator(cid,&ci)){
    fprintf(stderr, "ColumnAppend CreateColumnBeadIterator failed");
    assert(0);
  }
  while ( (nid = NextColumnBead(&ci)) != -1 ) {
    bead = GetBead(beadStore,nid);
    if ( bead->next != -1 && bead->next != bid) {
      AlignBead(column->lid,AppendGapBead(nid));
    }
  }
  column->ma_id =  prev->ma_id;
  column->ma_index =  prev->ma_index + 1;
  AddColumnToMANode(column->ma_id,*column);
  return column->lid;
}

int32 ColumnPrepend(int32 cid, int32 bid) {
  // bid is the offset of the Bead seeding the column
  Column *column;
  Column *prev;
  Column *next;
  ColumnBeadIterator ci;
  int32 nid;
  Bead *bead = GetBead(beadStore,bid);
  Bead *call,*prevcall,*nextcall;
  // make sure this bead exists before continuing
  assert(bead != NULL);
  column = CreateColumn(bid);
  // make sure this column exists before continuing
  assert(column != NULL);
  call = GetBead(beadStore,column->call);
  next = GetColumn(columnStore,cid);
  nextcall = GetBead(beadStore,next->call);
  column->prev = next->prev;
  column->next = cid;
  call->prev = nextcall->prev;
  call->next = nextcall->boffset;
  next->prev = column->lid;
  nextcall->prev = call->boffset; 
  if ( column->prev != -1 ) {
    prev = GetColumn(columnStore,column->prev);
    prev->next = column->lid;
  }
  if ( call->prev != -1 ) {
    prevcall = GetBead(beadStore,call->prev);
    prevcall->next = call->boffset;
  }
  if(! CreateColumnBeadIterator(cid,&ci)){
    fprintf(stderr, "ColumnPrepend CreateColumnBeadIterator failed");
    assert(0);
  }
  while ( (nid = NextColumnBead(&ci)) != -1 ) {
    bead = GetBead(beadStore,nid);
    if ( bead->prev != -1 && bead->prev != bid) {
      AlignBead(column->lid,PrependGapBead(nid));
    }
  }
  column->ma_id =  next->ma_id;
  column->ma_index =  next->ma_index - 1;
  AddColumnToMANode(column->ma_id,*column);
  if ( column->prev == -1 ) {
    GetMANode(manodeStore,column->ma_id)->first = column->lid;
  }
  return column->lid;
}

int32 FirstColumn(int32 mid, int32 bid) {
  // bid is the offset of the Bead seeding the column
  Column *column;
  Bead *bead = GetBead(beadStore,bid);
  AssertPtr(bead);
  assert(bead != NULL);
  column = CreateColumn(bid);
  assert(column != NULL);
  column->ma_id =  mid;
  column->ma_index =  0;
  AddColumnToMANode(mid,*column);
  return column->lid;
}


int MergeCompatible(int32 cid) {
  // test for Level 1 (neighbor) merge compatibility of cid with right neighbor
  // and merge if compatible
  Column *column, *merge_column;
  Bead *cbead, *mbead;
  char cchar, mchar;
  int32 mid; // id of bead to merge
  int mergeok = 1;
  column = GetColumn(columnStore,cid);
  assert(column != NULL);
  if (column->next == -1) return 0;
  merge_column = GetColumn(columnStore,column->next);
  assert(merge_column != NULL);
  cbead = GetBead(beadStore,column->call);
  while (mergeok && cbead->down != - 1) {
    cbead = GetBead(beadStore,cbead->down);
    mid = cbead->next;
    if ( mid == -1 ) continue;
    mbead =  GetBead(beadStore,mid);
    cchar = *Getchar(sequenceStore,cbead->soffset);
    mchar = *Getchar(sequenceStore,mbead->soffset);
    if ( ! ((cchar == '-') | ( mchar == '-')) ) mergeok = 0;
  }
  if ( mergeok ) { // go ahead and do merge (to left)
    cbead = GetBead(beadStore,column->call);
    while (cbead->down != - 1) {
      cbead = GetBead(beadStore,cbead->down);
      mid = cbead->next;
      if ( mid == -1 ) continue;
      mbead =  GetBead(beadStore,mid);
      cchar = *Getchar(sequenceStore,cbead->soffset);
      mchar = *Getchar(sequenceStore,mbead->soffset);
      if ( ((cchar == '-') && ( mchar != '-')) ) {
        mid = LeftEndShiftBead(cbead->boffset,mid);
        cbead = GetBead(beadStore,mid);
      }
    }
    //wrap up with any trailing guys from right that need to move left
    while ( GetDepth(merge_column) != GetBaseCount(&merge_column->base_count,'-') ) {
      mbead = GetBead(beadStore,merge_column->call);
      while ( mbead->down != -1 ) { 
        mbead = GetBead(beadStore,mbead->down);
        if ( *Getchar(sequenceStore,mbead->soffset) != '-' ) {
          UnAlignBead(mbead->boffset);
          AlignBead(cbead->column_index,mbead->boffset);
          cbead = mbead;
          break;
        }
        if(GetDepth(merge_column) <= 0){
          fprintf(stderr, "MergeCompatible empty column");
          assert(0);
        }
      }
    }
    return 1;
  } else { return 0; }
}

int AverageDepth(int32 bgn, int32 end) {
  int depth=0;
  int ncolumns = 0;
  ColumnIterator ci; 
  int nid;
  if ( !CreateColumnIterator(bgn,&ci)) {
    fprintf(stderr, "AverageDepth CreateColumnIterator failed");
    assert(0);
  }
  while ( (nid = NextColumn(&ci)) != -1 ) {
    if (nid == end) break;
    depth += GetDepth(GetColumn(columnStore,nid));
    ncolumns++;
  }
  if ( ncolumns == 0) return 0;
  return (int) depth/ncolumns;
}
      
void ShowColumn(int32 cid) {
  Column *column = GetColumn(columnStore,cid);
  Bead *call;
  Bead *bead;
  FragType type;
  UnitigType utype;
  ColumnBeadIterator ci;
  int32 bid;
  if(!CreateColumnBeadIterator(cid,&ci)){
    fprintf(stderr, "ShowColumn CreateColumnBeadIterator failed");
    assert(0);
  }
  call = GetBead(beadStore,column->call);
  fprintf(stderr,"\nstore_index: %-20d ( prev: %d next: %d)\n",column->lid,column->prev,column->next);
  fprintf(stderr,"ma_index:    %-20d\n",column->ma_index);
  fprintf(stderr,"------------------\n");
  fprintf(stderr,"composition:\n");
  while ( (bid = NextColumnBead(&ci)) != -1 ) {
    // do something here
    bead = GetBead(beadStore,bid);
    type = GetFragment(fragmentStore,bead->frag_index)->type;
    utype = GetFragment(fragmentStore,bead->frag_index)->utype;
    fprintf(stderr,"             %c /%c (%10d) <-- %d %c/%c\n",*Getchar(sequenceStore,bead->soffset),
            *Getchar(qualityStore,bead->soffset),
            bid,
            bead->frag_index,
            type,
            (type == AS_UNITIG)?utype:' '); 
  }
  fprintf(stderr,"------------------\n");
  fprintf(stderr,"call:        %c /%c\n",toupper(*Getchar(sequenceStore,call->soffset)),*Getchar(qualityStore,call->soffset));
}

char QVInRange(int q) {
  // Ensure that QV returned is within allowable range
  // outliers are rounded down or up to endpoints of range
  // [ CNS_MIN_QV , CNS_MAX_QV ]
  if (  q > CNS_MAX_QV ) {
    return (char) CNS_MAX_QV + '0';
  } else if (  q < CNS_MIN_QV ) {
    return (char) CNS_MIN_QV + '0';
  } else {
    return (char) q + '0';
  }
}

static int
Iid2ReadId(int32 iid, int32 *iids, int nr)
{
  int i;
  for (i=0; i<nr; i++)
    {
      if (iid == iids[i])
        return i;
    }
  return (-1);
}

static int
IsNewRead(int32 iid, int32 *iid_list, int nr)
{
  int i;
  for (i=0; i<nr; i++)
    {
      if (iid_list[i] == iid)
        return 0;
    }
  return 1;
}


//*********************************************************************************
// Function: BaseCalling
// Purpose: Calculate the consensus base for the given column
//*********************************************************************************
int
BaseCall(int32 cid, int quality, double *var, VarRegion  *vreg,
         int target_allele, char *cons_base, int verbose, int get_scores,
         CNS_Options *opp)
{
  /* NOTE: negative target_allele means the the alleles will be used */

  Column *column=GetColumn(columnStore,cid);
  Bead *call = GetBead(beadStore, column->call);
  Bead *bead;
  int best_read_base_count[CNS_NP]  = {0};
  int other_read_base_count[CNS_NP] = {0};
  int guide_base_count[CNS_NP]      = {0};
 
  char bases[CNS_NALPHABET] = {'-', 'A', 'C', 'G', 'T', 'N'};
  int best_read_qv_count[CNS_NP] = {0};
  int other_read_qv_count[CNS_NP] = {0};
  int highest_qv[CNS_NP] = {0};
  int highest2_qv[CNS_NP] = {0};

  int b_read_depth=0, o_read_depth=0, guide_depth=0;
  int score=0;
  int bi;
  int32 bid;
  int32 iid = 0;
  char cqv, cbase;
  int qv = 0;
  static  double cw[CNS_NP];      // "consensus weight" for a given base
  static  double tau[CNS_NP];
  FragType type;
  UnitigType utype;
  ColumnBeadIterator ci;
  int used_surrogate=0;
  int sum_qv_cbase=0, sum_qv_all=0;
  int k;

  vreg->nb = 0;

  //  Make sure that we have valid options here, we then reset the
  //  pointer to the freshly copied options, so that we can always
  //  assume opp is a valid pointer
  //
  CNS_Options  opp_private;
  if (opp == NULL) {
    opp_private.split_alleles   = CNS_OPTIONS_SPLIT_ALLELES_DEFAULT;
    opp_private.smooth_win      = CNS_OPTIONS_SMOOTH_WIN_DEFAULT;
    opp_private.max_num_alleles = CNS_OPTIONS_MAX_NUM_ALLELES;
    opp = &opp_private;
  }


  if(!CreateColumnBeadIterator(cid, &ci)){
    fprintf(stderr, "BaseCall CreateColumnBeadIterator failed");
    assert(0);
  }

  *var = 0.;
  if (quality > 0)
    {
      static int guides_alloc=0;
      static VarArrayBead  *guides;
      static VarArrayBead  *b_reads;
      static VarArrayBead  *o_reads;
      static VarArrayint16 *tied;
      uint32 bmask;
      int    num_b_reads, num_o_reads, num_guides;
      Bead  *gb;
      int    cind;
      double tmpqv;
      int16  bi;
      int    b_read_count = 0;
      int    frag_cov=0;
      int16  max_ind=0;
      double max_cw=0.0;   // max of "consensus weights" of all bases
      double normalize=0.;
      int    nr=0, max_nr=100;
      int32 *column_iid_list = (int32 *)safe_malloc(max_nr*sizeof(int32));

      if (!guides_alloc) {
        guides = CreateVA_Bead(16);
        b_reads  = CreateVA_Bead(16);
        o_reads  = CreateVA_Bead(16);
        tied   = CreateVA_int16(32);
        guides_alloc = 1;
      }
      else {
        ResetBead(guides);
        ResetBead(b_reads);
        ResetBead(o_reads);
        Resetint16(tied);
      }
      for (bi=0;bi<CNS_NP;bi++) {
        tau[bi] = 1.0;
      }

      // Scan a column of aligned bases (=beads).
      // Sort the beads into three groups:
      //      - those corresponding to the reads of the best allele,
      //      - those corresponding to the reads of the other allele and
      //      - those corresponding to non-read fragments (aka guides)
      while ( (bid = NextColumnBead(&ci)) != -1)
        {
          bead =  GetBead(beadStore,bid);
          cbase = *Getchar(sequenceStore,bead->soffset);    // current base
          qv = (int) ( *Getchar(qualityStore,bead->soffset)-'0');
          if ( cbase == 'N' ) {
            // skip 'N' base calls
            // fprintf(stderr,
            //    "encountered 'n' base in fragment data at column cid=%d\n",
            //    cid);
            continue;
          }
          bmask = AMASK[BaseToInt(cbase)];
          type  = GetFragment(fragmentStore,bead->frag_index)->type;
          iid   = GetFragment(fragmentStore,bead->frag_index)->iid;
          k     = Iid2ReadId(iid, vreg->iids, vreg->nr);

          if ((type == AS_READ)   ||
              (type == AS_EXTR)   ||
              (type == AS_TRNR))
            {
              // Filter out "duplicated" reads with the same iid
              if (!IsNewRead(iid, column_iid_list, nr))
                {
                  fprintf(stderr, "Read iid= %d occurs more than once ", iid);
                  fprintf(stderr, "in MSA for contig #%d at pos= %d\n",
                          contig_id, cid);
                  continue;
                }

              column_iid_list[nr] = iid;
              nr++;
              if (nr == max_nr)
                {
                  max_nr += 100;
                  column_iid_list = (int32 *)safe_realloc(column_iid_list,
                                                          max_nr*sizeof(int32));
                }

              // Will be used when detecting alleles
              if (target_allele < 0 && get_scores)
                {
                  vreg->curr_bases[vreg->nb] = cbase;
                  vreg->iids[vreg->nb]  = iid;
                  vreg->nb++;
                  if (vreg->nb == vreg->max_nr)
                    {
                      vreg->max_nr += INITIAL_NR;
                      vreg->curr_bases = (char *)safe_realloc(vreg->curr_bases,
                                                              vreg->max_nr*sizeof(char));
                      vreg->iids = (int32 *)safe_realloc(vreg->iids,
                                                         vreg->max_nr*sizeof(int32));
                    }
                }

              // Will be used when detecting variation
              if (((target_allele < 0)  ||   // use any allele
                   !opp->split_alleles  ||   // use any allele
                   (vreg->nr >  0  &&
                    vreg->reads[k].allele_id == target_allele))) // use the best allele
                {
                  best_read_base_count[BaseToInt(cbase)]++;
                  best_read_qv_count[BaseToInt(cbase)] += qv;
                  AppendBead(b_reads, bead);
                }
              else
                {
                  other_read_base_count[BaseToInt(cbase)]++;
                  other_read_qv_count[BaseToInt(cbase)] += qv;
                  AppendBead(o_reads, bead);
                }

              if (highest_qv[BaseToInt(cbase)] < qv)
                {
                  highest2_qv[BaseToInt(cbase)] = highest_qv[BaseToInt(cbase)];
                  highest_qv[BaseToInt(cbase)] = qv;
                }
              else if (highest_qv[BaseToInt(cbase)] >= qv &&
                       highest2_qv[BaseToInt(cbase)] < qv)
                {
                  highest2_qv[BaseToInt(cbase)] = qv;
                }
            }
          else
            {
              guide_base_count[BaseToInt(cbase)]++;
              AppendBead(guides, bead);
            }

          if ( type != AS_UNITIG ) {
            frag_cov++;
          }
        }

      b_read_depth = GetNumBeads(b_reads);
      o_read_depth = GetNumBeads(o_reads);
      guide_depth  = GetNumBeads(guides);

      // For each base, calculate tau
      // It will be used to calculate cw
      if (b_read_depth > 0)
        {
          for (cind = 0; cind < b_read_depth; cind++)
            {
              gb = GetBead(b_reads, cind);
              cbase = *Getchar(sequenceStore,gb->soffset);
              qv = (int) ( *Getchar(qualityStore,gb->soffset)-'0');
              if ( qv == 0 )
                qv += 5;
              bmask = AMASK[BaseToInt(cbase)];
              for (bi=0;bi<CNS_NP;bi++) {
                if ( (bmask>>bi) & 1 ) {
                  tau[bi]*= PROB[qv];
                } else {
                  tau[bi]*= (double) TAU_MISMATCH * EPROB[qv];
                }
              }
            }
        }
      else
        {
          for (cind = 0; cind < o_read_depth; cind++)
            {
              gb = GetBead(o_reads, cind);
              cbase = *Getchar(sequenceStore,gb->soffset);
              qv = (int) ( *Getchar(qualityStore,gb->soffset)-'0');
              if ( qv == 0 )
                qv += 5;
              bmask = AMASK[BaseToInt(cbase)];
              for (bi=0;bi<CNS_NP;bi++) {
                if ( (bmask>>bi) & 1 ) {
                  tau[bi]*= PROB[qv];
                } else {
                  tau[bi]*= (double) TAU_MISMATCH * EPROB[qv];
                }
              }
            }
        }

      // If there are no reads, use fragments of other types
      if (b_read_depth == 0 && o_read_depth == 0)
        {
          for (cind = 0; cind < guide_depth; cind++) {
            gb = GetBead(guides,cind);
            type  = GetFragment(fragmentStore,gb->frag_index)->type;
            utype = GetFragment(fragmentStore,gb->frag_index)->utype;

            if ( type == AS_UNITIG &&
                 ((utype != AS_STONE_UNITIG &&
                   utype != AS_PEBBLE_UNITIG &&
                   utype != AS_OTHER_UNITIG) || b_read_depth > 0))
              {
                continue;
              }
            used_surrogate=1;
            // only for surrogates, use their basecalls/quality in contig consensus
            cbase = *Getchar(sequenceStore,gb->soffset);
            qv = (int) ( *Getchar(qualityStore, gb->soffset)-'0');
            if ( qv == 0 )
              qv += 5;
            bmask = AMASK[BaseToInt(cbase)];
            for (bi=0; bi<CNS_NP; bi++) {
              if ( (bmask>>bi) & 1 ) {
                tau[bi] *= PROB[qv];
              } else {
                tau[bi] *= (double) TAU_MISMATCH * EPROB[qv];
              }
            }
          }
        }

      max_ind = 0;      // consensus is gap
      max_cw  = 0.0;

      //  This is gross.

      for (bi=0; bi<5; bi++) {
        cw[bi]     = tau[bi] * 0.2;
        normalize += cw[bi];
      }

      if (normalize)
        normalize = 1./normalize;

      // Calculate max_ind as {i | cw[i] -> max_cw}
      // Store all other indexes { i | cw[i] == max_cw } in VA Array tied
      for (bi=0; bi<5; bi++) {
        cw[bi] *= normalize;

        if (cw[bi] > max_cw + ZERO_PLUS) {
          max_ind = bi;
          max_cw = cw[bi];
          Resetint16(tied);
        } else if (DBL_EQ_DBL(cw[bi], max_cw)) {
          Appendint16(tied,&bi);
        }
      }

      // If max_cw == 0, then consensus base call will be a gap
      // (max_ind==0)
      //
      // Otherwise, it will be selected RANDOMLY (!!!) from all
      // {i|cw[i]==max_cw}

      if (DBL_EQ_DBL(max_cw, (double)0.0)) {
        max_ind = 0;      // consensus is gap
      } else {
        if (GetNumint16s(tied)> 0) {
          Appendint16(tied, &max_ind);
          max_ind = *Getint16(tied,1);
          max_cw = cw[max_ind];
        }
      }








      if ( verbose ) {
        fprintf(stdout,"calculated probabilities:\n");
        for (bi=0;bi<CNS_NP;bi++)
          fprintf(stdout,"%c = %16.8f %c\n", RALPHABET[bi],cw[bi], (bi == max_ind) ? '*' : ' ');
      }

      // Set the consensus base quality value
      cbase = RALPHABET[max_ind];      // consensus base
      if (DBL_EQ_DBL(max_cw, (double)1.0)) {
        cqv = CNS_MAX_QV+'0';
        Setchar(qualityStore, call->soffset, &cqv);
      }
      else
        {
          if ( frag_cov != 1 || used_surrogate)
            {
              tmpqv =  -10.0 * log10(1.0-max_cw);
              qv = DBL_TO_INT(tmpqv);
              if ((tmpqv - qv)>=.50)
                qv++;
            }
          cqv = QVInRange(qv);
        }


      *cons_base = cbase;
      if (target_allele <  0 || target_allele == vreg->alleles[0].id)
        {
          Setchar(sequenceStore, call->soffset, &cbase);
          Setchar(qualityStore, call->soffset, &cqv);
        }

      // Detecting variation
      for (bi=0; bi<CNS_NALPHABET-1; bi++)
        b_read_count += best_read_base_count[bi];

      for (bi=0; bi<CNS_NALPHABET-1; bi++) {
        // NALAPHBET-1 to exclude "n" base call
        bmask = AMASK[bi];  // mask for indicated base
        if ( ! ((bmask>>max_ind) & 1) ) {
          // penalize only if base in not represented in call
          score += best_read_base_count[bi] + other_read_base_count[bi]
            + guide_base_count[bi];
        }
        /* To be considered, varied base should be confirmed by another base 
         * and either overall quality should be >= 60
         * (Granger's suggestion - GD)
         */
        if (*cons_base != '-' && bases[bi] != '-' &&
            best_read_base_count[bi] >  1 &&
            (float)best_read_qv_count[bi]/(float)best_read_base_count[bi]
            >= MIN_AVE_QV_FOR_VARIATION)
          {
            sum_qv_all += best_read_qv_count[bi];
            if (IntToBase(bi) == cbase)
              sum_qv_cbase = best_read_qv_count[bi];
          }
        else if ((*cons_base == '-' || bases[bi] == '-') &&
                 best_read_base_count[bi] >  1 &&
                 highest_qv[bi] + highest2_qv[bi] >= MIN_SUM_QVS_FOR_VARIATION)
          {
            sum_qv_all += best_read_qv_count[bi];
            if (IntToBase(bi) == cbase)
              sum_qv_cbase = best_read_qv_count[bi];
          } 
      }
      if ((b_read_count <= 1 ) || (sum_qv_all == 0))
        {
          *var = 0.;
          if (opp->smooth_win > 0 && cbase == '-')
            *var = -2;
        }
      else
        {
          *var = 1. - (double)sum_qv_cbase / (double)sum_qv_all;
          if (opp->smooth_win > 0 && cbase == '-')
            {
              *var = - (*var);
            }
        }
      safe_free(column_iid_list);
      return score;
    }
  else if (quality == 0 )
    {
      int max_count=0,max_index=-1,tie_count=0;
      int tie_breaker, max_tie, i;

      if(!CreateColumnBeadIterator(cid,&ci)) {
        fprintf(stderr, "BaseCount CreateColumnBeadIterator failed");
        assert(0);
      }
      while ( (bid = NextColumnBead(&ci)) != -1 ) {
        bead = GetBead(beadStore,bid);
        cbase = *Getchar(sequenceStore,bead->soffset);
        qv = (int) ( *Getchar(qualityStore, bead->soffset)-'0');
        type = GetFragment(fragmentStore,bead->frag_index)->type;
        if (type  != AS_READ &&
            type  != AS_EXTR &&
            type  != AS_TRNR ) {
          guide_base_count[BaseToInt(cbase)]++;
        }
        else {
          best_read_base_count[BaseToInt(cbase)]++;
        }
      }
      for (i=0; i<CNS_NALPHABET; i++) {
        if (best_read_base_count[i]+guide_base_count[i] > max_count) {
          max_count = best_read_base_count[i] + guide_base_count[i];
          max_index = i;
        }
      }
      if ( best_read_base_count[max_index] + guide_base_count[max_index] >
           (b_read_depth                 + guide_depth)/2 )
        {
          tie_count = 0;
        }
      else
        {
          for (i=0;i<CNS_NALPHABET;i++)
            {
              if (best_read_base_count[i]+guide_base_count[i] == max_count)
                {
                  max_index = i;
                  tie_count++;
                }
            }
        }
      max_tie=-1;
      if ( tie_count > 1 ) {
        for (i=1;i<CNS_NALPHABET;i++)
          {     /* i starts at 1 to prevent ties */
            /* from being broken with '-'    */
            if ( best_read_base_count[i]+guide_base_count[i] == max_count )
              {
                /* Break unresolved ties with random numbers: */
                tie_breaker = random();
                if (tie_breaker > max_tie) {
                  max_tie = tie_breaker;
                  max_index = i;
                }
              }
          }
      }
      cbase=toupper(RALPHABET[max_index]);
      Setchar(sequenceStore, call->soffset, &cbase);
      cqv = 0 + '0';
      Setchar(qualityStore, call->soffset, &cqv);
      for (bi=0;bi<CNS_NALPHABET;bi++) {
        if (bi != BaseToInt(cbase))
          {
            score += best_read_base_count[bi]+guide_base_count[bi];
          }
      }
      return score;
    }
  else if (quality == -1 ) {
    // here, just promote the aligned fragment's seq and quality to the basecall
    char bqv;
    bid = NextColumnBead(&ci);
    bead =  GetBead(beadStore,bid);
    cbase = *Getchar(sequenceStore, bead->soffset);
    bqv  = *Getchar(qualityStore,bead->soffset);
    Setchar(sequenceStore, call->soffset, &cbase);
    Setchar(qualityStore, call->soffset, &bqv);
    return score;
  }
  return score;
}


static void
SetDefault(VarRegion  *vreg)
{
  vreg->nr = 0;
}

static void
SmoothenVariation(double *var, int len, int window)
{
  int i;
  double *y = (double *)safe_malloc(len * sizeof(double));

  if (window == 0)
    {
      safe_free(y);
      return;
    }
  for (i=0; i<len; i++)
    {
      int j, left_win=0, right_win=0;
      double sum_var = (var[i]>ZERO_MINUS)?var[i]:((var[i]<-1.)?0.:-var[i]);
      int max_left_win = window/2;
      int max_right_win = window - max_left_win;

      j = i-1;
      while (j>=0 && left_win <= max_left_win)
        {
          if (var[j] > ZERO_MINUS)   // consensus is not gap
            {
              left_win++;
              sum_var += var[j];
            }
          else if (var[j]>-1+ZERO_MINUS) // consensus is gap, var != 0
            sum_var -= var[j];

          j--;
        }
      j = i+1;

      while (j<len && right_win <= max_right_win)
        {
          if (var[j] > ZERO_MINUS)  // consensus is not gap
            {
              right_win++;
              sum_var += var[j];
            }
          else if (var[j]>-1+ZERO_MINUS) // consensus is gap, var != 0
            sum_var -= var[j];
          j++;
        }
      y[i] = (left_win+right_win > 0) ?
        sum_var/(double)(left_win+right_win) : var[i];
    }
  for (i=0; i<len; i++)
    {
      var[i] = y[i];
    }
  safe_free(y);
}

static void
GetReadIidsAndNumReads(int cid, VarRegion  *vreg)
{
  int      cind;
  int16    bi;
  Column  *column=GetColumn(columnStore,cid);
  Bead    *call = GetBead(beadStore, column->call);
  Bead    *bead;
  int32    bid;
  int32    iid;
  FragType type;
  ColumnBeadIterator ci;
  int nr=0, max_nr=100;
  int32 *column_iid_list = (int32 *)safe_malloc(max_nr*sizeof(int32));

  if(!CreateColumnBeadIterator(cid, &ci)){
    fprintf(stderr, "GetReadIidsAndNumReads CreateColumnBeadIterator failed");
    assert(0);
  }
  while ( (bid = NextColumnBead(&ci)) != -1 )
    {
      char base;
 
      bead =  GetBead(beadStore,bid);
      base = *Getchar(sequenceStore,bead->soffset);
      if ( base == 'N' ) 
        continue;
      type = GetFragment(fragmentStore,bead->frag_index)->type;
      iid  = GetFragment(fragmentStore,bead->frag_index)->iid;

      if ((type == AS_READ) ||
          (type == AS_EXTR) ||
          (type == AS_TRNR))
        {
          // Filter out "duplicated" reads with the same iid
          if (!IsNewRead(iid, column_iid_list, nr))
            continue;

          column_iid_list[nr] = iid;
          nr++;
          if (nr == max_nr)
            {
              max_nr += 100;
              column_iid_list = (int32 *)safe_realloc(column_iid_list,
                                                      max_nr*sizeof(int32));
            }

          if (IsNewRead(iid, vreg->iids, vreg->nr)) {
            if (vreg->nr == vreg->max_nr) {
              int l;
              vreg->max_nr += MIN_ALLOCATED_DEPTH;
              vreg->iids = (int32 *)safe_realloc(vreg->iids, 
                                                 vreg->max_nr*sizeof(int32));
              for (l=vreg->nr; l<vreg->max_nr; l++) {
                vreg->iids[l] = -1;
              }
            }
            vreg->iids[vreg->nr] = iid;
            vreg->nr++;  
          }
        }
    }
  safe_free(column_iid_list);
}

static void
AllocateDistMatrix(VarRegion  *vreg, int init)
{
  int j, k;

  vreg->dist_matrix = (int **)safe_calloc(vreg->nr, sizeof(int *));
  for (j=0; j<vreg->nr; j++)
    {
      vreg->dist_matrix[j] = (int *)safe_calloc(vreg->nr, sizeof(int));
      for (k=0; k<vreg->nr; k++)
        vreg->dist_matrix[j][k] = init;
    }
}

static void
OutputDistMatrix(FILE *fout, VarRegion  *vreg)
{
  int j, k;

  fprintf(fout, "Distance matrix=\n");
  for (j=0; j<vreg->nr; j++)
    {
      for (k=0; k<vreg->nr; k++)
        fprintf(fout, " %d", vreg->dist_matrix[j][k]);             
      fprintf(fout, "\n");
    }
}


/*******************************************************************************
 * Function: ClusterReads
 * Purpose:  detect allele and split reads between the alleles
 *******************************************************************************
 */
static void
ClusterReads(Read *reads, int nr, Allele *alleles, int32 *na, int32 *nca, 
             int **dist_matrix)
{
  int aid, anr, row, col;

  *na = 0;

  // Iniytialize alleles
   
  // Process zero elements first                                         
  for (row=0; row<nr; row++)
    {
      for (col=row+1; col<nr; col++)
        {
          if (dist_matrix[row][col]!=0)
            continue;

          if (reads[row].allele_id < 0 && 
              reads[col].allele_id < 0)
            {
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
            }
          else if (reads[row].allele_id < 0 && 
                   reads[col].allele_id >=0)
            {
              // Already existing allele
              aid = reads[col].allele_id;
              reads[row].allele_id = aid;
              alleles[aid].weight += ROUND(reads[row].ave_qv); 
              anr = alleles[aid].num_reads;
              alleles[aid].read_ids[anr] = row;
              alleles[aid].read_iids[anr] = reads[row].iid;
              alleles[aid].num_reads++;
            }
          else if (reads[row].allele_id >=0 &&
                   reads[col].allele_id < 0)
            {
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
  for (row=0; row<nr; row++)
    {
      if (reads[row].allele_id < 0)
        {
          // New allele
          reads[row].allele_id = *na;
          alleles[*na].weight = ROUND(reads[row].ave_qv);
          alleles[*na].uglen  = reads[row].uglen;
          alleles[*na].read_ids[0] = row;
          alleles[*na].read_iids[0] = reads[row].iid;
          alleles[*na].num_reads = 1;
          alleles[*na].id = *na;
          (*na)++;
        }
    }
}
// Reverse sort by weight
static void
SortAllelesByWeight(Allele *alleles, int32 num_alleles, Read *reads)
{
  int i, j, best_id;
  Allele temp;

  for (i=0; i<num_alleles; i++)
    {
      int best_weight = alleles[i].weight;
      best_id = -1;
      for (j=i+1; j<num_alleles; j++)
        {
          if (best_weight < alleles[j].weight)
            {
              best_weight = alleles[j].weight;
              best_id = j;
            }
        }
      if (best_id >= 0)
        {
          temp       = alleles[i];
          alleles[i] = alleles[best_id];
          alleles[best_id] = temp;
        }
    }
  // Update allele_id of reads
  for (i=0; i<num_alleles; i++)
    {
      for (j=0; j<alleles[i].num_reads; j++)
        {
          int read_id = alleles[i].read_ids[j];
          reads[read_id].allele_id = i;
        }
    } 
}

// Reverse sort confirmed alleles by ungapped length
static void
SortAllelesByLength(Allele *alleles, int32 num_alleles, Read *reads)
{
  int i, j, best_id;
  Allele temp;

  for (i=0; i<num_alleles; i++)
    {
      int best_uglen = alleles[i].uglen; 
      best_id = -1;
      for (j=i+1; j<num_alleles; j++)
        {
          if (best_uglen  < alleles[j].uglen )
            {
              best_uglen  = alleles[j].uglen ;
              best_id = j;
            }
        }
      if (best_id >= 0)
        {
          temp       = alleles[i];
          alleles[i] = alleles[best_id];
          alleles[best_id] = temp;
        }
    }
  // Update allele_id of reads
  for (i=0; i<num_alleles; i++)
    {
      for (j=0; j<alleles[i].num_reads; j++)
        {
          int read_id = alleles[i].read_ids[j];
          reads[read_id].allele_id = i;
        }
    }
}

// Sort confirmed alleles according to their mapping
// between two "phased" VAR records
static void
SortAllelesByMapping(Allele *alleles, int32 nca, Read *reads,
                     int *allele_map)
{
  int i, j, k;
  Allele temp;

  for (i=0; i<nca; i++)
    {
      // j is id of the allele that should be at i-th place
      for (j=0; j<nca; j++)
        if (allele_map[j] == i) break;

      for (k=i; k<nca; k++)
        {
          if (alleles[k].id == j)
            {
              temp       = alleles[i];
              alleles[i] = alleles[k];
              alleles[k] = temp;
              break;
            }
        }
    }

  // Update allele_ids
  for (i=0; i<nca; i++)
    {
      alleles[i].id = i;
      for (j=0; j<alleles[i].num_reads; j++)
        {
          int read_id = alleles[i].read_ids[j];
          reads[read_id].allele_id = i;
        }
    }
}

static int
is_good_base(char b)
{
  if (b == '-')             return 1;
  if (b == 'a' || b == 'A') return 1;
  if (b == 'c' || b == 'C') return 1;
  if (b == 'g' || b == 'G') return 1;
  if (b == 't' || b == 'T') return 1;
  if (b == 'n' || b == 'N') return 1;
  return 0;
}

static void
UpdateScoreNumRunsOfGaps(VarRegion vreg, int prev_nr, char *prev_bases,
                         int32 *prev_iids, int get_scores)
{
  int i, j;

  // Updating count of stretches of gaps
  for (i=0; i<prev_nr; i++) {
    if (prev_bases[i] == '-')
      continue;

    for (j=0; j<vreg.nb; j++) {
      if (vreg.curr_bases[j] != '-')
        continue;

      if (prev_iids[i] == vreg.iids[j])
        {
          if (get_scores == 1)
            NumRunsOfGapsInUnitigReads++;
          else if (get_scores == 2)
            NumRunsOfGapsInContigReads++;
        }
    }
  }
}

static void
UpdateScoreNumGaps(char cbase, int get_scores)
{
  if (cbase == '-')
    {
      if (get_scores == 1)
        NumGapsInUnitigs++;
      else if (get_scores == 2)
        NumGapsInContigs++;
    }
}

static void
UpdateScores(VarRegion vreg, char *cbase, int nca)
{
  int i, j;

  for (i=0; i<nca; i++)
    {
      for (j=i+1; j<nca; j++)
        {
          if (cbase[i] != cbase[j])
            NumAAMismatches++;
        }
    }
}


static void 
GetReadsForVARRecord(Read *reads, int32 *iids, int32 nvr,
                     int beg, int end, int32 *cids)
{
  int k;

  for (k=beg; k<=end; k++)
    {
      Column *column=GetColumn(columnStore,cids[k]);
      Bead  *bead;
      int32  bid;
      int32  iid;
      FragType type;
      ColumnBeadIterator ci;
      char  base, qv;
      int i, j;
      int    nr=0, max_nr=100;
      int32 *column_iid_list = (int32 *)safe_malloc(max_nr*sizeof(int32));

      if(!CreateColumnBeadIterator(cids[k], &ci)){
        fprintf(stderr, "GetReadsForVARRecord CreateColumnBeadIterator failed");
        assert(0);
      }

      // Collect bases and usids in the coluimn
      while ( (bid = NextColumnBead(&ci)) != -1 )
        {
          bead = GetBead(beadStore,bid);
          type = GetFragment(fragmentStore,bead->frag_index)->type;
          iid  = GetFragment(fragmentStore,bead->frag_index)->iid;

          if ((type == AS_READ)   
              //              || (type == AS_EXTR) 
              //              || (type == AS_TRNR)
              )
            {
              // Filter out "duplicated" reads with the same iid
              if (!IsNewRead(iid, column_iid_list, nr))
                continue;

              column_iid_list[nr] = iid;
              nr++;
              if (nr == max_nr)
                {
                  max_nr += 100;
                  column_iid_list = (int32 *)safe_realloc(column_iid_list,
                                                          max_nr*sizeof(int32));
                }

              base = *Getchar(sequenceStore,bead->soffset);
              if (base != '-')
                {
                  qv = (int)(*Getchar(qualityStore,bead->soffset)-'0');
                }
              else // set qvs of boundary gaps to qvs of adjacent bases
                {
                  Bead *prev_bead = GetBead(beadStore, bead->prev);
                  Bead *next_bead = GetBead(beadStore, bead->next);
                  qv = 0;
                  if (prev_bead != NULL)
                    {
                      char prev_base=*Getchar(sequenceStore, prev_bead->soffset);
                      int  prev_qv  =(int)(*Getchar(qualityStore,
                                                    prev_bead->soffset)-'0');
                      if (prev_base != '-') { qv = prev_qv; } 
                      // otherwise, it stays QV_FOR_MULTI_GAP
                    }
                  if (next_bead != NULL)
                    {
                      char next_base=*Getchar(sequenceStore, next_bead->soffset);
                      int  next_qv  =(int)(*Getchar(qualityStore,
                                                    next_bead->soffset)-'0');
                      if (next_base != '-' && 
                          (qv == 0 || qv > next_qv)) 
                        qv = next_qv; 
                    }
                }
              iid  =  GetFragment(fragmentStore,bead->frag_index)->iid;
              i    =  Iid2ReadId(iid, iids, nvr);
    
              if (i < 0 || i>=nvr) {
                continue;
              }
              reads[i].bases[k-beg] = base;
              reads[i].qvs[k-beg] = qv;
              reads[i].iid = iid;
            }
        }
      safe_free(column_iid_list);
    }

  // Reset qvs of internal gaps to MIN(qv_first_gap, qv_last_gap); 
  // Compute ave_qvs
  for (k=0; k<nvr; k++)
    {
      int i, j, m = end-beg+1;
      reads[k].uglen = 0;
      reads[k].ave_qv = 0.;
      for (i=0; i<m; i++)
        {
          if (reads[k].bases[i] != '-')
            {
              reads[k].ave_qv += (double)reads[k].qvs[i];
              reads[k].uglen++;
            }
          else   // gap
            {
              int first_gap = i;
              int first_qv  = reads[k].qvs[first_gap];
              int last_gap  = i;
              int last_qv;
              int min_qv;
              if (first_qv == 0 && i>0) first_qv = reads[k].qvs[i-1];
              while (last_gap<m && reads[k].bases[last_gap] == '-')
                last_gap++;
              if (last_gap == m || reads[k].bases[last_gap] != '-')
                last_gap--;
              last_qv = reads[k].qvs[last_gap];
              if (last_qv == 0 && last_gap<m-1) 
                last_qv = reads[k].qvs[last_gap+1];
              if (first_qv != 0 && last_qv  != 0) 
                min_qv = (first_qv < last_qv) ? first_qv : last_qv;
              else if (first_qv == 0 && last_qv != 0)
                min_qv = last_qv;
              else if (first_qv != 0 && last_qv == 0)
                min_qv = first_qv;
              else // both == 0               
                min_qv = QV_FOR_MULTI_GAP;
              for (j=first_gap; j<=last_gap; j++)
                {
                  reads[k].qvs[j] = min_qv;
                  reads[k].ave_qv += (double)min_qv;
                } 
              i = last_gap;
            } 
        }
      reads[k].ave_qv /= (double)m;
    }

#if 0
  fprintf(stderr, "In GetReads: ave_qvs= ");
  for (k=0; k<nvr; k++) 
    fprintf(stderr, "%3.1f ", reads[k].ave_qv);
  fprintf(stderr, "\n");
#endif
}

static int
GetDistanceBetweenReads(char *read1, char *read2, int len)
{
  int i, j, k, uglen1=0, uglen2=0, uglen;
  int dist, gapped_dist = 0, ungapped_dist = 0;
  char *ugread1 = (char*)safe_malloc(len*sizeof(char));
  char *ugread2 = (char*)safe_malloc(len*sizeof(char));

  // Compute gapped distance
  for (k=0; k<len; k++) {
    if (read1[k] != read2[k])
      gapped_dist++;

    if (read1[k] != '-')
      {
        ugread1[uglen1] = read1[k];
        uglen1++;
      }
    if (read2[k] != '-')
      {
        ugread2[uglen2] = read2[k];
        uglen2++;
      }
  }
    
  uglen = (uglen1<uglen2) ? uglen2:uglen1;
  for (k=0; k<uglen; k++)
    {
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

static void
PopulateDistMatrix(Read *reads, int len, VarRegion  *vreg)
{
  int i, j;

  // Update the matrix
  for (i=0; i<vreg->nr; i++) {
    for (j=i; j<vreg->nr; j++) {
      vreg->dist_matrix[i][j] =
        GetDistanceBetweenReads(reads[i].bases, reads[j].bases, len);
      vreg->dist_matrix[j][i] = vreg->dist_matrix[i][j];
    }
  }
}

static int
GetTheMostDistantRead(int curr_read_id, int32 nr, int32 **dist_matrix)
{
  int i, dist_read_id = -1;
  int max_dist = -1;
  for (i=0; i<nr; i++)
    {
      if (i == curr_read_id)
        continue;

      if (max_dist < dist_matrix[curr_read_id][i])
        {
          max_dist = dist_matrix[curr_read_id][i];
          dist_read_id = i;
        }
    }
  return dist_read_id; 
}

static void
OutputReads(FILE *fout, Read *reads, int32 nr, int32 width)
{
  int i, j;
  fprintf(fout, "\nReads =\n");

  for (i=0; i<nr; i++) {
    fprintf(fout, "%d   ", reads[i].allele_id);
    for (j=0; j<width; j++)
      fprintf(fout, "%c", reads[i].bases[j]);
    fprintf(fout, "\n");
  }
  fprintf(fout, "\n\n");
}

static void
OutputAlleles(FILE *fout, VarRegion *vreg)
{
  int i, j;
  fprintf(fout,   "Outputting alleles:\n");
  fprintf(fout,   "nr= %d na= %d nca= %d\n", vreg->nr, vreg->na, vreg->nca);
  fprintf(fout,   "Num_reads= ");
  for (i=0; i<vreg->na; i++)
    {
      fprintf(fout,   "%d ", vreg->alleles[i].num_reads);
    }
  fprintf(fout,   "\n");
  fprintf(fout,   "Weights= ");
  for (i=0; i<vreg->na; i++)
    {
      fprintf(fout,   "%d ", vreg->alleles[i].weight);
    }
  fprintf(fout,   "\n");
  fprintf(fout,   "Reads= \n");
  for (i=0; i<vreg->na; i++)
    {
      fprintf(fout,   "   Allele order= %d, id= %d:\n", i, vreg->alleles[i].id);
      for (j=0; j<vreg->alleles[i].num_reads; j++)
        {
          int k, read_id = vreg->alleles[i].read_ids[j];
          int len = vreg->end-vreg->beg+1;
          fprintf(fout,   "    %d   ", read_id);
          for (k=0; k<len; k++)
            fprintf(fout,   "%c", vreg->reads[read_id].bases[k]);
          fprintf(fout,   "   %d\n", vreg->alleles[i].read_iids[j]);
        }
    }
}

static void
PopulateVARRecord(int is_phased, int32 *cids, int32 *nvars, int32 *min_len_vlist,
                  IntMultiVar **v_list, VarRegion vreg, CNS_Options *opp, int get_scores,
                  int32 *conf_read_iids)
{
  double fict_var;
  int   m;
  int   num_reported_alleles= (vreg.nca < 2) ? 2 : vreg.nca;
  char  cbase;
  char *base = (char*)safe_calloc(num_reported_alleles,sizeof(char));
  char  buf[10000];
  int   tot_num_conf_reads = 0;

  for (m=0; m<vreg.nca; m++)
    tot_num_conf_reads += vreg.alleles[m].num_reads;

  if (!(*v_list)) {
    *v_list = (IntMultiVar *)safe_malloc(*min_len_vlist*
                                         sizeof(IntMultiVar));
  }
  if (*nvars == *min_len_vlist) {
    *min_len_vlist += 10;
    *v_list = (IntMultiVar *)safe_realloc(*v_list, *min_len_vlist*
                                          sizeof(IntMultiVar));
  }
  vreg_id++;
  (*v_list)[*nvars].position.bgn = vreg.beg;
  (*v_list)[*nvars].position.end = vreg.end+1;
  (*v_list)[*nvars].num_reads = (int32)vreg.nr;
  (*v_list)[*nvars].num_conf_alleles = vreg.nca;
  (*v_list)[*nvars].anchor_size = opp->smooth_win;
  (*v_list)[*nvars].var_length = vreg.end+1-vreg.beg;
  (*v_list)[*nvars].curr_var_id = vreg_id;
  (*v_list)[*nvars].phased_var_id = is_phased ? vreg_id-1 : -1;
  (*v_list)[*nvars].weights = (char*)safe_calloc(num_reported_alleles,
                                                 (5+2)* sizeof(char));
  (*v_list)[*nvars].nr_conf_alleles = (char*)safe_calloc(num_reported_alleles,
                                                         (2+2)* sizeof(char));
  (*v_list)[*nvars].var_seq = (char*)safe_malloc(num_reported_alleles*
                                                 (vreg.end-vreg.beg+2)* sizeof(char));
  (*v_list)[*nvars].conf_read_iids= 
    (char*)safe_calloc(vreg.nr, (15+2)*sizeof(char));
  NumVARRecords++;
  {
    int al, rd;
    int32 shift   = vreg.end-vreg.beg+2;
    int distant_read_id = -42, distant_allele_id = -42;
    if (vreg.nca < 2)
      {
        distant_read_id = GetTheMostDistantRead(vreg.alleles[0].read_ids[0],
                                                vreg.nr, vreg.dist_matrix);
        distant_allele_id = vreg.reads[distant_read_id].allele_id;
      }
#if DEBUG_VAR_RECORDS
    fprintf(stderr, "VAR beg= %d end= %d\n", vreg.beg, vreg.end);
    OutputReads(stderr, vreg.reads, vreg.nr, vreg.end-vreg.beg+1);
    OutputDistMatrix(stderr, &vreg);
    OutputAlleles(stderr, &vreg);
#endif

    for (m=0; m<vreg.end-vreg.beg+1; m++)
      {
        for (al=num_reported_alleles-1; al >=0; al--)
          {
            if (al == 0 || al < vreg.nca)
              {
                int read_id = vreg.alleles[al].read_ids[0];
                base[al] = vreg.reads[read_id].bases[m];
                if (al == 0)
                  {
                    int32 cid = cids[vreg.beg+m];
                    Column *column=GetColumn(columnStore,cid);
                    Bead *call = GetBead(beadStore, column->call);
#if DEBUG_CONSENSUS_CALL
                    // Check the consistency of a consensus call
                    BaseCall(cids[vreg.beg+m], 1, &fict_var, &vreg,
                             vreg.alleles[al].id, &cbase, 0, 0, opp);
                    if (cbase != base[al])
                      {
                        fprintf(stderr, "Error setting the consensus base #%d  %c/%c\n", m, cbase, base[al]);
                        OutputReads(stderr, vreg.reads, vreg.nr, vreg.end-vreg.beg+1);
                        OutputDistMatrix(stderr, &vreg);
                        OutputAlleles(stderr, &vreg);
                        exit(-1);
                      }
#endif

                    // Set the consensus quality
                    BaseCall(cids[vreg.beg+m], 1, &fict_var, &vreg,
                             -1, &cbase, 0, 0, opp);

                    // Set the consensus base
                    Setchar(sequenceStore, call->soffset, &base[al]);

                  }
              }
            else // vreg.nca < 2 and al == 1
              {
                base[al] = vreg.reads[distant_read_id].bases[m];
              }
            (*v_list)[*nvars].var_seq[m+al*shift] = base[al];
          }
        if (get_scores > 0)
          UpdateScores(vreg, base, num_reported_alleles);
      }

    sprintf((*v_list)[*nvars].weights, "");
    sprintf((*v_list)[*nvars].nr_conf_alleles, "");
    for (al=0; al < num_reported_alleles; al++)
      {
        int    weight = vreg.alleles[al].weight;
        int    num_reads = vreg.alleles[al].num_reads;
        char *format_weight = (al < num_reported_alleles-1) ?
          "%d/" : "%d\0";
        char *format_num_reads = (al < num_reported_alleles-1) ?
          "%d/" : "%d\0";

        (*v_list)[*nvars].var_seq[-1+(al+1)*shift] =
          (al < num_reported_alleles-1) ? '/' : '\0';

        sprintf(buf, format_weight, ROUND(weight));
        (*v_list)[*nvars].weights = strcat((*v_list)[*nvars].weights, buf);

        sprintf(buf, format_num_reads, num_reads);
        (*v_list)[*nvars].nr_conf_alleles = strcat((*v_list)[*nvars].nr_conf_alleles, buf);
      }
#if 0
    fprintf(stderr, "len= %d var_seq = %s\n", vreg.end-vreg.beg+1, (*v_list)[*nvars].var_seq);

#endif
    sprintf((*v_list)[*nvars].conf_read_iids, "");
    for (rd=0; rd < tot_num_conf_reads; rd++)
      {
        char *format_iids = (rd < tot_num_conf_reads-1) ?
          "%d/" : "%d\0";
        sprintf(buf, format_iids, conf_read_iids[rd]);
        (*v_list)[*nvars].conf_read_iids =
          strcat((*v_list)[*nvars].conf_read_iids, buf);
      }

    for (al=0; al < num_reported_alleles; al++)
      {
        if ((*v_list)[*nvars].var_seq[al*shift]             == '-' &&
            (*v_list)[*nvars].var_seq[al*shift + shift - 2] == '-')
          {
            NumVARStringsWithFlankingGaps++;
          }
      }
  }
  safe_free(base);
#if DEBUG_ABACUS
  fprintf(stderr, "VARiation= %s\n", (*v_list)[*nvars].var_seq);
#endif
  (*nvars)++;
}



// Allocate memrory for reads
static void
AllocateMemoryForReads(Read **reads, int32 nr, int32 len,
                       int default_qv)
{
  int i, j;    

  *reads = (Read *)safe_malloc(nr*sizeof(Read));
  for (i=0; i<nr; i++) 
    {
      (*reads)[i].allele_id = -1;
      (*reads)[i].ave_qv = 0.;
      (*reads)[i].bases = (char *)safe_malloc(len*sizeof(char));
      (*reads)[i].qvs   = (int  *)safe_malloc(len*sizeof(int ));
      for(j=0; j<len; j++)
        {
          (*reads)[i].bases[j] = '-';
          (*reads)[i].qvs[j] = default_qv;
        }
    }
}

// Allocate memrory for alleles
static void
AllocateMemoryForAlleles(Allele **alleles, int32 nr, int32 *na)
{
  int j;

  *na = 0;
  *alleles = (Allele *)safe_calloc(nr, sizeof(Allele));
  for (j=0; j<nr; j++)
    {
      (*alleles)[j].id = -1;
      (*alleles)[j].weight = 0;
      (*alleles)[j].read_ids = (int *)safe_calloc(nr, sizeof(int));
      (*alleles)[j].read_iids = (int32 *)safe_calloc(nr, sizeof(int32));
    }
}

static int
PhaseWithPrevVreg(int32 nca, Allele *alleles, Read *reads, int32 **allele_map,
                  int32 prev_nca, int32 *prev_nca_iid, int32 prev_nca_iid_max,
                  int32 prev_ncr, int32 *prev_ncr_iid, int32 prev_ncr_iid_max)
{
  int   i, j, k, l;
  int   is_phased = 0;
  int32 **allele_matrix;
  int32  *check;
  int   num_reads = 0;

  if (prev_nca == nca && nca >= 2)
    {
      *allele_map = (int32 *)safe_calloc(nca, sizeof(int32));
      check      = (int32 *)safe_calloc(nca, sizeof(int32));
      allele_matrix    = (int32 **)safe_calloc(prev_nca, sizeof(int32 *));
      allele_matrix[0] = (int32  *)safe_calloc(prev_nca * nca, sizeof(int32));
      for (i=1; i<prev_nca; i++)
        allele_matrix[i] = allele_matrix[i-1] + nca;

      assert(prev_ncr < prev_ncr_iid_max);
      assert(prev_nca < prev_nca_iid_max);

      /* Populate the allele_matrix:
       *   columns  = confirmed alleles in the previous VAR record
       *   rows     = confirmed alleles in the current  VAR record
       *   elements = # of reads in the corresponding allele
       */
      for (i=0; i<nca; i++)   // i = allele id in current VAR record
        {
          for (j=0; j<alleles[i].num_reads; j++)
            {
              int read_id=0;
              l = 0; // allele id in the prev VAR record

              for (k=0; k<prev_ncr; k++)
                {
                  assert (l < prev_nca);  //  data in prev_nca_iid only valid up to prev_nca!

                  if (read_id == prev_nca_iid[l])  // start of a new allele
                    {
                      l++;
                      read_id = 0;
                    }
                  read_id++;

                  if (prev_ncr_iid[k] == alleles[i].read_iids[j])
                    {
                      assert(l < prev_nca);  //  matrix only allocated up to prev_nca
                      allele_matrix[l][i]++;
                      num_reads++;
                    }
                }
            }
        }

      /* Check if    alleles of previous VAR record
       * map well on alleles of current  VAR record. They will do, if:
       * - maximal element in each row in allele_matrix is greater than half
       *   the sum of all elements, and
       * - maximal elements of different rows are located in different
       *   columns
       */
      for (i=0; i<nca; i++)
        {
          (*allele_map)[i] = check[i] = -1;
        }
      for (i=0; i<nca; i++)
        {
          int sum = 0, max = -1, j_best = -1;
          for (j=0; j<prev_nca; j++)
            {
              sum += allele_matrix[i][j];
              if (max < allele_matrix[i][j])
                {
                  max = allele_matrix[i][j];
                  j_best = j;
                }
            }
          if (2*max > sum)
            {
              (*allele_map)[i] = j_best;
              check[j_best] = 1;
            }
        }

      {
        int product = 1;
        for (i=0; i< nca; i++)
          product *= (check[i]+1);
        if (product > 0)
          is_phased = 1;
      }

      /* Check if    alleles of current  VAR record
       * map well on alleles of previous VAR record. They will do, if:
       * - maximal element in each col in allele_matrix is greater than half
       *   the sum of all elements, and
       * - maximal elements of different cols are located in different
       *   rows
       */
      if (!is_phased)
        {
          for (i=0; i<nca; i++)
            {
              (*allele_map)[i] = check[i] = -1;
            }
          for (j=0; j<nca; j++) // loop through all columns
            {
              int sum = 0, max = -1, i_best = -1;
              for (i=0; i<nca; i++)
                {
                  sum += allele_matrix[i][j];
                  if (max < allele_matrix[i][j])
                    {
                      max = allele_matrix[i][j];
                      i_best = i;
                    }
                }
              if (2*max > sum)
                {
                  (*allele_map)[i_best] = j;
                  check[i_best] = 1;
                }
            }
          {
            int product = 1;
            for (i=0; i< nca; i++)
              product *= (check[i]+1);
            if (product > 0)
              is_phased = 1;
          }
        }
      safe_free(allele_matrix[0]);
      safe_free(allele_matrix);
      safe_free(check);
    } /* if (prev_nca == nca) */ 

  return is_phased;
}

static void
show_confirmed_reads(VarRegion *vreg)
{
  int j, l;
  fprintf(stderr, "Confirmed reads=\n");
  for (j=0; j<vreg->nr; j++)
    {
      if (vreg->reads[j].allele_id >= vreg->nca)
        continue;

      for (l=0; l< vreg->end-vreg->beg+1; l++)
        {
          fprintf(stderr, "%c", vreg->reads[j].bases[l]);
        }
      fprintf(stderr, " allele= %d iid= %d\n",
              vreg->reads[j].allele_id, vreg->reads[j].iid);
    }
  fprintf(stderr, "\n");
}
static void
show_reads_of_confirmed_alleles(VarRegion *vreg)
{
  int j, l;
  fprintf(stderr, "Confirmed reads=\n");
  for (j=0; j<vreg->nca; j++)
    {
      for (l=0; l< vreg->alleles[j].num_reads; l++)
        {
          fprintf(stderr, "%d  %d \n", vreg->alleles[j].read_ids[l],
                  vreg->alleles[j].read_iids[l]);
        }
    }
  fprintf(stderr, "\n");
}


static void
show_reads(VarRegion *vreg)
{
  int j, l;

  fprintf(stderr, "Num_reads= %d vreg.beg= %d vreg.end= %d\n",
          vreg->nr, vreg->beg, vreg->end);

  fprintf(stderr, "Reads=\n");
  for (j=0; j<vreg->nr; j++)
    {
      for (l=0; l< vreg->end-vreg->beg+1; l++)
        {
          fprintf(stderr, "%c", vreg->reads[j].bases[l]);
        }
      fprintf(stderr, "\n");
    }
  fprintf(stderr, "Ave_qvs= \n");
  for (j=0; j<vreg->nr; j++)
    fprintf(stderr, "%f ", vreg->reads[j].ave_qv);

  fprintf(stderr, "\n");
}


//*********************************************************************************
// Basic MultiAlignmentNode (MANode) manipulation
//*********************************************************************************
int 
RefreshMANode(int32 mid, int quality, CNS_Options *opp, int32 *nvars, 
              IntMultiVar **v_list, int make_v_list, int get_scores)
{
  // refresh columns from cid to end
  // if quality == -1, don't recall the consensus base
  int     i=0, j=0, l=0, index=0, len_manode = MIN_SIZE_OF_MANODE;
  int32   cid=0, *cids=NULL, *prev_iids=NULL;
  int     window=0, vbeg=0, vend=0, max_prev_nr=INITIAL_NR,
    prev_nr=0;
  char    cbase=0, abase=0, *prev_bases = NULL;
  char   *var_seq=NULL;
  double *varf=NULL, *svarf=NULL;
  Column *column = NULL;
  VarRegion  vreg = {0};
  char  **reads = NULL;
  MANode *ma = GetMANode(manodeStore,mid);
  int32   min_len_vlist = 10;

  // Variables used to phase VAR records
  int32  prev_nca          = 0;     // valid size of array prev_nca_iid 
  int32  prev_ncr          = 0;     // valid size of array prev_ncr_iid 
  int32  prev_nca_iid_max  = 4096;  // allocated size of arrays; err on the large side
  int32  prev_ncr_iid_max  = 4096;  //   hopefully avoiding reallocs.
  int32 *prev_nca_iid      = NULL;  // number of reads in 10 first confirmed alleles
  int32 *prev_ncr_iid      = NULL;  // iids of the first 100 reads, rev. sorted by allele

  //  Make sure that we have valid options here, we then reset the
  //  pointer to the freshly copied options, so that we can always
  //  assume opp is a valid pointer
  //
  CNS_Options  opp_private;
  if (opp == NULL) {
    opp_private.split_alleles   = CNS_OPTIONS_SPLIT_ALLELES_DEFAULT;
    opp_private.smooth_win      = CNS_OPTIONS_SMOOTH_WIN_DEFAULT;
    opp_private.max_num_alleles = CNS_OPTIONS_MAX_NUM_ALLELES;
    opp = &opp_private;
  }

  window = opp->smooth_win;

  if (nvars)
    *nvars = 0;

  assert(ma != NULL);

  if ( ma->first == -1 ) 
    return 1;

  SetDefault(&vreg);
  vreg.max_nr = MIN_ALLOCATED_DEPTH;
  vreg.iids  = (int32 *)safe_calloc(vreg.max_nr, sizeof(int32));
  vreg.curr_bases =  (char *)safe_calloc(vreg.max_nr, sizeof(char));

  varf     = (double *)safe_calloc(len_manode, sizeof(double));
  cids     =  (int32 *)safe_calloc(len_manode, sizeof(int32));
  Resetint32(ma->columns);
  cid = ma->first;
  vreg.nr = -1;

  if (get_scores > 0) {
    prev_bases = (char  *)safe_malloc(max_prev_nr*sizeof(char ));
    prev_iids  = (int32 *)safe_malloc(max_prev_nr*sizeof(int32));
  }

  // Calculate variation as a function of position in MANode. 
  while ( cid  > -1 ) 
    {
      column = GetColumn(columnStore, cid);
      assert(column != NULL);
      if ( quality != -2 ) 
        {
          if (index >= len_manode)
            {
              len_manode += MIN_SIZE_OF_MANODE;
              varf  = (double *)safe_realloc(varf,  len_manode*sizeof(double));
              cids = (int32 *)safe_realloc(cids, len_manode*sizeof(int32));
            }
          // Call consensus using all alleles
          // The goal is to detect a variation at a given position
          BaseCall(cid, quality, &(varf[index]), &vreg, -1, &cbase, 0, get_scores, 
                   opp);
          cids[index] = cid;
        }
      column->ma_index = index;
      AppendVA_int32(ma->columns, &cid);
      // sanity check
      if (index>0) {
        int32 prev= *Getint32(ma->columns, index-1);
        Column *pcol= GetColumn(columnStore, prev);
        if( prev != column->prev ||  pcol->next != column->lid)
          {
            fprintf(stderr, "RefreshMANode column relationships violated");
            assert(0);
          }
      }

      if (get_scores> 0)
        {
#if 0
          fprintf(stderr, "vreg.nb=%d vreg.curr_bases=", vreg.nb);
          for (i=0; i<vreg.nb; i++) 
            fprintf(stderr, "%c", vreg.curr_bases[i]);
          fprintf(stderr, " prev_nr=%d prev_bases=", prev_nr);
          for (i=0; i<prev_nr; i++)
            fprintf(stderr, "%c", prev_bases[i]);
          fprintf(stderr, " NumRunsOfGaps=%d \nvreg.iids= ", NumRunsOfGaps);
          for (i=0; i<vreg.nb; i++)
            fprintf(stderr, "%d ", vreg.iids[i]);
          fprintf(stderr, "\n");
          fprintf(stderr, "prev_iids= ");
          for (i=0; i<prev_nr; i++)
            fprintf(stderr, "%d ", prev_iids[i]);
          fprintf(stderr, "\n");                
#endif
          UpdateScoreNumRunsOfGaps(vreg, prev_nr, prev_bases, prev_iids, 
                                   get_scores);
          UpdateScoreNumGaps(cbase, get_scores);
          if (vreg.nb > max_prev_nr) {
            max_prev_nr =  vreg.nb;
            prev_bases = (char *)safe_realloc(prev_bases,
                                              max_prev_nr*sizeof(char));
            prev_iids  = (int32 *)safe_realloc(prev_iids,
                                               max_prev_nr*sizeof(int32));
          }
          prev_nr = vreg.nb;
          for (i=0; i<vreg.nb; i++) {
            prev_bases[i] = vreg.curr_bases[i];
            prev_iids[i]  = vreg.iids[i];
          }
        }

      cid = column->next;
      index++;
    }

  if (get_scores == 1) {
    NumColumnsInUnitigs += index;
  }
  else if (get_scores == 2) {
    NumColumnsInContigs += index;
  }

  if (opp->split_alleles == 0 || quality <= 0 || make_v_list == 0)
    {
      safe_free(vreg.curr_bases);
      safe_free(vreg.iids);
      safe_free(varf);
      safe_free(cids);
      if (get_scores > 0) {
        safe_free(prev_bases);
        safe_free(prev_iids);
      }
      return 1;
    }

  assert(nvars  != NULL);
  assert(v_list != NULL);

  // Proceed further only if accurate base calls are needed
  // Smoothen variation 
  len_manode = index -1;
  svarf= (double *)safe_calloc(len_manode, sizeof(double));
  for (i=0; i<len_manode; i++) {
    svarf[i] = varf[i];
    if (varf[i] < ZERO_MINUS)
      varf[i] = (varf[i] < -1.) ? 0. : -varf[i];
#if 0
    if (varf[i] > 0)
      fprintf(stderr, "i= %d varf= %f\n", i, varf[i]);
#endif
  }
  SmoothenVariation(svarf, len_manode, window);

  prev_nca_iid = safe_calloc(prev_nca_iid_max, sizeof(int32));
  prev_ncr_iid = safe_calloc(prev_ncr_iid_max, sizeof(int32));

  for (i=0; i<len_manode; i++) 
    { 
      if (svarf[i] == 0) {
        continue;
      }
      else 
        {
          // Process a region of variation
          int32  is_phased = 0;
          int32 *conf_read_iids = NULL;
          double fict_var;
          int32 *allele_map;

          vreg.beg = vbeg = vend = i;

          /* Set beginning of unsmoothed VAR region */
          while (vreg.beg < len_manode - 1 && 
                 DBL_EQ_DBL(varf[vreg.beg], (double)0.0))
            vreg.beg++;

          /* Set end of smoothed VAR region */
          while ((vend < len_manode-1) && (svarf[vend] > ZERO_PLUS))
            vend++;

          vreg.end = vend;

          /* Set end of unsmoothed VAR region */
          while (vreg.end >0 && varf[vreg.end] < ZERO_PLUS)
            vreg.end--;

          // Store iids of all the reads in current region
          vreg.nr = 0;
          for(l=0; l<vreg.max_nr; l++)
            vreg.iids[l] = -1;

          // Get all the read iids
          // Calculate the total number of reads, vreg.nr (corresponding to any allele)
          for (j=vreg.beg; j<=vreg.end; j++)
            GetReadIidsAndNumReads(cids[j], &vreg);

          // Allocate memrory for reads
          AllocateMemoryForReads(&vreg.reads, vreg.nr, vreg.end - vreg.beg + 1,
                                 0);

          GetReadsForVARRecord(vreg.reads, vreg.iids, vreg.nr, vreg.beg,
                               vreg.end, cids);
          // Calculate a sum of qvs for each read within a variation region
          // Populate the distance matrix

          if (SHOW_READS) show_reads(&vreg);

          AllocateDistMatrix(&vreg, -1);
          PopulateDistMatrix(vreg.reads, vreg.end-vreg.beg+1, &vreg);

          if (OUTPUT_DISTANCE_MATRIX)
            OutputDistMatrix(stderr, &vreg);

          // Allocate memory for alleles
          AllocateMemoryForAlleles(&vreg.alleles, vreg.nr, &vreg.na);

          // Populate vreg.alleles array
          // Determine the best allele and the number of reads in this allele
          ClusterReads(vreg.reads, vreg.nr, vreg.alleles, &(vreg.na),
                       &(vreg.nca), vreg.dist_matrix);

          is_phased = PhaseWithPrevVreg(vreg.nca, vreg.alleles, vreg.reads, 
                                        &allele_map,
                                        prev_nca, prev_nca_iid, prev_nca_iid_max,
                                        prev_ncr, prev_ncr_iid, prev_ncr_iid_max);

          if (is_phased)
            {
              SortAllelesByMapping(vreg.alleles, vreg.nca, vreg.reads, allele_map);
              safe_free(allele_map);
            }
          else
            SortAllelesByWeight(vreg.alleles, vreg.na, vreg.reads);

          // Create a list of iids of confirmed reads
          {
            int start_num = 0, num_conf_reads = 0;

            for (i=0; i<vreg.nca; i++)
              num_conf_reads += vreg.alleles[i].num_reads;

            conf_read_iids = (int32 *)safe_calloc(num_conf_reads, sizeof(int32));

            for (i=0; i<vreg.nca; i++)
              {
                for (j=0; j<vreg.alleles[i].num_reads; j++)
                  {
                    conf_read_iids[start_num+j] = vreg.alleles[i].read_iids[j];
                  }
                start_num += vreg.alleles[i].num_reads;
              }
          }

            
          if (SHOW_CONFIRMED_READS) show_confirmed_reads(&vreg);
          if (SHOW_CONFIRMED_READS) show_reads_of_confirmed_alleles(&vreg);

          // Update the info about the previous VAR record :
          // prev_nca, prev_nca_iid and prev_ncr_iid
          {
            int iv, jv, kv = 0;

            //  we should have space enough for these, right?  If
            //  not, count and reallocate new space -- be sure to
            //  clear it!

            prev_nca = vreg.nca;
            prev_ncr = 0;
            for (iv = 0; iv < vreg.nca; iv++)
              prev_ncr += vreg.alleles[iv].num_reads;

            if (prev_nca_iid_max < prev_nca) {
              while (prev_nca_iid_max < prev_nca)
                prev_nca_iid_max *= 2;
              safe_free(prev_nca_iid);
              prev_nca_iid = safe_calloc(prev_nca_iid_max, sizeof(int32));
            }
            if (prev_ncr_iid_max < prev_ncr) {
              while (prev_ncr_iid_max < prev_ncr)
                prev_ncr_iid_max *= 2;
              safe_free(prev_ncr_iid);
              prev_ncr_iid = safe_calloc(prev_ncr_iid_max, sizeof(int32));
            }

            memset(prev_nca_iid, 0, sizeof(int32) * prev_nca_iid_max);
            memset(prev_ncr_iid, 0, sizeof(int32) * prev_ncr_iid_max);

            for (iv = 0; iv < vreg.nca; iv++) {
              prev_nca_iid[iv] = vreg.alleles[iv].num_reads;
              for (jv = 0; jv < vreg.alleles[iv].num_reads; jv++)
                prev_ncr_iid[kv++] = vreg.alleles[iv].read_iids[jv];
            }

            if (SHOW_ALLELES) OutputAlleles(stderr, &vreg);

            /* Store variations in a v_list */
            PopulateVARRecord(is_phased, cids, nvars, &min_len_vlist, v_list,
                              vreg, opp, get_scores, conf_read_iids);
          }

          i = vend;

          for (j=0; j<vreg.nr; j++)
            {
              safe_free(vreg.dist_matrix[j]);
              safe_free(vreg.reads[j].bases);
              safe_free(vreg.reads[j].qvs);
              safe_free(vreg.alleles[j].read_ids);
              safe_free(vreg.alleles[j].read_iids);
            }
          safe_free(vreg.dist_matrix);
          safe_free(vreg.reads);
          safe_free(vreg.alleles);
          vreg.nr = 0;
          safe_free(conf_read_iids);
        }
    }
  safe_free(vreg.curr_bases);
  safe_free(vreg.iids);
  safe_free(varf);
  safe_free(svarf);
  safe_free(cids);
  if (get_scores > 0) {
    safe_free(prev_bases);
    safe_free(prev_iids);
  }
  safe_free(prev_nca_iid); 
  safe_free(prev_ncr_iid); 
  return 1;
}

int SeedMAWithFragment(int32 mid, int32 fid, int quality, 
                       CNS_Options *opp) 
{
  MANode *ma = NULL;                             
  Fragment *fragment = GetFragment(fragmentStore,fid);
  FragmentBeadIterator fi;
  int32 cid;
  int32 bid;

  ma = GetMANode(manodeStore, mid);
  assert(ma != NULL);
  assert(fragment != NULL);
  if(!CreateFragmentBeadIterator(fid,&fi)){
    fprintf(stderr, "SeedMAWithFragment CreateFragmentBeadIterator failed");
    assert(0);
  }
  bid = NextFragmentBead(&fi);
  cid = FirstColumn(mid,bid);
  while ( (bid = NextFragmentBead(&fi)) != -1 ) {
    cid = ColumnAppend(cid, bid);
  }
  fragment->manode=mid;
  {
    IntMultiVar *vl = NULL;
    int32 nv=0;
    RefreshMANode(mid, quality, opp, &nv, &vl, 0, 0); 
    safe_free(vl);
  }
  return 1;
}

int InvertTrace(int alen, int blen, Overlap *O) {
  int aend=alen+2;
  int bend=blen+2;
  int n_dels=0;
  int32 *otrace=O->trace;
  int32 *t=otrace;
  int32 *s;
  int32 tmp;
  while ( *t != 0 ) {
    n_dels++; t++;
  }
  t=otrace;
  s=t+n_dels-1;
  while (  s - t > 0 ) {
    tmp = *t;
    if ( *s < 0 ) {
      *t = - (aend + *s);
    } else {
      *t = (bend - *s);
    }
    if ( tmp < 0 ) {
      *s = - (aend + tmp);
    } else {
      *s = (bend - tmp);
    }
    t++;s--;
  }
  if ( s == t ) {
    if ( *s < 0 ) {
      *s = - (aend + *s);
    } else {
      *s = (bend - *s);
    }
  }
  tmp =O->begpos;
  O->begpos = - O->endpos;
  O->endpos = - tmp;
  return 1;
}

int * UnpackTrace(int ahang, signed char *rdelta) {
  int32 apos, bpos, idel, i, count, rdel;
  int32 delta_pos=0;
  static int32 delta[AS_READ_MAX_LEN];
  
  apos = ahang;
  bpos=0;
  while ( apos < 0 ) {
    apos++;bpos++;
  }
  if ( rdelta == NULL || rdelta[0] == 0 ) {
    delta[0] = 0;
    return delta;
  } else {
    for (idel=0;rdelta[idel]!=0;idel++) {  
      rdel = rdelta[idel];
      if ( rdel == AS_LONG_DELTA_CODE ) {
        apos+=AS_LONGEST_DELTA;
        bpos+=AS_LONGEST_DELTA;
      } else if ( rdel == AS_POLY_DELTA_CODE) {
        idel++;
        rdel = rdelta[idel];
        count = (rdel > 0)?rdel:-rdel;
        if ( rdel < 0 ) { // add gaps to a (neg. trace)
          for (i=0;i<count;i++) {
            delta[delta_pos++] = -apos;
          }
        } else {
          for (i=0;i<count;i++) { // add gaps to a (pos. trace)
            delta[delta_pos++] = bpos;
          }
        }
      } else if ( rdel < 0 ) { // align |rdel-1| positions then gap a
        for (i=0;i< -rdel - 1;i++) {
          apos++;bpos++;
        } 
        delta[delta_pos++] = -apos;
        bpos++;
      } else { // align |rdel-1| positions then gap b
        for (i=0;i< rdel - 1;i++) {
          apos++;bpos++;
        } 
        delta[delta_pos++] = bpos;
        apos++;
      }
    }
  }
  delta[delta_pos]=0;
  return delta;
}


typedef enum { 
  CNS_ALN_NONE = 'N', 
  CNS_ALN_THIN_OLAP = 'T', 
  CNS_ALN_WIDE = 'W', 
  CNS_ALN_ORIENTATION = 'O', 
  CNS_ALN_HIGH_ERATE = 'E', 
  CNS_ALN_SWAP = 'S', 
  CNS_ALN_ORIENTATION_AND_SWAP = 'B', 
  CNS_ALN_REAL_WIDE = 'X', 
  CNS_ALN_SUPER_WIDE = 'Z', 
  CNS_ALN_SEARCH_ALL = 'A', 
  CNS_ALN_EXPLICIT_DP_COMPARE = 'D',
  CNS_ALN_END_GAPS = 'G'
} CNS_AlignTrick;

typedef struct CNS_AlignParams {
  int bandBgn;
  int bandEnd;
  int maxBegGap;
  int maxEndGap;
  int opposite;
  double erate;
  double thresh;
  int minlen;
  int what;
} CNS_AlignParams;


Overlap *Compare(char *a, int alen,char *b, int blen,Overlap *(*COMPARE_FUNC)(COMPARE_ARGS), CNS_AlignParams *params) {
  Overlap *O;
  int maxbegdef=MaxBegGap;
  int maxenddef=MaxEndGap;
  if ( params->bandBgn > alen ) {
    params->bandBgn = alen;
  }
  if ( params->bandEnd > alen ) {
    params->bandEnd = alen;
  }
  if ( params->bandEnd <-blen ) {
    params->bandEnd = -blen;
  }
  if ( params->bandBgn <-blen ) {
    params->bandBgn = -blen;
  }
  if ( params->erate > 0.30 ) {
    //  Do NOT, EVER, look for an alignment with more than 30% error.
    //  This frequently occurs when AS_CNS_ERROR_RATE is set to 0.08
    //  or higher.
    params->erate = 0.30;
  }
  MaxBegGap = params->maxBegGap;
  MaxEndGap = params->maxEndGap;
  O = (*COMPARE_FUNC)(a,b, params->bandBgn, params->bandEnd,params->opposite,
                      params->erate,params->thresh,params->minlen, params->what);
  MaxBegGap = maxbegdef;
  MaxEndGap = maxenddef;
  return O;
}

void ReportOverlap(FILE *fp, Overlap *(*COMPARE_FUNC)(COMPARE_ARGS), CNS_AlignParams params,
                   int32 aiid,char atype,int32 biid,char btype,Overlap *O,int expected_hang) { 
  FILE *se=stderr;
  // This writes the basic characteristics of the overlap to both stderr, AND 
  // fp if fp is non-null
  if (O == NULL) return;

  if (fp == NULL ) {
    ReportOverlap(se,COMPARE_FUNC,params,aiid,atype,biid,btype,O,expected_hang); 
    return;
  }
  fprintf(fp,"========================================================\n");
  if ( COMPARE_FUNC == DP_Compare ) {
    fprintf(fp,"DP_Compare ");
  } else if (COMPARE_FUNC == Local_Overlap_AS_forCNS ) {
    fprintf(fp,"Local_Overlap_AS_forCNS ");
  } else {
    fprintf(fp,"An alternate aligner ");
  }

  fprintf(fp,"found overlap between %d (%c) and %d (%c) ahang: %d, bhang: %d (expected hang was %d)\n", 
          aiid,atype,biid,btype,O->begpos,O->endpos,expected_hang); 
  fprintf(fp,"Alignment params: %d %d %d %d %d %5.2f %g %d %d\n", params.bandBgn, params.bandEnd,params.maxBegGap,params.maxEndGap,params.opposite,
          params.erate,params.thresh,params.minlen, params.what);
  if (O->begpos < 0 ) fprintf(fp,"Beware, encountered unexpected negative ahang!\n"); 
  fflush(fp);
  if (fp != se ){
    ReportOverlap(stderr,COMPARE_FUNC,params,aiid,atype,biid,btype,O,expected_hang); 
  }
  return; 
} 

void PrintOverlap(FILE *fp, char *a, char *b, Overlap *O) {
  FILE *se=stderr;
  if ( O==NULL ) return;
  if (fp==NULL) {
    PrintOverlap(se,a,b,O);
    return;
  }
  Print_Overlap(fp,a,b,O);
  if ( fp != se) {
    Print_Overlap(stderr,a,b,O);
  }
  return;
}

void PrintAlarm(FILE *fp, char *msg) {
  FILE *se=stderr;
  if ( msg==NULL ) return;
  if (fp==NULL) {
    PrintAlarm(se,msg);
    return;
  }
  fprintf(fp,msg);
  if ( fp != se) {
    PrintAlarm(stderr,msg);
  }
  return;
}

void ReportTrick(FILE *fp, CNS_AlignTrick trick) {
  FILE *se=stderr;

  if (fp == NULL ) {
    ReportTrick(se,trick);
    return;
  }
  fprintf(fp,"\n========================================================");
  switch (trick) {
    case  CNS_ALN_END_GAPS:
      fprintf(fp,"\nLarge LocalAligner endgaps were allowed");
      break;
    case  CNS_ALN_HIGH_ERATE:
      fprintf(fp,"\nHigh erate was used");
      break;
    case  CNS_ALN_ORIENTATION:
      fprintf(fp,"\nOrientation reversed");
      break;
    case  CNS_ALN_THIN_OLAP:
      fprintf(fp,"\nThin overlap was used");
      break;
    case CNS_ALN_WIDE:
      fprintf(fp,"\nWide band was used");
      break;
    case CNS_ALN_SWAP:
      fprintf(fp,"\nFragments were swapped");
      break;
    case CNS_ALN_ORIENTATION_AND_SWAP: 
      fprintf(fp,"\nOrientation reversed AND fragments were swapped");
      break;
    case CNS_ALN_REAL_WIDE:
      fprintf(fp,"\nExtra-wide band was used");
      break;
    case CNS_ALN_SUPER_WIDE:
      fprintf(fp,"\nSuper-wide band was used");
      break;
    case CNS_ALN_SEARCH_ALL:
      fprintf(fp,"\nWhole search space was explored");
      break;
    case CNS_ALN_EXPLICIT_DP_COMPARE:
      fprintf(fp,"\nDP_Compare was called explicitly");
      break;
    case CNS_ALN_NONE: 
      fprintf(fp,"\nDefaults were used");
      break;
    default:
      fprintf(fp,"\nUnrecognized trick code %d",trick);
      assert(FALSE);
  }
  fprintf(fp," to capture overlap\n");
  if ( fp != se ) {
    ReportTrick(stderr,trick);
  }
  return;
}

//*********************************************************************************
// Look for the required overlap between two fragments, and return trace
//*********************************************************************************

int GetAlignmentTrace(int32 afid, int32 aoffset, int32 bfid, int32 *ahang, 
                      int32 ovl, VA_TYPE(int32) *trace, OverlapType *otype, 
                      Overlap *(*COMPARE_FUNC)(COMPARE_ARGS),int show_olap, int allow_big_endgaps) 
{
  // create a pair of dummy fragments to feed DP_Compare interface
  // aoffset is going to be used to indicate position in consensus sequence from which to start alignment
  // this will be triggered when afid == -1
  Fragment *afrag = NULL, *bfrag = NULL;
  char *a, *b;
  int32 aiid,biid;
  int32 alen,blen;
  int32 ahang_input=*ahang;
  int32 ahang_tmp;
  int32 ahang_adj=0;
  char atype,btype;
  static char cnstmpseq[2*AS_READ_MAX_LEN+1];
  static VarArrayint32 *cns_trace=NULL;
  int slip;
  int *tmp;
  Overlap *O;
  Bead *call;
  double CNS_ERATE=AS_CNS_ERROR_RATE;
  CNS_AlignTrick trick=CNS_ALN_NONE;
  int align_to_consensus=0;
  CNS_AlignParams params;

  CNS_AlignParams LOCAL_DEFAULT_PARAMS = {0,0,0,0,0,AS_CNS_ERROR_RATE,CNS_DP_THRESH,CNS_DP_MINLEN,AS_FIND_ALIGN};

  int     CNS_TIGHTSEMIBANDWIDTH = 6;
  int     CNS_LOOSESEMIBANDWIDTH = 100;
  int     CNS_DP_THIN_MINLEN     = 10;

  assert((0.0 <= AS_CNS_ERROR_RATE) && (AS_CNS_ERROR_RATE <= AS_MAX_ERROR_RATE));

  if (AS_CNS_ERROR_RATE > 0.06)
    CNS_TIGHTSEMIBANDWIDTH = 100;

  if (afid < 0  ) { 
    // copy CNS sequence into cnstmpseq, then set a to start
    int ic; 
    char callchar;
    Column *col;
    Bead *cb=GetBead(beadStore,aoffset);
    align_to_consensus=1;
    if (cns_trace==NULL) {
      cns_trace = CreateVA_int32(256);
    } else {
      ResetVA_int32(cns_trace);
    }
    for (ic=0;ic<2*AS_READ_MAX_LEN;) {
      if ( cb != NULL ) {
        col = GetColumn(columnStore,cb->column_index);
        call = GetBead(beadStore,col->call);
        //callchar=GetMaxBaseCount(&col->base_count,0);
        callchar=*Getchar(sequenceStore,call->soffset);
        if ( callchar != '-') {
          cnstmpseq[ic++] = callchar;
        } else {
          AppendVA_int32(cns_trace,&ic);
        }  
        cb = GetBead(beadStore,cb->next);
      } else {
        cnstmpseq[ic] = '\0';
        break;
      }
    }
    a = &cnstmpseq[0];
    aiid=-1;
    atype = 'M'; // "multialignment consensus"
    CNS_ERATE = 2*CNS_ERATE;
  } else {
    afrag = GetFragment(fragmentStore,afid);
    assert(afrag!=NULL);
    a = Getchar(sequenceStore,afrag->sequence);
    aiid = afrag->iid;
    atype = afrag->type;
    if ( atype == AS_UNITIG || atype == AS_CONTIG ) CNS_ERATE = 2*CNS_ERATE;
  }
  bfrag = GetFragment(fragmentStore,bfid);
  assert(bfrag!=NULL);
  biid = bfrag->iid;
  btype = bfrag->type;
  b = Getchar(sequenceStore,bfrag->sequence);
  assert(a != NULL);
  assert(b != NULL);
  alen = strlen(a); 
  blen = strlen(b);
  LOCAL_DEFAULT_PARAMS.maxBegGap = MaxBegGap;
  LOCAL_DEFAULT_PARAMS.maxEndGap = MaxEndGap;

  if ( allow_big_endgaps > 0 ) {
    LOCAL_DEFAULT_PARAMS.maxBegGap = allow_big_endgaps;
    LOCAL_DEFAULT_PARAMS.maxEndGap = allow_big_endgaps;
    PrintAlarm(stderr,"GetAlignmentTrace()-- NOTE: Looking for local alignment with large endgaps.\n");
  }
  LOCAL_DEFAULT_PARAMS.bandBgn=ahang_input-CNS_TIGHTSEMIBANDWIDTH;
  LOCAL_DEFAULT_PARAMS.bandEnd=ahang_input+CNS_TIGHTSEMIBANDWIDTH;
  if ( bfrag->type == AS_UNITIG ) LOCAL_DEFAULT_PARAMS.erate = 2*AS_CNS_ERROR_RATE;

  // Compare with the default parameters:

  params = LOCAL_DEFAULT_PARAMS;
  O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);

  if ( O == NULL ) {
    // look for potentially narrower overlap
    params.minlen=CNS_DP_THIN_MINLEN;
    O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    if (O!=NULL) trick=CNS_ALN_THIN_OLAP;
    else params = LOCAL_DEFAULT_PARAMS;
  }
  if ( O == NULL && ( strchr(a,'N') != NULL || strchr(b,'N') != NULL || bfrag->type==AS_UNITIG )) {
    // there are N's in the sequence, or they are consensus sequences, loosen the erate to compensate
    params.erate=2*AS_CNS_ERROR_RATE;
    O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    if (O!=NULL) trick=CNS_ALN_HIGH_ERATE;
    else params = LOCAL_DEFAULT_PARAMS;
  }
  if ( O == NULL  && ( ALIGNMENT_CONTEXT == AS_MERGE || bfrag->type == AS_UNITIG) ) { 
    // broaden scope out, and look for potentially narrower overlap
    // don't do this for normal unitig or contig alignments, which should be accurately placed
    params.bandBgn=ahang_input-2*CNS_LOOSESEMIBANDWIDTH;
    params.bandEnd=ahang_input+2*CNS_LOOSESEMIBANDWIDTH;
    params.erate =2*AS_CNS_ERROR_RATE;
    O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    if (O!=NULL) trick=CNS_ALN_WIDE;
    else params=LOCAL_DEFAULT_PARAMS;
  }
  if ( O == NULL && ( ALIGNMENT_CONTEXT == AS_MERGE || bfrag->type == AS_UNITIG) ) {
    // broaden even more, loosen erate
    params.bandBgn=ahang_input-3*CNS_LOOSESEMIBANDWIDTH;
    params.bandEnd=ahang_input+3*CNS_LOOSESEMIBANDWIDTH;
    params.erate=2*CNS_ERATE;
    O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    if ( O != NULL ) {
      if ( O->diffs / O->length > CNS_ERATE ) O = NULL;
    }
    if (O!=NULL) trick=CNS_ALN_REAL_WIDE;
    else params=LOCAL_DEFAULT_PARAMS;
  }
  if ( O == NULL && ( ALIGNMENT_CONTEXT == AS_MERGE || bfrag->type == AS_UNITIG) ) {
    // broaden even more, loosen erate
    params.bandBgn=ahang_input-5*CNS_LOOSESEMIBANDWIDTH;
    params.bandEnd=ahang_input+5*CNS_LOOSESEMIBANDWIDTH;
    params.erate=2*CNS_ERATE;
    O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    if ( O != NULL ) {
      if ( O->diffs / O->length > CNS_ERATE ) O = NULL;
    }
    if (O!=NULL) trick=CNS_ALN_SUPER_WIDE;
    else params=LOCAL_DEFAULT_PARAMS;
  }
  if ( O == NULL && ( ALIGNMENT_CONTEXT == AS_MERGE || bfrag->type == AS_UNITIG) ) {
    // broaden even more, loosen erate
    params.bandBgn=ahang_input-2*CNS_LOOSESEMIBANDWIDTH;
    params.bandEnd=ahang_input+2*CNS_LOOSESEMIBANDWIDTH;
    params.erate=2*CNS_ERATE;
    params.minlen=CNS_DP_THIN_MINLEN;
    O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    if ( O != NULL ) {
      if ( O->diffs / O->length > CNS_ERATE ) O = NULL;
    }
    if (O!=NULL) trick=CNS_ALN_THIN_OLAP;
    else params=LOCAL_DEFAULT_PARAMS;
  }
  if ( O == NULL && ( ALIGNMENT_CONTEXT == AS_MERGE || bfrag->type == AS_UNITIG) ) {
    // broaden even more, loosen erate
    params.bandBgn=-blen;
    params.bandEnd=alen;
    params.erate=2*CNS_ERATE;
    params.minlen=CNS_DP_THIN_MINLEN;
    O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    if ( O != NULL ) {
      if ( O->diffs / O->length > CNS_ERATE ) O = NULL;
    }
    if (O!=NULL) trick=CNS_ALN_THIN_OLAP;
    else params=LOCAL_DEFAULT_PARAMS;
  }
 
  if ( O == NULL || (O->begpos < CNS_NEG_AHANG_CUTOFF &&  ! allow_neg_hang)) { // if begpos is negative, then this isn't the intended olap
    // perhaps a poor prefix is terminating the search, or causing an alternate
    // overlap to be found 
    // try from other end
    SequenceComplement(a,NULL);
    SequenceComplement(b,NULL);
    ahang_tmp = alen - ahang_input - blen; // calculate the hang if coming from the right instead
    // note: the preceding calc may be problematic: we really would like to have the bhang, and
    // the equation above gives exactly the bhang only if the number of gaps in A and B is the
    // same -- otherwise, we are off by the number of gaps

    params.bandBgn = ahang_tmp-CNS_TIGHTSEMIBANDWIDTH;
    params.bandEnd = ahang_tmp+CNS_TIGHTSEMIBANDWIDTH;
    O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    if ( O == NULL || O->endpos  > -CNS_NEG_AHANG_CUTOFF ) {
      params.bandBgn = ahang_tmp-2*CNS_LOOSESEMIBANDWIDTH;
      params.bandEnd = ahang_tmp+2*CNS_LOOSESEMIBANDWIDTH;
      //params.bandEnd = (ALIGNMENT_CONTEXT==AS_MERGE || bfrag->type==AS_UNITIG)?ahang_tmp+2*CNS_LOOSESEMIBANDWIDTH:0,
      if ( ALIGNMENT_CONTEXT == AS_MERGE || bfrag->type == AS_UNITIG ) params.erate=2*CNS_ERATE;
      O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    }
    if ( O == NULL || O->endpos  > -CNS_NEG_AHANG_CUTOFF ) {
      //try full length of fragments, due to troubles estimating the original bhang
      params.bandBgn = -blen;
      params.bandEnd = alen;
      if ( ALIGNMENT_CONTEXT == AS_MERGE || bfrag->type == AS_UNITIG ) params.erate=2*CNS_ERATE;
      O = Compare(a,alen,b,blen,COMPARE_FUNC,&params);
    }
    if ( O == NULL || O->endpos > - CNS_NEG_AHANG_CUTOFF) {
      // here, we'll try to swap the fragments too
      params.bandBgn = -ahang_tmp-2*CNS_LOOSESEMIBANDWIDTH;
      params.bandEnd = -ahang_tmp+2*CNS_LOOSESEMIBANDWIDTH;
      O = Compare(b,blen,a,alen,COMPARE_FUNC,&params);
      if (O != NULL ) {
        int i=0;
        while (O->trace[i]!=0){
          O->trace[i++]*=-1;
        }
        O->begpos*=-1;
        O->endpos*=-1;
        trick=CNS_ALN_ORIENTATION_AND_SWAP;
      }
    } else { // the orientation alone was enough to find the right overlap
      trick=CNS_ALN_ORIENTATION;
    }
    // restore input strings and ahang for sanity if second call is necessary
    SequenceComplement(a,NULL);
    SequenceComplement(b,NULL);

    if ( O != NULL ) InvertTrace(alen,blen,O);
    else params = LOCAL_DEFAULT_PARAMS;
  }
  if ( O == NULL || ( O->begpos < CNS_NEG_AHANG_CUTOFF && ! allow_neg_hang) ) {
    // this still isn't a good overlap
    // try to see whether just swapping the fragments to see if that locates the overlap
    int tmp;
    params = LOCAL_DEFAULT_PARAMS;
    params.bandBgn = ahang_input-3*CNS_LOOSESEMIBANDWIDTH;
    params.bandEnd = ahang_input+3*CNS_LOOSESEMIBANDWIDTH;
    params.bandEnd = alen-CNS_DP_MINLEN;

    tmp=params.bandBgn;
    params.bandBgn= -params.bandEnd;
    params.bandEnd= -tmp;

    O = Compare(b,blen,a,alen,COMPARE_FUNC,&params);

    if (O != NULL ) {
      int i=0;
      while (O->trace[i]!=0){
        O->trace[i++]*=-1;
      }
      O->begpos*=-1;
      O->endpos*=-1;
      trick=CNS_ALN_SWAP;
    } else {
      params = LOCAL_DEFAULT_PARAMS;
    }
  }

  if ( O == NULL ) {
    // Here, we're convinced there is NO acceptable overlap with this COMPARE_FUNC

    fprintf(stderr,"GetAlignmentTrace()-- Could not find overlap between %d (%c) and %d (%c) estimated ahang: %d\n",
            aiid,atype,biid,btype,ahang_input);

#if 0
    fprintf(stderr,">Afrag%d\n%s\n",aiid, a);
    fprintf(stderr,">Bfrag%d\n%s\n",biid, a);
#endif

    return 0;
  }

  // from this point on, we have an Overlap

  if ( O->begpos < 0 ) {  
    // this is an undesirable situation... by construction, we anticipate all 
    // ahangs to be non-negative
    if (show_olap) {
      ReportTrick(stderr,trick);
      ReportOverlap(stderr,COMPARE_FUNC,params,aiid,atype,biid,btype,O,ahang_input);
    }

    if ( O->begpos < CNS_NEG_AHANG_CUTOFF && ! allow_neg_hang)  {
      if (show_olap) {
        if (O->begpos > -12) 
          fprintf(stderr,"GetAlignmentTrace()-- DIAGNOSTIC: would have accepted bad olap with %d bp slip\n",ahang_input-O->begpos); // diagnostic - remove soon!
        PrintOverlap(stderr, a, b, O);
        PrintAlarm(stderr,"GetAlignmentTrace()-- NOTE: Negative ahang is unacceptably large. Will not use this overlap.\n");
      }
      if ( O->begpos < -10 ) //added to get lsat 3 human partitions through
        return 0;
    }
  }
  slip = O->begpos - ahang_input;
  if (slip < 0 ) slip *=-1;
  if ( ALIGNMENT_CONTEXT != AS_MERGE && bfrag->type != AS_UNITIG && slip > CNS_TIGHTSEMIBANDWIDTH && COMPARE_FUNC == DP_Compare ) {  
    if (show_olap) {
      ReportTrick(stderr,trick);
      ReportOverlap(stderr,COMPARE_FUNC,params,aiid,atype,biid,btype,O,ahang_input);
      PrintOverlap(stderr, a, b, O);
      PrintAlarm(stderr,"GetAlignmentTrace()-- NOTE: Slip is unacceptably large. Will not use this overlap.\n");
      fprintf(stderr,"GetAlignmentTrace()-- DIAGNOSTIC: would have accepted bad olap with %d bp slip\n",slip); // diagnostic - remove soon!
    }
    //     if (O->begpos < 0 && slip < 15 ) {} //added to get last 3 human partitions through
    //        else 
    return 0;
  }

  if (show_olap) {
    if (trick != CNS_ALN_NONE) {
      // write something to the logs to show that heroic efforts were made
      ReportTrick(stderr,trick);
      ReportOverlap(stderr,COMPARE_FUNC,params,aiid,atype,biid,btype,O,ahang_input);
      PrintOverlap(stderr, a, b, O);
    }
  }

  tmp = O->trace;
  ResetVA_int32(trace);
  *otype = (O->endpos<0)?AS_CONTAINMENT:AS_DOVETAIL;
  *ahang = O->begpos + ahang_adj; // approximate ahang is replaced with found ahang
  if ( ! align_to_consensus ) {
    while ( *tmp != 0) {
      if ( *tmp < 0 ) *tmp -= ahang_adj;
      AppendVA_int32(trace,tmp);
      tmp++;
    }
  } else { 
    // here, the ungapped consensus was used, and the trace needs to be adjusted to gapped
    int num_c_gaps=GetNumint32s(cns_trace);
    int i;
    int32 *ctrace;
    int agaps=0;
    int bgaps=0;
    int new_gap_in_b;
    int ahang_gaps=0;
    int cgaps=0;
    tmp = O->trace;
    //fprintf(stderr,"revised trace: ");
    for (i=0;i<num_c_gaps;i++) {
      ctrace=Getint32(cns_trace,i);
      if ( (*ctrace + 1) < O->begpos ) {
        ahang_gaps += 1;
        cgaps++;
      } else {
        while ( *tmp != 0 ) {
          if ( *tmp < 0 ) {
            if ( (*ctrace+1)  > -*tmp ) {
              *tmp -= cgaps;
              AppendVA_int32(trace, tmp);
              //fprintf(stderr,"%d ",*tmp);
              agaps++;
              tmp++;
              continue;
            }
          } else {
            if ( *tmp + bgaps + *ahang < (*ctrace+1) + agaps ) { 
              AppendVA_int32(trace, tmp);
              //fprintf(stderr,"%d ",*tmp);
              bgaps++;
              tmp++;
              continue;
            }
          }
          new_gap_in_b= (*ctrace+1) + agaps - *ahang - ahang_gaps - bgaps;
          AppendVA_int32(trace, &new_gap_in_b);
          cgaps++;
          //fprintf(stderr,"%d ",new_gap_in_b);
          break;
        } 
      }
    }
    while ( *tmp != 0 ) {
      if ( *tmp < 0 ) {
        *tmp -= ahang_gaps+cgaps;
        AppendVA_int32(trace, tmp);
        //fprintf(stderr,"%d ",*tmp);
        tmp++;
      } else {
        AppendVA_int32(trace, tmp);
        //fprintf(stderr,"%d ",*tmp);
        tmp++;
      }
    }
  }
  return 1;
}

int MarkAsContained(int32 fid) {
  Fragment *frag = GetFragment(fragmentStore,fid);
  assert(frag != NULL);
  frag->contained = 1;
  return 1;
}

int IsContained(int32 fid) {
  Fragment *frag = GetFragment(fragmentStore,fid);
  assert(frag != NULL);
  return frag->contained;
}

int32 ApplyIMPAlignment(int32 afid, int32 bfid, int32 ahang, int32 *trace) {
  /* We assume that the bfrag frag is contained in the afrag, as a fragment
     would be to a gapped multialignment consensus 
  */
  Fragment *afrag;
  Fragment *bfrag;
  int aboffset;
  int blen;
  int bboffset;
  int apos;
  int bpos;
  Bead *abead;
  int binsert;
  afrag= GetFragment(fragmentStore,afid);
  assert(afrag != NULL);
  bfrag= GetFragment(fragmentStore,bfid);
  assert(bfrag != NULL);
  aboffset = afrag->beads;
  blen = bfrag->length;
  bboffset = bfrag->beads;
  apos = aboffset+ahang;
  bpos = 0;
  while ( (NULL != trace) && *trace != 0 ) {
    // align  (*trace-bpos) positions
    while ( *trace-bpos>0) {
      abead = GetBead(beadStore,apos++);
      AlignBead(abead->column_index, bboffset+bpos++);
    } 
    abead = GetBead(beadStore,apos++);
    binsert = AlignBead(abead->column_index, AppendGapBead(bboffset+bpos-1));
    trace++;
  }
  // now, finish up aligning the rest of b
  while ( bpos<blen ) {
    abead = GetBead(beadStore,apos++);
    AlignBead(abead->column_index, bboffset+bpos++);
  }
  bfrag->manode=afrag->manode;
  return bpos;
}

//*********************************************************************************
// Align the b fragment to the previously aligned a fragment, applying the given
// ahang and trace
//*********************************************************************************

int32 ApplyAlignment(int32 afid, int32 aoffset,int32 bfid, int32 ahang, int32 *trace) {
  // aoffset is going to be used to indicate position in consensus sequence from which to start alignment
  // this will be triggered when afid == -1
  Fragment *afrag = NULL;
  Fragment *bfrag = NULL;
  int32 aboffset, bboffset; // offsets of first beads in fragments
  int32 apos, bpos; // local offsets as alignment progresses
  int32 alen, blen;
  int32 ovl_remaining, column_appends, column_index;
  int32 first_touched_column;
  int32 last_a_aligned,last_b_aligned;
  int32 next_to_align;
  int32 binsert;
  int32 *aindex;
  Bead *abead;
  Bead *gbead;
  int32 ipx, off;
  int align_to_consensus=0;
  if ( afid < 0 ) align_to_consensus = 1;
  if ( align_to_consensus) {
    aboffset = aoffset;
    alen =0;
    { Bead *ab=GetBead(beadStore,aboffset);
      while ( ab!=NULL && alen < 2*AS_READ_MAX_LEN ){
        alen++;
        ab = GetBead(beadStore,ab->next);
      }
    }
  } else {
    afrag= GetFragment(fragmentStore,afid);
    assert(afrag != NULL);
    alen = afrag->length;
    aboffset = afrag->beads;
  }
  aindex = (int32 *)safe_malloc(alen*sizeof(int32)); 
  { Bead *ab=GetBead(beadStore,aboffset);
    int ai;
    if ( align_to_consensus ) {
      for (ai=0;ai<alen;ai++) {
        aindex[ai] = ab->boffset;
        ab = GetBead(beadStore,ab->next);
      }
    } else { 
      for (ai=0;ai<alen;ai++) {
        aindex[ai] = aboffset+ai;
      }
    }
  }
  bfrag= GetFragment(fragmentStore,bfid);
  assert(bfrag != NULL);
  blen = bfrag->length;
  bboffset = bfrag->beads;
  last_a_aligned = -1;
  last_b_aligned = -1;
  apos = MAX(ahang,0);
  bpos = 0;

  if ( ahang == alen ) { // special case where fragments abutt
    abead = GetBead(beadStore,aindex[alen-1]);
  } else {
    abead = GetBead(beadStore,aindex[apos]);
  }

  first_touched_column = abead->column_index;
  if ( ahang < 0 ) {
    gbead = GetBead(beadStore,bboffset);
    while ( bpos < -ahang ) {
      ColumnPrepend(first_touched_column,bboffset+bpos);
      bpos++;
    }
    last_b_aligned = bboffset+bpos-1;
  }

  last_a_aligned = GetBead(beadStore,aindex[apos])->prev;

  while ( (NULL != trace) && *trace != 0 ) {
    if ( *trace < 0 ) {  // gap is in afrag
      // align ( - *trace - apos ) positions
      while ( apos < (- *trace - 1)) {
        abead = GetBead(beadStore,aindex[apos]);
        AlignBead(abead->column_index, bboffset+bpos);
        last_a_aligned = aindex[apos];
        last_b_aligned = bboffset+bpos;
        apos++; bpos++;
        binsert = bboffset+bpos-1;
        while ( abead->next > -1 && (abead = GetBead(beadStore,abead->next))->boffset != aindex[apos] ) {
          // remember bead offset in case AppendGapBead messes up the pointer (MP)
          int32 off = abead->boffset;
          // insert a gap bead in b and align to    
          binsert = AppendGapBead(binsert);
          abead = GetBead(beadStore, off);
          AlignBead(abead->column_index, binsert);
          last_a_aligned = abead->boffset;
          last_b_aligned = binsert;
        }
      }
      // insert a gap column to accommodate bpos "insert"
      //   via:
      // insert a gap in afrag; insert new column seeded with that
      // then align bpos to that column
      //                           apos
      //                           | | | | |
      //               | | | | | | a a a a a 
      //               a a a a a a     
      //                   b b b b   
      //                           b b b b
      //                           bpos
      //
      //                         |
      //                         V
      //                             apos
      //                           * | | | | |
      //               | | | | | | | a a a a a 
      //               a a a a a a -    
      //                   b b b b b 
      //                           * b b b
      //                             bpos
      //              * is new column
      abead = GetBead(beadStore,aindex[apos]);
      // in case the last aligned bead in a is not apos->prev
      //   (Because gap beads were inserted, for example)
      binsert = bboffset+bpos-1;
      while ( abead->prev != last_a_aligned ) {
        binsert = AppendGapBead(binsert);
        next_to_align = (GetBead(beadStore,last_a_aligned))->next;
        AlignBead( (GetBead(beadStore,next_to_align))->column_index, binsert);
        last_a_aligned = next_to_align;
        last_b_aligned = binsert;
      }
      ColumnAppend((GetColumn(columnStore,abead->column_index))->prev,bboffset+bpos);
      abead = GetBead(beadStore,aindex[apos]);
      last_a_aligned = abead->prev;
      last_b_aligned = bboffset+bpos;
      bpos++;
    } else { // gap is in bfrag
      // align ( *trace - bpos ) positions
      while ( bpos < (*trace - 1) ) {
        abead = GetBead(beadStore,aindex[apos]);
        AlignBead(abead->column_index, bboffset+bpos);
        last_a_aligned = aindex[apos];
        last_b_aligned = bboffset+bpos;
        apos++; bpos++;
        binsert = bboffset+bpos-1;
        while ( abead->next > -1 && (abead = GetBead(beadStore,abead->next))->boffset != aindex[apos] ) {
          // remember bead offset in case AppendGapBead messes up the pointer (MP)
          int32 off = abead->boffset;
          // insert a gap bead in b and align to    
          binsert = AppendGapBead(binsert);
          abead = GetBead(beadStore, off);
          AlignBead(abead->column_index, binsert);
          last_a_aligned = abead->boffset;
          last_b_aligned = binsert;
        }
      }
      // insert a gap bead at bpos to represent bpos "delete"
      // and align the gap position with abead
      //                           apos
      //                           | | | | |
      //               | | | | | | a a a a a 
      //               a a a a a a     
      //                   b b b b   
      //                           b b b b
      //                           bpos
      //
      //                         |
      //                         V
      //                             apos
      //                             | | | |
      //               | | | | | | | a a a a 
      //               a a a a a a a    
      //                   b b b b - 
      //                             b b b b
      //                             bpos
      //              (no new column is required)

      // Jason 24-Jan-2005. Modified this code block to be resistant to the
      // stale-pointer effect when AppendGapBead realloc's the array.
      // Similar to bug fix below by MP.
      // Compare back to the original if the assert ever fails.
      ipx = aindex[apos];
      abead = GetBead(beadStore,ipx);
      off = abead->boffset;
      assert (off == ipx); 
      binsert = AppendGapBead(last_b_aligned);
      abead = GetBead(beadStore, off);
      binsert = AlignBead(abead->column_index, binsert);
      last_a_aligned = ipx;
      last_b_aligned = binsert;
      apos++;

      while ( abead->next > -1 && (abead = GetBead(beadStore,abead->next))->boffset != aindex[apos] ) {
        // remember bead offset in case AppendGapBead messes up the pointer (MP)
        off = abead->boffset;
        // insert a gap bead in b and align to    
        binsert = AppendGapBead(binsert);
        abead = GetBead(beadStore, off);
        AlignBead(abead->column_index, binsert);
        last_a_aligned = abead->boffset;
        last_b_aligned = binsert;
      }
    }
    trace++;
  }
  // remaining alignment contains no indels
  ovl_remaining  = (blen-bpos < alen-apos)?blen-bpos:alen-apos;
  while ( ovl_remaining-- > 0 ) {

    abead = GetBead(beadStore,aindex[apos]);
    AlignBead(abead->column_index, bboffset+bpos);
    last_a_aligned = abead->boffset;
    last_b_aligned = bboffset+bpos;
    apos++;bpos++; 
    binsert = bboffset+bpos-1;
    while ( abead->next > -1 && apos < alen &&
            (abead = GetBead(beadStore,abead->next))->boffset != 
            aindex[apos] ) {
      // insert a gap bead in b and align to    
      // variables needed because realloc in AppendGapBead may invalidate
      // the abead pointer (MP)
      int32 abeadIndex = abead->column_index;
      int32 abeadOffset = abead->boffset;
      binsert = AppendGapBead(binsert);
      AlignBead(abeadIndex, binsert);
      abead = GetBead(beadStore, abeadOffset);
      last_a_aligned = abeadOffset;
      last_b_aligned = binsert;
    }
  }
  column_appends = blen-bpos;
  column_index = abead->column_index;

  if ( column_appends > 0 ) {
    // First, if there are any previously aligned columns to right of abead 
    // insert gaps into b to align with these columns
    Column *pcol = GetColumn(columnStore,column_index);
    while ( pcol->next != -1 )  {
      binsert = AppendGapBead(binsert);
      column_index = pcol->next;
      AlignBead(column_index, binsert);
      pcol = GetColumn(columnStore,column_index);
    }
    // then, add on trailing (dovetail) beads from b
    while (column_appends-- > 0 ) {
      column_index = ColumnAppend(column_index,bboffset+bpos);
      bpos++;
    }
  }
  safe_free(aindex);
  bfrag->manode=afrag->manode;
  return last_b_aligned;
}

//*********************************************************************************
// Utility functions for MultiAlignmentNodes (MANode) and Fragments
//*********************************************************************************

int GetMANodeConsensus(int32 mid, VA_TYPE(char) *sequence, VA_TYPE(char) *quality) {
  ConsensusBeadIterator bi;
  Bead *bead;
  int32 bid;
  int length=GetMANodeLength(mid);
  int i=0;
  if ( sequence == NULL ) {
    sequence = CreateVA_char(length+1);
  }
  ResetVA_char(sequence);
  EnableRangeVA_char(sequence,length+1);
  if ( quality == NULL ) {
    quality = CreateVA_char(length+1);
  }
  ResetVA_char(quality);
  EnableRangeVA_char(quality,length+1);
  if(!CreateConsensusBeadIterator(mid,&bi)){
    fprintf(stderr, "GetMANodeConsensus CreateConsensusBeadIterator failed");
    assert(0);
  }
  while ( (bid = NextConsensusBead(&bi)) != -1 ) {
    bead = GetBead(beadStore,bid);
    SetVA_char(sequence,i,Getchar(sequenceStore,bead->soffset));
    SetVA_char(quality,i,Getchar(qualityStore,bead->soffset));
    i++;
  }
  return length;
}

int32 *GetFragmentDeltas(int32 fid, VA_TYPE(int32) *deltas, int length) {
  int32 delta_count = GetNumint32s(deltas);
  int32 bid;
  int32 deltas_added=0;
  FragmentBeadIterator fi;
  int32 index=0;
  if(!CreateFragmentBeadIterator(fid,&fi)){
    fprintf(stderr, "GetFragmentDeltas CreateFragmentBeadIterator failed");
    assert(0);
  }
  while ( (bid = NextFragmentBead(&fi)) != -1 && index < length) { // the index < length eliminates any endgaps from the delta list KAR, 09/19/02
    if ( *Getchar(sequenceStore,GetBead(beadStore,bid)->soffset) == '-' ) {
      Appendint32(deltas,&index);
      deltas_added++;
    } else {
      index++;
    }
  } 
  return Getint32(deltas,delta_count);
}

int GetMANodePositions(int32 mid, int mesg_n_frags, IntMultiPos *imps, int mesg_n_unitigs, IntUnitigPos *iups, VA_TYPE(int32) *deltas) {
  MANode *ma = GetMANode(manodeStore,mid);
  Fragment *fragment;
  SeqInterval position;
  IntMultiPos *fimp;
  IntUnitigPos *fump;
  int ndeletes=0;
  int odlen=0;
  int32 n_frags=0,n_unitigs=0;
  int32 i,delta_pos,prev_num_deltas;
  int hash_rc;
  assert(ma != NULL);
  if ( deltas == NULL ) {
    deltas = CreateVA_int32(gaps_in_alignment);
  } else {
    ResetVA_int32(deltas);
  }
  for (i=0;i<GetNumFragments(fragmentStore);i++) {
    fragment = GetFragment(fragmentStore,i);
    if ( fragment->deleted || fragment->manode != mid) { 
      ndeletes++; 
      continue;
    }
    position.bgn = GetColumn(columnStore,(GetBead(beadStore,fragment->beads))->column_index)->ma_index;
    position.end = GetColumn(columnStore,
                             (GetBead(beadStore,fragment->beads+fragment->length-1))->column_index)->ma_index+1;
    if ( odlen > 0 ) {
      assert (iups[0].delta_length == odlen);
    }
    if ( fragment->type  == AS_UNITIG ) {
      //fprintf(stderr,"INDEX %d, UNITIG %d, id %d ",i,n_unitigs,fragment->iid);
      assert( n_unitigs<mesg_n_unitigs ); // don't overwrite end of iup list from protomsg.
      fump = &iups[n_unitigs++];
      if(fump->ident != fragment->iid){
        fprintf(stderr, "GetMANodePositions UnitigPos id mismatch");
        assert(0);
      }
      fump->position.bgn = (fragment->complement)?position.end:position.bgn;
      fump->position.end = (fragment->complement)?position.bgn:position.end;
      fump->delta = NULL;  // just for the moment; 
      prev_num_deltas = GetNumint32s(deltas);
      GetFragmentDeltas(i,deltas,fragment->length);
      fump->delta_length = GetNumint32s(deltas)-prev_num_deltas;
      if ( n_unitigs==1 ) odlen=fump->delta_length;
      //fprintf(stderr,"Unitig %d, delta_length = %d\n", fump->ident,fump->delta_length);
    } else {
      //fprintf(stderr,"INDEX %d, READ %d, id %d ",i,n_frags,fragment->iid);
      if (ExistsInHashTable_AS (fragmentMap, fragment->iid, 0)) {
        uint64  count = LookupValueInHashTable_AS (fragmentMap, fragment->iid, 0);
        assert((count == 1) || (count == 2));
        if (count == 1) {
          // indicates that the fragment appears in the contig's f_list;
          // and this is the first time it's been seen in the fragmentStore       
          // mark that it's been seen by adding a ref to it
          //fprintf(stderr,"Fragment %d (index %d) found in the contig fragmentMap\n",i,fragment->iid);
          ReplaceInHashTable_AS(fragmentMap, fragment->iid, 0, 2, 0);
        } else {
          //fprintf(stderr,"Fragment %d (index %d) already seen in the fragmentStore\n",i,fragment->iid);
          continue;
        }
      } else {
        //fprintf(stderr,"Fragment %d not in the contig's f_list\n",fragment->iid);
        continue; // this one is not in the contig's f_list (belongs to a surrogate unitig)
      }
      assert( n_frags<mesg_n_frags ); // don't overwrite end of imp list from protomsg.
      fimp = &imps[n_frags++];
      fimp->ident = fragment->iid;
      fimp->type = fragment->type;
      fimp->position.bgn = (fragment->complement)?position.end:position.bgn;
      fimp->position.end = (fragment->complement)?position.bgn:position.end;
      fimp->delta = NULL;  // just for the moment; 
      prev_num_deltas = GetNumint32s(deltas);
      GetFragmentDeltas(i,deltas,fragment->length);
      fimp->delta_length = GetNumint32s(deltas)-prev_num_deltas;
      //fprintf(stderr,"Fragment %d, delta_length = %d\n", fimp->ident,fimp->delta_length);
    }
  }
  // now, loop through again to asign pointers to delta in imps
  // have to do this at the end to ensure that deltas isn't realloced out from under references
  delta_pos=0;  
  n_frags = 0;n_unitigs = 0;
  for (i=0;i<GetNumFragments(fragmentStore);i++) {
    fragment = GetFragment(fragmentStore,i);
    if ( fragment->deleted || fragment->manode != mid) continue;
    //fprintf(stderr,"fragment type: %c ", (char) fragment->type);
    if ( fragment->type  == AS_UNITIG ) {
      //fprintf(stderr,"UNITIG %d, (dlen %d) ",n_unitigs,iups[n_unitigs].delta_length);
      iups[n_unitigs].delta = Getint32(deltas,delta_pos);
      if (iups[n_unitigs].delta == NULL) {
        assert(iups[n_unitigs].delta_length == 0 );
      }
      delta_pos+= iups[n_unitigs].delta_length;
      n_unitigs++;
    } else {
      if (ExistsInHashTable_AS(fragmentMap, fragment->iid, 0)) {
        // all of the contig's fragments should've had their value set to 2 in previous block
        assert(2 == LookupValueInHashTable_AS(fragmentMap, fragment->iid, 0));
        // now, remove this guy from the hash_table;
        DeleteFromHashTable_AS(fragmentMap, fragment->iid, 0);
      } else  { 
        continue;
      }
      imps[n_frags].delta = Getint32(deltas,delta_pos);
      if (imps[n_frags].delta == NULL) { 
        assert(imps[n_frags].delta_length == 0 );
      }
      delta_pos+= imps[n_frags].delta_length;
      n_frags++;
    }
    //fprintf(stderr,"index %d fragment iid = %d, delta_pos = %d\n", i,fragment->iid,delta_pos);
  }
  return n_frags;
}

#if 0
//  Used only in MultiAlignUnitig
int PrintFrags(FILE *out, int accession, IntMultiPos *all_frags, int num_frags, 
               GateKeeperStore *frag_store) {
  int i,lefti,righti;
  int isread,isforward;
  int num_matches;
  GenericMesg pmesg;  
  FragMesg fmesg;
  fragRecord fsread;

  for (i=0;i<num_frags;i++) {
    isread = (all_frags[i].type == AS_READ ||
              all_frags[i].type == AS_EXTR ||
              all_frags[i].type == AS_TRNR)?1:0;
    if (all_frags[i].position.bgn<all_frags[i].position.end) {
      lefti = all_frags[i].position.bgn;
      righti = all_frags[i].position.end;
      isforward = 1;
    } else {
      righti = all_frags[i].position.bgn;
      lefti = all_frags[i].position.end;
      isforward = 0;
    }

    getFrag(gkpStore,all_frags[i].ident,&fsread,FRAG_S_ALL);

    fmesg.sequence = getFragRecordSequence(&fsread);
    fmesg.quality  = getFragRecordQuality(&fsread);

    fmesg.clear_rng.bgn = getFragRecordClearRegionBegin(&fsread, clear_range_to_use);
    fmesg.clear_rng.end = getFragRecordClearRegionEnd  (&fsread, clear_range_to_use);

    fmesg.iaccession = all_frags[i].ident;
    fmesg.type = all_frags[i].type;
    fmesg.eaccession = getFragRecordUID(&fsread);
    fmesg.action = AS_ADD;
    fmesg.source = getFragRecordSource(&fsread);

    pmesg.t = MESG_IFG;
    pmesg.m = &fmesg;
    WriteProtoMesg_AS(out,&pmesg); // write out the Fragment message
  }

  return 1;
}
#endif


void PrintAlignment(FILE *print, int32 mid, int32 from, int32 to, CNS_PrintKey what) {
  /*
    Print the columns of MANode mid from column index "from" to column index "to"
    (use 0 and -1 to print all columns)
    here's the intent for the what values;
    CNS_QUIET      = (int)'Q', // quiet,  print nothing
    CNS_STATS_ONLY = (int)'S', // print only 1-line statistic summary
    CNS_ALIGNMENT  = (int)'A', // print the multialignment, sans CNS
    CNS_CONSENSUS  = (int)'C', // print the multialignment, with CNS
    CNS_DOTS       = (int)'D', // print the multialignment, dot format
    CNS_NODOTS     = (int)'N', // print the multialignment, nodot format
    CNS_EDIT_SCORE = (int)'E'  // print the edit score column by column
  */
  MANode *ma = GetMANode(manodeStore,mid);
  int32 ma_length=GetMANodeLength(mid);
  int32 i,num_frags;
#ifdef PRINTUIDS
  int64 *fids;
#else
  int32 *fids;
#endif
  char *types;
  int32 window_start, wi;
  VA_TYPE(char) *sequenceSpace,*qualitySpace;
  char *sequence, *quality;
  char pc;
  FragmentBeadIterator *read_it;
  int32 bid;
  Bead *bead;
  Fragment *fragment;
  SeqInterval *positions;
  int dots=0;

  if(what == CNS_VIEW_UNITIG)what=CNS_DOTS;
  if (what != CNS_CONSENSUS && what != CNS_DOTS && what != CNS_NODOTS && what != CNS_VERBOSE ) return;
  if (what == CNS_DOTS) dots = 1;
  if (what == CNS_NODOTS) dots = 2;
  if (to == -1 ) {
    to = ma_length;
  }
  if(from < 0 || from > to || to > ma_length){
    fprintf(stderr, "PrintAlignment column range invalid");
    assert(0);
  }
  // now, adjust from column so that start is divisible by 100
  // (purely for convenience in viewing)
  from = ((int) from/100)*100;
  if (((int) to/100)*100 != to ) {
    to = ((int) to/100 + 1)*100;
  } else { 
    to = ((int) to/100)*100;
  }

#undef SHOW_MSA_ON_ONE_PAGE
#ifdef SHOW_MSA_ON_ONE_PAGE
  ALNPAGEWIDTH=to-from+1;
#endif

  sequenceSpace = CreateVA_char(ma_length);
  qualitySpace = CreateVA_char(ma_length);
  GetMANodeConsensus(mid,sequenceSpace,qualitySpace);
  sequence = Getchar(sequenceSpace,0);
  quality = Getchar(qualitySpace,0);
  num_frags = GetNumFragments(fragmentStore);
  read_it = (FragmentBeadIterator *) safe_calloc(num_frags,sizeof(FragmentBeadIterator));
#ifdef PRINTUIDS
  fids = (int64 *) safe_calloc(num_frags,sizeof(int64));
#else
  fids = (int32 *) safe_calloc(num_frags,sizeof(int32));
#endif
  types = (char *) safe_calloc(num_frags,sizeof(char));
  positions = (SeqInterval *) safe_calloc(num_frags,sizeof(SeqInterval));
  for (i=0;i<num_frags;i++) {
    int bgn_column;
    int end_column;
    fragment = GetFragment(fragmentStore,i); 
    if ( fragment->deleted || fragment->manode != mid) {
      fids[i] = 0;
      continue;
    }
    bgn_column = (GetBead(beadStore,fragment->beads))->column_index;
    end_column = (GetBead(beadStore,fragment->beads+fragment->length-1))->column_index;
#ifdef PRINTUIDS
    if(fragment->type==AS_READ){
      fids[i] = fragment->uid;
    } else {
      fids[i] = fragment->iid;
    }
#else
    fids[i] = fragment->iid;
#endif
    types[i] = fragment->type;
    if ( bgn_column > -1 && end_column > -1 ) {
      positions[i].bgn = GetColumn(columnStore,bgn_column)->ma_index;
      positions[i].end = GetColumn(columnStore, end_column)->ma_index+1;
    }
    NullifyFragmentBeadIterator(&read_it[i]);
  }
  window_start = from;
  fprintf(print,"\n\n================  MultiAlignment ID %d ==================\n\n",ma->iid);
  while ( window_start < to ) {

    fprintf(print,"\n%d\n%-*.*s <<< consensus\n",window_start,ALNPAGEWIDTH,ALNPAGEWIDTH,&sequence[window_start]);
    fprintf(print,"%-*.*s <<< quality\n\n",ALNPAGEWIDTH,ALNPAGEWIDTH,&quality[window_start]);
    for (i=0;i<num_frags;i++) {
      if ( fids[i] == 0 ) continue;
      for (wi = window_start;wi<window_start+ALNPAGEWIDTH;wi++) { 
        if ( IsNULLIterator(&read_it[i]) ) {
          if ( positions[i].bgn < wi && positions[i].end > wi ) {
            if(!CreateFragmentBeadIterator(i,&read_it[i])){
              fprintf(stderr, "PrintAlignment CreateFragmentBeadIterator failed");
              assert(0);
            }
            bid = NextFragmentBead(&read_it[i]);
            while ( GetColumn(columnStore,(bead=GetBead(beadStore,bid))->column_index)->ma_index < wi ) {
              bid = NextFragmentBead(&read_it[i]);
            }
            if ( bid > -1 ) {
              pc = *Getchar(sequenceStore,(GetBead(beadStore,bid))->soffset);
              if (dots == 1) {
                // check whether matches consensus, and make it a dot if so
                if (pc == sequence[wi]) pc = '.';
              }
              if (dots == 2) {
                if (pc == sequence[wi]) pc = ' ';
              }
              fprintf(print,"%c",tolower(pc));
            } 
          } else if ( positions[i].bgn ==  wi ) {
            if(!CreateFragmentBeadIterator(i,&read_it[i])){
              fprintf(stderr, "PrintAlignment CreateFragmentBeadIterator failed");
              assert(0);
            }
          } else if ( positions[i].bgn > window_start &&  positions[i].bgn < window_start+ALNPAGEWIDTH) {
            fprintf(print," ");
          } else if ( positions[i].end >= window_start &&  positions[i].end < window_start+ALNPAGEWIDTH) {
            fprintf(print," ");
          } else {
            break;
          }
        } 
        if ( ! IsNULLIterator(&read_it[i]) ) {
          bid = NextFragmentBead(&read_it[i]); 
          if ( bid > -1 ) {
            pc = *Getchar(sequenceStore,(GetBead(beadStore,bid))->soffset);
            if (dots == 1 ) {
              // check whether matches consensus, and make it a dot if so
              if (pc == sequence[wi]) pc = '.';
            }
            if (dots == 2 ) {
              // check whether matches consensus, and make it a dot if so
              if (pc == sequence[wi]) pc = ' ';
            }
            fprintf(print,"%c",tolower(pc));
          } else {
            fprintf(print," ");
            NullifyFragmentBeadIterator(&read_it[i]);
          }
        }  
#ifdef PRINTUIDS
        if ( wi == window_start+ALNPAGEWIDTH - 1 ) fprintf(print," <<< %ld (%c)\n",fids[i],types[i]);
#else
        if ( wi == window_start+ALNPAGEWIDTH - 1 ) fprintf(print," <<< %d (%c)\n",fids[i],types[i]);
#endif
      }
    }
    window_start+=ALNPAGEWIDTH;
  }
  safe_free(read_it);
  safe_free(fids);
  safe_free(types);
  safe_free(positions);
}

int RemoveNullColumn(int32 nid) {
  Column *null_column = GetColumn(columnStore,nid);
  Bead *call;
  Bead *bead;
  
  assert(null_column != NULL);
  if(GetDepth(null_column) != GetBaseCount(&null_column->base_count,'-')){
    fprintf(stderr, "RemoveNullColumn depth(null_column)!=gap basecount");
    assert(0);
  }
  call = GetBead(beadStore,null_column->call);
  while ( call->down != -1 ) {
    bead = GetBead(beadStore,call->down);
    // heal wound left by lateral removal
    if (bead->prev != -1 ) GetBead(beadStore,bead->prev)->next = bead->next;
    if (bead->next != -1 ) GetBead(beadStore,bead->next)->prev = bead->prev;
    UnAlignBead(bead->boffset);
  }
  // heal wound left by lateral removal of call
  if (call->prev != -1 ) GetBead(beadStore,call->prev)->next = call->next;
  if (call->next != -1 ) GetBead(beadStore,call->next)->prev = call->prev;
  // now, reset column pointers to bypass the removed column
  if (null_column->prev > -1 ) GetColumn(columnStore,null_column->prev)->next = null_column->next;
  if (null_column->next > -1 ) GetColumn(columnStore,null_column->next)->prev = null_column->prev;
  return 1;
}

//*********************************************************************************
// Simple sweep through the MultiAlignment columns, looking for columns
// to merge and removing null columns
//*********************************************************************************

int32 MergeRefine(int32 mid, IntMultiVar **v_list, int32 *num_vars, 
                  CNS_Options *opp, int get_scores)
{
  MANode      *ma = NULL;
  int32 cid;
  int32 nid;
  int32 removed=0;
  int32 merged;
  Column *column,*next_column;

  ma = GetMANode(manodeStore,mid);
  assert(ma != NULL);
  for (cid=ma->first;cid!=-1;){
    column = GetColumn(columnStore,cid);
    merged = MergeCompatible(cid);
    if (merged) {
      nid = column->next;
      while (nid > -1 ) {
        next_column = GetColumn(columnStore,nid);
        if ( GetDepth(next_column) == GetBaseCount(&next_column->base_count,'-') ) {
          // remove this column and try MergeColumns again
          RemoveNullColumn(nid);
          MergeCompatible(cid);
          nid = GetColumn(columnStore,cid)->next;
        } else {
          break;
        }
      }
    }
    cid = column->next;
  }
  {
    IntMultiVar *vl=NULL;
    int32 nv=0;
    int i, make_v_list=0;
 
    if (v_list && num_vars)
      make_v_list = 1;
    RefreshMANode(mid, 1, opp, &nv, &vl, make_v_list, get_scores);
    if (make_v_list && num_vars)
      {
        if (nv > 0)
          {
            *v_list = (IntMultiVar *)safe_realloc(*v_list, nv * sizeof(IntMultiVar));
            *num_vars = nv;
            for (i=0; i<nv; i++)
              {
                (*v_list)[i] = vl[i];
              }
          }
        else
          {
            safe_free(*v_list);
            *num_vars = 0;
          }
      }
    safe_free(vl);
  }
  return removed;
}

//*********************************************************************************
// Simple sweep through the MultiAlignment columns, tabulating discrepencies by QV
//*********************************************************************************

int32 AlternateDiscriminator(int32 mid, int32 *allmismatches,int32 *hqmismatches, 
                             int32 *hqsum, int32 *basecount) 
{
  MANode *ma = GetMANode(manodeStore,mid);
  int32 cid;
  int32 nid;
  int32 beadcount=0;
  int i;
  char *call;
  Bead *bead;
  ColumnBeadIterator ci;
  Column *column;
  static int qvtab[60];
  int hqtab=0;
  int alltab=0;

  assert(ma != NULL);
  for (i=0;i<60;i++) qvtab[i] = 0;
  
  for (cid=ma->first;cid!=-1;){
    column = GetColumn(columnStore,cid);
    call = Getchar(sequenceStore, GetBead(beadStore,column->call)->soffset);
    if(! CreateColumnBeadIterator(cid,&ci)){
      fprintf(stderr, "AlternateDiscriminator CreateColumnBeadIterator failed");
      assert(0);
    }
    while ( (nid = NextColumnBead(&ci)) != -1 ) {
      beadcount++;
      bead = GetBead(beadStore,nid);
      if ( *Getchar(sequenceStore,bead->soffset) != *call ) {
        // discrepancy between fragment base and consensus call
        qvtab[(int) *Getchar(qualityStore,bead->soffset)-'0']++;
        // fprintf(stderr,"%d\n",(int) *Getchar(qualityStore,bead->soffset)-'0');
      } 
    }
    cid = column->next;
  }
  
  *hqsum=0;
  for (i=0;i<60;i++) {
    if ( i>= ALT_QV_THRESH ) {
      hqtab+=qvtab[i];
      *hqsum+=i*qvtab[i];
    }
    alltab+=qvtab[i];
  }
  *allmismatches = alltab;  
  *hqmismatches = hqtab;  
  *basecount = beadcount;
  return hqtab;
}

//*********************************************************************************
// Utility functions for Abacus
//*********************************************************************************

char *GetAbacus(Abacus *a, int32 i, int32 j) {
  return (a->beads+i*(a->columns+2)+j+1);
}

void SetAbacus(Abacus *a, int32 i, int32 j, char c) 
{
  int32 offset = i*(a->columns+2)+j+1;
  if(i<0 || i>a->rows-1){
    fprintf(stderr, "i=%d j=%d a->rows=%d\n", i, j, a->rows);
    fprintf(stderr, "SetAbacus attempt to write beyond row range");
    assert(0);
  }
  if(j<0 || j>a->columns-1){
    fprintf(stderr, "i=%d j=%d a->columns=%d\n", i, j, a->columns);
    fprintf(stderr, "SetAbacus attempt to write beyond column range");
    assert(0);
  }
  a->beads[offset] = c; 
}

int ResetCalls(Abacus *a) {
  int j;
  for (j=0;j<a->columns;j++) {
    a->calls[j] = 'n';
  }
  return 1;
}

int ResetIndex(VA_TYPE(int32) * indices, int n) {
  int32 value=0;
  int i;
  for( i=0;i<n;i++) {
    SetVA_int32(indices,i,&value);
  }
  return 1;
}

Abacus *CreateAbacus(int32 mid, int32 from, int32 end) 
{
  // from,and end are ids of the first and last columns in the columnStore
  Abacus             *abacus;
  int32               columns=1, rows=0, i, j, bid, orig_columns, set_column;
  Column             *column,*last;
  ColumnBeadIterator  bi;
  Bead               *bead;
  MANode             *ma;
#define  MAX_MID_COLUMN_NUM 50
  static int                 mid_column_points[MAX_MID_COLUMN_NUM] = { 75, 150};
  Column             *mid_column[MAX_MID_COLUMN_NUM] = { NULL, NULL };
  int                 next_mid_column=0;
  int                 max_mid_columns = 0;

  //  Macaque, using overlap based trimming, needed mid_column points
  //  at rather small intervals to pass.  Even without OBT, macaque
  //  needed another point at 63 to pass.
  //
  //  This change was tested on macaque, and did not change the
  //  results (except for allowing one partition to finish....).  You
  //  can revert to the original behavior by undef'ing the following.
  //
#define EXTRA_MID_COLUMNS 1
#ifdef EXTRA_MID_COLUMNS
  if(mid_column_points[0]==75){
    for (i=0; i<MAX_MID_COLUMN_NUM; i++){
      mid_column_points[i] = i * (AS_FRAG_MIN_LEN-1) + 30;
      mid_column[i]=NULL;
    }
  }
  max_mid_columns=MAX_MID_COLUMN_NUM;

#else
  max_mid_columns=2;
#endif

  ma = GetMANode(manodeStore, mid);
  column = GetColumn(columnStore, from);

  assert(ma != NULL);
  assert(column != NULL);
  assert(abacus_indices != NULL);

  ResetIndex(abacus_indices,GetNumFragments(fragmentStore));

#if 0
  fprintf(stderr, "column_ids= %lu", from);
#endif
  // first, just determine requires number of rows and columns for Abacus
  while( column->next != end  && column->next != -1) { // this test looks subtly wrong: most of the loop could be done
    // as long as column!=NULL?  -- ALH
    columns++;

#if 0
    fprintf(stderr, ",%lu", column->next);
#endif

    if( next_mid_column<MAX_MID_COLUMN_NUM && 
        columns == mid_column_points[next_mid_column]){
      mid_column[next_mid_column] = GetColumn(columnStore,column->lid);
      next_mid_column++;
    }
    column = GetColumn(columnStore,column->next);
    // GD: this is where base calling code should be called  
  }

#if 0
  fprintf(stderr, "\n");
#endif

  if(columns>MAX_MID_COLUMN_NUM*(AS_FRAG_MIN_LEN-1)){
    fprintf(stderr,"WARNING: CreateAbacus called with such a large window that small fragments might slip through the cracks...\n");
  }

  max_mid_columns=next_mid_column;

  orig_columns = columns;
  last = column;
  column = GetColumn(columnStore, from);

  if(!CreateColumnBeadIterator(column->lid,&bi)){
    fprintf(stderr, "CreateAbacus CreateColumnBeadIterator failed");
    assert(0);
  }
  while ( (bid = NextColumnBead(&bi)) != -1 ) {
    bead = GetBead(beadStore,bid);
    rows++;
    SetVA_int32(abacus_indices,bead->frag_index,&rows);
  }

  if(!CreateColumnBeadIterator(last->lid,&bi)){
    fprintf(stderr, "CreateAbacus CreateColumnBeadIterator failed");
    assert(0);
  }
  while ( (bid = NextColumnBead(&bi)) != -1 ) {
    bead = GetBead(beadStore,bid);
    if ( *Getint32(abacus_indices,bead->frag_index) == 0 ) {
      rows++;
      SetVA_int32(abacus_indices,bead->frag_index,&rows);
    }
  }

  // a little fragment may sneak in, have to ensure the abacus has a
  // row for it.  The introduction of late- and mid column was done
  // to eliminate a problem with a degenerate alignment consistenting
  // of essentially one long poly run.  (encountered in unitig
  // 1618966 of the NOV'01 human vanilla assembly) it happened that a
  // little fragment was caught even between the mid_column and end
  // column, so it's index wasn't in the index set...  which causes a
  // "SetAbacus" beyond row range error.  putting in two mid-columns
  // will hopefully catch all fragments in the abacus range.
  //
  for (i=0; i<max_mid_columns; i++) {
    if (mid_column[i] != NULL) {
      if(!CreateColumnBeadIterator(mid_column[i]->lid,&bi)) {
        fprintf(stderr, "CreateAbacus CreateColumnBeadIterator failed");
        assert(0);
      }
      while ((bid = NextColumnBead(&bi)) != -1) {
        bead = GetBead(beadStore,bid);
        if ( *Getint32(abacus_indices,bead->frag_index) == 0 ) {
          rows++;
          SetVA_int32(abacus_indices,bead->frag_index,&rows);
        }
      }
    }
  }

  abacus = (Abacus *) safe_malloc(sizeof(Abacus));
  abacus->start_column = from;
  abacus->end_column = last->lid;
  abacus->rows = rows;
  abacus->window_width = orig_columns;
  abacus->columns = 3*orig_columns;
  abacus->shift = UNSHIFTED;
  abacus->beads = (char *) safe_calloc(rows*(abacus->columns+2),sizeof(char)); 
  abacus->calls = (char *) safe_calloc((abacus->columns),sizeof(char)); 
  // two extra gap columns, plus "null" borders

  // now, fill the center third of abacus with chars from the columns

  for (i=0;i<rows*(abacus->columns+2);i++) {
    abacus->beads[i] = 'n'; // initialize to "null" code
  } 
  columns = 0;
  while( column->lid != end  && column->lid != -1) {
    if(!CreateColumnBeadIterator(column->lid,&bi)){
      fprintf(stderr, "CreateAbacus CreateColumnBeadIterator failed");
      assert(0);
    }
    set_column = columns+orig_columns;
    while ( (bid = NextColumnBead(&bi)) != -1 ) {
      bead = GetBead(beadStore,bid);
      SetAbacus(abacus, *Getint32(abacus_indices,bead->frag_index)-1,
                set_column, *Getchar(sequenceStore,bead->soffset));
    }
    columns++;
    column = GetColumn(columnStore,column->next);
  }

  for (i=0;i<rows;i++) {
    set_column = orig_columns;
    for (j=0;j<set_column;j++) {
      SetAbacus(abacus,i,j,'-');
    }
    set_column = 2*orig_columns-1;
    for (j=set_column+1;j<abacus->columns;j++) {
      SetAbacus(abacus,i,j,'-');
    }
  }
  ResetCalls(abacus);
  return abacus;
}

void DeleteAbacus(Abacus *abacus) {
  safe_free(abacus->beads);
  safe_free(abacus->calls);
  safe_free(abacus);
}

Abacus *CloneAbacus(Abacus *abacus) {
  Abacus *clone;
  int32 rows=abacus->rows;
  int32 columns=abacus->columns;
  clone = (Abacus *) safe_malloc(sizeof(Abacus));
  clone->beads = (char *) safe_calloc(rows*(columns+2),sizeof(char)); // 
  clone->calls = (char *) safe_calloc((columns),sizeof(char)); 
  clone->rows = rows;
  clone->window_width = abacus->window_width;
  clone->columns = columns;
  clone->start_column = abacus->start_column;
  clone->end_column = abacus->end_column;
  clone->shift = abacus->shift;
  memcpy(clone->beads,abacus->beads,rows*(columns+2)*sizeof(char));
  memcpy(clone->calls,abacus->calls,columns*sizeof(char));
  return clone;
}


void ShowAbacus(Abacus *abacus) {
  int32 i;
  char form[10];
  sprintf(form,"%%%d.%ds\n",abacus->columns,abacus->columns);
#if  0
  fprintf(stderr, "form = %s\n", form);
#endif
  fprintf(stderr,"\nstart column: %d\n",abacus->start_column);
  for (i=0;i<abacus->rows;i++) {
    fprintf(stderr,form,GetAbacus(abacus,i,0));
  }
  fprintf(stderr,"\n");
  fprintf(stderr,form,abacus->calls);
}  

int32 ScoreAbacus(Abacus *abacus, int *cols)  
{ 
  // cols is the number of "good" (non-null) columns found
  // GD: This function counts the total number of bases which
  //   - are different from column's "consensus" call and
  //   - are not 'n'
  //
  BaseCount *counts;
  int score=0;
  char b;
  int i,j;
  counts = (BaseCount *) safe_calloc(abacus->columns,sizeof(BaseCount));
  memset(counts,'\0',abacus->columns*sizeof(BaseCount));
  *cols=0;

  for (i=0;i<abacus->rows;i++) {
    for (j=0;j<abacus->columns;j++) {
      b = *GetAbacus(abacus,i,j);
      if ( b == '-' ) {
        if ( j>0 && j < abacus->columns-1) {
          if ((*GetAbacus(abacus,i,j-1) == 'n')  ||
              (*GetAbacus(abacus,i,j+1) == 'n') ) 
            {
              b = 'n';
            }
        }
      }
      IncBaseCount(&counts[j],b);
    }
  }
  // now, for each column, generate the majority call
  for (j=0;j<abacus->columns;j++) {
    if ( GetBaseCount(&counts[j],'-') + GetBaseCount(&counts[j],'n') == counts[j].depth ) {
      // null (all-gap) column. Flag with an 'n' basecall
      abacus->calls[j] = 'n';
    } 
    else {
      *cols=*cols+1;
      abacus->calls[j] = GetMaxBaseCount(&counts[j],0);
      // and then tally edit score
      score += counts[j].depth - counts[j].count[BaseToInt(abacus->calls[j])] -
        counts[j].count[CNS_NALPHABET-1]; // don't count 'n's
    }
  }

  safe_free(counts);
  return score;       
}

int32 AffineScoreAbacus(Abacus *abacus)
{
  // This simply counts the number of opened gaps, to be used in tie breaker
  //   of edit scores.
  int score=0;
  char b;
  int i,j;
  int start_column, end_column;

  if (abacus->shift == LEFT_SHIFT)
    {
      start_column = 0;
      end_column   = abacus->columns/3;
    }
  else if (abacus->shift == RIGHT_SHIFT)
    {
      start_column = 2*abacus->columns/3;
      end_column   =   abacus->columns;
    }
  else //  abacus->shift == UNSHIFTED
    {
      start_column =   abacus->columns/3;
      end_column   = 2*abacus->columns/3;
    }

  for (i=0;i<abacus->rows;i++)
    {
      int in_gap=0;
      for (j=start_column;j<end_column;j++)
        {
          b = *GetAbacus(abacus,i,j);
          //      if ( abacus->calls[j] != 'n')
          //      commented out in order to make gap_score
          //      of the orig_abacus non-zero - GD
          {// don't look at null columns
            if ( b != '-' )
              {
                in_gap=0;
              }
            else
              {
                // Size of a gap does not matter, their number in a row does - GD
                if ( ! in_gap )
                  {
                    in_gap = 1;
                    score++;
                  }
              }
          }
        }
    }
  return score;
}

int MergeAbacus(Abacus *abacus, int merge_dir)
{
  // sweep through abacus from left to right
  // testing for Level 1 (neighbor) merge compatibility of each column
  // with right neighbor and merge if compatible
  //
  //  GD: this code will merge practically any
  int  i, j, k;
  int  mergeok, next_column_good, curr_column_good;
  char prev, curr, next;
  int  last_non_null = abacus->columns-1;
  int  first_non_null = 0;
  int  columns_merged = 0;

  // determine the rightmost and leftmost columns
  // not totally composed of gaps
  for (j=abacus->columns-1;j>0;j--)
    {
      int null_column = 1;
      for (i=0; i<abacus->rows; i++) {
        curr = *GetAbacus(abacus,i,j);
        if (curr != '-') null_column = 0;
      }
      if (!null_column)
        break;
      last_non_null = j;
    }
  for (j=0; j<abacus->columns;j++)
    {
      int null_column = 1;
      for (i=0; i<abacus->rows; i++) {
        curr = *GetAbacus(abacus,i,j);
        if (curr != '-') null_column = 0;
      }
      if (!null_column)
        break;
      first_non_null = j;
    }
#if DEBUG_ABACUS
  fprintf(stderr, "abacus->columns=%d first_non_null = %d last_non_null= %d\n",
          abacus->columns, first_non_null, last_non_null);
#endif
  if (merge_dir < 0)
    {
      for (j=0;j<last_non_null;j++)
        {
          int num_gaps=0, num_ns=0;
          mergeok = 1;
          next_column_good = -1;
          for (i=0;i<abacus->rows;i++)
            {
              curr = *GetAbacus(abacus,i,j);
              next = *GetAbacus(abacus,i,j+1);
              // at least in one column there should be a gap
              // or, alternatively, both should be 'n'
              if (curr != '-' && next != '-') {
                if (curr != 'n' || next != 'n') {
                  mergeok = 0;
                  break;
                }
                else
                  num_ns++;
              }
              else
                num_gaps++;

              // next column should contain at least one good base - a, c, g or t
              if (next != '-' && next != 'n') {
                next_column_good = i;
              }
            }
          //fprintf(stderr, "column= %d mergeok= %d next_column_good= %d\n", j, mergeok, next_column_good);
          if (mergeok && next_column_good >= 0 && num_gaps > num_ns) 
            {
              columns_merged++;
              for (i=0;i<abacus->rows;i++) {
                curr = *GetAbacus(abacus,i,j  );
                next = *GetAbacus(abacus,i,j+1);
                if (curr == 'n' && next == 'n')
                  {
                    continue;
                  }
                if (next != '-' && next != 'n' )
                  {
                    SetAbacus(abacus, i, j  , next);
                    SetAbacus(abacus, i, j+1, curr);
                  }
              }
              // The entire j+1-th column now contains only gaps or n's
              // Remove it by shifting all the subsequent columns
              // one position to the left
              for (i=0;i<abacus->rows;i++)
                {
                  curr = *GetAbacus(abacus,i,j  );
                  next = *GetAbacus(abacus,i,j+1);
                  if (curr == 'n' && next == 'n')
                    continue;
                  for (k=j+1; k<last_non_null; k++)
                    {
                      next= *GetAbacus(abacus,i,k+1);
                      SetAbacus(abacus, i, k, next);
                    }
                  SetAbacus(abacus, i, last_non_null, '-');
                }
              // Return to the previous coljumn to see if it can be merged again
              j--;
            }
        }
    }
  else /* merge_dir > 0 */
    {
      for (j=last_non_null-1; j>first_non_null; j--)
        {
          int num_gaps=0, num_ns=0;
          mergeok = 1;
          curr_column_good = -1;
          for (i=0;i<abacus->rows;i++)
            {
              curr = *GetAbacus(abacus,i,j);
              next = *GetAbacus(abacus,i,j+1);
              // in at least one column there should be a gap
              // or, alternatively, both should be 'n'
              if (curr != '-' && next != '-') {
                if (curr != 'n' || next != 'n') {
                  mergeok = 0;
                  break;
                }
                else
                  num_ns++; 
              }
              else
                num_gaps++;
              // current column should contain at least one good base - a, c, g or t
              if (curr != '-' && curr != 'n')
                {
                  curr_column_good = i;
                }
            }
          //fprintf(stderr, "column= %d mergeok= %d next_column_good= %d\n", j, mergeok, next_column_good);
          if (mergeok && curr_column_good >= 0 && num_gaps > num_ns) 
            {
              columns_merged++;
              for (i=0;i<abacus->rows;i++) {
                curr = *GetAbacus(abacus,i,j  );
                next = *GetAbacus(abacus,i,j+1);
                if (curr == 'n' && next == 'n')
                  {
                    continue;
                  }
                if (curr != '-' && curr != 'n' ) {
                  SetAbacus(abacus, i, j  , next);
                  SetAbacus(abacus, i, j+1, curr);
                }
              }
              // The entire j-th column contains gaps
              // Remove it by shifting all the previous columns
              // one position to the right
              for (i=0;i<abacus->rows;i++)
                {
                  curr = *GetAbacus(abacus,i,j  );
                  next = *GetAbacus(abacus,i,j+1);
                  if (curr == 'n' && next == 'n')
                    continue;
                  for (k=j; k>first_non_null; k--)
                    {
                      prev = *GetAbacus(abacus,i,k-1);
                      SetAbacus(abacus, i, k, prev);
                    }
                  SetAbacus(abacus, i, first_non_null, '-');
                }
              // Return to the next column to see if it can be merged again
              j++;
            }
        }
    }
#if 0
  fprintf(stderr, "Columns merged=%d\n", columns_merged);
#endif
  return columns_merged;
}


/* Refine the original abacus:
 *    - make sure that reads that belong to the same (confirmed) allele
 *      are aligned identically
 *    - this alignment will be dictated by the last read in the allele
 */

int32
RefineOrigAbacus(Abacus *abacus, VarRegion vreg)
{
  int32 i, j, k, l;
  int32 il;
  char c, call;
  ResetCalls(abacus);

  for (j=abacus->window_width; j<2*abacus->window_width; j++)
    {
      // Look through confirmed alleles only
      for (k=0; k<vreg.nca; k++)
        {
          // Line number corresponding to the last read in allele
          il = vreg.alleles[k].read_ids[vreg.alleles[k].num_reads-1];
          c  = *GetAbacus(abacus, il, j);
          for (l=0; l<vreg.alleles[k].num_reads-1; l++)
            {
              i = vreg.alleles[k].read_ids[l];
              SetAbacus(abacus, i, j, c);
            }
        }
    }
}

int32 
LeftShift(Abacus *abacus, VarRegion vreg, int *lcols) 
{  
  // lcols is the number of non-null columns in result
  int32 i, j, k, l, ccol, pcol;
  char c, call;
  ResetCalls(abacus);
#if 0
  fprintf(stderr, "Abacus region:\n");
  fprintf(stderr, "nr= %d na= %d nca= %d\n", vreg.nr, vreg.na, vreg.nca);
  fprintf(stderr, "Order of left-shifting alleles:\n");
  for (i=0; i<vreg.na; i++)
    {
      fprintf(stderr, "Allele %d uglen=%d weight= %d\n", i,
              vreg.alleles[i].uglen, vreg.alleles[i].weight);
      fprintf(stderr, "   Reads:\n");
      for (j=0; j<vreg.alleles[i].num_reads; j++)
        {
          int k, read_id = vreg.alleles[i].read_ids[j];
          int len = vreg.end-vreg.beg+1;
          fprintf(stderr, "    %d   ", read_id);
          for (k=0; k<len; k++)
            fprintf(stderr, "%c", vreg.reads[read_id].bases[k]);
          fprintf(stderr, "\n");
        }
    }
#endif
  for (j=abacus->window_width; j<2*abacus->window_width; j++) 
    {
      for (k=0; k<vreg.na; k++)
        {
          for (l=0; l<vreg.alleles[k].num_reads; l++) 
            {
              i = vreg.alleles[k].read_ids[l];
              c = *GetAbacus(abacus, i, j);
              ccol = j;
              if ( c != '-' ) 
                {
                  //look to the left for a suitable placement
                  // will be safe on left since abacus has 'n' border
                  while (*GetAbacus(abacus, i, ccol-1) == '-') {
                    ccol--;
                  }
                  // from ccol back up to j, look for column with matching call
                  for (pcol = ccol; pcol<j; pcol++) 
                    {
                      call = abacus->calls[pcol]; 
                      if ( call != 'n' && call != c && c != 'n') {
                        // GD: consensus in a column == '-' ? 
                        continue;
                      }     
                      if ( call == 'n') {
                        // GD: 
                        // Found the leftmost column with non-gap consensus.
                        // Now, reset its consensus "dynamically" to the 
                        // current base
                        // Potential problem: the result will generally 
                        // depend on the order in which rows 
                        // are processed
                        abacus->calls[pcol] = c;
#if 0
                        fprintf(stderr, "j= %d i= %d calls[%d]= %c\n", j, i, pcol, c);
#endif
                      } 
                      if (abacus->calls[pcol] == c || c == 'n') {
                        // swap bases in columns pcol and j of row i
                        SetAbacus(abacus, i, j, '-');
                        SetAbacus(abacus, i, pcol, c);
                        break;
                      }
                    }
                  if (*GetAbacus(abacus, i, j) != '-') 
                    abacus->calls[j] = c;
                }
            }
        }
    }
#if 0
  fprintf(stderr, "Test calls=\n");
  for (j=0;j<abacus->columns;j++)
    fprintf(stderr, "%c", abacus->calls[j]);
  fprintf(stderr, "\n");
#endif
#if DEBUG_ABACUS
  fprintf(stderr, "Abacus after LeftShift before Merge:\n");
  ShowAbacus(abacus); 
#endif
  MergeAbacus(abacus, -1);
#if 0
  fprintf(stderr, "Abacus after Merge:\n");
  ShowAbacus(abacus);
#endif
  abacus->shift = LEFT_SHIFT;
  return ScoreAbacus(abacus,lcols);
}

int32 RightShift(Abacus *abacus, VarRegion vreg, int *rcols) 
{ // rcols is the number of non-null columns in result
  int32 i, j, k, l, ccol, pcol;
  char c, call;
  ResetCalls(abacus);
#if 0
  fprintf(stderr, "Abacus region:\n");
  fprintf(stderr, "nr= %d na= %d nca= %d\n", vreg.nr, vreg.na, vreg.nca);
  fprintf(stderr, "Order of left-shifting alleles:\n");
  for (i=0; i<vreg.na; i++)
    {
      fprintf(stderr, "Allele %d uglen=%d weight= %d\n", i,
              vreg.alleles[i].uglen, vreg.alleles[i].weight);
      fprintf(stderr, "   Reads:\n");
      for (j=0; j<vreg.alleles[i].num_reads; j++)
        {
          int k, read_id = vreg.alleles[i].read_ids[j];
          int len = vreg.end-vreg.beg+1;
          fprintf(stderr, "    %d   ", read_id);
          for (k=0; k<len; k++)
            fprintf(stderr, "%c", vreg.reads[read_id].bases[k]);
          fprintf(stderr, "\n");
        }
    }
#endif

  for (j=2*abacus->window_width-1;j>abacus->window_width-1;j--) 
    {
      for (k=0; k<vreg.na; k++)
        {
          for (l=0; l<vreg.alleles[k].num_reads; l++)
            {
              i = vreg.alleles[k].read_ids[l];
              c = *GetAbacus(abacus,i,j);
              ccol = j;
              if ( c != '-' ) 
                {
                  //look to the right for a suitable placement
                  // will be safe on right since abacus has 'n' border
                  while ( *GetAbacus(abacus,i,ccol+1) == '-' ) 
                    ccol++;
                  // now, from ccol back down to j, look for column with matching call
                  for ( pcol = ccol;pcol>j;pcol--) 
                    {
                      call = abacus->calls[pcol]; 
                      if ( call != 'n' && call != c && c != 'n' ) 
                        continue;
                      if ( call == 'n') 
                        abacus->calls[pcol] = c;
                      if (abacus->calls[pcol] == c || c == 'n' ) {
                        SetAbacus(abacus,i,j,'-');
                        SetAbacus(abacus,i,pcol,c);
                        break;
                      }
                    }
                  if ( *GetAbacus(abacus,i,j) != '-' ) 
                    abacus->calls[j] = c;
                }
            }
        }
    }
#if DEBUG_ABACUS
  fprintf(stderr, "Abacus after RightShift before Merge:\n");
  ShowAbacus(abacus);
#endif
  MergeAbacus(abacus, 1);
  abacus->shift = RIGHT_SHIFT;
  return ScoreAbacus(abacus,rcols);
}

int32 MixedShift(Abacus *abacus, int *mcols, VarRegion  vreg, int lpos, int rpos,
                 char *template, int long_allele, int short_allele)
{
  // lcols is the number of non-null columns in result
  int32 i, j, k, l, ccol, pcol;
  char c, call;
  ResetCalls(abacus);
  int32 window_beg, window_end;
  int shift =0;

  if (abacus->shift == LEFT_SHIFT)
    {
      window_beg = 0;
      window_end = abacus->window_width;
    }
  else if (abacus->shift == UNSHIFTED)
    {
      window_beg = abacus->window_width;
      window_end = 2* abacus->window_width;
    }
  else
    {
      window_beg = 2*abacus->window_width;
      window_end = 3*abacus->window_width;
    }

  /* Populate calls */
  for (j=window_beg; j<window_end; j++)
    abacus->calls[j] = template[j];

  /* Perform left shift */
  for (j=window_beg;j<=MIN(window_end, lpos);j++)
    {
      for (k=0; k<vreg.na; k++)
        {
          for (l=0; l<vreg.alleles[k].num_reads; l++)
            {
              i = vreg.alleles[k].read_ids[l];
              // Only reads from short allele shouls be shifted
              if (vreg.alleles[i].id != short_allele)
                continue;

              c = *GetAbacus(abacus,i,j);
              ccol = j;
              if ( c != '-' )
                {
                  //look to the left for a suitable placement
                  // will be safe on left since abacus has 'n' border
                  while (( *GetAbacus(abacus,i,ccol-1) == '-' ) &&
                         (ccol > window_beg)) {
                    ccol--;
                  }
                  // now, from ccol back up to j, look for column with matching call
                  for ( pcol = ccol;pcol<j;pcol++) {
                    call = abacus->calls[pcol];
                    if ( call != 'n' && call != c && c != 'n') 
                      // GD: consensus in a column == '-' ?
                      continue;
                        
                    if ( call == 'n') {
                      // GD: found the leftmost column with non-gap consensus =>
                      //     reset it consensus "dynamically" to the current base
                      //     Potential problem: this code is biased  in the sense that
                      //     the result will generally depend on the order in which
                      //     reads i(or rows) are processed
                      abacus->calls[pcol] = c;
                    }
                    if (abacus->calls[pcol] == c || c == 'n') {
                      // swap bases in columns pcol and j of row i
                      SetAbacus(abacus,i,j,'-');
                      SetAbacus(abacus,i,pcol,c);
                      break;
                    }
                  }
                  if ( *GetAbacus(abacus,i,j) != '-' ) {
                    abacus->calls[j] = c;
                  }
                }
            }
        }
    }
#if 0
  fprintf(stderr, "In MixedShift: window_beg=%d lpos=%d rpos=%d  window_end=%d\n",
          window_beg, lpos, rpos, window_end);
  fprintf(stderr, "Abacus calls=\n");
  for (i=window_beg; i<window_end; i++)
    fprintf(stderr, "%c", abacus->calls[i]);
  fprintf(stderr, "\n");
#endif

  /* Perform right shift */
  for (j=window_end-1;j>(rpos>0?rpos:window_end);j--)
    {
      for (k=0; k<vreg.na; k++)
        {
          for (l=0; l<vreg.alleles[k].num_reads; l++)
            {
              i = vreg.alleles[k].read_ids[l];
              // Only reads from short allele shouls be shifted
#if 0
              fprintf(stderr, "i=%d vreg.alleles[i]=%d short_allele=%d\n", 
                      i, vreg.alleles[i], short_allele);
#endif
              if (vreg.alleles[i].id != short_allele)
                continue;

              c = *GetAbacus(abacus,i,j);
              ccol = j;
              if ( c != '-' )
                {
                  //look to the right for a suitable placement
                  // will be safe on right since abacus has 'n' border
                  while (( *GetAbacus(abacus,i,ccol+1) == '-') &&
                         (ccol+1<window_end) ) 
                    ccol++;
#if 0
                  fprintf(stderr, "ccol=%d\n", ccol);
#endif
                  // now, from ccol back down to j, look for column with matching call
                  for ( pcol = ccol;pcol>j;pcol--) {
                    call = abacus->calls[pcol];
#if 0
                    fprintf(stderr, "i=%d j=%d c=%c pcol=%d call=%d \n", 
                            i, j, c, pcol, call);
#endif
                    if ( call != 'n' && call != c && c != 'n' ) {
                      continue;
                    }
                    if ( call == 'n') {
                      abacus->calls[pcol] = c;
                    }
#if 0
                    fprintf(stderr, "abacus->calls=%c c=%c\n", abacus->calls, c);
#endif
                    if (abacus->calls[pcol] == c || c == 'n' ) {
#if 0
                      fprintf(stderr, "Swapping elements (%d, %d)=%c  and (%d, %d)='-'\n", 
                              i, j, c, i, pcol);
#endif
                      SetAbacus(abacus,i,j,'-');
                      SetAbacus(abacus,i,pcol,c);
                      break;
                    }
                  }
                  if ( *GetAbacus(abacus,i,j) != '-' ) 
                    abacus->calls[j] = c;
                }
            }
        }
    }
  // MergeAbacus(abacus, 1);
  abacus->shift = MIXED_SHIFT;
  return ScoreAbacus(abacus,mcols);
}

void GetAbacusBaseCount(Abacus *a, BaseCount *b) {
  int j;
  ResetBaseCount(b);
  for (j=0;j<a->columns;j++) {
    IncBaseCount(b,a->calls[j]); 
  }
}

int ColumnMismatch(Column *c) {
  char maxchar =  GetMaxBaseCount(&c->base_count,0);
  return c->base_count.depth - c->base_count.count[BaseToInt(maxchar)];
}

char GetBase(int s) {
  return *Getchar(sequenceStore,s);
}

int ApplyAbacus(Abacus *a, CNS_Options *opp) 
{
  Column *column; 
  int columns=0;
  int32 bid,eid,i;
  char a_entry;
  double fict_var;   // variation is a column
  Bead *bead,*exch_bead;
  VarRegion  vreg;

  SetDefault(&vreg);
  if ( a->shift == LEFT_SHIFT) 
    {
      column = GetColumn(columnStore,a->start_column);
      assert(column != NULL);
      while (columns<a->window_width) 
        {
          char base;
          bid = GetBead(beadStore,column->call)->down;
          while ( bid != -1 ) 
            {
              // Update all beads in a given column
              bead = GetBead(beadStore,bid);
              i =  *Getint32(abacus_indices,bead->frag_index) - 1;
              a_entry = *GetAbacus(a,i,columns);
              if ( a_entry == 'n') 
                {
                  exch_bead = GetBead(beadStore,bead->up);
                  //fprintf(stderr,"Unaligning trailing gaps from %d.\n",bid);
                  UnAlignTrailingGapBeads(bid);
                } 
              else if ( a_entry != *Getchar(sequenceStore,bead->soffset)) 
                {
                  //  Look for matching bead in frag and exchange
                  exch_bead = GetBead(beadStore,bead->boffset);
                  if ( NULL == exch_bead ) {
                    //fprintf(stderr,"Uh-oh... out of beads in fragment. (LEFT_SHIFT)\n");
                    eid = AppendGapBead(bead->boffset);
                    //fprintf(stderr,"Adding gapbead %d\n",eid);
                    AlignBead(GetColumn(columnStore,bead->column_index)->next,eid);
                    exch_bead = GetBead(beadStore,eid);
                  }
                  while (  a_entry != *Getchar(sequenceStore,exch_bead->soffset)) 
                    {
                      if (exch_bead->next == -1 ) {
                        //fprintf(stderr,"Uh-oh... out of beads in fragment. (LEFT_SHIFT)\n");
                        eid = AppendGapBead(exch_bead->boffset);
                        //fprintf(stderr,"Adding gapbead %d\n",eid);
                        AlignBead(GetColumn(columnStore,exch_bead->column_index)->next,eid);
                      } else if (exch_bead->column_index == a->end_column) {
                        //fprintf(stderr,"Uh-oh... out of beads in window. (LEFT_SHIFT)\n");
                        eid = AppendGapBead(exch_bead->boffset);
                        //fprintf(stderr,"Adding gapbead %d\n",eid);
                        // ColumnAppend(exch_bead->column_index,eid);
                        { // mods (ALH) to handle reallocation of columnStore
                          int32 off   = 0;
                          int curridx = column->lid;

                          // remember bead offset in case realloc() moves the
                          // memory and exch_bead becomes a stale pointer
                          off = exch_bead->boffset; 
                          ColumnAppend(exch_bead->column_index,eid);
                          exch_bead = GetBead(beadStore, off);

                          column=GetColumn(columnStore,curridx);
                        } 
                      }
                      exch_bead = GetBead(beadStore,exch_bead->next);
                    }
                  /* 
                     fprintf(stderr,"LeftShifting bead %d (%c) with bead %d (%c).\n",
                     bid, *Getchar(sequenceStore,GetBead(beadStore,bid)->soffset),
                     exch_bead->boffset, *Getchar(sequenceStore,GetBead(beadStore,exch_bead->boffset)->soffset));
                  */ 
                  LeftEndShiftBead(bid,exch_bead->boffset);
                } else {
                exch_bead = bead; // no exchange necessary;
              }
              bid = exch_bead->down;
              /*
                fprintf(stderr,"New bid is %d (%c), from %d down\n",
                bid, (bid > -1)?*Getchar(sequenceStore,GetBead(beadStore,bid)->soffset):'n',
                exch_bead->boffset);
              */
            }
          BaseCall(column->lid, 1, &fict_var, &vreg, -1, &base, 0, 0, opp);
          column = GetColumn(columnStore,column->next);
          columns++;
        } 
    } 
  else if ( a->shift == RIGHT_SHIFT)
    {
      column = GetColumn(columnStore,a->end_column);
      assert(column != NULL);
      while (columns<a->window_width) {
        char base;
        bid = GetBead(beadStore,column->call)->down;
        while ( bid != -1 ) {
          bead = GetBead(beadStore,bid);
          i =  *Getint32(abacus_indices,bead->frag_index) - 1;
          a_entry = *GetAbacus(a,i,a->columns-columns-1);
          if ( a_entry == 'n' ) {
            exch_bead = GetBead(beadStore,bead->up);
            //fprintf(stderr,"Unaligning trailing gaps from %d.\n",bid);
            UnAlignTrailingGapBeads(bid);
          } else if ( a_entry != *Getchar(sequenceStore,bead->soffset)) {
            //  Look for matching bead in frag and exchange
            exch_bead = GetBead(beadStore,bead->boffset);
            if ( NULL == exch_bead ) {
              //fprintf(stderr,"Uh-oh... out of beads in fragment. (RIGHT_SHIFT)\n");
              eid = PrependGapBead(bead->boffset);
              //fprintf(stderr,"Adding gapbead %d\n",eid);
	     
              AlignBead(GetColumn(columnStore,bead->column_index)->prev,eid);
              exch_bead = GetBead(beadStore,eid);
            }

            while (  a_entry != *Getchar(sequenceStore,exch_bead->soffset)) {
              int exchbeadprev = exch_bead->prev;

              if (exch_bead->prev == -1 ) {
                //fprintf(stderr,"Uh-oh... out of beads in fragment. (RIGHT_SHIFT)\n");
                exchbeadprev = eid = PrependGapBead(exch_bead->boffset);
                //fprintf(stderr,"Adding gapbead %d\n",eid);
                AlignBead(GetColumn(columnStore,exch_bead->column_index)->prev,eid);
              } else if (exch_bead->column_index == a->start_column) {
                //fprintf(stderr,"Uh-oh... out of beads in window. (RIGHT_SHIFT)\n");
                exchbeadprev = eid = AppendGapBead(exch_bead->prev);
                //fprintf(stderr,"Adding gapbead %d\n",eid);

                {// ALH's change to fix reallocation of column store
                  int curridx = column->lid;
                  ColumnAppend(GetColumn(columnStore,exch_bead->column_index)->prev,eid);	         
                  column = GetColumn(columnStore, curridx);
                }
              }
              exch_bead = GetBead(beadStore, exchbeadprev);
            }
            /*
              fprintf(stderr,"RightShifting bead %d (%c) with bead %d (%c).\n",
              exch_bead->boffset, *Getchar(sequenceStore,GetBead(beadStore,exch_bead->boffset)->soffset),
              bid, *Getchar(sequenceStore,GetBead(beadStore,bid)->soffset));
            */
            RightEndShiftBead(exch_bead->boffset,bid);
          } else {
            exch_bead = bead; // no exchange necessary;
          }
          bid = exch_bead->down;
          /*
            fprintf(stderr,"New bid is %d (%c), from %d down\n",
            bid, (bid>-1)?*Getchar(sequenceStore,GetBead(beadStore,bid)->soffset):'n',
            exch_bead->boffset);
          */
        }
        BaseCall(column->lid, 1, &fict_var, &vreg, -1, &base, 0, 0, opp);
        column = GetColumn(columnStore,column->prev);
        columns++;
      } 
    }
  return 1;
}

int IdentifyWindow(Column **start_column, int *stab_bgn, CNS_RefineLevel level) 
{
  Column *stab;
  Column *pre_start;
  int win_length=1;
  int rc=0;
  int gap_count=0;
  char poly;
  *stab_bgn = (*start_column)->next; 
  stab = GetColumn(columnStore,*stab_bgn);
  switch (level) {
    case CNS_SMOOTH: 
      // in this case, we just look for a string of gaps in the consensus sequence
      if ( GetBase( GetBead(beadStore,(*start_column)->call)->soffset ) != '-' ) break;
      // here, there's a '-' in the consensus sequence, see if it expands 
      while( GetBase( GetBead(beadStore,stab->call)->soffset) == '-' )  {
        // move stab column ahead
        if ( stab->next != -1 ) {
          *stab_bgn = stab->next;
          stab = GetColumn(columnStore,*stab_bgn);
          win_length++;
        } else {
          break;
        }
      }
      if ( win_length > 1 ) rc = win_length;
      break;
    case CNS_POLYX:
      // here, we're looking for a string of the same character
      gap_count=GetColumnBaseCount(*start_column,'-');
      poly =  GetBase(GetBead(beadStore,(*start_column)->call)->soffset);
      if ( poly != '-' ) {
        char cb;
       
        while( (cb = GetBase(GetBead(beadStore,stab->call)->soffset)) == poly || cb == '-' )  {
          // move stab column ahead
          if ( stab->next != -1 ) {
            *stab_bgn = stab->next;
            gap_count+=GetColumnBaseCount(stab,'-');
            stab = GetColumn(columnStore,*stab_bgn);
            win_length++;
          } else {
            break;
          }
        }
        // capture trailing gap-called columns
        if ( win_length > 2 ) {
          while( GetBase(GetBead(beadStore,stab->call)->soffset) == '-' )  {
            if ( GetMaxBaseCount(&stab->base_count,1) != poly ) break;
            if ( stab->next != -1 ) {
              *stab_bgn = stab->next;
              gap_count+=GetColumnBaseCount(stab,'-');
              stab = GetColumn(columnStore,*stab_bgn);
              win_length++;
            } else {
              break;
            }
          }
          // now that a poly run with trailing gaps is established, look for leading gaps
          pre_start = *start_column;
          while ( pre_start->prev != -1 ) {
            char cb;
            pre_start = GetColumn(columnStore,pre_start->prev);
            if ( (cb = GetBase(GetBead(beadStore,pre_start->call)->soffset)) != '-' && cb != poly ) break;
            *start_column = pre_start; 
            gap_count+=GetColumnBaseCount(pre_start,'-');
            win_length++;
          }
        } else {
          break;
        }
      }
      if ( (*start_column)->prev != -1 && win_length > 2 && gap_count > 0) {
        //fprintf(stderr,"POLYX candidate (%c) at column %d, width %d, gapcount %d\n", poly,(*start_column)->ma_index,win_length,gap_count);
        rc = win_length;
      }
      break;
    case CNS_INDEL:
      /*
        in this case, we look for a string mismatches, indicating a poor alignment region
        which might benefit from Abacus refinement
        heuristics: 
        > stable border on either side of window of width:  STABWIDTH
        > fewer than STABMISMATCH in stable border
              
        _              __              ___
        SSSSS SSSSS    SSSSS .SSSS+    SSSSS  .SSSS+
        SSSSS SSSSS    SSSSS .SSSS+    SSSSS  .SSSS+
        SSSSS SSSSS => SSSSS .SSSS+ => SSSSS  .SSSS+
        SSSSS SSSSS    SSSSS .SSSS+    SSSSS  .SSSS+
        SSSSS_SSSSS    SSSSS_.SSSS+    SSSSS__.SSSS+
        |               \               \ 
        |\_______________\_______________\______ growing 'gappy' window 
        start_column
      */
      {
        int cum_mm=0;
        int stab_mm=0;
        int stab_gaps=0;
        int stab_width=0;
        int stab_bases=0;
        Column *stab_end;

        cum_mm = ColumnMismatch(*start_column);
        if ( cum_mm > 0 && GetColumnBaseCount(*start_column,'-') > 0) {
          stab = *start_column;
          stab = GetColumn(columnStore,(*start_column)->next);
          stab_end = stab;
          while ( stab_end->next != -1 && stab_width < STABWIDTH) {
            stab_mm+=ColumnMismatch(stab_end);
            stab_gaps+=GetColumnBaseCount(stab_end,'-');
            stab_bases+=GetDepth(stab_end);
            stab_end = GetColumn(columnStore,stab_end->next);
            stab_width++;
          }
          if ( stab_bases == 0 ) break;
          //  Floating point 'instability' here?
          while( (double)stab_mm/(double)stab_bases > 0.02 ||  //  CNS_SEQUENCING_ERROR_EST
                 (double)stab_gaps/(double)stab_bases > .25  ){
            int mm=ColumnMismatch(stab);
            int gp=GetColumnBaseCount(stab,'-');
            int bps=GetDepth(stab);
            // move stab column ahead
            if ( stab_end->next != -1 ) {
              stab_mm+=ColumnMismatch(stab_end);
              stab_bases+=GetDepth(stab_end);
              stab_gaps+=GetColumnBaseCount(stab_end,'-');
              stab_end = GetColumn(columnStore,stab_end->next);
              stab_mm-=mm;
              stab_gaps-=gp;
              stab_bases-=bps;
              cum_mm+=mm;
              stab = GetColumn(columnStore,stab->next);
              win_length++;
            } else {
              break;
            }
          }
          *stab_bgn = stab->lid;
        }
        if ( win_length > 1 ) rc = win_length;
      }
      break;
    default:
      break;
  }
  return rc;
}

void ShowCalls(Abacus *abacus)
{
  int j;
  //  fprintf(stderr, "Calls=\n");
  for (j=0;j<abacus->columns;j++) 
    fprintf(stderr, "%c", abacus->calls[j]);
  fprintf(stderr, "\n");
}

static void 
GetReadsForAbacus(Read *reads, Abacus *abacus)
{
  int i, j, shift=0;
  char base;

#if 0
  fprintf(stderr, "rows=%lu shift=%c window_width=%lu \nReads= \n",
          abacus->rows, (char)abacus->shift, abacus->window_width);
#endif
  if (abacus->shift == UNSHIFTED) 
    shift = abacus->columns;
  else if (abacus->shift == RIGHT_SHIFT)
    shift = 2*abacus->columns; 
  for (i=0; i<abacus->rows; i++) {
    reads[i].id = i;
    reads[i].allele_id = -1;
    reads[i].ave_qv = 20.;        // qvs are hardly available for abacus
    reads[i].uglen = 0;
    for (j=0; j<abacus->columns; j++) {
      base = *GetAbacus(abacus,i,j);
      if (is_good_base(base))
        reads[i].bases[j] = base;                        
      if (base != '-')
        reads[i].uglen++;
    }
  }
}

static int 
base2int(char b)
{
  if (b == '-')             return 0;
  if (b == 'a' || b == 'A') return 1;
  if (b == 'c' || b == 'C') return 2;
  if (b == 'g' || b == 'G') return 3;
  if (b == 't' || b == 'T') return 4;
  if (b == 'n' || b == 'N') return 5;
  fprintf(stderr, "base2int b out of range");
  assert(0);
}

static void
GetConsensusForAbacus(VarRegion  *vreg, Read *reads, Abacus *abacus, 
                      char ***consensus)
{
  int i, j;
  char bases[CNS_NALPHABET] = {'-', 'A', 'C', 'G', 'T', 'N'};
  // Allocate memory for consensus
  *consensus = (char **)safe_malloc(2 * sizeof(char *));
  for (i=0; i<2; i++) {
    (*consensus)[i] = (char *)safe_malloc(3*abacus->window_width * sizeof(char));
    for (j=0; j<3*abacus->window_width; j++)
      (*consensus)[i][j] = '-';
  }

  // Call consensus
  for (i=0; i<3*abacus->window_width; i++) 
    {
      int bcount0[6] = {0, 0, 0, 0, 0, 0};
      int bcount1[6] = {0, 0, 0, 0, 0, 0};
      int best_count0=0, second_best_count0=0;
      int best_count1=0, second_best_count1=0;
      char cbase0, cbase1;
      for (j=0; j<abacus->rows; j++) {
#if 0
        fprintf(stderr, " reads[%d][%d]= %c\n", j, i, reads[j].bases[i]);
#endif
        if (is_good_base(reads[j].bases[i]))
          {
            if   (vreg->alleles[j].id == 0) 
              bcount0[base2int(reads[j].bases[i])]++;         
            else                       
              bcount1[base2int(reads[j].bases[i])]++;
          }
      }
      for (j=0; j<CNS_NALPHABET; j++) {
        if (best_count0 < bcount0[j]) {
          second_best_count0 = best_count0;
          best_count0 = bcount0[j];
          cbase0 = bases[j];
        }
        else if (  best_count0 >= bcount0[j] && 
                   second_best_count0 <  bcount0[j]) {
          second_best_count0  = bcount0[j];
        }
      }
      for (j=0; j<CNS_NALPHABET; j++) {
        if (best_count1 < bcount1[j]) {
          second_best_count1 = best_count1;
          best_count1 = bcount1[j];
          cbase1 = bases[j];
        }
        else if (  best_count1 >= bcount1[j] &&
                   second_best_count1 <  bcount1[j]) {
          second_best_count1  = bcount1[j];
        }
      }
      if (best_count0 == second_best_count0)
        (*consensus)[0][i] = 'N';
      else
        (*consensus)[0][i] = cbase0;
      if (best_count1 == second_best_count1)
        (*consensus)[1][i] = 'N';
      else
        (*consensus)[1][i] = cbase1;
    }
}

/* Create ungapped consensus sequences and map them to gapped consensus sequences */
static void
MapConsensus(int ***imap, char **consensus,  char ***ugconsensus, 
             int len, int *uglen)
{
  int i, j, k;
  uglen[0] = uglen[1] = 0;
  *ugconsensus = (char **)safe_malloc(2*sizeof(char *));
  *imap        = (int  **)safe_malloc(2*sizeof(int  *));
  for (i=0; i<2; i++)
    {
      (*ugconsensus)[i] = (char *)safe_malloc(len*sizeof(char));    
      (*imap)[i]        = (int  *)safe_malloc(len*sizeof(int ));
      for (j=0; j<len; j++)
        (*imap)[i][j] = j;                
      k=0;
      for (j=0; j<len; j++)
        {
          if (consensus[i][j] != '-')
            {
              (*ugconsensus)[i][k] = consensus[i][j];
              (*imap)[i][k] = j;
              k++;
            }   
        }
      uglen[i] = k;
    }
}

/* Count gaps in the short and long consensus sequences */
static void
CountGaps(char **consensus, int len, int *gapcount)
{
  int i, j, first_base, last_base;
   
  for (i=0; i<2; i++)
    {
      last_base = len-1;
      while ((last_base > 0) && (consensus[i][last_base] == '-'))
        last_base--;

      gapcount[i] = 0;
      first_base = -1; 
      for (j=0; j<last_base + 1; j++)
        {
          if (consensus[i][j] != '-')
            first_base = j;

          if ((first_base >= 0) && (consensus[i][j] == '-'))
            gapcount[i]++;
        }
    }
}

/*
  Find an adjusted left boundary for long and short ungapped sequences (that is, 
  the leftmost position starting from which the sequences will match):
  6.1 select a size of a "probe" k-mer, say, 3
  6.2 scan the short ungapped consensus from the left to the right
  6.3 for each position in the short ungapped consensus, get k-mer starting at 
  this position
  6.4 find the leftmost position of this k-mer in the long ungapped sequence
  6.5 set adjusted left boundary of short ungapped consensus to the position of
  k-mer with leftmost occurrence in the long ungapped sequence
*/
static void
FindAdjustedLeftBounds(int *adjleft, char **ugconsensus, int *uglen,
                       int short_allele, int long_allele)
{
  int   s, l;
  char *ps, *pl;

  adjleft[short_allele] = uglen[short_allele]-1;
  adjleft[long_allele]  = uglen[long_allele]-1;
  for (s=0; s<uglen[short_allele] - MSTRING_SIZE; s++)
    {
      ps = ugconsensus[short_allele] + s;
      for (l=0; l<uglen[long_allele] - MSTRING_SIZE; l++)
        {
          pl = ugconsensus[long_allele] + l;    
          if (strncmp(pl, ps, MSTRING_SIZE) == 0)
            {
              if (adjleft[0] + adjleft[1] > s+l) 
                { 
                  adjleft[long_allele]  = l;
                  adjleft[short_allele] = s;
                }
            }          
        }
    }          
  if ((adjleft[long_allele]  == uglen[long_allele]-1) &&
      (adjleft[short_allele] == uglen[short_allele]-1))  
    {
      adjleft[short_allele] = 0;
      adjleft[long_allele]  = 0;
    } 
}

static void
FindAdjustedRightBounds(int *adjright,  char **ugconsensus, int *uglen,
                        int short_allele, int long_allele)
{
  int   s, l;
  char *ps, *pl;

  adjright[short_allele] = uglen[short_allele]-1;
  adjright[long_allele]  = uglen[ long_allele]-1;
  for (s=uglen[short_allele] - MSTRING_SIZE-1; s>= 0; s--)
    {
      ps = ugconsensus[short_allele] + s;
      for (l=uglen[long_allele] - MSTRING_SIZE-1; l>=0; l--)
        {
          pl = ugconsensus[long_allele] + l;
          if (strncmp(pl, ps, MSTRING_SIZE) == 0)
            {
              if (adjright[0] + adjright[1] > 
                  uglen[short_allele]-1 -(s+MSTRING_SIZE)
                  +uglen[long_allele] -1 -(l+MSTRING_SIZE))
                {
                  adjright[long_allele] = 
                    uglen[long_allele]-1-(l+MSTRING_SIZE);
                  adjright[short_allele] = 
                    uglen[short_allele]-1 -(s+MSTRING_SIZE);
                }
            }
        }
    }
  if ((adjright[long_allele]  == uglen[long_allele]-1) &&
      (adjright[short_allele] == uglen[short_allele]-1))
    {
      adjright[short_allele] = 0;
      adjright[long_allele]  = 0;
    }
}

static void
GetLeftScore(char **ugconsensus, int *uglen, int **imap, int *adjleft, 
             int short_allele, int long_allele, int *maxscore, int *maxpos)
{
  int i, score = 0;
    
  *maxscore = 0;
  *maxpos   = adjleft[short_allele];
  i = 0;
  while ((i < uglen[short_allele] - adjleft[short_allele]) &&
         (i < uglen[ long_allele] - adjleft[ long_allele]))
    {
      int lpos = i + adjleft[long_allele];
      int spos = i + adjleft[short_allele]; 
      if (ugconsensus[short_allele][spos] == ugconsensus[long_allele][lpos])
        score++;
      else
        score--;
      if (*maxscore < score)
        {
          *maxscore = score;
          *maxpos   = spos;
        }
      i++;
    }
  /* Position in ungapped consensus */
  *maxpos = imap[short_allele][*maxpos];
}

static void
GetRightScore(char **ugconsensus, int *uglen, int **imap, int *adjright, 
              int short_allele, int long_allele, int *maxscore, int *maxpos)
{
  int i, j, score = 0;

  *maxscore = 0;
  *maxpos   = uglen[short_allele]-1-adjright[short_allele];
  i = uglen[long_allele]-1;
  j = uglen[short_allele]-1;
  while ((j >= adjright[short_allele]) &&
         (i >= adjright[ long_allele]))
    {
      int lpos = i - adjright[long_allele];
      int spos = j - adjright[short_allele];
      if (ugconsensus[short_allele][spos] == ugconsensus[long_allele][lpos])
        score++;
      else
        score--;
      if (*maxscore < score)
        {
          *maxscore = score;
          *maxpos = spos;
        }
      i--;
      j--;
    }
  /* Position in ungapped consensus */
#if 0
  fprintf(stderr, "long_allele=%d  maxpos =%d \n", long_allele, *maxpos);
#endif
  *maxpos = imap[short_allele][*maxpos];
}

static void
AdjustShiftingInterfaces(int *lpos, int *rpos, int lscore, int rscore,
                         int *adjleft, int *adjright, int long_allele, int short_allele)
{
#if 0
  fprintf(stderr, "\nlpos=%d rpos=%d lscore=%d rscore=%d \n", *lpos, *rpos, lscore, rscore);
  fprintf(stderr, "adjleft = %d %d  adjright= %d %d \n", adjleft[0], adjleft[1],
          adjright[0], adjright[1]);
#endif
  if (adjleft[long_allele] > 5)
    {
      *lpos = -1;
      lscore = -1;
    }
  if (adjright[long_allele] > MAX_SIZE_OF_ADJUSTED_REGION)
    {
      *rpos = -1;
      rscore = -1;
    }

  /* Set teh posiytion of shifting interface */
  if (*lpos <= *rpos)
    {
    }
  else /* lpos >  rpos */
    {
      if ((lscore > 0) && (rscore > 0))
        {
          //         if (adjleft[short_allele] <= adjright[short_allele])
          if (lscore > rscore)
            *rpos = *lpos;
          else
            *lpos = *rpos;
        }
      else if ((lscore > 0) && (rscore <= 0))
        *rpos = -1;
      else
        *lpos = -1;
    }
}

void  GetTemplateForAbacus(char **template, char **consensus, int len,
                           char **ugconsensus, int *uglen, int lpos, int rpos, int **imap,
                           int *adjleft, int *adjright, int short_allele, int long_allele)
{
  int i, j;

  *template = (char *)safe_malloc(len*sizeof(char));
  for (i=0; i<len; i++)
    (*template)[i] = consensus[long_allele][i];

  /* Set Ns in the left part of the template */
  i = 0;
  while ((imap[long_allele][i] <= lpos) &&
         (i < uglen[short_allele] - adjleft[short_allele]) &&
         (i < uglen[ long_allele] - adjleft[ long_allele]))
    {
      int lpos = i + adjleft[long_allele];
      int spos = i + adjleft[short_allele];
      if ((ugconsensus[short_allele][spos] != ugconsensus[long_allele][lpos]) &&
          ((*template)[imap[long_allele][lpos]] != '-'))
        (*template)[imap[long_allele][lpos]] = 'n';
      i++;
    }

  /* template bases before adjusted left boundary should be 'N's */
  if (adjleft[long_allele] > 0 && lpos > 0)
    {
      i = imap[long_allele][adjleft[long_allele]-1];
      j = 0;
      while ((j < adjleft[short_allele]) && (i >= 0))
        {
          if (consensus[long_allele][i] != '-')
            {
              (*template)[i] = 'n';
              j++;
              i--;
            }
          else
            i--;
        }
    }

  /* Set Ns in the right part of the template */
  i = uglen[long_allele]-1-adjright[long_allele];
  j = uglen[short_allele]-1-adjright[short_allele];
  while ((i >= adjleft[long_allele]) &&
         (j >= adjleft[short_allele]) &&
         (imap[long_allele][i] > rpos))
    {
      if ((ugconsensus[short_allele][j] != ugconsensus[long_allele][i]) &&
          ((*template)[imap[long_allele][i]] != '-'))
        {
          (*template)[imap[long_allele][i]] =  'n';
        }
      i--;
      j--;
    }

  /* template bases after adjusted right boundary should be 'N's */
  if ((adjright[long_allele] > 0) && (rpos > 0))
    {
      for (i = uglen[long_allele]-adjright[long_allele];
           i < uglen[long_allele];
           i++)
        {
          j = imap[long_allele][i];
          if (consensus[long_allele][j] != '-')
            {
              (*template)[i] = 'n';
            }
        }
    }
}

int RefineWindow(MANode *ma, Column *start_column, int stab_bgn,
                 CNS_Options *opp ) 
{
  int orig_columns=0, left_columns=0, right_columns=0, best_columns=0;
  // Mismatch, gap and total scores:
  int32 orig_mm_score=0, left_mm_score=0, right_mm_score=0, best_mm_score=0;
  int32 orig_gap_score=0, left_gap_score=0, right_gap_score=0,
    best_gap_score = 0;
  int32 orig_total_score, left_total_score, right_total_score, 
    best_total_score;
  int32 max_element = 0, score_reduction;
  BaseCount abacus_count;
  Abacus *left_abacus, *orig_abacus, *right_abacus, *best_abacus;
  VarRegion  vreg;

  orig_abacus = CreateAbacus(ma->lid,start_column->lid,stab_bgn);
  SetDefault(&vreg);
  vreg.nr = orig_abacus->rows;
  //ShowAbacus(orig_abacus);
  orig_mm_score = ScoreAbacus(orig_abacus,&orig_columns);

  // Process reads of the original abacus
  AllocateDistMatrix(&vreg, 0);
  AllocateMemoryForReads(&vreg.reads, orig_abacus->rows, orig_abacus->columns,
                         QV_FOR_MULTI_GAP);
  vreg.beg = 0;
  vreg.end = orig_abacus->columns-1;
  GetReadsForAbacus(vreg.reads, orig_abacus);
#if 0
  OutputReads(stderr, vreg.reads, vreg.nr, orig_abacus->columns);
#endif
  PopulateDistMatrix(vreg.reads, orig_abacus->columns, &vreg);
#if DEBUG_ABACUS
  OutputDistMatrix(stderr, &vreg);
#endif
  AllocateMemoryForAlleles(&vreg.alleles, vreg.nr, &vreg.na);
  ClusterReads(vreg.reads, vreg.nr, vreg.alleles, &vreg.na,
               &vreg.nca, vreg.dist_matrix);
  SortAllelesByLength(vreg.alleles, vreg.nca, vreg.reads);

  RefineOrigAbacus(orig_abacus, vreg);
  //ShowAbacus(orig_abacus);
  orig_mm_score = ScoreAbacus(orig_abacus, &orig_columns);

#if 0
  fprintf(stderr, "vreg.alleles= ");
  for (i=0; i<vreg.nr; i++)
    fprintf(stderr, "%d", vreg.alleles[i]);

#endif

#if DEBUG_ABACUS
  fprintf(stderr, "\n\nOrigCalls=\n");
  ShowCalls(orig_abacus);
  fprintf(stderr, "Abacus=\n");
  ShowAbacus(orig_abacus);
  fprintf(stderr, "\n");
#endif
  //ShowAbacus(orig_abacus);
  left_abacus = CloneAbacus(orig_abacus);
  left_mm_score = LeftShift(left_abacus, vreg, &left_columns);
#if DEBUG_ABACUS
  fprintf(stderr, "\n\nLeftShiftCalls=\n");
  ShowCalls(left_abacus);
  fprintf(stderr, "Abacus=\n");
  ShowAbacus(left_abacus);
  fprintf(stderr, "\n");
#endif
  right_abacus = CloneAbacus(orig_abacus);
  right_mm_score = RightShift(right_abacus, vreg, &right_columns);
#if DEBUG_ABACUS
  fprintf(stderr, "\n\nRightShiftCalls=\n");
  ShowCalls(right_abacus);
  fprintf(stderr, "Abacus=\n");
  ShowAbacus(right_abacus);
  fprintf(stderr, "\n");
#endif
  //fprintf(stderr,"Abacus Report:\norig_mm_score: %d left_mm_score: %d right_mm_score: %d\n",
  //             orig_mm_score,left_mm_score,right_mm_score);
  //ShowAbacus(left_abacus);
  //ShowAbacus(right_abacus);
  // determine best score and apply abacus to real columns
  orig_gap_score  = AffineScoreAbacus(orig_abacus);
  left_gap_score  = AffineScoreAbacus(left_abacus);
  right_gap_score = AffineScoreAbacus(right_abacus);
  best_abacus     = orig_abacus;
  best_columns    = orig_columns;
  best_gap_score  = orig_gap_score;
  best_mm_score   = orig_mm_score;
  orig_total_score  = orig_mm_score  + orig_columns  + orig_gap_score;
  left_total_score  = left_mm_score  + left_columns  + left_gap_score; 
  right_total_score = right_mm_score + right_columns + right_gap_score;
  best_total_score  = orig_total_score;

#if DEBUG_ABACUS
  fprintf(stderr, "In RefineWindow: beg= %lu end= %d\n", 
          start_column->lid, stab_bgn);
  fprintf(stderr, "    abacus->columns= %d, abacus->rows= %d\n", 
          orig_abacus->columns, orig_abacus->rows);
  fprintf(stderr, "    w_width left= %d orig= %d right= %d\n", 
          left_abacus->window_width, orig_abacus->window_width,
          right_abacus->window_width);
  fprintf(stderr, "    mm_score left= %d orig= %d right= %d\n", 
          left_mm_score, orig_mm_score, right_mm_score);
  fprintf(stderr, "     columns left= %d orig= %d right= %d\n", 
          left_columns, orig_columns, right_columns);
  fprintf(stderr, "   gap_score left= %d orig= %d right= %d\n", 
          left_gap_score, orig_gap_score, right_gap_score);
  fprintf(stderr, " total_score left= %d orig= %d right= %d\n", 
          left_total_score, orig_total_score, right_total_score);
#endif

  // Use the total score to refine the abacus
  if ( left_total_score < orig_total_score || 
       right_total_score < orig_total_score )
    {
      if ( left_total_score <= right_total_score ) {
        score_reduction += orig_total_score - left_total_score;
        //fprintf(stderr,"\nTry to apply LEFT abacus:\n");
        //ShowAbacus(left_abacus);
        GetAbacusBaseCount(left_abacus, &abacus_count);
#if 0
        fprintf(stderr, " Applying left abacus\n");
#endif
        best_abacus      = left_abacus;
        best_mm_score    = left_mm_score;
        best_columns     = left_columns;
        best_gap_score   = left_gap_score;
        best_total_score = left_total_score;
      }
      else
        {
          score_reduction += orig_total_score - right_total_score;
          //fprintf(stderr,"\nTry to apply RIGHT abacus:\n");
          //ShowAbacus(right_abacus);
          GetAbacusBaseCount(right_abacus,&abacus_count);
#if 0
          fprintf(stderr, " Applying right abacus\n");
#endif
          best_abacus      = right_abacus;
          best_mm_score    = right_mm_score;
          best_columns     = right_columns;
          best_gap_score   = right_gap_score;
          best_total_score = right_total_score;
        }
    }

#if 0
  fprintf(stderr, "Best Abacus Before MixedShift=\n");
  ShowAbacus(best_abacus);
#endif
  { 
    int i, j;
    char  **consensus=NULL, **ugconsensus=NULL, *template=NULL;
    int   **imap=NULL, uglen[2]={0,0}, adjleft[2]={-1,-1}, adjright[2]={-1,-1};
    int     gapcount[2], short_allele=-1, long_allele=-1;
    int     lscore=0, rscore=0, lpos=-1, rpos=-1;
    int     mixed_columns=0;
    int32   mixed_mm_score=0, mixed_gap_score=0;
    Abacus *mixed_abacus=NULL;

    GetConsensusForAbacus(&vreg, vreg.reads, best_abacus, &consensus);
#if 0
    fprintf(stderr, "\nconsensus[0]=\n");
    for (i=0; i<3*best_abacus->window_width; i++)
      fprintf(stderr, "%c", consensus[0][i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "consensus[1]=\n");
    for (i=0; i<3*best_abacus->window_width; i++)
      fprintf(stderr, "%c", consensus[1][i]);
    fprintf(stderr, "\n\n");
#endif
    CountGaps(consensus, 3*best_abacus->window_width, gapcount);
    short_allele = (gapcount[0] >= gapcount[1]) ? 0 : 1;
    long_allele  = (gapcount[0] <  gapcount[1]) ? 0 : 1;
#if 0
    fprintf(stderr, "gapcounts[0, 1] = %d %d\n", gapcount[0], gapcount[1]);
#endif
    if (gapcount[short_allele] == 0)
      {
#if 0
        fprintf(stderr, "No MixedShift will be performed: gapcount[short_allele] = %d\n", gapcount[short_allele]);
#endif
        ApplyAbacus(best_abacus, opp);

        DeleteAbacus(orig_abacus);
        DeleteAbacus(left_abacus);
        DeleteAbacus(right_abacus);
        if (vreg.nr > 0)
          {
            int j;
            for (j=0; j<vreg.nr; j++)
              {
                safe_free(vreg.alleles[j].read_ids);
                safe_free(vreg.alleles[j].read_iids);
                safe_free(vreg.dist_matrix[j]);
                safe_free(vreg.reads[j].bases);
                safe_free(vreg.reads[j].qvs);
              }
            safe_free(vreg.reads);
            safe_free(vreg.alleles);
            safe_free(vreg.dist_matrix);
          }
        safe_free(consensus[0]);
        safe_free(consensus[1]);
        safe_free(consensus);
        return score_reduction;
      }

    /* Now try the mixed consensus */
    MapConsensus(&imap, consensus, &ugconsensus, 3*best_abacus->window_width,
                 uglen);
    if ((uglen[0] < MSTRING_SIZE) || (uglen[1] < MSTRING_SIZE))
      {
#if 0
        fprintf(stderr, "No MixedShift will be performed: uglen = %d %d\n", uglen[0], uglen[1]);
#endif
        ApplyAbacus(best_abacus, opp);

        DeleteAbacus(orig_abacus);
        DeleteAbacus(left_abacus);
        DeleteAbacus(right_abacus);
        if (vreg.nr > 0)
          {
            int j;
            for (j=0; j<vreg.nr; j++)
              {
                safe_free(vreg.alleles[j].read_ids);
                safe_free(vreg.alleles[j].read_iids);
                safe_free(vreg.dist_matrix[j]);
                safe_free(vreg.reads[j].bases);
                safe_free(vreg.reads[j].qvs);
              }
            safe_free(vreg.reads);
            safe_free(vreg.alleles);
            safe_free(vreg.dist_matrix);
          }
        safe_free(consensus[0]);
        safe_free(consensus[1]);
        safe_free(consensus);
        for (i=0; i<2; i++)
          {
            safe_free(ugconsensus[i]);
            safe_free(imap[i]);
          }
        safe_free(ugconsensus);
        safe_free(imap);
        return score_reduction;
      }

#if 0
    fprintf(stderr, "\nuglen[0]=%d ugconsensus[0] =\n", uglen[0]);
    for (i=0; i<uglen[0]; i++)
      fprintf(stderr, "%c", ugconsensus[0][i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "uglen[1]=%d ugconsensus[1] =\n", uglen[1]);
    for (i=0; i<uglen[1]; i++)
      fprintf(stderr, "%c", ugconsensus[1][i]);
    fprintf(stderr, "\n\n");
#endif
    FindAdjustedLeftBounds(adjleft, ugconsensus, uglen, short_allele,
                           long_allele);
    FindAdjustedRightBounds(adjright, ugconsensus, uglen, short_allele,
                            long_allele);
#if 0
    fprintf(stderr, "Adjusted left bounds 0, 1= %d %d \n", adjleft[0], adjleft[1]);
    fprintf(stderr, "Adjusted right bounds 0, 1= %d %d \n", adjright[0], adjright[1]);
#endif
    GetLeftScore(ugconsensus, uglen, imap, adjleft,
                 short_allele, long_allele, &lscore, &lpos);
    GetRightScore(ugconsensus, uglen, imap, adjright,
                  short_allele, long_allele, &rscore, &rpos);
    AdjustShiftingInterfaces(&lpos, &rpos, lscore, rscore,
                             adjleft, adjright, long_allele, short_allele);
    GetTemplateForAbacus(&template, consensus, 3*best_abacus->window_width,
                         ugconsensus, uglen, lpos, rpos, imap, adjleft, adjright,
                         short_allele, long_allele);

    mixed_abacus = CloneAbacus(orig_abacus);
#if 0
    {
      fprintf(stderr, "Template = \n");
      for (i=0; i<3*best_abacus->window_width; i++)
        fprintf(stderr, "%c", template[i]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "Start calling MixedShift\n\n");
    fprintf(stderr, "   Final lpos=%d rpos=%d window_width=%d long_allele=%d short_allele=%d\n",
            lpos, rpos, best_abacus->window_width, long_allele, short_allele);
#endif
    mixed_mm_score = MixedShift(mixed_abacus, &mixed_columns, vreg, lpos, rpos,
                                template, long_allele, short_allele);

#if 0
    fprintf(stderr, "Mixed abacus=\n");
    ShowAbacus(mixed_abacus);
    fprintf(stderr, "End calling MixedShift\n\n");
    fprintf(stderr, "mixed_mm_score=%d bast_score=%d\n", mixed_mm_score, best_mm_score);
    fprintf(stderr, "mixed_columns=%d best_columns=%d\n", mixed_columns, best_columns);
    fprintf(stderr, "mixed_gap_score=%d best_gap_score=%d\n", mixed_gap_score, best_gap_score);
#endif
    mixed_gap_score  = AffineScoreAbacus(mixed_abacus);
#if 0
    fprintf(stderr, "mixed_gap_score=%d best_gap_score=%d mixed_columns=%d best_columns=%d mixed_mm_score=%d best_mm_score=%d\n", mixed_gap_score, best_gap_score, mixed_columns, best_columns, mixed_mm_score, best_mm_score);
#endif
    if ( (mixed_gap_score <  best_gap_score) ||
         ((mixed_gap_score == best_gap_score) && (mixed_columns < best_columns))
         ||
         ((mixed_gap_score == best_gap_score) && (mixed_columns == best_columns) &&
          (mixed_mm_score < best_mm_score)))
      {
        best_abacus    = mixed_abacus;
        best_mm_score  = mixed_mm_score;
        best_columns   = mixed_columns;
        best_gap_score = mixed_gap_score;
      }
#if 0
    ShowCalls(best_abacus);
    fprintf(stderr, "Best Abacus after MixedShift =\n");
    ShowAbacus(best_abacus);
#endif

    //      OutputDistMatrix(stderr, &vreg);

#if 0
    {
      int j;
      fprintf(stderr, "Consensus0 =\n");
      for (j=0; j<3*best_abacus->window_width; j++)
        fprintf(stderr, "%c", consensus[0][j]);
      fprintf(stderr, "\n\n");
      fprintf(stderr, "Consensus1 =\n");
      for (j=0; j<3*best_abacus->window_width; j++)
        fprintf(stderr, "%c", consensus[1][j]);
      fprintf(stderr, "\n\n");
    }
#endif

    /* Otherwise, try to do a more sophisticated shift:
     * - only shifting reads of the shortest allele
     * - only within a subregion of abacus window where the alleles match
     */
#if 0
    fprintf(stderr, "Applying the Best abacus\n");
#endif
    ApplyAbacus(best_abacus, opp);

    //      fprintf(stderr, "vreg.nr = %d\n", vreg.nr);

    DeleteAbacus(orig_abacus);
    DeleteAbacus(left_abacus);
    DeleteAbacus(right_abacus);
    DeleteAbacus(mixed_abacus);
    {
      safe_free(consensus[0]);
      safe_free(consensus[1]);
      safe_free(consensus);
      safe_free(ugconsensus[0]);
      safe_free(ugconsensus[1]);
      safe_free(ugconsensus);
      safe_free(imap[0]);
      safe_free(imap[1]);
      safe_free(imap);
      safe_free(template);
    }
    if (vreg.nr > 0)
      {
        int j;
        for (j=0; j<vreg.nr; j++)
          {
            safe_free(vreg.alleles[j].read_ids);
            safe_free(vreg.alleles[j].read_iids);
            safe_free(vreg.dist_matrix[j]);
            safe_free(vreg.reads[j].bases);
            safe_free(vreg.reads[j].qvs);
          }
        safe_free(vreg.reads);
        safe_free(vreg.alleles);
        safe_free(vreg.dist_matrix);
      }
  }
  return score_reduction;
}


//*********************************************************************************
// Abacus Refinement:
//   AbacusRefine contains the logic for sweeping through the multialignment,
//   and identifying candidate windows for refinement.
//   Each window is cast into an abacus, which is left and right shifted.
//   The best resulting base arrangement (if different from input) is then
//   applied to window of the MultiAlignment
//*********************************************************************************

int AbacusRefine(MANode *ma, int32 from, int32 to, CNS_RefineLevel level,
                 CNS_Options *opp) 
{
  // from and to are in ma's column coordinates
  int32 sid, eid, stab_bgn;
  int32 ma_length = GetMANodeLength(ma->lid);
  int32 score_reduction=0;
  int32 orig_length = ma_length; 
  int32 refined_length = orig_length;
  Column *start_column;
  int i;
 
  if(from < 0 || from > ma_length-1){
    fprintf(stderr, "AbacusRefine range (from) invalid");
    assert(0);
  }
  if ( to == -1 ) to = ma_length-1;
  if(to <= from || to > ma_length-1){
    fprintf(stderr, "AbacusRefine range (to) invalid");
    assert(0);
  }

  ResetIndex(abacus_indices,GetNumFragments(fragmentStore));
  sid = *Getint32(ma->columns,from);   // id of the starting column
  eid = *Getint32(ma->columns,to);     // id of the ending column
  start_column = GetColumn(columnStore,sid);

  while (start_column->lid != eid) 
    {
      int window_width = IdentifyWindow(&start_column,&stab_bgn, level);
      // start_column stands as the candidate for first column in window 
      // look for window start and stop

      if (window_width > 0)
        {
#if DEBUG_ABACUS
          fprintf(stderr, "In AbacusRefine window_width= %d\n", window_width);
#endif
          //
          // refine in window
          if ( start_column->prev == -1 ) {
            // if start_column->prev == -1, insert a gap column for maneuvering room
            int32 newbead;
            Bead *firstbead;
            firstbead = GetBead(beadStore,GetBead(beadStore,start_column->call)->down);
            newbead   = AppendGapBead(firstbead->boffset);
            //  beadStore could have changed position, so grab
            //  firstbead again; we could instead have cached boffset
            //  and column_index
            firstbead = GetBead(beadStore,GetBead(beadStore,start_column->call)->down);
            fprintf(stderr,"Adding gapbead %d after %d to add abacus room for abacus abutting left of multialignment\n",
                    newbead, firstbead->boffset);
            ColumnAppend(firstbead->column_index,newbead);
          }

          //  if the window is too big, there's likely a polymorphism
          //  that won't respond well to abacus, so skip it
          //
          //  BPW saw crashes with large window_width's (1333, 3252,
          //  1858, 675, 855, 1563, 601, 1102).  The longest
          //  window_width that worked was 573.  Previous versions
          //  used 100 here.  Not sure what it should be.
          //
          if ( window_width < MAX_WINDOW_FOR_ABACUS_REFINE )
            score_reduction += RefineWindow(ma,start_column,stab_bgn, opp); 

          start_column = GetColumn(columnStore, stab_bgn);
        }
      start_column = GetColumn(columnStore, stab_bgn);
    }
  {
    int32 nv=0;
    IntMultiVar *vl=NULL;
    RefreshMANode(ma->lid, 1, opp, &nv, &vl, 0, 0);
    safe_free(vl);
  }
  refined_length = GetMANodeLength(ma->lid);
  if ( refined_length < orig_length ) 
    {
      //fprintf(stderr,"Column reduction = ", orig_length-GetMANodeLength(ma->lid));
    }
  return score_reduction;
}


int MANode2Array(MANode *ma, int *depth, char ***array, int ***id_array,
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
    if (max_depth > MAX_ALLOWED_MA_DEPTH )  return 0;
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
      fbgn = GetColumn(columnStore,(GetBead(beadStore,frag->beads))->column_index)->ma_index;
      fend = GetColumn(columnStore,
                       (GetBead(beadStore,frag->beads+frag->length-1))->column_index)->ma_index+1;
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
        fb = GetBead(beadStore,frag->beads);
        bcolumn =  GetColumn(columnStore,fb->column_index);
        if(!CreateFragmentBeadIterator(fid,&fi)){
          fprintf(stderr, "MANode2Array CreateFragmentBeadIterator failed");
          assert(0);
        }
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

int RealignToConsensus(int32 mid,
                       char *sequence,
                       char *quality,
                       int32 fid_bgn, 
                       int32 fid_end, 
                       Overlap *(*COMPARE_FUNC)(COMPARE_ARGS),
                       CNS_Options *opp) 
{
  // this is prototype code code in development
  static char cnstmpseq[2*AS_READ_MAX_LEN+1];
  static char cnstmpqlt[2*AS_READ_MAX_LEN+1];
  int i;
  MANode *ma_realigned=NULL;
  Fragment *afrag;
  Bead *afirst;
  Column *col;
  int cns_fid;
  char *stmp=cnstmpseq,*qtmp=cnstmpqlt;
  int32 ahang;
  int32 aoffset;
  int olap_success;
  int ovl=0;
  OverlapType otype;
  static VA_TYPE(int32) *trace=NULL;
  if ( trace == NULL ) {
    trace = CreateVA_int32(AS_READ_MAX_LEN);
  } else {
    ResetVA_int32(trace);
  }
  for (i=0;;i++){
    if ( sequence[i] != '\0' ) {
      if ( sequence[i] != '-' ) {
        *stmp++ = sequence[i];
        *qtmp++ = quality[i];
      }
    } else {
      *stmp=sequence[i];
      *qtmp=quality[i];
      break;
    }
  }
  if ( fid_end == -1 ) fid_end = GetNumFragments(fragmentStore);
  cns_fid = AppendArtificialFragToLocalStore( (FragType)'c',
                                              0, 
                                              0,
                                              0,
                                              (UnitigType) 'X',
                                              cnstmpseq,cnstmpqlt,stmp-cnstmpseq);
  ma_realigned = CreateMANode(GetNumMANodes(manodeStore));
  assert(ma_realigned!=NULL);
  SeedMAWithFragment(ma_realigned->lid, cns_fid, 0, opp);
  for (i=fid_bgn;i<fid_end;i++) {

    afrag = GetFragment(fragmentStore,i);
    afirst = GetBead(beadStore,afrag->beads);
    col = GetColumn(columnStore,afirst->column_index);
    aoffset = col->call;
    // starting consensus bead  is call of columne where first a basepair is aligned
    olap_success = GetAlignmentTrace(cns_fid, aoffset, i, &ahang, ovl, trace, &otype,DP_Compare,SHOW_OLAP,0);
    UnAlignFragment(i);
    ApplyAlignment(cns_fid,aoffset,i,ahang,Getint32(trace,0));
    afrag->deleted = 0;
    GetMultiAlignInStore(unitigStore, mid);
    {
      IntMultiVar *vl=NULL;
      int32 nv=0;
      RefreshMANode(mid, 0, opp, &nv, &vl, 0, 0);
    }
  }
 
  return 1;
}

int32
GetFragmentIndex(IntFragment_ID ident2, IntMultiPos *positions, int num_frags)
{
  int i;
  for (i=0; i<num_frags; i++)
    {
      if (ident2 == positions[i].ident)
        {
          return i;
        }
    }
  return -1;
}

int MultiAlignUnitig(IntUnitigMesg *unitig, 
                     GateKeeperStore *fragStore,
                     VA_TYPE(char) *sequence,
                     VA_TYPE(char) *quality, 
                     VA_TYPE(int32) *deltas, 
                     CNS_PrintKey printwhat, 
                     int mark_contains, 
                     Overlap *(*COMPARE_FUNC)(COMPARE_ARGS),
                     CNS_Options *opp) 
{
  // The function will return 0 if successful, and -1 if unsuccessful 
  // (due to overlap failure)
  int32 fid,i,align_to;
  int32 num_reads=0,num_columns=0;
#ifdef ALIGN_TO_CONSENSUS
  int32 aoffset;
#endif 

  int do_rez=1;
  // command line arg that is now obsolete
  // mark_contains is used in the case where post-unitigging processes 
  // (SplitUnitig, extendClearRange,e.g.) are used to re-align unitigs after 
  // fragments have been altered... With an extended clear range, 
  // the fragment may now "contain" another which it used to be a dovetail 
  // relationship with, or which used to contain it.  The mark_contains flag 
  // tells MultiAlignUnitig that there may be such a new relationship, and that 
  // it should be detected and marked as the alignment is being formed.
  // Without this marking, the multialignment is likely to have pieces 
  // which are not properly aligned, and which will appear as block indels 
  // (large gap-blocks) which will foil future overlaps involving
  // the consensus sequence of this "reformed" unitig

  int score_reduction;
  int complement;
  MANode *ma;
  SeqInterval *offsets;
  static VA_TYPE(int32) *trace=NULL;
  IntMultiPos *positions=unitig->f_list;
  int num_frags = unitig->num_frags;
  int unitig_forced = 0;
#ifdef NEW_UNITIGGER_INTERFACE
  int32 *is_pointed = NULL, *is_aligned = NULL, i_afrag;
#endif

  if ((num_frags == 1) && 
      (positions[0].position.bgn == positions[0].position.end))
    {
      fprintf(stderr, "MultiAlignUnitig()-- Warning: unitig %d contains a single fragment of length 0 !\n", unitig->iaccession);
      return EXIT_FAILURE;
    }

  //  Make sure that we have valid options here, we then reset the
  //  pointer to the freshly copied options, so that we can always
  //  assume opp is a valid pointer
  //

  ALIGNMENT_CONTEXT=AS_CONSENSUS;
  gkpStore=fragStore;

  RALPH_INIT = InitializeAlphTable();
  offsets = (SeqInterval *) safe_calloc(num_frags,sizeof(SeqInterval));
  for (i=0;i<num_frags;i++) {
    num_columns = ( positions[i].position.bgn>num_columns)? positions[i].position.bgn : num_columns;
    num_columns = ( positions[i].position.end>num_columns)? positions[i].position.end : num_columns;
  }
  ResetStores(num_frags,num_columns);

  fragmentMap = CreateScalarHashTable_AS(2*(num_frags));

  for (i=0;i<num_frags;i++) {
    complement = (positions[i].position.bgn<positions[i].position.end)?0:1;
    switch (positions[i].type) 
      {
        case AS_READ:
        case AS_EXTR:
        case AS_TRNR:
          {
            num_reads++;
            if (HASH_SUCCESS != InsertInHashTable_AS(fragmentMap,positions[i].ident, 0, 1, 0)) {
              fprintf(stderr,"MultiAlignUnitig()-- Failure to insert ident %d in hashtable\n", positions[i].ident); 
              assert(0);
            }

            fid = AppendFragToLocalStore(positions[i].type, 
                                         positions[i].ident, 
                                         complement,
                                         positions[i].contained,
                                         AS_OTHER_UNITIG, NULL);

            offsets[fid].bgn = complement ? positions[i].position.end : positions[i].position.bgn;
            offsets[fid].end = complement ? positions[i].position.bgn : positions[i].position.end;
            break;
          }
        case AS_UNITIG:
        default:
          {
            fprintf(stderr, "MultiAlignUnitig()-- Failed to determine the type of fragment %d in unitig %d\n",
                    i, unitig->iaccession);
            DeleteHashTable_AS(fragmentMap);
            DeleteMANode(ma->lid);
            return EXIT_FAILURE;
          }
      }
  }

  ma = CreateMANode(unitig->iaccession);
  assert(ma->lid == 0);
  if ( trace == NULL ) {
    trace = CreateVA_int32(AS_READ_MAX_LEN);
  } else {
    ResetVA_int32(trace);
  }

  SeedMAWithFragment(ma->lid, GetFragment(fragmentStore,0)->lid,0, opp);

  // Now, loop on remaining fragments, aligning to:
  //    a)  containing frag (if contained)
  // or b)  previously aligned frag

#ifdef NEW_UNITIGGER_INTERFACE
  is_pointed = safe_calloc(num_frags,sizeof(int32));
  is_aligned = safe_calloc(num_frags,sizeof(int32));
  is_aligned[0] = 1;
#endif

  for (i=1;i<num_frags;i++) {
    Fragment *bfrag = GetFragment(fragmentStore,i);
    Fragment *afrag = NULL;
    int ahang;
    int olap_success=0;
    int ovl=0;
    OverlapType otype;

#ifdef NEW_UNITIGGER_INTERFACE
    i = 0;
    afrag = GetFragment(fragmentStore,i);
    bfrag = NULL;
#endif

    // check whether contained, if so
    // align_to = containing
    // else 
    int frag_forced=0;
#ifdef NEW_UNITIGGER_INTERFACE
#if 0
    fprintf(stderr, "i=%d contained=%d align_to=%d num_frags=%d\n", 
            i, afrag->contained, align_to, num_frags);
#endif
    // Don't process the last fragment unless it is contained!
    if ((i==num_frags-1) && !afrag->contained)
      continue;
    align_to = GetFragmentIndex(positions[i].ident2, positions, num_frags);
    if (align_to < 0)
      continue;
    assert(align_to >= 0);
    ahang = positions[i].ahang;
    is_pointed[align_to] = 1;
    if (align_to < i)
      {
        // Redefine afrag and bfrag: afrag should be the upstream one!
        bfrag = afrag;
        afrag = GetFragment(fragmentStore, align_to);
      }
    else
      bfrag = GetFragment(fragmentStore, align_to);
#if 0
    fprintf(stderr, 
            "i= %d align_to= %d ident= %d ident2= %d a_contained=%d b_contained= %d ahang= %d bhang= %d num_frags=%d\n", 
            i, align_to, positions[i].ident, positions[i].ident2, afrag->contained, 
            bfrag->contained, ahang, positions[i].bhang, num_frags);
    fprintf(stderr, "Aligned: a=%d b=%d \n", is_aligned[i], is_aligned[align_to]);
#endif
#ifdef NEW_UNITIGGER_INTERFACE
    if (!is_aligned[i] && !is_aligned[align_to])
#endif
      {
        fprintf(stderr, "MultiAlignUnitig()-- Failed to process unitig %d: ", unitig->iaccession);
        fprintf(stderr, "MultiAlignUnitig()-- none of fragments a (i=%d, id=%d) and b (i=%d, id=%d) is pre-aligned\n",
                i, positions[i].ident, align_to, positions[align_to].ident);
        DeleteHashTable_AS(fragmentMap);
        DeleteMANode(ma->lid);
        return EXIT_FAILURE;        
      }
#else
    align_to = i-1;
    while (! olap_success) 
#endif
      {
#ifndef NEW_UNITIGGER_INTERFACE
        if (align_to < 0)
          {
            if (VERBOSE_MULTIALIGN_OUTPUT) {
              fprintf(stderr, "MultiAlignUnitig()-- hit the beginning of fragment list: no fragment upstream overlaps with current fragment %d\n",
                      bfrag->iid);
            }
            break;
          }
        afrag = GetFragment(fragmentStore, align_to);

        /* If bfrag is contained, then afrag should be its container.
         * If bfrag is not contained, then afrag should be just uncontained
         */
        if ( bfrag->contained ) {
          while ( align_to>-1 ) {
            if ( afrag->iid == bfrag->contained && afrag->contained != afrag->iid)
              break;
            align_to--;
            if ( align_to > -1)
              afrag = GetFragment(fragmentStore, align_to);
          }
        } else {
          while ( align_to>0 && afrag->contained ) {
            align_to--;
            if (align_to > -1)
              afrag = GetFragment(fragmentStore, align_to);
          }
        }
        if ( align_to < 0 ) {
          if (VERBOSE_MULTIALIGN_OUTPUT) {
            if (bfrag->contained )
              {
                fprintf(stderr, "MultiAlignUnitig()-- bfrag %d is contained, but no container is found upstream\n",
                        bfrag->iid);
              }
            else
              {
                fprintf(stderr, "MultiAlignUnitig()-- bfrag %d is not contained, and no uncontained afrag is found upstream\n",
                        bfrag->iid);
              }
          }
          break;
        }
        ahang = offsets[bfrag->lid].bgn - offsets[afrag->lid].bgn;

        /* Make sure ahang is above the cutoff value. 
         * If it's not, may need to sort fragments begfore processing
         */
        if (ahang < CNS_NEG_AHANG_CUTOFF && (! allow_neg_hang))
          {
            if (VERBOSE_MULTIALIGN_OUTPUT) {
              fprintf(stderr, "MultiAlignUnitig()-- too negative ahang is detected for afrag %d and bfrag %d; proceed to the next upstraem afrag\n",
                      afrag->iid,  bfrag->iid);
            }
            align_to--;
            continue;
          }
#endif 

        ovl = offsets[afrag->lid].end - offsets[bfrag->lid].bgn;

#if 0
        fprintf(stderr, "MultiAlignUnitig()-- Aligning frag #%d (iid %d, range %d,%d) to afrag iid %d range %d,%d -- ovl=%d ahang=%d\n", 
                positions[i].ident,
                bfrag->iid,
                offsets[bfrag->lid].bgn,
                offsets[bfrag->lid].end,
                afrag->iid,
                offsets[afrag->lid].bgn,
                offsets[afrag->lid].end,
                ovl,
                ahang);
#endif

        if (ovl < 0)
          {
            if (VERBOSE_MULTIALIGN_OUTPUT) {
              fprintf(stderr, "MultiAlignUnitig()-- positions of afrag %d and bfrag %d do not overlap; proceed to the next upstream afrag\n",
                      afrag->iid, bfrag->iid);
            }
            align_to--;
            continue;
          }
            
#ifdef ALIGN_TO_CONSENSUS
        { 
          Bead *afirst = GetBead(beadStore,afrag->beads);
          Column *col = GetColumn(columnStore,afirst->column_index);
          aoffset = col->call;
          // starting consensus bead  is call of columne where first a basepair 
          // is aligned
          olap_success = GetAlignmentTrace(-1, aoffset, bfrag->lid, &ahang, 
                                           ovl, trace, &otype,COMPARE_FUNC,DONT_SHOW_OLAP,0);
        }
#else
        olap_success = GetAlignmentTrace(afrag->lid, 0, bfrag->lid, &ahang, 
                                         ovl, trace, &otype,DP_Compare,DONT_SHOW_OLAP,0);

        if ( ! olap_success && COMPARE_FUNC != DP_Compare ) {
          // try again, perhaps with alternate overlapper
          olap_success = GetAlignmentTrace(afrag->lid, 0, bfrag->lid, &ahang, 
                                           ovl, trace, &otype,COMPARE_FUNC,DONT_SHOW_OLAP,0);
        }
#endif
        if ( !olap_success ) {
#ifndef NEW_UNITIGGER_INTERFACE       
          align_to--;
#endif
          fprintf(stderr, "MultiAlignUnitig()-- MultiAlignUnitig: positions of %d (%c) and %d (%c) overlap, but GetAlignmentTrace returns no overlap success estimated ahang: %d%s\n",
                  afrag->iid,afrag->type,bfrag->iid,bfrag->type,
                  ahang, (bfrag->contained)?" (reported contained)":"");
        }
      } /* ! olap_success */

    if ( ! olap_success ) {
      if (VERBOSE_MULTIALIGN_OUTPUT) {
        if ( bfrag->contained && afrag->iid != bfrag->contained ) {
          // report a more meaningful error in the case were overlap with
          //   a declared contained isn't successful
          fprintf(stderr, "MultiAlignUnitig()-- Could not find overlap between bfrag %d (%c) and its containing fragment, %d.\n",
                  bfrag->iid,bfrag->type, bfrag->contained);
        } else {
          fprintf(stderr, "MultiAlignUnitig()-- Could (really) not find overlap between afrag %d (%c) and bfrag %d (%c); estimated ahang: %d\n",
                  afrag->iid,afrag->type,bfrag->iid,bfrag->type,ahang);
        }

        //PrintFrags(stderr,0,&positions[i],1,gkpStore);
      }

      if ( allow_forced_frags ) {
        frag_forced = 1;
        unitig_forced = 1;
      } else {
        DeleteHashTable_AS(fragmentMap);
        DeleteMANode(ma->lid);
        return EXIT_FAILURE;
      }
    }
    if ( mark_contains && otype == AS_CONTAINMENT ) { 
      MarkAsContained(i);
    }
#ifdef NEW_UNITIGGER_INTERFACE
    i_afrag = (i<align_to) ? i : align_to;
    if (i_afrag && !is_pointed[i] && !afrag->contained && !bfrag->contained)
      {
        fprintf(stderr, "Failed to process unitig %d: no fragment pointed to uncontained fragment %d.\n",
                unitig->iaccession, positions[i].ident);
        DeleteHashTable_AS(fragmentMap);
        DeleteMANode(ma->lid);
        return EXIT_FAILURE;         
      }
#endif
    if ( frag_forced ) {
      ApplyAlignment(afrag->lid,0,bfrag->lid,ahang,Getint32(trace,0));
    } else {
#ifdef ALIGN_TO_CONSENSUS
      ApplyAlignment(-1,aoffset,bfrag->lid,ahang,Getint32(trace,0));
#else
      ApplyAlignment(afrag->lid,0,bfrag->lid,ahang,Getint32(trace,0));
#endif
    }
#ifdef NEW_UNITIGGER_INTERFACE
    is_aligned[i       ] = 1;
    is_aligned[align_to] = 1;
#endif
  } /* loop through all the unitigs */
#ifdef NEW_UNITIGGER_INTERFACE
  safe_free(is_pointed);
  safe_free(is_aligned);
#endif

  RefreshMANode(ma->lid, 0, opp, NULL, NULL, 0, 0);

  safe_free(offsets);

  if (printwhat == CNS_VERBOSE) 
    PrintAlignment(stderr,ma->lid,0,-1,printwhat);
  if ( ! unitig_forced ) {
    score_reduction = AbacusRefine(ma,0,-1,CNS_SMOOTH, opp);
    MergeRefine(ma->lid, NULL, NULL, opp, 0);
    AbacusRefine(ma,0,-1,CNS_POLYX, opp);
    MergeRefine(ma->lid, NULL, NULL, opp, 0);
    if (printwhat == CNS_VERBOSE)
      PrintAlignment(stderr,ma->lid,0,-1,printwhat);
    AbacusRefine(ma,0,-1,CNS_INDEL, opp);
    MergeRefine(ma->lid, NULL, NULL, opp, 1);

    if (printwhat != CNS_QUIET && printwhat != CNS_STATS_ONLY)
        PrintAlignment(stderr,ma->lid,0,-1,printwhat);
  }
  GetMANodeConsensus(ma->lid,sequence,quality);
  unitig->consensus = Getchar(sequence,0);
  unitig->quality = Getchar(quality,0);
  GetMANodePositions(ma->lid, num_frags,unitig->f_list, 0,NULL, deltas);
  unitig->length = GetNumchars(sequence)-1;
  if ( do_rez) {
    char **multia = NULL;
    int **id_array = NULL;
    int depth;
    int i;
    int rc;
    char srcadd[32];
    int addlen;
    double prob_value=0;

    rc = MANode2Array(ma, &depth, &multia, &id_array,0);
    if ( rc ) {
      prob_value = AS_REZ_MP_MicroHet_prob(multia,
                                           id_array,
                                           gkpStore,
                                           unitig->length,
                                           depth);
    } else {
      prob_value = 0;
    }

    addlen = sprintf(srcadd,"\nmhp:%e",prob_value); 
#ifdef AS_ENABLE_SOURCE
    if ( unitig->source != NULL ) {
      memcpy(&SRCBUFFER[0],unitig->source,strlen(unitig->source)+1);
      strcat(&SRCBUFFER[0],srcadd);
      unitig->source = &SRCBUFFER[0];
    } else {
      memcpy(&SRCBUFFER[0],srcadd,addlen+1);
      unitig->source = &SRCBUFFER[0];
    }
#endif
    if ( rc ) {
      for (i=0;i<depth;i++) {
        safe_free(multia[2*i]);
        safe_free(multia[2*i+1]);
        safe_free(id_array[i]);
      }
      safe_free(multia);
      safe_free(id_array);
    }
  }

  DeleteHashTable_AS(fragmentMap);
  DeleteMANode(ma->lid);
  return EXIT_SUCCESS;
}

int IsDovetail(SeqInterval a,SeqInterval b) {
  int ahang;
  int alen, blen;
  int acomplement=0,bcomplement=0;
  alen = a.end - a.bgn;
  blen = b.end - b.bgn;
  if ( alen < 0 ) {
    acomplement++;
    alen = -alen; 
  }
  if ( blen < 0 ) {
    bcomplement++;
    blen = -blen;
  }
  if ( acomplement && bcomplement) {
    ahang = b.end - a.end;
  } else if ( acomplement && !bcomplement ) {
    ahang = b.bgn - a.end;
  } else if ( ! acomplement && bcomplement ) {
    ahang = b.end -a.bgn; 
  } else {
    ahang = b.bgn - a.bgn;
  }
  if ( ahang >= alen ) return -1;
  return ahang;
}

int32 PlaceFragments(int32 fid, Overlap *(*COMPARE_FUNC)(COMPARE_ARGS),
                     CNS_Options *opp) 
{
  /*
    all of fid's component frags will be aligned to it
    (not to eachother)
          
    fcomplement==0                                fcomplement==1
       
    A)       fid                                  C)     fid 
    ------------------>                            <----------------
    --->                                                        <---
    bid (bcomplement==0)                                       bid
       
    B)       fid                                  D)     fid
    ------------------>                            <----------------
    <---                                                        --->
    bid (bcomplement==1)                                       bid
       
  */
  int n_frags,i,ahang,ovl,fcomplement,bcomplement;
  int32 blid,afid = NULLINDEX;
  OverlapType otype;
  static VA_TYPE(int32) *trace = NULL;
  Fragment *afrag=GetFragment(fragmentStore,fid);
  CNS_AlignedContigElement *bfrag=GetCNS_AlignedContigElement(fragment_positions,afrag->components);
  MultiAlignT *ma;

  fcomplement = afrag->complement;
  n_frags = 0;
  if ( afrag->n_components == 0 ) return 0;
  if ( trace== NULL ) {
    trace = CreateVA_int32(AS_READ_MAX_LEN);
  } else {
    ResetVA_int32(trace);
  }
  if ( USE_SDB ) {
    assert(sequenceDB != NULL);
    ma = loadMultiAlignTFromSequenceDB(sequenceDB, afrag->iid, TRUE);
  } else {
    ma = GetMultiAlignInStore(unitigStore,afrag->iid);
  }
  for (i=0;bfrag->frg_or_utg==CNS_ELEMENT_IS_FRAGMENT;i++,bfrag++,n_frags++) 
    {
      int align_failure=0;
      int containFound=0;

      if (!ExistsInHashTable_AS(fragmentMap, bfrag->idx.fragment.frgIdent, 0))
        continue;

      bcomplement = (bfrag->position.bgn < bfrag->position.end) ? 0 : 1;

      // next test to to see whether IUM's fragment is in the ICM
      //if (  Getint32(fragment_indices,bfrag->ident) == NULL  ) continue;
      //if ( ! *Getint32(fragment_indices,bfrag->ident)) continue;

      blid = AppendFragToLocalStore(bfrag->idx.fragment.frgType,
                                    bfrag->idx.fragment.frgIdent,
                                    (bcomplement != fcomplement),
                                    bfrag->idx.fragment.frgContained,
                                    AS_OTHER_UNITIG, NULL);

      afrag = GetFragment(fragmentStore,fid); 

#if 0
      {
        Fragment *tfrag=GetFragment(fragmentStore,blid);
        ReplaceInHashTable_AS(fragmentMap, tfrag->iid, 0, tfrag->lid, 0);
      }
#endif

      if ( bfrag->idx.fragment.frgContained > 0 ) {
        if (!ExistsInHashTable_AS(fragmentMap, bfrag->idx.fragment.frgContained, 0)) {

#define ALLOW_MISSING_CONTAINER_TO_HANDLE_SURROGATE_RESOLUTION
#ifndef ALLOW_MISSING_CONTAINER_TO_HANDLE_SURROGATE_RESOLUTION
          fprintf(stderr, "Could not find containing fragment %d in local store\n",
                  bfrag->idx.fragment.frgContained);
          return EXIT_FAILURE;
#else
          fprintf(stderr, "Could not find containing fragment %d in local store -- due to surrogate resolution?\n",
                  bfrag->idx.fragment.frgContained);
#endif
        } else {
          containFound=1;
          afid = fid;
        }
      }

      if(!containFound){ // either not contained or container not found
        // afrag = GetFragment(fragmentStore,blid-1); 
        // if ( afrag->type == AS_UNITIG ) {
        //   afrag = GetFragment(fragmentStore,fid); 
        // }
        afid = fid;
      }
      afrag = GetFragment(fragmentStore,afid); 
      if ( afrag == NULL ) {
        fprintf(stderr,
                "Lookup failure in CNS: attempting to align %d with %d, but aligned frag %d could not be found\n",
                bfrag->idx.fragment.frgIdent,afid,afid);
        assert(0);
      }
      ovl = GetFragment(fragmentStore,blid)->length;
      if ( fcomplement && bcomplement) {
        ahang = afrag->length - bfrag->position.bgn; /* Case D */
      } else if ( fcomplement && !bcomplement ) {
        ahang = afrag->length - bfrag->position.end; /* Case C */
      } else if ( ! fcomplement && bcomplement ) {
        ahang = bfrag->position.end;                 /* Case B */
      } else {
        ahang = bfrag->position.bgn;                 /* Case A */
      }

      if ( ! GetAlignmentTrace(afrag->lid, 0,blid, &ahang, ovl, trace, &otype, 
                               DP_Compare,DONT_SHOW_OLAP,0)  
           &&
           ! GetAlignmentTrace(afrag->lid, 0,blid, &ahang, ovl, trace, &otype, 
                               COMPARE_FUNC,SHOW_OLAP,0)      ) 
        {
          Bead *afirst = GetBead(beadStore,afrag->beads+ahang);
          Column *col = GetColumn(columnStore,afirst->column_index);
          MANode *manode = GetMANode(manodeStore,col->ma_id);
          int i;
          {
            IntMultiVar *vl = NULL;
            int32        nv  = 0;
            RefreshMANode(manode->lid, 0, opp, &nv, &vl, 0, 0);
          }
          afirst = GetBead(beadStore,afrag->beads+ahang);
          col = GetColumn(columnStore,afirst->column_index);

          fprintf(stderr, "Could (really) not find overlap between %d (%c) and %d (%c) estimated ahang: %d\n",
                  afrag->iid,
                  afrag->type,
                  bfrag->idx.fragment.frgIdent,
                  bfrag->idx.fragment.frgType,
                  ahang);
          fprintf(stderr, "Ejecting fragment %d from contig\n",
                  bfrag->idx.fragment.frgIdent);

          //PrintAlignment(cnslog, manode->lid, col->ma_index, col->ma_index+500,CNS_CONSENSUS);
          // From here, try some previously aligned fragments
          //CleanExit("",__LINE__,1);
          // okay, try to pitch this fragment
          { Fragment *b = GetFragment(fragmentStore,blid);
            b->deleted = 1;
            align_failure = 1;
          }
        }
      if ( ! align_failure ) {
        ApplyAlignment(afrag->lid,0,blid,ahang,Getint32(trace,0));
      }
    }
  return n_frags;
}


int MultiAlignContig(IntConConMesg *contig,
                     VA_TYPE(char) *sequence, VA_TYPE(char) *quality, 
                     VA_TYPE(int32) *deltas, CNS_PrintKey printwhat, 
                     Overlap *(*COMPARE_FUNC)(COMPARE_ARGS), CNS_Options *opp)   
{
  MANode *ma;
  int           num_unitigs,num_frags;
  int32         num_columns=0;
  int           complement;
  int           forced_contig=0;
  int32         fid,i,align_to;
  IntUnitigPos *upositions; 
  SeqInterval  *offsets;
  int           total_aligned_elements=0;
  static        VA_TYPE(int32) *trace=NULL;

  if (contig == NULL ) 
    {
      fprintf(stderr, "Null contig #%d detected\n", contig->iaccession);
      return EXIT_FAILURE;
    }

  num_unitigs = contig->num_unitigs;
  num_frags = contig->num_pieces;
  upositions = contig->unitigs;
  total_aligned_elements=num_frags+num_unitigs;

  RALPH_INIT = InitializeAlphTable();

  offsets = (SeqInterval *) safe_calloc(num_unitigs,sizeof(SeqInterval));
  for (i=0;i<num_unitigs;i++) {
    num_columns = ( upositions[i].position.bgn>num_columns)? upositions[i].position.bgn : num_columns;
    num_columns = ( upositions[i].position.end>num_columns)? upositions[i].position.end : num_columns;
  }

  ResetStores(num_unitigs,num_columns);

  {
    int hash_rc;

    fragmentMap = CreateScalarHashTable_AS(2*(num_frags+num_unitigs));
    for (i=0;i<num_frags;i++) {
      if (ExistsInHashTable_AS (fragmentMap, contig->pieces[i].ident, 0)) {
        // indicates that the fragment appears more than once in the f_list;
        fprintf(stderr,"Failure to insert ident %d in fragment hashtable, already present\n",contig->pieces[i].ident); 
        assert(0);
      }
      InsertInHashTable_AS(fragmentMap, contig->pieces[i].ident, 0, 1, 0);
    }
    contig_id = contig->iaccession;
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
        Fragment *f;
        char     *s;

        f = GetFragment(fragmentStore,i); 
        s = Getchar(sequenceStore,f->sequence);
        fprintf(stderr, ">unitig-%d\n%s\n", f->iid, s);
      }
    }

    ma = CreateMANode(contig->iaccession);
    if ( trace == NULL ) {
      trace = CreateVA_int32(AS_READ_MAX_LEN);
    } else {
      ResetVA_int32(trace);
    }
    
    // Seed multiAlignment with 1st fragment of 1st unitig 
    SeedMAWithFragment(ma->lid, GetFragment(fragmentStore,0)->lid,0, opp);
    PlaceFragments(GetFragment(fragmentStore,0)->lid,COMPARE_FUNC, opp);
     
    // Now, loop on remaining fragments, aligning to:
    //    a)  containing frag (if contained)
    // or b)  previously aligned frag
    for (i=1;i<num_unitigs;i++) 
      {
        int ahang,ovl;
        int32 alid,blid;
        int32 last_b_aligned;
        OverlapType otype;
        Fragment *afrag = NULL;
        Fragment *bfrag = GetFragment(fragmentStore,i); 
        int olap_success=0;
        int try_contained=0;
        Fragment *afrag_first = NULL;
        int ahang_first = -1;
        blid = bfrag->lid;
        // check whether contained, if so
        // align_to = containing
        // else 
        align_to = i-1;
        while (! olap_success) 
          {
            while ( align_to > 0 && ( (try_contained)?0:IsContained(align_to)) ) {
              align_to--;
            }
            if ( align_to < 0 ) {
              if (VERBOSE_MULTIALIGN_OUTPUT) {
                fprintf(stderr, "MultiAlignContig: hit the beginning of ");
                fprintf(stderr, "unitig list: no unitig upstream ");
                fprintf(stderr, "overlaps with current unitig %d\n",
                        bfrag->iid);
              }

              if ( align_to < 0 && ! try_contained ) {
                try_contained = 1;
                if (VERBOSE_MULTIALIGN_OUTPUT) {
                  fprintf(stderr, "MultiAligncontig: trying contained afrags ");
                  fprintf(stderr, "for bfrag %d\n", bfrag->iid);
                }
                align_to = i-1;
                continue;
              }

              break;
            }
            afrag = GetFragment(fragmentStore, align_to);
            alid = afrag->lid;
            ovl = offsets[alid].end - offsets[blid].bgn;
            if (ovl <= 0)
              {
                if (VERBOSE_MULTIALIGN_OUTPUT) {
                  fprintf(stderr, "MultiAlignContig: positions of afrag ");
                  fprintf(stderr, "%d and bfrag %d do not overlvreg. ",
                          afrag->iid, bfrag->iid);
                  fprintf(stderr, "Proceed to the next upstream afrag\n");
                }
                align_to--;
                continue;
              }
            else /* ovl > 0 */
              {
                ahang = offsets[blid].bgn - offsets[alid].bgn;

                if ( ahang_first == -1) {
                  ahang_first = ahang;
                  afrag_first = afrag;
                }
                ResetVA_int32(trace);

                if (VERBOSE_MULTIALIGN_OUTPUT)
                  fprintf(stderr,"Attemping alignment (DP_Compare) of afrag %d (%c) and bfrag %d (%c) estimated ahang: %d\n",
                          afrag->iid,afrag->type,bfrag->iid,bfrag->type,ahang);
                olap_success = GetAlignmentTrace(afrag->lid, 0,bfrag->lid, &ahang, ovl,
                                                 trace, &otype, DP_Compare,DONT_SHOW_OLAP,0);

                if ( !olap_success && COMPARE_FUNC != DP_Compare ) {
                  if (VERBOSE_MULTIALIGN_OUTPUT)
                    fprintf(stderr,"Attemping alignment (COMPARE_FUNC) of afrag %d (%c) and bfrag %d (%c) estimated ahang: %d\n",
                            afrag->iid,afrag->type,bfrag->iid,bfrag->type,ahang);
                  olap_success = GetAlignmentTrace(afrag->lid, 0,bfrag->lid, &ahang, ovl,
                                                   trace, &otype, COMPARE_FUNC,SHOW_OLAP,0);
                }

                // here, calculate the appropriate allowable endgvreg.
                //  ------------ afrag
                //     ----------- bfrag
                //         ------------------ nextfrag
                //  offsets[nextfrag].bgn - offsets[bfrag].bgn
             
                if ( !olap_success && COMPARE_FUNC != DP_Compare ) {
                  int max_gap=0;
                  int nlid=blid+1;
                  if ( nlid < num_unitigs )
                    max_gap = offsets[nlid].bgn - offsets[blid].bgn;
                  else
                    max_gap = 800;

                  if (VERBOSE_MULTIALIGN_OUTPUT)
                    fprintf(stderr, "Trying local aligner on unitigs %d and %d, allowing for large endgaps (up to %d)\n",
                            afrag->iid, bfrag->iid, max_gap);
                  olap_success = GetAlignmentTrace(afrag->lid, 0,bfrag->lid, &ahang, ovl,
                                                   trace, &otype, COMPARE_FUNC,SHOW_OLAP,max_gap);
                }
              }
            if ( ! olap_success ) {
              if (VERBOSE_MULTIALIGN_OUTPUT) {
                fprintf(stderr, "MultiAlignContig: positions of ");
                fprintf(stderr, "afrag %d (%c) and bfrag %d (%c) overlap, but ",
                        afrag->iid,afrag->type,bfrag->iid,bfrag->type);
                fprintf(stderr, "GetAlignmentTrace returns no overlap success ");
                fprintf(stderr, "estimated ahang: %d\n", ahang);
              }
              align_to--;
              if ( align_to < 0 && ! try_contained ) {
                try_contained = 1;
                if (VERBOSE_MULTIALIGN_OUTPUT) {
                  fprintf(stderr, "MultiAligncontig: trying contained afrags ");
                  fprintf(stderr, "for bfrag %d\n", bfrag->iid);
                }
                align_to = i-1;
              }
            }
          }

        if ( ! olap_success ) {
          fprintf(stderr,"Could (really) not find overlap between %d (%c) and %d (%c), estimated ahang %d\n", 
                  afrag->iid,afrag->type,bfrag->iid,bfrag->type, ahang);

          if (FORCE_UNITIG_ABUT == 0) {
            fprintf(stderr,"You can force these to abut with '-D forceunitigabut'\n");
            safe_free(offsets);
            return EXIT_FAILURE;
          }

          fprintf(stderr,"Forcing the alignment between %d (%c) and %d (%c)\n", 
                  afrag->iid,afrag->type,bfrag->iid,bfrag->type);

          forced_contig = 1; 

          //  BUG! (?)  If there is no overlap between the first and
          //  second frag, we might not have *_first set.
          //
          //  The following happened.  The original placement was:
          //
          //    unitig1    ------------------
          //    unitig2                     -----
          //    unitig3                       --------
          //
          //  BUT, that thin overlap wasn't found between #2 and #1,
          //  instead, we have this:
          //
          //    unitig1    ------------------
          //    unitig2         ------
          //    unitig3                       --------
          //
          //  #3 tried to align to #2, saw it was contained, skipped
          //  it.  There is now no overlap between #3 and #1, and so
          //  we never set the *_first variables.
          //  
          //  We're likely pretty screwed up here already, and since
          //  it took manual intervention to get here (being a
          //  recompilation), we're gonna force it to go together by
          //  picking the closest unitig (#1 in this case) as our
          //  first.
          //
          //  If you do see this again, bpw suggests yanking out
          //  unitig#2 (it's a rock, right?)  and forcing the
          //  alignment.
          //

          if (afrag_first) {
            afrag = afrag_first;
            ahang = ahang_first;
          } else {
            //  Dang, we're really screwed.  Nobody overlapped with us.
            //  Cross our fingers and find the closest end point.
            //
            int   maxOvl = -offsets[blid].bgn;

            if (VERBOSE_MULTIALIGN_OUTPUT)
              fprintf(stderr, "MultiAlignContig:  YIKES!  Your unitig doesn't overlap with anything!  Picking the closest thing!\n");

            align_to = i-1;

            while (align_to >= 0) {
              if ((try_contained) ? 0 : IsContained(align_to)) {
                //  NOP!  Found a contained frag, and we want to skip it.
              } else if (maxOvl < offsets[alid].end - offsets[blid].bgn) {
                afrag  = GetFragment(fragmentStore, align_to);
                alid   = afrag->lid;
                maxOvl = offsets[alid].end - offsets[blid].bgn;
                ahang  = offsets[blid].bgn - offsets[alid].bgn;

                fprintf(stderr, "MultiAlignContig:  RESET align_to=%d alid=%d maxOvl=%d ahang=%d\n", align_to, alid, maxOvl, ahang);
              }

              align_to--;
            }
          }

          if (VERBOSE_MULTIALIGN_OUTPUT)
            fprintf(stderr, "MultiAlignContig:  Forcing abut between afrag %d (%c) and bfrag %d (%c).\n",
                    afrag->iid, afrag->type,
                    bfrag->iid, bfrag->type);

          assert(afrag != NULL);

          //  If our ahang is too big, force a 20bp overlap.
          //
          if (ahang + 20 > afrag->length) {
            fprintf(stderr, "RESET: ahang from %d to %d\n", ahang, afrag->length-20);
            ahang = afrag->length - 20;
          }

          //  And now, after all this work, we simply die in the
          //  ApplyAlignment() ten lines below.  Too bad.....

          otype = AS_DOVETAIL;
        }  //  End of forcing the unitig to abut.

        //  Unitig is placed, or we just forced it to be placed.

        if ( otype == AS_CONTAINMENT ) { 
          MarkAsContained(i);
        }

        last_b_aligned=ApplyAlignment(afrag->lid,0,bfrag->lid,ahang,Getint32(trace,0));

        PlaceFragments(bfrag->lid,COMPARE_FUNC, opp);
      }
    //   contig->num_vars = 20;     // affect .cns/ICM
    {
      IntMultiVar  *vl=NULL;
      int32 nv=0;
      RefreshMANode(ma->lid, 0, opp, &nv, &vl, 0, 0);
    }
    // Now, must find fragments in regions of overlapping unitigs, and adjust 
    // their alignments as needed
     
    // DeleteVA_int32(trace);
    if (printwhat == CNS_VERBOSE) {
      fprintf(stderr,"Initial pairwise induced alignment\n");
      PrintAlignment(stderr,ma->lid,0,-1,printwhat);
    }
#if 0
    fprintf(stderr, "Contig = %lu\n\n", contig->iaccession);
#endif
    AbacusRefine(ma,0,-1,CNS_SMOOTH, opp);
    {
      IntMultiVar  *vl=NULL;
      int32 nv=0;
      int score_reduction;
      MergeRefine(ma->lid, NULL, NULL, opp, 0);
      AbacusRefine(ma,0,-1,CNS_POLYX, opp);
      if (printwhat == CNS_VERBOSE) {
        fprintf(stderr,"\nPOLYX refined alignment\n");
        PrintAlignment(stderr,ma->lid,0,-1,printwhat);
      }
      RefreshMANode(ma->lid, 0, opp, &nv, &vl, 0, 0);
      score_reduction = AbacusRefine(ma,0,-1,CNS_INDEL, opp);
      MergeRefine(ma->lid, &(contig->v_list), &(contig->num_vars), opp, 2);
    }

    if ((printwhat == CNS_VERBOSE || printwhat == CNS_VIEW_UNITIG)) { 
      fprintf(stderr,"\nFinal refined alignment\n");
      PrintAlignment(stderr,ma->lid,0,-1,printwhat);
    }
    if ( num_frags == 0 ) {
      PrintAlignment(stderr,ma->lid,0,-1,printwhat);
    }
    GetMANodeConsensus(ma->lid,sequence,quality);
    contig->consensus = Getchar(sequence,0);
    contig->quality  = Getchar(quality,0);
    contig->num_pieces = GetMANodePositions(ma->lid, num_frags, contig->pieces, 
                                            num_unitigs, contig->unitigs, deltas);
    contig->length = GetNumchars(sequence)-1;
    contig->forced = forced_contig;

    DeleteMANode(ma->lid);
    DeleteHashTable_AS(fragmentMap);
  }

  safe_free(offsets);
  return EXIT_SUCCESS; 
}



int MultiAlignContig_ReBasecall(MultiAlignT *cma, VA_TYPE(char) *sequence, VA_TYPE(char) *quality, 
                                CNS_Options *opp) 
{
  MANode *ma; // this is to build, for purposes of ascii printout or analysis
  MultiAlignStoreT *contigStore;
  int contigID=cma->maID;
  int num_unitigs,num_frags;
  int32 num_columns=0;
  int32 fid,i;
  IntMultiPos *fpositions; 
  IntUnitigPos *upositions; 
  IntMultiVar *vl;
  int32 nv;

  // static VA_TYPE(int32) *trace=NULL;
  static int32 *tracep=NULL;
   
  num_frags=GetNumIntMultiPoss(cma->f_list);
  num_unitigs=GetNumIntUnitigPoss(cma->u_list);
  fpositions=GetIntMultiPos(cma->f_list,0);
  upositions=GetIntUnitigPos(cma->u_list,0);

  RALPH_INIT = InitializeAlphTable();

  if ( tracep == NULL ) {
    tracep = safe_malloc(sizeof(int32)*(AS_READ_MAX_LEN+1));
  } 

  ResetStores(num_unitigs,num_columns);
  contigStore = CreateMultiAlignStoreT();
  SetMultiAlignInStore(contigStore,cma->maID,cma);

  ma = CreateMANode(contigID);
  fid = AppendFragToLocalStore(AS_CONTIG,
                               contigID,
                               0,
                               0,
                               AS_OTHER_UNITIG, contigStore);
  SeedMAWithFragment(ma->lid, GetFragment(fragmentStore,0)->lid,-1, opp);
     
  // Now, loop on the fragments, applying the computed alignment from the InMultiPos:
  for (i=0;i<num_frags;i++) {
    IntMultiPos *imp=fpositions +i;
    int ahang;
    int32 blid;
    Fragment *afrag=GetFragment(fragmentStore,0); // always align to the contig consensus
    int fcomplement=(imp->position.bgn<imp->position.end)?0:1;
        
    ahang=(fcomplement)?imp->position.end:imp->position.bgn;
    blid = AppendFragToLocalStore(imp->type,
                                  imp->ident,
                                  fcomplement,
                                  0,
                                  AS_OTHER_UNITIG, NULL);
    assert ( imp->delta_length < AS_READ_MAX_LEN );
    memcpy(tracep,imp->delta,imp->delta_length*sizeof(int32));
    tracep[imp->delta_length]=0;
    ApplyIMPAlignment(afrag->lid,blid,ahang,tracep);
  }
  for (i=0;i<num_unitigs;i++) {
    IntUnitigPos *iup=upositions +i;
    int ahang;
    int32 blid;
    Fragment *afrag=GetFragment(fragmentStore,0); // always align to the contig consensus
    int fcomplement=(iup->position.bgn<iup->position.end)?0:1;
        
    ahang=(fcomplement)?iup->position.end:iup->position.bgn;
    blid = AppendFragToLocalStore(AS_UNITIG,
                                  iup->ident,
                                  fcomplement,
                                  0,
                                  iup->type, unitigStore);
    assert ( iup->delta_length < AS_READ_MAX_LEN );
    memcpy(tracep,iup->delta,iup->delta_length*sizeof(int32));
    tracep[iup->delta_length]=0;
    ApplyIMPAlignment(afrag->lid,blid,ahang,tracep);
  }
  RefreshMANode(ma->lid, 1, opp, &nv, &vl, 1, 0);
  GetMANodeConsensus(ma->lid,sequence,quality);
  // DeleteVA_int32(trace);
  DeleteMANode(ma->lid);
  if ( contigStore ) DeleteMultiAlignStoreT(contigStore);
  return 0; 
}



int ExamineMANode(FILE *outFile,int32 sid, int32 mid, UnitigData *tigData,int num_unitigs,
                  CNS_Options *opp) 
{
  int index=0,ugindex=0;
  int32 cid;
  Column *column;
  Bead *cbead;
  int unitig_index=0;
  int tindex=0;
  UnitigData *tig;
  MANode *ma = GetMANode(manodeStore,mid);
  VarRegion  vreg;
  char base;

  SetDefault(&vreg);
  assert(ma != NULL);
  if ( ma->first == -1 ) return 1;
  cid = ma->first;
  while ( cid  > -1 ) {
    char base;
    char qv;
    int tig_depth=0;
    double var;

    column = GetColumn(columnStore,cid);
    assert(column != NULL);
    cbead = GetBead(beadStore,column->call); 
    base = *Getchar(sequenceStore,cbead->soffset);
    qv = *Getchar(qualityStore,cbead->soffset);
    fprintf(outFile,"%d\t%d\t%d\t%d\t%c\t%c\t" ,sid,ma->iid,index,ugindex,base,qv);
    ShowBaseCountPlain(outFile,&column->base_count);
    BaseCall(cid, 1, &var, &vreg, vreg.alleles[0].id, &base, 0, 0, opp); 
    // recall with quality on (and QV parameters set by user)
    fprintf(outFile,"%c\t%c\t", *Getchar(sequenceStore,cbead->soffset), 
            *Getchar(qualityStore,cbead->soffset));
    // restore original consensus basecall/quality
    Setchar(sequenceStore,cbead->soffset,&base);
    Setchar(qualityStore,cbead->soffset,&qv);
    tig=&tigData[unitig_index];
    while ( index >= tig->right && unitig_index < num_unitigs-1) { 
      unitig_index++; 
      tig++;
    }
    tindex=unitig_index;
    while ( tindex < num_unitigs && index  >= tig->left && index < tig->right ) {
      tig_depth++;
      fprintf(outFile,"%d\t%c\t%f\t%d\t",tig->ident,tig->type,tig->coverage_stat, tig->length);
      tindex++;
      tig++;
    }
   
    fprintf(outFile,"\n");
    if ( *Getchar(sequenceStore,cbead->soffset) != '-') ugindex++;
    index++;
    cid = column->next;
  }
  return 1;
}

// static int utl_counts[4]={0,0,0,0};


int ExamineConfirmedMMColumns(FILE *outFile,int32 sid, int32 mid, UnitigData *tigData,int num_unitigs) {
  int index=0,ugindex=0;
  int32 cid;
  Column *column;
  Column *last_mm=NULL;
  Column *frag_start_column;
  static VA_TYPE(Bead) *shared_left=NULL;
  static VA_TYPE(Bead) *shared_right=NULL;
  HashTable_AS *bhash=NULL;
  MANode *ma = GetMANode(manodeStore,mid);
  
  assert(ma != NULL);
  if ( ma->first == -1 ) return 1;
  if ( bhash==NULL ) bhash = CreateScalarHashTable_AS(5000);
  if ( shared_left== NULL ) {
    shared_left = CreateVA_Bead(100);
    shared_right = CreateVA_Bead(100);
  } else {
    ResetVA_Bead(shared_left);
    ResetVA_Bead(shared_right);
  }
  cid = ma->first;
  while ( cid  > -1 ) {
    Fragment *frag;
    Bead *cbead;
    Bead *fbead = NULL;
    char base;
    char qv;
    int bid;
    int depth=0;
    column = GetColumn(columnStore,cid);
    assert(column != NULL);
    cbead = GetBead(beadStore,column->call); 
    base = *Getchar(sequenceStore,cbead->soffset);
    qv = *Getchar(qualityStore,cbead->soffset);
    depth=GetDepth(column);
    // check to see whether there is a confirmed mismatch
    if ( depth > GetBaseCount(&column->base_count,base)+1) {
      // potential for a confirmed mismatch
      char mm=GetConfMM(&column->base_count,base);
      if ( mm != base ) { // this condition indicates a "positive" return from the preceding GetConfMM call 
        if ( last_mm == NULL ) {
          last_mm=column;
        } else {
          //check for compatibility with last confirmed mismatch
          ColumnBeadIterator bi;
          ResetVA_Bead(shared_left);
          ResetVA_Bead(shared_right);
          ResetHashTable_AS(bhash);
          if(!CreateColumnBeadIterator(column->lid,&bi)){
            fprintf(stderr, "CreateAbacus CreateColumnBeadIterator failed");
            assert(0);
          }
          while ( (bid = NextColumnBead(&bi)) != -1 ) {
            cbead = GetBead(beadStore,bid);
            frag = GetFragment(fragmentStore,fbead->frag_index);
            fbead = GetBead(beadStore,frag->beads);
            frag_start_column = GetColumn(columnStore,fbead->column_index);
            if ( column ->ma_index <= last_mm->ma_index ) {
              // shared frag; 
              InsertInHashTable_AS(bhash, frag->lid, 0, bid, 0);
            }
          }
          if ( GetNumBeads(shared_right) > 0 ) {
            if(!CreateColumnBeadIterator(last_mm->lid,&bi)){
              fprintf(stderr, "CreateAbacus CreateColumnBeadIterator failed");
              assert(0);
            }
            while ( (bid = NextColumnBead(&bi)) != -1 ) {
              cbead = GetBead(beadStore,bid);
              frag = GetFragment(fragmentStore,cbead->frag_index);
              if (ExistsInHashTable_AS(bhash, frag->lid, 0)) {
                AppendVA_Bead(shared_left,cbead); 
                cbead = GetBead(beadStore, LookupValueInHashTable_AS(bhash, frag->lid, 0));
                AppendVA_Bead(shared_right,cbead); 
              }
            } 
          }
          // now, look at the partitions of the shared_left and shared_right, and see whether they conflict
          if ( GetNumBeads(shared_left) > 3 ) {
            ShowColumn(last_mm->lid);
            ShowColumn(column->lid);
          }
        }
      }
    }

    
    if ( *Getchar(sequenceStore,cbead->soffset) != '-') ugindex++;
    index++;
    cid = column->next;
  }
  return 1;
}




MultiAlignT *ReplaceEndUnitigInContig( tSequenceDB *sequenceDBp,
                                       GateKeeperStore *frag_store,
                                       uint32 contig_iid, uint32 unitig_iid, int extendingLeft,
                                       Overlap *(*COMPARE_FUNC)(COMPARE_ARGS),
                                       CNS_Options *opp){
  int32 cid,tid; // local id of contig (cid), and unitig(tid)
  int32 aid,bid;  
  int i,num_unitigs;
  MultiAlignT *oma;
  MultiAlignT *cma;
  IntUnitigPos *u_list;
  IntMultiPos *f_list;
  IntMultiVar  *v_list;
  int append_left=0;
  int num_frags=0;
  int complement=0;
  MANode *ma;
  Fragment *cfrag; 
  Fragment *tfrag = NULL;
  static VA_TYPE(int32) *trace=NULL;

  //ALIGNMENT_CONTEXT=AS_CONSENSUS;
  ALIGNMENT_CONTEXT=AS_MERGE;

  //  We need to reset the global sequenceDB pointer -- if we call
  //  this from anything but consensus, the global pointer isn't set.
  //
  sequenceDB = sequenceDBp;

  USE_SDB    = 1;
  RALPH_INIT = InitializeAlphTable();

  gkpStore = frag_store;

  oma =  loadMultiAlignTFromSequenceDB(sequenceDBp, contig_iid, FALSE);

  ResetStores(2,GetNumchars(oma->consensus)+MAX_EXTEND_LENGTH);

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
                               AS_OTHER_UNITIG, NULL);

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
    trace = CreateVA_int32(AS_READ_MAX_LEN);
  ResetVA_int32(trace);

  {
    int ahang,ovl,pos_offset=0;  
    int tigs_adjusted_pos=0;
    OverlapType otype;
    int olap_success=0;
    cfrag=GetFragment(fragmentStore,cid);
    for(i=0;i<num_unitigs;i++) {
      uint32 id=u_list[i].ident;
      if ( id == unitig_iid ) {
        int bgn=u_list[i].position.bgn;
        int end=u_list[i].position.end;
        int complement_tmp=(bgn<end)?0:1;
        int left=(complement_tmp)?end:bgn;
        int right=(complement_tmp)?bgn:end;
        complement=complement_tmp;
        tid = AppendFragToLocalStore(AS_UNITIG,
                                     id,
                                     complement,
                                     0,
                                     AS_OTHER_UNITIG, NULL);
        tfrag=GetFragment(fragmentStore,tid);
        ovl = right-left;  // this is the size of the original (non-extended) unitig
        if ( extendingLeft ) {
          // need to set aid to unitig to preserve positive ahang
          append_left=1;
          aid=tid;
          bid=cid;
          // and ahang estimate is the diff in size between 
          // new unitig (GetFragment(fragmentStore,tid)->length) and old unitig (right-left) 
          ahang = GetFragment(fragmentStore,tid)->length - (right-left);
        }  else {
          aid=cid;
          bid=tid;
          ahang=left;
        }
        SeedMAWithFragment(ma->lid,aid,0, opp);

        //  do the alignment 
        olap_success = GetAlignmentTrace(aid, 0,bid,&ahang,ovl,trace,&otype, DP_Compare,DONT_SHOW_OLAP,0);
        if ( !olap_success && COMPARE_FUNC != DP_Compare )
          olap_success = GetAlignmentTrace(aid, 0,bid,&ahang,ovl,trace,&otype, COMPARE_FUNC,DONT_SHOW_OLAP,0);

        //  If the alignment fails -- usually because the ahang is
        //  negative -- return an empty alignment.  This causes
        //  extendClearRanges (the sole user of this function) to
        //  gracefully handle the failure.
        //
        if (olap_success == 0) {
          return(NULL);
          assert(olap_success);
        }

        ApplyAlignment(aid, 0, bid, ahang, Getint32(trace,0));

        {
          IntMultiVar *vl;
          int32 nv;
          RefreshMANode(ma->lid, 0, opp, &nv, &vl, 0, 0);
        }

        //PrintAlignment(stderr,ma->lid,0,-1,'C');

        break;
      }
    }
  }

  // Now, want to generate a new MultiAlignT which is an appropriate adjustment of original
  cma = CreateMultiAlignT();
  cma->consensus = CreateVA_char(GetMANodeLength(ma->lid)+1);
  cma->quality = CreateVA_char(GetMANodeLength(ma->lid)+1);

  GetMANodeConsensus(ma->lid, cma->consensus, cma->quality);
  // no deltas required at this stage 
  // merge the f_lists and u_lists by cloning and concating
  cma->f_list = Clone_VA(oma->f_list);
  cma->fdelta = CreateVA_int32(0);
  cma->u_list = Clone_VA(oma->u_list);
  cma->udelta = CreateVA_int32(0);
  cma->v_list = Clone_VA(oma->v_list);

  {
    CNS_AlignedContigElement *components;
    CNS_AlignedContigElement *tcomponents;
    CNS_AlignedContigElement *contig_component;
    CNS_AlignedContigElement *aligned_component;
    int ifrag=0;
    int iunitig=0;
    IntMultiPos *imp;
    IntUnitigPos *iup;
    Fragment *frag;
    int ci=0;
    int tc=0; //unitig component index
    int32 bgn,end,left,right,tmp;
    int range_bgn=0,range_end=0,new_tig=0;
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
      left = GetColumn(columnStore, 
                       GetBead(beadStore,frag->beads + left)->column_index)->ma_index;
      right= GetColumn(columnStore, 
                       GetBead(beadStore,frag->beads + right-1)->column_index)->ma_index + 1;
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
        imp->sourceInt = aligned_component->idx.fragment.frgSource;
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

MultiAlignT *MergeMultiAligns( tSequenceDB *, GateKeeperStore *, 
                               VA_TYPE(IntMultiPos) *, int , int , 
                               Overlap *(*COMPARE_FUNC)(COMPARE_ARGS), CNS_Options *);

MultiAlignT *MergeMultiAlignsFast_new( tSequenceDB *sequenceDBp,
                                       GateKeeperStore *frag_store, VA_TYPE(IntElementPos) *positions, 
                                       int quality, int verbose, Overlap *(*COMPARE_FUNC)(COMPARE_ARGS),
                                       CNS_Options *opp)
{
  // this is the functionality used in traditionl CGW contigging
  //     I'm now extending it so that "contained" contigs are handled appropriately,
  //     which is necessitated by "local unitigging" (a.k.a. "meta-unitigging")

  static VA_TYPE(IntMultiPos) *mpositions=NULL;
  static IntMultiPos mpos;
  IntElementPos *epos = GetIntElementPos(positions,0);
  int i;

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

  return(MergeMultiAligns(sequenceDBp, frag_store, mpositions, quality, verbose, COMPARE_FUNC, opp));
}

MultiAlignT *MergeMultiAligns( tSequenceDB *sequenceDBp,
                               GateKeeperStore *frag_store, 
                               VA_TYPE(IntMultiPos) *positions, 
                               int quality, 
                               int verbose, 
                               Overlap *(*COMPARE_FUNC)(COMPARE_ARGS),
                               CNS_Options *opp)
{
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
  ALIGNMENT_CONTEXT = AS_MERGE;
   
  RALPH_INIT = InitializeAlphTable();

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
  } else {
    for (i=0;i<num_contigs;i++) {
      complement = 
        (cpositions[i].position.bgn<cpositions[i].position.end)
        ? 0 : 1;
      fid = AppendFragToLocalStore(cpositions[i].type, 
                                   cpositions[i].ident, 
                                   complement,
                                   0,
                                   AS_OTHER_UNITIG, NULL);
      offsets[fid].bgn = complement?cpositions[i].position.end:cpositions[i].position.bgn;
      offsets[fid].end = complement?cpositions[i].position.bgn:cpositions[i].position.end;
#if 0
      if ( complement ) {
        fprintf(stderr,"%10d:  %12d <---- %12d\n",cpositions[i].ident,offsets[fid].bgn,offsets[fid].end);
      } else {
        fprintf(stderr,"%10d:  %12d ----> %12d\n",cpositions[i].ident,offsets[fid].bgn,offsets[fid].end);
      }
#endif
    }
     
    ma = CreateMANode(cpositions[0].ident);
    if ( trace == NULL ) {
      trace = CreateVA_int32(AS_READ_MAX_LEN);
    } else {
      ResetVA_int32(trace);
    }

    SeedMAWithFragment(ma->lid, GetFragment(fragmentStore,0)->lid,0,opp);
     
    // Now, loop on remaining fragments, aligning to:
    //    a)  containing frag (if contained)
    // or b)  previously aligned frag
    for (i=1;i<num_contigs;i++) 
      {
        int ahang,ovl;
        int32 alid,blid;
        OverlapType otype;
        int olap_success=0;
        int try_contained=0;
        Fragment *afrag = NULL;
        Fragment *bfrag = GetFragment(fragmentStore,i); 
        blid = bfrag->lid;
        // check whether contained, if so
        // align_to = containing
        // else 
        align_to = i-1;

        //  BPW is not sure why this is in a while loop.  If we fail to
        //  find an overlap, the last thing we do is break out of here
        //  (the last "if (!olap_success)").  Feel free to document...

        while (! olap_success) 
          {
            while ( align_to > 0 && IsContained(align_to) ) {
              align_to--;
            }
            if ( align_to < 0 ) {
              if (VERBOSE_MULTIALIGN_OUTPUT) {
                fprintf(stderr, "MergeMultiAligns: unable to find uncontained contig ");
                fprintf(stderr, "upstream from current contig %d\n",
                        bfrag->iid);
              }
              break;
            }
            afrag = GetFragment(fragmentStore, align_to);
            alid = afrag->lid;
            ovl = offsets[alid].end - offsets[blid].bgn;
            if( ovl <= 0 ){
              if (VERBOSE_MULTIALIGN_OUTPUT) {
                fprintf(stderr, "MergeMultiAligns: uncontained contig upstream is found, ");
                fprintf(stderr, "but positions indicate no overlap ");
                fprintf(stderr, "between contigs %d and %d bailing...", afrag->iid, bfrag->iid);
              }
              DeleteMANode(ma->lid);
              safe_free(offsets);
              return NULL;
            }
            if ( offsets[alid].end > offsets[blid].end ) { // containment
              /* GD: this is a containment, assuming that 
               * offsets[alid].beg < offsets[blid].beg
               */ 
              ahang = afrag->length - bfrag->length 
                - (offsets[alid].end-offsets[blid].end);
            } else {
              ahang = afrag->length - ovl;
            }
            olap_success = GetAlignmentTrace(afrag->lid, 0,bfrag->lid, &ahang, ovl, 
                                             trace, &otype,DP_Compare,DONT_SHOW_OLAP,0);
            if ( !olap_success && COMPARE_FUNC != DP_Compare ) {
              olap_success = GetAlignmentTrace(afrag->lid, 0,bfrag->lid, &ahang, ovl, 
                                               trace, &otype,COMPARE_FUNC,SHOW_OLAP,0);
            }
            if ( ! olap_success ) {
              if (VERBOSE_MULTIALIGN_OUTPUT) {
                fprintf(stderr, "MergeMultiAligns: positions of contigs %d and %d ",
                        afrag->iid, bfrag->iid);
                fprintf(stderr, "overlap, but GetAlignmentTrace does not return ");
                fprintf(stderr, "overlap success\n");
              }
              break; 
            }
          }
        if ( ! olap_success ) {
          fprintf(stderr, "MergeMultiAligns failed to find overlap between contigs %d and %d, bailing...\n",
                  afrag->iid, bfrag->iid);
          DeleteMANode(ma->lid);
          safe_free(offsets);
          return NULL;
        }
        if ( otype == AS_CONTAINMENT ) { 
          MarkAsContained(i);
        }
        ApplyAlignment(afrag->lid,0,bfrag->lid,ahang,Getint32(trace,0));
      } /* loop through all contigs */

    {
      IntMultiVar *vl;
      int32 nv;
      RefreshMANode(ma->lid, 0, opp, &nv, &vl, 0, 0);
    }
    // DeleteVA_int32(trace);
  }
  {

    // Now, want to generate a new MultiAlignT which merges the u_list and f_list of the contigs
    // merge the f_lists and u_lists by cloning and concating (or constructing dummy, when dealing with single read

    int ifrag;
    int iunitig;
    IntMultiPos *imp;
    IntUnitigPos *iup;

    cma = CreateMultiAlignT();
    cma->consensus = CreateVA_char(GetMANodeLength(ma->lid)+1);
    cma->quality = CreateVA_char(GetMANodeLength(ma->lid)+1);

    GetMANodeConsensus(ma->lid, cma->consensus, cma->quality);

    // no deltas required at this stage 
    cma->fdelta = CreateVA_int32(0);
    cma->udelta = CreateVA_int32(0);
  
    if( isChunk(cpositions[0].type) ){
      MultiAlignT *ma = loadMultiAlignTFromSequenceDB(sequenceDBp, cpositions[0].ident, cpositions[0].type == AS_UNITIG);
      cma->f_list = Clone_VA(ma->f_list);
      cma->v_list = Clone_VA(ma->v_list);
      cma->u_list = Clone_VA(ma->u_list);
    } else {
      assert(isRead(cpositions[0].type));
      cma->f_list = CreateVA_IntMultiPos(0);
      cma->v_list = CreateVA_IntMultiVar(0);
      cma->u_list = CreateVA_IntUnitigPos(0);
      AppendVA_IntMultiPos(cma->f_list,cpositions+0);
    }

    for (i=1;i<num_contigs;i++) {

      if( isChunk(cpositions[i].type) ){
        MultiAlignT *ma = loadMultiAlignTFromSequenceDB(sequenceDBp, cpositions[i].ident, cpositions[i].type == AS_UNITIG);
        ConcatVA_IntMultiPos(cma->f_list,ma->f_list);
        ConcatVA_IntMultiPos(cma->v_list,ma->v_list);
        ConcatVA_IntUnitigPos(cma->u_list,ma->u_list);
      } else {
        assert(isRead(cpositions[i].type));
        AppendVA_IntMultiPos(cma->f_list,cpositions+i);
      }
    }


    ifrag=0;
    iunitig=0;
    for (i=0;i<num_contigs;i++) 
      {

        Fragment *cfrag=GetFragment(fragmentStore,i);  /* contig pseudo-frag */

        if(isChunk(cfrag->type))
          {

            CNS_AlignedContigElement *components=GetCNS_AlignedContigElement
              (fragment_positions,cfrag->components);
            CNS_AlignedContigElement *compci;

            int ci=0;
            int32 bgn,end,left,right,tmp;
            // make adjustments to positions
            while (ci < cfrag->n_components) { 
              compci = &components[ci];
              if ( cfrag->complement ) {
                bgn = cfrag->length-compci->position.bgn;
                end = cfrag->length-compci->position.end;
              } else {
                bgn = compci->position.bgn;
                end = compci->position.end;
              }
              left = (bgn<end)?bgn:end;
              right = (bgn<end)?end:bgn;
              left = GetColumn(columnStore, 
                               GetBead(beadStore,cfrag->beads + left)
                               ->column_index)->ma_index;
              right   = GetColumn(columnStore, 
                                  GetBead(beadStore,cfrag->beads + right-1)
                                  ->column_index)->ma_index + 1;
              tmp = bgn;
              bgn = (bgn<end)?left:right;
              end = (tmp<end)?right:left;

              if (compci->frg_or_utg==CNS_ELEMENT_IS_UNITIG) {
                iup = GetIntUnitigPos(cma->u_list,iunitig);
                iup->position.bgn = bgn;
                iup->position.end = end;
                iup->delta_length = 0;
                iup->delta = NULL;
                ci++;iunitig++;
              } else {
                imp = GetIntMultiPos(cma->f_list,ifrag);
                imp->ident = compci->idx.fragment.frgIdent;
                imp->sourceInt = compci->idx.fragment.frgSource;
                imp->position.bgn = bgn;
                imp->position.end = end;
                imp->delta_length = 0;
                imp->delta = NULL;
#if 0
                fprintf(stderr,
                        "Placing " F_CID " at " F_COORD "," F_COORD 
                        " based on positions " F_COORD "," F_COORD
                        " (compl %d length %d within input parent)\n",
                        imp->ident, bgn,end,
                        compci->position.bgn,compci->position.end,
                        cfrag->complement, cfrag->length);
#endif
                ci++;ifrag++;
              }
            }
          } else {

          int32 bgn,end;

          assert(isRead(cfrag->type));

          // make adjustments to positions due to application of traces??

          bgn = GetBead(beadStore,cfrag->beads)->column_index;
          end = GetBead(beadStore,cfrag->beads + cfrag->length -1 )->column_index + 1;
          if(cfrag->complement){
            int32 tmp = bgn;
            bgn = end;
            end = tmp;
          }

          imp = GetIntMultiPos(cma->f_list,ifrag);
          imp->position.bgn = bgn;
          imp->position.end = end;

#if 0
          fprintf(stderr,
                  "Placing " F_CID " at " F_COORD "," F_COORD 
                  " based on positions " F_COORD "," F_COORD
                  " (compl %d length %d within input parent)\n",
                  imp->ident, bgn,end,
                  offsets[i].bgn, offsets[i].end,
                  cfrag->complement, cfrag->length);
#endif
          ifrag++;
        }
      }
  }

  DeleteMANode(ma->lid);
  safe_free(offsets);
  return cma; 
}
/* end of MergeMultiAlign */

int32 AppendArtificialFragToLocalStore(FragType type, int32 iid, int complement,int32 contained,
                                       UnitigType utype, char *seq, char *qlt, int len) {
  static char seqbuffer[AS_READ_MAX_LEN+1];
  static char qltbuffer[AS_READ_MAX_LEN+1];
  char *sequence=seqbuffer,*quality=qltbuffer;
  int i;
  Fragment fragment;
  
  if ( len > AS_READ_MAX_LEN ) {
    fprintf(stderr, "AppendArtificialFragToLocalStore: input too long for buffer");
    assert(0);
  }
  for (i=0;i<len;i++) {
    seqbuffer[i]=*seq++; 
    qltbuffer[i]=*qlt++;
  }
  seqbuffer[len] = '\0';
  qltbuffer[len] = '\0';
  fragment.uid=iid;
  fragment.source = NULL;
  fragment.length = len;
  fragment.n_components = 0;  // no component frags or unitigs
  fragment.components = -1;
  if (complement) {
    SequenceComplement(sequence, quality);
  }
  fragment.lid = GetNumFragments(fragmentStore);
  fragment.iid = iid;
  fragment.type = type;
  fragment.utype = utype;
  fragment.complement = complement;
  fragment.contained = contained; 
  fragment.deleted = 0; 
  fragment.sequence = GetNumchars(sequenceStore);
  fragment.quality = GetNumchars(qualityStore);
  fragment.beads = GetNumBeads(beadStore);
  AppendRangechar(sequenceStore, fragment.length + 1, sequence);
  AppendRangechar(qualityStore, fragment.length + 1, quality);
  {Bead bead;
    int32 boffset;
    int32 soffset;
    int32 foffset;
    boffset = fragment.beads;
    soffset = fragment.sequence;
    bead.up = -1;
    bead.down = -1;
    bead.frag_index = fragment.lid;
    bead.column_index = -1;
    for (foffset = 0; foffset < fragment.length; foffset++ ) {
      bead.foffset = foffset;
      bead.boffset = boffset+foffset;
      bead.soffset = soffset+foffset;
      bead.next = bead.boffset +1;
      bead.prev = bead.boffset -1;
      if ( foffset == fragment.length-1) bead.next = -1;
      if ( foffset == 0 ) bead.prev = -1;
      SetVA_Bead(beadStore,boffset+foffset,&bead);
    }
  }
  AppendVA_Fragment(fragmentStore,&fragment);
  return fragment.lid;
}

int SetupSingleColumn(char *sequence, char *quality,
                      char *frag_type, char *unitig_type, CNS_Options *opp) 
{
  // returns the columnd id in the columnStore
  int32 fid,i;
  MANode *ma;
  int column_depth=0;
 
  if (sequence != NULL ) column_depth = strlen(sequence);
  if ( column_depth==0 ) return -1;
  RALPH_INIT = InitializeAlphTable();

  for (i=0;i<column_depth;i++) {
    fid = AppendArtificialFragToLocalStore((FragType)frag_type[i],
                                           i, 
                                           0,
                                           0,
                                           (UnitigType) unitig_type[i],
                                           &sequence[i],&quality[i],1);
  }


  ma = CreateMANode(GetNumMANodes(manodeStore));
  assert(ma->lid == 0);

  SeedMAWithFragment(ma->lid, GetFragment(fragmentStore,0)->lid,0, opp);
  for (i=1;i<column_depth;i++) {
    ApplyAlignment(i-1,0,i,0,NULL);
  }

  return GetMANode(manodeStore,ma->lid)->first;
}
