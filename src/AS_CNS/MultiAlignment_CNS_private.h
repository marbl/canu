
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

#ifndef MULTIALIGNMENT_CNS_PRIVATE_H
#define MULTIALIGNMENT_CNS_PRIVATE_H

static const char *rcsid_MULTIALIGNMENT_CNS_PRIVATE_H = "$Id: MultiAlignment_CNS_private.h,v 1.21 2011-11-17 08:19:59 brianwalenz Exp $";

#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapStore.h"

//  These are used ONLY IN MultiAlignment_CNS.c.

#define ALT_QV_THRESH                      30
#define DONT_SHOW_OLAP                      0
#define SHOW_OLAP                           1
#define MIN_AVE_QV_FOR_VARIATION           21
#define MIN_SUM_QVS_FOR_VARIATION          60
#define QV_FOR_MULTI_GAP                   14

#define CNS_DP_THRESH                       1e-6
#define CNS_DP_MINLEN                      30

#define CNS_NEG_AHANG_CUTOFF               -5

#define INITIAL_NR                        100
#define MAX_WINDOW_FOR_ABACUS_REFINE      100
#define STABWIDTH                           6

#undef DEBUG_ABACUS
#undef DEBUG_ABACUS_ALIGN
#undef DEBUG_VAR_RECORDS  //  BROKEN!
#undef DEBUG_GET_ALIGNMENT_TRACE

#define MSTRING_SIZE                        3
#define MAX_SIZE_OF_ADJUSTED_REGION         5

//  GetAlignmentTrace alignment contexts.
//
//  Old ones were AS_CONSENSUS=0, AS_MERGE=1
//
typedef enum {
  GETALIGNTRACE_UNITIG    = 'u',  //  unitig consensus: Fragment to fragment
  GETALIGNTRACE_CONTIGU   = 'U',  //  contig consensus: Unitig to unitig
  GETALIGNTRACE_CONTIGF   = 'F',  //  contig consensus: Fragment to unitig
  GETALIGNTRACE_MERGE     = 'm'   //  cgw scaffold merging
} GetAlignmentTraceContext;

#define CNS_MIN_QV 0
#define CNS_MAX_QV 60

#define RINDEXMAX 128

#define CNS_NALPHABET 7
#define CNS_NP 32


typedef struct {
  int32    id;
  int32    iid;
  char    *bases;      // gapped sequence
  int32   *qvs;        // quality values
  double   ave_qv;
  int32    allele_id;
  int32    uglen;      // ungapped length
} Read;

typedef struct {
  int32  id;
  int32  num_reads;
  int32 *read_ids;
  int32 *read_iids;
  int32  weight;
  int32  uglen;      // ungapped length
} Allele;

typedef struct {
  /*  This structure is used when recalling consensus bases
   *  to use only one of two alleles
   */
  int32    beg;         // position of the left boundary
  int32    end;         // position of the right boundary
  int32    nr;          // number of reads in the region of variation
  int32    max_nr;
  int32    nb;          // number of "current" bases
  int32    na;          // total number of detected alleles
  int32    nca;         // number of confirmed alleles
  char    *curr_bases;  // dim = nr
  char    *types;       // dim = nr
  int32   *iids;        // iids of the reads
  Read    *reads;
  Allele  *alleles;
  int32  **dist_matrix; // nr x nr matrix of cross-distances between reads
} VarRegion;

// CNS_AlignedContigElement holds temporary summaries of
// components of a contig.
// Each component may be a fragment or a unitig of fragments.
// The flag frg_or_utg distinguishes these cases.
// This structure is local to MultiAlignment_CNS.

typedef struct {
  FragType             frgType;
  IntFragment_ID       frgIdent;
  IntFragment_ID       frgContained;
  IntUnitig_ID         frgInUnitig;
} CNS_FragmentContigElement;

typedef struct {
  UnitigType           utgType;
  IntUnitig_ID         utgIdent;
} CNS_UnitigContigElement;


#define CNS_ELEMENT_IS_FRAGMENT 'F'
#define CNS_ELEMENT_IS_UNITIG   'U'

typedef struct {
  union {
    CNS_FragmentContigElement fragment;
    CNS_UnitigContigElement unitig;
  } idx;
  char                 frg_or_utg; // use CNS_ELEMENT_IS_FRAGMENT
  SeqInterval          position;
  int32                delta_length;
  int32               *delta;
} CNS_AlignedContigElement;

VA_DEF(CNS_AlignedContigElement)




//  This class is designed to prevent spurious assignment to anything except other beadIdx.  A
//  method is provided to explicitly set the value, but doing "beadIdx b = 4" will not work.
//
//  operator size_t -- used for automatically converting to an int for use in GetBead().
//  set()           -- force a setting
//  isInvalid()     -- test if the bead has a valid index
//
//  Ideally, 'operator size_t' would go away, since it could be hiding some bad usage,
//  like 'int32 bid = someOtherBeadIdx', but there are 192 instances of GetBead() and I
//  don't currently feel like tracking those down -- it'll be easier to replace the
//  beadStore with a real class.
//
class beadIdx {
public:
  beadIdx()  { _idx = 0xffffffff; };
  ~beadIdx() {                            };

  //operator const size_t () { return(_idx); };

  void   set(size_t forcedIdx) { _idx = forcedIdx; };
  uint32 get(void)             { return(_idx); };

  bool   isInvalid(void) { return(_idx == 0xffffffff); };
  bool   isValid(void)   { return(_idx != 0xffffffff); };

  bool   operator!=(beadIdx that) const  { return(_idx != that._idx); };
  bool   operator==(beadIdx that) const  { return(_idx == that._idx); };
private:
  uint32   _idx;
};

class seqIdx {
public:
  seqIdx()  { _idx = 0xffffffff; };
  ~seqIdx() {                            };

  operator const size_t () { return(_idx); };

  void   set(size_t forcedIdx) { _idx = forcedIdx; };
  uint32 get(void)             { return(_idx); };

  bool   isInvalid(void) { return(_idx == 0xffffffff); };
  bool   isValid(void)   { return(_idx != 0xffffffff); };

  bool   operator!=(seqIdx that) const  { return(_idx != that._idx); };
  bool   operator==(seqIdx that) const  { return(_idx == that._idx); };
private:
  uint32   _idx;
};





//VA_DEF(beadIdx)

typedef struct {
  beadIdx boffset; // Location in BeadStore
  seqIdx  soffset; // Location in sequence/qualityStores
  int32  foffset; // Location in Fragment sequence
  beadIdx prev;
  beadIdx next;
  beadIdx up;
  beadIdx down;  // navigation in multialignment (global offsets)
  int32  frag_index; // Location of containing fragment in fragmentStore
  int32  column_index; // Location of alignment column in columnStore
} Bead;

VA_DEF(Bead)

inline
Bead*
GetBead(VA_TYPE(Bead) *beadStore, beadIdx bid) {
  return(GetBead(beadStore, bid.get()));
}


typedef struct {
  FragType type;
  UnitigType utype;
  uint32 iid;
#ifdef PRINTUIDS
  uint64 uid;
#endif
  int32   lid;            // index in sequence/quality/fragment store
  int32   length;
  int32   complement;
  int32   container_iid;    // if non-zero, the iid of our container
  int32   is_contained;     // if non-zero, consensus detected this fragment is contained
  int32   deleted;
  int32   manode;
  seqIdx   sequence;       // global index of first sequence character
  beadIdx firstbead;      // global index of first "bead"
  int32   n_components;   // number of component frags (in case of "unitig" Fragments)
  int32   components;     // global index of first component frag
  char   *source;         // consensus just carried this through - no mods
} Fragment;

VA_DEF(Fragment)


typedef struct {
  int32 count[CNS_NALPHABET];
  int32 depth;
} BaseCount;

typedef struct {
  int32     lid;  // index in columnStore
  beadIdx   call; // global offset in beadStore;
  int32     next;
  int32     prev; // navigation in columnStore;
  int32     ma_id;     // MANode membership;
  int32     ma_index;  // index in MANode; // refreshed only periodically // seems to also be gapped position in the align
  BaseCount base_count;
} Column;

VA_DEF(Column)


//  This is the basic multialignment atom: A collection (possibly
//  empty) of columns given by their offsets in the global columnStore
typedef struct {
  int32  lid;      // MANode id in the manodeStore
  int32  iid;      // MANode's iid
  int32  first;
  int32  last;
  VA_TYPE(int32) *columnList;  //  Used in AbacusRefine to get random access to specific column
} MANode;

VA_DEF(MANode)



typedef struct {
  Column column;
  beadIdx bead;
} ColumnBeadIterator;

typedef struct {
  Fragment fragment;
  beadIdx   bead;
  bool     isNull;
} FragmentBeadIterator;

typedef struct {
  int32 manode_id;
  beadIdx bead;
} ConsensusBeadIterator;


typedef enum {
  LEFT_SHIFT  = (int) 'L', // Left Shifted
  RIGHT_SHIFT = (int) 'R', // Right Shifted
  UNSHIFTED   = (int) 'U', // Unshifted
  MIXED_SHIFT = (int) 'M'  // shifted in different directions
} ShiftStatus;

typedef struct {
  int32 start_column;
  int32 end_column;
  int32 rows;
  int32 columns;
  int32 window_width;
  ShiftStatus shift;
  char *beads;
  char *calls;
} AbacusDataStructure;



extern gkStore               *gkpStore;
extern OverlapStore          *ovlStore;
extern MultiAlignStore       *tigStore;

extern HashTable_AS          *fragmentMap;

extern VA_TYPE(char) *sequenceStore;
extern VA_TYPE(char) *qualityStore;
extern VA_TYPE(Bead) *beadStore;

extern VA_TYPE(Fragment) *fragmentStore;
extern VA_TYPE(Column)   *columnStore;
extern VA_TYPE(MANode)   *manodeStore;

extern VA_TYPE(int32) *fragment_indices;
extern VA_TYPE(int32) *abacus_indices;

extern VA_TYPE(CNS_AlignedContigElement) *fragment_positions;

extern double EPROB[CNS_MAX_QV-CNS_MIN_QV+1];
extern double PROB[CNS_MAX_QV-CNS_MIN_QV+1];
extern int32  RINDEX[RINDEXMAX];
extern char   ALPHABET[6];
extern char   RALPHABET[CNS_NP];
extern char   RALPHABETC[CNS_NP];
extern double TAU_MISMATCH;
extern uint32 AMASK[5];

extern int32 thisIsConsensus;

extern int32 NumColumnsInUnitigs;
extern int32 NumRunsOfGapsInUnitigReads;
extern int32 NumGapsInUnitigs;
extern int32 NumColumnsInContigs;
extern int32 NumRunsOfGapsInContigReads;
extern int32 NumGapsInContigs;
extern int32 NumAAMismatches;
extern int32 NumVARRecords;
extern int32 NumVARStringsWithFlankingGaps;
extern int32 NumUnitigRetrySuccess;

extern int32 DUMP_UNITIGS_IN_MULTIALIGNCONTIG;
extern int32 VERBOSE_MULTIALIGN_OUTPUT;
extern int32 FORCE_UNITIG_ABUT;

#define MIN_STAIRCASE_TRACE_SIZE 2

extern ssize_t previousStaircaseTraceEntry;
extern size_t previousStaircaseSize;
extern int32 previousStaircaseFragmentId;

//  Functions used by lots of pieces internally to AS_CNS.  Defined in
//  MultiAlgnment_CNS.c.

int
IncBaseCount(BaseCount *b, char c);
int
DecBaseCount(BaseCount *b, char c);
int
GetBaseCount(BaseCount *b, char c);
int
GetColumnBaseCount(Column *b, char c);
int
GetDepth(Column *c);
void
ResetBaseCount(BaseCount *b);
void
ShowBaseCount(BaseCount *b);
void
ShowBaseCountPlain(FILE *out,BaseCount *b);
char
GetMaxBaseCount(BaseCount *b, int32 start_index);
void
CheckColumnBaseCount(Column *c);

MANode *
CreateMANode(int32 iid);
void
DeleteMANode(int32 iid);
int32
GetMANodeLength(int32 mid);
void
SeedMAWithFragment(int32 mid,
                   int32 fid,
                   int32 quality,
                   CNS_Options *opp);
int
GetMANodeConsensus(int32 mid, VA_TYPE(char) *sequence, VA_TYPE(char) *quality);
int
GetMANodePositions(int32 mid, MultiAlignT *ma);

void
CreateColumnBeadIterator(int32 cid,ColumnBeadIterator *bi);
beadIdx
NextColumnBead(ColumnBeadIterator *bi);
void
NullifyFragmentBeadIterator(FragmentBeadIterator *bi);
int
IsNULLIterator(FragmentBeadIterator *bi);
void
CreateFragmentBeadIterator(int32 fid,FragmentBeadIterator *bi);
beadIdx
NextFragmentBead(FragmentBeadIterator *bi);
void
CreateConsensusBeadIterator(int32 mid,ConsensusBeadIterator *bi);
beadIdx
NextConsensusBead(ConsensusBeadIterator *bi);

void
ClearBead(beadIdx bid);
void
AlignBeadToColumn(int32 cid, beadIdx bid, char *label);
beadIdx
UnAlignBeadFromColumn(beadIdx bid);
beadIdx
UnAlignTrailingGapBeads(beadIdx bid);
void
LateralExchangeBead(beadIdx lid, beadIdx rid);
beadIdx
AppendGapBead(beadIdx bid);
beadIdx
PrependGapBead(beadIdx bid);

Column *
CreateColumn(beadIdx bid);
void
AddColumnToMANode(int32 ma, Column column);
int32
ColumnAppend(int32 cid, beadIdx bid);
void
ShowColumn(int32 cid);

void
ResetStores(int32 num_bases, int32 num_frags, int32 num_columns);
int32
AppendFragToLocalStore(FragType          type,
                       int32             iid,
                       int32             complement,
                       int32             contained,
                       UnitigType        utype);

void
AllocateDistMatrix(VarRegion  *vreg, int32 init);
void
OutputDistMatrix(FILE *fout, VarRegion  *vreg);
void
PopulateDistMatrix(Read *reads, int32 len, VarRegion  *vreg);
void
OutputReads(FILE *fout, Read *reads, int32 nr, int32 width);
void
OutputAlleles(FILE *fout, VarRegion *vreg);
void
AllocateMemoryForReads(Read **reads, int32 nr, int32 len,
                       int32 default_qv);
void
AllocateMemoryForAlleles(Allele **alleles, int32 nr, int32 *na);
void
SortAllelesByLength(Allele *alleles, int32 num_alleles, Read *reads);
void
SortAllelesByWeight(Allele *alleles, int32 num_alleles, Read *reads);
void
SortAllelesByMapping(Allele *alleles, int32 nca, Read *reads, int32 *allele_map);
void
ClusterReads(Read *reads, int32 nr, Allele *alleles, int32 *na, int32 *nca, int32 **dist_matrix);


//
//  Main blocks of functionality.  All are in files named after the function. 
//
int
AbacusRefine(MANode *ma, int32 from, int32 to, CNS_RefineLevel level,
             CNS_Options *opp);

int
RefreshMANode(int32 mid, int32 quality, CNS_Options *opp, int32 *nvars,
              IntMultiVar **v_list, int32 make_v_list, int32 get_scores);

void
ApplyAlignment(int32 afid,
               int32 alen, beadIdx *aindex,
               int32 bfid,
               int32 ahang,
               int32 *trace);

bool IsStaircaseAlignment(int32 fragmentId, int32* trace);

void
PrintAlignment(FILE *print, int32 mid, int32 from, int32 to, CNS_PrintKey what);

void
MergeRefine(int32 mid, VA_TYPE(IntMultiVar) *v_list,
            int32 utg_alleles, CNS_Options *opp, int32 get_scores);


int
GetAlignmentTrace(int32                      afid,
                  char                      *aseq_input,
                  int32                      bfid,
                  int32                     *ahang,
                  int32                     *bhang,
                  int32                      expected_length,
                  VA_TYPE(int32)            *trace,
                  OverlapType               *otype,
                  AS_ALN_Aligner            *alignFunction,
                  int32                      show_olap,
                  int32                      allow_big_endgaps,
                  GetAlignmentTraceContext   alignment_context,
                  double                     input_erate);

int
GetAlignmentTraceDriver(Fragment                    *afrag,
                        char                        *aseq,
                        Fragment                    *bfrag,
                        int32                       *ahang,
                        int32                       *bhang,
                        int32                        expected_length,
                        VA_TYPE(int32)              *trace,
                        OverlapType                 *otype,
                        GetAlignmentTraceContext     alignment_context,
                        int32                        max_gap);

int
BaseCall(int32 cid, int32 quality, double *var, VarRegion  *vreg,
         int32 target_allele, char *cons_base, int32 verbose, int32 get_scores,
         CNS_Options *opp);

#endif
