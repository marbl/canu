
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
#ifndef MULTIALIGNMENT_CNS_INCLUDE
#define MULTIALIGNMENT_CNS_INCLUDE

#include "AS_global.h"
#include "MultiAlignStore_CNS.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "AS_ALN_forcns.h"
#include "dpc_CNS.h"
#include "AS_MSG_pmesg.h"
#include "AS_SDB_SequenceDB.h"
#define CNS_MIN_QV 0
#define CNS_MAX_QV 60
#define CNS_NALPHABET 6
#define CNS_NP 32 
#define COMPARE_ARGS char *aseq, char *bseq, int beg, int end, int opposite, \
                    double erate, double thresh, int minlen, \
				 CompareOptions what

#define MIN_SIZE_OF_MANODE 10000
#define BC_MAX(a,b)  (((a)>(b))?(a):(b))
#define BC_MIN(a,b)  (((a)<(b))?(a):(b))
#define FREE(x) if((x)!=NULL) {free((char *)(x));(x)=NULL;}
#define MIN_ALLOCATED_DEPTH 100

extern int DUMP_UNITIGS_IN_MULTIALIGNCONTIG;
extern int VERBOSE_MULTIALIGN_OUTPUT;

#define CNS_OPTIONS_SPLIT_ALLELES_DEFAULT  1
#define CNS_OPTIONS_SMOOTH_WIN_DEFAULT    10
#define CNS_OPTIONS_MAX_NUM_ALLELES        2

typedef struct {
   int split_alleles;
   int smooth_win;
   int max_num_alleles;
} CNS_Options;

typedef struct {
int      id;
char    *bases;      // gapped sequence
int     *qvs;        // quality values
double   ave_qv; 
int      allele_id;
int      uglen;      // ungapped length
} Read;

typedef struct {
int    id;
int    num_reads;
int   *read_ids;
int    weight;
int    uglen;      // ungapped length
} Allele;

typedef struct {
/*  This structure is used when recalling consensus bases
 *  to use only one of two alleles
 */
int32    beg;         // position of the left boundary
int32    end;         // position of the right boundary
int32    nr;          // number of reads in the region of variation
int32    max_nr;
int      nb;          // number of "current" bases
int32    na;          // total number of detected alleles
int32    nca;         // number of confirmed alleles
char    *curr_bases;  // dim = nr
char    *types;       // dim = nr
int32   *iids;        // iids of the reads
Read    *reads;
Allele  *alleles;
int32  **dist_matrix; // nr x nr matrix of cross-distances between reads
} VarRegion;

// -----------------------------------
// Jason introduced this new structure to address previous
// ambiguous use of IntMultiPos to hold both Fragments and Unitigs.

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
    int32		         frgSource;    
} CNS_FragmentContigElement;

typedef struct {
    UnitigType           utgType;
    IntUnitig_ID         utgIdent;
    int32                utgFirst; // index of this unitig's first fragment in fragment_positions
    int32                utgLast; // index of this unitig's last fragment in fragment_positions
} CNS_UnitigContigElement;

#define CNS_ELEMENT_IS_FRAGMENT 'F' 
#define CNS_ELEMENT_IS_UNITIG 'U' 
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

// Macro expand the Variable-Array-of-Type data structure 
// and its associated functions.
VA_DEF(CNS_AlignedContigElement)

// -----------------------------------

typedef struct {
int32 boffset; // Location in BeadStore
int32 soffset; // Location in sequence/qualityStores
int32 foffset; // Location in Fragment sequence
int32 prev;    
int32 next;
int32 up;
int32 down;  // navigation in multialignment (global offsets)
int32 frag_index; // Location of containing fragment in fragmentStore
int32 column_index; // Location of alignment column in columnStore
} Bead;

// Macro expand the Variable-Array-of-Type data structure 
// and its associated functions.
VA_DEF(Bead)

// -----------------------------------

typedef struct {
FragType type;
UnitigType utype;
uint32 iid;
uint64 uid;
int32 lid;  // index in sequence/quality/fragment store
int32 length;
int complement;
int contained;
int deleted;
int manode;
int32 sequence;  // global index of first sequence character
int32 quality;  // global index of first quality character
int32 beads;  // global index of first "bead"
int32 n_components; // number of component frags (in case of "unitig" Fragments)
int32 components;  // global index of first component frag 
int32 bactig;  // global index of bactig backbone of unitig pseudo frag 
char *source;  // consensus just carried this through - no mods
} Fragment;

// Macro expand the Variable-Array-of-Type data structure 
// and its associated functions.
VA_DEF(Fragment)

// -----------------------------------

typedef struct {
   int32 count[CNS_NALPHABET];
   int32 depth;
} BaseCount;

typedef struct {
int32 lid;  // index in columnStore
int32 call; // global offset in beadStore;
int32 next; 
int32 prev; // navigation in columnStore;
int32 ma_id;     // MANode membership;
int32 ma_index;  // index in MANode; // refreshed only periodically
BaseCount base_count;
} Column;

// Macro expand the Variable-Array-of-Type data structure 
// and its associated functions.
VA_DEF(Column)

// -----------------------------------

typedef struct {
//  This is the basic multialignment atom: 
//  A collection (possibly empty) of columns
//  Given by their offsets in the global columnStore
int32 lid;      // MANode id in the manodeStore
int32 iid;      // MANode's iid
int32 first;
int32 last;
VA_TYPE(int32) *columns;
} MANode; 

// Macro expand the Variable-Array-of-Type data structure 
// and its associated functions.
VA_DEF(MANode)

// -----------------------------------

MANode * CreateMANode(int32 iid);

static char ALPHABET[] = {'-','a','c','g','t','n'};

static int RALPH_INIT=0;
static char RALPHABET[CNS_NP] = {'-','A','C','G','T',
                               'a','c','g','t',   // -A, -C, -G, -T
                                   'M','R','W',   //     AC, AG, AT
                                       'S','Y',   //         CG, CT
                                           'K',   //             GT
                                   'm','r','w',   //    -AC,-AG,-AT
                                       's','y',   //        -CG,-CT
                                           'k',   //            -GT
                               'V','H','D','B',   //ACG,ACT,AGT,CGT
                               'v','h','d','b',   //-ACG,-ACT,-AGT,-CGT
                                           'X','x','N'};// ACGT, -ACGT, ??
static char RALPHABETC[CNS_NP] = {'-','T','G','C','A',
                               't','g','c','a',   // -A, -C, -G, -T
                                   'K','Y','W',   //     AC, AG, AT
                                       'S','R',   //         CG, CT
                                           'M',   //             GT
                                   'k','y','w',   //    -AC,-AG,-AT
                                       's','r',   //        -CG,-CT
                                           'm',   //            -GT
                               'B','D','H','V',   //ACG,ACT,AGT,CGT
                               'b','d','h','v',   //-ACG,-ACT,-AGT,-CGT
                                       'X','x','N'};// ACGT,-ACGT, ??

static float COMP_BIAS[CNS_NP] = {0.2,.2,.2,.2,.2,
                                     0.0, 0.0, 0.0, 0.0,  // here will be SNP rate
                                          0.0, 0.0, 0.0,  
                                               0.0, 0.0,  
                                                    0.0,  
                                         0.,   0.,   0.,   // initially only 2 haplotypes
                                               0.,   0.,   // so "tri- quad- morphisms" not possible
                                                     0.,  
                                   0.,   0.,   0.,   0.,  
                                   0.,   0.,   0.,   0.,  
                                               0.,   0.,   0.};
static double TAU_MISMATCH = (double)1./(5. - 1.); 
static uint32 AMASK[] = {   
    013607700741, // -
    015670707042, // a
    016733131104, // c
    017355252210, // g
    017566464420}; 

int InitializeAlphTable(void);
char BaseComplement(char c);
void SequenceComplement(char *sequence, char *quality);

typedef enum {
  CNS_QUIET      = (int)'Q', // quiet,  print nothing
  CNS_STATS_ONLY = (int)'S', // print only 1-line statistic summary
  CNS_ALIGNMENT  = (int)'A', // print the multialignment, sans CNS
  CNS_CONSENSUS  = (int)'C', // print the multialignment, with CNS
  CNS_DOTS       = (int)'D', // print the multialignment, dot format
  CNS_NODOTS     = (int)'N', // print the multialignment, "nodot" format
  CNS_EDIT_SCORE = (int)'E', // print the edit score column by column
  CNS_VIEW_UNITIG = (int)'U',  // show the unitigs in the contig alignment
  CNS_VERBOSE = (int)'V'  // verbose pre-post refinment output
} CNS_PrintKey;   // determine the format for PrintAlignment

typedef enum {
  CNS_SMOOTH = 1, // only eliminate pairwise construction artifacts
  CNS_POLYX = 2, // align poly-X regions
  CNS_INDEL = 4 // push apart mushed block indels
}  CNS_RefineLevel;

MultiAlignT * MANodeToMultiAlignT(MANode *ma);
int32 GetMANodeLength(int32 mid);

typedef struct {
   Column column;
   int32 bead;
} ColumnBeadIterator;

typedef struct {
   int32 cid;
} ColumnIterator;

typedef struct {
   Fragment fragment;
   int32 bead;
} FragmentBeadIterator;

typedef struct {
   int32 manode_id;
   int32 bead;
} ConsensusBeadIterator;

int GetMANodeConsensus(int32 mid, VA_TYPE(char) *, VA_TYPE(char) *);
int GetMANodePositions(int32 , int, IntMultiPos *, int, IntUnitigPos *, 
    VA_TYPE(int32) *);
void PrintAlignment(FILE *, int32, int32, int32, CNS_PrintKey );
    int32 MergeRefine(int32 , IntMultiVar **, int32 *, CNS_Options *, int);

typedef enum {
  LEFT_SHIFT  = (int) 'L', // Left Shifted
  RIGHT_SHIFT = (int) 'R', // Right Shifted
  UNSHIFTED   = (int) 'U', // Unshifted
  MIXED_SHIFT = (int) 'M'  // shifted in different directions
} ShiftStatus;

typedef struct {
  int32 start_column, end_column, rows, columns, window_width;
  ShiftStatus shift;
  char *beads;
  char *calls;
} Abacus;

typedef struct {
  int32 ident;
  int32 length;
  float32 coverage_stat;
  int32 left;
  int32 right;
  UnitigType type;
  // UnitigStatus_t status;
  //double prob;
} UnitigData;

VA_DEF(UnitigData)

typedef struct {
  int32 ident;
  int32 length;
  int32 num_contig_pairs;
  int32 contig_pairs;
} ScaffoldData;
VA_DEF(ScaffoldData)

// MergeMultiAlignsFast_new is the original CGW/CNS interface for contigging
 // and is now a wrapper around the more modern MergeMultiAligns which allows
 // "contained" relationships among the input contigs
MultiAlignT *MergeMultiAlignsFast_new( tSequenceDB *, FragStoreHandle, 
    VA_TYPE(IntElementPos) *, int, int, Overlap *(*)(COMPARE_ARGS),
    CNS_Options *opp);

// MergeMultiAligns 
MultiAlignT *MergeMultiAligns( tSequenceDB *, FragStoreHandle, 
      VA_TYPE(IntMultiPos) *, int, int, Overlap *(*)(COMPARE_ARGS), CNS_Options *opp);


MultiAlignT *ReplaceEndUnitigInContig( tSequenceDB *, FragStoreHandle ,
      uint32 , uint32 , int , Overlap *(*)(COMPARE_ARGS), CNS_Options *opp);

void ResetStores(int32 num_frags, int32 num_columns);
int SetupSingleColumn(char *, char *, char *, char *, CNS_Options *opp);
int BaseCall(int32 , int , double *, VarRegion  *, int, char *, int, int,
   CNS_Options *);
void ShowColumn(int32 cid);

int MultiAlignUnitig(IntUnitigMesg *, FragStoreHandle, VA_TYPE(char) *,
    VA_TYPE(char) *, VA_TYPE(int32) *, CNS_PrintKey , int, 
    Overlap *(*COMPARE_FUNC)(COMPARE_ARGS), CNS_Options *opp);

int ExamineMANode(FILE *, int32 , int32 , UnitigData *, int , CNS_Options *opp);
void ResetBaseCount(BaseCount *b);
int IncBaseCount(BaseCount *b,char c);
char GetConfMM(BaseCount *b,int mask);
int BaseToInt(char c);
int32 AppendFragToLocalStore(FragType, int32, int , int32 ,
    char *, UnitigType , MultiAlignStoreT *);
int GetAlignmentTrace(int32 , int32 , int32 , int32 *, int32 , 
                      VA_TYPE(int32) *, OverlapType *,
                      Overlap *(*COMPARE_FUNC)(COMPARE_ARGS),
                      int , int );
int MultiAlignContig_NoCompute(FILE *, int ,MultiAlignT *,
    tSequenceDB *, VA_TYPE(UnitigData) *, CNS_Options *opp);

#endif
