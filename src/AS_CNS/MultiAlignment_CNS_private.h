
//  These are used ONLY IN MultiAlignment_CNS.c.

typedef struct {
  int      id;
  int32    iid;
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
  int   *read_iids;
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

VA_DEF(CNS_AlignedContigElement)


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

VA_DEF(Bead)


typedef struct {
  FragType type;
  UnitigType utype;
  uint32 iid;
#ifdef PRINTUIDS
  uint64 uid;
#endif
  int32 lid;            // index in sequence/quality/fragment store
  int32 length;
  int complement;
  int container_iid;    // if non-zero, the iid of our container
  int is_contained;     // if non-zero, consensus detected this fragment is contained
  int deleted;
  int manode;
  int32 sequence;       // global index of first sequence character
  int32 quality;        // global index of first quality character
  int32 firstbead;      // global index of first "bead"
  int32 n_components;   // number of component frags (in case of "unitig" Fragments)
  int32 components;     // global index of first component frag
  char *source;         // consensus just carried this through - no mods
} Fragment;

VA_DEF(Fragment)


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

VA_DEF(Column)


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

VA_DEF(MANode)


static char ALPHABET[] = {'-','a','c','g','t','n'};

static int RALPH_INIT=0;
static char RALPHABET[CNS_NP] = {'-','A','C','G','T','N',
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
static char RALPHABETC[CNS_NP] = {'-','T','G','C','A','N',
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

static double TAU_MISMATCH = (double)1./(5. - 1.);
static uint32 AMASK[] = {
  013607700741, // -
  015670707042, // a
  016733131104, // c
  017355252210, // g
  017566464420};



typedef struct {
  Column column;
  int32 bead;
} ColumnBeadIterator;

typedef struct {
  Fragment fragment;
  int32 bead;
} FragmentBeadIterator;

typedef struct {
  int32 manode_id;
  int32 bead;
} ConsensusBeadIterator;


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
  float   coverage_stat;
  int32 left;
  int32 right;
  UnitigType type;
} UnitigData;

VA_DEF(UnitigData)

typedef struct {
  int32 ident;
  int32 length;
  int32 num_contig_pairs;
  int32 contig_pairs;
} ScaffoldData;

VA_DEF(ScaffoldData)



extern GateKeeperStore       *gkpStore;
extern tSequenceDB           *sequenceDB;
extern HashTable_AS          *fragmentMap;
extern MultiAlignStoreT      *unitigStore;

extern VA_TYPE(char) *sequenceStore;
extern VA_TYPE(char) *qualityStore;
extern VA_TYPE(Bead) *beadStore;

extern VA_TYPE(Fragment) *fragmentStore;
extern VA_TYPE(Column)   *columnStore;
extern VA_TYPE(MANode)   *manodeStore;

extern int USE_SDB;

//extern int allow_forced_frags;
//extern int allow_neg_hang;

extern int NumColumnsInUnitigs;
extern int NumRunsOfGapsInUnitigReads;
extern int NumGapsInUnitigs;
extern int NumColumnsInContigs;
extern int NumRunsOfGapsInContigReads;
extern int NumGapsInContigs;
extern int NumAAMismatches;
extern int NumVARRecords;
extern int NumVARStringsWithFlankingGaps;
extern int NumUnitigRetrySuccess;

extern int DUMP_UNITIGS_IN_MULTIALIGNCONTIG;
extern int VERBOSE_MULTIALIGN_OUTPUT;
extern int FORCE_UNITIG_ABUT;
extern int clear_range_to_use;
