
#ifndef ABACUS_H
#define ABACUS_H

static const char *rcsid_ABACUS_H = "$Id:  $";

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"

//  Probably can't change these

#define CNS_MIN_QV 0
#define CNS_MAX_QV 60

#define CNS_NALPHABET 7
#define CNS_NP 32

//  Probably can change these

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



//  Tables and other static data - needs to be global so subclasses can access it

extern bool    DATAINITIALIZED;
extern double  EPROB[CNS_MAX_QV - CNS_MIN_QV + 1];  // prob of error for each quality value
extern double  PROB [CNS_MAX_QV - CNS_MIN_QV + 1];  // prob of correct call for each quality value (1-eprob)
extern int32   RINDEX[256];
extern char    ALPHABET[6];
extern char    RALPHABET[CNS_NP];
extern char    RALPHABETC[CNS_NP];
extern double  TAU_MISMATCH;
extern uint32  AMASK[5];




class abSequence;
class abBase;
class abBead;
class abColumn;
class abMultiAlign;
class abAbacus;
class abSeqBeadIterator;
class abColBeadIterator;

#include "abIDs.H"
#include "abSequence.H"
#include "abBase.H"
#include "abBead.H"
#include "abColumn.H"
//#include "abIterators.H"  After abAbacus
#include "abVariants.H"
#include "abMultiAlign.H"





class abAbacus {
public:
  abAbacus(gkStore *gkpStore);
  ~abAbacus();

public:

  //
  //  Accessors to lists
  //

  abSequence   *getSequence(abSeqID sid) {
    if (sid.isValid() == false)
      return(NULL);
    assert(sid.get() < _sequencesLen);
    return(_sequences + sid.get());
  };


  char          getBase(abBaseID bid) {
    if (bid.isValid() == false)
      return(0);
    assert(bid.get() < _basesLen);
    return(_bases[bid.get()]);
  };

  char          getBase(abBeadID bid) {
    return(getBase(getBead(bid)->baseIdx()));  //  Crashes if bid is not valid.
  };

  char         *getBases(abSequence *seq) {
    return(_bases + seq->firstBase().get());
  };


  char          getQual(abBaseID bid) {
    if (bid.isValid() == false)
      return(0);
    assert(bid.get() < _basesLen);
    return(_quals[bid.get()]);
  };

  char          getQual(abBeadID bid) {
    return(getQual(getBead(bid)->baseIdx()));  //  Crashes if bid is not valid.
  };

  char         *getQuals(abSequence *seq) {
    return(_quals + seq->firstBase().get());
  };


  abBead       *getBead(abBeadID bid) {
    if (bid.isValid() == false)
      return(NULL);
    assert(bid.get() < _beadsLen);
    return(_beads + bid.get());
  };


  abColumn     *getColumn(abColID cid) {
    if (cid.isValid() == false)
      return(NULL);
    assert(cid.get() < _columnsLen);
    return(_columns + cid.get());
  };

  abColumn     *getColumn(abBeadID bid) {
    if (bid.isValid() == false)
      return(NULL);

    return(getColumn(getBead(bid)->colIdx()));
  };


  abMultiAlign *getMultiAlign(abMultiAlignID mid) {
    if (mid.isValid() == false)
      return(NULL);
    assert(mid.get() < _multiAlignsLen);
    return(_multiAligns + mid.get());
  };

  //
  //  Setters to lists - lists with no setter shouldn't be set except during construction
  //

  void          setBase(abBaseID bid, char b) {
    assert(bid.isValid() == true);
    assert(bid.get() < _basesLen);
    _bases[bid.get()] = b;
  };

  void          setQual(abBaseID bid, char q) {
    assert(bid.isValid() == true);
    assert(bid.get() < _basesLen);
    _quals[bid.get()] = q;
  };

  //
  //  Constructors
  //

  abSeqID         addRead  (gkStore *gkpStore, uint32 readID, bool complemented);  //  Adds gkpStore read 'readID' to the abacus; former AppendFragToLocalStore
  abSeqID         addUnitig(gkStore *gkpStore, uint32 readID, bool complemented);  //  NOT SUPPORTED
  
  abBeadID        addBead(char base=0, char qual=0);

  abBeadID        appendGapBead(abBeadID bid);
  abBeadID        prependGapBead(abBeadID bid);

  //  Former ColumnAppend
  abColID         appendColumn(abColID cid, abBeadID bid);
  abColID         prependColumn(abColID cid, abBeadID bid);

  //  Adds a column, with a single bead
  //  Former CreateColumn()
  //
  abColID         addColumn(abMultiAlignID mid, abBeadID bid);

  void            alignBeadToColumn(abColID cid, abBeadID bid, char *label);
  abBeadID        unalignBeadFromColumn(abBeadID bid);
  abBeadID        unalignTrailingGapBeads(abBeadID bid);

  void            checkColumnBaseCounts();

  void            lateralExchangeBead(abBeadID lid, abBeadID rid);


  //  After all reads are added, create a new multialign by specifying the first sequence
  //  Former SeedMAWithFragment
  //
  abMultiAlignID  addMultiAlign(abSeqID sid);


  //  Options:
  //    split_alleles - CNS_OPTIONS_SPLIT_ALLELES_DEFAULT = 1
  //    smooth_win    - CNS_OPTIONS_MIN_ANCHOR_DEFAULT    = 11
  //    do_phasing    - CNS_OPTIONS_DO_PHASING_DEFAULT    = 1
private:
  void   getReadsForVAR(abVarRegion &region, abColID *cids);

public:
  void   refreshMultiAlign(abMultiAlignID  mid,
                           uint32          quality       = 0,       //  0 (contig consensus), -2, 1 (unitig consensus)
                           uint32          splitAlleles  = 1,
                           uint32          smoothWindow  = 11,
                           uint32          doPhasing     = 1,
                           uint32         *nvars         = NULL,
                           abVarRead     **v_list        = NULL,
                           uint32          make_v_list   = 0,        //  0, 1 (set to majority) or 2 (make a var record, for contigs)
                           bool            getScores     = false);   //  WAS: 0, 1 (unitig) or 2 (contig)

private:
  void            baseCallMajority(abColID cid);
  void            baseCallQuality(abVarRegion   &vreg,
                                  abColID        cid,
                                  double        &var,
                                  int32          target_allele,
                                  bool           getScores,
                                  int32          split_alleles,
                                  int32          smooth_win);
public:
  char            baseCall(abVarRegion &vreg,
                           abColID      cid,
                           bool         highQuality,
                           double      &var,
                           int32        target_allele,
                           bool         getScores,
                           int32        split_alleles,    //  these need to be set in the class
                           int32        smooth_win);      //

  void            applyAlignment(abSeqID  afid,
                                 int32    alen, abBeadID *aindex,
                                 abSeqID  bfid,
                                 int32    ahang,
                                 int32   *trace);



  //
  //  Major operations
  //
#if 0
int      AbacusRefine(MANode *ma,
                        int32 from,
                        int32 to,
                        CNS_RefineLevel level,
                        CNS_Options *opp);



  void     PrintAlignment(FILE *print,
                          int32 mid,
                          int32 from,
                          int32 to);

#endif



public:
  //  The iterators make calls to abAbacus, and to make them inline, they need to be defined after abAbacus.
  //  That means these constructors cannot be inlined.
  abSeqBeadIterator  *createSeqBeadIterator(abSeqID sid);
  abColBeadIterator  *createColBeadIterator(abColID cid);

  //
  //  Objects that we're trying to build multialigns for
  //

private:
  uint32            _sequencesLen;
  uint32            _sequencesMax;
  abSequence       *_sequences;  //  VA_TYPE(Fragment) *fragmentStore; - indices into bases/quals/beads
public:
  uint32             numberOfSequences(void) { return(_sequencesLen); };

private:
  uint32            _basesLen;
  uint32            _basesMax;
  //abBase           *_bases;      //  VA_TYPE(char) *sequenceStore, *qualityStore
  char             *_bases;      //  Used to be a tuple-type, but then we couldn't grab a pointer to the
  char             *_quals;      //  read for use in alignment - we would have had to copy it first

  //
  //  Work space for the multialigns
  //

  uint32            _beadsLen;
  uint32            _beadsMax;
  abBead           *_beads;      //  VA_TYPE(Bead) *beadStore;

  uint32            _columnsLen;
  uint32            _columnsMax;
  abColumn         *_columns;    //  VA_TYPE(Column)   *columnStore;   - basecall, profile, index in alignment, indexed pointers to prev/next cols

  //
  //  Multialigns - more than one, so we can handle splitting/joining without creating a new abacus
  //

  uint32            _multiAlignsLen;
  uint32            _multiAlignsMax;
  abMultiAlign     *_multiAligns;      //  VA_TYPE(MANode)   *manodeStore;   - first and last cols, profile, index into multialign, VA of component cols

  //
  //
  //

  vector<int32>       _readIdx;      //  VA_TYPE(int32) *fragment_indices; - convenience array
  vector<int32>       _abacusIdx;    //  VA_TYPE(int32) *abacus_indices;   - convenience array

  //  unused?
  //vector<alignedContig>  //  VA_TYPE(CNS_AlignedContigElement) *fragment_positions;

  //  Only used for remembering surrogates in contigs
  //HashTable_AS          *fragmentMap;

  //  Parameters
  uint32  smoothingWindow;
  uint32  splitAlleles;


  //  Statistics
private:
  uint32 NumColumnsInUnitigs;
  uint32 NumRunsOfGapsInUnitigReads;
  uint32 NumGapsInUnitigs;

  uint32 NumColumnsInContigs;
  uint32 NumRunsOfGapsInContigReads;
  uint32 NumGapsInContigs;

  uint32 NumAAMismatches;

  uint32 NumVARRecords;
  uint32 NumVARStringsWithFlankingGaps;

  uint32 NumUnitigRetrySuccess;


private:
  //  Options
  uint32 VERBOSE_MULTIALIGN_OUTPUT;
  uint32 FORCE_UNITIG_ABUT;

  uint32 MULTIALIGN_PRINT_WIDTH;
  uint32 MULTIALIGN_PRINT_SPACING;


  //  Data for this multialignment
private:
#if 0
  int32         start_column;
  int32         end_column;
  int32         rows;
  int32         columns;
  int32         window_width;
  ShiftStatus   shift;
  char         *beads;
  char         *calls;
#endif
};


#include "abIterators.H"

#endif
