


typedef enum {
  abAbacus_Smooth  = 1, // only eliminate pairwise construction artifacts
  abAbacus_Poly_X  = 2, // align poly-X regions
  abAbacus_Indel   = 4  // push apart mushed block indels
} abAbacusRefineLevel;



//  This is the basic multialignment atom: A collection (possibly empty) of columns given by their
//  offsets in the global columnStore.
//
//  Former MANode
//
class abMultiAlign {
public:
  abMultiAlign() {
    iid = 0;
  };
  ~abMultiAlign() {
  };

  abMultiAlignID    ident(void)       { return(lid); };

  abColID           firstColumn(void) { return(first); };
  abColID           lastColumn(void)  { return(last); };

  //  GetMANodeLength - also guarded against null nodes and returned -1
  uint32            length(void) {
    return(columnList.size());
  };

  vector<abColID>  &columns(void) {  //  For unitigConsensus, to rebuild
    return(columnList);
  }

  void              getConsensus(abAbacus *abacus, char *bases, char *quals, uint32 &len, uint32 max);
  void              getConsensus(abAbacus *abacus, char *bases, char *quals);
  void              getConsensus(abAbacus *abacus, tgTig *tig);
  void              getPositions(abAbacus *abacus, tgTig *tig);

  void              mergeRefine(abAbacus *abacus);
  int32             refineWindow(abAbacus *abacus, abColumn *start_column, abColID stab_bgn);
  int32             refine(abAbacus            *abacus,
                           abAbacusRefineLevel  level,
                           uint32               from = 0,
                           uint32               to   = UINT32_MAX);

  //  There are two multiAlign displays; this one, and one in tgTig.
  void              display(abAbacus  *abacus,
                            FILE      *F    = stdout,
                            uint32     from = 0,
                            uint32     to   = UINT32_MAX);

private:
  abMultiAlignID         lid;         //  MANode id in the manodeStore
  uint32                 iid;         //  tigStore tigID
  abColID                first;       //  first column in the multialignment
  abColID                last;
  vector<abColID>        columnList;  //  Used in AbacusRefine to get random access to specific column

  friend class abAbacus;
};



