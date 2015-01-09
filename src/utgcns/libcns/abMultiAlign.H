


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
  abMultiAlign();
  ~abMultiAlign();

  abMultiAlignID ident(void) { return(lid); };

  //  GetMANodeLength - also guarded against null nodes and returned -1
  uint32         getLength(void) {
    return(columnList.size());
  };

  uint32         getConsensus(abAbacus *abacus,
                              vector<char> &bases,
                              vector<char> &quals);

  void            mergeRefine(abAbacus *abacus);
  int32           refineWindow(abAbacus *abacus, abColumn *start_column, abColID stab_bgn);
  int32           refine(abAbacus            *abacus,
                         abAbacusRefineLevel  level,
                         uint32               from = 0,
                         uint32               to   = UINT32_MAX);

private:
  abMultiAlignID         lid;         //  MANode id in the manodeStore
  uint32                 iid;         //  tigStore id
  abColID                first;       //  first column in the multialignment
  abColID                last;
  vector<abColID>        columnList;  //  Used in AbacusRefine to get random access to specific column

  friend class abAbacus;
};



