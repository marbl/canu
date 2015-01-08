

//  This is the basic multialignment atom: A collection (possibly empty) of columns given by their
//  offsets in the global columnStore.
//
//  Former MANode
//
class abMultiAlign {
public:
  abMultiAlign();
  ~abMultiAlign();

  //  GetMANodeLength - also guarded against null nodes and returned -1
  uint32         getLength(void) {
    return(columnList.size());
  };

  uint32         getConsensus(abAbacus *ab,
                              vector<char> &bases,
                              vector<char> &quals);

  void            mergeRefine(abAbacus *ab);


private:
  abMultiAlignID         lid;         //  MANode id in the manodeStore
  uint32                 iid;         //  tigStore id
  abColID                first;       //  first column in the multialignment
  abColID                last;
  vector<abColID>        columnList;  //  Used in AbacusRefine to get random access to specific column

  friend class abAbacus;
};



