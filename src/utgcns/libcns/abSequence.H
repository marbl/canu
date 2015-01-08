
//  Former Fragment, a sequence loaded into the multialign
//
class abSequence {
public:
  abSequence() {
  };
  ~abSequence() {
  };

  abBeadID        firstBead(void) {
    return(firstbead);
  };

  bool             isRead(void) {
    return(is_read);
  };


private:
  //FragType    type;
  //UnitigType  utype;
  uint32           iid;            //  external, aka gatekeeper, ID;
  abBaseID         lid;            //  index in sequence/quality/fragment store

  bool             is_read;

  uint32           length;
  bool             complement;
  uint32           container_iid;  //  if non-zero, the iid of our container
  bool             is_contained;   //  if non-zero, consensus detected this fragment is contained
  bool             deleted;
  abMultiAlignID   manode;
  abSeqID          sequence;       //  global index of first sequence character
  abBeadID         firstbead;      //  global index of first "bead"
  uint32           n_components;   //  number of component frags (in case of "unitig" Fragments)
  uint32           components;     //  global index of first component frag

  friend class abAbacus;
};

