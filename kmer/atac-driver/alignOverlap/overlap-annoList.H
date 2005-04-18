

//  List of the annotation.  Used for classifying each piece of the
//  annotation, e.g., U followed by 1 followed by U means that
//  somebody really did map something uniquely, where Y followed by 1
//  is probably just an extension.
//
//  This only works if assemblyA is the reference!
//
class annoList {
public:
  char    type;
  u32bit  iid1,  pos1,  len1;   //  The position on the reference axis
  u32bit  iid2a, pos2a, len2a;  //  The position on mapping 1
  u32bit  iid2b, pos2b, len2b;  //  The position on mapping 2

  void add(char type_,
           u32bit  iid1_,  u32bit pos1_,  u32bit len1_,
           u32bit  match1, match_t *m1,
           u32bit  match2, match_t *m2) {
    type  = type_;
    iid1  = iid1_;
    pos1  = pos1_;
    len1  = len1_;

    iid2a = match1;
    pos2a = 0;
    len2a = 0;
    if (m1) {
      pos2a = m1->pos2;
      len2a = m1->len2;
    }

    iid2b = match2;
    pos2b = 0;
    len2b = 0;
    if (m2) {
      pos2b = m2->pos2;
      len2b = m2->len2;
    }
  }
};
