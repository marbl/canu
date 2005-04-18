//  Given two ATAC-format match files, compute the amount they agree / disagree.

class match_t {
public:
  u32bit  matchiid;  //  internal id, usually pointing to an entry in matchList
  u32bit  matchuid;  //  external id
  u32bit  iid1, pos1, len1, ori1;
  u32bit  iid2, pos2, len2, ori2;
};


