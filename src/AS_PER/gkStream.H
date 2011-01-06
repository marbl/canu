
////////////////////////////////////////
//
//  gkStream -- instead of random access to get a gkFrag, this will do
//  (slightly) optimized sequential access.
//
class gkStream {
public:
  gkStream(gkStore *gkp, AS_IID start, AS_IID end, uint32 flags);
  ~gkStream();

  void    reset(AS_IID start, AS_IID end);
  int     next(gkFragment *fr);

private:
  gkStore           *gkp;

  int                flags;

  AS_IID             bgnIID;
  AS_IID             curIID;
  AS_IID             endIID;

  StreamStruct      *fpk;
  StreamStruct      *qpk;

  StreamStruct      *fnm;
  StreamStruct      *snm;
  StreamStruct      *qnm;

  StreamStruct      *fsb;
  StreamStruct      *ssb;
  StreamStruct      *qsb;

  friend class gkStore;
};


