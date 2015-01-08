






class unitigConsensus {
public:
  unitigConsensus(MultiAlignT *ma_, CNS_Options *opp_) {
    ma       = ma_;
    numfrags = GetNumIntMultiPoss(ma->f_list);
    fraglist = GetVA_IntMultiPos(ma->f_list, 0);
    fragback = NULL;
    opp      = opp_;
    trace    = NULL;
    manode   = NULL;
    utgpos   = NULL;
    cnspos   = NULL;
    tiid     = 0;
    piid     = -1;

    frankensteinLen = 0;
    frankensteinMax = 0;
    frankenstein    = NULL;
    frankensteinBof = NULL;
  };

  ~unitigConsensus() {
    DeleteVA_int32(trace);
    if (manode)
      DeleteMANode(manode->lid);

    safe_free(fragback);
    safe_free(utgpos);
    safe_free(cnspos);
    safe_free(frankenstein);
    safe_free(frankensteinBof);
  };

  int32  initialize(int32 *failed); 

  void   reportStartingWork(void);
  void   reportFailure(int32 *failed);
  void   reportSuccess(int32 *failed);

  int32  moreFragments(void)  { tiid++;  return (tiid < numfrags); };

  int32  computePositionFromParent(bool doContained);
  int32  computePositionFromLayout(void);
  int32  computePositionFromAlignment(void);

  void   rebuild(bool recomputeFullConsensus);

  bool   rejectAlignment(bool allowBhang, bool allowAhang, ALNoverlap *O);

  int32  alignFragment(void);
  void   applyAlignment(int32 frag_aiid=-1, int32 frag_ahang=0, int32 *frag_trace=NULL);

  void   generateConsensus(void);
  void   restoreUnitig(void);

private:
  MultiAlignT    *ma;
  int32           numfrags;
  IntMultiPos    *fraglist;
  IntMultiPos    *fragback;

  CNS_Options    *opp;

  VA_TYPE(int32) *trace;
  int32           traceBgn;

  MANode         *manode;
  SeqInterval    *utgpos;   //  Original unitigger location, DO NOT MODIFY
  SeqInterval    *cnspos;   //  Actual location in frankenstein.

  //int32           ovl;    //  Expected overlap in bases to the frankenstein
  //int32           ahang;  //  Expected hangs to the frankenstein
  //int32           bhang;

  int32           tiid;   //  This frag IID
  int32           piid;   //  Parent frag IID - if -1, not valid

  int32           frankensteinLen;
  int32           frankensteinMax;
  char           *frankenstein;
  beadIdx         *frankensteinBof;
};
