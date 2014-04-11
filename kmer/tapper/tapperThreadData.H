class tapperThreadData {
public:
  tapperThreadData(tapperGlobalData *g) {
    posn1fMax  = 256 * 1024;
    posn1fLen  = 0;
    posn1f     = new uint64 [posn1fMax];

    posn1rMax  = 256 * 1024;
    posn1rLen  = 0;
    posn1r     = new uint64 [posn1rMax];

    posn2fMax  = 256 * 1024;
    posn2fLen  = 0;
    posn2f     = new uint64 [posn2fMax];

    posn2rMax  = 256 * 1024;
    posn2rLen  = 0;
    posn2r     = new uint64 [posn2rMax];

    numHappiesMax = 256 * 1024;
    tag1happies   = new uint32 [numHappiesMax];
    tag1mate      = new uint32 [numHappiesMax];
    tag1tangled   = new uint32 [numHappiesMax];
    tag2happies   = new uint32 [numHappiesMax];
    tag2mate      = new uint32 [numHappiesMax];
    tag2tangled   = new uint32 [numHappiesMax];
    tangle        = 0L;
  };

  ~tapperThreadData() {
    delete [] posn1f;
    delete [] posn1r;
    delete [] posn2f;
    delete [] posn2r;
    delete [] tag1happies;
    delete [] tag1mate;
    delete [] tag1tangled;
    delete [] tag2happies;
    delete [] tag2mate;
    delete [] tag2tangled;
    delete [] tangle;
  };

public:
  uint64  posn1fMax;
  uint64  posn1fLen;
  uint64 *posn1f;

  uint64  posn1rMax;
  uint64  posn1rLen;
  uint64 *posn1r;

  uint64  posn2fMax;
  uint64  posn2fLen;
  uint64 *posn2f;

  uint64  posn2rMax;
  uint64  posn2rLen;
  uint64 *posn2r;

  uint32  numHappiesMax;

  uint32 *tag1happies;
  uint32 *tag1mate;
  uint32 *tag1tangled;

  uint32 *tag2happies;
  uint32 *tag2mate;
  uint32 *tag2tangled;

  intervalList *tangle;
};

