class tapperThreadData {
public:
  tapperThreadData(tapperGlobalData *g) {
    posn1fMax  = 256 * 1024;
    posn1fLen  = 0;
    posn1f     = new u64bit [posn1fMax];

    posn1rMax  = 256 * 1024;
    posn1rLen  = 0;
    posn1r     = new u64bit [posn1rMax];

    posn2fMax  = 256 * 1024;
    posn2fLen  = 0;
    posn2f     = new u64bit [posn2fMax];

    posn2rMax  = 256 * 1024;
    posn2rLen  = 0;
    posn2r     = new u64bit [posn2rMax];

    numHappiesMax = 256 * 1024;
    tag1happies   = new u32bit [numHappiesMax];
    tag1mate      = new u32bit [numHappiesMax];
    tag2happies   = new u32bit [numHappiesMax];
    tag2mate      = new u32bit [numHappiesMax];

    tangle        = 0L;
  };

  ~tapperThreadData() {
    delete [] posn1f;
    delete [] posn1r;
    delete [] posn2f;
    delete [] posn2r;
    delete [] tag1happies;
    delete [] tag1mate;
    delete [] tag2happies;
    delete [] tag2mate;
    delete [] tangle;
  };

public:
  u64bit  posn1fMax;
  u64bit  posn1fLen;
  u64bit *posn1f;

  u64bit  posn1rMax;
  u64bit  posn1rLen;
  u64bit *posn1r;

  u64bit  posn2fMax;
  u64bit  posn2fLen;
  u64bit *posn2f;

  u64bit  posn2rMax;
  u64bit  posn2rLen;
  u64bit *posn2r;

  u32bit  numHappiesMax;
  u32bit *tag1happies;
  u32bit *tag1mate;
  u32bit *tag2happies;
  u32bit *tag2mate;

  intervalList *tangle;
};

