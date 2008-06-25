class tapperThreadData {
public:
  tapperThreadData(tapperGlobalData *g) {
    posn1fMax = 1024;
    posn1fLen = 0;
    posn1f    = new u64bit [posn1fMax];

    posn1rMax = 1024;
    posn1rLen = 0;
    posn1r    = new u64bit [posn1rMax];

    posn2fMax = 1024;
    posn2fLen = 0;
    posn2f    = new u64bit [posn2fMax];

    posn2rMax = 1024;
    posn2rLen = 0;
    posn2r    = new u64bit [posn2rMax];

    hits1Max  = 1048576;
    hits1Len  = 0;
    hits1     = new u64bit [hits1Max];

    hits2Max  = 1048576;
    hits2Len  = 0;
    hits2     = new u64bit [hits2Max];

    numMers = 0;

    startTime = getTime();
  };
  ~tapperThreadData() {
    delete [] posn1f;
    delete [] posn1r;
    delete [] posn2f;
    delete [] posn2r;
    delete [] hits1;
    delete [] hits2;
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

  u64bit  hits1Max;
  u64bit  hits1Len;
  u64bit *hits1;

  u64bit  hits2Max;
  u64bit  hits2Len;
  u64bit *hits2;

  char    themer[256];

  u64bit  numMers;

  double  startTime;
};

