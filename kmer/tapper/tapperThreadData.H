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
  };

  ~tapperThreadData() {
    delete [] posn1f;
    delete [] posn1r;
    delete [] posn2f;
    delete [] posn2r;
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
};

