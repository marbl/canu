
class searcherState {
public:
  u64bit         posnMax;
  u64bit         posnLen;
  u64bit        *posn;

#ifdef __APPLE__
  u32bit         pad;
#endif

  double         encodeTime;
  double         maskTime;
  double         searchTime;
  double         filterTime;

  searcherState() {
    posnMax = 16384;
    posnLen = 0;
    posn    = new u64bit [ posnMax ];

    encodeTime = 0.0;
    maskTime   = 0.0;
    searchTime = 0.0;
    filterTime = 0.0;
  };

  ~searcherState() {
    delete [] posn;
  };
};
