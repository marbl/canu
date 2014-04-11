
class searcherState {
public:
  uint64         posnMax;
  uint64         posnLen;
  uint64        *posn;

#ifdef __APPLE__
  uint32         pad;
#endif

  double         encodeTime;
  double         maskTime;
  double         searchTime;
  double         filterTime;

  searcherState() {
    posnMax = 16384;
    posnLen = 0;
    posn    = new uint64 [ posnMax ];

    encodeTime = 0.0;
    maskTime   = 0.0;
    searchTime = 0.0;
    filterTime = 0.0;
  };

  ~searcherState() {
    delete [] posn;
  };
};
