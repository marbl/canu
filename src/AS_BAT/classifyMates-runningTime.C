



class onlineMeanStdDev {
public:
  onlineMeanStdDev() {
    hMin   = UINT32_MAX;
    hMax   = 0;
    hAlc   = 0;
    hDat   = new uint32 [hAlc];

    memset(hDat, 0, sizeof(uint32) * hAlc);

    hCnt   = 0;

    computedMean   = 0.0;
    computedStdDev = 0.0;
  };
  ~onlineMeanStdDev() {
    delete [] hDat;
  };

  uint32    min(void)     { return(hMin);           };
  uint32    max(void)     { return(hMax);           };

  double    mean(void)    { return(computedMean);   };
  double    stddev(void)  { return(computedStdDev); };

  uint32    numData(void) { return(hCnt);           };

  void      addDataPoint(uint32 dp) {

    if (dp >= hAlc) {
      while (hAlc <= dp)
        hAlc++;

      uint32  *hNew = new uint32 [hAlc];

      memset(hNew, 0,    sizeof(uint32) * (hAlc));
      memcpy(hNew, hDat, sizeof(uint32) * (hMax + 1));

      delete [] hDat;
      hDat = hNew;
    }

    if (hMin > dp)  hMin = dp;
    if (hMax < dp)  hMax = dp;

    hDat[dp]++;
    
    hCnt++;
  };

  void      recompute(void) {
    double m = 0;  //  Scratch copies; we update the real ones quickly at the end.
    double s = 0;

    uint32 c = 0;

    for (uint32 dp=hMin; dp<=hMax; dp++) {
      c += hDat[dp];
      m += hDat[dp] * dp;
    }

    if (c != hCnt)
      fprintf(stderr, "WARNING: c=%u hCnt=%u\n", c, hCnt);
    assert(c == hCnt);

    m /= hCnt;

    for (uint32 dp=hMin; dp<=hMax; dp++)
      s += hDat[dp] * (dp - m) * (dp - m);

    s /= hCnt - 1;
    s  = sqrt(s);

    computedMean   = m;
    computedStdDev = s;
  };

private:
  uint32    hMin;  //  Minimum value we've seen.
  uint32    hMax;  //  Maximum value we've seen.
  uint32    hAlc;  //  Maximum value we've allocated space for in hDat.
  uint32   *hDat;  //  Histogram of values.

  uint32    hCnt;  //  Number of datapoints - sum of hDat.

  double    computedMean;
  double    computedStdDev;
};
