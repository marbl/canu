class tapperGlobalData {
public:
  char   *genName;
  char   *colName;
  char   *qryName;

  u32bit  numThreads;
  bool    beVerbose;

  u32bit  tagSize;
  u32bit  maxError;

  seqFile    *QF;

  seqStream  *SS;
  merStream  *MS;
  positionDB *PS;


  tapperGlobalData() {
    genName   = 0L;
    colName   = 0L;
    qryName   = 0L;

    numThreads = 2;
    beVerbose  = true;

    tagSize   = 0;
    maxError  = 3;

    QF = 0L;

    SS = 0L;
    MS = 0L;
    PS = 0L;
  };

  ~tapperGlobalData() {
    delete QF;
    delete PS;
    delete MS;
    delete SS;
  }

  void   initialize(void) {
    rewriteFileAsColorACGT(genName, colName);

    QF = openSeqFile(qryName);

    //  Read the first sequence to get the tagSize we should use.  We
    //  don't support changing tagSizes during a run, and will filter
    //  out any sequences that aren't this length.  See the comments
    //  in the loader as to why there is a -2 here.
    //
    {
      seqInCore *s = QF->getSequenceInCore();
      tagSize = s->sequenceLength() - 2;
      delete s;

      fprintf(stderr, "Discovered tag size of "u32bitFMT".\n", tagSize);

      QF->rewind();
    }

    fprintf(stderr, "Building seqStream\n");
    SS = new seqStream(colName, true);

    fprintf(stderr, "Building merStream\n");
    MS = new merStream(new kMerBuilder(tagSize), SS);

    fprintf(stderr, "Building positionDB\n");
    PS = new positionDB(MS, tagSize, 0, 0L, 0L, 0L, 0, 0, beVerbose, true);
  };

};
