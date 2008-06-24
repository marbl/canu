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

    if (maxError == 0) {
      fprintf(stderr, "tapper: initialize()- maxError is 0 -- use snapper!\n");
      exit(1);
    }

    PS->setUpMismatchMatcher(maxError);
  };

private:

  //  Inplace converts an acgt sequence to a color-space sequence.
  void  convertACGTtoColor(char *color, char *acgt, u32bit len) {
    char       l = 'n';  //  We always start the color encoding assuming the -1 letter is a gap
    char       n = 0;

    for (u32bit i=0; i<len; i++) {
      n        = acgt[i];
      color[i] = baseToColor[l][n];
      l        = n;
    }
  };

  void  rewriteFileAsColorACGT(char *acgtname, char  *colorname) {
    seqFile    *F  = openSeqFile(acgtname);

    if (fileExists(colorname)) {
      fprintf(stderr, "ColorFastA '%s' exists.  NOT recreating.\n", colorname);
      return;
    }

    fprintf(stderr, "Rewriting '%s' as ColorFastA '%s'.\n", acgtname, colorname);

    errno = 0;
    FILE       *CF = fopen(colorname, "w");
    if (errno)
      fprintf(stderr, "Failed to create '%s': %s\n", colorname, strerror(errno)), exit(1);

    while (F->eof() == false) {
      seqInCore *f = F->getSequenceInCore();
      convertACGTtoColor(f->sequence(), f->sequence(), f->sequenceLength());
      fprintf(CF, "%s\n%s\n", f->header(), f->sequence());
      delete f;
    }

    delete F;
    fclose(CF);
  };

};
