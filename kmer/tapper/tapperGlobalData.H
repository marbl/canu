

class tapperGlobalData {
public:
  tapperGlobalData();
  ~tapperGlobalData();

  void   initialize(void);

private:
  void  convertACGTtoColor(char *color, char *acgt, u32bit len);
  void  rewriteFileAsColorACGT(char *acgtname, char  *colorname);

public:
  char             *genName;
  char             *colName;
  char             *qryName;
  char             *outName;

  recordFile       *outIndex;
  recordFile       *outFragment;
  recordFile       *outMated;
  recordFile       *outSingleton;
  recordFile       *outTangled;

  u32bit            numThreads;
  bool              beVerbose;

  u32bit            tagSize;

  u32bit            maxColorError;
  u32bit            maxBaseError;

  tapperTagFile    *TF;

  seqStream        *SS;
  merStream        *MS;
  positionDB       *PS;

  FastACache       *GS;
};





tapperGlobalData::tapperGlobalData() {
  genName   = 0L;
  colName   = 0L;
  qryName   = 0L;
  outName   = 0L;

  outIndex     = 0L;
  outFragment  = 0L;
  outMated     = 0L;
  outSingleton = 0L;
  outTangled   = 0L;

  numThreads = 2;
  beVerbose  = false;

  maxColorError  = 3;
  maxBaseError   = 5;

  TF = 0L;

  SS = 0L;
  MS = 0L;
  PS = 0L;

  GS = 0L;
}

tapperGlobalData::~tapperGlobalData() {
  delete outIndex;
  delete outFragment;
  delete outMated;
  delete outSingleton;
  delete outTangled;

  delete TF;
  delete PS;
  delete MS;
  delete SS;
  delete GS;
}



void
tapperGlobalData::initialize(void) {
  rewriteFileAsColorACGT(genName, colName);

  TF = new tapperTagFile(qryName);

  //  See the comments in the loader about the -1.
  tagSize = TF->metaData()->tagSize() - 1;

  if (tagSize > 32) {
    fprintf(stderr, "tag size too big for this implementation.\n");
    exit(1);
  }

  fprintf(stderr, "Building seqStream\n");
  SS = new seqStream(colName, true);

  fprintf(stderr, "Building merStream\n");
  MS = new merStream(new kMerBuilder(tagSize), SS);  //  MEMORY LEAK kMerBuilder

  fprintf(stderr, "Building positionDB\n");
  PS = new positionDB(MS, tagSize, 0, 0L, 0L, 0L, 0, 0, beVerbose, true);

  PS->setUpMismatchMatcher(maxColorError);

  GS = new FastACache(genName, 0, true, false);

  {
    char filename[FILENAME_MAX];

    sprintf(filename, "%s.tapperMappedIndex", outName);
    outIndex = new recordFile(filename, 0, sizeof(tapperResult), true);

    sprintf(filename, "%s.tapperMappedFragment", outName);
    outFragment = new recordFile(filename, 0, sizeof(tapperResultFragment), true);

    sprintf(filename, "%s.tapperMappedMated", outName);
    outMated = new recordFile(filename, 0, sizeof(tapperResultMated), true);

    sprintf(filename, "%s.tapperMappedSingleton", outName);
    outSingleton = new recordFile(filename, 0, sizeof(tapperResultSingleton), true);

    sprintf(filename, "%s.tapperMappedTangled", outName);
    outTangled = new recordFile(filename, 0, sizeof(tapperResultTangled), true);
  }
}



//  Inplace converts an acgt sequence to a color-space sequence.
void
tapperGlobalData::convertACGTtoColor(char *color, char *acgt, u32bit len) {
  char       l = 'n';  //  We always start the color encoding assuming the -1 letter is a gap
  char       n = 0;

  for (u32bit i=0; i<len; i++) {
    n        = acgt[i];
    color[i] = baseToColor[l][n];
    l        = n;
  }
}



void
tapperGlobalData::rewriteFileAsColorACGT(char *acgtname, char  *colorname) {
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

  fclose(CF);

  delete F;
}
