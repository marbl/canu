#include "positionDB.H"
#include "seqCache.H"

#if defined (__SVR4) && defined (__sun) 
// Solaris defines SS and GS in sys/regset.h
#undef GS
#undef SS
#endif

class tapperGlobalData {
public:
  tapperGlobalData();
  ~tapperGlobalData();

  void   initialize(void);

private:
  void     convertACGTtoColor(char *color, char *acgt, uint32 len);
  void     rewriteFileAsColorACGT(char *acgtname, char  *colorname);

public:
  char             *genName;
  char             *qryName;
  char             *outName;

  uint32            bgnRead;
  uint32            endRead;

  uint32            thisPartition;
  uint32            numPartitions;

  uint32            repeatThreshold;

  uint32            maxMemory;
  uint32            numThreads;
  bool              beVerbose;

  uint32            tagSize;

  uint32            maxColorError;
  uint32            maxBaseError;

  tapperTagFile        *TF;
  tapperResultFile     *TA;

  seqStream        *SS;
  merStream        *MS;
  positionDB       *PS;

  seqCache         *GS;
};





tapperGlobalData::tapperGlobalData() {
  genName   = 0L;
  qryName   = 0L;
  outName   = 0L;

  bgnRead =  uint32ZERO;
  endRead = ~uint32ZERO;

  thisPartition = 0;
  numPartitions = 1;

  repeatThreshold = 500;

  maxMemory  = 0;
  numThreads = 2;
  beVerbose  = false;

  maxColorError  = 3;
  maxBaseError   = 5;

  TF = 0L;
  TA = 0L;

  SS = 0L;
  MS = 0L;
  PS = 0L;

  GS = 0L;
}

tapperGlobalData::~tapperGlobalData() {
  delete TF;
  delete TA;
  delete PS;
  delete MS;
  delete SS;
  delete GS;
}



void
tapperGlobalData::initialize(void) {
  char  colName[FILENAME_MAX];

  sprintf(colName, "%s.colorspace", genName);
  rewriteFileAsColorACGT(genName, colName);

  TF = new tapperTagFile(qryName, 'r');

  if (numPartitions > 1) {
    if (thisPartition >= numPartitions) {
      fprintf(stderr, "ERROR: invalid partition n="uint32FMT" m="uint32FMT".\n", thisPartition, numPartitions);
      exit(1);
    }

    //  File has either fragment tags OR mate pairs, never both.

    uint32  numTags = (TF->numberOfFragmentTags() + TF->numberOfMatePairs()) / numPartitions + 1;

    bgnRead = numTags * thisPartition;
    endRead = numTags * thisPartition + numTags;

    fprintf(stderr, "Set partition for "uint64FMT" frags or "uint64FMT" mates: -begin "uint32FMT" -end "uint32FMT"\n",
            TF->numberOfFragmentTags(), TF->numberOfMatePairs(), bgnRead, endRead);
  }

  //  Set ranges that we want to compute.
  TF->setBegin(bgnRead);
  TF->setEnd(endRead);

  //  See the comments in the loader about the -1.
  tagSize = TF->metaData()->tagSize() - 1;

  if (tagSize > 32) {
    fprintf(stderr, "tag size too big for this implementation.\n");
    exit(1);
  }

  fprintf(stderr, "Building seqStream\n");
  SS = new seqStream(colName);

  fprintf(stderr, "Building merStream\n");
  MS = new merStream(new kMerBuilder(tagSize),
                     SS,
                     true, false);

  sprintf(colName, "%s.ms"uint32FMT".ce"uint32FMT".posDB", genName, tagSize, maxColorError);

  if (fileExists(colName)) {
    fprintf(stderr, "Loading positionDB\n");
    PS = new positionDB(colName, tagSize, 0, maxColorError);
  } else {
    fprintf(stderr, "Building positionDB\n");
    PS = new positionDB(MS, tagSize, 0, 0L, 0L, 0L, 0, 0, maxColorError, maxMemory, beVerbose);

    PS->saveState(colName);
  }

  delete MS;
  MS = 0L;

  GS = new seqCache(genName, 0, false);
  GS->loadAllSequences();

  TA = new tapperResultFile(outName, 'w');

  //  We get races unless we prebuild the AQI stuff.  I don't want to
  //  make this a requirement of the constructor, since only
  //  multithreaded codes have this problem, and it is perfectly valid
  //  for a file to have alignments with different max error
  //  rates....while tapper will only write with these two maximums.
  //
  TA->AQIlength(maxBaseError, maxColorError);
}



//  Inplace converts an acgt sequence to a color-space sequence.
void
tapperGlobalData::convertACGTtoColor(char *color, char *acgt, uint32 len) {
  char       l = 'n';  //  We always start the color encoding assuming the -1 letter is a gap
  char       n = 0;

  for (uint32 i=0; i<len; i++) {
    n        = acgt[i];
    color[i] = baseToColor[l][n];
    l        = n;
  }
}



void
tapperGlobalData::rewriteFileAsColorACGT(char *acgtname, char  *colorname) {
  seqCache   *F  = new seqCache(acgtname);

  if (fileExists(colorname)) {
    fprintf(stderr, "ColorFastA '%s' exists.  NOT recreating.\n", colorname);
    return;
  }

  fprintf(stderr, "Rewriting '%s' as ColorFastA '%s'.\n", acgtname, colorname);

  errno = 0;
  FILE       *CF = fopen(colorname, "w");
  if (errno)
    fprintf(stderr, "Failed to create '%s': %s\n", colorname, strerror(errno)), exit(1);

  seqInCore *f = F->getSequenceInCore();
  while (f) {
    convertACGTtoColor(f->sequence(), f->sequence(), f->sequenceLength());
    fprintf(CF, "%s\n%s\n", f->header(), f->sequence());
    delete f;
    f = F->getSequenceInCore();
  }

  fclose(CF);

  delete F;
}
