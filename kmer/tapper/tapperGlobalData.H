#include "positionDB.H"


class tapperGlobalData {
public:
  tapperGlobalData();
  ~tapperGlobalData();

  void   initialize(void);

  u32bit   getQualityHistogramIndex(u32bit numBaseMismatch, u32bit numColorMismatch, u32bit numColorError) {
    return(histogramIndices[numBaseMismatch * 256 + numColorMismatch * 16 + numColorError]);
  };

private:
  void     convertACGTtoColor(char *color, char *acgt, u32bit len);
  void     rewriteFileAsColorACGT(char *acgtname, char  *colorname);

public:
  char             *genName;
  char             *qryName;
  char             *outName;

  recordFile       *outIndex;
  recordFile       *outFragment;
  recordFile       *outMated;
  recordFile       *outSingleton;
  recordFile       *outTangled;
  recordFile       *outAlignQual;

  u32bit            bgnRead;
  u32bit            endRead;

  u32bit            repeatThreshold;

  u32bit            maxMemory;
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

  u32bit           *histogramIndices;
  u32bit            histogramIndicesLen;
};





tapperGlobalData::tapperGlobalData() {
  genName   = 0L;
  qryName   = 0L;
  outName   = 0L;

  outIndex     = 0L;
  outFragment  = 0L;
  outMated     = 0L;
  outSingleton = 0L;
  outTangled   = 0L;
  outAlignQual = 0L;

  bgnRead =  u32bitZERO;
  endRead = ~u32bitZERO;

  repeatThreshold = 500;

  maxMemory  = 0;
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
  delete outAlignQual;

  delete TF;
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
  SS = new seqStream(colName, true);

  fprintf(stderr, "Building merStream\n");
  MS = new merStream(new kMerBuilder(tagSize), SS);  //  MEMORY LEAK kMerBuilder

  sprintf(colName, "%s.posDB", genName);

  if (fileExists(colName)) {
    fprintf(stderr, "Loading positionDB\n");
    PS = new positionDB(colName);
  } else {
    fprintf(stderr, "Building positionDB\n");
    PS = new positionDB(MS, tagSize, 0, 0L, 0L, 0L, 0, 0, maxColorError, maxMemory, beVerbose);

    PS->saveState(colName);
  }

  GS = new FastACache(genName, 0, true, false);

  {
    char filename[FILENAME_MAX];

    sprintf(filename, "%s.tapperMappedIndex", outName);
    outIndex = new recordFile(filename, 0, sizeof(tapperResultIndex), 'w');

    sprintf(filename, "%s.tapperMappedFragment", outName);
    outFragment = new recordFile(filename, 0, sizeof(tapperResultFragment), 'w');

    sprintf(filename, "%s.tapperMappedMated", outName);
    outMated = new recordFile(filename, 0, sizeof(tapperResultMated), 'w');

    sprintf(filename, "%s.tapperMappedSingleton", outName);
    outSingleton = new recordFile(filename, 0, sizeof(tapperResultSingleton), 'w');

    sprintf(filename, "%s.tapperMappedTangled", outName);
    outTangled = new recordFile(filename, 0, sizeof(tapperResultTangled), 'w');

    sprintf(filename, "%s.tapperMappedAlignQual", outName);
    outAlignQual = new recordFile(filename, 0, sizeof(u32bit), 'w');
  }

  {
    //  min base mismatches for i color mismatches - the min is (I
    //  think always) the sum of the mins for the prime decomposition.
    //    9 - 3,3,3   -> min 6 base mismatches
    //    9 - 2,3,4   -> min 5 base mismatches
    //    9 - 2,2,5   -> min 5 base mismatches
    //    9 - 2,2,2,3 -> min 5 base mismatches
    //
    //  max base mismatches is, for the most part, used defined, but 0
    //  and 1 color mismatches are forced to 0 color mismatches.
    //
    //  finally, it is impossible to have just one color mismatch.

    //                cm 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    u32bit bmmin[16] = { 0, 0, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8 };

    u32bit maxc = maxColorError;  //  max color errors given to the mapper
    u32bit maxb = maxBaseError;   //  max base  errors given to the mapper

    u32bit index = 0;

    histogramIndices = new u32bit [4096];

    for (u32bit ii=0; ii < 4096; ii++)
      histogramIndices[ii] = ~u32bitZERO;

    //  A special case for 0.
    for (u32bit ce=0; ce <= maxc; ce++) {
      fprintf(stderr, "histogramIndices["u32bitFMT"/"u32bitFMT"/"u32bitFMT"] = "u32bitFMT"\n", 0, 0, ce, index);
      histogramIndices[0 * 256 + 0 * 16 + ce] = index++;
    }
    
    for (u32bit cm=2; cm <= maxc; cm++)
      for (u32bit ce=0; ce <= maxc - cm; ce++)
        for (u32bit bm=bmmin[cm]; bm <= maxb; bm++) {
          fprintf(stderr, "histogramIndices["u32bitFMT"/"u32bitFMT"/"u32bitFMT"] = "u32bitFMT"\n", bm, cm, ce, index);
          histogramIndices[bm * 256 + cm * 16 + ce] = index++;
        }

    histogramIndicesLen = index;
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
