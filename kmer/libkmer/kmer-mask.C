#include <stdio.h>
#include <stdlib.h>

#include "util++.H"
#include "bio++.H"
#include "libmeryl.H"
#include "existDB.H"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

#include "sweatShop.H"




class fastqRecord {
public:
  fastqRecord(uint32 ml) {
    maxLength = ml;
    alloc     = new char [maxLength * 8];

    a2        = alloc + 0 * maxLength;
    af        = alloc + 1 * maxLength;
    am        = alloc + 2 * maxLength;
    a4        = alloc + 3 * maxLength;

    a1[0]     = 0;
    a2[0]     = 0;
    af[0]     = 0;
    am[0]     = 0;
    a3[0]     = 0;
    a4[0]     = 0;

    aLength   = 0;
    aRetained = 0.0;
    aLabel    = 0;

    b2        = alloc + 4 * maxLength;
    bf        = alloc + 5 * maxLength;
    bm        = alloc + 6 * maxLength;
    b4        = alloc + 7 * maxLength;

    b1[0]     = 0;
    b2[0]     = 0;
    bf[0]     = 0;
    bm[0]     = 0;
    b3[0]     = 0;
    b4[0]     = 0;

    bLength   = 0;
    bRetained = 0.0;
    bLabel    = 0;
  };

  ~fastqRecord() {
    delete [] alloc;
  };


  bool          load(FILE *FASTQ1, FILE *FASTQ2) {
    bool  tooShort = false;

    a2[maxLength - 2] = 0;
    b2[maxLength - 2] = 0;

    if (FASTQ1) {
      fgets(a1,      1024, FASTQ1);  chomp(a1);
      fgets(a2, maxLength, FASTQ1);  chomp(a2);
      fgets(a3,      1024, FASTQ1);  chomp(a3);
      fgets(a4, maxLength, FASTQ1);  chomp(a4);

      aLength   = strlen(a2);
      aRetained = 0.0;
      aLabel    = 0;

      if (a2[maxLength - 2] != 0)
        tooShort = true;
    }

    if (FASTQ2) {
      fgets(b1,      1024, FASTQ2);  chomp(b1);
      fgets(b2, maxLength, FASTQ2);  chomp(b2);
      fgets(b3,      1024, FASTQ2);  chomp(b3);
      fgets(b4, maxLength, FASTQ2);  chomp(b4);

      bLength   = strlen(b2);
      bRetained = 0.0;
      bLabel    = 0;

      if (b2[maxLength - 2] != 0)
        tooShort = true;
    }

    if (tooShort) {
      fprintf(stderr, "ERROR: -l too small for reads:\n");
      fprintf(stderr, "       a = '%s'\n", a1);
      fprintf(stderr, "       b = '%s'\n", b1);
      exit(1);
    }

    return(!feof(FASTQ1));
  };


  void          write(FILE *FASTQ1, FILE *FASTQ2) {

    if (FASTQ1)
      fprintf(FASTQ1, "%s fractionRetained=%.3f\n%s\n%s\n%s\n", a1, aRetained, am, a3, a4);

    if (FASTQ2)
      fprintf(FASTQ2, "%s fractionRetained=%.3f\n%s\n%s\n%s\n", b1, bRetained, bm, b3, b4);
  };



public:
  uint32        maxLength;
  char         *alloc;

  char          a1[1024];
  char         *a2;
  char         *af;
  char         *am;
  char          a3[1024];
  char         *a4;
  
  uint32        aLength;
  double        aRetained;
  uint32        aLabel;

  char          b1[1024];
  char         *b2;
  char         *bf;
  char         *bm;
  char          b3[1024];
  char         *b4;

  uint32        bLength;
  double        bRetained;
  uint32        bLabel;
};




class maskGlobal {
public:
  maskGlobal() {
    merName         = NULL;

    seq1Name        = NULL;
    seq2Name        = NULL;

    outPrefix       = NULL;

    merSize         = 0;
    maxLength       = 512;

    existName       = NULL;
    minSize         = 0;
    extend          = 0;
    keepNovel       = false;
    keepConfirmed   = false;

    demote          = false;
    promote         = false;
    discard         = true;

    lowThreshold    = 1. / 3.;
    highThreshold   = 2. / 3.;

    for (uint32 ii=0; ii<1001; ii++)
      scoreHistogram[ii] = 0;

    for (uint32 ii=0; ii<4; ii++)
      for (uint32 jj=0; jj<4; jj++)
        thresholdCounts[ii][jj] = 0.0;

    outputHistogram = NULL;

    exist = NULL;

    FASTQ1      = NULL;
    FASTQ1pipe  = false;

    FASTQ2      = NULL;
    FASTQ2pipe  = false;

    OUTPUT1[0]  = NULL;
    OUTPUT1[1]  = NULL;
    OUTPUT1[2]  = NULL;

    OUTPUT2[0]  = NULL;
    OUTPUT2[1]  = NULL;
    OUTPUT2[2]  = NULL;
  };

  ~maskGlobal() {
  };


public:
  char         *merName;

  char         *seq1Name;
  char         *seq2Name;

  char         *outPrefix;

  uint32        merSize;
  uint32        maxLength;

  char         *existName;
  uint32        minSize;
  uint32        extend;
  bool          keepNovel;
  bool          keepConfirmed;

  bool          demote;
  bool          promote;
  bool          discard;

  double        lowThreshold;
  double        highThreshold;

  uint32        scoreHistogram[1001];
  uint64        thresholdCounts[4][4];

  char         *outputHistogram;

  existDB      *exist;

  FILE         *FASTQ1;
  bool          FASTQ1pipe;

  FILE         *FASTQ2;
  bool          FASTQ2pipe;

  FILE         *OUTPUT1[4];
  FILE         *OUTPUT2[4];
};















//  Masks mers present in the database from the input sequences.  Chains together
//  across small bits of missing mer.

void
printBits(char *S, uint32 Slen, char *found, char *display, const char *label) {
  for (uint32 i=0; i<Slen; i++)
    display[i] = (found[i]) ? '1' : '0';

  display[Slen] = 0;

  fprintf(stdout, "%s\n%s\n", label, display);
}




//  Scan the read for kmers that exist in the DB.  Set a bit for each kmer that exists.
void
buildMask(char *S, uint32 Slen, char *found, bool keepNovel, existDB *exist, uint32 merSize) {
  merStream    MS(new kMerBuilder(merSize),
                  new seqStream(S, Slen),
                  true, true);

  for (uint32 i=0; i<Slen; i++)
    found[i] = false;

  while (MS.nextMer())
    if (exist->exists(MS.theFMer()) || exist->exists(MS.theRMer()))
      found[MS.thePositionInSequence()] = true;
}




//  Searched for isolated 'true' bits, and removes them.  Isolated means fewer than minSize true
//  bits are adjacent.
//
void
removeIsolatedMers(char *S, uint32 Slen, char *found, uint32 minSize) {
  uint32  bgn   = 0;
  uint32  end   = 0;
  bool    inRun = false;

  for (uint32 ii=0; ii<Slen; ii++) {

    //  Start of a run of 'true'.
    if ((found[ii] == true) && (inRun == false)) {
      bgn = ii;
      end = ii;
      inRun = true;
    }

    //  End of run of 'true'.  If small, destroy it.
    if ((found[ii] == false) && (inRun == true)) {
      end = ii;

      if (end - bgn < minSize)
        for (uint32 jj=bgn; jj<end; jj++)
          found[jj] = false;

      inRun = false;
    }
  }
}



//  Convert the mer-start-based mask to a base-covering mask, allowing an extra uncovered
//  'extension' bases in between two blocks to join.
void
convertToBases(char *S, uint32 Slen, char *found, uint32 merSize, uint32 extend) {
  uint32   isMasking = 0;

  for (uint32 ii=0; ii<Slen; ii++) {

    if (found[ii])
      isMasking = merSize;

    //  If the last mer we've found, see if we can extend over bases.
    if ((isMasking == 1) &&
        (extend > 0)) {
      for (uint32 jj=ii; (jj<Slen) && (jj <= ii + extend + 1); jj++) {
        if (found[jj] == true)
          isMasking = jj - ii + 2;
      }
    }

    if (isMasking > 0) {
      found[ii] = true;
      isMasking--;
    }
  }
}



//  Assumes the found[] array represents base-based masking.
//  Returns the fraction of the sequence that is not masked.
double
maskSequence(char *S, uint32 Slen, char *found, bool keepNovel, char *display) {
  uint32  saved  = 0;

  for (uint32 ii=0; ii<Slen; ii++) {
    if (found[ii] == keepNovel) {
      display[ii] = 'n';
    } else {
      display[ii] = S[ii];
      saved++;
    }
  }

  display[Slen] = 0;

  return((double)saved  / Slen);
}








void *
fastqLoader(void *G) {
  maskGlobal   *g = (maskGlobal *)G;
  fastqRecord  *s = new fastqRecord(g->maxLength);

  if (s->load(g->FASTQ1, g->FASTQ2) == false) {
    delete s;
    s = NULL;
  }

  return(s);
}


void
maskWorker(void *G, void *T, void *S) {
  maskGlobal   *g = (maskGlobal *)G;
  //maskThread  *t = (maskThread *)T;
  fastqRecord  *s = (fastqRecord *)S;

  buildMask(s->a2, s->aLength, s->af, g->keepNovel, g->exist, g->merSize);
  buildMask(s->b2, s->bLength, s->bf, g->keepNovel, g->exist, g->merSize);
  //printBits(S, found, display , "INITIAL");

  removeIsolatedMers(s->a2, s->aLength, s->af, g->minSize);
  removeIsolatedMers(s->b2, s->bLength, s->bf, g->minSize);
  //printBits(S, found, display, "ISOLATED REMOVAL");

  convertToBases(s->a2, s->aLength, s->af, g->merSize, g->extend);
  convertToBases(s->b2, s->bLength, s->bf, g->merSize, g->extend);
  //printBits(S, found, display, "BASE COVERAGE");

  s->aRetained = maskSequence(s->a2, s->aLength, s->af, g->keepNovel, s->am);
  s->bRetained = maskSequence(s->b2, s->bLength, s->bf, g->keepNovel, s->bm);

  s->aLabel = (s->aRetained < g->lowThreshold) ? 0 : ((s->aRetained < g->highThreshold) ? 1 : 2);
  s->bLabel = (s->bRetained < g->lowThreshold) ? 0 : ((s->bRetained < g->highThreshold) ? 1 : 2);

  if ((s->aLabel != s->bLabel) && (g->demote)) {
    s->aLabel = MIN(s->aLabel, s->bLabel);
    s->bLabel = MIN(s->aLabel, s->bLabel);
  }

  if ((s->aLabel != s->bLabel) && (g->promote)) {
    s->aLabel = MAX(s->aLabel, s->bLabel);
    s->bLabel = MAX(s->aLabel, s->bLabel);
  }

  if ((s->aLabel != s->bLabel) && (g->discard)) {
    s->aLabel = 3;
    s->bLabel = 3;
  }
}


void
fastqWriter(void *G, void *S) {
  maskGlobal  *g = (maskGlobal *)G;
  fastqRecord *s = (fastqRecord *)S;

  s->write(g->OUTPUT1[s->aLabel], g->OUTPUT2[s->bLabel]);

  g->thresholdCounts[s->aLabel][s->bLabel]++;

  g->scoreHistogram[(uint32)(1000 * s->aRetained)]++;
  g->scoreHistogram[(uint32)(1000 * s->bRetained)]++;

  delete s;
}






FILE *
openInput(char *filename, bool &P) {
  char  C[2 * FILENAME_MAX];
  FILE *F = NULL;
  int32 L = strlen(filename);

  if ((L > 6) && (strcmp(filename + L - 6, ".fastq") == 0)) {
    F = fopen(filename, "r");
    P = false;
  }

  if ((L > 3) && (strcmp(filename + L - 3, ".gz") == 0)) {
    sprintf(C, "gzip -dc %s", filename);
    F = popen(C, "r");
    P = true;
  }

  if ((L > 4) && (strcmp(filename + L - 4, ".bz2") == 0)) {
    sprintf(C, "bzip2 -dc %s", filename);
    F = popen(C, "r");
    P = true;
  }

  if ((L > 3) && (strcmp(filename + L - 3, ".xz") == 0)) {
    sprintf(C, "xz -dc %s", filename);
    F = popen(C, "r");
    P = true;
  }

  return(F);
}

void
closeInput(FILE *F, char *filename, bool P) {
  if (F)
    if (P)
      pclose(F);
    else
      fclose(F);
}



FILE *
openOutput(char *prefix, const char *extension) {
  char  N[FILENAME_MAX];
  FILE *F = NULL;

  sprintf(N, "%s.%s.fastq", prefix, extension);
  F = fopen(N, "w");
  if (errno)
    fprintf(stderr, "ERROR: failed to open '%s': %s\n", N, strerror(errno)), exit(1);

  return(F);
}

void
closeOutput(FILE *F, char *prefix, const char *extension) {
  if (F)
    fclose(F);
}



int
main(int argc, char **argv) {
  maskGlobal   *g = new maskGlobal();

  uint32        numWorkers = 1;
  bool          beVerbose  = false;

  int32 arg=1;
  int32 err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-mdb") == 0) {
      g->merName = argv[++arg];

    } else if (strcmp(argv[arg], "-ms") == 0) {
      g->merSize = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      g->maxLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-edb") == 0) {
      g->existName = argv[++arg];

    } else if (strcmp(argv[arg], "-1") == 0) {
      g->seq1Name = argv[++arg];
    } else if (strcmp(argv[arg], "-2") == 0) {
      g->seq2Name = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      g->outPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      g->minSize = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      g->extend = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      numWorkers = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;

    } else if (strcmp(argv[arg], "-novel") == 0) {
      g->keepNovel = true;
    } else if (strcmp(argv[arg], "-confirmed") == 0) {
      g->keepConfirmed = true;

    } else if (strcmp(argv[arg], "-demote") == 0) {
      g->demote = true;
    } else if (strcmp(argv[arg], "-promote") == 0) {
      g->promote = true;
    } else if (strcmp(argv[arg], "-discard") == 0) {
      g->discard = true;

    } else if (strncmp(argv[arg], "-lowthreshold", 3) == 0) {
      g->lowThreshold = atof(argv[++arg]);
    } else if (strncmp(argv[arg], "-highthreshold", 3) == 0) {
      g->highThreshold = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-h") == 0) {
      g->outputHistogram = argv[++arg];
      //} else if (strcmp(argv[arg], "-o") == 0) {
      //  outputSequence = atoi(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }
  if ((g->keepNovel == false) && (g->keepConfirmed == false))
    err++;
  if (err) {
    fprintf(stderr, "usage: %s [-novel | -confirmed] ...\n", argv[0]);
    fprintf(stderr, "  -mdb mer-database      load masking kmers from meryl 'mer-database'\n");
    fprintf(stderr, "  -ms  mer-size          \n");
    fprintf(stderr, "  -edb exist-database    save masking kmers to 'exist-database' for faster restarts\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -1 in.1.fastq          input reads - fastq, fastq.gz, fastq.bz2 or fastq.xz\n");
    fprintf(stderr, "  -2 in.2.fastq                      - (optional, but if not present, messes up the output classification)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o  out                output reads:\n");
    fprintf(stderr, "                            out.fullymasked.[12].fastq      - reads with below 'lowthreshold' bases retained\n");
    fprintf(stderr, "                            out.partiallymasked.[12].fastq  - reads in between\n");
    fprintf(stderr, "                            out.retained.[12].fastq         - reads with more than 'hightreshold' bases retained\n");
    fprintf(stderr, "                            out.discarded.[12].fastq        - reads with conflicting status\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -m min-size            ignore database hits below this many consecutive kmers (%d)\n", g->minSize);
    fprintf(stderr, "  -e extend-size         extend database hits across this many missing kmers (%d)\n", g->extend);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -novel                 RETAIN novel sequence not present in the database\n");
    fprintf(stderr, "  -confirmed             RETAIN confirmed sequence present in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -promote               promote the lesser RETAINED read to the status of the more RETAINED read\n");
    fprintf(stderr, "                           read1=fullymasked and read2=partiallymasked -> both are partiallymasked\n");
    fprintf(stderr, "  -demote                demote the more RETAINED read to the status of the lesser RETAINED read\n");
    fprintf(stderr, "                           read1=fullymasked and read2=partiallymasked -> both are fullymasked\n");
    fprintf(stderr, "  -discard               discard pairs with conflicting status (DEFAULT)\n");
    fprintf(stderr, "                           read1=fullymasked and read2=partiallymasked -> both are discarded\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "stats on stderr, number of sequences with amount RETAINED:\n");
    fprintf(stderr, "  -lowthreshold t        (%.4f)\n", g->lowThreshold);
    fprintf(stderr, "  -highthreshold t       (%.4f)\n", g->highThreshold);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -h histogram           write a histogram of the amount of sequence RETAINED\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t t                   use 't' compute threads\n");
    fprintf(stderr, "  -v                     show progress\n");

    if ((g->keepNovel == false) && (g->keepConfirmed == false))
      fprintf(stderr, "ERROR: exactly one of -novel and -confirmed must be supplied.\n");

    exit(1);
  }

  //  Open inputs

  g->FASTQ1     = openInput(g->seq1Name, g->FASTQ1pipe);
  g->FASTQ2     = openInput(g->seq2Name, g->FASTQ2pipe);

  g->OUTPUT1[0] = openOutput(g->outPrefix, "fullymasked.1");
  g->OUTPUT1[1] = openOutput(g->outPrefix, "partiallymasked.1");
  g->OUTPUT1[2] = openOutput(g->outPrefix, "retained.1");
  g->OUTPUT1[3] = openOutput(g->outPrefix, "discarded.1");

  g->OUTPUT2[0] = openOutput(g->outPrefix, "fullymasked.2");
  g->OUTPUT2[1] = openOutput(g->outPrefix, "partiallymasked.2");
  g->OUTPUT2[2] = openOutput(g->outPrefix, "retained.2");
  g->OUTPUT2[3] = openOutput(g->outPrefix, "discarded.2");


  //  Load data

  if ((g->existName != NULL) && (fileExists(g->existName))) {
    if (beVerbose)
      fprintf(stderr, "Load existDB existName='%s'.\n", g->existName);
    g->exist = new existDB(g->existName);

  } else {
    if (beVerbose)
      fprintf(stderr, "Build existDB from merName='%s'.\n", g->merName);
    g->exist = new existDB(g->merName, g->merSize, existDBnoFlags, 0, ~uint32ZERO);

    if (g->existName != NULL) {
      if (beVerbose)
        fprintf(stderr, "Save existDB into existName='%s'.\n", g->existName);
      g->exist->saveState(g->existName);
    }
  }

  //  Process!

  sweatShop *ss = new sweatShop(fastqLoader, maskWorker, fastqWriter);

  ss->setNumberOfWorkers(numWorkers);

  ss->setWorkerBatchSize(1024);

  ss->setLoaderQueueSize(numWorkers * 81920);
  ss->setWriterQueueSize(numWorkers * 81920);

  ss->run(g, beVerbose);

  closeInput(g->FASTQ1, g->seq1Name, g->FASTQ1pipe);
  closeInput(g->FASTQ2, g->seq1Name, g->FASTQ2pipe);

  closeOutput(g->OUTPUT1[0], g->outPrefix, "fulymasked.1");
  closeOutput(g->OUTPUT1[1], g->outPrefix, "partiallymasked.1");
  closeOutput(g->OUTPUT1[2], g->outPrefix, "retained.1");
  closeOutput(g->OUTPUT1[3], g->outPrefix, "discarded.1");

  closeOutput(g->OUTPUT2[0], g->outPrefix, "fulymasked.2");
  closeOutput(g->OUTPUT2[1], g->outPrefix, "partiallymasked.2");
  closeOutput(g->OUTPUT2[2], g->outPrefix, "retained.2");
  closeOutput(g->OUTPUT2[3], g->outPrefix, "discarded.2");

  fprintf(stderr, "            bBelow    bNormal   bHigh     bDiscarded\n");
  fprintf(stderr, "aBelow      %8lu  %8lu  %8lu  %8lu\n", g->thresholdCounts[0][0], g->thresholdCounts[0][1], g->thresholdCounts[0][2], g->thresholdCounts[0][3]);
  fprintf(stderr, "aNormal     %8lu  %8lu  %8lu  %8lu\n", g->thresholdCounts[1][0], g->thresholdCounts[1][1], g->thresholdCounts[1][2], g->thresholdCounts[1][3]);
  fprintf(stderr, "aHigh       %8lu  %8lu  %8lu  %8lu\n", g->thresholdCounts[2][0], g->thresholdCounts[2][1], g->thresholdCounts[2][2], g->thresholdCounts[2][3]);
  fprintf(stderr, "aDiscarded  %8lu  %8lu  %8lu  %8lu\n", g->thresholdCounts[3][0], g->thresholdCounts[3][1], g->thresholdCounts[3][2], g->thresholdCounts[3][3]);

  if (g->outputHistogram != NULL) {
    FILE *H = fopen(g->outputHistogram, "w");

    fprintf(H, "# amount of sequence retained\n");
    for (uint32 i=0; i<1001; i++)
      if (g->scoreHistogram[i] > 0)
        fprintf(H, "%.4f\t%u\n", i / 1000.0, g->scoreHistogram[i]);

    fclose(H);
  }

  delete g;

  exit(0);
}
