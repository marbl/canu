#include <stdio.h>
#include <stdlib.h>

#include "util++.H"
#include "bio++.H"
#include "libmeryl.H"
#include "existDB.H"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

//  Masks mers present in the database from the input sequences.  Chains together
//  across small bits of missing mer.

void
printBits(char *S, uint32 Slen, bool *found, char *display, const char *label) {
  for (uint32 i=0; i<Slen; i++)
    display[i] = (found[i]) ? '1' : '0';

  display[Slen] = 0;

  fprintf(stdout, "%s\n%s\n", label, display);
}




//  Scan the read for kmers that exist in the DB.  Set a bit for each kmer that exists.
void
buildMask(char *S, uint32 Slen, bool *found, bool keepNovel, existDB *exist, uint32 merSize) {
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
removeIsolatedMers(char *S, uint32 Slen, bool *found, uint32 minSize) {
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
convertToBases(char *S, uint32 Slen, bool *found, uint32 merSize, uint32 extend) {
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
maskSequence(char *S, uint32 Slen, bool *found, bool keepNovel, char *display) {
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






int
main(int argc, char **argv) {
  char         *merName         = NULL;

  char         *seq1Name        = NULL;
  char         *seq2Name        = NULL;

  char         *outPrefix       = NULL;

  uint32        merSize         = 0;
  char         *existName       = NULL;
  uint32        minSize         = 0;
  uint32        extend          = 0;
  bool          keepNovel       = false;
  bool          keepConfirmed   = false;

  double        lowThreshold    = 1. / 3.;
  double        highThreshold   = 2. / 3.;

  //  Save threshold stats.  For each read, we can label (0) below threshold, (1) in between, (2) above threshold
  //
  uint64        thresholdCounts[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };

  char         *outputHistogram = NULL;

  int32 arg=1;
  int32 err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-mdb") == 0) {
      merName = argv[++arg];

    } else if (strcmp(argv[arg], "-ms") == 0) {
      merSize = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-edb") == 0) {
      existName = argv[++arg];

    } else if (strcmp(argv[arg], "-1") == 0) {
      seq1Name = argv[++arg];
    } else if (strcmp(argv[arg], "-2") == 0) {
      seq2Name = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      minSize = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      extend = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-novel") == 0) {
      //  Retains kmers that do NOT exist in the DB.
      keepNovel = true;

    } else if (strcmp(argv[arg], "-confirmed") == 0) {
      //  Retains kmers that exist in the DB.
      keepConfirmed = true;

    } else if (strncmp(argv[arg], "-lowthreshold", 3) == 0) {
      lowThreshold = atof(argv[++arg]);
    } else if (strncmp(argv[arg], "-highthreshold", 3) == 0) {
      highThreshold = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-h") == 0) {
      outputHistogram = argv[++arg];
    //} else if (strcmp(argv[arg], "-o") == 0) {
    //  outputSequence = atoi(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }
  if ((keepNovel == false) && (keepConfirmed == false))
    err++;
  if (err) {
    fprintf(stderr, "usage: %s [-novel | -confirmed] ...\n", argv[0]);
    fprintf(stderr, "  -mdb mer-database      load masking kmers from meryl 'mer-database'\n");
    fprintf(stderr, "  -ms  mer-size          \n");
    fprintf(stderr, "  -edb exist-database    save masking kmers to 'exist-database' for faster restarts\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -1 in.1.fastq.gz       input reads - must be .gz!\n");
    fprintf(stderr, "  -2 in.2.fastq.gz                   - (optional)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o  out                output reads:\n");
    fprintf(stderr, "                            out.fullymasked.[12].fastq      - reads with below 'lowthreshold' bases retained\n");
    fprintf(stderr, "                            out.partiallymasked.[12].fastq  - reads in between\n");
    fprintf(stderr, "                            out.retained.[12].fastq         - reads with more than 'hightreshold' bases retained\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -m min-size            ignore database hits below this many consecutive kmers (%d)\n", minSize);
    fprintf(stderr, "  -e extend-size         extend database hits across this many missing kmers (%d)\n", extend);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -novel                 RETAIN novel sequence not present in the database\n");
    fprintf(stderr, "  -confirmed             RETAIN confirmed sequence present in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "stats on stderr, number of sequences with amount RETAINED:\n");
    fprintf(stderr, "  -lowthreshold t        (%.4f)\n", lowThreshold);
    fprintf(stderr, "  -highthreshold t       (%.4f)\n", highThreshold);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -h histogram           write a histogram of the amount of sequence RETAINED\n");
    fprintf(stderr, "\n");

    if ((keepNovel == false) && (keepConfirmed == false))
      fprintf(stderr, "ERROR: exactly one of -novel and -confirmed must be supplied.\n");

    exit(1);
  }

  //  Open inputs

  
  FILE *FASTQ1      = NULL;
  bool  FASTQ1pipe  = false;

  FILE *FASTQ2      = NULL;
  bool  FASTQ2pipe  = false;

  FILE *OUTPUT1[3]  = { NULL, NULL, NULL };
  FILE *OUTPUT2[3]  = { NULL, NULL, NULL };


  if (seq1Name) {
    char *CMD   = new char [2 * FILENAME_MAX];
    int32 LEN   = strlen(seq1Name);

    sprintf(CMD, "gzip -dc %s", seq1Name);
    FASTQ1 = popen(CMD, "r");
    if (errno)
      fprintf(stderr, "ERROR: failed to open -1 '%s': %s\n", seq1Name, strerror(errno));

    sprintf(CMD, "%s.fullymasked.1.fastq", outPrefix);
    OUTPUT1[0] = fopen(CMD, "w");
    if (errno)
      fprintf(stderr, "ERROR: failed to open 1[0] '%s': %s\n", CMD, strerror(errno));

    sprintf(CMD, "%s.partiallymasked.1.fastq", outPrefix);
    OUTPUT1[1] = fopen(CMD, "w");
    if (errno)
      fprintf(stderr, "ERROR: failed to open 1[1] '%s': %s\n", CMD, strerror(errno));

    sprintf(CMD, "%s.retained.1.fastq", outPrefix);
    OUTPUT1[2] = fopen(CMD, "w");
    if (errno)
      fprintf(stderr, "ERROR: failed to open 1[2] '%s': %s\n", CMD, strerror(errno));
  }

  if (seq2Name) {
    char *CMD   = new char [2 * FILENAME_MAX];
    int32 LEN   = strlen(seq1Name);

#if 0
    if ((LEN > 6) && (strcmp(seq2Name + LEN - 6, ".fastq") == 0)) {
      FASTQ2     = fopen(seq2Name, "r");
      FASTQ2pipe = false;
    }
    if ((LEN > 3) && (strcmp(seq2Name + LEN - 3, ".gz") == 0)) {
      sprintf(CMD, "gzip -dc %s", seq2Name);
      FASTQ2     = popen(CMD, "r");
      FASTQ2pipe = true;
    }
    if ((LEN > 4) && (strcmp(seq2Name + LEN - 4, ".bz2") == 0)) {
      sprintf(CMD, "bzip2 -dc %s", seq2Name);
      FASTQ2     = popen(CMD, "r");
      FASTQ2pipe = true;
    }
    if ((LEN > 3) && (strcmp(seq2Name + LEN - 3, ".xz") == 0)) {
      sprintf(CMD, "xz -dc %s", seq2Name);
      FASTQ2     = popen(CMD, "r");
      FASTQ2pipe = true;
    }
#endif

    sprintf(CMD, "gzip -dc %s", seq2Name);
    FASTQ2 = popen(CMD, "r");
    if (errno)
      fprintf(stderr, "ERROR: failed to open -2 '%s': %s\n", seq2Name, strerror(errno));

    sprintf(CMD, "%s.fullymasked.2.fastq", outPrefix);
    OUTPUT2[0] = fopen(CMD, "w");
    if (errno)
      fprintf(stderr, "ERROR: failed to open 2[[0] '%s': %s\n", CMD, strerror(errno));

    sprintf(CMD, "%s.partiallymasked.2.fastq", outPrefix);
    OUTPUT2[1] = fopen(CMD, "w");
    if (errno)
      fprintf(stderr, "ERROR: failed to open 2[1] '%s': %s\n", CMD, strerror(errno));

    sprintf(CMD, "%s.retained.2.fastq", outPrefix);
    OUTPUT2[2] = fopen(CMD, "w");
    if (errno)
      fprintf(stderr, "ERROR: failed to open 2[2] '%s': %s\n", CMD, strerror(errno));
  }

  //  Load data

  existDB *exist = 0L;

  if ((existName != NULL) && (fileExists(existName))) {
    fprintf(stderr, "Load existDB existName='%s'.\n", existName);
    exist = new existDB(existName);

  } else {
    fprintf(stderr, "Build existDB from merName='%s'.\n", merName);
    fprintf(stderr, "Save existDB into existName='%s'.\n", existName);
    exist = new existDB(merName, merSize, existDBnoFlags, 0, ~uint32ZERO);

    if (existName != NULL)
      exist->saveState(existName);
  }

  uint32        allocLen = 1048576;

  bool         *found1    = new bool [allocLen];
  char         *display1  = new char [allocLen];

  char         *a1 = new char [1024];
  char         *a2 = new char [allocLen];
  char         *a3 = new char [1024];
  char         *a4 = new char [allocLen];

  bool         *found2    = new bool [allocLen];
  char         *display2  = new char [allocLen];

  char         *b1 = new char [1024];
  char         *b2 = new char [allocLen];
  char         *b3 = new char [1024];
  char         *b4 = new char [allocLen];

  uint32        al = 0;
  uint32        bl = 0;

  uint32        scoreHistogram[1001] = { 0 };

  speedCounter  C("    %7.2f reads -- %5.2f reads/second\r", 1.0, 0x1ffff, true);

  if (FASTQ1) {
    fgets(a1,     1024, FASTQ1);  chomp(a1);
    fgets(a2, allocLen, FASTQ1);  chomp(a2);
    fgets(a3,     1024, FASTQ1);  chomp(a3);
    fgets(a4, allocLen, FASTQ1);  chomp(a4);

    al = strlen(a2);
  }

  if (FASTQ2) {
    fgets(b1,     1024, FASTQ2);  chomp(b1);
    fgets(b2, allocLen, FASTQ2);  chomp(b2);
    fgets(b3,     1024, FASTQ2);  chomp(b3);
    fgets(b4, allocLen, FASTQ2);  chomp(b4);

    bl = strlen(b2);
  }

  fprintf(stderr, "Begin.\n");

  while (((FASTQ1 != NULL) || (FASTQ2 != NULL)) &&
         ((FASTQ1 == NULL) || (feof(FASTQ1) == false)) &&
         ((FASTQ2 == NULL) || (feof(FASTQ2) == false))) {

    buildMask(a2, al, found1, keepNovel, exist, merSize);
    buildMask(b2, bl, found2, keepNovel, exist, merSize);
    //printBits(S, found, display , "INITIAL");

    removeIsolatedMers(a2, al, found1, minSize);
    removeIsolatedMers(b2, bl, found2, minSize);
    //printBits(S, found, display, "ISOLATED REMOVAL");

    convertToBases(a2, al, found1, merSize, extend);
    convertToBases(b2, bl, found2, merSize, extend);
    //printBits(S, found, display, "BASE COVERAGE");

    double fractionRetained1 = maskSequence(a2, al, found1, keepNovel, display1);
    double fractionRetained2 = maskSequence(b2, bl, found2, keepNovel, display2);

    uint32  aLabel = (fractionRetained1 < lowThreshold) ? 0 : ((fractionRetained1 < highThreshold) ? 1 : 2);
    uint32  bLabel = (fractionRetained2 < lowThreshold) ? 0 : ((fractionRetained2 < highThreshold) ? 1 : 2);

    if (OUTPUT1[aLabel])
      fprintf(OUTPUT1[aLabel], "%s fractionRetained=%.3f\n%s\n%s\n%s\n", a1, fractionRetained1, display1, a3, a4);

    if (OUTPUT2[bLabel])
      fprintf(OUTPUT2[bLabel], "%s fractionRetained=%.3f\n%s\n%s\n%s\n", b1, fractionRetained2, display2, b3, b4);

    thresholdCounts[aLabel][bLabel]++;

    scoreHistogram[(uint32)(1000 * fractionRetained1)]++;
    scoreHistogram[(uint32)(1000 * fractionRetained2)]++;

    //  Load next.

    if (FASTQ1) {
      fgets(a1,     1024, FASTQ1);  chomp(a1);
      fgets(a2, allocLen, FASTQ1);  chomp(a2);
      fgets(a3,     1024, FASTQ1);  chomp(a3);
      fgets(a4, allocLen, FASTQ1);  chomp(a4);

      al = strlen(a2);
    }

    if (FASTQ2) {
      fgets(b1,     1024, FASTQ2);  chomp(b1);
      fgets(b2, allocLen, FASTQ2);  chomp(b2);
      fgets(b3,     1024, FASTQ2);  chomp(b3);
      fgets(b4, allocLen, FASTQ2);  chomp(b4);

      bl = strlen(b2);
    }

    C.tick();
  }

  C.finish();

  pclose(FASTQ1);
  pclose(FASTQ2);

  fclose(OUTPUT1[0]);
  fclose(OUTPUT1[1]);
  fclose(OUTPUT1[2]);

  fclose(OUTPUT2[0]);
  fclose(OUTPUT2[1]);
  fclose(OUTPUT2[2]);

  fprintf(stderr, "         bBelow    bNormal   bHigh\n");
  fprintf(stderr, "aBelow   %8lu  %8lu  %8lu\n", thresholdCounts[0][0], thresholdCounts[0][1], thresholdCounts[0][2]);
  fprintf(stderr, "aNormal  %8lu  %8lu  %8lu\n", thresholdCounts[1][0], thresholdCounts[1][1], thresholdCounts[1][2]);
  fprintf(stderr, "aHigh    %8lu  %8lu  %8lu\n", thresholdCounts[2][0], thresholdCounts[2][1], thresholdCounts[2][2]);
         

  delete [] display1;
  delete [] display2;

  delete [] found1;
  delete [] found2;

  if (outputHistogram != NULL) {
    FILE *H = fopen(outputHistogram, "w");

    fprintf(H, "# amount of sequence retained\n");
    for (uint32 i=0; i<1001; i++)
      if (scoreHistogram[i] > 0)
        fprintf(H, "%.4f\t%u\n", i / 1000.0, scoreHistogram[i]);

    fclose(H);
  }

  exit(0);
}
