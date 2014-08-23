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
  char         *seqName         = NULL;
  uint32        merSize         = 0;
  char         *existName       = NULL;
  uint32        minSize         = 0;
  uint32        extend          = 0;
  bool          keepNovel       = false;

  uint32        belowLow        = 0;
  double        lowThreshold    = 1. / 3.;

  uint32        aboveHigh       = 0;
  double        highThreshold   = 2. / 3.;

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

    } else if (strcmp(argv[arg], "-s") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      minSize = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      extend = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-novel") == 0) {
      //  Retains kmers that do NOT exist in the DB.
      keepNovel = true;

    } else if (strcmp(argv[arg], "-confirmed") == 0) {
      //  Retains kmers that exist in the DB.
      keepNovel = false;

    } else if (strcmp(argv[arg], "-lowthreshold") == 0) {
      lowThreshold = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-highthreshold") == 0) {
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
  if (err) {
  }

  //  Open inputs

  char *CMD   = new char [2 * FILENAME_MAX];
  sprintf(CMD, "gzip -dc %s", seqName);
  FILE *FASTQ = popen(CMD, "r");
  if (errno)
    fprintf(stderr, "ERROR: failed to open '%s': %s\n", seqName, strerror(errno));

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

  //seqCache     *F     = new seqCache(seqName);
  //seqInCore    *S     = 0L;

  uint32    allocLen = 1048576;
  bool     *found    = new bool [allocLen];
  char     *display  = new char [allocLen];

  char     *a1 = new char [1024];
  char     *a2 = new char [allocLen];
  char     *a3 = new char [1024];
  char     *a4 = new char [allocLen];

  uint32    al = strlen(a2);

  uint32   scoreHistogram[1001] = { 0 };

  fgets(a1,     1024, FASTQ);  chomp(a1);
  fgets(a2, allocLen, FASTQ);  chomp(a2);
  fgets(a3,     1024, FASTQ);  chomp(a3);
  fgets(a4, allocLen, FASTQ);  chomp(a4);

  al = strlen(a2);

  fprintf(stderr, "Begin.\n");

  while (!feof(FASTQ)) {
    //(S = F->getSequenceInCore()) != 0L)

    buildMask(a2, al, found, keepNovel, exist, merSize);
    //printBits(S, found, display , "INITIAL");

    removeIsolatedMers(a2, al, found, minSize);
    //printBits(S, found, display, "ISOLATED REMOVAL");

    convertToBases(a2, al, found, merSize, extend);
    //printBits(S, found, display, "BASE COVERAGE");

    double fractionRetained = maskSequence(a2, al, found, keepNovel, display);

    //fprintf(stdout, "%.3f\t%s\n", fractionRetained, S->header());
    fprintf(stdout, "%s fractionRetained=%.3f\n%s\n%s\n%s\n",
            a1, fractionRetained,
            display,
            a3,
            a4);

    if (fractionRetained < lowThreshold)
      belowLow++;

    if (fractionRetained > highThreshold)
      aboveHigh++;

    scoreHistogram[(uint32)(1000 * fractionRetained)]++;

    //delete S;

    fgets(a1,     1024, FASTQ);  chomp(a1);
    fgets(a2, allocLen, FASTQ);  chomp(a2);
    fgets(a3,     1024, FASTQ);  chomp(a3);
    fgets(a4, allocLen, FASTQ);  chomp(a4);

    al = strlen(a2);
  }

  fprintf(stderr, "below %8.6f:   %u\n", lowThreshold, belowLow);
  fprintf(stderr, "above %8.6f:   %u\n", highThreshold, aboveHigh);

  delete [] display;
  delete [] found;

  if (outputHistogram != NULL) {
    FILE *H = fopen(outputHistogram, "w");

    fprintf(H, "# amount of sequence not masked\n");
    for (uint32 i=0; i<1001; i++)
      if (scoreHistogram[i] > 0)
        fprintf(H, "%.4f\t%u\n", i / 1000.0, scoreHistogram[i]);

    fclose(H);
  }

  exit(0);
}
