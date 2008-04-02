#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bio++.H"
#include "libmeryl.H"

#define MAX_COVERAGE 51

class merMaskedSequence {
public:
  merMaskedSequence(char *fastaName, char *merylName) {
    _numSeq = 0;
    _seqLen = 0L;
    _masking = 0L;
    _merSize = 0;

    strcpy(_fastaName, fastaName);
    strcpy(_merylName, merylName);

    strcpy(_maskMersName, merylName);
    strcat(_maskMersName, ".maskMers");

    if (fileExists(_maskMersName))
      loadMasking();
    else
      buildMasking();
  };
  ~merMaskedSequence() {
    delete [] _seqLen;
    for (u32bit i=0; i<_numSeq; i++)
      delete [] _masking[i];
    delete [] _masking;
  };

public:
  u32bit     numSeq(void)                { return(_numSeq); };
  u32bit     seqLen(u32bit i)            { return(_seqLen[i]); };
  char       masking(u32bit s, u32bit p) { return(_masking[s][p]); };

  char      *merylName(void)             { return(_merylName); };

  u32bit     merSize(void)               { return(_merSize); };

private:
  void       loadMasking(void);      //  Read the masking from the saved file
  void       saveMasking(void);      //  Write the masking to a file
  void       buildMasking(void);     //  Read the mers to build the masking

  u32bit    _numSeq;
  u32bit   *_seqLen;
  char    **_masking;

  u32bit    _merSize;

  char      _fastaName[FILENAME_MAX];
  char      _merylName[FILENAME_MAX];
  char      _maskMersName[FILENAME_MAX];
};


void
merMaskedSequence::loadMasking(void) {
  FILE  *maskMersFile = fopen(_maskMersName, "r");

  fread(&_numSeq,  sizeof(u32bit), 1, maskMersFile);
  fread(&_merSize, sizeof(u32bit), 1, maskMersFile);

  _seqLen = new u32bit [_numSeq];
  _masking = new char * [_numSeq];

  fprintf(stderr, u32bitFMT" sequences in '%s'\n", _numSeq, _fastaName);

  fread( _seqLen,   sizeof(u32bit), _numSeq, maskMersFile);

  for (u32bit i=0; i<_numSeq; i++) {
    _masking[i] = new char [_seqLen[i]];
    memset(_masking[i], 'g', sizeof(char) * _seqLen[i]);
  }

  for (u32bit i=0; i<_numSeq; i++)
    fread(_masking[i], sizeof(char), _seqLen[i], maskMersFile);

  fclose(maskMersFile);
}


void
merMaskedSequence::saveMasking(void) {
  FILE  *maskMersFile = fopen(_maskMersName, "w");

  fwrite(&_numSeq,  sizeof(u32bit), 1,       maskMersFile);
  fwrite(&_merSize, sizeof(u32bit), 1,       maskMersFile);
  fwrite( _seqLen,  sizeof(u32bit), _numSeq, maskMersFile);

  for (u32bit i=0; i<_numSeq; i++)
    fwrite(_masking[i], sizeof(char), _seqLen[i], maskMersFile);

  fclose(maskMersFile);
}


void
merMaskedSequence::buildMasking(void) {
  seqFile       *F   = openSeqFile(_fastaName);
  seqStream     *STR = new seqStream(F, true);
  //seqStore      *STO = new seqStore(fastaName, STR);

  _numSeq  = STR->numberOfSequences();
  _seqLen  = new u32bit [_numSeq];
  _masking = new char * [_numSeq];
  _merSize = 0;

  fprintf(stderr, u32bitFMT" sequences in '%s'\n", _numSeq, _fastaName);

  for (u32bit i=0; i<_numSeq; i++) {
    _seqLen[i] = STR->lengthOf(i);
    _masking[i] = new char [_seqLen[i]];
    memset(_masking[i], 'g', sizeof(char) * _seqLen[i]);
  }

  //  g -> gap in sequence
  //  u -> unique mer
  //  r -> repeat mer

  merylStreamReader *MS = new merylStreamReader(_merylName);
  speedCounter      *CT = new speedCounter(" Masking mers in sequence: %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);

  _merSize = MS->merSize();

  while (MS->nextMer()) {
    //fprintf(stderr, "mer count="u64bitFMT" pos="u32bitFMT"\n", MS->theCount(), MS->getPosition(0));

    if (MS->theCount() == 1) {
      u32bit p = MS->getPosition(0);
      u32bit s = STR->sequenceNumberOfPosition(p);
      p -= STR->startOf(s);

      _masking[s][p] = 'u';
    } else {
      for (u32bit i=0; i<MS->theCount(); i++) {
        u32bit p = MS->getPosition(i);
        u32bit s = STR->sequenceNumberOfPosition(p);
        p -= STR->startOf(s);

        _masking[s][p] = 'r';
      }
    }

    CT->tick();
  }

  delete CT;

  delete MS;

  //delete STO;
  delete STR;
  delete F;

  saveMasking();
}


void
computeDensity(merMaskedSequence *S) {
  char    outputName[FILENAME_MAX];
  FILE   *outputFile;
  u32bit  windowSizeMax = 1000;

  for (u32bit s=0; s<S->numSeq(); s++) {
    sprintf(outputName, "%s.density.seq"u32bitFMTW(02), S->merylName(), s);
    outputFile = fopen(outputName, "w");

    fprintf(outputFile, "#window\tunique\trepeat\tgaps\n");

    //  Not the most efficient, but good enough for us right now.

    for (u32bit p=0; p<S->seqLen(s); ) {
      u32bit  windowSize    = 0;
      u32bit  uniqueSum     = 0;
      u32bit  repeatSum     = 0;
      u32bit  gapSum        = 0;

      while ((windowSize < windowSizeMax) &&
             (p < S->seqLen(s))) {
        char m = S->masking(s, p);

        if (m == 'u')  uniqueSum++;
        if (m == 'g')  gapSum++;
        if (m == 'r')  repeatSum++;

        windowSize++;
        p++;
      }

      fprintf(outputFile, u32bitFMT"\t%f\t%f\t%f\n",
              p - windowSize,
              (double)uniqueSum / windowSize,
              (double)repeatSum / windowSize,
              (double)gapSum    / windowSize);
    }

    fclose(outputFile);
  }
}


void
computeMateRescue(merMaskedSequence *S) {
  char    outputName[FILENAME_MAX];
  FILE   *outputFile;

  s32bit  mean      = 3000;
  s32bit  stddev    = 500;

  assert(mean > 3 * stddev);

  //  +1 because a value of 1.0 == histogramLen
  //
  u32bit   histogramLen = 1000;
  u32bit  *histogram    = new u32bit [histogramLen + 1];

  //  For each 'r' mer, compute the number of 'u' mers
  //  that are within some mean +- stddev range.
  //
  //  We count for two blocks:
  //
  //            |   <- mean ->   |  <- mean ->    |
  //  ---[block1]---------------mer---------------[block2]---
  //
  //  Once we know that, we can compute the probability that
  //  a repeat mer can be rescued.
  //
  //  p1 = uniq/total   -- for 1 X coverage
  //  pn = 1 - (1-p1)^n -- for n X coverage
  //

  double *normal     = 0L;
  s32bit  normalZero = 0;

  {
    double  a = 1.0 / (stddev * sqrt(2 * M_PI));
    double  c = 2 * stddev * stddev;

    s32bit  b1l = (s32bit)floor(-3 * stddev);
    s32bit  b1h = (s32bit)ceil ( 3 * stddev);

    normal     = new double [b1h - b1l + 1];
    normalZero = -b1l;

    for (s32bit l=0; l<b1h - b1l + 1; l++)
      normal[l] = 0.0;

    for (s32bit l=b1l; l<b1h; l++)
      normal[l + normalZero] = a * exp(- l*l / c);
  }



  for (u32bit s=0; s<S->numSeq(); s++) {
    double  numRR[MAX_COVERAGE] = {0};  //  num repeats rescued (expected) for [] X coverage

    memset(histogram, 0, sizeof(u32bit) * (histogramLen + 1));

    //speedCounter *CT = new speedCounter(" Examining repeats: %7.2f Krepeats -- %5.2f Krepeats/second\r", 1000.0, 0x1ffff, true);

    u32bit  numRT = 0;  //  num repeats total

    for (s32bit p=0; p<S->seqLen(s); p++) {
      if (S->masking(s, p) == 'r') {
        numRT++;

        s32bit b1l = (s32bit)floor(p - mean - 3 * stddev);
        s32bit b1h = (s32bit)ceil (p - mean + 3 * stddev);

        s32bit b2l = (s32bit)floor(p + mean - 3 * stddev);
        s32bit b2h = (s32bit)ceil (p + mean + 3 * stddev);

        if (b1l < 0)            b1l = 0;
        if (b1h < 0)            b1h = 0;
        if (b1h > S->seqLen(s)) b1h = S->seqLen(s);

        if (b2l < 0)            b2l = 0;
        if (b2h > S->seqLen(s)) b2h = S->seqLen(s);
        if (b2l > S->seqLen(s)) b2l = S->seqLen(s);

        double  pRescue = 0.0;
        double  pRepeat = 0.0;

        //fprintf(stderr, "b1: %d-%d  b2:%d-%d\n", b1l, b1h, b2l, b2h);

        for (s32bit l=b1l; l<b1h; l++) {
          if (S->masking(s, l) == 'u')
            pRescue += normal[l - p + mean + normalZero];
        }

        for (s32bit l=b2l; l<b2h; l++) {
          if (S->masking(s, l) == 'u')
            pRescue += normal[l - p - mean + normalZero];
        }

        //  We're summing over two distributions.
        pRescue /= 2.0;
        pRepeat  = 1.0 - pRescue;

        for (u32bit x=1; x<MAX_COVERAGE; x++) {
          numRR[x] += 1 - pRepeat;
          pRepeat *= (1.0 - pRescue);
        }

        //  Histogram of the probability of rescuing a repeat
        histogram[(s32bit)floor(histogramLen * pRescue)]++;

        //CT->tick();
      }
    }

    //delete CT;

    sprintf(outputName, "%s.mateRescue.seq"u32bitFMTW(02)",out", S->merylName(), s);
    outputFile = fopen(outputName, "w");

    fprintf(outputFile, "seqIID\tmerSize\t#totalRepeats\texpectedRescue\tXcoverage\n");

    for (u32bit x=1; x<MAX_COVERAGE; x++)
      fprintf(outputFile, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%.0f\t"u32bitFMT"\n",
              s, S->merSize(), numRT, numRR[x], x);

    fclose(outputFile);

    sprintf(outputName, "%s.mateRescue.seq"u32bitFMTW(02)".histogram", S->merylName(), s);
    outputFile = fopen(outputName, "w");

    fprintf(outputFile, "#pRescue\tnumRepeats\tfraction_repeats_higher_probability\n");

    u32bit  cumulative = 0;

    for (u32bit i=0; i<=histogramLen; i++) {
      cumulative += histogram[i];
      fprintf(outputFile, "%f\t"u32bitFMT"\t%f\n",
              i / (double)histogramLen,
              histogram[i],
              1.0 - cumulative / (double)numRT);
    }

    fclose(outputFile);
  }

  delete [] histogram;
}



int
main(int argc, char **argv) {
  char     *merylName  = 0L;
  char     *fastaName  = 0L;

  bool      doDensity  = false;
  bool      doRescue   = false;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-mers") == 0) {
      merylName = argv[++arg];

    } else if (strcmp(argv[arg], "-seq") == 0) {
      fastaName = argv[++arg];

    } else if (strcmp(argv[arg], "-d") == 0) {
      doDensity = true;

    } else if (strcmp(argv[arg], "-r") == 0) {
      doRescue = true;

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err) || (merylName == 0L) || (fastaName == 0L)) {
    fprintf(stderr, "usage: %s -mers mers -seq fasta > output\n", argv[0]);
    exit(1);
  }

  merMaskedSequence *S = new merMaskedSequence(fastaName, merylName);

  if (doDensity)
    computeDensity(S);

  if (doRescue)
    computeMateRescue(S);

  return(0);
}
