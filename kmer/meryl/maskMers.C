#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bio++.H"
#include "seqStream.H"
#include "libmeryl.H"

#include <algorithm>

#define MAX_COVERAGE 51

class mateRescueData {
public:
  mateRescueData() {
    _mean       = 0;
    _stddev     = 0;
    _coverage   = 0;
    _normal     = 0L;
    _normalZero = 0;
  };

  void init(s32bit mean_, s32bit stddev_, u32bit coverage_) {
    _mean      = mean_;
    _stddev    = stddev_;
    _coverage  = coverage_;

    assert(_mean > 3 * _stddev);

    double  a = 1.0 / (_stddev * sqrt(2 * M_PI));
    double  c = 2 * _stddev * _stddev;

    s32bit  b1l = (s32bit)floor(-3 * _stddev);
    s32bit  b1h = (s32bit)ceil ( 3 * _stddev);

    _normal     = new double [b1h - b1l + 1];
    _normalZero = -b1l;

    for (s32bit l=0; l<b1h - b1l + 1; l++)
      _normal[l] = 0.0;

    for (s32bit l=b1l; l<b1h; l++)
      _normal[l + _normalZero] = a * exp(- l*l / c);
  };
  ~mateRescueData() {
  };

  s32bit   mean(void)       { return(_mean); };
  s32bit   stddev(void)     { return(_stddev); };
  u32bit   coverage(void)   { return(_coverage); };

  double   normal(s32bit p) { return(_normal[p + _normalZero]); };

private:
  s32bit  _mean;
  s32bit  _stddev;
  u32bit  _coverage;

  double *_normal;
  s32bit  _normalZero;
};


class merMaskedSequence {
public:
  merMaskedSequence(char *fastaName_, char *merylName_, u32bit onlySeqIID_=~u32bitZERO) {
    _numSeq    = 0;
    _seqLen    = 0L;
    _masking   = 0L;
    _repeatID  = 0L;
    _merSize   = 0;

    strcpy(_fastaName, fastaName_);
    strcpy(_merylName, merylName_);

    strcpy(_maskMersName, _merylName);
    strcat(_maskMersName, ".maskMers");

    if (fileExists(_maskMersName))
      loadMasking(onlySeqIID_);
    else
      buildMasking();
  };
  ~merMaskedSequence() {
    delete [] _seqLen;
    for (u32bit i=0; i<_numSeq; i++) {
      delete [] _masking[i];
      delete [] _repeatID[i];
    }
    delete [] _masking;
    delete [] _repeatID;
  };

public:
  u32bit     numSeq(void)                 { return(_numSeq); };
  s32bit     seqLen(u32bit i)             { return(_seqLen[i]); };
  char       masking(u32bit s, u32bit p)  { return(_masking[s][p]); };
  u32bit     repeatID(u32bit s, u32bit p) { return(_repeatID[s][p]); };

  u32bit     merSize(void)                { return(_merSize); };

private:
  void       loadMasking(u32bit onlySeqIID_=~u32bitZERO);  //  Read the masking from the saved file
  void       saveMasking(void);                            //  Write the masking to a file
  void       buildMasking(void);                           //  Read the mers to build the masking

  u32bit    _numSeq;
  s32bit   *_seqLen;   //  signed just for convenience later (positions are signed for same reason)
  char    **_masking;
  u32bit  **_repeatID;

  u32bit    _merSize;

  char      _fastaName[FILENAME_MAX];
  char      _merylName[FILENAME_MAX];
  char      _maskMersName[FILENAME_MAX];
};


void
merMaskedSequence::loadMasking(u32bit onlySeqIID_) {
  FILE  *maskMersFile = fopen(_maskMersName, "r");

  fread(&_numSeq,  sizeof(u32bit), 1, maskMersFile);
  fread(&_merSize, sizeof(u32bit), 1, maskMersFile);

  _seqLen   = new s32bit   [_numSeq];
  _masking  = new char   * [_numSeq];
  _repeatID = new u32bit * [_numSeq];

  fprintf(stderr, u32bitFMT" sequences in '%s'\n", _numSeq, _fastaName);

  fread( _seqLen,   sizeof(u32bit), _numSeq, maskMersFile);

  for (u32bit i=0; i<_numSeq; i++) {
    _masking[i]  = 0L;
    _repeatID[i] = 0L;

    if ((onlySeqIID_ >= _numSeq) || (onlySeqIID_ == i)) {
      fprintf(stderr, "Loading sequence "u32bitFMT" of length "u32bitFMT"\n", i, _seqLen[i]);

      _masking[i]  = new char   [_seqLen[i]];
      _repeatID[i] = new u32bit [_seqLen[i]];

      //memset(_masking[i],  'g', sizeof(char)   * _seqLen[i]);
      //memset(_repeatID[i],   0, sizeof(u32bit) * _seqLen[i]);

      fread(_masking[i],  sizeof(char),   _seqLen[i], maskMersFile);
      fread(_repeatID[i], sizeof(u32bit), _seqLen[i], maskMersFile);
    } else {
      fseek(maskMersFile, sizeof(char)   * _seqLen[i], SEEK_CUR);
      fseek(maskMersFile, sizeof(u32bit) * _seqLen[i], SEEK_CUR);
      _seqLen[i] = 0;
    }
  }

  fclose(maskMersFile);
}


void
merMaskedSequence::saveMasking(void) {
  FILE  *maskMersFile = fopen(_maskMersName, "w");

  fwrite(&_numSeq,  sizeof(u32bit), 1,       maskMersFile);
  fwrite(&_merSize, sizeof(u32bit), 1,       maskMersFile);
  fwrite( _seqLen,  sizeof(u32bit), _numSeq, maskMersFile);

  for (u32bit i=0; i<_numSeq; i++) {
    fwrite(_masking[i],  sizeof(char),   _seqLen[i], maskMersFile);
    fwrite(_repeatID[i], sizeof(u32bit), _seqLen[i], maskMersFile);
  }

  fclose(maskMersFile);
}


void
merMaskedSequence::buildMasking(void) {
  seqStream     *STR = new seqStream(_fastaName);

  _numSeq   = STR->numberOfSequences();

  _seqLen   = new s32bit   [_numSeq];
  _masking  = new char   * [_numSeq];
  _repeatID = new u32bit * [_numSeq];

  _merSize  = 0;

  fprintf(stderr, u32bitFMT" sequences in '%s'\n", _numSeq, _fastaName);

  for (u32bit i=0; i<_numSeq; i++) {
    _seqLen[i]   = STR->lengthOf(i);

    _masking[i]  = new char   [_seqLen[i]];
    _repeatID[i] = new u32bit [_seqLen[i]];

    memset(_masking[i],  'g', sizeof(char)   * _seqLen[i]);
    memset(_repeatID[i],   0, sizeof(u32bit) * _seqLen[i]);
  }

  //  g -> gap in sequence
  //  u -> unique mer
  //  r -> repeat mer
  //
  //  For all the r's we also need to remember the other locations
  //  that repeat is at.  We annotate the map with a repeat id, set if
  //  another copy of the repeat is nearby.

  merylStreamReader *MS = new merylStreamReader(_merylName);
  speedCounter      *CT = new speedCounter(" Masking mers in sequence: %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);

  u32bit             rid = 0;

  _merSize = MS->merSize();

  while (MS->nextMer()) {
    //fprintf(stderr, "mer count="u64bitFMT" pos="u32bitFMT"\n", MS->theCount(), MS->getPosition(0));

    if (MS->theCount() == 1) {
      u32bit p = MS->getPosition(0);
      u32bit s = STR->sequenceNumberOfPosition(p);
      p -= STR->startOf(s);

      _masking[s][p] = 'u';
    } else {
      std::sort(MS->thePositions(), MS->thePositions() + MS->theCount());

      u32bit  lastS = ~u32bitZERO;
      u32bit  lastP = 0;

      rid++;

      for (u32bit i=0; i<MS->theCount(); i++) {
        u32bit p = MS->getPosition(i);
        u32bit s = STR->sequenceNumberOfPosition(p);
        p -= STR->startOf(s);

        //  Always set the masking.
        _masking[s][p] = 'r';

        //  If there is a repeat close by, set the repeat ID.
        if ((s == lastS) && (lastP + 40000 > p)) {
          _repeatID[s][lastP] = rid;
          _repeatID[s][p]     = rid;
        }

        lastS = s;
        lastP = p;
      }
    }

    CT->tick();
  }

  delete CT;

  delete MS;

  delete STR;

  saveMasking();
}


void
computeDensity(merMaskedSequence *S, char *outputPrefix) {
  char    outputName[FILENAME_MAX];
  FILE   *outputFile;
  u32bit  windowSizeMax = 10000;

  for (u32bit s=0; s<S->numSeq(); s++) {

    //  seqLen == 0 iff that sequence is not loaded.
    if (S->seqLen(s) == 0)
      continue;

    sprintf(outputName, "%s.density.seq"u32bitFMTW(02), outputPrefix, s);
    outputFile = fopen(outputName, "w");

    fprintf(stderr, "Starting '%s'\n", outputName);

    fprintf(outputFile, "#window\tunique\trepeat\tgaps\n");

    //  Not the most efficient, but good enough for us right now.

    for (s32bit p=0; p<S->seqLen(s); ) {
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


void
computeMateRescue(merMaskedSequence *S, char *outputPrefix, mateRescueData *lib, u32bit libLen) {
  char    outputName[FILENAME_MAX];
  FILE   *outputFile;
  FILE   *outputData;

  u32bit  closeRepeatsLen = 0;
  u32bit  closeRepeatsMax = 80000;
  s32bit *closeRepeats    = new s32bit [closeRepeatsMax];

  speedCounter *CT = new speedCounter(" Examining repeats: %7.2f Kbases -- %5.2f Kbases/second\r", 1000.0, 0x1ffff, true);

  u32bit  totalDepth = 0;
  for (u32bit l=0; l<libLen; l++)
    totalDepth += lib[l].coverage();

  for (u32bit s=0; s<S->numSeq(); s++) {

    //  seqLen == 0 iff that sequence is not loaded.
    if (S->seqLen(s) == 0)
      continue;

    fprintf(stderr, "Starting sequence "u32bitFMT"\n", s);

    sprintf(outputName, "%s.mateRescue.seq"u32bitFMTW(02)".out", outputPrefix, s);
    outputFile = fopen(outputName, "w");

    sprintf(outputName, "%s.mateRescue.seq"u32bitFMTW(02)".dat", outputPrefix, s);
    outputData = fopen(outputName, "w");

    double  numRR[MAX_COVERAGE] = {0};  //  num repeats rescued (expected) for [] X coverage
    double  numNR[MAX_COVERAGE] = {0};  //  num repeats nonrescuable (expected) for [] X coverage

    u32bit  numRT = 0;  //  num repeats total

    for (s32bit p=0; p<S->seqLen(s); p++) {
      CT->tick();

      double pRtot = 0.0;
      double pFtot = 0.0;

      if ((S->masking(s, p) != 'g') &&
          (S->masking(s, p) != 'u') &&
          (S->masking(s, p) != 'r'))
        fprintf(stderr, "INVALID MASKING - got %d = %c\n", S->masking(s, p), S->masking(s, p));


      if (S->masking(s, p) == 'r') {
        numRT++;

        //  Index over x-coverage in libraries.  MUST BE 1.
        u32bit  ridx = 1;

        for (u32bit l=0; l<libLen; l++) {
          s32bit mean   = lib[l].mean();
          s32bit stddev = lib[l].stddev();
          
          //  Build a list of the same repeat close to this guy.
          closeRepeatsLen = 0;

          if (S->repeatID(s, p) > 0) {
            s32bit pl = (s32bit)floor(p - 3 * stddev);
            s32bit ph = (s32bit)ceil (p + 3 * stddev);

            if (pl < 0)             pl = 0;
            if (ph > S->seqLen(s))  ph = S->seqLen(s);

            for (s32bit pi=pl; pi<ph; pi++)
              if ((S->repeatID(s, pi) == S->repeatID(s, p)) && (pi != p))
                closeRepeats[closeRepeatsLen++] = pi;
          }


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

          //fprintf(stderr, "b1: %d-%d  b2:%d-%d\n", b1l, b1h, b2l, b2h);

          //  probability we can rescue this repeat with this mate pair
          double  pRescue = 0.0;
          double  pFailed = 0.0;

          if (closeRepeatsLen == 0) {
            //  No close repeats, use the fast method.
            for (s32bit b=b1l; b<b1h; b++) {
              if (S->masking(s, b) == 'u')
                pRescue += lib[l].normal(b - p + mean);
            }

            for (s32bit b=b2l; b<b2h; b++) {
              if (S->masking(s, b) == 'u')
                pRescue += lib[l].normal(b - p - mean);
            }
          } else {
            //  Close repeats, gotta be slow.
            for (s32bit b=b1l; b<b1h; b++) {
              if (S->masking(s, b) == 'u') {
                s32bit  mrl = b + mean - 3 * stddev;
                s32bit  mrh = b + mean + 3 * stddev;

                bool    rescuable = true;

                for (u32bit cri=0; rescuable && cri<closeRepeatsLen; cri++)
                  if ((mrl <= closeRepeats[cri]) && (closeRepeats[cri] <= mrh))
                    rescuable = false;

                if (rescuable)
                  pRescue += lib[l].normal(b - p + mean);
                else
                  pFailed += lib[l].normal(b - p + mean);
              }
            }

            for (s32bit b=b2l; b<b2h; b++) {
              if (S->masking(s, b) == 'u') {
                s32bit  mrl = b - mean - 3 * stddev;
                s32bit  mrh = b - mean + 3 * stddev;

                bool    rescuable = true;

                for (u32bit cri=0; rescuable && cri<closeRepeatsLen; cri++)
                  if ((mrl <= closeRepeats[cri]) && (closeRepeats[cri] <= mrh))
                    rescuable = false;

                if (rescuable)
                  pRescue += lib[l].normal(b - p - mean);
                else
                  pFailed += lib[l].normal(b - p - mean);
              }
            }
          }

          //  We're summing over two distributions.
          pRescue /= 2.0;
          pFailed /= 2.0;

          //  Compute probability of rescuing with libraries we've
          //  seen already, and the expected number of repeats
          //  rescued.
          //
          //  We keep track of the probability we rescue this repeat
          //  with additional coverage of libraries.  First 1x of the
          //  first lib, then 2x of the first, etc, etc.
          //
          {
            double  pR = 1.0;
            double  pF = 1.0;
            for (u32bit x=0; x<lib[l].coverage(); x++) {
              //  Makes it here.  pRescue != 1.0
              pR *= (1.0 - pRescue);
              numRR[ridx] += 1 - pR;
              pRtot       += 1 - pR;

              pF *= (1.0 - pFailed);
              numNR[ridx] += 1 - pF;
              pFtot       += 1 - pF;

              ridx++;
            }
          }
        }  // over all libs

        fprintf(outputData, s32bitFMT"\t%f\t%f\n", p, pRtot / totalDepth, pFtot / totalDepth);

      }  //  if masking is r
    }  // over all positions

    fprintf(outputFile, "seqIID\tmerSize\ttRepeat\teRescue\teFailed\tXcov\tmean\tstddev\n");

    for (u32bit x=1, l=0, n=0; l<libLen; x++) {
      fprintf(outputFile, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%.0f\t%.0f\t"u32bitFMT"\t"s32bitFMT"\t"s32bitFMT"\n",
              s, S->merSize(), numRT, numRR[x], numNR[x], x, lib[l].mean(), lib[l].stddev());
      n++;
      if (n >= lib[l].coverage()) {
        l++;
        n = 0;
      }
    }

    fclose(outputFile);
    fclose(outputData);
  }

  delete CT;
}



int
main(int argc, char **argv) {
  char     *merylName    = 0L;
  char     *fastaName    = 0L;
  char     *outputPrefix = 0L;

  u32bit    onlySeqIID = ~u32bitZERO;

  bool      doDensity  = false;
  bool      doRescue   = false;

  mateRescueData  lib[MAX_COVERAGE];
  u32bit          libLen  = 0;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-mers") == 0) {
      merylName = argv[++arg];

    } else if (strcmp(argv[arg], "-seq") == 0) {
      fastaName = argv[++arg];

    } else if (strcmp(argv[arg], "-only") == 0) {
      onlySeqIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-output") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-d") == 0) {
      doDensity = true;

    } else if (strcmp(argv[arg], "-r") == 0) {
      if (atoi(argv[arg+3]) > 0) {
        doRescue = true;
        lib[libLen++].init(atoi(argv[arg+1]), atoi(argv[arg+2]), atoi(argv[arg+3]));
      }
      arg += 3;

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err) || (merylName == 0L) || (fastaName == 0L) || (outputPrefix == 0L)) {
    fprintf(stderr, "usage: %s -mers mers -seq fasta -output prefix [-d] [-r mean stddev coverage]\n", argv[0]);
    exit(1);
  }

  merMaskedSequence *S = new merMaskedSequence(fastaName, merylName, onlySeqIID);

  if (doDensity)
    computeDensity(S, outputPrefix);

  if (doRescue)
    computeMateRescue(S, outputPrefix, lib, libLen);

  return(0);
}
