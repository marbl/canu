
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    kmer/meryl/maskMers.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2008-MAR-31 to 2014-APR-11
 *      are Copyright 2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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

  void init(int32 mean_, int32 stddev_, uint32 coverage_) {
    _mean      = mean_;
    _stddev    = stddev_;
    _coverage  = coverage_;

    assert(_mean > 3 * _stddev);

    double  a = 1.0 / (_stddev * sqrt(2 * M_PI));
    double  c = 2 * _stddev * _stddev;

    int32  b1l = (int32)floor(-3 * _stddev);
    int32  b1h = (int32)ceil ( 3 * _stddev);

    _normal     = new double [b1h - b1l + 1];
    _normalZero = -b1l;

    for (int32 l=0; l<b1h - b1l + 1; l++)
      _normal[l] = 0.0;

    for (int32 l=b1l; l<b1h; l++)
      _normal[l + _normalZero] = a * exp(- l*l / c);
  };
  ~mateRescueData() {
  };

  int32   mean(void)       { return(_mean); };
  int32   stddev(void)     { return(_stddev); };
  uint32   coverage(void)   { return(_coverage); };

  double   normal(int32 p) { return(_normal[p + _normalZero]); };

private:
  int32  _mean;
  int32  _stddev;
  uint32  _coverage;

  double *_normal;
  int32  _normalZero;
};


class merMaskedSequence {
public:
  merMaskedSequence(char *fastaName_, char *merylName_, uint32 onlySeqIID_=~uint32ZERO) {
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
    for (uint32 i=0; i<_numSeq; i++) {
      delete [] _masking[i];
      delete [] _repeatID[i];
    }
    delete [] _masking;
    delete [] _repeatID;
  };

public:
  uint32     numSeq(void)                 { return(_numSeq); };
  int32     seqLen(uint32 i)             { return(_seqLen[i]); };
  char       masking(uint32 s, uint32 p)  { return(_masking[s][p]); };
  uint32     repeatID(uint32 s, uint32 p) { return(_repeatID[s][p]); };

  uint32     merSize(void)                { return(_merSize); };

private:
  void       loadMasking(uint32 onlySeqIID_=~uint32ZERO);  //  Read the masking from the saved file
  void       saveMasking(void);                            //  Write the masking to a file
  void       buildMasking(void);                           //  Read the mers to build the masking

  uint32    _numSeq;
  int32   *_seqLen;   //  signed just for convenience later (positions are signed for same reason)
  char    **_masking;
  uint32  **_repeatID;

  uint32    _merSize;

  char      _fastaName[FILENAME_MAX];
  char      _merylName[FILENAME_MAX];
  char      _maskMersName[FILENAME_MAX];
};


void
merMaskedSequence::loadMasking(uint32 onlySeqIID_) {
  FILE  *maskMersFile = fopen(_maskMersName, "r");

  fread(&_numSeq,  sizeof(uint32), 1, maskMersFile);
  fread(&_merSize, sizeof(uint32), 1, maskMersFile);

  _seqLen   = new int32   [_numSeq];
  _masking  = new char   * [_numSeq];
  _repeatID = new uint32 * [_numSeq];

  fprintf(stderr, uint32FMT" sequences in '%s'\n", _numSeq, _fastaName);

  fread( _seqLen,   sizeof(uint32), _numSeq, maskMersFile);

  for (uint32 i=0; i<_numSeq; i++) {
    _masking[i]  = 0L;
    _repeatID[i] = 0L;

    if ((onlySeqIID_ >= _numSeq) || (onlySeqIID_ == i)) {
      fprintf(stderr, "Loading sequence "uint32FMT" of length "uint32FMT"\n", i, _seqLen[i]);

      _masking[i]  = new char   [_seqLen[i]];
      _repeatID[i] = new uint32 [_seqLen[i]];

      //memset(_masking[i],  'g', sizeof(char)   * _seqLen[i]);
      //memset(_repeatID[i],   0, sizeof(uint32) * _seqLen[i]);

      fread(_masking[i],  sizeof(char),   _seqLen[i], maskMersFile);
      fread(_repeatID[i], sizeof(uint32), _seqLen[i], maskMersFile);
    } else {
      fseek(maskMersFile, sizeof(char)   * _seqLen[i], SEEK_CUR);
      fseek(maskMersFile, sizeof(uint32) * _seqLen[i], SEEK_CUR);
      _seqLen[i] = 0;
    }
  }

  fclose(maskMersFile);
}


void
merMaskedSequence::saveMasking(void) {
  FILE  *maskMersFile = fopen(_maskMersName, "w");

  fwrite(&_numSeq,  sizeof(uint32), 1,       maskMersFile);
  fwrite(&_merSize, sizeof(uint32), 1,       maskMersFile);
  fwrite( _seqLen,  sizeof(uint32), _numSeq, maskMersFile);

  for (uint32 i=0; i<_numSeq; i++) {
    fwrite(_masking[i],  sizeof(char),   _seqLen[i], maskMersFile);
    fwrite(_repeatID[i], sizeof(uint32), _seqLen[i], maskMersFile);
  }

  fclose(maskMersFile);
}


void
merMaskedSequence::buildMasking(void) {
  seqStream     *STR = new seqStream(_fastaName);

  _numSeq   = STR->numberOfSequences();

  _seqLen   = new int32   [_numSeq];
  _masking  = new char   * [_numSeq];
  _repeatID = new uint32 * [_numSeq];

  _merSize  = 0;

  fprintf(stderr, uint32FMT" sequences in '%s'\n", _numSeq, _fastaName);

  for (uint32 i=0; i<_numSeq; i++) {
    _seqLen[i]   = STR->lengthOf(i);

    _masking[i]  = new char   [_seqLen[i]];
    _repeatID[i] = new uint32 [_seqLen[i]];

    memset(_masking[i],  'g', sizeof(char)   * _seqLen[i]);
    memset(_repeatID[i],   0, sizeof(uint32) * _seqLen[i]);
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

  uint32             rid = 0;

  _merSize = MS->merSize();

  while (MS->nextMer()) {
    //fprintf(stderr, "mer count="uint64FMT" pos="uint32FMT"\n", MS->theCount(), MS->getPosition(0));

    if (MS->theCount() == 1) {
      uint32 p = MS->getPosition(0);
      uint32 s = STR->sequenceNumberOfPosition(p);
      p -= STR->startOf(s);

      _masking[s][p] = 'u';
    } else {
      std::sort(MS->thePositions(), MS->thePositions() + MS->theCount());

      uint32  lastS = ~uint32ZERO;
      uint32  lastP = 0;

      rid++;

      for (uint32 i=0; i<MS->theCount(); i++) {
        uint32 p = MS->getPosition(i);
        uint32 s = STR->sequenceNumberOfPosition(p);
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
  uint32  windowSizeMax = 10000;

  for (uint32 s=0; s<S->numSeq(); s++) {

    //  seqLen == 0 iff that sequence is not loaded.
    if (S->seqLen(s) == 0)
      continue;

    sprintf(outputName, "%s.density.seq"uint32FMTW(02), outputPrefix, s);
    outputFile = fopen(outputName, "w");

    fprintf(stderr, "Starting '%s'\n", outputName);

    fprintf(outputFile, "#window\tunique\trepeat\tgaps\n");

    //  Not the most efficient, but good enough for us right now.

    for (int32 p=0; p<S->seqLen(s); ) {
      uint32  windowSize    = 0;
      uint32  uniqueSum     = 0;
      uint32  repeatSum     = 0;
      uint32  gapSum        = 0;

      while ((windowSize < windowSizeMax) &&
             (p < S->seqLen(s))) {
        char m = S->masking(s, p);

        if (m == 'u')  uniqueSum++;
        if (m == 'g')  gapSum++;
        if (m == 'r')  repeatSum++;

        windowSize++;
        p++;
      }

      fprintf(outputFile, uint32FMT"\t%f\t%f\t%f\n",
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
computeMateRescue(merMaskedSequence *S, char *outputPrefix, mateRescueData *lib, uint32 libLen) {
  char    outputName[FILENAME_MAX];
  FILE   *outputFile;
  FILE   *outputData;

  uint32  closeRepeatsLen = 0;
  uint32  closeRepeatsMax = 80000;
  int32 *closeRepeats    = new int32 [closeRepeatsMax];

  speedCounter *CT = new speedCounter(" Examining repeats: %7.2f Kbases -- %5.2f Kbases/second\r", 1000.0, 0x1ffff, true);

  uint32  totalDepth = 0;
  for (uint32 l=0; l<libLen; l++)
    totalDepth += lib[l].coverage();

  for (uint32 s=0; s<S->numSeq(); s++) {

    //  seqLen == 0 iff that sequence is not loaded.
    if (S->seqLen(s) == 0)
      continue;

    fprintf(stderr, "Starting sequence "uint32FMT"\n", s);

    sprintf(outputName, "%s.mateRescue.seq"uint32FMTW(02)".out", outputPrefix, s);
    outputFile = fopen(outputName, "w");

    sprintf(outputName, "%s.mateRescue.seq"uint32FMTW(02)".dat", outputPrefix, s);
    outputData = fopen(outputName, "w");

    double  numRR[MAX_COVERAGE] = {0};  //  num repeats rescued (expected) for [] X coverage
    double  numNR[MAX_COVERAGE] = {0};  //  num repeats nonrescuable (expected) for [] X coverage

    uint32  numRT = 0;  //  num repeats total

    for (int32 p=0; p<S->seqLen(s); p++) {
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
        uint32  ridx = 1;

        for (uint32 l=0; l<libLen; l++) {
          int32 mean   = lib[l].mean();
          int32 stddev = lib[l].stddev();

          //  Build a list of the same repeat close to this guy.
          closeRepeatsLen = 0;

          if (S->repeatID(s, p) > 0) {
            int32 pl = (int32)floor(p - 3 * stddev);
            int32 ph = (int32)ceil (p + 3 * stddev);

            if (pl < 0)             pl = 0;
            if (ph > S->seqLen(s))  ph = S->seqLen(s);

            for (int32 pi=pl; pi<ph; pi++)
              if ((S->repeatID(s, pi) == S->repeatID(s, p)) && (pi != p))
                closeRepeats[closeRepeatsLen++] = pi;
          }


          int32 b1l = (int32)floor(p - mean - 3 * stddev);
          int32 b1h = (int32)ceil (p - mean + 3 * stddev);

          int32 b2l = (int32)floor(p + mean - 3 * stddev);
          int32 b2h = (int32)ceil (p + mean + 3 * stddev);

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
            for (int32 b=b1l; b<b1h; b++) {
              if (S->masking(s, b) == 'u')
                pRescue += lib[l].normal(b - p + mean);
            }

            for (int32 b=b2l; b<b2h; b++) {
              if (S->masking(s, b) == 'u')
                pRescue += lib[l].normal(b - p - mean);
            }
          } else {
            //  Close repeats, gotta be slow.
            for (int32 b=b1l; b<b1h; b++) {
              if (S->masking(s, b) == 'u') {
                int32  mrl = b + mean - 3 * stddev;
                int32  mrh = b + mean + 3 * stddev;

                bool    rescuable = true;

                for (uint32 cri=0; rescuable && cri<closeRepeatsLen; cri++)
                  if ((mrl <= closeRepeats[cri]) && (closeRepeats[cri] <= mrh))
                    rescuable = false;

                if (rescuable)
                  pRescue += lib[l].normal(b - p + mean);
                else
                  pFailed += lib[l].normal(b - p + mean);
              }
            }

            for (int32 b=b2l; b<b2h; b++) {
              if (S->masking(s, b) == 'u') {
                int32  mrl = b - mean - 3 * stddev;
                int32  mrh = b - mean + 3 * stddev;

                bool    rescuable = true;

                for (uint32 cri=0; rescuable && cri<closeRepeatsLen; cri++)
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
            for (uint32 x=0; x<lib[l].coverage(); x++) {
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

        fprintf(outputData, int32FMT"\t%f\t%f\n", p, pRtot / totalDepth, pFtot / totalDepth);

      }  //  if masking is r
    }  // over all positions

    fprintf(outputFile, "seqIID\tmerSize\ttRepeat\teRescue\teFailed\tXcov\tmean\tstddev\n");

    for (uint32 x=1, l=0, n=0; l<libLen; x++) {
      fprintf(outputFile, uint32FMT"\t"uint32FMT"\t"uint32FMT"\t%.0f\t%.0f\t"uint32FMT"\t"int32FMT"\t"int32FMT"\n",
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

  uint32    onlySeqIID = ~uint32ZERO;

  bool      doDensity  = false;
  bool      doRescue   = false;

  mateRescueData  lib[MAX_COVERAGE];
  uint32          libLen  = 0;

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
