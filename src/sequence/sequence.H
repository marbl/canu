
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
 *    src/sequence/sequence.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "strings.H"
#include "files.H"
#include "sampledDistribution.H"

#include <vector>
#include <set>
#include <algorithm>

using namespace std;



class summarizeParameters {
public:
  summarizeParameters() {
    genomeSize   = 0;

    limitTo1x    = false;
    breakAtN     = false;

    asSimple     = false;
    asLength     = false;

    asSequences  = true;
    asBases      = false;
  };

  ~summarizeParameters() {
  };


  void      finalize(void) {
  }


  uint64    genomeSize;

  bool      limitTo1x;
  bool      breakAtN;

  bool      asSimple;
  bool      asLength;

  bool      asSequences;
  bool      asBases;
};



class extractParameters {
public:
  extractParameters() {
    asReverse      = false;
    asComplement   = false;
    asUpperCase    = false;
    asLowerCase    = false;
    doMasking      = false;
    maskWithN      = true;
  };

  ~extractParameters() {
  };


  void      finalize(void) {

    //  If no base range specified, output all bases.

    if (baseBgn.size() == 0) {
      baseBgn.push_back(0);
      baseEnd.push_back(UINT64_MAX);
    }

    //  If no sequence range specified, output all sequences.

    if (seqsBgn.size() == 0) {
      seqsBgn.push_back(1);
      seqsEnd.push_back(UINT64_MAX);
    }

    //  If no length restriction, output all lengths.

    if (lensBgn.size() == 0) {
      lensBgn.push_back(0);
      lensEnd.push_back(UINT64_MAX);
    }

    //  Check and adjust the sequence ranges.
    //
    //  To the user, sequences begin at ONE, not ZERO.
    //  To us, sequences begin at zero.

    for (uint32 si=0; si<seqsBgn.size(); si++) {
      if (seqsBgn[si] == 0) {
        fprintf(stderr, "ERROR: sequences begin at 1, not zero.\n");
        exit(1);
      }

      seqsBgn[si] -= 1;
    }

    //  Check and adjust the base ranges.  These are space based.  A quirk in the
    //  command line parsing results in bgn == end if a single number is supplied;
    //  we interpret that to mean 'output the base at space N'.

    for (uint32 bi=0; bi<baseBgn.size(); bi++) {
      if (baseBgn[bi] == baseEnd[bi])
        baseEnd[bi] += 1;
    }
  };


  vector<uint64>  baseBgn;    //  Base ranges to print
  vector<uint64>  baseEnd;    //

  vector<uint64>  seqsBgn;    //  Sequence ranges to print
  vector<uint64>  seqsEnd;    //

  vector<uint64>  lensBgn;    //  Length ranges to print
  vector<uint64>  lensEnd;    //

  bool          asReverse;
  bool          asComplement;
  bool          asUpperCase;
  bool          asLowerCase;

  bool          doMasking;    //  Mask out any base not in baseBgn/baseEnd with 'N'
  bool          maskWithN;    //  Mask with lowercase sequence instead of 'N'
};



class generateParameters {
public:
  generateParameters() {
    minLength             = 0;
    maxLength             = 10000;

    nSeqs                 = 0;
    nBases                = 0;

    useGaussian           = true;
    gMean                 = 0;
    gStdDev               = 0;

    useMirror             = false;
    mirrorInput           = NULL;
    mirrorDistribution    = 0.0;
    mirrorDistributionLen = 0;
    mirrorDistributionMax = 0;
    mirrorDistributionSum = 0;

    aFreq                 = 0.25;
    cFreq                 = 0.25;
    gFreq                 = 0.25;
    tFreq                 = 0.25;
  };

  ~generateParameters() {
  };


  void      finalize(void) {

    //  Check for invalid.  If not set up, just return.

    if ((nSeqs == 0) && (nBases == 0))
      return;

    if (minLength > maxLength)
      fprintf(stderr, "ERROR:  Told to generate sequences with min length larger than max length.\n"), exit(1);

    //  Unlimit any unset limit.

    fprintf(stderr, "Generating up to " F_U64 " sequences and up to " F_U64 " bases.\n", nSeqs, nBases);

    if (nSeqs  == 0)      nSeqs = UINT64_MAX;
    if (nBases == 0)     nBases = UINT64_MAX;

    //  Set Gaussian mean and standard deviation

    gMean   = (minLength + maxLength) / 2.0;
    gStdDev = (maxLength - minLength) / 6.0;

    //  Load lengths to mirror.

    //  Set base frequencies.

    double  fSum = aFreq + cFreq + gFreq + tFreq;

    fprintf(stderr, "Using base frequencies:\n");
    fprintf(stderr, "  A = %7.5f / %7.5f = %7.5f\n", aFreq, fSum, aFreq / fSum);
    fprintf(stderr, "  C = %7.5f / %7.5f = %7.5f\n", cFreq, fSum, cFreq / fSum);
    fprintf(stderr, "  G = %7.5f / %7.5f = %7.5f\n", gFreq, fSum, gFreq / fSum);
    fprintf(stderr, "  T = %7.5f / %7.5f = %7.5f\n", tFreq, fSum, tFreq / fSum);

    aFreq /= fSum;
    cFreq /= fSum;
    gFreq /= fSum;
    tFreq /= fSum;
  };


  uint64    minLength;
  uint64    maxLength;

  uint64    nSeqs;
  uint64    nBases;

  bool      useGaussian;
  double    gMean;
  double    gStdDev;

  bool      useExponential;

  bool      useMirror;
  char     *mirrorInput;
  double    mirrorDistribution;
  uint64    mirrorDistributionLen;
  uint64    mirrorDistributionMax;
  uint64    mirrorDistributionSum;

  double    aFreq;
  double    cFreq;
  double    gFreq;
  double    tFreq;
};



class simulateParameters {
public:
  simulateParameters() {
    genomeSize      = 0;
    circular        = false;

    desiredCoverage = 0.0;
    desiredNumReads = 0;
    desiredNumBases = 0;

    desiredMinLength = 0;
    desiredMaxLength = UINT32_MAX;

    memset(genomeName,  0, FILENAME_MAX+1);
    memset(distribName, 0, FILENAME_MAX+1);
    memset(outputName,  0, FILENAME_MAX+1);
  }

  ~simulateParameters() {
  }


  void      initialize(void) {
  };


  void      finalize(void) {
    dist.loadDistribution(distribName);
  };


  uint64  genomeSize;
  bool    circular;

  double  desiredCoverage;
  uint64  desiredNumReads;
  uint64  desiredNumBases;

  uint32  desiredMinLength;
  uint32  desiredMaxLength;

  sampledDistribution  dist;

  char    genomeName[FILENAME_MAX+1];
  char    distribName[FILENAME_MAX+1];
  char    outputName[FILENAME_MAX+1];
};



class sampleParameters {
public:
  sampleParameters() {
    isPaired        = false;

    desiredCoverage = 0.0;
    genomeSize      = 0;

    desiredNumReads = 0;
    desiredNumBases = 0;

    desiredFraction = 0.0;

    memset(output1, 0, FILENAME_MAX+1);
    memset(output2, 0, FILENAME_MAX+1);
  }

  ~sampleParameters() {
  }


  void      initialize(void) {
  };


  void      finalize(void) {
  }


  bool    isPaired;

  double  desiredCoverage;
  uint64  genomeSize;

  uint64  desiredNumReads;
  uint64  desiredNumBases;

  double  desiredFraction;

  char    output1[FILENAME_MAX+1];
  char    output2[FILENAME_MAX+1];
};



class shiftRegisterParameters {
public:
  shiftRegisterParameters() {
    len   = 0;
    sr[0] = 0;
    sv[0] = 0;
  };
  ~shiftRegisterParameters() {
  };

  void    initialize(void) {
    uint32 srLen = strlen(sr);
    uint32 svLen = strlen(sv);

    if (srLen != svLen)
      fprintf(stderr, "ERROR: srLen %u (%s) svLen %u (%s)\n",
              srLen, sr,
              svLen, sv);
    assert(srLen == svLen);

    len = srLen;
  };

  uint32  len;

  char    sr[17];
  char    sv[17];
};



class mutateParameters {
public:
  mutateParameters() {
    for (uint32 ii=0; ii<256; ii++) {
      for (uint32 jj=0; jj<256; jj++)
        pS[ii][jj] = 0.0;
      pI[ii] = 0.0;
      pD[ii] = 0.0;

      pSubstitute[ii] = 0.0;
    }

    pInsert = 0.0;
  };

  ~mutateParameters() {
  };

public:
  void      setProbabilitySubstititue(double p, char baseToReplace, char replacementBase) {
    fprintf(stderr, "sub %c -> %c %f\n", baseToReplace, replacementBase, p);
    pS[baseToReplace][replacementBase] = p;
  };

  void      setProbabilityInsert(double p, char baseToInsert) {
    fprintf(stderr, "ins %c %f\n", baseToInsert, p);
    pI[baseToInsert] = p;
  };

  void      setProbabilityDelete(double p, char baseToDelete) {
    fprintf(stderr, "del %c %f\n", baseToDelete, p);
    pD[baseToDelete] = p;
  };

  void      finalize(void) {

    for (uint32 ii=0; ii<256; ii++)
      pSubstitute[ii] = 0.0;

    pInsert = 0.0;

    for (uint32 ii=0; ii<256; ii++) {
      for (uint32 jj=0; jj<256; jj++)
        pSubstitute[ii] += pS[ii][jj];

      pInsert           += pI[ii];
      pDelete           += pD[ii];
    }
  };


  //  probability of substituting base [a] with [b]    (base based)
  //  probability of inserting    an A, C, G, T, N     (space based)
  //  probability of deleting     base [a]             (base based)
  //
  //  extensions:
  //  probability of inserting ACGTN between bases a and b.
  //
  double  pS[256][256];
  double  pI[256];
  double  pD[256];

  double  pSubstitute[256];   //  Probability of substituting 'a' for anything.
  double  pInsert;            //  Probability of inserting any base.
  double  pDelete;
};




void doSummarize    (vector<char *> &inputs, summarizeParameters     &sumPar);
void doExtract      (vector<char *> &inputs, extractParameters       &extPar);
void doGenerate     (                        generateParameters      &genPar);
void doSimulate     (vector<char *> &inputs, simulateParameters      &simPar);
void doSample       (vector<char *> &inputs, sampleParameters        &samPar);
void doShiftRegister(                        shiftRegisterParameters &srPar);
void doMutate       (vector<char *> &inputs, mutateParameters        &mutPar);
