
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sequence.H"

#include "mt19937ar.H"
#include "files.H"
#include "strings.H"

#include <vector>
#include <set>
#include <algorithm>

using namespace std;


enum opMode {
  modeSummarize,
  modeExtract,
  modeGenerate,
  modeSimulate,
  modeSample,
  modeShift,
  modeUnset
};


class summarizeParameters {
public:
  summarizeParameters() {
    genomeSize   = 0;

    limitTo1x    = false;
    breakAtN     = false;

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





bool
doSummarize_loadSequence(dnaSeqFile  *sf,
                         bool         asSequences,
                         char       *&name,   uint32    &nameMax,
                         char       *&seq,
                         uint8      *&qlt,    uint64    &seqMax,
                         uint64      &seqLen) {

  if (asSequences)
    return(sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen));

  //  Otherwise, piece it together from multiple calls to get bases.
  //  loadBases() returns true if bases were loaded, and sets endOfSeq
  //  if the block returned is to the end of a sequence.

  uint64   bufferMax = 23;
  uint64   bufferLen = 0;
  char    *buffer    = new char [bufferMax];
  bool     endOfSeq  = false;

  resizeArray(name, 0, nameMax, (uint32)1024);
  resizeArrayPair(seq, qlt, 0, seqMax, seqLen+1, resizeArray_doNothing);

  name[0] = 0;
  seq[0]  = 0;
  qlt[0]  = 0;

  seqLen = 0;

  while (sf->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
    if (seqLen + bufferLen >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 2 * (seqLen + bufferLen + 1));

    assert(seqLen + bufferLen + 1 < seqMax);

    memcpy(seq + seqLen, buffer, sizeof(char) * bufferLen);
    seqLen += bufferLen;

    seq[seqLen] = 0;

    if (endOfSeq)
      break;
  }

  //  We get here in two ways:
  //    loadBases() immediately hit EOF, it will return false and set endOfSeq = false.
  //    Otherwise, it found a sequence and endOfSeq (here) _must_ be true.
  //  So, we can just return endOfSeq to mean 'there might be another sequence to load'.

  delete [] buffer;

  return(endOfSeq);
}



void
doSummarize(vector<char *>       &inputs,
            summarizeParameters  &sumPar) {

  vector<uint64>  lengths;

  uint64          nSeqs  = 0;
  uint64          nBases = 0;

  uint32          mer = 0;

  uint64          mn[4]     = {0};
  uint64          dn[4*4]   = {0};
  uint64          tn[4*4*4] = {0};

  double          nmn = 0;
  double          ndn = 0;
  double          ntn = 0;

  uint32          nameMax = 0;
  char           *name    = NULL;
  uint64          seqMax  = 0;
  char           *seq     = NULL;
  uint8          *qlt     = NULL;
  uint64          seqLen  = 0;

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf = new dnaSeqFile(inputs[ff]);

    //  If sequence,
    //    Count mono-, di- and tri-nucleotides.
    //    Count number of mono-, di- and tri-nucleotides.
    //    Count number of sequences and total bases.
    //    Save the lengths of sequences.

    while (doSummarize_loadSequence(sf, sumPar.asSequences, name, nameMax, seq, qlt, seqMax, seqLen) == true) {
      uint64  pos = 0;
      uint64  bgn = 0;

      if (pos < seqLen) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
      }

      if (pos < seqLen) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
        dn[mer & 0x0f]++;
      }

      while (pos < seqLen) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
        dn[mer & 0x0f]++;
        tn[mer & 0x3f]++;
      }

      nmn +=                    (seqLen-0);
      ndn += (seqLen < 2) ? 0 : (seqLen-1);
      ntn += (seqLen < 3) ? 0 : (seqLen-2);

      //  If we're NOT splitting on N, add one sequence of the given length.

      if (sumPar.breakAtN == false) {
        nSeqs  += 1;
        nBases += seqLen;

        lengths.push_back(seqLen);
        continue;
      }

      //  But if we ARE splitting on N, add multiple sequences.

      pos = 0;
      bgn = 0;

      while (pos < seqLen) {

        //  Skip any N's.
        while ((pos < seqLen) && ((seq[pos] == 'n') ||
                                  (seq[pos] == 'N')))
          pos++;

        //  Remember our start position.
        bgn = pos;

        //  Move ahead until the end of sequence or an N.
        while ((pos < seqLen) && ((seq[pos] != 'n') &&
                                  (seq[pos] != 'N')))
          pos++;

        //  If a sequence, increment stuff.
        if (pos - bgn > 0) {
          nSeqs  += 1;
          nBases += pos - bgn;

          lengths.push_back(pos - bgn);
        }
      }
    }

    //  All done!

    delete sf;
  }

  delete [] name;
  delete [] seq;
  delete [] qlt;

  //  Finalize.

  sort(lengths.begin(), lengths.end(), greater<uint64>());

  if (sumPar.genomeSize == 0)
    sumPar.genomeSize = nBases;

  uint64   lSum  = 0;                                 //  Sum of the lengths we've encountered so far

  uint32   nStep = 10;                                //  Step of each N report.
  uint32   nVal  = nStep;                             //  Index of the threshold we're next printing.
  uint64   nThr  = sumPar.genomeSize * nVal / 100;    //  Threshold lenth; if sum is bigger, emit and move to the next threshold

  ////////////////////////////////////////

  uint32   nCols    = 63;   //  Magic number to make the histogram the same width as the trinucleotide list
  uint32   nRows    = 0;    //  Height of the histogram; dynamically set.
  uint32   nRowsMin = 50;   //  Nothing really magic, just fits on the screen.

  //  Count the number of lines we expect to get in the NG table.

  for (uint32 ii=0; ii<nSeqs; ii++) {
    lSum += lengths[ii];

    while (lSum >= nThr) {
      nRows++;

      if      (nVal  <   200)  nVal += nStep;
      else if (nVal  <  2000)  nVal += nStep * 10;
      else if (nVal  < 20000)  nVal += nStep * 100;
      else                     nVal += nStep * 1000;

      nThr  = sumPar.genomeSize * nVal / 100;
    }
  }

  uint64   minLength = lengths[nSeqs-1];
  uint64   maxLength = lengths[0];

  if (nRows < nRowsMin)                    //  If there are too few lines, make it some minimum value,
    nRows = nRowsMin;                      //  otherwise, make the length histogram plot the same size.

  double   bucketSized = (double)(maxLength - minLength) / nRows;
  uint32   bucketSize  = (uint32)ceil(bucketSized);

  if (bucketSize == 0)
    bucketSize = 1;

  nRows = (maxLength - minLength) / bucketSize;

  if (nRows > nRowsMin)
    nRows = nRowsMin;

  lSum  = 0;                               //  Reset for actually generating the length histogram.
  nStep = 10;
  nVal  = nStep;
  nThr  = sumPar.genomeSize * nVal / 100;

  //  Generate the length histogram.

  uint32  *nSeqPerLen = new uint32 [nRows + 1];

  for (uint32 rr=0; rr<nRows+1; rr++)                   //  Clear the histogram.
    nSeqPerLen[rr] = 0;

  for (uint32 ii=0; ii<nSeqs; ii++) {                   //  Count number of sequences per size range.
    uint32 r = (lengths[ii] - minLength) / bucketSize;

    nSeqPerLen[r]++;
  }

  uint64  maxCount = 1;                                 //  Avoids divide-by-zero, even if zero can never actually occur.

  for (uint32 rr=0; rr<nRows+1; rr++)                   //  Find the maximum number of sequences in a size range.
    if (maxCount < nSeqPerLen[rr])
      maxCount = nSeqPerLen[rr];

  char **histPlot = new char * [nRows + 1];

  for (uint32 rr=0; rr<nRows+1; rr++)                   //  28 is a magic number based on the histPlot[] format string below.
    histPlot[rr] = new char [28 + nCols + 1];           //  28 = 9 + 1 + 9 + 1 + 7 + 1

  for (uint32 rr=0; rr<nRows+1; rr++) {                 //  Generate histogram in text.
    uint32  nn = (uint32)ceil(nSeqPerLen[rr] * nCols / (double)maxCount);
    uint64  lo = (rr+0) * bucketSize + minLength;
    uint64  hi = (rr+1) * bucketSize + minLength - 1;

    if (lo == hi)
      sprintf(histPlot[rr], "%9" F_U64P "           %7" F_U32P "|",
              lo, nSeqPerLen[rr]);
    else
      sprintf(histPlot[rr], "%9" F_U64P "-%-9" F_U64P " %7" F_U32P "|",
              lo, hi, nSeqPerLen[rr]);

    for (uint32 cc=0; cc<nn; cc++)
      histPlot[rr][28 + cc] = '-';

    histPlot[rr][28 + nn]   = 0;
  }

  //  Output N table, with length histogram appended at the end of each line.

  uint32  hp = 0;

  fprintf(stdout, "\n");
  fprintf(stdout, "G=%-12" F_U64P "                     sum of  ||               length     num\n", sumPar.genomeSize);
  fprintf(stdout,   "NG         length     index       lengths  ||                range    seqs\n");
  fprintf(stdout,   "----- ------------ --------- ------------  ||  ------------------- -------\n");

  //  Write lines if we're showing all data, or if we're below 1x coverage.

  for (uint32 ii=0; ii<nSeqs; ii++) {
    lSum += lengths[ii];

    while (lSum >= nThr) {
      if ((sumPar.limitTo1x == false) ||
          (nVal <= 100)) {
        if (hp <= nRows)
          fprintf(stdout, "%05"    F_U32P " %12" F_U64P " %9" F_U32P " %12" F_U64P "  ||  %s\n",
                  nVal, lengths[ii], ii, lSum,
                  histPlot[hp++]);
        else
          fprintf(stdout, "%05"    F_U32P " %12" F_U64P " %9" F_U32P " %12" F_U64P "  ||\n",
                  nVal, lengths[ii], ii, lSum);
      }

      if      (nVal <   200)   nVal += nStep;
      else if (nVal <  2000)   nVal += nStep * 10;
      else if (nVal < 20000)   nVal += nStep * 100;
      else                     nVal += nStep * 1000;

      nThr  = sumPar.genomeSize * nVal / 100;
    }
  }

  //  If we're displaying exactly 1x, write empty lines to get to there.

  if (sumPar.limitTo1x == true) {
    while (nVal <= 100) {
      if (hp <= nRows)
        fprintf(stdout, "%05"    F_U32P " %12s %9s %12s  ||  %s\n",
                nVal, "-", "-", "-",
                histPlot[hp++]);
      else
        fprintf(stdout, "%05"    F_U32P " %12s %9s %12s  ||\n",
                nVal, "-", "-", "-");

      nVal += nStep;
    }
  }

  //  Now the final summary line.

  if (sumPar.genomeSize == 0)
    fprintf(stdout, "%07.3fx           %9" F_U64P " %12" F_U64P "  ||  %s\n", 0.0, nSeqs, lSum, histPlot[hp++]);   //  Occurs if only empty sequences in the input!
  else if (hp <= nRows)
    fprintf(stdout, "%07.3fx           %9" F_U64P " %12" F_U64P "  ||  %s\n", (double)lSum / sumPar.genomeSize, nSeqs, lSum, histPlot[hp++]);
  else
    fprintf(stdout, "%07.3fx           %9" F_U64P " %12" F_U64P "  ||\n",     (double)lSum / sumPar.genomeSize, nSeqs, lSum);

  //fprintf(stdout,    "                                           ||  %s\n", histPlot[hp++]);
  //fprintf(stdout,   "                 genome-size %12" F_U64P "  ||  %s\n", sumPar.genomeSize, histPlot[hp++]);

  while (hp <= nRows)
    fprintf(stdout, "                                           ||  %s\n", histPlot[hp++]);

  fprintf(stdout, "\n");

  //  Report.  To get alphabetic ordering correct, hardcoded.

#define FMT "%12" F_U64P " %6.4f"
#define GC "%05.02f%%"

  if (nmn == 0)  nmn = 1;   //  Avoid divide by zero.
  if (ndn == 0)  ndn = 1;
  if (ntn == 0)  ntn = 1;

  double gc = 100.0 * (mn[0x01] + mn[0x03]) / (mn[0x00] + mn[0x01] + mn[0x03] + mn[0x02]);

  fprintf(stdout, "--------------------- --------------------- ----------------------------------------------------------------------------------------------\n");
  fprintf(stdout, "       mononucleotide          dinucleotide                                                                                  trinucleotide\n");
  fprintf(stdout, "--------------------- --------------------- ----------------------------------------------------------------------------------------------\n");
  fprintf(stdout, ""             FMT " A" FMT " AA" FMT " AAA " FMT " AAC " FMT " AAG " FMT " AAT\n", mn[0x00], mn[0x00] / nmn, dn[0x00], dn[0x00] / ndn, tn[0x00], tn[0x00] / ntn, tn[0x01], tn[0x01] / ntn, tn[0x03], tn[0x03] / ntn, tn[0x02], tn[0x02] / ntn);
  fprintf(stdout, ""             FMT " C" FMT " AC" FMT " ACA " FMT " ACC " FMT " ACG " FMT " ACT\n", mn[0x01], mn[0x01] / nmn, dn[0x01], dn[0x01] / ndn, tn[0x04], tn[0x04] / ntn, tn[0x05], tn[0x05] / ntn, tn[0x07], tn[0x07] / ntn, tn[0x06], tn[0x06] / ntn);
  fprintf(stdout, ""             FMT " G" FMT " AG" FMT " AGA " FMT " AGC " FMT " AGG " FMT " AGT\n", mn[0x03], mn[0x03] / nmn, dn[0x03], dn[0x03] / ndn, tn[0x0c], tn[0x0c] / ntn, tn[0x0d], tn[0x0d] / ntn, tn[0x0f], tn[0x0f] / ntn, tn[0x0e], tn[0x0e] / ntn);
  fprintf(stdout, ""             FMT " T" FMT " AT" FMT " ATA " FMT " ATC " FMT " ATG " FMT " ATT\n", mn[0x02], mn[0x02] / nmn, dn[0x02], dn[0x02] / ndn, tn[0x08], tn[0x08] / ntn, tn[0x09], tn[0x09] / ntn, tn[0x0b], tn[0x0b] / ntn, tn[0x0a], tn[0x0a] / ntn);
  fprintf(stdout, "                     " FMT " CA" FMT " CAA " FMT " CAC " FMT " CAG " FMT " CAT\n",                           dn[0x04], dn[0x04] / ndn, tn[0x10], tn[0x10] / ntn, tn[0x11], tn[0x11] / ntn, tn[0x13], tn[0x13] / ntn, tn[0x12], tn[0x12] / ntn);
  fprintf(stdout, "      --GC--  --AT-- " FMT " CC" FMT " CCA " FMT " CCC " FMT " CCG " FMT " CCT\n",                           dn[0x05], dn[0x05] / ndn, tn[0x14], tn[0x14] / ntn, tn[0x15], tn[0x15] / ntn, tn[0x17], tn[0x17] / ntn, tn[0x16], tn[0x16] / ntn);
  fprintf(stdout, "      " GC "  " GC " " FMT " CG" FMT " CGA " FMT " CGC " FMT " CGG " FMT " CGT\n", gc, 100 - gc,             dn[0x07], dn[0x07] / ndn, tn[0x1c], tn[0x1c] / ntn, tn[0x1d], tn[0x1d] / ntn, tn[0x1f], tn[0x1f] / ntn, tn[0x1e], tn[0x1e] / ntn);
  fprintf(stdout, "                     " FMT " CT" FMT " CTA " FMT " CTC " FMT " CTG " FMT " CTT\n",                           dn[0x06], dn[0x06] / ndn, tn[0x18], tn[0x18] / ntn, tn[0x19], tn[0x19] / ntn, tn[0x1b], tn[0x1b] / ntn, tn[0x1a], tn[0x1a] / ntn);
  fprintf(stdout, "                     " FMT " GA" FMT " GAA " FMT " GAC " FMT " GAG " FMT " GAT\n",                           dn[0x0c], dn[0x0c] / ndn, tn[0x30], tn[0x30] / ntn, tn[0x31], tn[0x31] / ntn, tn[0x33], tn[0x33] / ntn, tn[0x32], tn[0x32] / ntn);
  fprintf(stdout, "                     " FMT " GC" FMT " GCA " FMT " GCC " FMT " GCG " FMT " GCT\n",                           dn[0x0d], dn[0x0d] / ndn, tn[0x34], tn[0x34] / ntn, tn[0x35], tn[0x35] / ntn, tn[0x37], tn[0x37] / ntn, tn[0x36], tn[0x36] / ntn);
  fprintf(stdout, "                     " FMT " GG" FMT " GGA " FMT " GGC " FMT " GGG " FMT " GGT\n",                           dn[0x0f], dn[0x0f] / ndn, tn[0x3c], tn[0x3c] / ntn, tn[0x3d], tn[0x3d] / ntn, tn[0x3f], tn[0x3f] / ntn, tn[0x3e], tn[0x3e] / ntn);
  fprintf(stdout, "                     " FMT " GT" FMT " GTA " FMT " GTC " FMT " GTG " FMT " GTT\n",                           dn[0x0e], dn[0x0e] / ndn, tn[0x38], tn[0x38] / ntn, tn[0x39], tn[0x39] / ntn, tn[0x3b], tn[0x3b] / ntn, tn[0x3a], tn[0x3a] / ntn);
  fprintf(stdout, "                     " FMT " TA" FMT " TAA " FMT " TAC " FMT " TAG " FMT " TAT\n",                           dn[0x08], dn[0x08] / ndn, tn[0x20], tn[0x20] / ntn, tn[0x21], tn[0x21] / ntn, tn[0x23], tn[0x23] / ntn, tn[0x22], tn[0x22] / ntn);
  fprintf(stdout, "                     " FMT " TC" FMT " TCA " FMT " TCC " FMT " TCG " FMT " TCT\n",                           dn[0x09], dn[0x09] / ndn, tn[0x24], tn[0x24] / ntn, tn[0x25], tn[0x25] / ntn, tn[0x27], tn[0x27] / ntn, tn[0x26], tn[0x26] / ntn);
  fprintf(stdout, "                     " FMT " TG" FMT " TGA " FMT " TGC " FMT " TGG " FMT " TGT\n",                           dn[0x0b], dn[0x0b] / ndn, tn[0x2c], tn[0x2c] / ntn, tn[0x2d], tn[0x2d] / ntn, tn[0x2f], tn[0x2f] / ntn, tn[0x2e], tn[0x2e] / ntn);
  fprintf(stdout, "                     " FMT " TT" FMT " TTA " FMT " TTC " FMT " TTG " FMT " TTT\n",                           dn[0x0a], dn[0x0a] / ndn, tn[0x28], tn[0x28] / ntn, tn[0x29], tn[0x29] / ntn, tn[0x2b], tn[0x2b] / ntn, tn[0x2a], tn[0x2a] / ntn);
  fprintf(stdout, "\n");

  //  Cleanup.

  for (uint32 rr=0; rr<nRows+1; rr++)
    delete [] histPlot[rr];

  delete [] histPlot;
  delete [] nSeqPerLen;

}



void
doExtract(vector<char *>    &inputs,
          extractParameters &extPar) {

  char            C[256] = {0};
  char            U[256] = {0};
  char            L[256] = {0};

  uint32          nameMax = 0;
  char           *name    = NULL;
  uint64          seqMax  = 0;
  char           *seq     = NULL;
  uint8          *qlt     = NULL;
  uint64          seqLen  = 0;

  uint64  outputStringLen = 0;
  uint64  outputStringMax = 0;
  char   *outputString    = NULL;

  //  Initialize complement. toUpper and toLower arrays.

  C['a'] = 't';  U['a'] = 'A';  L['a'] = 'a';
  C['c'] = 'g';  U['c'] = 'C';  L['c'] = 'c';
  C['g'] = 'c';  U['g'] = 'G';  L['g'] = 'g';
  C['t'] = 'a';  U['t'] = 'T';  L['t'] = 't';

  C['A'] = 'T';  U['A'] = 'A';  L['A'] = 'a';
  C['C'] = 'G';  U['C'] = 'C';  L['C'] = 'c';
  C['G'] = 'C';  U['G'] = 'G';  L['G'] = 'g';
  C['T'] = 'A';  U['T'] = 'T';  L['T'] = 't';



  for (uint32 fi=0; fi<inputs.size(); fi++) {
    dnaSeqFile  *sf   = new dnaSeqFile(inputs[fi], true);

    //  Allocate a string big enough to hold the largest output.
    //
    //  Later, maybe, we can analyze the bases to output and make this exactly the correct size.

    uint64  maxStringLength = 0;

    for (uint32 ss=0; ss<sf->numberOfSequences(); ss++)
      maxStringLength = max(maxStringLength, sf->sequenceLength(ss));

    resizeArray(outputString, 0, outputStringMax, maxStringLength + 1);

    //fprintf(stderr, "seqs - length %u first %u %u\n", extPar.seqsBgn.size(), extPar.seqsBgn[0], extPar.seqsEnd[0]);

    for (uint32 si=0; si<extPar.seqsBgn.size(); si++) {
      uint64  sbgn = extPar.seqsBgn[si];
      uint64  send = extPar.seqsEnd[si];

      sbgn = min(sbgn, sf->numberOfSequences());
      send = min(send, sf->numberOfSequences());

      //fprintf(stderr, "sbgn %u send %u\n", sbgn, send);

      for (uint32 ss=sbgn; ss<send; ss++) {
        uint64  seqLen = sf->sequenceLength(ss);

        //fprintf(stderr, "lens - length %u first %u %u\n", extPar.lensBgn.size(), extPar.lensBgn[0], extPar.lensEnd[0]);

        for (uint32 li=0; li<extPar.lensBgn.size(); li++) {
          uint64  lmin = extPar.lensBgn[li];
          uint64  lmax = extPar.lensEnd[li];

          if ((seqLen < lmin) ||
              (lmax < seqLen))
            seqLen = UINT64_MAX;
        }

        if (seqLen == UINT64_MAX)
          continue;

        if (sf->findSequence(ss) == false) {
          //fprintf(stderr, "Failed to find sequence #%u in file '%s'\n", ss, inputs[fi]);
          continue;
        }

        if (sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen) == false) {
          //fprintf(stderr, "Failed to load sequence #%u in file '%s'\n", ss, inputs[fi]);
          continue;
        }

        //fprintf(stderr, "base - length %u first %u %u\n", extPar.baseBgn.size(), extPar.baseBgn[0], extPar.baseEnd[0]);

        outputStringLen = 0;

        for (uint32 bi=0; bi<extPar.baseBgn.size(); bi++) {
          uint64  bbgn = extPar.baseBgn[bi];
          uint64  bend = extPar.baseEnd[bi];

          bbgn = min(bbgn, seqLen);
          bend = min(bend, seqLen);

          //fprintf(stderr, "base - seqLen %u base[%u] %u %u limited %u %u\n", seqLen, bi, extPar.baseBgn[bi], extPar.baseEnd[bi], bbgn, bend);

          if (bbgn == bend)
            continue;

          memcpy(outputString + outputStringLen, seq + bbgn, bend - bbgn);

          outputStringLen += bend - bbgn;
        }

        outputString[outputStringLen] = 0;

        if (extPar.asReverse)
          reverse(outputString, outputString + outputStringLen);

        if (extPar.asComplement)
          for (uint32 ii=0; ii<outputStringLen; ii++)
            outputString[ii] = C[outputString[ii]];

        if (extPar.asUpperCase)
          for (uint32 ii=0; ii<outputStringLen; ii++)
            outputString[ii] = U[outputString[ii]];

        if (extPar.asLowerCase)
          for (uint32 ii=0; ii<outputStringLen; ii++)
            outputString[ii] = L[outputString[ii]];

        fprintf(stdout, ">%s\n%s\n", name, outputString);
      }
    }

    //  Done with this file.  Get rid of it.

    delete sf;
  }

  //  Cleanup.

  delete [] name;
  delete [] seq;
  delete [] qlt;
  delete [] outputString;
}







void
doGenerate(generateParameters &genPar) {
  mtRandom   MT;

  uint64  nSeqs  = 0;
  uint64  nBases = 0;

  uint64  seqLen = 0;
  uint64  seqMax = 65536;
  char   *seq    = new char  [seqMax + 1];
  uint8  *qlt    = new uint8 [seqMax + 1];

  double Athresh = genPar.aFreq;
  double Cthresh = genPar.aFreq + genPar.cFreq;
  double Gthresh = genPar.aFreq + genPar.cFreq + genPar.gFreq;
  double Tthresh = genPar.aFreq + genPar.cFreq + genPar.gFreq + genPar.tFreq;

  while ((nSeqs  < genPar.nSeqs) &&
         (nBases < genPar.nBases)) {
    double   len = MT.mtRandomGaussian(genPar.gMean, genPar.gStdDev);

    while (len < -0.5)
      len = MT.mtRandomGaussian(genPar.gMean, genPar.gStdDev);

    seqLen = (uint64)round(len);

    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, 0, seqMax, seqLen+1, resizeArray_doNothing);

    for (uint64 ii=0; ii<seqLen; ii++) {
      double  bp = MT.mtRandomRealOpen();

      if        (bp < Athresh) {
        seq[ii] = 'A';
        qlt[ii] = 0;

      } else if (bp < Cthresh) {
        seq[ii] = 'C';
        qlt[ii] = 0;

      } else if (bp < Gthresh) {
        seq[ii] = 'G';
        qlt[ii] = 0;

      } else {
        seq[ii] = 'T';
        qlt[ii] = 0;
      }
    }

    seq[seqLen] = 0;
    qlt[seqLen] = 0;

    AS_UTL_writeFastA(stdout,
                      seq, seqLen, 0,
                      ">random%08 " F_U64P "\n", nSeqs);

    nSeqs  += 1;
    nBases += seqLen;
  }

  delete [] seq;
  delete [] qlt;
}



void
doSimulate(vector<char *> &inputs) {
}




class seqEntry {
public:
  seqEntry(mtRandom &MT, uint64 pos_) {
    pos = pos_;
    rnd = MT.mtRandomRealOpen();
  };

  uint64   pos;    //  Position in the seqLengths vector.
  double   rnd;    //  A random number.
};


bool seqOrderRandom(const seqEntry &a, const seqEntry &b) {
  return(a.rnd < b.rnd);
}

bool seqOrderNormal(const seqEntry &a, const seqEntry &b) {
  return(a.pos < b.pos);
}



void
doSample_sample(sampleParameters &samPar,
                uint64            numSeqsTotal,
                uint64            numBasesTotal,
                vector<uint64>   &seqLengths,
                vector<seqEntry> &seqOrder) {

  //  Randomize the sequences.

  sort(seqOrder.begin(), seqOrder.end(), seqOrderRandom);

  //  Do some math to figure out what sequences to report.

  if (samPar.desiredCoverage > 0.0) {
    samPar.desiredNumBases = (uint64)ceil(samPar.desiredCoverage * samPar.genomeSize);
  }

  if (samPar.desiredNumReads > 0) {
    fprintf(stderr, "Emitting " F_U64 " reads.\n",
            samPar.desiredNumReads);

    for (uint64 ii=0; ii<numSeqsTotal; ii++)
      if (ii < samPar.desiredNumReads)
        ;
      else
        seqLengths[seqOrder[ii].pos] = 0;
  }

  if (samPar.desiredNumBases > 0) {
    if (samPar.desiredCoverage > 0)
      fprintf(stderr, "Emitting %.3fx coverage; " F_U64 " bases.\n",
              samPar.desiredCoverage,
              samPar.desiredNumBases);
    else
      fprintf(stderr, "Emitting " F_U64 " bases.\n",
              samPar.desiredNumBases);

    for (uint64 nbe=0, ii=0; ii<numSeqsTotal; ii++) {
      if (nbe < samPar.desiredNumBases)
        nbe += seqLengths[seqOrder[ii].pos];
      else
        seqLengths[seqOrder[ii].pos] = 0;
    }
  }

  if (samPar.desiredFraction > 0.0) {
    fprintf(stderr, "Emitting %.4f fraction of the reads.\n",
            samPar.desiredFraction);

    for (uint64 ii=0; ii<numSeqsTotal; ii++) {
      if (seqOrder[ii].rnd < samPar.desiredFraction)
        ;
      else
        seqLengths[seqOrder[ii].pos] = 0;
    }
  }

  //  Unrandomize the sequences.  Not needed for 2-pass, but needed for 1-pass.

  sort(seqOrder.begin(), seqOrder.end(), seqOrderNormal);
}



void
doSample_paired(vector<char *> &inputs, sampleParameters &samPar) {

  samPar.initialize();

  vector<uint64>    numSeqsPerFile;
  vector<uint64>    seqLengths;
  vector<seqEntry>  seqOrder;

  uint64            numSeqsTotal  = 0;
  uint64            numBasesTotal = 0;

  vector<char *>    names;
  vector<char *>    sequences;
  vector<char *>    qualities;

  mtRandom          MT;

  //  Open output files. If paired, replace #'s in the output names with 1 or 2.

  FILE *outFile1 = NULL;
  FILE *outFile2 = NULL;

  {
    char  *a = strrchr(samPar.output1, '#');
    char  *b = strrchr(samPar.output2, '#');

    if (a == NULL)
      fprintf(stderr, "ERROR: Failed to find '#' in output name '%s'\n", samPar.output1), exit(1);
    if (b == NULL)
      fprintf(stderr, "ERROR: Failed to find '#' in output name '%s'\n", samPar.output2), exit(1);

    *a = '1';
    *b = '2';

    outFile1 = AS_UTL_openOutputFile(samPar.output1);
    outFile2 = AS_UTL_openOutputFile(samPar.output2);
  }

  //  Scan the inputs, saving the number of sequences in each and the length of each sequence.

  dnaSeq   seq1;
  dnaSeq   seq2;

  for (uint32 ff=0; ff<inputs.size(); ff += 2) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff+0]);
    dnaSeqFile  *sf2 = new dnaSeqFile(inputs[ff+1]);
    uint64       num = 0;

    bool   sf1more = sf1->loadSequence(seq1);
    bool   sf2more = sf2->loadSequence(seq2);

    while ((sf1more == true) &&
           (sf2more == true)) {
      seqLengths.push_back(seq1.length() + seq2.length());
      seqOrder.push_back(seqEntry(MT, numSeqsTotal));

      numSeqsTotal  += 1;
      numBasesTotal += seq1.length() + seq2.length();

      num += 1;

      sf1more = sf1->loadSequence(seq1);
      sf2more = sf2->loadSequence(seq2);
    }

    numSeqsPerFile.push_back(num);

    delete sf1;
    delete sf2;
  }

  //  Figure out what to output.

  doSample_sample(samPar, numSeqsTotal, numBasesTotal, seqLengths, seqOrder);

  //  Scan the inputs again, this time emitting sequences if their saved length isn't zero.

  for (uint32 ff=0; ff<inputs.size(); ff += 2) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff+0]);
    dnaSeqFile  *sf2 = new dnaSeqFile(inputs[ff+1]);
    uint64       num = 0;

    bool   sf1more = sf1->loadSequence(seq1);
    bool   sf2more = sf2->loadSequence(seq2);

    while ((sf1more == true) &&
           (sf2more == true)) {
      if (seqLengths[num] > 0) {
        AS_UTL_writeFastA(outFile1, seq1.bases(), seq1.length(), 0, ">%s\n", seq1.name());
        AS_UTL_writeFastA(outFile2, seq2.bases(), seq2.length(), 0, ">%s\n", seq2.name());
      }

      num += 1;

      sf1more = sf1->loadSequence(seq1);
      sf2more = sf2->loadSequence(seq2);
    }

    delete sf1;
    delete sf2;
  }

  AS_UTL_closeFile(outFile1, samPar.output1);
  AS_UTL_closeFile(outFile2, samPar.output2);
}



void
doSample_single(vector<char *> &inputs, sampleParameters &samPar) {

  samPar.initialize();

  vector<uint64>    numSeqsPerFile;
  vector<uint64>    seqLengths;
  vector<seqEntry>  seqOrder;

  uint64            numSeqsTotal  = 0;
  uint64            numBasesTotal = 0;

  vector<char *>    names;
  vector<char *>    sequences;
  vector<char *>    qualities;

  mtRandom          MT;

  //  Open output files. If paired, replace #'s in the output names with 1 or 2.

  FILE *outFile1 = AS_UTL_openOutputFile(samPar.output1);

  //  Scan the inputs, saving the number of sequences in each and the length of each sequence.

  dnaSeq   seq1;

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff]);
    uint64       num = 0;

    while (sf1->loadSequence(seq1)) {
      seqLengths.push_back(seq1.length());
      seqOrder.push_back(seqEntry(MT, numSeqsTotal));

      numSeqsTotal  += 1;
      numBasesTotal += seq1.length();

      num += 1;
    }

    numSeqsPerFile.push_back(num);

    delete sf1;
  }

  //  Figure out what to output.

  doSample_sample(samPar, numSeqsTotal, numBasesTotal, seqLengths, seqOrder);

  //  Scan the inputs again, this time emitting sequences if their saved length isn't zero.

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff]);
    uint64       num = 0;

    while (sf1->loadSequence(seq1)) {
      if (seqLengths[num] > 0)
        AS_UTL_writeFastA(outFile1, seq1.bases(), seq1.length(), 0, ">%s\n", seq1.name());

      num += 1;
    }

    delete sf1;
  }

  AS_UTL_closeFile(outFile1, samPar.output1);
}



void
doSample(vector<char *> &inputs, sampleParameters &samPar) {

  if (samPar.isPaired == false)
    doSample_single(inputs, samPar);
  else
    doSample_paired(inputs, samPar);
}



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


void
doShiftRegister(shiftRegisterParameters &srPar) {
  uint32 len = 0;
  char   sr[17];
  uint8  sv[17];

  srPar.initialize();

  len = srPar.len;
  len = 8;

  for (uint32 ii=0; ii<srPar.len; ii++) {
    sr[ii] =  srPar.sr[ii];
    sv[ii] = (srPar.sv[ii] == '1') ? 1 : 0;
  }

  sr[0] = 'A';  sv[0] = 1;   //  Oldest value
  sr[1] = 'A';  sv[1] = 0;
  sr[2] = 'A';  sv[2] = 0;
  sr[3] = 'A';  sv[3] = 0;
  sr[4] = 'A';  sv[4] = 0;
  sr[5] = 'A';  sv[5] = 0;
  sr[6] = 'A';  sv[6] = 0;
  sr[7] = 'G';  sv[7] = 0;  //  Newest value

  sr[srPar.len] = 0;
  sv[srPar.len] = 0;

  uint32   kmer   = 0;
  uint32  *detect = new uint32 [65536];

  while ((sv[7] != 0) ||
         (sv[6] != 0) ||
         (sv[5] != 0) ||
         (sv[4] != 0) ||
         (sv[3] != 0) ||
         (sv[2] != 0) ||
         (sv[1] != 0) ||
         (sv[0] != 0)) {

    sr[0] = 'A';
    sr[1] = 'A';
    sr[2] = 'A';
    sr[3] = 'A';
    sr[4] = 'A';
    sr[5] = 'A';
    sr[6] = 'A';
    sr[7] = 'G';

    kmer = 0x0003;

#define STRING
#undef  KMER

#ifdef STRING
      fprintf(stdout, "%s", sr);
#endif

    memset(detect, 0, sizeof(uint32) * 65536);
    detect[kmer] = 1;

    for (uint32 ii=0; ii<65536; ii++) {
#ifdef KMER
      fprintf(stdout, "%06u %s 0x%04x %u\n", ii, sr, kmer, detect[kmer]);
#endif

      uint32  next = 0;

      for (uint32 kk=0; kk<len; kk++)
        if (sv[kk])
          next += (sr[kk] >> 1) & 0x03;

      if      ((next & 0x03) == 0x00)  next = 'A';
      else if ((next & 0x03) == 0x01)  next = 'C';
      else if ((next & 0x03) == 0x03)  next = 'G';
      else if ((next & 0x03) == 0x02)  next = 'T';

      for (uint32 kk=0; kk<len-1; kk++)
        sr[kk] = sr[kk+1];

      sr[len-1] = next;

#ifdef STRING
      fprintf(stdout, "%c", next);
#endif

      kmer <<= 2;
      kmer  |= ((sr[len-1] >> 1) & 0x3);
      kmer &= 0xffff;

      detect[kmer]++;

      if (detect[kmer] == 2) {
        if (ii > 0) {
#ifdef KMER
          fprintf(stdout, "%06u %s 0x%04x %u\n", ii, sr, kmer, detect[kmer]);
#endif
#ifdef STRING
          fprintf(stdout, "\n");
#endif
          fprintf(stderr, "cycle at ii=%5u  %d%d%d%d%d%d%d%d\n",
                  ii, sv[0], sv[1], sv[2], sv[3], sv[4], sv[5], sv[6], sv[7]);
        }
        break;
      }
    }

    sv[7]++;

    if (sv[7] > 1)  { sv[7] = 0;  sv[6]++; }
    if (sv[6] > 1)  { sv[6] = 0;  sv[5]++; }
    if (sv[5] > 1)  { sv[5] = 0;  sv[4]++; }
    if (sv[4] > 1)  { sv[4] = 0;  sv[3]++; }
    if (sv[3] > 1)  { sv[3] = 0;  sv[2]++; }
    if (sv[2] > 1)  { sv[2] = 0;  sv[1]++; }
    if (sv[1] > 1)  { sv[1] = 0;  sv[0]++; }
    if (sv[0] > 1)  { break; }
  }

  delete [] detect;
}



int
main(int argc, char **argv) {
  vector<char *>              inputs;

  opMode                      mode = modeUnset;

  summarizeParameters         sumPar;
  extractParameters           extPar;
  generateParameters          genPar;
  sampleParameters            samPar;
  shiftRegisterParameters     srPar;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {

    //  SUMMARIZE

    if      (strcmp(argv[arg], "summarize") == 0) {
      mode = modeSummarize;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-size") == 0)) {
      sumPar.genomeSize = strtoull(argv[++arg], NULL, 10);
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-1x") == 0)) {
      sumPar.limitTo1x = true;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-split-n") == 0)) {
      sumPar.breakAtN = true;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-assequences") == 0)) {
      sumPar.asSequences = true;
      sumPar.asBases     = false;
    }

    else if ((mode == modeSummarize) && (strcmp(argv[arg], "-asbases") == 0)) {
      sumPar.asSequences = false;
      sumPar.asBases     = true;
    }

    //  EXTRACT

    else if (strcmp(argv[arg], "extract") == 0) {
      mode = modeExtract;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-bases") == 0)) {
      decodeRange(argv[++arg], extPar.baseBgn, extPar.baseEnd);
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-sequences") == 0)) {
      decodeRange(argv[++arg], extPar.seqsBgn, extPar.seqsEnd);
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-reverse") == 0)) {
      extPar.asReverse = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-complement") == 0)) {
      extPar.asComplement = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-rc") == 0)) {
      extPar.asReverse = true;
      extPar.asComplement = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-upper") == 0)) {
      extPar.asUpperCase = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-lower") == 0)) {
      extPar.asLowerCase = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-length") == 0)) {
      decodeRange(argv[++arg], extPar.lensBgn, extPar.lensEnd);
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-lowermask") == 0)) {
      extPar.doMasking = true;
      extPar.maskWithN = false;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "-nmask") == 0)) {
      extPar.doMasking = true;
      extPar.maskWithN = true;
    }

    else if ((mode == modeExtract) && (strcmp(argv[arg], "") == 0)) {
    }


    //  GENERATE

    else if (strcmp(argv[arg], "generate") == 0) {
      mode = modeGenerate;
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-min") == 0)) {
      genPar.minLength = strtouint64(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-max") == 0)) {
      genPar.maxLength = strtouint64(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-sequences") == 0)) {
      genPar.nSeqs = strtouint64(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-bases") == 0)) {
      genPar.nBases = strtouint64(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-guassian") == 0)) {
      genPar.useGaussian = true;
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-mirror") == 0)) {
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-gc") == 0)) {
      double  gc = strtodouble(argv[++arg]);
      double  at = 1.0 - gc;

      genPar.gFreq = genPar.cFreq = gc / 2.0;
      genPar.aFreq = genPar.tFreq = at / 2.0;
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-at") == 0)) {
      double  at = strtodouble(argv[++arg]);
      double  gc = 1.0 - at;

      genPar.gFreq = genPar.cFreq = gc / 2.0;
      genPar.aFreq = genPar.tFreq = at / 2.0;
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-a") == 0) ){
      genPar.aFreq = strtodouble(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-c") == 0)) {
      genPar.cFreq = strtodouble(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-g") == 0)) {
      genPar.gFreq = strtodouble(argv[++arg]);
    }

    else if ((mode == modeGenerate) && (strcmp(argv[arg], "-t") == 0)) {
      genPar.tFreq = strtodouble(argv[++arg]);
    }

    //  SIMULATE

    else if (strcmp(argv[arg], "simulate") == 0) {
      mode = modeSimulate;
    }

    //  SAMPLE

    else if (strcmp(argv[arg], "sample") == 0) {
      mode = modeSample;
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-paired") == 0)) {
      samPar.isPaired = true;
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-output") == 0)) {
      strncpy(samPar.output1, argv[++arg], FILENAME_MAX);  //  #'s in the name will be replaced
      strncpy(samPar.output2, argv[  arg], FILENAME_MAX);  //  by '1' or '2' later.
    }


    else if ((mode == modeSample) && (strcmp(argv[arg], "-coverage") == 0)) {      //  Sample reads up to some coverage C
      samPar.desiredCoverage = strtodouble(argv[++arg]);
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-genomesize") == 0)) {
      samPar.genomeSize = strtouint64(argv[++arg]);
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-bases") == 0)) {         //  Sample B bases
      samPar.desiredNumBases = strtouint64(argv[++arg]);
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-reads") == 0)) {         //  Sample N reads
      samPar.desiredNumReads = strtouint64(argv[++arg]);
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-pairs") == 0)) {         //  Sample N pairs of reads
      samPar.desiredNumReads = strtouint64(argv[++arg]) * 2;
    }

    else if ((mode == modeSample) && (strcmp(argv[arg], "-fraction") == 0)) {      //  Sample F fraction
      samPar.desiredFraction = strtodouble(argv[++arg]);
    }

    //  SHIFT

    else if (strcmp(argv[arg], "shift") == 0) {
      mode = modeShift;
    }

    else if ((mode == modeShift) && (strcmp(argv[arg], "-len") == 0)) {
    }

    else if ((mode == modeShift) && (strcmp(argv[arg], "-init") == 0)) {
      strcpy(srPar.sr, argv[++arg]);
    }

    else if ((mode == modeShift) && (strcmp(argv[arg], "-map") == 0)) {
      strcpy(srPar.sv, argv[++arg]);
    }

    else if ((mode == modeShift) && (strcmp(argv[arg], "") == 0)) {
    }

    //  INPUTS

    else if (fileExists(argv[arg]) == true) {
      inputs.push_back(argv[arg]);
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown parameter '%s'\n", argv[arg]);
      err.push_back(s);
    }
    arg++;
  }

  if  (mode == modeUnset)
    err.push_back("ERROR:  No mode (summarize, extract, generate or simulate) specified.\n");
  if (inputs.size() == 0)
    err.push_back("ERROR:  No sequence files supplied.\n");

  if (inputs.size() == 0)
    err.push_back("ERROR:  No input files supplied.\n");



  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [mode] [options] [sequence_file ...]\n", argv[0]);
    fprintf(stderr, "\n");

    if (mode == modeUnset) {
      fprintf(stderr, "MODES:\n");
      fprintf(stderr, "  summarize      report N50, length histogram, mono-, di- and tri-nucleotide frequencies\n");
      fprintf(stderr, "  extract        extract the specified sequences\n");
      fprintf(stderr, "  sample         emit existing sequences randomly\n");
      fprintf(stderr, "  generate       generate random sequences\n");
      fprintf(stderr, "  simulate       errors in existing sequences\n");
      fprintf(stderr, "\n");
    }

    if ((mode == modeUnset) || (mode == modeSummarize)) {
      fprintf(stderr, "OPTIONS for summarize mode:\n");
      fprintf(stderr, "  -size          base size to use for N50 statistics\n");
      fprintf(stderr, "  -1x            limit NG table to 1x coverage\n");
      fprintf(stderr, "  -assequences   load data as complete sequences (for testing)\n");
      fprintf(stderr, "  -asbases       load data as blocks of bases    (for testing)\n");
      fprintf(stderr, "\n");
    }

    if ((mode == modeUnset) || (mode == modeExtract)) {
      fprintf(stderr, "OPTIONS for extract mode:\n");
      fprintf(stderr, "  -bases     baselist extract bases as specified in the 'list' from each sequence\n");
      fprintf(stderr, "  -sequences seqlist  extract ordinal sequences as specified in the 'list'\n");
      fprintf(stderr, "  -reverse            reverse the bases in the sequence\n");
      fprintf(stderr, "  -complement         complement the bases in the sequence\n");
      fprintf(stderr, "  -rc                 alias for -reverse -complement\n");
      fprintf(stderr, "  -upcase\n");
      fprintf(stderr, "  -downcase\n");
      fprintf(stderr, "  -length min max     print sequence if it is at least 'min' bases and at most 'max' bases long\n");
      fprintf(stderr, "  \n");
      fprintf(stderr, "                      a 'baselist' is a set of integers formed from any combination\n");
      fprintf(stderr, "                      of the following, seperated by a comma:\n");
      fprintf(stderr, "                           num       a single number\n");
      fprintf(stderr, "                           bgn-end   a range of numbers:  bgn <= end\n");
      fprintf(stderr, "                      bases are spaced-based; -bases 0-2,4 will print the bases between\n");
      fprintf(stderr, "                      the first two spaces (the first two bases) and the base after the\n");
      fprintf(stderr, "                      fourth space (the fifth base).\n");
      fprintf(stderr, "  \n");
      fprintf(stderr, "                      a 'seqlist' is a set of integers formed from any combination\n");
      fprintf(stderr, "                      of the following, seperated by a comma:\n");
      fprintf(stderr, "                           num       a single number\n");
      fprintf(stderr, "                           bgn-end   a range of numbers:  bgn <= end\n");
      fprintf(stderr, "                      sequences are 1-based; -sequences 1,3-5 will print the first, third,\n");
      fprintf(stderr, "                      fourth and fifth sequences.\n");
      fprintf(stderr, "  \n");
    }

    if ((mode == modeUnset) || (mode == modeSample)) {
      fprintf(stderr, "OPTIONS for sample mode:\n");
      fprintf(stderr, "  -paired             treat inputs as paired sequences; the first two files form the\n");
      fprintf(stderr, "                      first pair, and so on.\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -output O           write output sequences to file O.  If paired, two files must be supplied.\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -coverage C         output C coverage of sequences, based on genome size G.\n");
      fprintf(stderr, "  -genomesize G       \n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -bases B            output B bases.\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -reads R            output R reads.\n");
      fprintf(stderr, "  -pairs P            output P pairs (only if -paired).\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -fraction F         output fraction F of the input bases.\n");
      fprintf(stderr, "\n");
    }

    if ((mode == modeUnset) || (mode == modeGenerate)) {
      fprintf(stderr, "OPTIONS for generate mode:\n");
      fprintf(stderr, "  -min M         minimum sequence length\n");
      fprintf(stderr, "  -max M         maximum sequence length\n");
      fprintf(stderr, "  -sequences N   generate N sequences\n");
      fprintf(stderr, "  -bases B       generate at least B bases, no more than B+maxLength-1 bases.\n");
      fprintf(stderr, "  -gaussian      99.73%% of the reads (3 standard deviations) will be between min and max\n");
      fprintf(stderr, "  -mirror F      \n");
      fprintf(stderr, "  -gc bias       sets GC/AT composition (default 0.50)\n");
      fprintf(stderr, "  -at bias       sets GC/AT composition (default 0.50)\n");
      fprintf(stderr, "  -a freq        sets frequency of A bases (default 0.25)\n");
      fprintf(stderr, "  -c freq        sets frequency of C bases (default 0.25)\n");
      fprintf(stderr, "  -g freq        sets frequency of G bases (default 0.25)\n");
      fprintf(stderr, "  -t freq        sets frequency of T bases (default 0.25)\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "The -gc option is a shortcut for setting all four base frequencies at once.  Order matters!\n");
      fprintf(stderr, "  -gc 0.6 -a 0.1 -t 0.3 -- sets G = C = 0.3, A = 0.1, T = 0.3\n");
      fprintf(stderr, "  -a 0.1 -t 0.3 -gc 0.6 -- sets G = C = 0.3, A = T = 0.15\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "Base frequencies are scaled to sum to 1.0.\n");
      fprintf(stderr, "  -a 1.25 -- results in a sum of 2.0 (1.25 + 0.25 + 0.25 + 0.25) so final frequencies will be:\n");
      fprintf(stderr, "             A =         1.25/2 = 0.625\n");
      fprintf(stderr, "             C = G = T = 0.25/2 = 0.125.\n");
      fprintf(stderr, "  -gc 0.8 -a 1.0 -t 0.2 -- sum is also 2.0, final frequencies will be:\n");
      fprintf(stderr, "             A =         1.00/2 = 0.5\n");
      fprintf(stderr, "             C = G =     0.40/2 = 0.2\n");
      fprintf(stderr, "             T =         0.20/2 = 0.1\n");
      fprintf(stderr, "\n");
    }

    if ((mode == modeUnset) || (mode == modeSimulate)) {
      fprintf(stderr, "OPTIONS for simulate mode:\n");
      fprintf(stderr, "\n");
    }

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  sumPar.finalize();
  genPar.finalize();
  extPar.finalize();

  switch (mode) {
    case modeSummarize:
      doSummarize(inputs, sumPar);
      break;
    case modeExtract:
      doExtract(inputs, extPar);
      break;
    case modeGenerate:
      doGenerate(genPar);
      break;
    case modeSimulate:
      doSimulate(inputs);
      break;
    case modeSample:
      doSample(inputs, samPar);
      break;
    case modeShift:
      doShiftRegister(srPar);
      break;
    default:
      break;
  }

  return(0);
}
