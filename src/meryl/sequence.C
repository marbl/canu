
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
 *    Brian P. Walenz beginning on 2018-JAN-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "libmeryl.H"
#include "libsequence.H"

#include "mt19937ar.H"
#include "AS_UTL_fasta.H"

#include <vector>
#include <algorithm>

using namespace std;


enum opMode {
  modeSummarize,
  modeExtract,
  modeGenerate,
  modeSimulate,
  modeSample,
  modeUnset
};


class generateParameters {
public:
  generateParameters() {
    minLength      = 0;
    maxLength      = 10000;

    nSeqs          = 0;
    nBases         = 0;

    useGaussian    = true;
    gMean          = 0;
    gStdDev        = 0;

    useMirror      = false;
    mirrorInput    = NULL;

    aFreq          = 0.25;
    cFreq          = 0.25;
    gFreq          = 0.25;
    tFreq          = 0.25;
  };

  ~generateParameters() {
  };


  void      finalize(void) {

    //  Check for invalid.

    if ((nSeqs == 0) && (nBases == 0))
      fprintf(stderr, "ERROR:  Told to generate infinite sequences and infinite bases.\n"), exit(1);

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




void
doSummarize(vector<char *> &inputs,
            uint64          genomeSize) {

  vector<uint64>  lengths;

  uint64          nSeqs  = 0;
  uint64          nBases = 0;

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

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf = new dnaSeqFile(inputs[ff]);

    uint64  len = sf->loadSequence(name, nameMax, seq, qlt, seqMax);
    uint32  mer = 0;

    //  If sequence,
    //    Count mono-, di- and tri-nucleotides.
    //    Count number of mono-, di- and tri-nucleotides.
    //    Count number of sequences and total bases.
    //    Save the lengths of sequences.

    while (len > 0) {
      uint64  pos = 0;

      if (pos < len) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
      }

      if (pos < len) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
        dn[mer & 0x0f]++;
      }

      while (pos < len) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
        dn[mer & 0x0f]++;
        tn[mer & 0x3f]++;
      }

      nmn +=                 (len-0);
      ndn += (len < 2) ? 0 : (len-1);
      ntn += (len < 3) ? 0 : (len-2);

      nSeqs  += 1;
      nBases += len;

      lengths.push_back(len);

      len = sf->loadSequence(name, nameMax, seq, qlt, seqMax);
    }

    //  All done!

    delete sf;
  }

  delete [] name;
  delete [] seq;
  delete [] qlt;

  //  Finalize.

  sort(lengths.begin(), lengths.end(), greater<uint64>());

  if (genomeSize == 0)
    genomeSize = nBases;

  uint64   lSum  = 0;                          //  Sum of the lengths we've encountered so far

  uint32   nStep = 10;                         //  Step of each N report.
  uint32   nVal  = nStep;                      //  Index of the threshold we're next printing.
  uint64   nThr  = genomeSize * nVal / 100;    //  Threshold lenth; if sum is bigger, emit and move to the next threshold

  ////////////////////////////////////////

  uint32   nCols    = 63;   //  Magic number to make the histogram the same width as the trinucleotide list
  uint32   nRows    = 0;
  uint32   nRowsMin = 50;   //  Nothing really magic, just fits on the screen.

  //  Count the number of rows we expect to get in the N histogram.

  for (uint32 ii=0; ii<nSeqs; ii++) {
    lSum += lengths[ii];

    while (lSum >= nThr) {
      nRows++;

      if      (nVal <   200)   nVal += nStep;
      else if (nVal <  2000)   nVal += nStep * 10;
      else if (nVal < 20000)   nVal += nStep * 100;
      else                     nVal += nStep * 1000;

      nThr  = genomeSize * nVal / 100;
    }
  }

  if (nRows < nRowsMin)
    nRows = nRowsMin;

  lSum  = 0;
  nStep = 10;
  nVal  = nStep;
  nThr  = genomeSize * nVal / 100;

  //  Generate the length histogram.

  char **histPlot = new char * [nRows + 1];

  for (uint32 rr=0; rr<nRows+1; rr++)           //  28 is a magic number based on the histPlot[] format string below.
    histPlot[rr] = new char [28 + nCols + 1];   //  28 = 9 + 1 + 9 + 1 + 7 + 1

  uint64   maxLength = lengths[0];
  uint32   maxCount  = 0;

  uint32  *nSeqPerLen = new uint32 [nRows + 1];

  for (uint32 rr=0; rr<nRows+1; rr++)                   //  Clear the histogram.
    nSeqPerLen[rr] = 0;

  for (uint32 ii=0; ii<nSeqs; ii++)                     //  Count number of sequences per column.
    nSeqPerLen[nRows * lengths[ii] / maxLength]++;

  for (uint32 rr=0; rr<nRows+1; rr++)                   //  Find the maximum count.
    if (maxCount < nSeqPerLen[rr])
      maxCount = nSeqPerLen[rr];

  for (uint32 rr=0; rr<nRows+1; rr++) {                 //  Generate histogram in text.
    uint32  nn = (uint32)ceil(nSeqPerLen[rr] * nCols / (double)maxCount);

    sprintf(histPlot[rr], "%9" F_U64P "-%-9" F_U64P " %7" F_U32P "|",
            (rr + 0) * maxLength / nRows,
            (rr + 1) * maxLength / nRows - 1,
            nSeqPerLen[rr]);

    for (uint32 cc=0; cc<nn; cc++)
      histPlot[rr][28 + cc] = '-';

    histPlot[rr][28 + nn]   = 0;
  }

  //  Output N table, with length histogram appended at the end of each line.

  uint32  hp = 0;

  fprintf(stdout, "\n");
  fprintf(stdout, "G=%-12" F_U64P "                     sum of  ||               length     num\n", genomeSize);
  fprintf(stdout,   "NG         length     index       lengths  ||                range    seqs\n");
  fprintf(stdout,   "----- ------------ --------- ------------  ||  ------------------- -------\n");

  for (uint32 ii=0; ii<nSeqs; ii++) {
    lSum += lengths[ii];

    while (lSum >= nThr) {
      fprintf(stdout, "%05"    F_U32P " %12" F_U64P " %9" F_U32P " %12" F_U64P "  ||  %s\n",
              nVal, lengths[ii], ii, lSum,
              histPlot[hp++]);

      if      (nVal <   200)   nVal += nStep;
      else if (nVal <  2000)   nVal += nStep * 10;
      else if (nVal < 20000)   nVal += nStep * 100;
      else                     nVal += nStep * 1000;

      nThr  = genomeSize * nVal / 100;
    }
  }

  fprintf(stdout, "%07.3fx           %9" F_U64P " %12" F_U64P "  ||  %s\n", (double)lSum / genomeSize, nSeqs, lSum, histPlot[hp++]);
  //fprintf(stdout,    "                                           ||  %s\n", histPlot[hp++]);
  //fprintf(stdout,   "                 genome-size %12" F_U64P "  ||  %s\n", genomeSize, histPlot[hp++]);

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

  fprintf(stdout, "--------------------- --------------------- ---------------------- ----------------------- ----------------------- -----------------------\n");
  fprintf(stdout, "       mononucleotide          dinucleotide                                                                                  trinucleotide\n");
  fprintf(stdout, "--------------------- --------------------- ---------------------- ----------------------- ----------------------- -----------------------\n");
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
doExtract(vector<char *> &inputs) {
}







void
doGenerate(generateParameters &parameters) {

  parameters.finalize();

  mtRandom   MT;

  uint64  nSeqs  = 0;
  uint64  nBases = 0;

  uint64  seqLen = 0;
  uint64  seqMax = 65536;
  char   *seq    = new char  [seqMax + 1];
  uint8  *qlt    = new uint8 [seqMax + 1];

  double Athresh = parameters.aFreq;
  double Cthresh = parameters.aFreq + parameters.cFreq;
  double Gthresh = parameters.aFreq + parameters.cFreq + parameters.gFreq;
  double Tthresh = parameters.aFreq + parameters.cFreq + parameters.gFreq + parameters.tFreq;

  while ((nSeqs  < parameters.nSeqs) &&
         (nBases < parameters.nBases)) {
    double   len = MT.mtRandomGaussian(parameters.gMean, parameters.gStdDev);

    while (len < 0.5)
      len = MT.mtRandomGaussian(parameters.gMean, parameters.gStdDev);

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
                      ">random%08 "F_U64P"\n", nSeqs);

    nSeqs  += 1;
    nBases += seqLen;
  }
}



void
doSimulate(vector<char *> &inputs) {
}



void
doSample(vector<char *> &inputs) {
}




int
main(int argc, char **argv) {
  vector<char *>      inputs;

  uint64              genomeSize = 0;

  opMode              mode = modeUnset;

  generateParameters  gPar;


  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {

    //  SUMMARIZE

    if      (strcmp(argv[arg], "summarize") == 0) {
      mode = modeSummarize;
    }

    else if (strcmp(argv[arg], "-gs") == 0) {
      genomeSize = strtoull(argv[++arg], NULL, 10);
    }

    //  EXTRACT

    else if (strcmp(argv[arg], "extract") == 0) {
      mode = modeExtract;
    }

    //  GENERATE

    else if (strcmp(argv[arg], "generate") == 0) {
      mode = modeGenerate;
    }

    else if (strcmp(argv[arg], "-min") == 0) {
      gPar.minLength = strtouint64(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-max") == 0) {
      gPar.maxLength = strtouint64(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-sequences") == 0) {
      gPar.nSeqs = strtouint64(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-bases") == 0) {
      gPar.nBases = strtouint64(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-guassian") == 0) {
      gPar.useGaussian = true;
    }

    else if (strcmp(argv[arg], "-mirror") == 0) {
    }

    else if (strcmp(argv[arg], "-gc") == 0) {
      double  gc = strtodouble(argv[++arg]);
      double  at = 1.0 - gc;

      gPar.gFreq = gPar.cFreq = gc / 2.0;
      gPar.aFreq = gPar.tFreq = at / 2.0;
    }

    else if (strcmp(argv[arg], "-at") == 0) {
      double  at = strtodouble(argv[++arg]);
      double  gc = 1.0 - at;

      gPar.gFreq = gPar.cFreq = gc / 2.0;
      gPar.aFreq = gPar.tFreq = at / 2.0;
    }

    else if (strcmp(argv[arg], "-a") == 0) { 
      gPar.aFreq = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-c") == 0) {
      gPar.cFreq = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-g") == 0) {
      gPar.gFreq = strtodouble(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-t") == 0) {
      gPar.tFreq = strtodouble(argv[++arg]);
    }

    //  SIMULATE

    else if (strcmp(argv[arg], "simulate") == 0) {
      mode = modeSimulate;
    }

    else if (strcmp(argv[arg], "sample") == 0) {
      mode = modeSimulate;
    }


    else if (strcmp(argv[arg], "") == 0) {
    }

    else if (strcmp(argv[arg], "") == 0) {
    }

    else if (AS_UTL_fileExists(argv[arg]) == true) {
      inputs.push_back(argv[arg]);
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "ERROR:  Unknown parameter '%s'\n", argv[arg]);
      err.push_back(s);
    }
    arg++;
  }

  if (mode == modeUnset)
    err.push_back("ERROR:  No mode (summarize, extract, generate or simulate) specified.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [mode] [options] [sequence_file ...]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "MODES:\n");
    fprintf(stderr, "  summarize      report N50, length histogram, mono-, di- and tri-nucleotide frequencies\n");
    fprintf(stderr, "  extract        extract the specified sequences\n");
    fprintf(stderr, "  generate       generate random sequences\n");
    fprintf(stderr, "  simulate       errors in existing sequences\n");
    fprintf(stderr, "  sample         emit existing sequences randomly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONS for summarize mode:\n");
    fprintf(stderr, "  -gs            genome size to use for N50 denominator\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONS for extract mode:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
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
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONS for simulate mode:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);
  }

  switch (mode) {
    case modeSummarize:
      doSummarize(inputs, genomeSize);
      break;
    case modeExtract:
      doExtract(inputs);
      break;
    case modeGenerate:
      doGenerate(gPar);
      break;
    case modeSimulate:
      doSimulate(inputs);
      break;
    case modeSample:
      doSample(inputs);
      break;
    default:
      break;
  }

  return(0);
}
