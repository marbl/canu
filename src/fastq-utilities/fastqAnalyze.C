
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
 *    src/AS_GKP/fastqAnalyze.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-FEB-24 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2015-JAN-13
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-03
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include <vector>
#include <algorithm>

using namespace std;


//  Map a letter to an index in the freq arrays.
uint32   baseToIndex[256];
char     indexToBase[8];

#define  FREQ_A    0
#define  FREQ_C    1
#define  FREQ_G    2
#define  FREQ_T    3
#define  FREQ_N    4  //  Just N's
#define  FREQ_g    5  //  Just -'s
#define  FREQ_Z    6  //  Everything else
#define  FREQ_NUM  7

#define  MAX_READ_LEN   1024 * 1024

class nucFreq {
public:
  nucFreq() {
    memset(mono, 0, sizeof(uint64) * FREQ_NUM);
    memset(di,   0, sizeof(uint64) * FREQ_NUM * FREQ_NUM);
    memset(tri,  0, sizeof(uint64) * FREQ_NUM * FREQ_NUM * FREQ_NUM);
  };

  uint64    mono[FREQ_NUM];
  uint64    di[FREQ_NUM][FREQ_NUM];
  uint64    tri[FREQ_NUM][FREQ_NUM][FREQ_NUM];
};

class nucOut {
public:
  nucOut(char a, char b, char c, uint64 cnt, double frq) {
    label[0] = a;
    label[1] = b;
    label[2] = c;
    label[3] = 0;

    count    = cnt;
    freq     = frq;
  };

  bool operator<(nucOut const &that) const {
    return(freq > that.freq);
  };

  char    label[4];
  uint64  count;
  double  freq;
};






void
doStats(char *inName,
        char *otName) {

  uint64           totSeqs  = 0;
  uint64           totBases = 0;
  vector<uint32>   seqLen;
  nucFreq         *freq = new nucFreq;

  char   A[MAX_READ_LEN];
  char   B[MAX_READ_LEN];
  char   C[MAX_READ_LEN];
  char   D[MAX_READ_LEN];

  errno = 0;
  FILE *F = fopen(inName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", inName, strerror(errno)), exit(1);

  //errno = 0;
  //FILE *O = fopen(otName, "w");
  //if (errno)
  //  fprintf(stderr, "Failed to open '%s' for writing: %s\n", otName, strerror(errno)), exit(1);

  while (!feof(F)) {
    fgets(A, MAX_READ_LEN, F);
    fgets(B, MAX_READ_LEN, F);  chomp(B);
    fgets(C, MAX_READ_LEN, F);
    fgets(D, MAX_READ_LEN, F);  chomp(D);

    if ((A[0] != '@') || (C[0] != '+')) {
      fprintf(stderr, "WARNING:  sequence isn't fastq.\n");
      fprintf(stderr, "WARNING:  %s",   A);
      fprintf(stderr, "WARNING:  %s\n", B);
      fprintf(stderr, "WARNING:  %s",   C);
      fprintf(stderr, "WARNING:  %s\n", D);
    }

    uint32  a   = 0;
    uint32  b   = baseToIndex[B[0]];
    uint32  c   = baseToIndex[B[1]];
    uint32  ii;

    freq->mono[b]++;
    freq->mono[c]++;

    freq->di[b][c]++;

    for (ii=2; B[ii]; ii++) {
      a = b;
      b = c;
      c = baseToIndex[B[ii]];

      freq->mono[c]++;
      freq->di[b][c]++;
      freq->tri[a][b][c]++;
    }

    ii--;

    seqLen.push_back(ii);

    totSeqs++;
    totBases += ii;

    if ((totSeqs % 10000) == 0)
      fprintf(stderr, "Reading "F_U64"\r", totSeqs);
  }

  fprintf(stderr, "Read    "F_U64"\n", totSeqs);

  fprintf(stdout, "%s\n", inName);
  fprintf(stdout, "\n");

  fprintf(stdout, "sequences\t"F_U64"\n", totSeqs);
  fprintf(stdout, "bases\t"F_U64"\n",     totBases);
  fprintf(stdout, "\n");
  fprintf(stdout, "average\t"F_U64"\n", totBases / totSeqs);
  fprintf(stdout, "\n");

  //sort(seqLen.begin(), seqLen.end());

  uint64   min        = UINT64_MAX;
  uint64   max        = 0;

  for (uint32 ii=0; ii<seqLen.size(); ii++) {
    min = MIN(min, seqLen[ii]);
    max = MAX(max, seqLen[ii]);
  }

  uint64  *histogram  = new uint64 [max + 1];

  for (uint32 ii=0; ii<=max; ii++)
    histogram[ii] = 0;

  for (uint32 ii=0; ii<seqLen.size(); ii++)
    histogram[seqLen[ii]]++;

  for (uint32 ii=min; ii<=max; ii++)
    fprintf(stdout, F_U32"\t"F_U64"\n", ii, histogram[ii]);

  delete [] histogram;

  vector<nucOut>    output;

  fprintf(stdout, "\n");
  fprintf(stdout, "mononucleotide\n");
  fprintf(stdout, "\n");
  for (uint32 ii=0; ii<FREQ_NUM; ii++)
    if (freq->mono[ii] > 0)
      output.push_back(nucOut(indexToBase[ii], 0, 0,
                              freq->mono[ii],
                              freq->mono[ii] * 100.0 / totBases));
  sort(output.begin(), output.end());
  for (uint32 ii=0; ii<output.size(); ii++)
    fprintf(stdout, "%s\t"F_U64"\t%.4f%%\n", output[ii].label, output[ii].count, output[ii].freq);
  output.clear();

  fprintf(stdout, "\n");
  fprintf(stdout, "dinucleotide\n");
  fprintf(stdout, "\n");
  for (uint32 ii=0; ii<FREQ_NUM; ii++)
    for (uint32 jj=0; jj<FREQ_NUM; jj++)
      if (freq->di[ii][jj] > 0)
        output.push_back(nucOut(indexToBase[ii], indexToBase[jj], 0,
                                freq->di[ii][jj],
                                freq->di[ii][jj] * 100.0 / totBases));
  sort(output.begin(), output.end());
  for (uint32 ii=0; ii<output.size(); ii++)
    fprintf(stdout, "%s\t"F_U64"\t%.4f%%\n", output[ii].label, output[ii].count, output[ii].freq);
  output.clear();

  fprintf(stdout, "\n");
  fprintf(stdout, "trinucleotide\n");
  fprintf(stdout, "\n");
  for (uint32 ii=0; ii<FREQ_NUM; ii++)
    for (uint32 jj=0; jj<FREQ_NUM; jj++)
      for (uint32 kk=0; kk<FREQ_NUM; kk++)
        if (freq->tri[ii][jj][kk] > 0)
          output.push_back(nucOut(indexToBase[ii], indexToBase[jj], indexToBase[kk],
                                  freq->tri[ii][jj][kk],
                                  freq->tri[ii][jj][kk] * 100.0 / totBases));
  sort(output.begin(), output.end());
  for (uint32 ii=0; ii<output.size(); ii++)
    fprintf(stdout, "%s\t"F_U64"\t%.4f%%\n", output[ii].label, output[ii].count, output[ii].freq);
  output.clear();

  //fclose(O);
  fclose(F);

  delete freq;
}



void
doAnalyzeQV(char *inName,
            bool &originalIsSolexa,
            bool &originalIsIllumina,
            bool &originalIsSanger) {

  if ((originalIsSolexa   == true) ||
      (originalIsIllumina == true) ||
      (originalIsSanger   == true))
    return;

  uint32 numValid  = 5;
  uint32 numTrials = 100000;

  char   A[MAX_READ_LEN];
  char   B[MAX_READ_LEN];
  char   C[MAX_READ_LEN];
  char   D[MAX_READ_LEN];

  errno = 0;
  FILE *F = fopen(inName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", inName, strerror(errno)), exit(1);

  //  Initially, it could be any of these.
  //
  bool   isNotSanger    = false;
  bool   isNotSolexa    = false;
  bool   isNotIllumina3 = false;
  bool   isNotIllumina5 = false;
  bool   isNotIllumina8 = false;

  uint32 qvCounts[256]  = {0};

  while ((!feof(F)) &&
         (numValid != 1) &&
         (numTrials > 0)) {
    fgets(A, MAX_READ_LEN, F);
    fgets(B, MAX_READ_LEN, F);
    fgets(C, MAX_READ_LEN, F);
    fgets(D, MAX_READ_LEN, F);  chomp(D);

    for (uint32 x=0; D[x] != 0; x++) {
      if (D[x] < '!')  isNotSanger    = true;
      if (D[x] < ';')  isNotSolexa    = true;
      if (D[x] < '@')  isNotIllumina3 = true;  //  Illumina 1.3
      if (D[x] < 'B')  isNotIllumina5 = true;  //  Illumina 1.5
      if (D[x] < '!')  isNotIllumina8 = true;  //  Illumina 1.5

      if ('I' < D[x])  isNotSanger    = true;
      if ('h' < D[x])  isNotSolexa    = true;
      if ('h' < D[x])  isNotIllumina3 = true;  //  Illumina 1.3
      if ('h' < D[x])  isNotIllumina5 = true;  //  Illumina 1.5
      if ('J' < D[x])  isNotIllumina8 = true;  //  Illumina 1.5

      qvCounts[D[x]]++;

      //fprintf(stderr, "%d%d%d%d%d %c %d\n",
      //        isNotSanger, isNotSolexa, isNotIllumina3, isNotIllumina5, isNotIllumina8, D[x], x);
    }

    //fprintf(stderr, "%d%d%d%d%d '%s'\n",
    //        isNotSanger, isNotSolexa, isNotIllumina3, isNotIllumina5, isNotIllumina8, D);

    numValid = 0;

    if (isNotSanger    == false)  numValid++;
    if (isNotSolexa    == false)  numValid++;
    if (isNotIllumina3 == false)  numValid++;
    if (isNotIllumina5 == false)  numValid++;
    if (isNotIllumina8 == false)  numValid++;

    numTrials--;
  }

  fclose(F);

  fprintf(stdout, "%s --", inName);

  if (isNotSanger    == false)  fprintf(stdout, " SANGER");
  if (isNotSolexa    == false)  fprintf(stdout, " SOLEXA");
  if (isNotIllumina3 == false)  fprintf(stdout, " ILLUMINA_1.3+");
  if (isNotIllumina5 == false)  fprintf(stdout, " ILLUMINA_1.5+");
  if (isNotIllumina8 == false)  fprintf(stdout, " ILLUMINA_1.8+");
  if (numValid       == 0)      fprintf(stdout, " NO_VALID_ENCODING");

  fprintf(stdout, "\n");

  if (numValid == 0) {
    fprintf(stdout, "QV histogram:\n");

   for (uint32 c=0, i=0; i<12; i++) {
      fprintf(stdout, "%3d: ", i * 10);

      for (uint32 j=0; j<10; j++, c++)
        fprintf(stdout, " %8d/%c", qvCounts[c], isprint(c) ? c : '.');

      fprintf(stdout, "\n");
    }
  }

  if (isNotSanger    == false)  originalIsSanger   = true;
  if (isNotSolexa    == false)  originalIsSolexa   = true;
  if (isNotIllumina3 == false)  originalIsIllumina = true;
  if (isNotIllumina5 == false)  originalIsIllumina = true;
  if (isNotIllumina8 == false)  originalIsSanger   = true;
}



void
doTransformQV(char *inName,
              char *otName,
              bool  originalIsSolexa,
              bool  originalIsIllumina,
              bool  originalIsSanger) {

  uint32 numValid = 0;

  char   A[MAX_READ_LEN];
  char   B[MAX_READ_LEN];
  char   C[MAX_READ_LEN];
  char   D[MAX_READ_LEN];

  if (originalIsSolexa    == true)  numValid++;
  if (originalIsIllumina  == true)  numValid++;
  if (originalIsSanger    == true)  numValid++;

  if (numValid == 0)
    fprintf(stderr, "No QV decision made.  No valid encoding found.  Specify a QV encoding to convert from.\n"), exit(0);

  if (numValid > 1)
    fprintf(stderr, "No QV decision made.  Multiple valid encodings found.  Specify a QV encoding to convert from.\n"), exit(0);


  if (originalIsSanger == true)
    fprintf(stderr, "No QV changes needed; original is in sanger format already.\n"), exit(0);

  errno = 0;
  FILE *F = fopen(inName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", inName, strerror(errno)), exit(1);

  errno = 0;
  FILE *O = fopen(otName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", otName, strerror(errno)), exit(1);

  while (!feof(F)) {
    fgets(A, MAX_READ_LEN, F);
    fgets(B, MAX_READ_LEN, F);
    fgets(C, MAX_READ_LEN, F);
    fgets(D, MAX_READ_LEN, F);  chomp(D);

    if (feof(F))
      break;

    for (uint32 x=0; D[x] != 0; x++) {
      if (originalIsSolexa) {
        double qs  = D[x] - '@';
        qs /= 10.0;
        qs  = 10.0 * log10(pow(10.0, qs) + 1);
        D[x] = lround(qs) + '0';
      }

      if (originalIsIllumina) {
        D[x] -= '@';
        D[x] += '!';
      }
    }

    fprintf(O, "%s%s%s%s\n", A, B, C, D);
  }

  fclose(F);
  fclose(O);
}




int
main(int argc, char **argv) {
  char   *inName = NULL;
  char   *otName = NULL;

  bool    originalIsSolexa   = false;
  bool    originalIsIllumina = false;
  bool    originalIsSanger   = false;

  bool    convertToSanger    = false;

  bool    computeStats       = false;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      otName          = argv[++arg];
      convertToSanger = true;

    } else if (strcmp(argv[arg], "-solexa") == 0) {
      originalIsSolexa = true;

    } else if (strcmp(argv[arg], "-illumina") == 0) {
      originalIsIllumina = true;

    } else if (strcmp(argv[arg], "-sanger") == 0) {
      originalIsSanger = true;

    } else if (strcmp(argv[arg], "-stats") == 0) {
      computeStats = true;

    } else if (inName == NULL) {
      inName = argv[arg];

    } else {
      err++;
    }

    arg++;
  }
  if ((err) || (inName == NULL)) {
    fprintf(stderr, "usage: %s [-stats] [-o output.fastq] input.fastq\n", argv[0]);
    fprintf(stderr, "  If no options are given, input.fastq is analyzed and a best guess for the\n");
    fprintf(stderr, "  QV encoding is output.  Otherwise, the QV encoding is converted to Sanger-style\n");
    fprintf(stderr, "  using this guess.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  In some cases, the encoding cannot be determined.  When this occurs, no guess is\n");
    fprintf(stderr, "  output.  For conversion, you can force the input QV type with:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -solexa     input QV is solexa\n");
    fprintf(stderr, "  -illumina   input QV is illumina\n");
    fprintf(stderr, "  -sanger     input QV is sanger\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o          sanger-style-output.fastq\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  If -stats is supplied, no QV analysis or conversion is performed, but some simple\n");
    fprintf(stderr, "  statistics are computed and output to stdout.\n");
    fprintf(stderr, "\n");

    exit(1);
  }

  for (uint32 ii=0; ii<256; ii++)
    baseToIndex[ii] = FREQ_Z;

  baseToIndex['a'] = FREQ_A;
  baseToIndex['A'] = FREQ_A;
  baseToIndex['c'] = FREQ_C;
  baseToIndex['C'] = FREQ_C;
  baseToIndex['g'] = FREQ_G;
  baseToIndex['G'] = FREQ_G;
  baseToIndex['t'] = FREQ_T;
  baseToIndex['T'] = FREQ_T;
  baseToIndex['n'] = FREQ_N;
  baseToIndex['N'] = FREQ_N;
  baseToIndex['-'] = FREQ_g;
  baseToIndex['-'] = FREQ_g;

  indexToBase[FREQ_A] = 'A';
  indexToBase[FREQ_C] = 'C';
  indexToBase[FREQ_G] = 'G';
  indexToBase[FREQ_T] = 'T';
  indexToBase[FREQ_N] = 'N';
  indexToBase[FREQ_g] = '-';
  indexToBase[FREQ_Z] = '?';


  if (computeStats)
    doStats(inName, otName), exit(0);

  doAnalyzeQV(inName,
              originalIsSolexa,
              originalIsIllumina,
              originalIsSanger);

  if (convertToSanger)
    doTransformQV(inName,
                  otName,
                  originalIsSolexa,
                  originalIsIllumina,
                  originalIsSanger);

  exit(0);
}
