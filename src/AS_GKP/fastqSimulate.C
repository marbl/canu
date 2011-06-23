
/**************************************************************************
 * Copyright (C) 2010, J Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id: fastqSimulate.C,v 1.8 2011-06-23 09:21:03 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"


static char reverseComplement[256];
static char errorBase[256][3];


#define QV_BASE  '!'


//  Returns random int in range bgn <= x < end.
//
int32
randomUniform(int32 bgn, int32 end) {
  return((end - bgn) * drand48() + bgn);
}



//  Generate a random gaussian using the Marsaglia polar method.
//
double
randomGaussian(double mean, double stddev) {
  double  u = 0.0;
  double  v = 0.0;
  double  r = 0.0;

  do {
    u = 2.0 * drand48() - 1.0;
    v = 2.0 * drand48() - 1.0;
    r = u * u + v * v;
  } while (r >= 1.0);

  if (r < 1e-10)
    r = 1e-10;

  r = sqrt(-2 * log(r) / r);

  //  Uniform gaussian is u*r and v*r.  We only use one of these.

  return(mean + u * r * stddev);
}





void
makeSequences(char    *frag,
              int32    fragLen,
              int32    readLen,
              char    *s1,
              char    *q1,
              char    *s2,
              char    *q2) {
  double perr = 0.01;

  assert(fragLen > 0);
  assert(readLen > 0);
  assert(fragLen > readLen);

  //  Build the reads.

  //fprintf(stderr, "makeSequences()-- fragLen=%d readLen=%d\n", fragLen, readLen);

  for (int32 p=0, i=0; p<readLen; p++, i++) {
    s1[p] = frag[i];
    q1[p] = QV_BASE + 39;

    if (drand48() < perr) {
      s1[p] = errorBase[s1[p]][randomUniform(0, 3)];
      q1[p] = QV_BASE + 8;
    }

    assert(s1[p] != '*');
  }

  for (int32 p=0, i=fragLen-1; p<readLen; p++, i--) {
    s2[p] = reverseComplement[frag[i]];
    q2[p] = QV_BASE + 39;

    if (drand48() < perr) {
      s2[p] = errorBase[s2[p]][randomUniform(0, 3)];
      q2[p] = QV_BASE + 8;
    }

    assert(s2[p] != '*');
  }  

  s1[readLen] = 0;
  q1[readLen] = 0;
  s2[readLen] = 0;
  q2[readLen] = 0;
}






void
makePE(char   *seq,
       int32   seqLen,
       FILE   *output,
       FILE   *output1,
       FILE   *output2,
       int32   readLen,
       int32   readPairs,
       int32   peShearSize,
       int32   peShearStdDev) {
  char   *s1 = new char [readLen + 1];
  char   *q1 = new char [readLen + 1];
  char   *s2 = new char [readLen + 1];
  char   *q2 = new char [readLen + 1];

  for (int32 np=0; np<readPairs; np++) {
    int32   len = randomGaussian(peShearSize, peShearStdDev);
    int32   bgn = randomUniform(0, seqLen - len);

    //  Scan the sequence, if we spanned a sequence break, don't use this pair

    for (int32 i=bgn; i<bgn+len; i++)
      if (seq[i] == '>')
        bgn = len = 0;

    if (len == 0) {
      //  Not a valid shearing.
      np--;
      continue;
    }

    if (len <= readLen) {
      np--;
      continue;
    }

    //  Read sequences from the ends.

    makeSequences(seq + bgn, len, readLen, s1, q1, s2, q2);

    //  Output sequences, with a descriptive ID

    fprintf(output, "@PE_%d_%d-%d/1\n", np, bgn, bgn+len);
    fprintf(output, "%s\n", s1);
    fprintf(output, "+\n");
    fprintf(output, "%s\n", q1);

    fprintf(output, "@PE_%d_%d-%d/2\n", np, bgn, bgn+len);
    fprintf(output, "%s\n", s2);
    fprintf(output, "+\n");
    fprintf(output, "%s\n", q2);

    fprintf(output1, "@PE_%d_%d-%d/1\n", np, bgn, bgn+len);
    fprintf(output1, "%s\n", s1);
    fprintf(output1, "+\n");
    fprintf(output1, "%s\n", q1);

    fprintf(output2, "@PE_%d_%d-%d/2\n", np, bgn, bgn+len);
    fprintf(output2, "%s\n", s2);
    fprintf(output2, "+\n");
    fprintf(output2, "%s\n", q2);
  }

  delete [] s1;
  delete [] q1;
  delete [] s2;
  delete [] q2;
}





void
makeMP(char   *seq,
       int32   seqLen,
       FILE   *output,
       FILE   *output1,
       FILE   *output2,
       int32   readLen,
       int32   readPairs,
       int32   mpInsertSize,
       int32   mpInsertStdDev,
       int32   mpShearSize,
       int32   mpShearStdDev,
       double  mpEnrichment,
       double  mpJunction) {
  char   *s1 = new char [readLen + 1];
  char   *q1 = new char [readLen + 1];
  char   *s2 = new char [readLen + 1];
  char   *q2 = new char [readLen + 1];
  char   *sh = new char [1048576];

  for (int32 np=0; np<readPairs; np++) {
    int32   len = randomGaussian(mpInsertSize, mpInsertStdDev);
    int32   bgn = randomUniform(0, seqLen - len);

    //  Scan the sequence, if we spanned a sequence break, don't use this pair

    for (int32 i=bgn; i<bgn+len; i++)
      if (seq[i] == '>')
        bgn = len = 0;

    if (len == 0) {
      //  Not a valid shearing.
      np--;
      continue;
    }

    //  Compute the size of the shearing

    int32   slen = randomGaussian(mpShearSize, mpShearStdDev);  //  shear size

    if (slen <= readLen) {
      np--;
      continue;
    }

    //  If we fail the mpEnrichment test, pick a random shearing and return PE reads.
    //  Otherwise, rotate the sequence to circularize and return MP reads.

    if (mpEnrichment < drand48()) {
      //  Failed to wash away non-biotin marked sequence, make PE
      int32  sbgn = bgn + randomUniform(0, len - slen);

      makeSequences(seq + sbgn, slen, readLen, s1, q1, s2, q2);

      //  Output sequences, with a descriptive ID

      fprintf(output, "@fPE_%d_%d-%d/1\n", np, sbgn, sbgn+slen);
      fprintf(output, "%s\n", s1);
      fprintf(output, "+\n");
      fprintf(output, "%s\n", q1);

      fprintf(output, "@fPE_%d_%d-%d/2\n", np, sbgn, sbgn+slen);
      fprintf(output, "%s\n", s2);
      fprintf(output, "+\n");
      fprintf(output, "%s\n", q2);

      fprintf(output1, "@fPE_%d_%d-%d/1\n", np, sbgn, sbgn+slen);
      fprintf(output1, "%s\n", s1);
      fprintf(output1, "+\n");
      fprintf(output1, "%s\n", q1);

      fprintf(output2, "@fPE_%d_%d-%d/2\n", np, sbgn, sbgn+slen);
      fprintf(output2, "%s\n", s2);
      fprintf(output2, "+\n");
      fprintf(output2, "%s\n", q2);

    } else {
      //  Successfully washed away non-biotin marked sequences, make MP
      //    Shift the fragment by randomGaussian() with mean slen / 2

      int32 shift = randomGaussian(slen / 2, slen / (2 * mpJunction));

      if (shift < 0)
        shift = 0;

      if (shift > slen)
        shift = slen;

      //  Put 'shift' bases from the end of the insert on the start of sh[],
      //  and then fill the remaining of sh[] with the beginning of the insert.
      //
      //  sh[] == [------>END] [BGN--------->]

      int32  pInsert = bgn;
      int32  pShift  = shift;

      while (pShift < slen)
        //  Copy BGN--->
        sh[pShift++] = seq[pInsert++];

      pInsert = bgn + len - shift;
      pShift  = 0;

      while (pInsert < bgn + len)
        //  Copy --->END
        sh[pShift++] = seq[pInsert++];

      assert(pShift == shift);

      sh[slen] = 0;

      makeSequences(sh, slen, readLen, s1, q1, s2, q2);

      //  Add a marker for the chimeric point.  This unfortunately includes some
      //  knowledge of makeSequences(); the second sequence is reverse complemented.
      //
      //  In r2, the junction is at position 'shift - (slen - readLen)', but the read
      //  is reverse complemented, and then the junction covers that base and the one
      //  previous.
      //
      if ((shift > 0) && (shift < readLen)) {
        q1[shift-1] = QV_BASE + 10;
        q1[shift-0] = QV_BASE + 10;
      }
      if ((shift > slen - readLen) && (shift < slen)) {
        q2[readLen - (shift + readLen - slen) - 1] = QV_BASE + 10;
        q2[readLen - (shift + readLen - slen) - 0] = QV_BASE + 10;

        assert((readLen - (shift + readLen - slen)) > 0);
        assert((readLen - (shift + readLen - slen)) < readLen);
      }

      char  type;

      type = 't';
      if (shift  < readLen)         type = 'a';
      if (shift == 0)               type = 'A';
      if (shift  > slen - readLen)  type = 'b';
      if (shift == slen)            type = 'B';

      fprintf(output, "@%cMP_%d_%d-%d_%d/%d/%d/1\n", type, np, bgn, bgn+len, shift, slen, bgn+len-shift);
      fprintf(output, "%s\n", s1);
      fprintf(output, "+\n");
      fprintf(output, "%s\n", q1);

      fprintf(output, "@%cMP_%d_%d-%d_%d/%d/%d/2\n", type, np, bgn, bgn+len, shift, slen, bgn+len-shift);
      fprintf(output, "%s\n", s2);
      fprintf(output, "+\n");
      fprintf(output, "%s\n", q2);

      fprintf(output1, "@%cMP_%d_%d-%d_%d/%d/%d/1\n", type, np, bgn, bgn+len, shift, slen, bgn+len-shift);
      fprintf(output1, "%s\n", s1);
      fprintf(output1, "+\n");
      fprintf(output1, "%s\n", q1);

      fprintf(output2, "@%cMP_%d_%d-%d_%d/%d/%d/2\n", type, np, bgn, bgn+len, shift, slen, bgn+len-shift);
      fprintf(output2, "%s\n", s2);
      fprintf(output2, "+\n");
      fprintf(output2, "%s\n", q2);
    }
  }
}





int
main(int argc, char **argv) {
  char      *fastaName = NULL;
  FILE      *fastaFile = NULL;

  int32      seqMax = 0;
  int32      seqLen = 0;
  char      *seq    = NULL;

  int32      readLen        = 0;    //  Length of read to generate
  int32      readPairs      = 0;    //  Number of pairs to generate, constant
  double     readCoverage   = 0.0;  //  Number of pairs to generate, based on length of sequence

  bool       peEnable       = false;
  int32      peShearSize    = 0;
  int32      peShearStdDev  = 0;

  bool       mpEnable       = false;
  int32      mpInsertSize   = 0;
  int32      mpInsertStdDev = 0;
  int32      mpShearSize    = 0;
  int32      mpShearStdDev  = 0;
  double     mpEnrichment   = 1.0;   //  success rate of washing away paired-end fragments
  double     mpJunction     = 3.0;   //  normally distributed junction location

  char      *outputPrefix   = NULL;
  char       outputName[FILENAME_MAX];
  FILE      *output         = NULL;
  FILE      *output1        = NULL;
  FILE      *output2        = NULL;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if       ((strcmp(argv[arg], "-f") == 0) ||
              (strcmp(argv[arg], "-F") == 0)) {
      fastaName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-l") == 0) {
      readLen = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-n") == 0) {
      readPairs = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-x") == 0) {
      readCoverage = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-pe") == 0) {
      if (arg + 2 >= argc) {
        fprintf(stderr, "Not enough args to -pe.\n");
        err++;
      } else {
        peEnable      = true;
        peShearSize   = atoi(argv[++arg]);
        peShearStdDev = atoi(argv[++arg]);
      }

    } else if (strcmp(argv[arg], "-mp") == 0) {
      if (arg + 6 >= argc) {
        fprintf(stderr, "Not enough args to -pe.\n");
        err++;
      } else {
        mpEnable       = true;
        mpInsertSize   = atoi(argv[++arg]);
        mpInsertStdDev = atoi(argv[++arg]);
        mpShearSize    = atoi(argv[++arg]);
        mpShearStdDev  = atoi(argv[++arg]);
        mpEnrichment   = atof(argv[++arg]);
        mpJunction     = atof(argv[++arg]);
      }
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) ||
      (fastaName == NULL) ||
      (outputPrefix == NULL) ||
      ((peEnable == false) &&
       (mpEnable == false))) {
    fprintf(stderr, "usage: %s -f reference.fasta -o output-prefix -l read-length ....\n", argv[0]);
    fprintf(stderr, "  -f ref.fasta    Use sequences in ref.fasta as the genome.\n");
    fprintf(stderr, "  -o name         Create outputs name.1.fastq and name.2.fastq (and maybe others).\n");
    fprintf(stderr, "  -l len          Create reads of length 'len' bases.\n");
    fprintf(stderr, "  -n np           Create 'np' pairs of reads.\n");
    fprintf(stderr, "  -x cov          Set 'np' to create reads that sample the genome to 'cov' coverage.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -pe shearSize shearStdDev\n");
    fprintf(stderr, "                  Create paired-end reads, from fragments of size 'shearSize +- shearStdDev'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -mp insertSize insertStdDev shearSize shearStdDev enrichment junction\n");
    fprintf(stderr, "                  Create mate-pair reads.  The pairs will be 'insertSize +- insertStdDev'\n");
    fprintf(stderr, "                  apart.  The circularized insert is then sheared into fragments of size\n");
    fprintf(stderr, "                  'shearSize +- shearStdDev'.  With probability 'enrichment' the fragment\n");
    fprintf(stderr, "                  containing the junction is used to form the pair of reads.  The junction\n");
    fprintf(stderr, "                  location is normally distributed through this fragment, with mean 'shearSize/2'\n");
    fprintf(stderr, "                  and std.dev 'shearSize/2/junction'.  With a 500bp fragment, and 100bp reads,\n");
    fprintf(stderr, "                  junction=3 will give about 6%% chimeric reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output QV's are the Sanger spec.\n");
    fprintf(stderr, "\n");

    if (fastaName == NULL)
      fprintf(stderr, "ERROR:  No fasta file (-f) supplied.\n");
    if (outputPrefix == NULL)
      fprintf(stderr, "ERROR:  No output prefix (-o) supplied.\n");
    if ((peEnable == false) && (mpEnable == false))
      fprintf(stderr, "ERROR:  No type (-pe or -mp) selected.\n");

    exit(1);
  }

  //
  //  Initialize
  //

  srand48(1);

  memset(reverseComplement, '&', sizeof(char) * 256);

  reverseComplement['A'] = 'T';
  reverseComplement['C'] = 'G';
  reverseComplement['G'] = 'C';
  reverseComplement['T'] = 'A';

  memset(errorBase, '*', sizeof(char) * 256 * 3);

  errorBase['A'][0] = 'C';  errorBase['A'][1] = 'G';  errorBase['A'][2] = 'T';
  errorBase['C'][0] = 'A';  errorBase['C'][1] = 'G';  errorBase['C'][2] = 'T';
  errorBase['G'][0] = 'A';  errorBase['G'][1] = 'C';  errorBase['G'][2] = 'T';
  errorBase['T'][0] = 'A';  errorBase['T'][1] = 'C';  errorBase['T'][2] = 'G';


  //
  //  Open output files, failing quickly.
  //

  errno = 0;

  sprintf(outputName, "%s.fastq", outputPrefix);
  output = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", outputName, strerror(errno)), exit(1);

  sprintf(outputName, "%s.1.fastq", outputPrefix);
  output1 = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", outputName, strerror(errno)), exit(1);

  sprintf(outputName, "%s.2.fastq", outputPrefix);
  output2 = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", outputName, strerror(errno)), exit(1);

  //
  //  Load all reference sequences into a single string.  Seperate different sequences with a '>', we'll not make
  //  fragments that span these markers (inefficiently, sigh).
  //

  fastaFile = fopen(fastaName, "r");
  if (errno)
    fprintf(stderr, "Failed to open fasta file '%s': %s\n", fastaName, strerror(errno)), exit(1);

  seqMax = AS_UTL_sizeOfFile(fastaName);
  seqLen = 0;
  seq    = new char [seqMax + 1];

  memset(seq, '.', sizeof(char) * seqMax);


  while (!feof(fastaFile)) {
    fgets(seq + seqLen, seqMax - seqLen, fastaFile);

    if (seq[seqLen] == '>')
      continue;

    for (;
         ((seq[seqLen] != '\n') && (seq[seqLen] != '\r') && (seq[seqLen] != 0));
         seqLen++)
      seq[seqLen] = toupper(seq[seqLen]);

    seq[seqLen++] = '>';

    assert(seqLen < seqMax);
  }

  fclose(fastaFile);

  seq[seqLen--] = 0;  //  Get rid of that last '>'

  fprintf(stderr, "READ sequence of length %d\n", seqLen);

  //
  //  If requested, compute the number of pairs to get a desired X of coverage
  //

  if (readCoverage > 0) {
    readPairs  = readCoverage * seqLen / readLen;
    readPairs /= 2;

    fprintf(stderr, "For %.2f X coverage of a %dbp genome, generate %d pairs of %dbp reads.\n",
            readCoverage, seqLen, readPairs, readLen);
  }

  //
  //
  //

  if (peEnable)
    makePE(seq, seqLen, output, output1, output2, readLen, readPairs, peShearSize, peShearStdDev);

  if (mpEnable)
    makeMP(seq, seqLen, output, output1, output2, readLen, readPairs, mpInsertSize, mpInsertStdDev, mpShearSize, mpShearStdDev, mpEnrichment, mpJunction);

  //
  //
  //

  fclose(output);
  fclose(output1);
  fclose(output2);

  delete [] seq;

  exit(0);
}
