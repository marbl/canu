
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

const char *mainid = "$Id: fastqSimulate.C,v 1.14 2011-08-30 02:57:17 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_reverseComplement.h"

#include <vector>
using namespace std;

vector<int32> seqStartPositions;

static char revComp[256];
static char errorBase[256][3];
static char validBase[256];

double readErrorRate = 0.01;  //  Fraction error
bool   allowGaps     = false;

const uint32 mpJunctionsNone   = 0;
const uint32 mpJunctionsNormal = 1;
const uint32 mpJunctionsAlways = 2;


#define QV_BASE  '!'


//  Returns random int in range bgn <= x < end.
//
int32
randomUniform(int32 bgn, int32 end) {
  if (bgn >= end)
    fprintf(stderr, "randomUniform()-- ERROR:  invalid range bgn=%d end=%d\n", bgn, end);
  assert(bgn < end);
  return((int32)floor((end - bgn) * drand48() + bgn));
}



//  Generate a random gaussian using the Marsaglia polar method.
//
int32
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

  return((int32)(mean + u * r * stddev));
}




int32
findSequenceIndex(int32 pos) {
  int32  seqIdx = 0;

  for (seqIdx=0; seqIdx<seqStartPositions.size(); seqIdx++)
    if (pos < seqStartPositions[seqIdx])
      break;

  assert(seqIdx > 0);
  seqIdx--;

  return(seqIdx);
}



void
makeSequences(char    *frag,
              int32    fragLen,
              int32    readLen,
              char    *s1,
              char    *q1,
              char    *s2,
              char    *q2) {

  for (int32 p=0, i=0; p<readLen; p++, i++) {
    s1[p] = frag[i];
    q1[p] = (validBase[s1[p]]) ? QV_BASE + 39 : QV_BASE + 2;

    if (drand48() < readErrorRate) {
      s1[p] = errorBase[s1[p]][randomUniform(0, 3)];
      q1[p] = (validBase[s1[p]]) ? QV_BASE + 8 : QV_BASE + 2;
    }

    assert(s1[p] != '*');
  }

  s1[readLen] = 0;
  q1[readLen] = 0;

  if ((fragLen == 0) || (s2 == NULL) || (q2 == NULL))
    return;

  for (int32 p=0, i=fragLen-1; p<readLen; p++, i--) {
    s2[p] = revComp[frag[i]];
    q2[p] = (validBase[s2[p]]) ? QV_BASE + 39 : QV_BASE + 2;

    if (drand48() < readErrorRate) {
      s2[p] = errorBase[s2[p]][randomUniform(0, 3)];
      q2[p] = (validBase[s2[p]]) ? QV_BASE + 8 : QV_BASE + 2;
    }

    assert(s2[p] != '*');
  }  

  s2[readLen] = 0;
  q2[readLen] = 0;
}




void
makeSE(char   *seq,
       int32   seqLen,
       FILE   *outputI,
       FILE   *outputC,
       int32   readLen,
       int32   numReads) {
  char   *s1 = new char [readLen + 1];
  char   *q1 = new char [readLen + 1];

  for (int32 nr=0; nr<numReads; nr++) {
  trySEagain:
    int32   len = readLen;
    int32   bgn = randomUniform(1, seqLen - len);
    int32   idx = findSequenceIndex(bgn);
    int32   zer = seqStartPositions[idx];

    //  Scan the sequence, if we spanned a sequence break or encounter a block of Ns, don't use this pair

    for (int32 i=bgn; i<bgn+len; i++)
      if ((seq[i] == '>') ||
          ((allowGaps == false) && (seq[i] == 'N')))
        goto trySEagain;

    //  Generate the sequence.

    makeSequences(seq + bgn, 0, readLen, s1, q1, NULL, NULL);

    //  Make sure the read doesn't contain N's (redundant in this particular case)

    for (int32 i=0; i<readLen; i++)
      if (s1[i] == 'N')
        goto trySEagain;

    //  Output sequence, with a descriptive ID

    fprintf(outputI, "@SE_%d_%d@%d-%d/1\n", nr, idx, bgn-zer, bgn+len-zer);
    fprintf(outputI, "%s\n", s1);
    fprintf(outputI, "+\n");
    fprintf(outputI, "%s\n", q1);
  }

  delete [] s1;
  delete [] q1;
}


void
makePE(char   *seq,
       int32   seqLen,
       FILE   *outputI,
       FILE   *outputC,
       FILE   *output1,
       FILE   *output2,
       int32   readLen,
       int32   numPairs,
       int32   peShearSize,
       int32   peShearStdDev) {
  char   *s1 = new char [readLen + 1];
  char   *q1 = new char [readLen + 1];
  char   *s2 = new char [readLen + 1];
  char   *q2 = new char [readLen + 1];

  for (int32 np=0; np<numPairs; np++) {
  tryPEagain:
    int32   len = randomGaussian(peShearSize, peShearStdDev);
    int32   bgn = randomUniform(1, seqLen - len);
    int32   idx = findSequenceIndex(bgn);
    int32   zer = seqStartPositions[idx];

    if (len <= readLen)
      goto tryPEagain;

    //  Scan the sequence, if we spanned a sequence break, don't use this pair

    for (int32 i=bgn; i<bgn+len; i++)
      if ((seq[i] == '>') ||
          ((allowGaps == false) && (seq[i] == 'N')))
        goto tryPEagain;


    //  Read sequences from the ends.

    makeSequences(seq + bgn, len, readLen, s1, q1, s2, q2);

    //  Make sure the reads don't contain N's

    for (int32 i=0; i<readLen; i++)
      if ((s1[i] == 'N') || (s2[i] == 'N'))
        goto tryPEagain;

    //  Output sequences, with a descriptive ID

    fprintf(outputI, "@PE_%d_%d@%d-%d/1\n", np, idx, bgn-zer, bgn+len-zer);
    fprintf(outputI, "%s\n", s1);
    fprintf(outputI, "+\n");
    fprintf(outputI, "%s\n", q1);

    fprintf(outputI, "@PE_%d_%d@%d-%d/2\n", np, idx, bgn-zer, bgn+len-zer);
    fprintf(outputI, "%s\n", s2);
    fprintf(outputI, "+\n");
    fprintf(outputI, "%s\n", q2);

    fprintf(output1, "@PE_%d_%d@%d-%d/1\n", np, idx, bgn-zer, bgn+len-zer);
    fprintf(output1, "%s\n", s1);
    fprintf(output1, "+\n");
    fprintf(output1, "%s\n", q1);

    fprintf(output2, "@PE_%d_%d@%d-%d/2\n", np, idx, bgn-zer, bgn+len-zer);
    fprintf(output2, "%s\n", s2);
    fprintf(output2, "+\n");
    fprintf(output2, "%s\n", q2);

    reverseComplement(s1, q1, readLen);
    reverseComplement(s2, q2, readLen);

    fprintf(outputC, "@PE_%d_%d@%d-%d/1\n", np, idx, bgn+len-zer, bgn-zer);
    fprintf(outputC, "%s\n", s1);
    fprintf(outputC, "+\n");
    fprintf(outputC, "%s\n", q1);

    fprintf(outputC, "@PE_%d_%d@%d-%d/2\n", np, idx, bgn+len-zer, bgn-zer);
    fprintf(outputC, "%s\n", s2);
    fprintf(outputC, "+\n");
    fprintf(outputC, "%s\n", q2);
  }

  delete [] s1;
  delete [] q1;
  delete [] s2;
  delete [] q2;
}





void
makeMP(char   *seq,
       int32   seqLen,
       FILE   *outputI,
       FILE   *outputC,
       FILE   *output1,
       FILE   *output2,
       int32   readLen,
       int32   numPairs,
       int32   mpInsertSize,
       int32   mpInsertStdDev,
       int32   mpShearSize,
       int32   mpShearStdDev,
       double  mpEnrichment,
       uint32  mpJunctions) {
  char   *s1 = new char [readLen + 1];
  char   *q1 = new char [readLen + 1];
  char   *s2 = new char [readLen + 1];
  char   *q2 = new char [readLen + 1];
  char   *sh = new char [1048576];

  for (int32 np=0; np<numPairs; np++) {
  tryMPagain:
    int32   len = randomGaussian(mpInsertSize, mpInsertStdDev);
    int32   bgn = randomUniform(1, seqLen - len);
    int32   idx = findSequenceIndex(bgn);
    int32   zer = seqStartPositions[idx];

    int32   slen = randomGaussian(mpShearSize, mpShearStdDev);  //  shear size

    if ((len  <= readLen) ||
        (slen <= readLen))
      goto tryMPagain;

    //  Scan the sequence, if we spanned a sequence break, don't use this pair

    for (int32 i=bgn; i<bgn+len; i++)
      if ((seq[i] == '>') ||
          ((allowGaps == false) && (seq[i] == 'N')))
        goto tryMPagain;


    //  If we fail the mpEnrichment test, pick a random shearing and return PE reads.
    //  Otherwise, rotate the sequence to circularize and return MP reads.

    if (mpEnrichment < drand48()) {
      //  Failed to wash away non-biotin marked sequence, make PE
      int32  sbgn = bgn + randomUniform(0, len - slen);

      makeSequences(seq + sbgn, slen, readLen, s1, q1, s2, q2);

      //  Make sure the reads don't contain N's

      for (int32 i=0; i<readLen; i++)
        if ((s1[i] == 'N') || (s2[i] == 'N'))
          goto tryMPagain;

      //  Output sequences, with a descriptive ID

      fprintf(outputI, "@fPE_%d_%d@%d-%d/1\n", np, idx, sbgn-zer, sbgn+slen-zer);
      fprintf(outputI, "%s\n", s1);
      fprintf(outputI, "+\n");
      fprintf(outputI, "%s\n", q1);

      fprintf(outputI, "@fPE_%d_%d@%d-%d/2\n", np, idx, sbgn-zer, sbgn+slen-zer);
      fprintf(outputI, "%s\n", s2);
      fprintf(outputI, "+\n");
      fprintf(outputI, "%s\n", q2);

      fprintf(output1, "@fPE_%d_%d@%d-%d/1\n", np, idx, sbgn-zer, sbgn+slen-zer);
      fprintf(output1, "%s\n", s1);
      fprintf(output1, "+\n");
      fprintf(output1, "%s\n", q1);

      fprintf(output2, "@fPE_%d_%d@%d-%d/2\n", np, idx, sbgn-zer, sbgn+slen-zer);
      fprintf(output2, "%s\n", s2);
      fprintf(output2, "+\n");
      fprintf(output2, "%s\n", q2);

      reverseComplement(s1, q1, readLen);
      reverseComplement(s2, q2, readLen);

      fprintf(outputC, "@fPE_%d_%d@%d-%d/1\n", np, idx, sbgn+slen-zer, sbgn-zer);
      fprintf(outputC, "%s\n", s1);
      fprintf(outputC, "+\n");
      fprintf(outputC, "%s\n", q1);

      fprintf(outputC, "@fPE_%d_%d@%d-%d/2\n", np, idx, sbgn+slen-zer, sbgn-zer);
      fprintf(outputC, "%s\n", s2);
      fprintf(outputC, "+\n");
      fprintf(outputC, "%s\n", q2);

    } else {
      //  Successfully washed away non-biotin marked sequences, make MP.  Shift the fragment by a
      //  random amount.  If we are not allowed to make junction reads, the shift must allow at
      //  least a read-length of sequence on each end.
      //
      //  If the shear size is less than two read lengths there is no way we can make a perfect MP
      //  pair.  We allow this case for normal junctions (both reads might be junction reads) but
      //  disallow it for the other two cases (no junction reads and all junction reads).

      int32 shift = 0;

      if (mpJunctions == mpJunctionsNormal) {
        shift = randomUniform(0, slen);

      } else if (mpJunctions == mpJunctionsNone) {
        if (slen <= 2 * readLen)
          goto tryMPagain;

        shift = randomUniform(readLen, slen - readLen);

      } else if (mpJunctions == mpJunctionsAlways) {
        if (slen <= 2 * readLen)
          goto tryMPagain;

        if (randomUniform(0, 100) < 50)
          shift = randomUniform(0, readLen);
        else
          shift = randomUniform(slen - readLen, slen);
      }

      if ((shift < 0) || (shift > slen))
        fprintf(stderr, "ERROR:  invalid shift %d.\n", shift);
      assert(shift >= 0);
      assert(shift <= slen);


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

      //  Make sure the reads don't contain N's

      for (int32 i=0; i<readLen; i++)
        if ((s1[i] == 'N') || (s2[i] == 'N'))
          goto tryMPagain;

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
      if (shift >= slen - readLen)  type = 'b';

      if (mpJunctions == mpJunctionsNone)
        assert(type == 't');
      if (mpJunctions == mpJunctionsAlways)
        assert(type != 't');

      fprintf(outputI, "@%cMP_%d_%d@%d-%d_%d/%d/%d/1\n", type, np, idx, bgn, bgn+len, shift, slen, bgn+len-shift);
      fprintf(outputI, "%s\n", s1);
      fprintf(outputI, "+\n");
      fprintf(outputI, "%s\n", q1);

      fprintf(outputI, "@%cMP_%d_%d@%d-%d_%d/%d/%d/2\n", type, np, idx, bgn, bgn+len, shift, slen, bgn+len-shift);
      fprintf(outputI, "%s\n", s2);
      fprintf(outputI, "+\n");
      fprintf(outputI, "%s\n", q2);

      fprintf(output1, "@%cMP_%d_%d@%d-%d_%d/%d/%d/1\n", type, np, idx, bgn, bgn+len, shift, slen, bgn+len-shift);
      fprintf(output1, "%s\n", s1);
      fprintf(output1, "+\n");
      fprintf(output1, "%s\n", q1);

      fprintf(output2, "@%cMP_%d_%d@%d-%d_%d/%d/%d/2\n", type, np, idx, bgn, bgn+len, shift, slen, bgn+len-shift);
      fprintf(output2, "%s\n", s2);
      fprintf(output2, "+\n");
      fprintf(output2, "%s\n", q2);

      reverseComplement(s1, q1, readLen);
      reverseComplement(s2, q2, readLen);

      fprintf(outputC, "@%cMP_%d_%d@%d-%d_%d/%d/%d/1\n", type, np, idx, bgn+len, bgn, shift, slen, bgn+len-shift);
      fprintf(outputC, "%s\n", s1);
      fprintf(outputC, "+\n");
      fprintf(outputC, "%s\n", q1);

      fprintf(outputC, "@%cMP_%d_%d@%d-%d_%d/%d/%d/2\n", type, np, idx, bgn+len, bgn, shift, slen, bgn+len-shift);
      fprintf(outputC, "%s\n", s2);
      fprintf(outputC, "+\n");
      fprintf(outputC, "%s\n", q2);
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
  int32      numReads       = 0;    //  Number of reads to generate, constant
  int32      numPairs       = 0;    //  Number of pairs to generate, constant (= numReads / 2)
  double     readCoverage   = 0.0;  //  Number of pairs to generate, based on length of sequence

  bool       seEnable       = false;

  bool       peEnable       = false;
  int32      peShearSize    = 0;
  int32      peShearStdDev  = 0;

  bool       mpEnable       = false;
  int32      mpInsertSize   = 0;
  int32      mpInsertStdDev = 0;
  int32      mpShearSize    = 0;
  int32      mpShearStdDev  = 0;
  double     mpEnrichment   = 1.0;   //  success rate of washing away paired-end fragments
  uint32     mpJunctions    = mpJunctionsNormal;

  char      *outputPrefix   = NULL;
  char       outputName[FILENAME_MAX];
  FILE      *outputI        = NULL;  //  Interleaved output
  FILE      *outputC        = NULL;  //  Interleaved output, reverse complemented
  FILE      *output1        = NULL;  //  A read output
  FILE      *output2        = NULL;  //  B read output

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
      numReads = numPairs = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-x") == 0) {
      readCoverage = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      readErrorRate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-allowgaps") == 0) {
      allowGaps = true;

    } else if (strcmp(argv[arg], "-nojunction") == 0) {
      mpJunctions = mpJunctionsNone;

    } else if (strcmp(argv[arg], "-alljunction") == 0) {
      mpJunctions = mpJunctionsAlways;

    } else if (strcmp(argv[arg], "-se") == 0) {
      seEnable = true;

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
      if (arg + 5 >= argc) {
        fprintf(stderr, "Not enough args to -mp.\n");
        err++;
      } else {
        mpEnable       = true;
        mpInsertSize   = atoi(argv[++arg]);
        mpInsertStdDev = atoi(argv[++arg]);
        mpShearSize    = atoi(argv[++arg]);
        mpShearStdDev  = atoi(argv[++arg]);
        mpEnrichment   = atof(argv[++arg]);
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
      ((seEnable == false) &&
       (peEnable == false) &&
       (mpEnable == false))) {
    fprintf(stderr, "usage: %s -f reference.fasta -o output-prefix -l read-length ....\n", argv[0]);
    fprintf(stderr, "  -f ref.fasta    Use sequences in ref.fasta as the genome.\n");
    fprintf(stderr, "  -o name         Create outputs name.1.fastq and name.2.fastq (and maybe others).\n");
    fprintf(stderr, "  -l len          Create reads of length 'len' bases.\n");
    fprintf(stderr, "  -n n            Create 'n' reads (for -se) or 'n' pairs of reads (for -pe and -mp).\n");
    fprintf(stderr, "  -x cov          Set 'np' to create reads that sample the genome to 'cov' coverage.\n");
    fprintf(stderr, "  -e err          Reads will contain fraction error 'e' (0.01 == 1%% error).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -allowgaps      Allow pairs to span N regions in the reference.  By default, pairs\n");
    fprintf(stderr, "                  are not allowed to span a gap.  Reads are never allowed to cover N's.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nojunction     For -mp, do not create chimeric junction reads.  Create only fully PE or\n");
    fprintf(stderr, "                  fully MP reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -se\n");
    fprintf(stderr, "                  Create single-end reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -pe shearSize shearStdDev\n");
    fprintf(stderr, "                  Create paired-end reads, from fragments of size 'shearSize +- shearStdDev'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -mp insertSize insertStdDev shearSize shearStdDev enrichment\n");
    fprintf(stderr, "                  Create mate-pair reads.  The pairs will be 'insertSize +- insertStdDev'\n");
    fprintf(stderr, "                  apart.  The circularized insert is then sheared into fragments of size\n");
    fprintf(stderr, "                  'shearSize +- shearStdDev'.  With probability 'enrichment' the fragment\n");
    fprintf(stderr, "                  containing the junction is used to form the pair of reads.  The junction\n");
    fprintf(stderr, "                  location is uniformly distributed through this fragment.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output QV's are the Sanger spec.\n");
    fprintf(stderr, "\n");

    if (fastaName == NULL)
      fprintf(stderr, "ERROR:  No fasta file (-f) supplied.\n");
    if (outputPrefix == NULL)
      fprintf(stderr, "ERROR:  No output prefix (-o) supplied.\n");
    if ((seEnable == false) && (peEnable == false) && (mpEnable == false))
      fprintf(stderr, "ERROR:  No type (-se or -pe or -mp) selected.\n");

    exit(1);
  }

  //
  //  Initialize
  //

  srand48(time(NULL));

  memset(revComp, '&', sizeof(char) * 256);

  revComp['A'] = 'T';
  revComp['C'] = 'G';
  revComp['G'] = 'C';
  revComp['T'] = 'A';
  revComp['N'] = 'N';

  memset(errorBase, '*', sizeof(char) * 256 * 3);

  errorBase['A'][0] = 'C';  errorBase['A'][1] = 'G';  errorBase['A'][2] = 'T';
  errorBase['C'][0] = 'A';  errorBase['C'][1] = 'G';  errorBase['C'][2] = 'T';
  errorBase['G'][0] = 'A';  errorBase['G'][1] = 'C';  errorBase['G'][2] = 'T';
  errorBase['T'][0] = 'A';  errorBase['T'][1] = 'C';  errorBase['T'][2] = 'G';
  errorBase['N'][0] = 'N';  errorBase['N'][1] = 'N';  errorBase['N'][2] = 'N';

  memset(validBase, 0, sizeof(char) * 256);

  validBase['A'] = 1;
  validBase['C'] = 1;
  validBase['G'] = 1;
  validBase['T'] = 1;


  //
  //  Open output files, failing quickly.
  //

  errno = 0;

  if (seEnable == true) {
    sprintf(outputName, "%s.s.fastq", outputPrefix);
    outputI = fopen(outputName, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", outputName, strerror(errno)), exit(1);
  }

  if (seEnable == false) {
    sprintf(outputName, "%s.i.fastq", outputPrefix);
    outputI = fopen(outputName, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", outputName, strerror(errno)), exit(1);

    sprintf(outputName, "%s.c.fastq", outputPrefix);
    outputC = fopen(outputName, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", outputName, strerror(errno)), exit(1);
  }

  if (peEnable || mpEnable) {
    sprintf(outputName, "%s.1.fastq", outputPrefix);
    output1 = fopen(outputName, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", outputName, strerror(errno)), exit(1);

    sprintf(outputName, "%s.2.fastq", outputPrefix);
    output2 = fopen(outputName, "w");
    if (errno)
      fprintf(stderr, "Failed to open output file '%s': %s\n", outputName, strerror(errno)), exit(1);
  }

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

    if (seq[seqLen] == '>') {
      seqLen++;
      seqStartPositions.push_back(seqLen);
      continue;
    }

    for (;
         ((seq[seqLen] != '\n') && (seq[seqLen] != '\r') && (seq[seqLen] != 0));
         seqLen++)
      seq[seqLen] = toupper(seq[seqLen]);

    assert(seqLen < seqMax);
  }

  fclose(fastaFile);

  seq[seqLen] = 0;

  fprintf(stderr, "READ sequence of length %d\n", seqLen);

  //
  //  If requested, compute the number of pairs to get a desired X of coverage
  //

  if (readCoverage > 0) {
    numReads = (int32)floor(readCoverage * seqLen / readLen);
    numPairs = numReads / 2;

    if (seEnable)
      fprintf(stderr, "For %.2f X coverage of a %dbp genome, generate %d %dbp reads.\n",
              readCoverage, seqLen, numReads, readLen);
    else
      fprintf(stderr, "For %.2f X coverage of a %dbp genome, generate %d pairs of %dbp reads.\n",
              readCoverage, seqLen, numPairs, readLen);
  }

  //
  //
  //

  if (seEnable)
    makeSE(seq, seqLen, outputI, outputC, readLen, numReads);

  if (peEnable)
    makePE(seq, seqLen, outputI, outputC, output1, output2, readLen, numPairs, peShearSize, peShearStdDev);

  if (mpEnable)
    makeMP(seq, seqLen, outputI, outputC, output1, output2, readLen, numPairs, mpInsertSize, mpInsertStdDev, mpShearSize, mpShearStdDev, mpEnrichment, mpJunctions);

  //
  //
  //

  fclose(outputI);
  fclose(outputC);

  if (peEnable || mpEnable) {
    fclose(output1);
    fclose(output2);
  }

  delete [] seq;

  exit(0);
}
