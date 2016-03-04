
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
 *    src/AS_GKP/fastqSimulate-sort.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2013-MAR-21 to 2013-AUG-01
 *      are Copyright 2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2015-JAN-13
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "AS_global.H"
#include "AS_UTL_fileIO.H"
#include "AS_UTL_reverseComplement.H"

#include <vector>
#include <algorithm>
using namespace std;

class pairedRead {
public:

  bool operator<(const pairedRead &that) const {
    return((seqId < that.seqId) ||
           ((seqId == that.seqId) && (seqPos < that.seqPos)));
  }


  uint32  seqId;
  uint32  seqPos;

  char   *readA;
  char   *readB;
};


uint32  hLen = 1024;
uint32  sLen = 1024 * 1024 * 16;

char   *a = NULL;
char   *b = NULL;
char   *c = NULL;
char   *d = NULL;

char *
readRead(FILE *inFile, uint32 &seq, uint32 &bgn, uint32 &end) {

  seq = 0;
  bgn = 0;
  end = 0;

  if (inFile == NULL)
    return(NULL);

  if (a == NULL) {
    a = new char [hLen];
    b = new char [sLen];
    c = new char [hLen];
    d = new char [sLen];
  }

  uint32  al = 0;
  uint32  bl = 0;
  uint32  cl = 0;
  uint32  dl = 0;

  fgets(a, hLen, inFile);  al = strlen(a);
  fgets(b, sLen, inFile);  bl = strlen(b);
  fgets(c, hLen, inFile);  cl = strlen(c);
  fgets(d, sLen, inFile);  dl = strlen(d);

  assert(a[0] == '@');
  assert(c[0] == '+');

  uint32  p=0;

  //  @tMP_0_0@2569105-2572074_239/553/2571835/1

  while (a[p] != '_')
    p++;
  p++;
  while (a[p] != '_')
    p++;
  p++;

  seq = strtoul(a+p, NULL, 10);

  while (a[p] != '@')
    p++;
  p++;

  bgn = strtoul(a+p, NULL, 10);

  while (a[p] != '-')
    p++;
  p++;

  end = strtoul(a+p, NULL, 10);

  //fprintf(stderr, "seq="F_U32" bgn="F_U32" end="F_U32" line %s",
  //        seq, bgn, end, a);

  char *retstr = new char [al + bl + cl + dl + 1];

  memcpy(retstr,                a, al);
  memcpy(retstr + al,           b, bl);
  memcpy(retstr + al + bl,      c, cl);
  memcpy(retstr + al + bl + cl, d, dl);

  retstr[al + bl + cl + dl] = 0;

  return(retstr);
}


int
main(int argc, char **argv) {
  char  *inName1 = NULL;
  char  *inName2 = NULL;
  char  *otName1 = NULL;
  char  *otName2 = NULL;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-i1") == 0) {
      inName1 = argv[++arg];

    } else if (strcmp(argv[arg], "-i2") == 0) {
      inName2 = argv[++arg];

    } else if (strcmp(argv[arg], "-o1") == 0) {
      otName1 = argv[++arg];

    } else if (strcmp(argv[arg], "-o2") == 0) {
      otName2 = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }

  if (inName1 == NULL)
    err++;
  if (otName1 == NULL)
    err++;
  if ((inName2 == NULL) != (otName2 == NULL))
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -i1 in.1.fastq [-i2 in.2.fastq] -o1 out.1.fastq [-o2 out.2.fastq]\n", argv[0]);

    if (inName1 == NULL)
      fprintf(stderr, "ERROR:  No in.1.fastq supplied with -i1.\n");
    if (otName1 == NULL)
      fprintf(stderr, "ERROR:  No out.1.fastq supplied with -i1.\n");
    if ((inName2 == NULL) != (otName2 == NULL))
      fprintf(stderr, "ERROR:  Matedness of input and output don't agree (neither or both -i2 and -o2 must be used).\n");

    exit(1);
  }

  errno = 0;
  FILE *inFile1 = (inName1 == NULL) ? NULL : fopen(inName1, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", inName1, strerror(errno)), exit(1);

  errno = 0;
  FILE *inFile2 = (inName2 == NULL) ? NULL : fopen(inName2, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", inName2, strerror(errno)), exit(1);

  errno = 0;
  FILE *otFile1 = (otName1 == NULL) ? NULL : fopen(otName1, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", otName1, strerror(errno)), exit(1);

  errno = 0;
  FILE *otFile2 = (otName2 == NULL) ? NULL : fopen(otName2, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", otName2, strerror(errno)), exit(1);

  //  Load reads.

  vector<pairedRead>  reads;
  pairedRead          pr;

  while (!feof(inFile1)) {
    uint32  seq1=0, bgn1=0, end1=0;
    uint32  seq2=0, bgn2=0, end2=0;

    char *a = readRead(inFile1, seq1, bgn1, end1);
    char *b = readRead(inFile2, seq2, bgn2, end2);

    if (feof(inFile1))
      break;

    if ((a != NULL) && (b != NULL)) {
      pr.seqId  = (seq1 < seq2) ? seq1 : seq2;
      pr.seqPos = (bgn1 < end1) ? bgn1 : end1;

      assert(bgn1 == bgn2);
      assert(end1 == end2);

    } else if (a != NULL) {
      pr.seqId  = seq1;
      pr.seqPos = bgn1;

    } else {
      assert(0);
    }

    pr.readA = a;
    pr.readB = b;

    reads.push_back(pr);
  }

  fprintf(stderr, "Loaded "F_U64" mated reads.\n", reads.size());

  sort(reads.begin(), reads.end());

  for (uint32 i=0; i<reads.size(); i++) {
    if (reads[i].readA)
      fputs(reads[i].readA, otFile1);
    if (reads[i].readB)
      fputs(reads[i].readB, otFile2);
  }

  if (inFile1)  fclose(inFile1);
  if (inFile2)  fclose(inFile2);
  if (otFile1)  fclose(otFile1);
  if (otFile2)  fclose(otFile2);

  return(0);
}
