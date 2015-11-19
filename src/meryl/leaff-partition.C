
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
 *    kmer/leaff/partition.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2009-FEB-07 to 2014-APR-11
 *      are Copyright 2009-2010,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "seqCache.H"

#include <math.h>

struct partition_s {
  uint32  length;
  uint32  index;
  uint32  partition;
};


static
int
partition_s_compare(const void *A, const void *B) {
  const partition_s *a = (const partition_s *)A;
  const partition_s *b = (const partition_s *)B;
  if (a->length < b->length)
    return(1);
  if (a->length > b->length)
    return(-1);
  return(0);
}


static
partition_s *
loadPartition(seqCache *F) {
  uint32        n = F->getNumberOfSequences();
  partition_s  *p = new partition_s [n];

  for (uint32 i=0; i<n; i++) {
    p[i].length    = F->getSequenceLength(i);
    p[i].index     = i;
    p[i].partition = 0;
  }

  qsort(p, n, sizeof(partition_s), partition_s_compare);

  return(p);
}


static
void
outputPartition(seqCache *F,
                char *prefix,
                partition_s *p, uint32 openP, uint32 n) {
  char  filename[1024];

  //  Check that everything has been partitioned
  //
  for (uint32 i=0; i<n; i++)
    if (p[i].partition == 0)
      fprintf(stderr, "ERROR: Failed to partition "F_U32"\n", i);

  if (prefix) {

    //  This rewrites the source fasta file into partitioned fasta files
    //
    for (uint32 o=1; o<=openP; o++) {
      sprintf(filename, "%s-%03"F_U32P".fasta", prefix, o);

      errno = 0;
      FILE *file = fopen(filename, "w");
      if (errno)
        fprintf(stderr, "Couldn't open '%s' for write: %s\n", filename, strerror(errno));

      for (uint32 i=0; i<n; i++)
        if (p[i].partition == o) {
          seqInCore *S = F->getSequenceInCore(p[i].index);
          fprintf(file, ">%s\n", S->header());
          fwrite(S->sequence(), sizeof(char), S->sequenceLength(), file);
          fprintf(file, "\n");

          if (S->sequenceLength() != p[i].length) {
            fprintf(stderr, "Huh?  '%s' "F_U32" != "F_U32"\n", S->header(), S->sequenceLength(), p[i].length);
          }

          delete S;
        }

      fclose(file);
    }

  } else {

    //  This dumps the partition information to stdout.
    //
    fprintf(stdout, F_U32"\n", openP);
    for (uint32 o=1; o<=openP; o++) {
      uint32  sizeP = 0;
      for (uint32 i=0; i<n; i++)
        if (p[i].partition == o)
          sizeP += p[i].length;
      fprintf(stdout, F_U32"]("F_U32")", o, sizeP);
      for (uint32 i=0; i<n; i++)
        if (p[i].partition == o)
          fprintf(stdout, " "F_U32"("F_U32")", p[i].index, p[i].length);
      fprintf(stdout, "\n");
    }

  }
}


void
partitionBySize(char *prefix, uint64 partitionSize, char *filename) {
  seqCache     *F = new seqCache(filename);
  uint32        n = F->getNumberOfSequences();
  partition_s  *p = loadPartition(F);

  uint32  openP = 1;  //  Currently open partition
  uint32  sizeP = 0;  //  Size of open partition
  uint32  seqsP = n;  //  Number of sequences to partition

  //  For any sequences larger than partitionSize, create
  //  partitions containing just one sequence
  //
  for (uint32 i=0; i<n; i++) {
    if (p[i].length > partitionSize) {
      p[i].partition = openP++;
      seqsP--;
    }
  }

  //  For the remaining, iterate through the list,
  //  greedily placing the longest sequence that fits
  //  into the open partition
  //
  while (seqsP > 0) {
    for (uint32 i=0; i<n; i++) {
      if ((p[i].partition == 0) &&
          (p[i].length + sizeP < partitionSize)) {
        p[i].partition = openP;
        sizeP += p[i].length;
        seqsP--;
      }
    }

    openP++;
    sizeP = 0;
  }

  outputPartition(F, prefix, p, openP-1, n);

  delete [] p;
  delete    F;
}


void
partitionByBucket(char *prefix, uint64 partitionSize, char *filename) {
  seqCache     *F = new seqCache(filename);
  uint32        n = F->getNumberOfSequences();
  partition_s  *p = loadPartition(F);

  if (partitionSize > n)
    partitionSize = n;

  //  The size, in bases, of each partition
  //
  uint32       *s = new uint32 [partitionSize];
  for (uint32 i=0; i<partitionSize; i++)
    s[i] = 0;

  //  For each sequence
  //
  for (uint32 nextS=0; nextS<n; nextS++) {

    //  find the smallest partition
    //
    uint32 openP = 0;
    for (uint32 i=0; i<partitionSize; i++)
      if (s[i] < s[openP])
        openP = i;

    //  add the next largest sequence to the open partition
    //
    s[openP] += p[nextS].length;
    p[nextS].partition = openP+1;
  }

  outputPartition(F, prefix, p, (uint32)partitionSize, n);

  delete [] p;
  delete    F;
}


void
partitionBySegment(char *prefix, uint64 numSegments, char *filename) {
  seqCache     *F = new seqCache(filename);
  uint32        n = F->getNumberOfSequences();
  partition_s  *p = new partition_s [n];
  uint32        numSeqPerPart = (uint32)ceil(n / (double)numSegments);

  for (uint32 i=0; i<n; i++) {
    p[i].length    = F->getSequenceLength(i);
    p[i].index     = i;
    p[i].partition = i / numSeqPerPart + 1;
  }

  outputPartition(F, prefix, p, numSegments, n);

  delete [] p;
  delete    F;
}
