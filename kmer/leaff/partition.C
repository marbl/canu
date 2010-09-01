#include "bio++.H"
#include "seqCache.H"

#include <math.h>

struct partition_s {
  u32bit  length;
  u32bit  index;
  u32bit  partition;
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
  u32bit        n = F->getNumberOfSequences();
  partition_s  *p = new partition_s [n];

  for (u32bit i=0; i<n; i++) {
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
                partition_s *p, u32bit openP, u32bit n) {
  char  filename[1024];

  //  Check that everything has been partitioned
  //
  for (u32bit i=0; i<n; i++)
    if (p[i].partition == 0)
      fprintf(stderr, "ERROR: Failed to partition "u32bitFMT"\n", i);

  if (prefix) {

    //  This rewrites the source fasta file into partitioned fasta files
    //
    for (u32bit o=1; o<=openP; o++) {
      sprintf(filename, "%s-"u32bitFMTW(03)".fasta", prefix, o);

      errno = 0;
      FILE *file = fopen(filename, "w");
      if (errno)
        fprintf(stderr, "Couldn't open '%s' for write: %s\n", filename, strerror(errno));

      for (u32bit i=0; i<n; i++)
        if (p[i].partition == o) {
          seqInCore *S = F->getSequenceInCore(p[i].index);
          fprintf(file, ">%s\n", S->header());
          fwrite(S->sequence(), sizeof(char), S->sequenceLength(), file);
          fprintf(file, "\n");

          if (S->sequenceLength() != p[i].length) {
            fprintf(stderr, "Huh?  '%s' "u32bitFMT" != "u32bitFMT"\n", S->header(), S->sequenceLength(), p[i].length);
          }

          delete S;
        }

      fclose(file);
    }

  } else {

    //  This dumps the partition information to stdout.
    //
    fprintf(stdout, u32bitFMT"\n", openP);
    for (u32bit o=1; o<=openP; o++) {
      u32bit  sizeP = 0;
      for (u32bit i=0; i<n; i++)
        if (p[i].partition == o)
          sizeP += p[i].length;
      fprintf(stdout, u32bitFMT"]("u32bitFMT")", o, sizeP);
      for (u32bit i=0; i<n; i++)
        if (p[i].partition == o)
          fprintf(stdout, " "u32bitFMT"("u32bitFMT")", p[i].index, p[i].length);
      fprintf(stdout, "\n");
    }

  }
}


void
partitionBySize(char *prefix, u64bit partitionSize, char *filename) {
  seqCache     *F = new seqCache(filename);
  u32bit        n = F->getNumberOfSequences();
  partition_s  *p = loadPartition(F);

  u32bit  openP = 1;  //  Currently open partition
  u32bit  sizeP = 0;  //  Size of open partition
  u32bit  seqsP = n;  //  Number of sequences to partition

  //  For any sequences larger than partitionSize, create
  //  partitions containing just one sequence
  //
  for (u32bit i=0; i<n; i++) {
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
    for (u32bit i=0; i<n; i++) {
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
partitionByBucket(char *prefix, u64bit partitionSize, char *filename) {
  seqCache     *F = new seqCache(filename);
  u32bit        n = F->getNumberOfSequences();
  partition_s  *p = loadPartition(F);

  if (partitionSize > n)
    partitionSize = n;

  //  The size, in bases, of each partition
  //
  u32bit       *s = new u32bit [partitionSize];
  for (u32bit i=0; i<partitionSize; i++)
    s[i] = 0;

  //  For each sequence
  //
  for (u32bit nextS=0; nextS<n; nextS++) {

    //  find the smallest partition
    //
    u32bit openP = 0;
    for (u32bit i=0; i<partitionSize; i++)
      if (s[i] < s[openP])
        openP = i;

    //  add the next largest sequence to the open partition
    //
    s[openP] += p[nextS].length;
    p[nextS].partition = openP+1;
  }

  outputPartition(F, prefix, p, (u32bit)partitionSize, n);

  delete [] p;
  delete    F;
}


void
partitionBySegment(char *prefix, u64bit numSegments, char *filename) {
  seqCache     *F = new seqCache(filename);
  u32bit        n = F->getNumberOfSequences();
  partition_s  *p = new partition_s [n];
  u32bit        numSeqPerPart = (u32bit)ceil(n / (double)numSegments);

  for (u32bit i=0; i<n; i++) {
    p[i].length    = F->getSequenceLength(i);
    p[i].index     = i;
    p[i].partition = i / numSeqPerPart + 1;
  }

  outputPartition(F, prefix, p, numSegments, n);

  delete [] p;
  delete    F;
}
