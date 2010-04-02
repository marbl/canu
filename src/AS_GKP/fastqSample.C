
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2010, J. Craig Venter Institute.
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

const char *mainid = "$Id: fastqSample.C,v 1.2 2010-04-02 06:14:58 brianwalenz Exp $";

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <errno.h>

#include "AS_global.h"

class aRead {
public:
  aRead() {};
  ~aRead() {};

  void  read(FILE *F) {
    fgets(a, 1024, F);
    fgets(b, 1024, F);
    fgets(c, 1024, F);
    fgets(d, 1024, F);
    if ((a[0] != '@') || (c[0] != '+')) {
      fprintf(stderr, "ERROR:  Not FastQ.  Read lines:\n");
      fprintf(stderr, "  %s", a);
      fprintf(stderr, "  %s", b);
      fprintf(stderr, "  %s", c);
      fprintf(stderr, "  %s", d);
      exit(1);
    }
  };
  void  write(FILE *F) {
    fputs(a, F);
    fputs(b, F);
    fputs(c, F);
    fputs(d, F);
  };

private:
  char  a[1024];
  char  b[1024];
  char  c[1024];
  char  d[1024];
};


int
main(int argc, char **argv) {
  aRead     Ar,        Br;
  FILE     *Ai = NULL, *Bi = NULL;
  FILE     *Ao = NULL, *Bo = NULL;

  uint32   *ids  = NULL;
  char     *save = NULL;

  int32  NUMREADS    = 0;        //  Number of mated reads in the input
  int32  GENOMESIZE  = 0;        //  Size of the genome in bp
  int32  READLENGTH  = 0;        //  For mated reads, 2x read size
  int32  COVERAGE    = 0;        //  Desired coverage in output
  char  *NAME        = NULL;     //  Prefix name

  char   path1[FILENAME_MAX];
  char   path2[FILENAME_MAX];

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-t") == 0) {
      //  Total reads in input
      NUMREADS = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-g") == 0) {
      //  Genome size
      GENOMESIZE = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      //  Length of read
      READLENGTH = atoi(argv[++arg]) * 2;

    } else if (strcmp(argv[arg], "-c") == 0) {
      //  Desired coverage
      COVERAGE = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-n") == 0) {
      NAME = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if (NAME == NULL)
    err++;
  if (GENOMESIZE == 0)
    err++;
  if (READLENGTH == 0)
    err++;
  if (COVERAGE == 0)
    err++;
  if (err) {
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "  -t T     total number of mate pairs in the input\n");
    fprintf(stderr, "  -g G     genome size\n");
    fprintf(stderr, "  -l L     length of a single read\n");
    fprintf(stderr, "  -c C     desired coverage in the output frags\n");
    fprintf(stderr, "  -n N     name (prefix) of the reads\n");
    fprintf(stderr, "\n");
    if (NAME == NULL)
      fprintf(stderr, "ERROR: no name supplied with -n.\n");
    if (GENOMESIZE == 0)
      fprintf(stderr, "ERROR: no genome size supplied with -g\n");
    if (READLENGTH == 0)
      fprintf(stderr, "ERROR: no read length supplied with -l\n");
    if (COVERAGE == 0)
      fprintf(stderr, "ERROR: no desired coverage supplied with -c\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  if (NUMREADS == 0) {
    uint32  Ac = 0;
    uint32  Bc = 0;

    fprintf(stderr, "Counting the number of reads in the input.\n");

    sprintf(path1, "%s.1.fastq", NAME);
    sprintf(path2, "%s.2.fastq", NAME);

    errno = 0;
    Ai = fopen(path1, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", path1, strerror(errno)), exit(1);

    errno = 0;
    Bi = fopen(path2, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", path2, strerror(errno)), exit(1);

    for (Ac=0; !feof(Ai); Ac++)
      Ar.read(Ai);
    for (Bc=0; !feof(Bi); Bc++)
      Ar.read(Bi);

    fclose(Ai);
    fclose(Bi);

    fprintf(stderr, "Found "F_U64" reads in '%s'\n", Ac, path1);
    fprintf(stderr, "Found "F_U64" reads in '%s'\n", Bc, path2);

    if (Ac != Bc) {
      fprintf(stderr, "ERROR:  Number of reads must be the same.\n");
      exit(1);
    }

    NUMREADS = Ac;
  }

  uint32    I = COVERAGE * GENOMESIZE / READLENGTH;

  if (I > NUMREADS) {
    fprintf(stderr, "ERROR: not enough reads, %d in input, %d needed, for desired coverage %d.\n",
            NUMREADS, I, COVERAGE);
    exit(1);
  }

  ids  = new uint32   [NUMREADS];
  save = new char     [NUMREADS];

  srand48(time(NULL));

  //  Make a list of IDs, and clear the array that tells us to emit a read.
  for (int i=0; i<NUMREADS; i++) {
    ids[i]  = i;
    save[i] = 0;
  }

  //  Randomize the ID list.
  for (int i=0; i<NUMREADS; i++) {
    int p = lrand48() % NUMREADS;
    int a = ids[p];

    ids[p] = ids[i];
    ids[i] = a;
  }

  //  Take the first I things from the randomized ID list, mark them as "to be output".
  for (uint32 i=0; i<I; i++)
    save[ids[i]] = 1;

  sprintf(path1, "%s.1.fastq", NAME);
  sprintf(path2, "%s.2.fastq", NAME);

  errno = 0;
  Ai = fopen(path1, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path1, strerror(errno)), exit(1);

  errno = 0;
  Bi = fopen(path2, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path2, strerror(errno)), exit(1);

  sprintf(path1, "%s.%03dx.1.fastq", NAME, COVERAGE);
  sprintf(path2, "%s.%03dx.2.fastq", NAME, COVERAGE);

  errno = 0;
  Ao = fopen(path1, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path1, strerror(errno)), exit(1);

  errno = 0;
  Bo = fopen(path2, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path2, strerror(errno)), exit(1);

  fprintf(stderr, "Extracting %d mate pairs into %s and %s\n",
          I, path1, path2);

  int i=0;
  int s=0;
  for (; !feof(Ai); i++) {
    Ar.read(Ai);
    Br.read(Bi);

    if (save[i]) {
      Ar.write(Ao);
      Br.write(Bo);
      s++;
    }
  }

  if (i > NUMREADS) {
    fprintf(stderr, "WARNING:  There are %d mates in the input; you claimed there are %d (-t option) mates.\n",
            i, NUMREADS);
    fprintf(stderr, "WARNING:  Result is not a random sample of the input file.\n");
  }

  if (i < NUMREADS) {
    fprintf(stderr, "WARNING:  There are %d mates in the input; you claimed there are %d (-t option) mates.\n",
            i, NUMREADS);
    fprintf(stderr, "WARNING:  Result is only %f X coverage.\n", (double)s * READLENGTH / GENOMESIZE);
  }

  fclose(Ai);
  fclose(Bi);

  fclose(Ao);
  fclose(Bo);

  return(0);
}
