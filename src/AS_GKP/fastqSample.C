
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

const char *mainid = "$Id: fastqSample.C,v 1.6 2011-07-29 01:56:11 brianwalenz Exp $";

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

  uint64    NUMINPUT    = 0;        //  Number of pairs in the input
  char     *NAME        = NULL;     //  Prefix name

  uint64    GENOMESIZE  = 0;        //  Size of the genome in bp
  uint64    READLENGTH  = 0;        //  For mated reads, 2x read size
  uint64    COVERAGE    = 0;        //  Desired coverage in output

  uint64    NUMOUTPUT   = 0;        //  Number of pairs to output

  char      path1[FILENAME_MAX];
  char      path2[FILENAME_MAX];

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-n") == 0) {
      NAME = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      //  Total reads in input
      NUMINPUT = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-g") == 0) {
      //  Genome size
      GENOMESIZE = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      //  Length of read
      READLENGTH = atoi(argv[++arg]) * 2;

    } else if (strcmp(argv[arg], "-c") == 0) {
      //  Desired coverage
      COVERAGE = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-p") == 0) {
      //  Desired coverage
      NUMOUTPUT = atoi(argv[++arg]);


    } else {
      err++;
    }

    arg++;
  }
  if (NAME == NULL)
    err++;
  if ((GENOMESIZE == 0) && (NUMOUTPUT == 0))
    err++;
  if ((READLENGTH == 0) && (NUMOUTPUT == 0))
    err++;
  if ((COVERAGE == 0) && (NUMOUTPUT == 0))
    err++;
  if (err) {
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "    -n N     name (prefix) of the reads\n");
    fprintf(stderr, "    -t T     total number of mate pairs in the input\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Method 1: specify desired coverage:\n");
    fprintf(stderr, "    -g G     genome size\n");
    fprintf(stderr, "    -l L     length of a single read\n");
    fprintf(stderr, "    -c C     desired coverage in the output reads (INTEGER ONLY)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Method 2: specify desired number of output pairs\n");
    fprintf(stderr, "    -p N     generate 2N reads, or N pairs of reads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Program will compute the number of pairs of reads to sample to generate\n");
    fprintf(stderr, "a set with the desired coverage.\n");
    fprintf(stderr, "\n");
    if (NAME == NULL)
      fprintf(stderr, "ERROR: no name supplied with -n.\n");
    if ((GENOMESIZE == 0) && (NUMOUTPUT == 0))
      fprintf(stderr, "ERROR: no genome size supplied with -g\n");
    if ((READLENGTH == 0) && (NUMOUTPUT == 0))
      fprintf(stderr, "ERROR: no read length supplied with -l\n");
    if ((COVERAGE == 0) && (NUMOUTPUT == 0))
      fprintf(stderr, "ERROR: no desired coverage supplied with -c\n");
    if (NUMOUTPUT == 0)
      fprintf(stderr, "ERROR: no desired number of pairs to output supplied with -p\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  if (NUMINPUT == 0) {
    uint64  Ac = 0;
    uint64  Bc = 0;

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

    NUMINPUT = Ac;
  }

  if (NUMOUTPUT == 0)
    NUMOUTPUT = COVERAGE * GENOMESIZE / READLENGTH;

  if (NUMOUTPUT > NUMINPUT) {
    fprintf(stderr, "ERROR: not enough reads, "F_U64" in input, "F_U64" needed, for desired coverage "F_U64".\n",
            NUMINPUT, NUMOUTPUT, COVERAGE);
    exit(1);
  }

  ids  = new uint32   [NUMINPUT];
  save = new char     [NUMINPUT];

  srand48(time(NULL));

  //  Make a list of IDs, and clear the array that tells us to emit a read.
  for (uint64 i=0; i<NUMINPUT; i++) {
    ids[i]  = i;
    save[i] = 0;
  }

  //  Randomize the ID list.
  for (uint64 i=0; i<NUMINPUT; i++) {
    uint64 p = (uint64)lrand48() % NUMINPUT;
    uint32 a = ids[p];

    ids[p] = ids[i];
    ids[i] = a;
  }

  //  Take the first NUMOUTPUT things from the randomized ID list, mark them as "to be output".
  for (uint64 i=0; i<NUMOUTPUT; i++)
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

  if (COVERAGE > 0) {
    sprintf(path1, "%s.%03"F_U64P"x.1.fastq", NAME, COVERAGE);
    sprintf(path2, "%s.%03"F_U64P"x.2.fastq", NAME, COVERAGE);
  } else {
    sprintf(path1, "%s.n="F_U64".1.fastq", NAME, NUMOUTPUT);
    sprintf(path2, "%s.n="F_U64".2.fastq", NAME, NUMOUTPUT);
  }

  errno = 0;
  Ao = fopen(path1, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path1, strerror(errno)), exit(1);

  errno = 0;
  Bo = fopen(path2, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", path2, strerror(errno)), exit(1);

  fprintf(stderr, "Extracting "F_U64" mate pairs into %s and %s\n",
          NUMOUTPUT, path1, path2);

  uint64 i=0;
  uint64 s=0;

  Ar.read(Ai);
  Br.read(Bi);

  for (; !feof(Ai); i++) {
    if (save[i]) {
      Ar.write(Ao);
      Br.write(Bo);
      s++;
    }

    Ar.read(Ai);
    Br.read(Bi);
  }

  if (i > NUMINPUT) {
    fprintf(stderr, "WARNING:  There are "F_U64" mates in the input; you claimed there are "F_U64" (-t option) mates.\n",
            i, NUMINPUT);
    fprintf(stderr, "WARNING:  Result is not a random sample of the input file.\n");
  }

  if (i < NUMINPUT) {
    fprintf(stderr, "WARNING:  There are "F_U64" mates in the input; you claimed there are "F_U64" (-t option) mates.\n",
            i, NUMINPUT);
    fprintf(stderr, "WARNING:  Result is only %f X coverage.\n", (double)s * READLENGTH / GENOMESIZE);
  }

  fclose(Ai);
  fclose(Bi);

  fclose(Ao);
  fclose(Bo);

  return(0);
}
