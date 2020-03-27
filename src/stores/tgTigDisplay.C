
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"

#include "sqStore.H"
#include "tgStore.H"


int
main(int argc, char **argv) {
  tgTig  tig;
  char  *seqName     = NULL;
  char  *tigFileName = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigFileName = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if (seqName == NULL)
    err++;
  if (tigFileName == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -S seqStore -t tigFile\n", argv[0]);
    exit(1);
  }

  sqStore  *seqStore = new sqStore(seqName);

  FILE *F = fopen(tigFileName, "r");

  tig.loadFromStream(F);

  AS_UTL_closeFile(F, tigFileName);

  uint32  displayWidth    = 250;
  uint32  displaySpacing  = 10;
  bool    withQV          = false;
  bool    withDots        = true;

  tig.display(stdout, seqStore, displayWidth, displaySpacing, withQV, withDots);

  exit(0);
}
