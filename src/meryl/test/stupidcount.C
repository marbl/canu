
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "bio++.H"

//  Reads a sequence file, outputs a list of the mers in it.  You can
//  then pipe this to unix sort and uniq to do a mercount.  You
//  probably don't want to count large things this way...

int
main(int argc, char **argv) {
  char    *seqName = 0L;
  uint32   merSize = 20;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-s") == 0) {
      seqName = argv[++arg];
    } else if (strcmp(argv[arg], "-m") == 0) {
      merSize = strtouint32(argv[++arg], 0L);
    }
    arg++;
  }

  if (seqName == 0L) {
    fprintf(stderr, "usage: %s [-m mersize] -s seqfile.fasta\n", argv[0]);
    exit(1);
  }

  seqStream       *CS = new seqStream(seqName, true);
  merStream       *MS = new merStream(new kMerBuilder(merSize), CS);
  char             str[1024];

  while (MS->nextMer())
    fprintf(stdout, "%s\n", MS->theFMer().merToString(str));

  delete MS;
  delete CS;

  exit(0);
}
