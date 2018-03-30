
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
 *    Brian P. Walenz beginning on 2017-SEP-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "seqCache.H"
#include "merStream.H"

#include "libmeryl.H"
#include "existDB.H"


#define  OP_NONE    0
#define  OP_STATS   1
#define  OP_REGIONS 2
#define  OP_DETAILS 3

int
main(int argc, char **argv) {
  uint32    merSize    = 16;
  char     *merylFile  = NULL;
  char     *existFile  = NULL;
  char     *fastaFile  = NULL;

  uint32    loCount    = 0;
  uint32    hiCount    = UINT32_MAX;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      merSize = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      merylFile = argv[++arg];
    } else if (strcmp(argv[arg], "-e") == 0) {
      existFile = argv[++arg];

    } else if (strcmp(argv[arg], "-l") == 0) {
      loCount = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-h") == 0) {
      hiCount = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-f") == 0) {
      fastaFile = argv[++arg];

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((existFile == NULL) && ((merylFile == NULL) || (fastaFile == NULL))) {
    fprintf(stderr, "usage: %s -m mersize [-s mers] -e dump-db [-l int] [-h int] [-f fasta] > output\n", argv[0]);
    fprintf(stderr, "\t-e: dump (create a new db) from the given -s mers to have random access\n");
    fprintf(stderr, "\t-l: kmers with counts >= l will be reported\n");
    fprintf(stderr, "\t-h: kmers with counts <= h will be reported\n");
    fprintf(stderr, "\t-f: fasta file, or -f - for stdin (assuming also being fasta)\n");
    fprintf(stderr, "\t    *If no fasta file is porvided, this dump creates a query file (-e) and exits.\n");
    fprintf(stderr, "\toutput: readid\ttotal-mers-in-read\ttotal-mers-in-s(0 if called from -e)\tfound-mers-in-read\n");
    exit(1);
  }

  existDB       *E = NULL;

  //  Load existence mers from the saved existDB, or create a new existDB and maybe save to disk.

  if ((existFile) && (AS_UTL_fileExists(existFile))) {
    E = new existDB(existFile);
  }

  else {
    E = new existDB(merylFile, merSize, existDBcounts, loCount, hiCount);

    if (existFile)
      E->saveState(existFile);
    
    if (fastaFile == NULL)
      exit(0);
  }

  //  Open the fasta/fastq input for random access, though we only read sequentially.

  fprintf(stderr, "..Start cashing\n");
  seqCache *F = new seqCache(fastaFile);
  fprintf(stderr, "..End cashing\n");

  for (uint32 Sid=0; Sid < F->getNumberOfSequences(); Sid++) {
    seqInCore  *S  = F->getSequenceInCore(Sid);
    merStream  *MS = new merStream(new kMerBuilder(merSize),
                                   new seqStream(S->sequence(), S->sequenceLength()),
                                   true, true);

    uint32   totalMers = 0;
    uint32   foundMers = 0;

    //  Stream the kmers in the sequence by the existDB table, counting how many mers
    //  are in the sequence and found in the table.

    while (MS->nextMer()) {
      totalMers++;

      //  merStream can also tell the position of the kmer in the sequence
      //pos = MS->thePositionInSequence();

      //  merStream can print the kmer as ASCII.
      //char merString[32];  //  Big enough to hold the mer
      //fprintf(stdout, "mer=%s\n", MS->theFMer().merToString(merString));

      if (E->count(MS->theFMer()) + E->count(MS->theRMer()) > 0)
        foundMers++;
    }

    fprintf(stdout, "%s\t%u\t%lu\t%u\n", S->header(), totalMers, E->numberOfMers(), foundMers);

    delete MS;
    delete S;
  }


  delete F;
  delete E;
}
