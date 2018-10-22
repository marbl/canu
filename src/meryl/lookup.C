
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "kmers.H"
#include "bits.H"

int
main(int argc, char **argv) {
  char   *inputDBname  = NULL;
  uint32  minV         = 0;
  uint32  maxV         = UINT32_MAX;
  char   *queryDBname = NULL;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-M") == 0) {   //  INPUT READS and RANGE TO PROCESS
      inputDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-min") == 0) {
      minV = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      maxV = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-Q") == 0) {
      queryDBname = argv[++arg];

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (inputDBname == NULL)
    err.push_back("No input meryl database (-M) supplied.\n");
  if (queryDBname == NULL)
    err.push_back("No query meryl database (-Q) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "(no help yet)\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }



  kmerCountFileReader   *rr = new kmerCountFileReader(inputDBname);
  kmerCountExactLookup  *ll = new kmerCountExactLookup(rr, minV, maxV);

  delete rr;

  rr = new kmerCountFileReader(queryDBname);

  uint64  tested    = 0;
  uint64  kmerFound = 0;
  uint64  valuFound = 0;

  while (rr->nextMer()) {
    kmer    k = rr->theFMer();
    uint32  c = rr->theCount();
    uint32  v = ll->value(k);

    tested++;

    if (v > 0)
      kmerFound++;

    if (c == v)
      valuFound++;
  }

  delete ll;
  delete rr;

  fprintf(stderr, "Tested %12lu kmers.\n", tested);
  fprintf(stderr, "Found  %12lu kmers.                     %12lu missed.\n", kmerFound, tested-kmerFound);
  fprintf(stderr, "Found  %12lu kmers with correct value.  %12lu missed.\n", valuFound, tested-valuFound);

  exit(0);
}
