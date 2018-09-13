
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
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
