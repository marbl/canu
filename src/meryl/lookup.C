
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

  kmerCountFileReader   *rr = new kmerCountFileReader(argv[1], false, true);
  kmerCountExactLookup  *ll = new kmerCountExactLookup(rr);

  delete rr;

  rr = new kmerCountFileReader(argv[2], false, true);

  uint64  tested = 0;
  uint64  found  = 0;

  while (rr->nextMer()) {
    kmer    k = rr->theFMer();
    uint32  c = rr->theCount();

    tested++;

    if (ll->value(k))
      found++;

    if ((tested % 100000) == 0)
      fprintf(stderr, "Tested %lu kmers, found %lu.\n", tested, found);
  }

  delete ll;
  delete rr;

  fprintf(stderr, "Tested %lu kmers, found %lu.\n", tested, found);

  exit(0);
}
