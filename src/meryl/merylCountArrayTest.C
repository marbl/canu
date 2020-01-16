
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
 *    Brian P. Walenz beginning on 2017-SEP-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "meryl.H"
#include "strings.H"
//#include "system.H"

#include "mt19937ar.H"

mtRandom  *mt = NULL;


void
display(char *l, kmdata s) {
  uint64 a = (s >> 64);
  uint64 b =  s;

  fprintf(stderr, "%s 0x%016lx 0x%016lx\n", l, a, b);
}


kmdata
setWord(uint32 w, uint64 t) {
  kmdata s;

  s   = t;
  s <<= 64;
  s  |= t;

  s <<= (128 - w);
  s >>= (128 - w);

  return(s);
}


kmdata
setWord(uint32 w, uint64 a, uint64 b) {
  kmdata s;

  s   = a;
  s <<= 64;
  s  |= b;

  s <<= (128 - w);
  s >>= (128 - w);

  return(s);
}


kmdata
setRandomWord(uint32 w) {
  kmdata s;

  s   = mt->mtRandom64();
  s <<= 64;
  s  |= mt->mtRandom64();

  s <<= (128 - w);
  s >>= (128 - w);

  return(s);
}





int
main(int argc, char **argv) {
  uint32 seed  = 0;
  uint32 iters = 0;
  uint32 words = 0;

  uint32 widthMin = 0;
  uint32 widthMax = 0;

  int err=0;
  int arg=1;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-seed") == 0) {
      seed = strtouint32(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-iters") == 0) {
      iters = strtouint32(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-words") == 0) {
      words = strtouint32(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-width") == 0) {
      decodeRange(argv[++arg], widthMin, widthMax);
    }

    else if (strcmp(argv[arg], "-iter") == 0) {
    }

    else if (strcmp(argv[arg], "-iter") == 0) {
    }

    else if (strcmp(argv[arg], "-iter") == 0) {
    }

    else {
      fprintf(stderr, "usage: ...\n");
    }

    arg++;
  }

  mt = new mtRandom(strtouint32(argv[2]));

  for (uint32 w=widthMin; w<=widthMax; w++) {
    merylCountArray   *A = new merylCountArray;

    A->initializeForTesting(w, words);

    //  Generate 'iters' random words.

    fprintf(stderr, "\n");
    fprintf(stderr, "Generate %u values of width %u.\n", words, w);

    kmdata *vals = new kmdata [iters];

    for (uint32 ii=0; ii<iters; ii++)
      vals[ii] = setRandomWord(w);

    //  Insert them into the table.  Check that they insert corretly.

    fprintf(stderr, "Insert and check each.\n");

    for (uint32 ii=0; ii<iters; ii++) {
      A->add(vals[ii]);
      //A->dumpData();

#if 1
      kmdata t = A->getSimple(ii);
      kmdata f = A->get(ii);

      //fprintf(stderr, "\n");
      //fprintf(stderr, "insert [%u] ", ii);
      //display("", vals[ii]);

      if (t != vals[ii]) {
        fprintf(stderr, "FAILED at iter ii %u\n", ii);
        display("val", vals[ii]);
        display("t  ", t);
      }
      assert(t == vals[ii]);

      if (t != f) {
        fprintf(stderr, "FAILED at iter ii %u\n", ii);
        display("t", t);
        display("f", f);
      }
      assert(t == f);
#endif
    }

    //  Check that all are still correct.

    fprintf(stderr, "Check all.\n");
#if 1
    for (uint32 ii=0; ii<iters; ii++) {
      kmdata t = A->get(ii);

      if (t != vals[ii]) {
        fprintf(stderr, "FAILED at iter ii %u\n", ii);
        display("val", vals[ii]);
        display("t  ", t);
      }
      assert(t == vals[ii]);
    }
#endif
    //  Cleanup for next loop.

    A->dumpStats();

    delete    A;
    delete [] vals;
  }

  return(0);
}
