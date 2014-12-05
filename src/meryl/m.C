#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util++.H"
#include "bio++.H"
#include "meryl.H"

#include "m-heap.H"

int
main(int argc, char **argv) {
  bool     beVerbose   = false;
  uint64   merSize     = 20;
  uint64   memLimit    = 768;
  char    *inName      = 0L;
  char    *outName     = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-verbose", 2) == 0) {
      beVerbose = true;
    } else if (strncmp(argv[arg], "-mersize", 4) == 0) {
      merSize = strtouint64(argv[++arg], 0L);
    } else if (strncmp(argv[arg], "-memory", 4) == 0) {
      memLimit = strtouint64(argv[++arg], 0L) * 1024 * 1024;
    } else if (strncmp(argv[arg], "-sequence", 2) == 0) {
      inName = argv[++arg];
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (inName == 0L) {
    fprintf(stderr, "usage: %s [-v] [-m mersize] [-memory Nmb] [-s seq.fasta]\n", argv[0]);
    exit(1);
  }

  outName = new char [strlen(inName) + 1];
  strcpy(outName, inName);

  seqStream  *seqstr = new seqStream(inName);
  seqStore   *seqsto = new seqStore(outName, seqstr);

  uint64      memUsed = seqsto->loadStoreInCore();
  uint64      numMers = seqsto->numberOfACGT();

#warning needed exact number of mers here

  fprintf(stderr, "Found "uint64FMT" mers in file of size "uint64FMT"\n", numMers, memUsed);

  if (memUsed > memLimit) {
    fprintf(stderr, "ERROR:  two-bit encoded sequence file is bigger than allowed memory usage.\n");
    exit(1);
  }

  //  Allocate a heap to fill up the rest of space

  //  Allocate a bitPackedHeap to store N merSize*2 integeers.
  //  N = (memLimit - memUsed) * 8 / (merSize * 2)
  //
  //  The bitPackedHeap doesn't care about the maximum size, only
  //  about the block size.
  //
  uint64   pointerWidth = logBaseTwo64(numMers);
  bitPackedMerHeap  *heap = new bitPackedMerHeap(seqsto, pointerWidth, 8 * 1024);

  speedCounter *S;

  uint64 N = (memLimit - memUsed) * 8 / pointerWidth;
  uint64 M = 0;

  fprintf(stderr, "Can store "uint64FMT" mer pointers of size "uint64FMT" in the heap.\n", N, pointerWidth);

  kMer mer;

  if (N > numMers)
    N = numMers;

  //  Initialize the heap with some numbers
  //
  S = new speedCounter(" Loading heap: %7.2f Mmers -- %8.1f mers/second\r", 1.0, 0x1ffff, beVerbose);
  while (M < N) {

#if 0
    heap->add(M);
    heap->get(mer);
    fprintf(stdout, "ADD "uint64FMT" -- %s\n", M, mer.merToString(str));
#endif

    heap->add(M++);
    S->tick();
  }
  delete S;

  //  Until we run out of mers, write things out of the heap.
  //
  S = new speedCounter(" Cycling heap: %7.2f Mmers -- %8.1f Mmers/second\r", 1.0, 0x1fff, beVerbose);
  while (M < numMers) {
    heap->add(M++);
    heap->get(mer);
    //fprintf(stdout, "GOT "uint64FMT" -- %s\n", M, mer.merToString(str));
    S->tick();
  }
  delete S;

  //  And finally, flush the heap.
  //
  S = new speedCounter(" Dumping heap: %7.2f Mmers -- %8.1f Mmers/second\r", 1.0, 0x1fff, beVerbose);
  uint64 idx = heap->get(mer);
  while (idx != ~uint64ZERO) {
    //fprintf(stdout, "OUT "uint64FMT" -- %s\n", idx, mer.merToString(str));
    idx = heap->get(mer);
    S->tick();
  }
  delete S;
}
