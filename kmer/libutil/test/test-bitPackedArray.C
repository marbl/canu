#include <unistd.h>
#include <time.h>
#include <math.h>

#include "util++.H"

u32bit wordSize = 41;
u32bit testSize =  1 * 1024 * 1024;
u32bit arrySize =  1 * 1024 * 1024;

int
u64bitcompare(const void *a, const void *b) {
  const u64bit   A = *(const u64bit *)a;
  const u64bit   B = *(const u64bit *)b;
  if (A<B) return(-1);
  if (A>B) return(1);
  return(0);
}

int
main(int argc, char **argv) {

  mt_s *mtctx = mtInit(time(NULL));

  //  Test the bitPackedArray by writing a bunch of random gibberish
  //  to it, and see if it's the same.

  u32bit  *pos = new u32bit [testSize];
  u64bit  *val = new u64bit [testSize];
  u64bit  *ans = new u64bit [arrySize];

  bitPackedArray *ARR  = new bitPackedArray(wordSize, 16);
  u32bit          fail = u32bitZERO;

#if 1
  fprintf(stderr, "Touching the end of the array and clearing.\n");
  //ARR->set(arrySize, 0);
  //ARR->clear();

  fprintf(stderr, "Generating random test data.\n");

  //  Hit every element first, just to do it
  for (u32bit i=0; i<arrySize; i++) {
    pos[i]       = i;
    val[i]       = mtRandom64(mtctx);
    val[i]      &= u64bitMASK(wordSize);
    ans[pos[i]]  = val[i];
  }

  //  Then hit random elements, with replacement, looking for bugs
  for (u32bit i=arrySize; i<testSize; i++) {
    pos[i]       = mtRandom32(mtctx) % arrySize;
    val[i]       = mtRandom64(mtctx);
    val[i]      &= u64bitMASK(wordSize);
    ans[pos[i]]  = val[i];
  }

  fprintf(stderr, "Filling array.\n");

  for (u32bit i=0; i<testSize; i++)
    ARR->set(pos[i], val[i]);

  fprintf(stderr, "Validating array.\n");

  for (u32bit i=0; i<arrySize; i++)
    if (ARR->get(i) != ans[i]) {
      fprintf(stderr, "FAIL at i="u32bitFMT"\n", i);
      fail++;

      if (fail > 1024) {
        fprintf(stderr, "bitPackedArray has errors, aborting!\n");
        return(1);
      }
    }

  if (fail) {
    fprintf(stderr, "bitPackedArray had "u32bitFMT" errors.\n", fail);
    return(1);
  }

  fprintf(stderr, "OK!\n");
#endif

  delete    ARR;
  delete [] pos;
  delete [] val;
  delete [] ans;

  //
  //
  //

  for (u32bit testNum=0; testNum<32; testNum++) {
    u32bit  thisTestSize = 0;
    u32bit  thisWordSize = 0;

    //  Test a BIG heap the first iteration.
    if (testNum == 0) {
      thisTestSize = 857353; //23987153;
      thisWordSize = 63;

      fprintf(stderr, "Building heap "u32bitFMT" (wordsize="u32bitFMT" testsize="u32bitFMT").\n",
              testNum, thisWordSize, thisTestSize);
    } else {
      thisTestSize = (mtRandom64(mtctx) % (2 * testNum)) * 1024 + 1024;
      thisWordSize = (mtRandom64(mtctx) % 63) + 1;
    }

    u32bit  blockSize = mtRandom64(mtctx) % 32 + 1;
    bitPackedHeap  *HEAP = new bitPackedHeap(thisWordSize, blockSize);

    val = new u64bit [thisTestSize];
    for (u32bit i=0; i<thisTestSize; i++) {
      val[i]       = mtRandom64(mtctx);
      val[i]      &= u64bitMASK(thisWordSize);
      HEAP->add(val[i]);
    }

    fprintf(stderr, "Testing heap "u32bitFMT" (wordsize="u32bitFMT" testsize="u32bitFMT").\n",
            testNum, thisWordSize, thisTestSize);

    qsort(val, thisTestSize, sizeof(u64bit), u64bitcompare);

    for (u32bit i=0; i<thisTestSize; i++) {
      u64bit  h = HEAP->get();

      //fprintf(stderr, "val["u32bitFMT"]="u64bitFMT" -- HEAP="u64bitFMT"\n", i, val[i], h);

      if (val[i] != h) {
        fprintf(stderr, "val["u32bitFMT"]="u64bitFMT" !! HEAP="u64bitFMT"\n", i, val[i], h);
        fail++;
        if (fail > 25) {
          fprintf(stderr, "bitPackedHeap has errors, aborting!\n");
          return(1);
        }
      }
    }
    
    if (fail) {
      fprintf(stderr, "bitPackedHeap had "u32bitFMT" errors.!\n", fail);
      return(1);
    }

    delete    HEAP;
    delete [] val;
  }

  fprintf(stderr, "OK!\n");

  return(fail);
}

