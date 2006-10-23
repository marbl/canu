#include <unistd.h>
#include <time.h>
#include <math.h>

#include "util++.H"

u32bit wordSize = 41;
u32bit testSize = 48 * 1024 * 1024;
u32bit arrySize =  8 * 1024 * 1024;

int
main(int argc, char **argv) {

  //  Test the bitPackedArray by writing a bunch of random gibberish
  //  to it, and see if it's the same.

  u32bit  *pos = new u32bit [testSize];
  u64bit  *val = new u64bit [testSize];
  u64bit  *ans = new u64bit [arrySize];

  bitPackedArray *ARR  = new bitPackedArray(wordSize);
  u32bit          fail = u32bitZERO;

  fprintf(stderr, "Touching the end of the array and clearing.\n");
  ARR->set(arrySize, 0);
  ARR->clear();

  fprintf(stderr, "Generating random test data.\n");

  mt_s *mtctx = mtInit(time(NULL));

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

      if (fail > 20)
        return(1);
    }

  delete [] pos;
  delete [] val;
  delete [] ans;

  return(fail);
}
