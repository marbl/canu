#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#include "util.h"

//  Does a quick test of memory mapped files.  First, writes a small
//  file, then it reads it back, checking the data.
//
//  Takes one optional argument, the size in MB of the file to map.

int
main(int argc, char **argv) {
  size_t   lw;
  u32bit  *ww = 0L;
  u32bit   idx = 0;
  u32bit   err = 0;
  FILE    *out;
  u32bit   blockSize = 1048576;
  u32bit   numBlocks = 32;

  if (argc == 2)
    numBlocks = strtou32bit(argv[1], 0L);

  //  The file must exist, and it must be large enough to contain all
  //  that we want to write.  So, we create the file and fill it with
  //  junk.
  //
  ww  = (u32bit *)malloc(sizeof(u32bit) * blockSize);
  if (ww == NULL) {
    fprintf(stderr, "can't allocate %d u32bit's for clearing the file.\n", blockSize);
    exit(1);
  }
  errno = 0;
  out = fopen("mmap.test.junk", "w");
  if (errno) {
    fprintf(stderr, "can't open 'mmap.test.junk' to fill with junk: %s\n", strerror(errno));
    exit(1);
  }
  for (idx=0; idx<numBlocks; idx++) {
    fprintf(stderr, "Writing initial blocks: "u32bitFMT"/"u32bitFMT"\r", idx, numBlocks), fflush(stderr);
    fwrite(ww, sizeof(u32bit), 1048576, out);
    if (errno) {
      fprintf(stderr, "can't write to 'mmap.test.junk': %s\n", strerror(errno));
      exit(1);
    }
  }
  fclose(out);
  free(ww);
  fprintf(stderr, "\n");

  //  Now, map it, and fill it with real data.
  //
  ww = (u32bit *)mapFile("mmap.test.junk", &lw, 'w');
  for (idx=0; idx<numBlocks * blockSize; idx++) {
    if ((idx & 0xfff) == 0)
      fprintf(stderr, "Writing: "u32bitFMT"/"u32bitFMT"\r", idx, numBlocks * blockSize), fflush(stderr);
    ww[idx] = idx;
  }
  unmapFile(ww, lw);
  fprintf(stderr, "\n");

  //  Map again, and check the data.
  //
  ww = mapFile("mmap.test.junk", &lw, 'r');
  for (idx=0; idx<numBlocks * blockSize; idx++) {
    if ((idx & 0xfff) == 0)
      fprintf(stderr, "Verifying: "u32bitFMT"/"u32bitFMT"\r", idx, numBlocks * blockSize), fflush(stderr);
    if (ww[idx] != idx)
      err++;
  }
  unmapFile(ww, lw);
  fprintf(stderr, "\n");

  unlink("mmap.test.junk");

  return (err != 0);  
}
