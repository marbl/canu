#include <stdio.h>
#include "../readBuffer.H"


int
main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "usage: %s file-to-read\n", argv[0]);
    exit(1);
  }

#if 0
  readBuffer B(argv[1], 999);
  for ( ; B.eof() == false; B.next())
    putchar(B.get());
  B.seek(0);
  for ( ; B.eof() == false; B.next())
    putchar(B.get());
#endif

#if 1
  readBuffer B(argv[1], 0);
  for ( ; B.eof() == false; B.next())
    putchar(B.get());
  fflush(stdout);
  B.seek(0);
  for ( ; B.eof() == false; B.next())
    putchar(B.get());
#endif

  //

#if 0
  readBuffer C(argv[1], 100);
  char       c[10000];
  int        l;

  while ((l = C.read(c, 10000))) {
    //fprintf(stderr, "READ %d bytes\n", l);
    for (int i=0; i<l; i++)
      putchar(c[i]);
  }
#endif

  fflush(stdout);
  fprintf(stderr, "All done!\n");
  fflush(stderr);
}

