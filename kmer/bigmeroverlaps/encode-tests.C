#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>     //  gettimeofday()
#include <sys/utsname.h>  //  uname()
#include <sys/resource.h> //  getrusage()

#include "gmp.h"

//  If you ever happen to want to link against installed libraries
//  in a given directory, LIBDIR, you must either use libtool, and
//  specify the full pathname of the library, or use the `-LLIBDIR'
//  flag during linking and do at least one of the following:
//     - add LIBDIR to the `LD_LIBRARY_PATH' environment variable
//       during execution
//     - add LIBDIR to the `LD_RUN_PATH' environment variable
//       during linking
//     - use the `-Wl,--rpath -Wl,LIBDIR' linker flag


double
getTime(void) {
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return(tp.tv_sec + (double)tp.tv_usec / 1000000.0);
}


int
main(int argc, char **argv) {
  mpz_t   a;
  mpz_t   a1;
  mpz_t   a2;
  mpz_t   b;
  mpz_t   r;

  fprintf(stderr, "sizeof(mpz_t) = %d\n", sizeof(mpz_t));

  mpz_init(a);
  mpz_init(a1);
  mpz_init(a2);
  mpz_init2(b, 200);
  mpz_init2(r, 200);

  mpz_set_ui(a, 0x12345678);
  mpz_init_set_str(a, "abcdef1234567890111abcdef1234567890abcdef1234567890", 16);

  unsigned long long   dataout[32] = {0};
  size_t               countout    =  0;

  mpz_export(dataout, &countout, 1, 8, 0, 0, a);

  fprintf(stderr, "countout = %ld\n", countout);

  for (size_t i=0; i<countout; i++)
    fprintf(stderr, "0x%08lx\n", dataout[i]);


  countout = 4;
  dataout[0] = 0x1020304050607080llu;
  dataout[1] = 0x8090a0b0c0d0e0f0llu;
  dataout[2] = 0xaaabbbcccdddeeefllu;
  dataout[3] = 0x857237157f7ae742llu;

  double  tbeg, tend;
  int     imax;


  //  import:
  //  1)  mpz number
  //  2)  number of elements to import
  //  3)  1 - most sig word first, -1 otherwise
  //  4)  size of word in bytes
  //  5)  endian - 1=msb, -1=lsb, 0=native
  //  6)  nails - most sig nail bits are skipped in each word, 0 for full words
  //  7)  array of data values
  //
  tbeg = getTime();
  imax = 100000000;
  for (int i=imax; --i; )
    mpz_import(a, countout, 1, 8, 0, 0, dataout);
  tend = getTime();

  fprintf(stderr, "%f us each (single four word import)\n", 1000000 * (tend - tbeg) / imax);
  fprintf(stderr, "%f / second\n", imax / (tend - tbeg));

  //  mpz_import seems to take 0.022499 us each



  countout = 4;
  dataout[0] = 0x0000000000607080llu;
  dataout[1] = 0x8090a0b0c0d0e0f0llu;
  dataout[2] = 0xaaabbbcccdddeeefllu;
  dataout[3] = 0x8572370000000000llu;

  mpz_import(b, countout, 1, 8, 0, 0, dataout);

  tbeg = getTime();
  imax = 100000000;
  for (int i=imax; --i; ) {
    countout = 4;
    dataout[0] = 0x1020304050607080llu;
    dataout[1] = 0x8090a0b0c0d0e0f0llu;
    dataout[2] = 0xaaabbbcccdddeeefllu;
    dataout[3] = 0x857237157f7ae742llu;

    mpz_import(a, countout, 1, 8, 0, 0, dataout);

    mpz_mod(r, a, b);
  }
  tend = getTime();

  fprintf(stderr, "%f us each (2x four word import, mod)\n", 1000000 * (tend - tbeg) / imax);
  fprintf(stderr, "%f / second\n", imax / (tend - tbeg));

  //  mpz_import + mpz_mod is 0.273291 us, so
  //  mpz_mod is 0.25 us (the same, even if we remove the import from the loop).

  char  outnumber[1024];

  mpz_get_str(outnumber, 16, a);
  fprintf(stdout, "%s\n", outnumber);
}
