#include <stdio.h>
#include <stdlib.h>

#include "libbri.H"

#define NUMMERS  (1024*1024)

int
main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s <1.fasta> <2.fasta>\n", argv[0]);
    fprintf(stderr, "Builds an array of the mers in 1, and brute-force\n");
    fprintf(stderr, "counts the number of times those occur in 2.\n");
    exit(1);
  }

  merStream *O = new merStream(20, argv[1]);
  merStream *T = new merStream(20, argv[2]);

  u64bit    *F = new u64bit [1024 * 1024];
  u64bit    *R = new u64bit [1024 * 1024];

  u32bit     l;

  for (l=0; O->nextMer() && l<NUMMERS; l++) {
    F[l] = O->theFMer();
    R[l] = O->theRMer();
  }

  fprintf(stderr, "Found %lu mers in the first file.\n", l);

  if (l == NUMMERS) {
    fprintf(stderr, "WARNING:  Input chopped to %lu mers!\n", l);
    fflush(stderr);
  }

  speedCounter *cnt = new speedCounter(" %8.5f Mmers (%8.5f Mmers/sec)\r",
                                     1000000,
                                     100000);

  u64bit     f, r;

  while (T->nextMer()) {
    cnt->tick();

    f = T->theFMer();
    r = T->theRMer();

    for (u32bit j=0; j<l; j++) {
      if (F[j] == f) {
        fprintf(stdout, "Got a F match for mer=%s at onePos=%lu (in mers), twoPos=%llu (in bases)\n",
                T->theFMerString(), j, T->thePosition());
        fflush(stdout);
        fprintf(stderr, "Got a F match for mer=%s at onePos=%lu (in mers), twoPos=%llu (in bases)\n",
                T->theFMerString(), j, T->thePosition());
        fflush(stderr);
      }

      if (R[j] == r) {
        fprintf(stdout, "Got a R match for mer=%s at onePos=%lu (in mers), twoPos=%llu (in bases)\n",
                T->theFMerString(), j, T->thePosition());
        fflush(stdout);
        fprintf(stderr, "Got a R match for mer=%s at onePos=%lu (in mers), twoPos=%llu (in bases)\n",
                T->theFMerString(), j, T->thePosition());
        fflush(stderr);
      }
    }
  }
}
