#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"

#define MERSIZE 5

int
main(int argc, char **argv) {
  int     e     = 0;
  u32bit  nfull = 0;
  u32bit  npart = 0;
  char    str[1024];
  u32bit  i;

  merStream  *MS2 = new merStream(MERSIZE, new seqStream(argv[1], true));

#if 0
  while (MS2->nextMer()) {
    nfull++;
    fprintf(stdout, "nfull: %s\n", MS2->theFMer().merToString(str));
  }
  fprintf(stdout, "nfull: "u32bitFMT"\n", nfull);
#endif

  seqStore   *SS   = new seqStore(argv[1], new seqStream(argv[1], true));
  merStream   *MS1 = new merStream(MERSIZE, SS);

#if 0
  MS1->setRange(0, 200);
  while (MS1->nextMer()) {
    npart++;
    fprintf(stdout, "npart1: %s\n", MS1->theFMer().merToString(str));
  }
  fprintf(stdout, "npart: "u32bitFMT"\n", npart);
  
  MS1->setRange(10, 99999999);
  while (MS1->nextMer()) {
    npart++;
    fprintf(stdout, "npart2: %s\n", MS1->theFMer().merToString(str));
  }
  fprintf(stdout, "npart: "u32bitFMT"\n", npart);
#endif

  for (i=0; i<10000; i++) {
    MS1->setRange(i, i+2);
    while (MS1->nextMer())
      fprintf(stdout, "%d %s\n", i, MS1->theFMer().merToString(str));
  }


  return(e);
}

