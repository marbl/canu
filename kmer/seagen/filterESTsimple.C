#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"

//  A very simple filter.
//
//  Output the top 50 hits or all hits above 0.2, whichever is _smaller_.

#include "hitReader.H"

#define UNIQ_THRESH 50
#define QUAL_THRESH  0.2

int
main(int argc, char **argv) {

  if (argc == 1) {
    fprintf(stderr, "ESTmapper utility function -- not for human use.\n");
    exit(1);
  }

  hitReader    HR(argc);

  //  takes no args

  int arg = 1;
  while (arg < argc) {
    HR.addInputFile(argv[arg]);
    arg++;
  }

  while (HR.loadHits()) {

    HR.sortByCoverage();

    //  Output top 'UNIQ_THRESH' hits

    u32bit  max = UNIQ_THRESH;

    if (max >= HR.numHits())
      max = HR.numHits();

    for (u32bit i=0; i<max; i++)
      ahit_printASCII(&HR[i].a, stdout);

#if 0
    u32bit  count = 0;

    for (u32bit i=0; i < HR.numHits(); i++)
      if (QUAL_THRESH <= HR[i].coverage)
        count++;

    if ((count > 0) && (count < UNIQ_THRESH)) {
      //  Output all hits above QUAL_THRESH
      for (u32bit i=0; i < HR.numHits(); i++)
        if (QUAL_THRESH <= HR[i].coverage)
          ahit_printASCII(&HR[i].a, stdout);
    } else {
      //  Output top 'UNIQ_THRESH' hits
      for (u32bit i=0; i < UNIQ_THRESH; i++)
        ahit_printASCII(&HR[i].a, stdout);
    }
#endif

  }

  return(0);
}
