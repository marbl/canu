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


int
main(int argc, char **argv) {

  if (argc == 1) {
    fprintf(stderr, "ESTmapper utility function -- not for human use.\n");
    exit(1);
  }

  hitReader    HR(argc);

  u32bit       uniqThresh         = 50;
  double       qualityThresh      = 0.2;

  //  takes no args

  int arg = 1;
  while (arg < argc) {
    HR.addInputFile(argv[arg]);
    arg++;
  }

  while (HR.loadHits()) {
    u32bit  count = 0;

    HR.sortByCoverage();

    //  Output top 'uniqThresh' hits

    u32bit  max = uniqThresh;

    if (max >= HR.numHits())
      max = HR.numHits();

    for (u32bit i=0; i<max; i++)
      ahit_printASCII(&HR[i].a, stdout);

#if 0
    for (u32bit i=0; i < HR.numHits(); i++)
      if (qualityThresh <= HR[i].coverage)
        count++;

    if ((count > 0) && (count < uniqThresh)) {
      //  Output all hits above qualityThresh
      for (u32bit i=0; i < HR.numHits(); i++)
        if (qualityThresh <= HR[i].coverage)
          ahit_printASCII(&HR[i].a, stdout);
    } else {
      //  Output top 'uniqThresh' hits
      for (u32bit i=0; i < uniqThresh; i++)
        ahit_printASCII(&HR[i].a, stdout);
    }
#endif

  }

  return(0);
}
