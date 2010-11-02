#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"
#include "hitReader.H"


int
main(int argc, char **argv) {

  if (argc < 2)
    fprintf(stderr, "ESTmapper utility function -- not for human use.\n"), exit(1);

  hitReader    HR(argc);
  double       L  = 0.2;
  double       H  = 0.6;
  double       V  = 0.7;
  double       M  = 0.3;
  double       MC = 0.2;
  u32bit       ML = 150;
  bool         beVerbose = false;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
    } else if (strcmp(argv[arg], "-l") == 0) {
      L = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-h") == 0) {
      H = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-v") == 0) {
      V = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-m") == 0) {
      M = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-mc") == 0) {
      MC = atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-ml") == 0) {
      ML = atoi(argv[++arg]);
    } else {
      HR.addInputFile(argv[arg]);
    }

    arg++;
  }

  if (beVerbose) {
    fprintf(stderr, "Filtering with:\n");
    fprintf(stderr, "  score difference of %4.2f or less -> 100.0%% of best score\n", L);
    fprintf(stderr, "  score difference of %4.2f or more -> %5.1f%% of best score\n", H, 100*V);
    fprintf(stderr, "  scores at least %4.2f are always output\n", M);
    fprintf(stderr, "  scores at least %4.2f AND at least "u32bitFMT" bases covered are always output\n", MC, ML);
  }

  while (HR.loadHits()) {
    HR.sortByCoverage();

    double  hiCov = HR[0].coverage;
    double  loCov = HR[0].coverage;
    for (u32bit i=0; i < HR.numHits(); i++)
      if ((HR[i].a._merged == false) && (loCov > HR[i].coverage))
        loCov = HR[i].coverage;
        
    double h = hiCov - loCov;
    double p = 0.0;

    if (h <= L)    p = 1.0;
    if (h >= H)    p = V;
    if (p == 0.0)  p = 1.0 - (1.0 - V) * (h - L) / (H - L);

    //  check p; it should be between V and 1.0
    if ((p > 1.0) || (p < V))
      fprintf(stderr, "error in p; p=%f\n", p);

    //  Output the top p% hits, by score.

    double cutL = HR[0].coverage - p * h;
    if (cutL > M)
      cutL = M;

    //  Save the hit if it has good coverage and it's either above
    //  the minimum coverage or long.  Also blindly save merged
    //  hits.
    //
    for (u32bit i=0; i < HR.numHits(); i++)
      if (((cutL <= HR[i].coverage) && ((MC <= HR[i].coverage) ||
                                        (ML <= HR[i].a._covered))) ||
          (HR[i].a._merged))
        ahit_printASCII(&HR[i].a, stdout);
  }

  return(0);
}
