#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"

//  $Id$

#define MAX_ESTS    (16 * 1024 * 1024)

typedef struct {
  bool     stillMore;
  bool     reload;
  FILE    *file;
  char     b[1024];
  aHit     a;
  bool     isBINARY;
} hitFile_s;

typedef struct {
  aHit     a;
  u32bit   estid;
  float    coverage;
  float    multiplicity;
} hit_s;


void
loadHit(hitFile_s *HF) {

  if (HF->stillMore && HF->reload) {
    HF->reload = false;

    if (HF->isBINARY) {
      ahit_readBinary(&HF->a, HF->file);
    } else {
      fgets(HF->b, 1024, HF->file);
      ahit_parseString(&HF->a, HF->b);
    }

    if (feof(HF->file))
      HF->stillMore = false;
  }
}


int
hitCompare(const void *a, const void *b) {
  const hit_s  *A = (const hit_s *)a;
  const hit_s  *B = (const hit_s *)b;

  if (A->coverage > B->coverage)
    return(-1);
  else
    return(A->coverage < B->coverage);
}


FILE *openOutputFile(char *n) {
  FILE *f = 0L;
  if (n) {
    errno = 0;
    while ((f = fopen(n, "w")) == 0L) {
      fprintf(stderr, "filterMRNA: ERROR: couldn't open '%s' for writing.\n%s\n", n, strerror(errno));
      sleep(1);
    }
  }
  return(f);
}


int
main(int argc, char **argv) {

  if (argc < 2) {
    fprintf(stderr, "ESTmapper utility function -- not for human use.\n");
    exit(1);
  }

  u32bit      *hitCounts    = new u32bit [MAX_ESTS];
  u32bit       numESTs      = 0;
  hitFile_s   *HF           = new hitFile_s [argc];
  u32bit       numFiles     = 0;

  double L  = 0.2;
  double H  = 0.6;
  double V  = 0.7;
  double M  = 0.3;
  double MC = 0.2;
  u32bit ML = 150;

  bool beVerbose = false;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-counts", 2) == 0) {
      errno = 0;
      FILE *F = fopen(argv[++arg], "r");
      if (F == 0L) {
        fprintf(stderr, "filterMRNA: ERROR: couldn't open '%s' for reading.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
      while (!feof(F)) {
#ifdef TRUE64BIT
        fscanf(F, " %u", hitCounts + numESTs);
#else
        fscanf(F, " %lu", hitCounts + numESTs);
#endif
        if (!feof(F))
          numESTs++;
      }
      fclose(F);
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
    } else if (strncmp(argv[arg], "-verbose", 2) == 0) {
      beVerbose = true;
    } else {
      errno = 0;

      HF[numFiles].stillMore = true;
      HF[numFiles].reload    = true;
      HF[numFiles].file      = fopen(argv[arg], "r");

      if (HF[numFiles].file == 0L) {
        fprintf(stderr, "filterMRNA: ERROR: couldn't open '%s' for reading.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }

      //  Binary or ASCII input?
      //
      char x = (char)fgetc(HF[numFiles].file);
      ungetc(x, HF[numFiles].file);

      HF[numFiles].isBINARY = (x != '-');

      if (HF[numFiles].file == 0L) {
        fprintf(stderr, "filterMRNA: ERROR opening '%s'\n", argv[arg]);
      } else {
        numFiles++;
      }
    }

    arg++;
  }

  //
  //  Check the args
  //
  if (beVerbose) {
    fprintf(stderr, "Filtering with:\n");
    fprintf(stderr, "  score difference of %4.2f or less -> 100.0%% of best score\n", L);
    fprintf(stderr, "  score difference of %4.2f or more -> %5.1f%% of best score\n", H, 100*V);
    fprintf(stderr, "  scores at least %4.2f are always output\n", M);
#ifdef TRUE64BIT
    fprintf(stderr, "  scores at least %4.2f AND at least %u bases covered are always output\n", MC, ML);
#else
    fprintf(stderr, "  scores at least %4.2f AND at least %lu bases covered are always output\n", MC, ML);
#endif
  }

  //  While there is input, merge
  //
  bool  keepGoing = true;

  while (keepGoing) {
    keepGoing = false;

    //  Load hits, if needed.
    //
    for (u32bit i=0; i<numFiles; i++)
      loadHit(HF+i);

    //  See if we're done.
    //
    for (u32bit i=0; i<numFiles; i++)
      keepGoing |= HF[i].stillMore;

    if (keepGoing) {

      //  Find the lowest ESTid
      //
      u32bit estOfInterest = 1 << 30;
      for (u32bit i=0; i<numFiles; i++)
        if ((HF[i].stillMore) && (estOfInterest > HF[i].a._qsIdx))
          estOfInterest = HF[i].a._qsIdx;

      //  Create a list of hits.  We know the size it's supposed to be,
      //  so we can allocate exactly that space.
      //
      //  For each file, save all hits with the estOfInterest to a list.
      //
      u32bit       listLen = 0;
      u32bit       listMax = hitCounts[estOfInterest];
      hit_s       *list    = new hit_s [listMax];

      for (u32bit i=0; i<numFiles; i++) {
        while ((HF[i].stillMore) && (HF[i].a._qsIdx == estOfInterest)) {
          if (listLen >= listMax) {
            fprintf(stderr, "ERROR:  Too many hits -- hitCounts might be invalid!\n");
            exit(1);
          }

          memcpy(&list[listLen].a, &HF[i].a, sizeof(aHit));

          list[listLen].coverage     = (float)HF[i].a._covered / (float)HF[i].a._numMers;
          list[listLen].multiplicity = (float)HF[i].a._matched / (float)HF[i].a._covered;

          //  aHit->_covered is in bases, but aHit->_numMers is the
          //  number of mers.  Possible for coverage to be > 1.0.
          //
          if (list[listLen].coverage > 1.0)
            list[listLen].coverage = 1.0;

          listLen++;
          HF[i].reload = true;
          loadHit(HF+i);
        }
      }

      //  Sort the hits by score (coverage), in decreasing order.
      //
      qsort(list, listLen, sizeof(hit_s), hitCompare);

      double h = list[0].coverage - list[listLen-1].coverage;
      double p = 0.0;

      if (h <= L)    p = 1.0;
      if (h >= H)    p = V;
      if (p == 0.0)  p = 1.0 - (1.0 - V) * (h - L) / (H - L);

      //  check p; it should be between V and 1.0
      if ((p > 1.0) || (p < V)) {
        fprintf(stderr, "error in p; p=%f\n", p);
      }

      //  Output the top p% hits, by score.
      //
      double cutL = list[0].coverage - p * h;

      if (cutL > M)
        cutL = M;

      for (u32bit i=0; i < listLen; i++) {
        if ((cutL <= list[i].coverage) &&
            ((MC <= list[i].coverage) || (ML <= list[i].a._covered))) {
          ahit_printASCII(&list[i].a, stdout);
        }
      }

      delete [] list;
    }
  }

  return(0);
}
