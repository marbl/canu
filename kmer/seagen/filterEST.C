#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"
#include "hitReader.H"

//  XXX: Thread the filter!  Really cool!  Pretty neat hack!  Usual
//  thing, a thread to read hits, n threads to do filtering, and a
//  thread to write filtered hits.  Not trivial, but maybe a win.

//  Global statistics
//
u32bit  hitsSaved     = 0;
u32bit  hitsFiltered  = 0;
u32bit  hitsUnknown   = 0;
u32bit  hitsTotal     = 0;

u32bit  seqsMapped    = 0;  //  Sequences that we mapped
u32bit  seqsPartial   = 0;  //  Sequences that we mapped, but missed a few good matches
u32bit  seqsMissed    = 0;  //  Sequences that we failed to map, but should have

u32bit  filterTP      = 0;
u32bit  filterFP      = 0;
u32bit  filterFNfilt  = 0;  //  false negatives from filtering
u32bit  filterFNunk   = 0;  //  false negatives from our failure to classify
u32bit  filterTN      = 0;

u32bit  goodPercentID = 94;
u32bit  goodCoverage  = 50;

//  Command line options
//
u32bit       uniqThresh         = 200;  //  Used to be 100
u32bit       reptThresh         = 200;  //  Used to be 100

FILE        *logFile            = 0L;

//  Filter results -- thread unsafe!
//
//  bool   decided    -- true if the filter could decide on how to filter the hits
//  char   label      -- if decided, the name of the decider
//  u32bit hitsToSave -- the number of hits to save
//  double qualToSave -- the quality threshold to filter at
//
bool         decided;
char        *label;
u32bit       hitsToSave;
double       qualToSave;




void
report(u32bit iid,
#ifdef WITH_ANSWERS
       u32bit filterTP,
       u32bit filterFP,
       u32bit filterFNfilt,
       u32bit filterFNunk,
       u32bit filterTN,
       u32bit seqsMapped,
       u32bit seqsPartial,
       u32bit seqsMissed,
#endif
       u32bit hitsSaved,
       u32bit hitsTotal,
       double perc) {

  fprintf(stderr,
          u32bitFMTW(9)"]"
#ifdef WITH_ANSWERS
          "  tp="u32bitFMTW(7)" fp="u32bitFMTW(7)" fnfilt="u32bitFMTW(7)" fnunkn="u32bitFMTW(7)" tn="u32bitFMTW(7)
          "  yea:"u32bitFMTW(7)" may:"u32bitFMTW(7)" nay:"u32bitFMTW(7)
#endif
          "  hits saved:"u32bitFMTW(8)"/"u32bitFMTW(8)" = %6.3f%%\r",
          iid,
#ifdef WITH_ANSWERS
          filterTP, filterFP, filterFNfilt, filterFNunk, filterTN,
          seqsMapped, seqsPartial, seqsMissed,
#endif
          hitsSaved, hitsTotal,
          perc);
}



void
complicatedFilter(hitReader &HR);



//  The simple filter just returns the top uniqThresh hits
//
void
simpleFilter(hitReader &HR) {

  decided    = true;
  label      = "simple";
  qualToSave = 0.0;
  hitsToSave = HR.numHits();

  if (HR.numHits() > uniqThresh)
    hitsToSave = uniqThresh;
}




int
main(int argc, char **argv) {

  if (argc == 1) {
    fprintf(stderr, "ESTmapper utility function -- not for human use.\n");
    exit(1);
  }

  hitReader    HR(argc);

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-uniquethreshold", 2) == 0) {
      uniqThresh  = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-repeatthreshold", 2) == 0) {
      reptThresh  = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-log", 2) == 0) {
      ++arg;
      errno = 0;
      logFile = fopen(argv[arg], "w");
      if (errno) {
        fprintf(stderr, "filterEST: ERROR: couldn't open logFile '%s' for writing.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
    } else {
      HR.addInputFile(argv[arg]);
    }

    arg++;
  }

  while (HR.loadHits()) {

    //  Not every filter we can think of needs the hits sorted, but
    //  it's nice to guarantee they are sorted.
    //
    HR.sortByCoverage();


    //simpleFilter(HR);
    complicatedFilter(HR);


#ifdef WITH_ANSWERS
    int tp=0, tn=0, fn=0, fp=0;
#endif

    //  If we still haven't figured out what to do, then the EST is
    //  labeled a repeat.  Otherwise, write the (filtered) hits to the
    //  file.
    //
    if (!decided) {
      hitsUnknown += HR.numHits();

#ifdef WITH_ANSWERS
      //  We've failed to classify all these hits, so anythig that looks good is a false negative
      //
      for (u32bit i=0; i < HR.numHits(); i++)
        if ((HR[i].mappedCoverage >= goodCoverage) && (HR[i].mappedIdentity >= goodPercentID)) {
          if (logFile) {
            fprintf(logFile, "FAILUNKN hit="u32bitFMTW(3)" id=%2d cv=%2d COV=%5.3f MUL=%5.3f: ", i, HR[i].mappedIdentity, HR[i].mappedCoverage, HR[i].coverage, HR[i].multiplicity);
            ahit_printASCII(&HR[i].a, logFile);
          }
          filterFNunk++;
          fn++;
        } else {
          tn++;
        }
#endif

    } else {
      for (u32bit i=0; i < HR.numHits(); i++) {
        if ((i < hitsToSave) && (qualToSave <= HR[i].coverage)) {
          hitsSaved++;
          ahit_printASCII(&HR[i].a, stdout);

#ifdef WITH_ANSWERS
          if ((HR[i].mappedCoverage >= goodCoverage) && (HR[i].mappedIdentity >= goodPercentID)) {
            tp++;
          } else {
            fp++;
          }
#endif
        } else if (HR[i].a._merged) {
          //  We merged this hit, so scores are incorrect.  Give it
          //  the benefit of the doubt and report it.
          //
          hitsSaved++;
          ahit_printASCII(&HR[i].a, stdout);
        } else {
          hitsFiltered++;
#ifdef WITH_ANSWERS
          //  Report hits that are false negatives
          //
          if ((HR[i].mappedCoverage >= goodCoverage) && (HR[i].mappedIdentity >= goodPercentID)) {
            if (logFile) {
              fprintf(logFile, "FAILFILT hit="u32bitFMTW(3)" id=%2d cv=%2d COV=%5.3f MUL=%5.3f: ", i, HR[i].mappedIdentity, HR[i].mappedCoverage, HR[i].coverage, HR[i].multiplicity);
              ahit_printASCII(&HR[i].a, logFile);
            }
            filterFNfilt++;
            fn++;
          } else {
            tn++;
          }
#endif
        }
      }
    }

#ifdef WITH_ANSWERS
    filterTP += tp;
    filterTN += tn;
    filterFP += fp;

    if (tp > 0)
      seqsMapped++;

    if (fn > 0)
      seqsPartial++;

    if ((tp == 0) && (fn > 0))
      seqsMissed++;
#endif

    hitsTotal += HR.numHits();

#ifdef WITH_ANSWERS
    //  Report if we saw falsenegatives (we should have printed FAIL into the log, too)
    //
    if (fn > 0)
      if (logFile)
        fprintf(logFile, u32bitFMT"] %sFALSENEGATIVE %10.10s  tp="u32bitFMTW(7)" fp="u32bitFMTW(7)" fn="u32bitFMTW(7)" tn="u32bitFMTW(7)"\n",
                HR.iid(),
                (tp > 0) ? "partial" : "fatal",
                label,
                tp, fp, fn, tn);
#endif


    if ((HR.iid() % 500) == 0) {
      report(HR.iid(),
#ifdef WITH_ANSWERS
             filterTP, filterFP, filterFNfilt, filterFNunk, filterTN,
             seqsMapped, seqsPartial, seqsMissed,
#endif
             hitsSaved, hitsTotal,
             100.0 * hitsSaved / hitsTotal);
      fflush(stderr);
    }
  }

  if (logFile)
    fclose(logFile);


  report(HR.iid(),
#ifdef WITH_ANSWERS
         filterTP, filterFP, filterFNfilt, filterFNunk, filterTN,
         seqsMapped, seqsPartial, seqsMissed,
#endif
         hitsSaved, hitsTotal,
         100.0 * hitsSaved / hitsTotal);
  fprintf(stderr, "\n");

  return(0);
}
