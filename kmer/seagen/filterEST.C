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
uint32  hitsSaved     = 0;
uint32  hitsFiltered  = 0;
uint32  hitsUnknown   = 0;
uint32  hitsTotal     = 0;

uint32  seqsMapped    = 0;  //  Sequences that we mapped
uint32  seqsPartial   = 0;  //  Sequences that we mapped, but missed a few good matches
uint32  seqsMissed    = 0;  //  Sequences that we failed to map, but should have

uint32  filterTP      = 0;
uint32  filterFP      = 0;
uint32  filterFNfilt  = 0;  //  false negatives from filtering
uint32  filterFNunk   = 0;  //  false negatives from our failure to classify
uint32  filterTN      = 0;

uint32  goodPercentID = 94;
uint32  goodCoverage  = 50;

//  Command line options
//
uint32       uniqThresh         = 200;  //  Used to be 100
uint32       reptThresh         = 200;  //  Used to be 100

FILE        *logFile            = 0L;

//  Filter results -- thread unsafe!
//
//  bool   decided    -- true if the filter could decide on how to filter the hits
//  char   label      -- if decided, the name of the decider
//  uint32 hitsToSave -- the number of hits to save
//  double qualToSave -- the quality threshold to filter at
//
bool         decided;
const char  *label;
uint32       hitsToSave;
double       qualToSave;




void
report(uint32 iid
#ifdef WITH_ANSWERS
       ,
       uint32 filterTP,
       uint32 filterFP,
       uint32 filterFNfilt,
       uint32 filterFNunk,
       uint32 filterTN,
       uint32 seqsMapped,
       uint32 seqsPartial,
       uint32 seqsMissed
#endif
       ) {

  fprintf(stderr,
          uint32FMTW(9)"]"
#ifdef WITH_ANSWERS
          "  tp="uint32FMTW(7)" fp="uint32FMTW(7)" fnfilt="uint32FMTW(7)" fnunkn="uint32FMTW(7)" tn="uint32FMTW(7)
          "  yea:"uint32FMTW(7)" may:"uint32FMTW(7)" nay:"uint32FMTW(7)
#endif
          "  hits saved:"uint32FMTW(8)"/"uint32FMTW(8)" = %6.3f%%\r",
          iid,
#ifdef WITH_ANSWERS
          filterTP, filterFP, filterFNfilt, filterFNunk, filterTN,
          seqsMapped, seqsPartial, seqsMissed,
#endif
          hitsSaved, hitsTotal,
          100.0 * hitsSaved / hitsTotal);
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
      for (uint32 i=0; i < HR.numHits(); i++)
        if ((HR[i].mappedCoverage >= goodCoverage) && (HR[i].mappedIdentity >= goodPercentID)) {
          if (logFile) {
            fprintf(logFile, "FAILUNKN hit="uint32FMTW(3)" id=%2d cv=%2d COV=%5.3f MUL=%5.3f: ", i, HR[i].mappedIdentity, HR[i].mappedCoverage, HR[i].coverage, HR[i].multiplicity);
            ahit_printASCII(&HR[i].a, logFile);
          }
          filterFNunk++;
          fn++;
        } else {
          tn++;
        }
#endif

    } else {
      for (uint32 i=0; i < HR.numHits(); i++) {
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
              fprintf(logFile, "FAILFILT hit="uint32FMTW(3)" id=%2d cv=%2d COV=%5.3f MUL=%5.3f: ", i, HR[i].mappedIdentity, HR[i].mappedCoverage, HR[i].coverage, HR[i].multiplicity);
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
        fprintf(logFile, uint32FMT"] %sFALSENEGATIVE %10.10s  tp="uint32FMTW(7)" fp="uint32FMTW(7)" fn="uint32FMTW(7)" tn="uint32FMTW(7)"\n",
                HR.iid(),
                (tp > 0) ? "partial" : "fatal",
                label,
                tp, fp, fn, tn);
#endif


    if ((HR.iid() % 500) == 0) {
      report(HR.iid()
#ifdef WITH_ANSWERS
             ,
             filterTP, filterFP, filterFNfilt, filterFNunk, filterTN,
             seqsMapped, seqsPartial, seqsMissed
#endif
             );
      fflush(stderr);
    }
  }

  if (logFile)
    fclose(logFile);


  report(HR.iid()
#ifdef WITH_ANSWERS
         ,
         filterTP, filterFP, filterFNfilt, filterFNunk, filterTN,
         seqsMapped, seqsPartial, seqsMissed
#endif
         );
  fprintf(stderr, "\n");

  return(0);
}
