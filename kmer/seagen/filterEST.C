#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"

//  Define this if your hits have answers and you're curious about how
//  well the filter is performing.  The GOOD_* values decide if a
//  match is a true positive or false negative.
//
//#define WITH_ANSWERS
//#define GOOD_PERCENTID   95
//#define GOOD_COVERAGE    80

#include "hitReader.H"

//  XXX: Thread the filter!  Really cool!  Pretty neat hack!  Usual
//  thing, a thread to read hits, n threads to do filtering, and a
//  thread to write filtered hits.  Not trivial, but maybe a win.

//  Global statistics
//
u64bit  hitsSaved     = 0;
u64bit  hitsFiltered  = 0;
u64bit  hitsUnknown   = 0;
u64bit  hitsDiscarded = 0;
u64bit  hitsTotal     = 0;

u64bit  seqsSaved     = 0;
u64bit  seqsUnknown   = 0;
u64bit  seqsDiscarded = 0;
u64bit  seqsTotal     = 0;

int     filterTP      = 0;
int     filterFP      = 0;
int     filterFNfilt  = 0;  //  false negatives from filtering
int     filterFNunk   = 0;  //  false negatives from our failure to classify
int     filterFNdisc  = 0;  //  false negatives from lack of signal knee
int     filterTN      = 0;

int
main(int argc, char **argv) {

  if (argc == 1) {
    fprintf(stderr, "ESTmapper utility function -- not for human use.\n");
    exit(1);
  }

  hitReader    HR(argc);

  u32bit       uniqThresh         = 200;  //  Used to be 100
  u32bit       reptThresh         = 200;  //  Used to be 100
  double       qualityThresh      = 0.2;

  bool         conservativeFilter = false;

  FILE        *logFile            = stderr;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-uniquethreshold", 2) == 0) {
      uniqThresh  = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-repeatthreshold", 2) == 0) {
      reptThresh  = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-qualitythreshold", 2) == 0) {
      qualityThresh = atof(argv[++arg]);
    } else if (strncmp(argv[arg], "-conservative", 2) == 0) {
      conservativeFilter = true;
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

    //  Pick a quality threshold based on the input scores.  This is
    //  only used if there are more than uniqThresh hits for this
    //  EST.
    //
    //  The initial threshold is 20% of the range of hits, but we'll
    //  decrease this if it is bigger than the threshold supplied by
    //  the user.

    double qualityThreshUse = 0.2 * (HR.bestScore() - HR.worstScore()) + HR.worstScore();
    if (qualityThreshUse > qualityThresh)
      qualityThreshUse = qualityThresh;


    //  Pick a number of hits to save.  If the number of hits to
    //  save is zero, the sequence is flagged repetitive and no hits
    //  are output.

    bool    decided           = false;
    u32bit  hitsToSave        = 0;
    u32bit  count             = 0;
    double  largestDifference = 0.0;

    ////////////////////////////////////////
    //
    //  If we've got fewer hits than the repeat threshold, save ALL
    //  the hits.
    //
    if (HR.numHits() <= uniqThresh) {
      decided          = true;
      hitsToSave       = HR.numHits();
      qualityThreshUse = 0.0;

      //  Try being a little more aggressive.  Search for the last point
      //  where the score difference is more than 0.1 and use that for a
      //  limit.
      //
      if (!conservativeFilter) {
        HR.sortByCoverage();

        u32bit i = HR.numHits() - 1;
        while (i>9 && (HR[i-10].coverage - HR[i].coverage) < 0.1)
          i--;

        hitsToSave       = i + 1;
        qualityThreshUse = HR[i].coverage;
      }


      if (logFile) {
        if (!conservativeFilter)
          fprintf(logFile, u32bitFMT"] all "u32bitFMT" hits saved (aggressively filtered to "u32bitFMT" hits).\n", HR.iid(), HR.numHits(), hitsToSave);
        else
          fprintf(logFile, u32bitFMT"] all "u32bitFMT" hits saved.\n", HR.iid(), HR.numHits());
        fflush(logFile);
      }
    } else {


      ////////////////////////////////////////
      //
      //  Apply the first filter.  Count the number of hits with score
      //  more than the quality threshold picked above.  If the number
      //  to save is less than the repeat threshold, save all the hits
      //  above the threshold.
      //
      if (!conservativeFilter) {
        HR.sortByCoverage();

        //  Find the last hit that is above our hits to save
        //
        u32bit i = HR.numHits() - 1;

        while (i>0 && HR[i].coverage>qualityThreshUse)
          i--;

        while (i>9 && (HR[i-10].coverage - HR[i].coverage) < 0.1)
          i--;

        hitsToSave       = i + 1;
        qualityThreshUse = HR[i].coverage;

        if ((hitsToSave > 0) && (hitsToSave < reptThresh)) {
          decided    = true;

          if (logFile) {
            fprintf(logFile, ""u32bitFMT"] aggressively filtered "u32bitFMT" hits down to "u32bitFMT" hits using threshold %f\n", HR.iid(), HR.numHits(), hitsToSave, qualityThreshUse);
            fflush(logFile);
          }
        } else {
          hitsToSave = 0;
        }

      } else {
        for (u32bit i=0; i < HR.numHits(); i++)
          if (qualityThreshUse <= HR[i].coverage)
            count++;

        if (count <= reptThresh) {
          decided    = true;
          hitsToSave = HR.numHits();

          if (logFile) {
            fprintf(logFile, ""u32bitFMT"] filtered "u32bitFMT" hits down to "u32bitFMT" hits using threshold %f\n", HR.iid(), HR.numHits(), count, qualityThreshUse);
            fflush(logFile);
          }
        }
      }
    }


    ////////////////////////////////////////
    //
    //  If we're not decided here, the EST had too many "good" hits to
    //  be filtered by the threshold method.  Try a more sophisticated
    //  (confusing) method.
    //
    if (!decided) {

      //  Do some bookkeeping -- sort the hits and find the largest
      //  score difference.
      //
      HR.sortByCoverage();

      for (u32bit i=1; i < HR.numHits(); i++)
        if (largestDifference < (HR[i-1].coverage - HR[i].coverage))
          largestDifference = HR[i-1].coverage - HR[i].coverage;


      //////////////////////////////////////////////////
      //
      //  If all scores are about the same, it's either a repeat or a
      //  lot of spurious matches, depending on the level of signal.
      //
      if ((HR.bestScore() - HR.worstScore()) < 0.1) {
        decided    = true;
        hitsToSave = 0;

        if (logFile) {
          fprintf(logFile, ""u32bitFMT"] has uniform signal strength, no hits saved: "u32bitFMT" hits, best=%f, worst=%f, largestdiff=%f\n",
                  HR.iid(), HR.numHits(), HR.bestScore(), HR.worstScore(), largestDifference);
        }
      } else if (largestDifference < 0.1) {

        ////////////////////////////////////////
        //
        //  If the largest difference is below 10% coverage, then it's not
        //  clear how to pick a threshold.
        //
        decided    = true;
        hitsToSave = 0;

        if (logFile) {
          fprintf(logFile, ""u32bitFMT"] has no clear signal knee, no hits saved: "u32bitFMT" hits, best=%f, worst=%f, largestdiff=%f\n",
                  HR.iid(), HR.numHits(), HR.bestScore(), HR.worstScore(), largestDifference);
        }
      } else {

        ////////////////////////////////////////
        //
        //  Identify any spike near the start.  If we see a spike,
        //  save the first uniqThresh hits.
        //
        //  If the largest difference (which we guarantee to be >= 10%
        //  coverage here) is in the first uniqThresh hits, then we
        //  have a spike and we output uniqThresh hits.
        //
        //  To narrow the range more, we find the last spot where the
        //  difference in scores over 10 hits is > 0.1.  This is a
        //  generous heuristic.
        //

        u32bit  spikeFound = 0;

        for (u32bit i=1; i < uniqThresh; i++)
          if ((HR[i-1].coverage - HR[i].coverage) > 0.1) {
            decided    = true;
            spikeFound = i;
          }

        //  If we have found a spike, start at hit[uniqThresh], search
        //  backwards for the first point where the difference in
        //  scores across 10 hits is larger than 0.1

        if (decided) {
          hitsToSave = uniqThresh;

          for (u32bit i=uniqThresh-1; i > 9; i--)
            if ((HR[i-10].coverage - HR[i].coverage) > 0.1) {
              hitsToSave = i + 1;
              break;
            }

          qualityThreshUse = HR[hitsToSave].coverage;

          if (logFile) {
            fprintf(logFile, ""u32bitFMT"] spike at "u32bitFMT", "u32bitFMT" hits saved:  thresh=%f, "u32bitFMT" hits, best=%f, worst=%f, largestdiff=%f\n",
                    HR.iid(), spikeFound, hitsToSave, qualityThreshUse, HR.numHits(), HR.bestScore(), HR.worstScore(), largestDifference);
          }
        } else {
          if (logFile) {
            fprintf(logFile, ""u32bitFMT"] is an unclassified signal, no hits saved, "u32bitFMT" hits, best=%f, worst=%f, largestdiff=%f\n",
                    HR.iid(), HR.numHits(), HR.bestScore(), HR.worstScore(), largestDifference);
          }
        }
      }
    }

#ifdef WITH_ANSWERS
    int tp=0, tn=0, fn=0, fp=0;
#endif

    //  If we still haven't figure out what to do, then the EST
    //  is labeled a repeat.  Otherwise, write the (filtered)
    //  hits to the file.
    //
    if (!decided) {
      hitsUnknown += HR.numHits();
      seqsUnknown++;

#ifdef WITH_ANSWERS
      //  We've failed to classify all these hits, so anythig that looks good is a false negative
      //
      for (u32bit i=0; i < HR.numHits(); i++)
        if ((HR[i].mappedCoverage >= GOOD_COVERAGE) && (HR[i].mappedIdentity >= GOOD_PERCENTID)) {
          fprintf(logFile, "FAILUNKN hit=%3d id=%2d cv=%2d COV=%5.3f MUL=%5.3f: ", i, HR[i].mappedIdentity, HR[i].mappedCoverage, HR[i].coverage, HR[i].multiplicity);
          ahit_printASCII(&HR[i].a, logFile);
          filterFNunk++;
          fn++;
        } else {
          tn++;
        }
#endif

    } else if (hitsToSave == 0) {
      hitsDiscarded += HR.numHits();
      seqsDiscarded++;

#ifdef WITH_ANSWERS
      //  We've discarded all these hits, so anythig that looks good is a false negative
      //
      for (u32bit i=0; i < HR.numHits(); i++)
        if ((HR[i].mappedCoverage >= GOOD_COVERAGE) && (HR[i].mappedIdentity >= GOOD_PERCENTID)) {
          fprintf(logFile, "FAILDISC hit=%3d id=%2d cv=%2d COV=%5.3f MUL=%5.3f: ", i, HR[i].mappedIdentity, HR[i].mappedCoverage, HR[i].coverage, HR[i].multiplicity);
          ahit_printASCII(&HR[i].a, logFile);
          filterFNdisc++;
          fn++;
        } else {
          tn++;
        }
#endif

    } else {
      seqsSaved++;

      for (u32bit i=0; i < hitsToSave; i++) {
        if (qualityThreshUse <= HR[i].coverage) {
          hitsSaved++;
          ahit_printASCII(&HR[i].a, stdout);

#ifdef WITH_ANSWERS
          if ((HR[i].mappedCoverage >= 50) && (HR[i].mappedIdentity >= 95)) {
            tp++;
          } else {
            fp++;
          }
#endif
        } else {
          hitsFiltered++;
#ifdef WITH_ANSWERS
          //  Report hits that are false negatives -- these are probably an algorithmic error
          //
          if ((HR[i].mappedCoverage >= GOOD_COVERAGE) && (HR[i].mappedIdentity >= GOOD_PERCENTID)) {
            fprintf(logFile, "FAILFILT hit=%3d id=%2d cv=%2d COV=%5.3f MUL=%5.3f: ", i, HR[i].mappedIdentity, HR[i].mappedCoverage, HR[i].coverage, HR[i].multiplicity);
            ahit_printASCII(&HR[i].a, logFile);
            filterFNfilt++;
            fn++;
          } else {
            tn++;
          }
#endif
        }
      }

#ifdef WITH_ANSWERS
      //
      //  Look at the unfiltered stuff -- did we miss anything?
      //
      for (u32bit i=hitsToSave; i < HR.numHits(); i++) {
        //  Report hits that are false negatives
        //
        if ((HR[i].mappedCoverage >= GOOD_COVERAGE) && (HR[i].mappedIdentity >= GOOD_PERCENTID)) {
          fprintf(logFile, "FAILFILT hit=%3d id=%2d cv=%2d COV=%5.3f MUL=%5.3f: ", i, HR[i].mappedIdentity, HR[i].mappedCoverage, HR[i].coverage, HR[i].multiplicity);
          ahit_printASCII(&HR[i].a, logFile);
          filterFNfilt++;
          fn++;
        } else {
          tn++;
        }
      }
#endif

    }

#ifdef WITH_ANSWERS
    filterTP += tp;
    filterTN += tn;
    filterFP += fp;
#endif

    hitsTotal += HR.numHits();
    seqsTotal++;

#ifdef WITH_ANSWERS
    if (fn > 0) {
      fprintf(logFile, ""u32bitFMT"] FALSENEGATIVE  tp=%7d fp=%7d fn=%7d tn=%7d   hits saved:"u32bitFMTW(6) filtered:"u32bitFMTW(6) unknown:"u64bitFMTW(6)" discarded:"u64bitFMTW(6)" total:"u64bitFMTW(6)"\n",
              HR.iid(),
              tp, fp, fn, tn,
              hitsSaved,
              hitsFiltered,
              hitsUnknown,
              hitsDiscarded,
              hitsTotal);
    }
#endif


#ifdef WITH_ANSWERS
    if ((HR.iid() % 500) == 0) {
      fprintf(stderr, "%9u] tp=%7d fp=%7d fnfilt=%7d fnunkn=%7d fndisc=%7d tn=%7d  seqs saved:%6u unknown:"u32bitFMT" discarded:"u32bitFMTW(6) total:"u32bitFMT@(6)  hits saved "u64bitFMTW(8)"/"u64bitFMTW(8)" = %8.3f%%\r",
              HR.iid(),
              filterTP, filterFP, filterFNfilt, filterFNunk, filterFNdisc, filterTN,
              seqsSaved,
              seqsUnknown,
              seqsDiscarded,
              seqsTotal,
              hitsSaved, hitsTotal,
              100.0 * hitsSaved / hitsTotal);
      fflush(stderr);
    }
#else
    if ((HR.iid() % 500) == 0) {
      fprintf(stderr, "%9u] seqs saved:%6u unknown:"u32bitFMT" discarded:%6u total:%6u  hits saved "u64bitFMTW(8)"/"u64bitFMTW(8)" = %8.3f%%\r",
              HR.iid(),
              seqsSaved,
              seqsUnknown,
              seqsDiscarded,
              seqsTotal,
              hitsSaved, hitsTotal,
              100.0 * hitsSaved / hitsTotal);
      fflush(stderr);
    }
#endif
  }

  fclose(logFile);

  fprintf(stderr, "\n");

  return(0);
}
