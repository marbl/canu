#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"
#include "hitReader.H"


extern u32bit       uniqThresh;
extern u32bit       reptThresh;
extern FILE        *logFile;
extern bool         decided;
extern char        *label;
extern u32bit       hitsToSave;
extern double       qualToSave;


void
complicatedFilter_1_unique(hitReader &HR) {

  if (HR.numHits() <= uniqThresh) {
    decided          = true;
    label            = "unique";
    hitsToSave       = HR.numHits();
    qualToSave       = 0.0;

    //  Try being a little more aggressive.  Search for the last
    //  point where the score difference across 10 hits is more than
    //  0.1 and use that for a limit.


    //  On the 100k test set #1 (ESTmapper paper, 1 Oct 2004) this modification
    //  results in:
    //    tp=106564 fp=1255487 fn=56705 tn=52816595
    //
    //  compared to saving all hits:
    //    tp=106579 fp=1914659 fn=56690 tn=52157423
    //
    //  That is, we lost 15 true matches and didn't polish 660,000
    //  matches -- 1.21% of the total, but 50% of what we actually
    //  need to polish.
    //
    u32bit i = HR.numHits() - 1;
    while ((i >= 10) && ((HR[i-10].coverage - HR[i].coverage) < 0.1))
      i--;

    hitsToSave = HR.numHits();
    qualToSave = HR[i].coverage;

    //  Save all hits with this coverage score!  This isn't really needed, but it
    //  makes the log message correct.
    //
    while ((i < HR.numHits()) && (qualToSave == HR[i].coverage))
      i++;

    if (logFile)
      fprintf(logFile, u32bitFMT"] unique: aggressively filtered to "u32bitFMT" hits out of "u32bitFMT" hits.\n",
              HR.iid(), i, HR.numHits());
  }
}





void
complicatedFilter_2_knee(hitReader &HR) {

  decided    = false;
  hitsToSave = 0;
  qualToSave = 0.0;

  //  Apply the same filter as used in #1 (the aggressive part), and accept
  //  it if the number of hits saved is below some threshold.

  u32bit i = HR.numHits() - 1;
  while ((i >= 10) && ((HR[i-10].coverage - HR[i].coverage) < 0.1))
    i--;

  //  If i==9, then we failed to find a knee, and we fail this filter
  //
  if (i < 10)
    return;

  hitsToSave = HR.numHits();
  qualToSave = HR[i].coverage;

  //  Save all hits with this coverage score!
  //
  while ((i < HR.numHits()) && (qualToSave == HR[i].coverage))
    i++;

  if (i <= reptThresh) {
    decided = true;
    label   = "knee";

    if (logFile)
      fprintf(logFile, u32bitFMT"] knee: filtered "u32bitFMT" hits down to "u32bitFMT" hits using threshold %f\n",
              HR.iid(), HR.numHits(), i, qualToSave);
  }
}



//  If all scores are about the same, it's either a repeat or a
//  lot of spurious matches, depending on the level of signal.
//
void
complicatedFilter_3_uniform(hitReader &HR) {

  decided    = false;
  hitsToSave = 0;
  qualToSave = 0.0;

  if ((HR.bestScore() - HR.worstScore()) < 0.1) {
    decided    = true;
    label      = "uniform";
    hitsToSave = reptThresh;
    qualToSave = 0.0;

    if (logFile) {
      fprintf(logFile, u32bitFMT"] uniform: uniform signal strength, saving the first "u32bitFMT" hits out of "u32bitFMT" hits, best=%f, worst=%f\n",
              HR.iid(), hitsToSave, HR.numHits(), HR.bestScore(), HR.worstScore());
    }
  }
}




//  If we're not decided here, the EST had too many "good" hits to
//  be filtered by the threshold method.  Try a more sophisticated
//  (confusing) method.
//
void
complicatedFilter_4_largestdifference(hitReader &HR) {

  decided    = false;
  hitsToSave = 0;
  qualToSave = 0.0;

  double  largestDifference = 0.0;
  for (u32bit i=1; i < HR.numHits(); i++)
    if (largestDifference < (HR[i-1].coverage - HR[i].coverage))
      largestDifference = HR[i-1].coverage - HR[i].coverage;

  //  If the largest difference is below 10% coverage, then it's not
  //  clear how to pick a threshold and we just save a bunch of hits.
  //
  if (largestDifference < 0.1) {
    decided    = true;
    label      = "diff";
    hitsToSave = reptThresh;
    qualToSave = 0.0;

    if (logFile)
      fprintf(logFile, u32bitFMT"] diff: has no clear signal knee, saving the first "u32bitFMT" hits out of "u32bitFMT" hits, best=%f, worst=%f, largestdiff=%f\n",
              HR.iid(), hitsToSave, HR.numHits(), HR.bestScore(), HR.worstScore(), largestDifference);
  }
}



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
void
complicatedFilter_5_spikes(hitReader &HR) {

  decided    = false;
  hitsToSave = 0;
  qualToSave = 0.0;

  u32bit  spikeFound = 0;
  for (u32bit i=1; i < uniqThresh; i++)
    if ((HR[i-1].coverage - HR[i].coverage) > 0.1)
      spikeFound = i;

  //  If we have found a spike, start at hit[uniqThresh], search
  //  backwards for the first point where the difference in
  //  scores across 10 hits is larger than 0.1

  if (spikeFound) {
    decided    = true;
    label      = "spike";
    hitsToSave = uniqThresh;
    qualToSave = 0.0;

    for (u32bit i=uniqThresh-1; i > 9; i--)
      if ((HR[i-10].coverage - HR[i].coverage) > 0.1) {
        hitsToSave = i + 1;
        break;
      }

    qualToSave = HR[hitsToSave].coverage;

    if (logFile)
      fprintf(logFile, u32bitFMT"] spike: at "u32bitFMT", "u32bitFMT" hits saved:  thresh=%f, "u32bitFMT" hits, best=%f, worst=%f\n",
              HR.iid(), spikeFound, hitsToSave, qualToSave, HR.numHits(), HR.bestScore(), HR.worstScore());
  }
}




void
complicatedFilter(hitReader &HR) {
  decided    = false;
  label      = "NOLABELERROR";
  qualToSave = 1.0;
  hitsToSave = 0;

  complicatedFilter_1_unique(HR);

  if (decided)
    return;

  complicatedFilter_2_knee(HR);

  if (decided)
    return;

  complicatedFilter_3_uniform(HR);

  if (decided)
    return;

  complicatedFilter_4_largestdifference(HR);

  if (decided)
    return;

  complicatedFilter_5_spikes(HR);

  if (decided)
    return;

  decided    = true;
  label      = "unknown";
  hitsToSave = reptThresh;
  qualToSave = 0.0;

  if (hitsToSave > HR.numHits())
    hitsToSave = HR.numHits();

  if (logFile)
    fprintf(logFile, u32bitFMT"] is an unclassified signal, "u32bitFMT" hits saved out of "u32bitFMT" hits, best=%f, worst=%f\n",
            HR.iid(), hitsToSave, HR.numHits(), HR.bestScore(), HR.worstScore());
}
