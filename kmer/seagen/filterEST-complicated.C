#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"
#include "hitReader.H"


extern uint32       uniqThresh;
extern uint32       reptThresh;
extern FILE        *logFile;
extern bool         decided;
extern const char  *label;
extern uint32       hitsToSave;
extern double       qualToSave;

double difference = 0.1;

void
complicatedFilter_1_unique(hitReader &HR) {

  if (HR.numHits() <= uniqThresh) {
    decided          = true;
    label            = "unique";
    hitsToSave       = HR.numHits();
    qualToSave       = 0.0;

    //  Try being a little more aggressive.  Search for the last
    //  point where the score difference across 10 hits is more than
    //  difference and use that for a limit.


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
    uint32 i = HR.numHits() - 1;
    while ((i >= 10) && ((HR[i-10].coverage - HR[i].coverage) < difference))
      i--;

    hitsToSave = HR.numHits();
    qualToSave = HR[i].coverage;

#if 0
    //  Take the middle hit, not the end.  This doesn't hurt too much
    //  (20 matches out of 100,000 ESTs, and we missed one EST
    //  completely) but only gains us 0.08% additional filtering.
    if (i >= 15)
      qualToSave = HR[i-5].coverage;
#endif

    //  Save all hits with this coverage score!  This isn't really needed, but it
    //  makes the log message correct.
    //
    while ((i < HR.numHits()) && (qualToSave == HR[i].coverage))
      i++;

    if (logFile)
      fprintf(logFile, uint32FMT"] unique: aggressively filtered to "uint32FMT" hits out of "uint32FMT" hits.\n",
              HR.iid(), i, HR.numHits());
  }
}





void
complicatedFilter_2_knee(hitReader &HR) {

  decided    = false;
  hitsToSave = 0;
  qualToSave = 0.0;

#if 0
  decided    = true;
  hitsToSave = 0;
  qualToSave = 1.1;
  return;
#endif

  //  Apply the same filter as used in #1 (the aggressive part), and accept
  //  it if the number of hits saved is below some threshold.

  uint32 i = HR.numHits() - 1;
  while ((i >= 10) && ((HR[i-10].coverage - HR[i].coverage) < difference))
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

  if (i <= uniqThresh) {
    decided = true;
    label   = "knee";

    if (logFile)
      fprintf(logFile, uint32FMT"] knee: filtered "uint32FMT" hits down to "uint32FMT" hits using threshold %f\n",
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

  if ((HR.bestScore() - HR.worstScore()) < difference) {
    decided    = true;
    label      = "uniform";
    hitsToSave = reptThresh;
    qualToSave = 0.0;

    if (logFile) {
      fprintf(logFile, uint32FMT"] uniform: uniform signal strength, saving the first "uint32FMT" hits out of "uint32FMT" hits, best=%f, worst=%f\n",
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
  for (uint32 i=1; i < HR.numHits(); i++)
    if (largestDifference < (HR[i-1].coverage - HR[i].coverage))
      largestDifference = HR[i-1].coverage - HR[i].coverage;

  //  If the largest difference is below 10% coverage, then it's not
  //  clear how to pick a threshold and we just save a bunch of hits.
  //
  if (largestDifference < difference) {
    decided    = true;
    label      = "diff";
    hitsToSave = reptThresh;
    qualToSave = 0.0;

    if (logFile)
      fprintf(logFile, uint32FMT"] diff: has no clear signal knee, saving the first "uint32FMT" hits out of "uint32FMT" hits, best=%f, worst=%f, largestdiff=%f\n",
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
//  difference in scores over 10 hits is > difference.  This is a
//  generous heuristic.
//
void
complicatedFilter_5_spikes(hitReader &HR) {

  decided    = false;
  hitsToSave = 0;
  qualToSave = 0.0;

  uint32  spikeFound = 0;
  for (uint32 i=1; i < uniqThresh; i++)
    if ((HR[i-1].coverage - HR[i].coverage) > difference)
      spikeFound = i;

  //  If we have found a spike, start at hit[uniqThresh], search
  //  backwards for the first point where the difference in
  //  scores across 10 hits is larger than difference
  //
  //  Seems like a NOP, but it loosens things up a bit.  Consider a
  //  spike between hits 3 and 4, but 1=2=3 and 4=5=6=7=8=9=10=11.  We
  //  find a spike, then find a nice place to cut it.  If we never
  //  find a nice place, we save the top uniqThresh hits.

  if (spikeFound) {
    decided    = true;
    label      = "spike";
    hitsToSave = uniqThresh;
    qualToSave = 0.0;

    for (uint32 i=uniqThresh-1; i > 9; i--)
      if ((HR[i-10].coverage - HR[i].coverage) > difference) {
        hitsToSave = i + 1;
        break;
      }

    qualToSave = HR[hitsToSave].coverage;

    if (logFile)
      fprintf(logFile, uint32FMT"] spike: at "uint32FMT", "uint32FMT" hits saved:  thresh=%f, "uint32FMT" hits, best=%f, worst=%f\n",
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
    fprintf(logFile, uint32FMT"] is an unclassified signal, "uint32FMT" hits saved out of "uint32FMT" hits, best=%f, worst=%f\n",
            HR.iid(), hitsToSave, HR.numHits(), HR.bestScore(), HR.worstScore());
}
