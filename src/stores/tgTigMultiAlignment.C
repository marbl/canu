

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"
#include "aligners.H"
#include "unitigConsensus.H"


bool
generateMultiAlignment(tgTig     *tig,
                       gkStore   *gkpStore,
                       uint32    *failed,
                       double     errorRate,     //  AS_CNS_ERROR_RATE
                       double     errorRateMax,
                       uint32     minOverlap) {  //  AS_OVERLAP_MIN_LEN
  //double             origErate          = errorRate;
  //uint32             origLen            = minOverlap;
  bool               failuresToFix      = false;
  unitigConsensus   *uc                 = NULL;

  //origErate          = errorRate;
  uc                 = new unitigConsensus(tig, errorRate, minOverlap);

  if (uc->initialize(gkpStore, failed) == FALSE) {
    fprintf(stderr, "generateMultiAlignment()--  Failed to initialize for tig %u with %u children\n", tig->tigID(), tig->numberOfChildren());
    goto returnFailure;
  }

  while (uc->moreFragments()) {

#ifdef SKIP_CONTAINS
    if (uc->fraglist[uc->tiid].contained != 0) {
      //  skip contained for now
      fprintf(stderr, "SKIP CONTAINED %u\n", uc->fraglist[uc->tiid].ident);
      continue;
    }
#endif

    uc->reportStartingWork();

#if 0
    //  Attempt at increasing quality for high error, didn't help.
    if (errorRate > 0.25) {
      if (uc->showAlgorithm())
        fprintf(stderr, "generateMultiAlignment()-- high error, decrease allowed error rate from %f to %f\n", errorRate, errorRate * 2 / 3);

      uc->setErrorRate(errorRate * 2 / 3);

      if (uc->computePositionFromAnchor()    && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromLayout()    && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromAlignment() && uc->alignFragment())  goto applyAlignment;

      uc->setErrorRate(errorRate);
    }
#endif

    //  First attempt, all default parameters

    if (uc->computePositionFromAnchor()    && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()    && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment() && uc->alignFragment())  goto applyAlignment;

    //  Second attempt, higher error rate.

    if (uc->showAlgorithm())
      fprintf(stderr, "generateMultiAlignment()-- increase allowed error rate from %f to %f\n", errorRate, MIN(errorRateMax, 2.0 * errorRate));

    uc->setErrorRate(MIN(errorRateMax, 2.0 * errorRate));

    if (uc->computePositionFromAnchor()    && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()    && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment() && uc->alignFragment())  goto applyAlignment;

    uc->setErrorRate(errorRate);

    //  Third attempt, thinner overlaps.  These come from bogart repeat splitting, it apanchorly
    //  doesn't enforce the minimum overlap length in those unitugs.

    while (minOverlap > 40) {
      if (uc->showAlgorithm())
        fprintf(stderr, "generateMultiAlignment()-- decrease minimum overlap from %d to %d\n", minOverlap, MAX(40, minOverlap / 2));

      uc->setMinOverlap(MAX(40, minOverlap / 2));

      if (uc->computePositionFromAnchor()    && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromLayout()    && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromAlignment() && uc->alignFragment())  goto applyAlignment;
    }

    uc->setMinOverlap(minOverlap);

    //  Fourth attempt, default parameters after recomputing consensus sequence.

    if (uc->showAlgorithm())
      fprintf(stderr, "generateMultiAlignment()-- recompute full consensus\n");

    uc->rebuild(true);

    if (uc->computePositionFromAnchor()    && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()    && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment() && uc->alignFragment())  goto applyAlignment;

    //  Final attempt, higher error rate.

    if (uc->showAlgorithm())
      fprintf(stderr, "generateMultiAlignment()-- increase allowed error rate from %f to %f\n", errorRate, MIN(errorRateMax, 4.0 * errorRate));

    uc->setErrorRate(MIN(errorRateMax, 4.0 * errorRate));

    if (uc->computePositionFromAnchor()    && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()    && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment() && uc->alignFragment())  goto applyAlignment;

    //  Failed to align the fragment.  Dang.

    uc->setErrorRate(errorRate);

#ifdef FAILURE_IS_FATAL
    fprintf(stderr, "FAILED TO ALIGN FRAG.  DIE.\n");
    assert(0);
#endif

    uc->reportFailure(failed);
    failuresToFix = true;
    continue;

  applyAlignment:
    uc->setErrorRate(errorRate);
    uc->setMinOverlap(minOverlap);

    uc->reportSuccess(failed);
    uc->applyAlignment();
    uc->rebuild(false);
  }

  if (failuresToFix)
    goto returnFailure;

  uc->generateConsensus();
  uc->exportToTig();

  //  NEED TO UPDATE THE tgTig

  delete uc;
  return(true);

 returnFailure:
  fprintf(stderr, "generateMultiAlignment()-- unitig %d FAILED.\n", tig->tigID());

  //  tgTig should have no changes.

  delete uc;
  return(false);
}
