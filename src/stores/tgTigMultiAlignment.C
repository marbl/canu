

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"
#include "aligners.H"
#include "unitigConsensus.H"



bool
generateMultiAlignment(tgTig     *tig,
                       gkStore   *gkpStore,
                       uint32    *failed) {
  double             origErate          = AS_CNS_ERROR_RATE;
  uint32             origLen            = AS_OVERLAP_MIN_LEN;
  bool               failuresToFix      = false;
  unitigConsensus   *uc                 = NULL;

  origErate          = AS_CNS_ERROR_RATE;
  uc                 = new unitigConsensus(tig);

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
    if (AS_CNS_ERROR_RATE > 0.25) {
      if (uc->showAlgorithm())
        fprintf(stderr, "generateMultiAlignment()-- high error, decrease allowed error rate from %f to %f\n", AS_CNS_ERROR_RATE, AS_CNS_ERROR_RATE * 2 / 3);

      AS_CNS_ERROR_RATE = AS_CNS_ERROR_RATE * 2 / 3;

      if (uc->computePositionFromAnchor()      && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

      AS_CNS_ERROR_RATE = origErate;
    }
#endif

    //  First attempt, all default parameters

    if (uc->computePositionFromAnchor()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    //  Second attempt, higher error rate.

    if (uc->showAlgorithm())
      fprintf(stderr, "generateMultiAlignment()-- increase allowed error rate from %f to %f\n", AS_CNS_ERROR_RATE, MIN(AS_MAX_ERROR_RATE, 2.0 * AS_CNS_ERROR_RATE));

    AS_CNS_ERROR_RATE = MIN(AS_MAX_ERROR_RATE, 2.0 * AS_CNS_ERROR_RATE);

    if (uc->computePositionFromAnchor()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    AS_CNS_ERROR_RATE = origErate;

    //  Third attempt, thinner overlaps.  These come from bogart repeat splitting, it apanchorly
    //  doesn't enforce the minimum overlap length in those unitugs.

    while (AS_OVERLAP_MIN_LEN > 40) {
      if (uc->showAlgorithm())
        fprintf(stderr, "generateMultiAlignment()-- decrease minimum overlap from %d to %d\n", AS_OVERLAP_MIN_LEN, MAX(40, AS_OVERLAP_MIN_LEN / 2));

      AS_OVERLAP_MIN_LEN = MAX(40, AS_OVERLAP_MIN_LEN / 2);

      if (uc->computePositionFromAnchor()      && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;
    }

    AS_OVERLAP_MIN_LEN = origLen;

    //  Fourth attempt, default parameters after recomputing consensus sequence.

    if (uc->showAlgorithm())
      fprintf(stderr, "generateMultiAlignment()-- recompute full consensus\n");

    uc->rebuild(true);

    if (uc->computePositionFromAnchor()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    //  Final attempt, higher error rate.

    if (uc->showAlgorithm())
      fprintf(stderr, "generateMultiAlignment()-- increase allowed error rate from %f to %f\n", AS_CNS_ERROR_RATE, MIN(AS_MAX_ERROR_RATE, 4.0 * AS_CNS_ERROR_RATE));

    AS_CNS_ERROR_RATE = MIN(AS_MAX_ERROR_RATE, 4.0 * AS_CNS_ERROR_RATE);

    if (uc->computePositionFromAnchor()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    //  Failed to align the fragment.  Dang.

    AS_CNS_ERROR_RATE = origErate;

#ifdef FAILURE_IS_FATAL
    fprintf(stderr, "FAILED TO ALIGN FRAG.  DIE.\n");
    assert(0);
#endif

    uc->reportFailure(failed);
    failuresToFix = true;
    continue;

  applyAlignment:
    AS_CNS_ERROR_RATE = origErate;

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
