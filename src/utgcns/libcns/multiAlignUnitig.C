




bool
multiAlignUnitig(MultiAlignT     *ma,
                 gkStore         *fragStore,
                 CNS_Options     *opp,
                 int32           *failed) {
  double             origErate          = AS_CNS_ERROR_RATE;
  uint32             origLen            = AS_OVERLAP_MIN_LEN;
  bool               failuresToFix      = false;
  unitigConsensus   *uc                 = NULL;

  origErate          = AS_CNS_ERROR_RATE;
  uc                 = new unitigConsensus(ma, opp);

  if (uc->initialize(failed) == FALSE)
    goto returnFailure;

  while (uc->moreFragments()) {


#ifdef SKIP_CONTAINS
    if (uc->fraglist[uc->tiid].contained != 0) {
      //  skip contained for now
      fprintf(stderr, "SKIP CONTAINED %u\n", uc->fraglist[uc->tiid].ident);
      continue;
    }
#endif

    if (VERBOSE_MULTIALIGN_OUTPUT)
      uc->reportStartingWork();

#if 0
    //  Attempt at increasing quality for high error, didn't help.
    if (AS_CNS_ERROR_RATE > 0.25) {
      if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
        fprintf(stderr, "MultiAlignUnitig()-- high error, decrease allowed error rate from %f to %f\n", AS_CNS_ERROR_RATE, AS_CNS_ERROR_RATE * 2 / 3);

      AS_CNS_ERROR_RATE = AS_CNS_ERROR_RATE * 2 / 3;

      if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

      AS_CNS_ERROR_RATE = origErate;
    }
#endif

    //  First attempt, all default parameters

    if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    //  Second attempt, higher error rate.

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "MultiAlignUnitig()-- increase allowed error rate from %f to %f\n", AS_CNS_ERROR_RATE, MIN(AS_MAX_ERROR_RATE, 2.0 * AS_CNS_ERROR_RATE));

    AS_CNS_ERROR_RATE = MIN(AS_MAX_ERROR_RATE, 2.0 * AS_CNS_ERROR_RATE);

    if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    AS_CNS_ERROR_RATE = origErate;

    //  Third attempt, thinner overlaps.  These come from bogart repeat splitting, it apparently
    //  doesn't enforce the minimum overlap length in those unitugs.

    while (AS_OVERLAP_MIN_LEN > 40) {
      if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
        fprintf(stderr, "MultiAlignUnitig()-- decrease minimum overlap from %d to %d\n", AS_OVERLAP_MIN_LEN, MAX(40, AS_OVERLAP_MIN_LEN / 2));

      AS_OVERLAP_MIN_LEN = MAX(40, AS_OVERLAP_MIN_LEN / 2);

      if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
      if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;
    }

    AS_OVERLAP_MIN_LEN = origLen;

    //  Fourth attempt, default parameters after recomputing consensus sequence.

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "MultiAlignUnitig()-- recompute full consensus\n");

    uc->rebuild(true);

    if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromLayout()      && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromAlignment()   && uc->alignFragment())  goto applyAlignment;

    //  Final attempt, higher error rate.

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ALGORITHM)
      fprintf(stderr, "MultiAlignUnitig()-- increase allowed error rate from %f to %f\n", AS_CNS_ERROR_RATE, MIN(AS_MAX_ERROR_RATE, 4.0 * AS_CNS_ERROR_RATE));

    AS_CNS_ERROR_RATE = MIN(AS_MAX_ERROR_RATE, 4.0 * AS_CNS_ERROR_RATE);

    if (uc->computePositionFromParent(false) && uc->alignFragment())  goto applyAlignment;
    if (uc->computePositionFromParent(true)  && uc->alignFragment())  goto applyAlignment;
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

  delete uc;
  return(true);

 returnFailure:
  fprintf(stderr, "MultiAlignUnitig()-- unitig %d FAILED.\n", ma->maID);

  uc->restoreUnitig();

  delete uc;
  return(false);
}
