#include "snapper2.H"



char*
doPolishS4(searcherState       *state,
           seqInCore           *seq,
           aHit               *&theHits,
           u32bit              &theHitsLen,
           logMsg              *theLog) {
  double   startTime = getTime();
  u32bit   outLen    = 0;
  u32bit   outMax    = 2 * 1024 * theHitsLen;
  char    *out       = 0L;

  //  For the autofilter
  u64bit   successes    = u64bitZERO;
  u64bit   successMask  = u64bitMASK(config._afLength);
  u32bit   attempts     = 0;

  if (theHitsLen == 0) {
    out    = new char [8];
    out[0] = 0;
    return(out);
  }

  try {
    out = new char [outMax];
  } catch (...) {
    fprintf(stderr, "doPolish()-- Can't allocate space for the output string ("u32bitFMT" bytes) in thread "u64bitFMT"\n", outMax, state->threadID);
    abort();
  }

  out[0] = 0;

  for (u32bit h=0; h<theHitsLen; h++) {

    //  If the hit was discarded, move along.
    //
    if (theHits[h]._status & AHIT_DISCARDED) {
#ifdef SHOW_HIT_DISCARDING
      theLog->add("Hit %u out of %u (%u -> %u[%u-%u]) cov=%u matched=%u numMers=%u DISCARDED\n",
                  h, theHitsLen,
                  seq->getIID(),
                  theHits[h]._dsIdx,
                  theHits[h]._dsLo,
                  theHits[h]._dsHi,
                  theHits[h]._covered,
                  theHits[h]._matched,
                  theHits[h]._numMers);
#endif
      continue;
    }


    //  If the hit was filtered out, move along.
    //
    if ((config._doValidation == false) &&
        ((theHits[h]._status & AHIT_POLISHABLE) == 0) &&
        ((theHits[h]._status & AHIT_HAS_UNIQUE) == 0))
      continue;


    //  If our recent success rate is pretty terrible, continue.
    //
    if (config._afEnabled) {

      if (attempts > config._afInit) {
        double  rat = countNumberOfSetBits64(successes) / (double)((attempts < config._afLength) ? attempts : config._afLength);

#if 0
        fprintf(stderr, "autofilter: hit "u32bitFMT" out of "u32bitFMT" (attempts="u32bitFMT") with rate %f\n",
                h, theHitsLen, attempts, rat);
#endif

        //  If we've hit the end of the good polishes, give up.  But
        //  still do all the stuff with unique mers in them.
        //
        if (((theHits[h]._status & AHIT_HAS_UNIQUE) == 0) &&
            (rat < config._afThreshold))
          continue;
      }

      attempts++;
    }



    //
    //  Polish it up!
    //


    seqInCore            *ESTseq = seq;
    seqInCore            *GENseq = genome->getSequenceInCore(theHits[h]._dsIdx);
    u32bit                GENlo  = theHits[h]._dsLo;
    u32bit                GENhi  = theHits[h]._dsHi;

    if (GENhi > GENseq->sequenceLength())
      GENhi = GENseq->sequenceLength();

    assert(GENlo < GENhi);

    bool    doForward =  theHits[h]._status & AHIT_DIRECTION_MASK;
    bool    doReverse = !doForward;

#ifdef SHOW_POLISHING
    theLog->add("Hit %u out of %u (%u -> %u[%u-%u]) dir=%c cov=%u matched=%u numMers=%u\n",
                h, theHitsLen,
                ESTseq->getIID(),
                theHits[h]._dsIdx,
                theHits[h]._dsLo,
                theHits[h]._dsHi,
                doForward ? 'F' : 'R',
                theHits[h]._covered,
                theHits[h]._matched,
                theHits[h]._numMers);
#endif


#ifdef SHOW_POLISHING_EXPENSIVE
    double startTime = getTime();
#endif

    sim4command     *P4 = new sim4command(ESTseq,
                                          GENseq,
                                          GENlo,
                                          GENhi,
                                          doForward,
                                          doReverse);



    ////////////////////////////////////////
    //
    //  Add hits to the command
    //
    //  addSeed() expects base-based, of the last position in
    //  the seed.  We have space-based, first position.  Adding
    //  the size of a mer fixes both.
    //
    if (doForward) {
      for (u32bit i=0, x, y; theHits[h]._ML->getMer(i, x, y); i++) {
#ifdef SHOW_HITS_ADDED
#ifdef SHOW_HITS_ADDED_AFTER_QUERY
        if (ESTseq->getIID() > SHOW_HITS_ADDED_AFTER_QUERY)
#endif
          theLog->add("FORWARDHIT GEN: hi:"u32bitFMT"-lo:"u32bitFMT" pos:"u32bitFMT" EST: len:"u32bitFMT" pos:"u32bitFMT"\n",
                      GENhi, GENlo, y, (u32bit)ESTseq->sequenceLength(), x);
#endif
        assert(y + config._KBmerSize >= GENlo);

        P4->addSeed(y - GENlo + config._KBmerSize,
                    x         + config._KBmerSize,
                    config._KBmerSize);
      }
    } else {
      for (u32bit i=0, x, y; theHits[h]._ML->getMer(i, x, y); i++) {
#ifdef SHOW_HITS_ADDED
#ifdef SHOW_HITS_ADDED_AFTER_QUERY
        if (ESTseq->getIID() > SHOW_HITS_ADDED_AFTER_QUERY)
#endif
          theLog->add("REVERSEHIT GEN: hi:"u32bitFMT"-lo:"u32bitFMT" pos:"u32bitFMT" EST: len:"u32bitFMT" pos:"u32bitFMT"\n",
                      GENhi, GENlo, y, (u32bit)ESTseq->sequenceLength(), x);
#endif
        //  Original form was (GENhi-GENlo) - (y-GENlo), which
        //  reduces to the below.  By reversing, we no longer need
        //  to add in the mersize, we're representing the end of
        //  the mer now!
        //
        assert(GENhi                    >= y);
        assert(ESTseq->sequenceLength() >= x);

        P4->addSeed(GENhi                    - y,
                    ESTseq->sequenceLength() - x,
                    config._KBmerSize);
      }
    }



    //  The main loop deletes the hits, but we take care of deleting _ML here.
    //  Maybe it should go in the destructor for the hits??
    //
    delete theHits[h]._ML;
    theHits[h]._ML = 0L;


    Sim4            *S4 = new Sim4(&sim4params);
    sim4polishList  *l4 = S4->run(P4);
    sim4polishList  &L4 = *l4;


    //  Clean up the matches -- remove small exons from the match,
    //  split things with big gaps into two matches.

    for (u32bit i=0; L4[i]; i++) {

#ifdef SHOW_MATCH_SPLITTING
      theLog->add("  match "u32bitFMT" has "u32bitFMT" exons.\n",
                  i, L4[i]->numExons);
      for (u32bit j=L4[i]->numExons; j--; )
        theLog->add("    exon "u32bitFMT" query:"u32bitFMT"-"u32bitFMT" genome:"u32bitFMT"-"u32bitFMT" id:%d nm:%d\n",
                    j,
                    L4[i]->exons[j].estFrom,
                    L4[i]->exons[j].estTo,
                    L4[i]->exons[j].genFrom,
                    L4[i]->exons[j].genTo,
                    L4[i]->exons[j].percentIdentity,
                    L4[i]->exons[j].numMatches);

#endif

      for (u32bit j=L4[i]->numExons; j--; ) {
        if (((L4[i]->exons[j].estTo - L4[i]->exons[j].estFrom) < config._discardExonLength)  ||
            (L4[i]->exons[j].percentIdentity < config._discardExonQuality)) {
#ifdef SHOW_MATCH_SPLITTING
          theLog->add("    Deleting exon "u32bitFMT" from query:"u32bitFMT"-"u32bitFMT" genome:"u32bitFMT"-"u32bitFMT"\n",
                      j,
                      L4[i]->exons[j].estFrom,
                      L4[i]->exons[j].estTo,
                      L4[i]->exons[j].genFrom,
                      L4[i]->exons[j].genTo);
#endif
          s4p_deleteExon(L4[i], j);
        }
      }

      //  Copy each exon into a new match ("split things with big gaps")

      while (L4[i]->numExons > 1) {
#ifdef SHOW_MATCH_SPLITTING
        theLog->add("    Saving exon "u32bitFMT" from query:"u32bitFMT"-"u32bitFMT" genome:"u32bitFMT"-"u32bitFMT"\n",
                    L4[i]->numExons-1,
                    L4[i]->exons[L4[i]->numExons-1].estFrom,
                    L4[i]->exons[L4[i]->numExons-1].estTo,
                    L4[i]->exons[L4[i]->numExons-1].genFrom,
                    L4[i]->exons[L4[i]->numExons-1].genTo);
#endif

        sim4polish *n = s4p_copyPolish_OneExon(L4[i], L4[i]->numExons-1);
        L4.push(n);
        s4p_deleteExon(L4[i], L4[i]->numExons-1);
      }

      //  Rebuild the stats on this guy -- we now have one exon, so just copy
      //  the exon stats to the global stats.

      if (L4[i]->numExons > 0) {
#ifdef SHOW_MATCH_SPLITTING
        theLog->add("    Saving exon "u32bitFMT" from query:"u32bitFMT"-"u32bitFMT" genome:"u32bitFMT"-"u32bitFMT"\n",
                    0,
                    L4[i]->exons[0].estFrom,
                    L4[i]->exons[0].estTo,
                    L4[i]->exons[0].genFrom,
                    L4[i]->exons[0].genTo);
#endif

        L4[i]->numMatches       = L4[i]->exons[0].numMatches;
        L4[i]->numMatchesN      = L4[i]->exons[0].numMatchesN;
        L4[i]->numCovered       = L4[i]->exons[0].genTo - L4[i]->exons[0].genFrom + 1;
        L4[i]->percentIdentity  = L4[i]->exons[0].percentIdentity;
        L4[i]->querySeqIdentity = s4p_percentCoverageApprox(L4[i]);
      } else {
#ifdef SHOW_MATCH_SPLITTING
        theLog->add("    All exons removed!\n");
#endif
        L4.remove(i);
        i--;
      }
    }


    //  Even though we don't expect multiple polishes, we still have to deal with
    //  them.  :-(

    //  Clear the 'match' flag and set qualities to zero.  XXX:
    //  Again, this should be already done, but we need to guarantee
    //  it.
    //
    //theHits[h]._status &= 0x00000003;
    //  (I guess we don't _need_ to do it....)

    u32bit  pi = 0;
    u32bit  pc = 0;

    for (u32bit i=0; L4[i]; i++) {

      //  We need to remember the best pair of percent
      //  identity/coverage.  These wil be stored in the hit after
      //  we process all matches.
      //
      if ((L4[i]->percentIdentity >= pi) &&
          (L4[i]->querySeqIdentity >= pc)) {
        pi = L4[i]->percentIdentity;
        pc = L4[i]->querySeqIdentity;
      }

#ifdef SHOW_POLISHING
      theLog->add("  match["u32bitFMT"] query:"u32bitFMT"-"u32bitFMT" genome:"u32bitFMT"-"u32bitFMT" id=%u cv=%d nm=%u\n",
                  i,
                  L4[i]->exons[0].estFrom,
                  L4[i]->exons[0].estTo,
                  L4[i]->exons[0].genFrom,
                  L4[i]->exons[0].genTo,
                  L4[i]->percentIdentity,
                  L4[i]->querySeqIdentity,
                  L4[i]->exons[0].numMatches);
#endif

      //  If we have a real hit, set the flag and save the output
      //
      if ((L4[i]->percentIdentity  >= config._minMatchIdentity) &&
          (L4[i]->querySeqIdentity >= config._minMatchCoverage)) {

        theHits[h]._status |= AHIT_VERIFIED;

        char *pstr = s4p_polishToString(L4[i]);

        u32bit l = (u32bit)strlen(pstr);
        if (outLen + l + 1 >= outMax) {
          outMax = outMax + outMax + l;
          char *o = 0L;
          try {
            o = new char [outMax];
          } catch (...) {
            fprintf(stderr, "doPolish()-- Can't reallocate space for the output string ("u32bitFMT" bytes) in thread "u64bitFMT"\n", outMax, state->threadID);
            abort();
          }
          memcpy(o, out, sizeof(char) * outLen);
          delete [] out;
          out = o;
        }

        memcpy(out + outLen, pstr, sizeof(char) * l);
        outLen += l;
        out[outLen] = 0;

        free(pstr);
      }
    }

    //  Save the best scores
    //
    theHits[h]._status |= pi << 16;
    theHits[h]._status |= pc << 24;

    successes <<= 1;
    if ((pi  >= config._minMatchIdentity) &&
        (pc >= config._minMatchCoverage)) {
      //fprintf(stderr, "GOOD "u32bitFMT" "u32bitFMT"\n", pi, pc);
      successes |= u64bitONE;
    } else {
      //fprintf(stderr, "BAD  "u32bitFMT" "u32bitFMT"\n", pi, pc);
      successes |= u64bitZERO;
    }
    successes  &= successMask;

    delete l4;
    delete S4;
    delete P4;

#ifdef SHOW_POLISHING_EXPENSIVE
    double elapsedTime = getTime() - startTime;
    if (elapsedTime >= SHOW_POLISHING_EXPENSIVE) {
      theLog->add("Hit %u out of %u (%u -> %u[%u-%u]) took %f seconds ().\n",
                  h, theHitsLen,
                  ESTseq->getIID(), GENseq->getIID(), theHits[h]._dsLo, theHits[h]._dsHi,
                  elapsedTime);
    }
#endif

    state->polished++;
  }  //  over all hits

  state->polishTime += getTime() - startTime;
  return(out);
}




