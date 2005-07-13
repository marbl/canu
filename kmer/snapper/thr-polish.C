#include "thr.H"



char*
doPolish(searcherState       *state,
         FastASequenceInCore *seq,
         aHit               *&theHits,
         u32bit              &theHitsLen,
         logMsg              *theLog) {
  double   startTime = getTime();
  u32bit   outputLen = 0;
  u32bit   outputMax = 2 * 1024 * theHitsLen;
  char    *output    = 0L;

  if (theHitsLen == 0) {
    output    = new char [8];
    output[0] = 0;
    return(output);
  }


  try {
    output = new char [outputMax];
  } catch (...) {
    fprintf(stderr, "doPolish()-- Can't allocate space for the output string ("u32bitFMT" bytes) in thread "u64bitFMT"\n", outputMax, state->threadID);
    abort();
  }

  output[0] = 0;

  for (u32bit h=0; h<theHitsLen; h++) {
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

    } else {

      if ((config._doValidation) ||
          (theHits[h]._status & AHIT_POLISHABLE)) {
        FastASequenceInCore  *ESTseq = seq;
        FastASequenceInCore  *GENseq = cache->getSequence(theHits[h]._dsIdx);
        u32bit                GENlo  = theHits[h]._dsLo;
        u32bit                GENhi  = theHits[h]._dsHi;

        if (GENhi > GENseq->sequenceLength())
          GENhi = GENseq->sequenceLength();

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
              theLog->add(stderr, "FORWARDHIT GEN: hi:"u32bitFMT"-lo:"u32bitFMT" pos:"u32bitFMT" EST: len:"u32bitFMT" pos:"u32bitFMT"\n",
                          GENhi, GENlo, y, (u32bit)ESTseq->sequenceLength(), x);
#endif
            P4->addSeed(y - GENlo + config._merSize,
                        x         + config._merSize,
                        config._merSize);
          }
        } else {
          for (u32bit i=0, x, y; theHits[h]._ML->getMer(i, x, y); i++) {
#ifdef SHOW_HITS_ADDED
#ifdef SHOW_HITS_ADDED_AFTER_QUERY
            if (ESTseq->getIID() > SHOW_HITS_ADDED_AFTER_QUERY)
#endif
              theLog->add(stderr, "REVERSEHIT GEN: hi:"u32bitFMT"-lo:"u32bitFMT" pos:"u32bitFMT" EST: len:"u32bitFMT" pos:"u32bitFMT"\n",
                          GENhi, GENlo, y, (u32bit)ESTseq->sequenceLength(), x);
#endif
            //  Original form was (GENhi-GENlo) - (y-GENlo), which
            //  reduces to the below.  By reversing, we no longer need
            //  to add in the mersize, we're representing the end of
            //  the mer now!
            //
            P4->addSeed(GENhi                    - y,
                        ESTseq->sequenceLength() - x,
                        config._merSize);
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
          if ((L4[i]->percentIdentity >= pi) && (L4[i]->querySeqIdentity >= pc)) {
            pi = L4[i]->percentIdentity;
            pc = L4[i]->querySeqIdentity;
          }

#ifdef SHOW_POLISHING
          theLog->add("  match["u32bitFMT"] : id=%u cov=%u\n", i, L4[i]->percentIdentity, L4[i]->querySeqIdentity);
#endif

          //  If we have a real hit, set the flag and save the output
          //
          if ((L4[i]->percentIdentity  >= config._minMatchIdentity) &&
              (L4[i]->querySeqIdentity >= config._minMatchCoverage)) {

            theHits[h]._status |= AHIT_VERIFIED;

            char *pstr = s4p_polishToString(L4[i]);

            u32bit l = (u32bit)strlen(pstr);
            if (outputLen + l + 1 >= outputMax) {
              outputMax = outputMax + outputMax + l;
              char *o = 0L;
              try {
                o = new char [outputMax];
              } catch (...) {
                fprintf(stderr, "doPolish()-- Can't reallocate space for the output string ("u32bitFMT" bytes) in thread "u64bitFMT"\n", outputMax, state->threadID);
                abort();
              }
              memcpy(o, output, sizeof(char) * outputLen);
              delete [] output;
              output = o;
            }

            memcpy(output + outputLen, pstr, sizeof(char) * l);
            outputLen += l;
            output[outputLen] = 0;

            free(pstr);
          }
        }

        //  Save the best scores
        //
        theHits[h]._status |= pi << 16;
        theHits[h]._status |= pc << 24;

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
      }
    }
  }

  state->polishTime += getTime() - startTime;

  return(output);
}




