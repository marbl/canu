#include "sim4.H"
#include "sim4polishBuilder.H"

//#define SHOW_OVERLAPPING_EXONS


static void
add_offset_exons(Exon *exons, int offset) {
  if (!offset || !exons)
    return;

  for (; exons; exons = exons->next_exon) {
    if (exons->toGEN) {
      exons->frEST += offset;
      exons->toEST += offset;
    }
  }
}


#if 0
static void
add_offset_aligns(edit_script_list *aligns, int offset) {
  if (!offset || !aligns)
    return;

  for (; aligns; aligns = aligns->next_script)
    aligns->offset2 += offset;
}
#endif


void
Sim4::maskExonsFromSeeds(sim4command *cmd,
                         Exon *theExon) {

  while (theExon) {
    if (theExon->toGEN) {
      for (u32bit x=0; x<cmd->numberOfExternalSeeds(); x++) {
        u32bit pos = cmd->externalSeedGENPosition(x);

        if (((u32bit)theExon->frGEN <= pos + 1) &&
            (pos <= (u32bit)theExon->toGEN + cmd->externalSeedLength(x)))
          cmd->maskExternalSeed(x);
      }
    }
    theExon = theExon->next_exon;
  }
}

void
Sim4::maskExonsFromGenomic(Exon *theExon,
                           char *f,
                           char *r,
                           int l) {

  while (theExon) {
    if (theExon->toGEN) {
      for (int i=theExon->frGEN-1; i<theExon->toGEN; i++)
        f[i] = 'N';
      for (int i=l-theExon->frGEN; i>=l-theExon->toGEN; i--)
        r[i] = 'N';
    }
    theExon = theExon->next_exon;
  }
}




sim4polishList*
Sim4::run(sim4command *cmd) {
  sim4polishBuilder   B;
  sim4polishList     *L = new sim4polishList;


  int    dist, match_ori;
  int    g_pA=0, f_pA=0, r_pA=0;
  int    g_pT=0, f_pT=0, r_pT=0;

  Exon   *fExons = NULL;
  Exon   *rExons = NULL; 

  edit_script_list *fAligns = NULL;
  edit_script_list *rAligns = NULL;

  int     matchesPrinted = 0;

  char    touppercache[256];

  for (int i=0; i<256; i++)
    touppercache[i] = (char)toupper(i);

  cmd->finalize();

  u32bit  dblen       = cmd->getGENhi() - cmd->getGENlo();
  char   *dbseq       = 0L;
  char   *dbrev       = 0L;
  char   *dbseqorig   = cmd->getGENsequence();

  int     estlen     = 0;
  char   *estseq     = 0L;
  char   *estrev     = 0L;
  char   *estseqorig = 0L;


  //  Allocate space for temporary sequence storage.  We need
  //  to allocate space for two copies of the database, and space
  //  for the longest EST (in case we need to print it out
  //  reverse complemented).
  //
  char   *seqStorage     = 0L;
  u32bit  seqStorageSize = 0;

  seqStorageSize  = 2 * dblen + 2 * cmd->getESTlength() + 8;
  seqStorage      = new char [seqStorageSize];

  //  Original, forward, reverse, cdna
  //
  dbseq  = seqStorage;
  dbrev  = seqStorage + dblen + 2;
  estseq = seqStorage + dblen + 2 + dblen + 2;
  estrev = seqStorage + dblen + 2 + dblen + 2 + cmd->getESTlength() + 2;


  //  Prepare the database sequence
  //
  //  Trimming to the correct range
  //  Convert to uppercase
  //  Reverse complement
  //
  for (u32bit i=0, j=cmd->getGENlo(), k=dblen-1; j<cmd->getGENhi(); i++, j++, k--) {
    dbseq[i] = touppercache[(int)dbseqorig[j]];
    dbrev[k] = complementSymbol[(int)dbseq[i]];
  }
  dbseq[dblen] = 0;
  dbrev[dblen] = 0;

  sim4_stats_t   st, rev_st;

  estseqorig = cmd->getESTsequence();
  estlen     = cmd->getESTlength();

  for (int i=0; i<estlen; i++)
    estseq[i] = touppercache[(int)estseqorig[i]];
  estseq[estlen] = 0;

  g_pT = g_pA = 0;

  if (globalParams->_ignorePolyTails) {
    get_polyAT(estseq, estlen, &g_pT, &g_pA);
  }


  //  GRRR!  XXXXX  This needs to be defined outside the loop, and before the goto's
  bool  pleaseContinueComputing = false;


  if (estlen - g_pA - g_pT <= 0)
    goto abort;


  matchesPrinted = 0;

  do {
    //fprintf(stderr, "sim4string::main loop begins!\n");

    int     nmatches  = 0;
    double  coverage  = 0;
    int     percentid = 0;

    pleaseContinueComputing = false;

    B.create(cmd->getESTidx(), estlen,
             cmd->getGENidx(), cmd->getGENlo(), cmd->getGENhi());

    if (globalParams->_includeDefLine) {
      B.setESTdefline(cmd->getESTheader());
      B.setGENdefline(cmd->getGENheader());
    }

    memset(&st,     0, sizeof(sim4_stats_t));
    memset(&rev_st, 0, sizeof(sim4_stats_t));

    if (cmd->externalSeedsExist() == false)
      bld_table(estseq - 1 + g_pT, estlen - g_pA - g_pT, wordSize, INIT);

    if (cmd->doForward()) {

      //  Initialize the sequences and lengths
      //
      //  genSeq was seq1
      //  estSeq was seq2
      //
      _genSeq = dbseq;
      _estSeq = estseq + g_pT;
      _genLen = dblen;
      _estLen = estlen - g_pT - g_pA;

      //  This should be in a better spot.
      _mspManager.setLength(_estLen);
      _mspManager.clearDiagonal(_genLen, _estLen);
      _mspManager.setScoreThreshold(mspThreshold1, globalParams->_interspecies);

#ifdef SHOW_EXTERNAL_SEEDING
      fprintf(stderr, "FWD: estLen = %d  genLen = %d\n", _estLen, _genLen);
#endif

      //  Find the seeds.
      //
      if (cmd->externalSeedsExist() == false) {
        exon_cores(_genSeq-1, _estSeq-1, _genLen, _estLen, 1, 1, 0, wordSize, mspThreshold1, PERM);
      } else {
#ifdef SHOW_EXTERNAL_SEEDING
        fprintf(stderr, "FWD: Using external seeds -- adding "u32bitFMT" seeds to sim4.\n", cmd->numberOfExternalSeeds());
#endif

        cmd->sortExternalSeeds();

        for (u32bit x=0; x<cmd->numberOfExternalSeeds(); x++)
          if (cmd->externalSeedLength(x) > 0)
            _mspManager.addHit(_genSeq-1, _estSeq-1,
                               _genLen, _estLen,
                               cmd->externalSeedGENPosition(x),
                               cmd->externalSeedESTPosition(x),
                               cmd->externalSeedLength(x));

        exon_list = _mspManager.doLinking(DEFAULT_WEIGHT, DEFAULT_DRANGE,
                                          1, 1,
                                          0,
                                          false,
                                          _genSeq, _estSeq);

#ifdef SHOW_EXTERNAL_SEEDING
        fprintf(stderr, "FWD: Added and chained, starting SIM4() run.\n");
#endif
      }

      fAligns = SIM4(&dist,
                     &fExons,
                     &f_pA,
                     &f_pT,
                     &st);

      //  Continued from util.C :: slide_intron()
      //
      //  If we are forcing the strand prediction, and we are still unknown,
      //  set the strand prediction to the match orientation. Since this 
      //  will be reversed later on, set it to FWD here.
      //
      if ((globalParams->_forceStrandPrediction) && (st.orientation == BOTH))
        st.orientation = FWD;

      //  If the match was deemed expensive, report
      //
      if (st.tooManyMSPs) {
        B.setNumberOfMatches(0, 0);
        B.setPercentIdentity(0);
        B.setMatchOrientation(SIM4_MATCH_FORWARD);
        B.setStrandOrientation(SIM4_STRAND_INTRACTABLE);
        B.addExon(1, estlen,
                  1, cmd->getGENhi() - cmd->getGENlo(),
                  st.numberOfMatches, 0, 0,
                  SIM4_INTRON_NONE);
        goto fail;
      }
    }

    if (cmd->doReverse()) {
      //  Initialize the sequences and lengths
      //
      //  genSeq was seq1
      //  estSeq was seq2
      //
      _genSeq = dbrev;
      _estSeq = estseq + g_pT;
      _genLen = dblen;
      _estLen = estlen - g_pT - g_pA;

      //  This should be in a better spot.
      _mspManager.setLength(_estLen);
      _mspManager.clearDiagonal(_genLen, _estLen);
      _mspManager.setScoreThreshold(mspThreshold1, globalParams->_interspecies);

#ifdef SHOW_EXTERNAL_SEEDING
      fprintf(stderr, "BWD: estLen = %d  genLen = %d g_pT=%d g_pA=%d\n", _estLen, _genLen, g_pT, g_pA);
#endif

      //  Find the seeds.
      //
      if (cmd->externalSeedsExist() == false) {
        exon_cores(_genSeq-1, _estSeq-1, _genLen, _estLen, 1, 1, 0, wordSize, mspThreshold1, PERM);
      } else {
#ifdef SHOW_EXTERNAL_SEEDING
        fprintf(stderr, "BWD: Using external seeds -- adding "u32bitFMT" seeds to sim4.\n", cmd->numberOfExternalSeeds());
#endif

        cmd->sortExternalSeeds();

        //  We have sorted the seeds in incresing genomic position,
        //  but we need to reverse everything.  We can do this by just
        //  adding the seeds backwards!
        //
        for (u32bit x=cmd->numberOfExternalSeeds(); x--; )
          if (cmd->externalSeedLength(x) > 0)
            _mspManager.addHit(_genSeq-1, _estSeq-1,
                               _genLen, _estLen,
                               cmd->externalSeedGENPosition(x),
                               cmd->externalSeedESTPosition(x),
                               cmd->externalSeedLength(x));

        exon_list = _mspManager.doLinking(DEFAULT_WEIGHT, DEFAULT_DRANGE,
                                          1, 1,
                                          0,
                                          false,
                                          _genSeq, _estSeq);
#ifdef SHOW_EXTERNAL_SEEDING
        fprintf(stderr, "BWD: Added and chained, starting SIM4() run.\n");
#endif
      }

      rAligns = SIM4(&dist,
                     &rExons,
                     &r_pA,
                     &r_pT,
                     &rev_st);

      //  Continued from util.C :: slide_intron()
      //
      //  If we are forcing the strand prediction, and we are still unknown,
      //  set the strand prediction to the match orientation. 
      //
      if ((globalParams->_forceStrandPrediction) && (rev_st.orientation == BOTH))
        rev_st.orientation = FWD;

      //  If the match was deemed expensive, report
      if (rev_st.tooManyMSPs) { 
        B.setNumberOfMatches(0, 0);
        B.setPercentIdentity(0);
        B.setMatchOrientation(SIM4_MATCH_COMPLEMENT);
        B.setStrandOrientation(SIM4_STRAND_INTRACTABLE);
        B.addExon(1, estlen,
                  1, cmd->getGENhi() - cmd->getGENlo(),
                  rev_st.numberOfMatches, 0, 0,
                  SIM4_INTRON_NONE);
        goto fail;
      }
    }


    if (st.numberOfMatches >= rev_st.numberOfMatches) {
      match_ori = FWD;

      if (globalParams->_ignorePolyTails) {
        add_offset_exons(fExons, g_pT);  
        
        //add_offset_aligns(fAligns, g_pT);
        for (edit_script_list *aligns = fAligns; aligns; aligns = aligns->next_script)
          aligns->offset2 += g_pT;
      }

      B.setPolyTails(g_pA + f_pA, g_pT + f_pT);

      if (fExons) {
        //  We used to mask the seeds down with the masking of the
        //  genomic, but reverse exons are flipped here, and we need
        //  unflipped exons to mask.
        //
        if (cmd->externalSeedsExist() && globalParams->_findAllExons)
          maskExonsFromSeeds(cmd, fExons);

        if (checkExonsForOverlaps(fExons)) {
#ifdef SHOW_OVERLAPPING_EXONS
          B.setNumberOfMatches(0, 0);
          B.setPercentIdentity(0);
          B.setMatchOrientation(SIM4_MATCH_FORWARD);
          B.setStrandOrientation(SIM4_STRAND_FAILED);

          //  XXX:  result contains the exons and alignments
          //B.addExon(1, estlen, 1, cmd->getGENhi() - cmd->getGENlo(), rev_st.numberOfMatches, 0, SIM4_INTRON_NONE);
#endif
          goto fail;
        }
      }
    } else {
      match_ori = BWD;

      if (globalParams->_ignorePolyTails) {
        add_offset_exons(rExons, g_pT); 

        //add_offset_aligns(rAligns, g_pT);
        for (edit_script_list *aligns = rAligns; aligns; aligns = aligns->next_script)
          aligns->offset2 += g_pT;
      }

      B.setPolyTails(g_pA + r_pA, g_pT + r_pT);

      if (rAligns && rAligns->next_script)
        script_flip_list(&rAligns);

      if (rExons) {
        if (cmd->externalSeedsExist() && globalParams->_findAllExons)
          maskExonsFromSeeds(cmd, rExons);

        //  This used to be right before appendExons() in
        //  the reverse match section, but we need it
        //  before we test for overlapping exons
        //
        complement_exons(&rExons, dblen, estlen);

        if (checkExonsForOverlaps(rExons)) {
#ifdef SHOW_OVERLAPPING_EXONS
          B.setNumberOfMatches(0, 0);
          B.setPercentIdentity(0);
          B.setMatchOrientation(SIM4_MATCH_COMPLEMENT);
          B.setStrandOrientation(SIM4_STRAND_FAILED);

          //  XXX:  result contains the exons and alignments
          //B.addExon(1, estlen, 1, cmd->getGENhi() - cmd->getGENlo(), rev_st.numberOfMatches, 0, SIM4_INTRON_NONE);
#endif
          goto fail;
        }
      }
    }

    if (match_ori == FWD) {
      nmatches  = st.numberOfMatches;
      percentid = st.percentID;
    } else {
      nmatches  = rev_st.numberOfMatches;
      percentid = rev_st.percentID;
    }

    coverage  = (double)nmatches / (double)estlen;


    //  Is this match decent?
    //
    pleaseContinueComputing = ((coverage  >= globalParams->_minCoverage) &&
                               (percentid >= globalParams->_minPercentExonIdentity) &&
                               (nmatches  >= globalParams->_minCoverageLength) &&
                               (nmatches  > 0));

    //  If we're supposed to print at least _alwaysReport things,
    //  and we found a match, keep going.
    //
    if ((matchesPrinted < globalParams->_alwaysReport) && (nmatches > 0))
      pleaseContinueComputing = true;

    //  However, if we have printed enough stuff, and the last one is
    //  below the thresholds, stop.
    //
    if ((matchesPrinted >= globalParams->_alwaysReport) &&
        ((coverage  < globalParams->_minCoverage) ||
         (percentid < globalParams->_minPercentExonIdentity)))
      pleaseContinueComputing = false;


    if (pleaseContinueComputing) {
      matchesPrinted++;

      if (match_ori == FWD) {
        B.setNumberOfMatches(st.numberOfMatches, st.numberOfNs);
        B.setPercentIdentity(st.percentID);
        B.setMatchOrientation(SIM4_MATCH_FORWARD);

        switch (st.orientation) {
          case FWD:
            B.setStrandOrientation(SIM4_STRAND_POSITIVE);
            break;
          case BWD:
            B.setStrandOrientation(SIM4_STRAND_NEGATIVE);
            break;
          default:
            B.setStrandOrientation(SIM4_STRAND_UNKNOWN);
            break;
        }
      } else {
        B.setNumberOfMatches(rev_st.numberOfMatches, rev_st.numberOfNs);
        B.setPercentIdentity(rev_st.percentID);
        B.setMatchOrientation(SIM4_MATCH_COMPLEMENT);
        B.setStrandOrientation(SIM4_STRAND_FAILED);

        switch (rev_st.orientation) {
          case FWD:
            B.setStrandOrientation(SIM4_STRAND_NEGATIVE);
            break;
          case BWD:
            B.setStrandOrientation(SIM4_STRAND_POSITIVE);
            break;
          default:
            B.setStrandOrientation(SIM4_STRAND_UNKNOWN);
            break;
        }
      }


      //  If we have external seeds, we need to mask out seeds that we
      //  used BEFORE we print alignments -- printing reverse
      //  alignments also switches from reverse-complemented genomic
      //  to reverse-complemented EST, and then we can't (easily) mask
      //  seeds!
      //
      //  Likewise, we can't do the normal masking before we print the
      //  alignments, else we'd just print out N's for the genome.
      //
      if (match_ori == FWD) {
        appendExons(B, fExons);

        if (globalParams->_printAlignments) {
          appendAlignments(B,
                           estseq, dbseq, estlen, dblen,
                           fAligns, fExons,
                           FWD);
        }

        if (globalParams->_findAllExons)
          maskExonsFromGenomic(fExons, dbseq, dbrev, dblen);
      } else {
        appendExons(B, rExons);

        if (globalParams->_printAlignments) {
          for (int i=0, k=estlen-1; i<estlen; i++, k--)
            estrev[k] = complementSymbol[(int)estseq[i]];
          estrev[estlen] = 0;

          appendAlignments(B,
                           estrev, dbseq, estlen, dblen,
                           rAligns, rExons,
                           BWD);
        }

        if (globalParams->_findAllExons)
          maskExonsFromGenomic(rExons, dbseq, dbrev, dblen);
      }
    }

  fail:

    //  These are now garbage collected
    //if (fAligns)  free_align(fAligns);
    //if (rAligns)  free_align(rAligns);
    //freeExonList(fExons);
    //freeExonList(rExons);

    fAligns = rAligns = 0L;
    fExons  = rExons  = 0L;

    L->push(B.release());
  } while (globalParams->_findAllExons && pleaseContinueComputing);

 abort:

  delete [] seqStorage;

  return(L);
}






////////////////////////////////////////////////////////////
//
//  Exons
//
////////////////////////////////////////////////////////////


bool
Sim4::checkExonsForOverlaps(Exon *theExons) {
  Exon          *a = theExons;
  Exon          *b = theExons->next_exon;

  while (b && b->toGEN) {
    if ((b->frGEN <= a->toGEN) ||
        (b->frEST <= a->toEST)) {
      return(true);
    }

    a = b;
    b = b->next_exon;
  }

  return(false);
}



void
Sim4::appendExons(sim4polishBuilder &B, Exon *theExons) {
  Exon          *theExon = theExons;

  while (theExon) {
    if (theExon->toGEN) {

#ifdef SPLSCORE
      //  Save the splice score (theExon->splScore);
      //    "%d-%d (%d-%d) <%d-%d-%d> %1.2f %s"
#error I do not know how to save the splice score!
#endif

      char ori = SIM4_INTRON_NONE;

      if ((theExon->next_exon) && (theExon->next_exon->toGEN)) {
        switch (theExon->ori) {
          case 'C':  //  <-
            ori = SIM4_INTRON_NEGATIVE;
            break;
          case 'E':  //  ==
            ori = SIM4_INTRON_GAP;
            break;
          case 'G':  //  ->
            ori = SIM4_INTRON_POSITIVE;
            break;
          case 'N':  //  --
            ori = SIM4_INTRON_AMBIGUOUS;
            break;
          default:
            ori = SIM4_INTRON_ERROR;
            break;
        }
      }

      B.addExon(theExon->frEST, theExon->toEST,
                theExon->frGEN, theExon->toGEN,
                theExon->numMatches,
                theExon->numNs,
                theExon->percentID,
                ori);
    }

    theExon = theExon->next_exon;
  }
}


////////////////////////////////////////////////////////////
//
//  Alignments
//
////////////////////////////////////////////////////////////




void
Sim4::IDISPLAY(sim4polishBuilder &builder,
               char *aString,
               char *bString,
               char *A,
               char *B,
               int   M,
               int   N,
               int  *S,
               int   AP,
               int   BP,
               int   est_strand,
               Exon *exons) {
  Exon *t0;
  register int    i,  j, op;
  int   starti, is_intron=0;

  if ((exons==NULL) || (!exons->toGEN && (exons->next_exon==NULL))) {
    builder.addExonAlignment("Empty exon list; no alignment possible!",
                             "Empty exon list; no alignment possible!");
    return;
  }

  /* find the starting exon for this alignment */
  t0 = exons;
  while (t0 && (((est_strand==2) && ((t0->frGEN!=AP) || (t0->frEST!=BP))) ||
                ((est_strand==1) && ((t0->frGEN!=BP) || (t0->frEST!=AP))))) {
    t0 = t0->next_exon;
  }

  if (!t0) {
    builder.addExonAlignment("Alignment fragment not found; no alignment possible!",
                             "Alignment fragment not found; no alignment possible!");
    return;
  }

  i = j = op = 0;

  starti = (t0->next_exon && t0->next_exon->toGEN) ? (t0->toGEN+1):-1;

  char *a = aString;
  char *b = bString;

#if 0
  fprintf(stderr, "M=%d N=%d\n", M, N);
  fprintf(stderr, "aString=0x%p\nbString=0x%p\n", aString, bString);
#endif

  while (i < M || j < N) {
    *a = *b = 0;
#if 0
    fprintf(stderr, "i=%d < M=%d and j=%d < N=%d\n", i, M, j, N);
    fprintf(stderr, "a=%s\n", aString);
    fprintf(stderr, "b=%s\n", bString);
#endif

    if (op == 0 && *S == 0) {
      op = *S++;
      i++;
      j++;
      if (A[i] == B[j]) {
        *a++ = (char)(A[i] + 'a' - 'A');
        *b++ = (char)(B[j] + 'a' - 'A');
      } else {
        *a++ = A[i];
        *b++ = B[j];
      }
    } else {
      if (op == 0)
        op = *S++; 

      if (op > 0) {
        if (est_strand==2) {
          *a++ = '-';
          *b++ = B[++j];
          op--;
        } else {
          if (j+BP==starti) {
            /* detected intron */
            t0 = t0->next_exon;
            starti=(t0->next_exon && t0->next_exon->toGEN)?(t0->toGEN+1):-1;
            /* print entire exon */
            is_intron = 1;
            j += op;
            op = 0;
          } else {
            *a++ = '-';
            *b++ = B[++j];
            op--;
          }
        }
      } else {
        if (est_strand==1) {
          *a++ = A[++i];
          *b++ = '-';
          op++;
        } else {
          if (i+AP==starti) {
            /* detected intron */
            t0 = t0->next_exon;
            starti=(t0->next_exon && t0->next_exon->toGEN)?(t0->toGEN+1):-1;
            is_intron = 1;
            i += -op;
            op = 0;
          } else {
            *a++ = A[++i];
            *b++ = '-';
            op++;
          }
        }
      }
    }

    if (is_intron || ((i >= M) && (j >= N))) {
      *a   = 0;
      *b   = 0;

      builder.addExonAlignment(aString, bString);

      a = aString;
      b = bString;

      is_intron = 0;
    }
  }
}




void
Sim4::S2A(edit_script *head, int *S) {
  edit_script *tp;
  int *lastS, i;

  tp = head;
  lastS = S;
  while (tp != NULL) {
    if (tp->op_type == SUBSTITUTE) {
      for (i=0; i<tp->num; ++i)
        *lastS++ = 0;
    } else if (tp->op_type == INSERT) {
      *lastS++ = -tp->num;
    } else {     /* DELETE */
      *lastS++ =  tp->num;
    }
    tp = tp->next;
  }
  *(S-1) = (int)(lastS - S);
}




void
Sim4::appendAlignments(sim4polishBuilder &builder,
                       char *s1,
                       char *s2,
                       int l1,
                       int l2, 
                       edit_script_list *Aligns,
                       Exon *Exons, 
                       int match_ori) {

  if (Aligns==NULL)
    return;

  //  Detemine the maximum length of an alignment by finding the
  //  longest exon.
  //
  int     maxAlignmentLength = 0;
  Exon   *theExon = Exons;

  while (theExon) {
    if (theExon->toGEN) {
      if (maxAlignmentLength < (theExon->toGEN - theExon->frGEN + theExon->toEST - theExon->frEST))
        maxAlignmentLength = theExon->toGEN - theExon->frGEN + theExon->toEST - theExon->frEST;
    }

    theExon = theExon->next_exon;
  }

  char *aString = new char [maxAlignmentLength + 4];
  char *bString = new char [maxAlignmentLength + 4];

  for(edit_script_list *aligns = Aligns; aligns; aligns = aligns->next_script) {
    int *S = (int *)ckalloc((2 * aligns->len2 + 1 + 1) * sizeof(int));
    S++;
    S2A(aligns->script, S);
    
    if (match_ori==FWD) {
      IDISPLAY(builder,
               aString,
               bString,
               s1 + aligns->offset2 - 1 - 1,
               s2 + aligns->offset1 - 1 - 1,
               aligns->len2,
               aligns->len1,
               S,
               aligns->offset2,
               aligns->offset1,
               1,
               Exons);
    } else {
      align_reverse(S);
      IDISPLAY(builder,
               aString,
               bString,
               s1 + l1 + 1 - (aligns->offset2 + aligns->len2 - 1) - 1 - 1,
               s2 + l2 + 1 - (aligns->offset1 + aligns->len1 - 1) - 1 - 1,
               aligns->len2,
               aligns->len1,
               S,
               l1 + 1 - (aligns->offset2+aligns->len2 - 1),
               l2 + 1 - (aligns->offset1+aligns->len1 - 1),
               1,
               Exons);
    }
    ckfree(S-1);
  }

  delete [] aString;
  delete [] bString;
}

