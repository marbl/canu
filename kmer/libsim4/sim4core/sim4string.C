#include "sim4.H"
#include "sim4db.H"
#include "libbri.H"

//#define MASKTEST
//#define BUGS
//#define SHOW_OVERLAPPING_EXONS

static void
add_offset_exons(Exon *exons, int offset) {
  Exon *t;

  if (!offset || !(exons)) return;

  t = exons;
  while (t) {
    if (t->toGEN) { t->frEST += offset; t->toEST += offset; }
    t = t->next_exon;
  }
}

static void
add_offset_aligns(edit_script_list *aligns, int offset) {
  edit_script_list *head;
  
  if (!offset || !aligns) return;

  head = aligns;
  while (head) { head->offset2 += offset; head = head->next_script; }
}


char *
appendOutput(char *outstring, char const *appendstring) {
  int   ol=0, al=0;
  char *newstring;
  char *t, *o;

  if (outstring) {
    ol = 0;
    while (outstring[ol])
      ol++;
  }

  if (appendstring) {
    al = 0;
    while (appendstring[al])
      al++;
  }

  t = newstring = new char [ol + al + 1];
  o = outstring;

  if (outstring) {
    while (*o)
      *(t++) = *(o++);
    delete [] outstring;
  }

  if (appendstring)
    while (*appendstring)
      *(t++) = *(appendstring++);

  *t = 0;

  return(newstring);
}


void
Sim4::maskExonsFromGenomic(Exon *theExons,
                           unsigned char *f, unsigned char *r,
                           int l) {
  Exon          *theExon = theExons;

#ifdef MASKTEST
  fprintf(stderr, "BEFORE masking:\n");
  fprintf(stderr, "f='%s'\n", f);
  fprintf(stderr, "r='%s'\n", r);
#endif

  while (theExon) {
    if (theExon->toGEN) {
      for (int i=theExon->frGEN-1; i<theExon->toGEN; i++)
        f[i] = 'N';
      for (int i=l-theExon->frGEN; i>=l-theExon->toGEN; i--)
        r[i] = 'N';
    }

    theExon = theExon->next_exon;
  }

#ifdef MASKTEST
  fprintf(stderr, "AFTER masking:\n");
  fprintf(stderr, "f='%s'\n", f);
  fprintf(stderr, "r='%s'\n", r);
#endif
}





char *
Sim4::run(sim4parameters *s4Params) {
  unsigned int  lineBufferSize = 4096;
  char         *lineBuffer     = new char [lineBufferSize];

  int    dist, match_ori;
  int    pA=0, g_pA=0, f_pA=0, r_pA=0, pT=0, g_pT=0, f_pT=0, r_pT=0;

  Exon   *Exons     = NULL;
  Exon   *rev_Exons = NULL; 

  edit_script_list *Aligns     = NULL;
  edit_script_list *rev_Aligns = NULL;

  int     matchesPrinted = 0;

  //  Make absolutely sure that the database sequence start and end
  //  positions are within the valid range.  Ideally, this should be
  //  checked by whatever generates the input, but it probably isn't.
  //
  //  If the end position is too big, make it the same as the sequence
  //  length.
  //
  //  If the start position is bigger than the (corrected) end
  //  position, make it 100K less than the end position.
  //
  if (s4Params->_dbHi > s4Params->getDBlength(s4Params->_dbIdx))
    s4Params->_dbHi = s4Params->getDBlength(s4Params->_dbIdx);

  if (s4Params->_dbLo > s4Params->_dbHi)
    if (s4Params->_dbHi > 100000)
      s4Params->_dbLo = s4Params->_dbHi - 100000;
    else
      s4Params->_dbLo = 0;


  int              dblen     = s4Params->_dbHi - s4Params->_dbLo;
  unsigned char   *dbseq;
  unsigned char   *dbrev;
  unsigned char   *dbseqorig   = s4Params->getDBsequence(s4Params->_dbIdx);
  bool             dbWasMasked = false;

  int              estlen;
  unsigned char   *estseq;
  unsigned char   *estrev;

  unsigned int     pleaseContinueComputing = false;

  char            *outputString    = 0L;
  char const      *strandIndicator = 0L;


  //  Allocate space for temporary sequence storage.  We need
  //  to allocate space for two copies of the database, and space
  //  for the longest EST (in case we need to print it out
  //  reverse complemented).
  //
  unsigned char   *seqStorage     = 0L;
  unsigned int     seqStorageSize = 0;

  for (unsigned int e=0; e<s4Params->_numESTs; e++)
    if (seqStorageSize < s4Params->getESTlength(s4Params->_ESTlist[e]))
      seqStorageSize = s4Params->getESTlength(s4Params->_ESTlist[e]);

  seqStorageSize += dblen + dblen + 6;
  seqStorage      = new unsigned char [seqStorageSize];

  //  Original, forward, reverse, cdna
  //
  //fprintf(stderr, "Allocating %d bytes for sequence storage\n", seqStorageSize);

  dbseq  = seqStorage;
  dbrev  = seqStorage + dblen + 2;
  estrev = seqStorage + dblen + 2 + dblen + 2;


  //  Prepare the database sequence
  //
  //  Trimming to the correct range
  //  Convert to uppercase
  //  Reverse complement
  //
  for (unsigned int i=0, j=s4Params->_dbLo, k=dblen-1; j<s4Params->_dbHi; i++, j++, k--) {
    dbseq[i] = dbseqorig[j];
    dbrev[k] = complementSymbol[dbseq[i]];
  }
  dbseq[dblen] = 0;
  dbrev[dblen] = 0;

  for (unsigned int e=0; e<s4Params->_numESTs; e++) {
    sim4_stats_t   st, rev_st;

    estseq = s4Params->getESTsequence(s4Params->_ESTlist[e]);
    estlen = s4Params->getESTlength(s4Params->_ESTlist[e]);

    g_pT = g_pA = 0;

    if (dbParams._ignorePolyTails) {
      get_polyAT(estseq, estlen, &g_pT, &g_pA);
      //fprintf(stderr, "Masked: %d %d\n", g_pT, g_pA);
    }

    if (estlen - g_pA - g_pT <= 0)
      continue;

    //  If the database was masked, restore it to the original state.
    //
    if (dbWasMasked) {
      dbWasMasked = false;
      for (unsigned int i=0, j=s4Params->_dbLo, k=dblen-1; j<s4Params->_dbHi; i++, j++, k--) {
        dbseq[i] = dbseqorig[j];
        dbrev[k] = complementSymbol[dbseq[i]];
      }
      dbseq[dblen] = 0;
      dbrev[dblen] = 0;
    }

    matchesPrinted = 0;

    do {
      memset(&st,     0, sizeof(sim4_stats_t));
      memset(&rev_st, 0, sizeof(sim4_stats_t));

      bld_table(estseq - 1 + g_pT, estlen - g_pA - g_pT, DEFAULT_W, INIT);

      if (s4Params->_doForward) {
        Aligns = SIM4(dbseq, estseq + g_pT,
                      dblen, estlen - g_pT - g_pA,
                      &dist,
                      &Exons,
                      &f_pA,
                      &f_pT,
                      &st);

        //  Continued from util.C :: slide_intron()
        //
        //  If we are forcing the strand prediction, and we are still unknown,
        //  set the strand prediction to the match orientation. Since this 
        //  will be reversed later on, set it to FWD here.
        //
        if ((dbParams._forceStrandPrediction) && (st.orientation == BOTH))
          st.orientation = FWD;

#if ABORT_EXPENSIVE
        if (st.tooManyMSPs) { 
          unsigned int    lbs = (strlen((char *)s4Params->getESTheader(s4Params->_ESTlist[e])) +
                                 strlen((char *)s4Params->getDBheader(s4Params->_dbIdx)) +
                                 2048);

          if (lineBufferSize < lbs) {
            delete [] lineBuffer;
            lineBufferSize  = lbs;
            lineBuffer      = new char [lineBufferSize];
          }

          sprintf(lineBuffer, "sim4begin\n%d[%d-0-0] %d[%d-%d] <0-0-0-forward-intractable>\n",
                  s4Params->_ESTlist[e],
                  estlen,
                  s4Params->_dbIdx,
                  s4Params->_dbLo,
                  s4Params->_dbHi);
          outputString = appendOutput(outputString, lineBuffer); 

          if (dbParams._includeDefLine) {
            sprintf(lineBuffer, "edef=%s\nddef=%s\n",
                    s4Params->getESTheader(s4Params->_ESTlist[e]),
                    s4Params->getDBheader(s4Params->_dbIdx));
            outputString = appendOutput(outputString, lineBuffer); 
          }

          sprintf(lineBuffer, "1-%d (1-%d) <%d-0-0>\nsim4end\n",
                  estlen,
                  s4Params->_dbHi - s4Params->_dbLo,
                  st.numberOfMatches);
          outputString = appendOutput(outputString, lineBuffer); 

          goto abort;
        }
#endif
      }

      if (s4Params->_doReverse) {
        rev_Aligns = SIM4(dbrev, estseq + g_pT,
                          dblen, estlen - g_pT - g_pA,
                          &dist,
                          &rev_Exons,
                          &r_pA,
                          &r_pT,
                          &rev_st);

        //  Continued from util.C :: slide_intron()
        //
        //  If we are forcing the strand prediction, and we are still unknown,
        //  set the strand prediction to the match orientation. 
        //
        if ((dbParams._forceStrandPrediction) && (rev_st.orientation == BOTH))
          rev_st.orientation = FWD;

#if ABORT_EXPENSIVE
        if (rev_st.tooManyMSPs) { 
          unsigned int    lbs = (strlen((char *)s4Params->getESTheader(s4Params->_ESTlist[e])) +
                                 strlen((char *)s4Params->getDBheader(s4Params->_dbIdx)) +
                                 2048);

          if (lineBufferSize < lbs) {
            delete [] lineBuffer;
            lineBufferSize  = lbs;
            lineBuffer      = new char [lineBufferSize];
          }

          sprintf(lineBuffer, "sim4begin\n%d[%d-0-0] %d[%d-%d] <0-0-0-complement-intractable>\n",
                  s4Params->_ESTlist[e],
                  estlen,
                  s4Params->_dbIdx,
                  s4Params->_dbLo,
                  s4Params->_dbHi);
          outputString = appendOutput(outputString, lineBuffer); 

          if (dbParams._includeDefLine) {
            sprintf(lineBuffer, "edef=%s\nddef=%s\n",
                    s4Params->getESTheader(s4Params->_ESTlist[e]),
                    s4Params->getDBheader(s4Params->_dbIdx));
            outputString = appendOutput(outputString, lineBuffer); 
          }

          sprintf(lineBuffer, "1-%d (1-%d) <%d-0-0>\nsim4end\n",
                  estlen,
                  s4Params->_dbHi - s4Params->_dbLo,
                  rev_st.numberOfMatches);
          outputString = appendOutput(outputString, lineBuffer); 

          goto abort;
        }
#endif
      }


      if (st.numberOfMatches >= rev_st.numberOfMatches) {
        match_ori = FWD;

        if (s4Params->_strandIndicator == 0L) {
          switch (st.orientation) {
          case FWD:
            strandIndicator = "forward";
            break;
          case BWD:
            strandIndicator = "reverse";
            break;
          default:
            strandIndicator = "unknown";
            break;
          }
        } else {
          strandIndicator = s4Params->_strandIndicator;
        }

        if (dbParams._ignorePolyTails) {
          add_offset_exons(Exons, g_pT);  
          add_offset_aligns(Aligns, g_pT);
        }

        if (rev_Exons) { free_list(rev_Exons); rev_Exons = NULL; }
        if (rev_Aligns) { free_align(rev_Aligns); rev_Aligns = NULL; }

        pT = g_pT + f_pT;
        pA = g_pA + f_pA;

        if (Exons) {
          char *result = checkExonsForOverlaps(Exons);
          if (result) {
#ifdef SHOW_OVERLAPPING_EXONS
            unsigned int    lbs = (strlen((char *)s4Params->getESTheader(s4Params->_ESTlist[e])) +
                                   strlen((char *)s4Params->getDBheader(s4Params->_dbIdx)) +
                                   2048);

            if (lineBufferSize < lbs) {
              delete [] lineBuffer;
              lineBufferSize  = lbs;
              lineBuffer      = new char [lineBufferSize];
            }

            if (dbParams._includeDefLine) {
              sprintf(lineBuffer, "sim4begin\n%d[%d-0-0] %d[%d-%d] <0-0-0-forward-failed>\nedef=%s\nddef=%s\n%ssim4end\n",
                      s4Params->_ESTlist[e],
                      estlen,
                      s4Params->_dbIdx,
                      s4Params->_dbLo,
                      s4Params->_dbHi,
                      s4Params->getESTheader(s4Params->_ESTlist[e]),
                      s4Params->getDBheader(s4Params->_dbIdx),
                      result);
            } else {
              sprintf(lineBuffer, "sim4begin\n%d[%d-0-0] %d[%d-%d] <0-0-0-forward-failed>\n%ssim4end\n",
                      s4Params->_ESTlist[e],
                      estlen,
                      s4Params->_dbIdx,
                      s4Params->_dbLo,
                      s4Params->_dbHi,
                      result);
            }

            outputString = appendOutput(outputString, lineBuffer);
#endif
            delete [] result;
            goto abort;
          }
        }
      } else {
        match_ori = BWD;

        if (s4Params->_strandIndicator == 0L) {
          switch (rev_st.orientation) {
            case FWD:
              strandIndicator = "reverse";
              break;
            case BWD:
              strandIndicator = "forward";
              break;
            default:
              strandIndicator = "unknown";
              break;
          }
        } else {
          strandIndicator = s4Params->_strandIndicator;
        }

        if (dbParams._ignorePolyTails) {
          add_offset_exons(rev_Exons, g_pT); 
          add_offset_aligns(rev_Aligns, g_pT);
        }

        if (rev_Aligns && rev_Aligns->next_script)
          script_flip_list(&rev_Aligns);

        //  This used to be right before appendExons() in
        //  the reverse match section, but we need it
        //  before we test for overlapping exons
        //
        complement_exons(&rev_Exons, dblen, estlen);

        if (Exons) { free_list(Exons);  Exons = NULL; }
        if (Aligns) { free_align(Aligns); Aligns = NULL; }

        pT = g_pT + r_pT;
        pA = g_pA + r_pA;

        if (rev_Exons) {
          char *result = checkExonsForOverlaps(rev_Exons);
          if (result) {
#ifdef SHOW_OVERLAPPING_EXONS
            unsigned int    lbs = (strlen((char *)s4Params->getESTheader(s4Params->_ESTlist[e])) +
                                   strlen((char *)s4Params->getDBheader(s4Params->_dbIdx)) +
                                   2048);

            if (lineBufferSize < lbs) {
              delete [] lineBuffer;
              lineBufferSize  = lbs;
              lineBuffer      = new char [lineBufferSize];
            }


            if (dbParams._includeDefLine) {
              sprintf(lineBuffer, "sim4begin\n%d[%d-0-0] %d[%d-%d] <0-0-0-reverse-failed>\nedef=%s\nddef=%s\n%ssim4end\n",
                      s4Params->_ESTlist[e],
                      estlen,
                      s4Params->_dbIdx,
                      s4Params->_dbLo,
                      s4Params->_dbHi,
                      s4Params->getESTheader(s4Params->_ESTlist[e]),
                      s4Params->getDBheader(s4Params->_dbIdx),
                      result);
            } else {
              sprintf(lineBuffer, "sim4begin\n%d[%d-0-0] %d[%d-%d] <0-0-0-reverse-failed>\n%ssim4end\n",
                      s4Params->_ESTlist[e],
                      estlen,
                      s4Params->_dbIdx,
                      s4Params->_dbLo,
                      s4Params->_dbHi,
                      result);
            }

            outputString = appendOutput(outputString, lineBuffer);
#endif
            delete [] result;
            goto abort;
          }
        }

      }







      int     nmatches  = (match_ori==FWD) ? st.numberOfMatches  : rev_st.numberOfMatches;
      double  coverage  = (double)nmatches / (double)estlen;
      int     percentid = (match_ori==FWD) ? st.percentID        : rev_st.percentID;


      //  Is this match decent?
      //
      pleaseContinueComputing = ((coverage  >= dbParams._minCoverage) &&
                                 (percentid >= dbParams._minPercentExonIdentity) &&
                                 (nmatches  >= dbParams._minCoverageLength) &&
                                 (nmatches  > 0));

      //  If we're supposed to print at least _alwaysReport things,
      //  and we found a match, keep going.
      //
      if ((matchesPrinted < dbParams._alwaysReport) && (nmatches > 0))
        pleaseContinueComputing = true;

      //  However, if we have printed enough stuff, and the last one is
      //  below the thresholds, stop.
      //
      if ((matchesPrinted >= dbParams._alwaysReport) &&
          ((coverage  < dbParams._minCoverage) ||
           (percentid < dbParams._minPercentExonIdentity)))
        pleaseContinueComputing = false;


      if (pleaseContinueComputing) {
        matchesPrinted++;

        outputString = appendOutput(outputString, "sim4begin\n");
        if (dbParams._includeDefLine) {
          unsigned int    lbs = (strlen((char *)s4Params->getESTheader(s4Params->_ESTlist[e])) +
                                 strlen((char *)s4Params->getDBheader(s4Params->_dbIdx)) +
                                 2048);

          if (lineBufferSize < lbs) {
            delete [] lineBuffer;
            lineBufferSize  = lbs;
            lineBuffer      = new char [lineBufferSize];
          }
          
          sprintf(lineBuffer, "%d[%d-%d-%d] %d[%d-%d] <%d-%d-%d-%s-%s>\nedef=%s\nddef=%s\n",
                  s4Params->_ESTlist[e],
                  estlen,
                  pA,
                  pT,

                  s4Params->_dbIdx,
                  s4Params->_dbLo,
                  s4Params->_dbHi,

                  (match_ori==FWD) ? st.numberOfMatches  : rev_st.numberOfMatches,
                  (match_ori==FWD) ? st.numberOfNs       : rev_st.numberOfNs,
                  (match_ori==FWD) ? st.percentID        : rev_st.percentID,
                  (match_ori==FWD) ? "forward"           : "complement",
                  strandIndicator,

                  s4Params->getESTheader(s4Params->_ESTlist[e]),
                  s4Params->getDBheader(s4Params->_dbIdx));
        } else {
          sprintf(lineBuffer, "%d[%d-%d-%d] %d[%d-%d] <%d-%d-%d-%s-%s>\n",
                  s4Params->_ESTlist[e],
                  estlen,
                  pA,
                  pT,

                  s4Params->_dbIdx,
                  s4Params->_dbLo,
                  s4Params->_dbHi,

                  (match_ori==FWD) ? st.numberOfMatches  : rev_st.numberOfMatches,
                  (match_ori==FWD) ? st.numberOfNs       : rev_st.numberOfNs,
                  (match_ori==FWD) ? st.percentID        : rev_st.percentID,
                  (match_ori==FWD) ? "forward"           : "complement",
                  strandIndicator);
        }
        outputString = appendOutput(outputString, lineBuffer);
        
        if (match_ori == FWD) {
          outputString = appendExons(outputString, Exons);
          if (dbParams._printAlignments) {
            outputString = appendAlignments(outputString,
                                            estseq, dbseq, estlen, dblen,
                                            &Aligns, Exons,
                                            FWD);
          }

          if (dbParams._findAllExons) {
            maskExonsFromGenomic(Exons, dbseq, dbrev, dblen);
            dbWasMasked = true;
          }
        } else {
          outputString = appendExons(outputString, rev_Exons);

          if (dbParams._printAlignments) {
            for (int i=0, k=estlen-1; i<estlen; i++, k--)
              estrev[k] = complementSymbol[estseq[i]];
            estrev[estlen] = 0;

            outputString = appendAlignments(outputString,
                                            estrev, dbseq, estlen, dblen,
                                            &rev_Aligns, rev_Exons,
                                            BWD);
          }

          if (dbParams._findAllExons) {
            maskExonsFromGenomic(rev_Exons, dbseq, dbrev, dblen);
            dbWasMasked = true;
          }
        }

        outputString = appendOutput(outputString, "sim4end\n");
      }

      if (Aligns)     { free_align(Aligns);     Aligns = NULL; }
      if (rev_Aligns) { free_align(rev_Aligns); rev_Aligns = NULL; }
      if (Exons)      { free_list(Exons);       Exons = NULL; }
      if (rev_Exons)  { free_list(rev_Exons);   rev_Exons = NULL; }
    } while (dbParams._findAllExons && pleaseContinueComputing);
  }

#if ABORT_EXPENSIVE
 abort:
#endif

  delete [] seqStorage;
  delete [] lineBuffer;

  if (outputString == 0L) {
    outputString = new char [1];
    outputString[0] = 0;
  }

  return(outputString);
}



char *
Sim4::checkExonsForOverlaps(Exon *theExons) {
  Exon          *a = theExons;
  Exon          *b = theExons->next_exon;

  while (b && b->toGEN) {
    if ((b->frGEN <= a->toGEN) ||
        (b->frEST <= a->toEST)) {
      return(appendExons(0L, theExons));
    }

    a = b;
    b = b->next_exon;
  }

  return(0L);
}


char *
Sim4::appendExons(char *outstring, Exon *theExons) {
  Exon          *theExon = theExons;
  char          *exonString = new char [1024 * 1024];
  char          *exonPrint  = exonString;

  *exonPrint = 0;

  while (theExon) {
    if (theExon->toGEN) {
#ifdef SPLSCORE
      if ((theExon->next_exon) && (theExon->next_exon->toGEN)) {
        sprintf(exonPrint, "%d-%d (%d-%d) <%d-%d-%d> %1.2f",
                theExon->frEST, theExon->toEST,
                theExon->frGEN, theExon->toGEN,
                theExon->numMatches,
                theExon->numNs,
                theExon->percentID,
                theExon->splScore);
      } else {
        sprintf(exonPrint, "%d-%d (%d-%d) <%d-%d-%d>",
                theExon->frEST, theExon->toEST,
                theExon->frGEN, theExon->toGEN,
                theExon->numMatches,
                theExon->numNs,
                theExon->percentID);
      }
#else
      sprintf(exonPrint, "%d-%d (%d-%d) <%d-%d-%d>",
              theExon->frEST, theExon->toEST,
              theExon->frGEN, theExon->toGEN,
              theExon->numMatches,
              theExon->numNs,
              theExon->percentID);
#endif

      while (*exonPrint)
        exonPrint++;

      if ((theExon->next_exon) && (theExon->next_exon->toGEN)) {
        switch (theExon->ori) {
        case 'C':  //  <-
          *(exonPrint++) = ' ';
          *(exonPrint++) = '<';
          *(exonPrint++) = '-';
          break;
        case 'E':  //  ==
          *(exonPrint++) = ' ';
          *(exonPrint++) = '=';
          *(exonPrint++) = '=';
          break;
        case 'G':  //  ->
          *(exonPrint++) = ' ';
          *(exonPrint++) = '-';
          *(exonPrint++) = '>';
          break;
        case 'N':  //  --
          *(exonPrint++) = ' ';
          *(exonPrint++) = '-';
          *(exonPrint++) = '-';
          break;
        default:
          sprintf(exonPrint, " appendExon: Inconsistent exon orientation '%c'.", theExon->ori);
          while (*exonPrint)
            exonPrint++;
        }  
      }

      *(exonPrint++) = '\n';
      *(exonPrint) = 0;
    }

    theExon = theExon->next_exon;
  }

  outstring = appendOutput(outstring, exonString);

#ifdef BUGS
  fprintf(stderr, "----------------------------------------\n");
  fprintf(stderr, exonString);
  fprintf(stderr, "----------");
  fprintf(stderr, outstring);
#endif

  delete [] exonString;

  return(outstring);
}








char*
Sim4::IDISPLAY(char *outputstring,
               char *aString,
               char *bString,
               uchar A[],
               uchar B[],
               int M, int N,
               int S[],
               int AP, int BP,
               int est_strand,
               Exon *exons) {
  Exon *t0;
  register int    i,  j, op;
  int   starti, is_intron=0;

#ifdef BUGS
  fprintf(stderr, "Helo I'm IDISPLAY\n");
#endif

  if ((exons==NULL) || (!exons->toGEN && (exons->next_exon==NULL))) {
    outputstring = appendOutput(outputstring, "Empty exon list?\n");
    return(outputstring);
  }

  /* find the starting exon for this alignment */
  t0 = exons;
#ifdef BUGS
  fprintf(stderr, "t0=0x%016lx %d-%d (%d-%d)\n", t0, t0->frEST, t0->toEST, t0->frGEN, t0->toGEN);
#endif
  while (t0 && (((est_strand==2) && ((t0->frGEN!=AP) || (t0->frEST!=BP))) ||
                ((est_strand==1) && ((t0->frGEN!=BP) || (t0->frEST!=AP))))) {
    t0 = t0->next_exon;
#ifdef BUGS
    fprintf(stderr, "t0=0x%016lx %d-%d (%d-%d)\n", t0, t0->frEST, t0->toEST, t0->frGEN, t0->toGEN);
#endif
  }

  if (!t0) {
    outputstring = appendOutput(outputstring, "Alignment fragment not found?\n");
    return(outputstring);
  }

#ifdef BUGS
  fprintf(stderr, "t0=0x%016lx AP=%d - BP=%d\n", t0, AP, BP);
#endif

  i = j = op = 0;

  starti = (t0->next_exon && t0->next_exon->toGEN) ? (t0->toGEN+1):-1;

  char *a = aString;
  char *b = bString;

#if 0
  fprintf(stdout, "A=%s\n", A);
  fprintf(stdout, "B=%s\n", B);
#endif

  while (i < M || j < N) {
#ifdef BUGS
    fprintf(stderr, "%8d < %8d || %8d < %8d  :  ", i, M, j, N);
#endif

    if (op == 0 && *S == 0) {
      op = *S++;
      i++;
      j++;
      if (A[i] == B[j]) {
        *a++ = A[i] + 'a' - 'A';
        *b++ = B[j] + 'a' - 'A';
      } else {
        *a++ = A[i];
        *b++ = B[j];
      }
    } else {
      if (op == 0) op = *S++; 
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

#ifdef BUGS
    fprintf(stderr, "%c - %c\n", *(a-1), *(b-1));
#endif

    if (is_intron || ((i >= M) && (j >= N))) {
      *a++ = '\n';
      *a   = 0;

#ifdef BUGS
      fprintf(stdout, "aString=(%5d)%s", strlen(aString), aString);
#endif

      *b++ = '\n';
      *b   = 0;

#ifdef BUGS
      fprintf(stdout, "bString=(%5d)%s", strlen(bString), bString);
#endif

      outputstring = appendOutput(outputstring, aString);
      outputstring = appendOutput(outputstring, bString);

      a = aString;
      b = bString;

      is_intron = 0;
    }
  }

  return(outputstring);
}




void
Sim4::S2A(edit_script *head, int *S)
{
  edit_script *tp;
  int *lastS, i;

  tp = head;
  lastS = S;
  while (tp != NULL) {
    if (tp->op_type == SUBSTITUTE) {
      for (i=0; i<tp->num; ++i) *lastS++ = 0;
    } else if (tp->op_type == INSERT) {
      *lastS++ = -tp->num;
    } else {     /* DELETE */
      *lastS++ =  tp->num;
    }
    tp = tp->next;
  }
  *(S-1) = lastS - S;
}




char*
Sim4::appendAlignments(char *outputstring,
                       uchar *s1,
                       uchar *s2,
                       int l1,
                       int l2, 
                       edit_script_list **Aligns,
                       Exon *Exons, 
                       int match_ori) {
  int *S;
  edit_script_list *head, *aligns;

  if (*Aligns==NULL)
    return(outputstring);

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


  aligns = *Aligns;
  while (aligns!=NULL) {
    head = aligns;
    aligns = aligns->next_script; 

    S = (int *)ckalloc((2*head->len2+1+1)*sizeof(int));
    S++;            
    S2A(head->script, S);
    Free_script(head->script);
    
    if (match_ori==FWD) {
      outputstring = IDISPLAY(outputstring,
                              aString,
                              bString,
                              s1 + head->offset2 - 1 - 1,
                              s2 + head->offset1 - 1 - 1,
                              head->len2,
                              head->len1,
                              S,
                              head->offset2,
                              head->offset1,
                              1,
                              Exons);
    } else {
      align_reverse(S);
      outputstring = IDISPLAY(outputstring,
                              aString,
                              bString,
                              s1 + l1 + 1 - (head->offset2 + head->len2 - 1) - 1 - 1,
                              s2 + l2 + 1 - (head->offset1 + head->len1 - 1) - 1 - 1,
                              head->len2,
                              head->len1,
                              S,
                              l1 + 1 - (head->offset2+head->len2 - 1),
                              l2 + 1 - (head->offset1+head->len1 - 1),
                              1,
                              Exons);
    }
    free(S-1);
    free(head);
  }
  *Aligns = NULL;

  //fprintf(stderr, "All done print aligns\n");

  delete [] aString;
  delete [] bString;

  return(outputstring);
}

