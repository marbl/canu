#include "snapper2.H"


#define MATCH            0
#define GAPA             1
#define GAPB             2
#define STOP             3

#define MATCHSCORE       3
#define GAPSCORE        -2
#define MISMATCHSCORE   -2


void
reverse(char *a, char *b, int len) {
  char   c=0;
  char  *s=a,  *S=a+len-1;
  char  *q=b,  *Q=b+len-1;

  while (s < S) {
    c    = *s;
    *s++ =  *S;
    *S-- =  c;

    c    = *q;
    *q++ = *Q;
    *Q-- =  c;
  }
}



class dpMatch {
public:
  dpMatch() {
    matches  = 0;
    alignLen = 0;

    begI = begJ = 0;
    endI = endJ = 0;
    lenA = lenB = 0;
  };

  int      matches;
  int      alignLen;

  int      begI, begJ;
  int      endI, endJ;
  int      lenA, lenB;

  char    *alignA;
  char    *alignB;
};


class dpMatrix {
private:
  typedef struct {
    unsigned int  score  : 30;
    unsigned int  action : 2;
  } dpCell;

public:
  dpMatrix() {
    aMax = 0;
    bMax = 0;

    alignA = 0L;
    alignB = 0L;
    matrix = 0L;
  };

  ~dpMatrix() {
    delete [] alignA;
    delete [] alignB;
    delete [] matrix;
  };

  void   dpMatrixInit(int lenA, int lenB) {

    if ((aMax <= lenA) || (bMax <= lenB)) {
      delete [] alignA;
      delete [] alignB;
      delete [] matrix;

      aMax = lenA + 1;
      bMax = lenB + 1;

      fprintf(stderr, "dpMatrix-- reallocate to "u32bitFMT" x "u32bitFMT"\n", aMax, bMax);

      alignA = new char   [aMax + bMax + 1];
      alignB = new char   [bMax + bMax + 1];
      matrix = new dpCell [aMax * bMax];
    }

    int i, j, p = 0;

    for (i=0; i<lenA+1; i++) {
      matrix[p].score = 1 << 29;
      matrix[p].action = STOP;
      p += bMax;
    }

    p = 0;
    for (j=0; j<lenB+1; j++) {
      matrix[p].score = 1 << 29;
      matrix[p].action = STOP;
      p++;
    }
  };

  int    dpMatrixCellGetScore(int a, int b) {
    return(matrix[a * bMax + b].score);
  };

  int    dpMatrixCellGetAction(int a, int b) {
    return(matrix[a * bMax + b].action);
  };

  void   dpMatrixCellSet(int a, int b, int score, int action) {
    dpCell  x;
    x.score  = score;
    x.action = action;
    matrix[a * bMax + b] = x;
  };

  dpMatch *dpAlign(char *stringA, int lenA,
                   char *stringB, int lenB,
                   dpMatch *match);
private:
  int      aMax;
  int      bMax;

  char    *alignA;
  char    *alignB;
  dpCell  *matrix;
};


dpMatch *
dpMatrix::dpAlign(char     *stringA,  int lenA,
                  char     *stringB,  int lenB,
                  dpMatch  *match) {

  int   i, j;

  dpMatrixInit(lenA, lenB);

  int   scoreMax  = 0;

  int   begI=0, endI=0, curI=0;
  int   begJ=0, endJ=0, curJ=0;

  for (i=1; i<=lenA; i++){
    for (j=1; j<=lenB; j++){

      //  Pick the max of these

      int ul = dpMatrixCellGetScore(i-1, j-1) + ((stringA[i-1] == stringB[j-1]) ? MATCHSCORE : MISMATCHSCORE);
      int lf = dpMatrixCellGetScore(i-1, j) + GAPSCORE;
      int up = dpMatrixCellGetScore(i, j-1) + GAPSCORE;

      // (i,j) is the beginning of a subsequence, our default behavior
      int sc = 1 << 29;
      int ac = STOP;

      if (sc < ul) {
        sc = ul;
        ac = MATCH;
      }

      if (sc < lf) {
        sc = lf;
        ac = GAPB;
      }

      if (sc < up) {
        sc = up;
        ac = GAPA;
      }

      dpMatrixCellSet(i, j, sc, ac);

      if (scoreMax < sc) {
        scoreMax  = sc;
        endI = curI = i;
        endJ = curJ = j;
      }
    }
  }

  int  alignLen  = 0;
  int  matches   = 0;
  int  terminate = 0;

  while (terminate == 0) {
    switch (dpMatrixCellGetAction(curI, curJ)) {
      case STOP:
        terminate = 1;
        break;
      case MATCH:
        alignA[alignLen] = stringA[curI-1];
        alignB[alignLen] = stringB[curJ-1];

        if (alignA[alignLen] == alignB[alignLen]) {
          alignA[alignLen] = tolower(alignA[alignLen]);
          alignB[alignLen] = tolower(alignB[alignLen]);
          matches++;
        } else {
          alignA[alignLen] = toupper(alignA[alignLen]);
          alignB[alignLen] = toupper(alignB[alignLen]);
        }

        curI--;
        curJ--;
        alignLen++;
        break;
      case GAPA:
        alignA[alignLen] = '-';
        alignB[alignLen] = stringB[curJ-1];
        curJ--;
        alignLen++;
        break;
      case GAPB:
        alignA[alignLen] = stringA[curI-1];
        alignB[alignLen] = '-';
        curI--;
        alignLen++;
        break;
    }
  }

  begI = curI;
  begJ = curJ;

  alignA[alignLen] = 0;
  alignB[alignLen] = 0;

  reverse(alignA, alignB, alignLen);

  match->matches  = matches;
  match->alignLen = alignLen;
  match->begI     = begI;
  match->begJ     = begJ;
  match->endI     = endI;
  match->endJ     = endJ;
  match->lenA     = lenA;
  match->lenB     = lenB;

#warning alignA and alignB aliases to dpMatrix
  match->alignA   = alignA;
  match->alignB   = alignB;

  return(match);
}








char*
doPolishDP(searcherState       *state,
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

  //  Move these to searcherState!
  dpMatch   match;
  dpMatrix  matrix;

  for (u32bit h=0; h<theHitsLen; h++) {

    //  If the hit was discarded, move along.
    //
    if (theHits[h]._status & AHIT_DISCARDED)
      continue;

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

    seqInCore            *QRYseq = seq;
    seqInCore            *GENseq = genome->getSequenceInCore(theHits[h]._dsIdx);
    u32bit                GENlo  = theHits[h]._dsLo;
    u32bit                GENhi  = theHits[h]._dsHi;

    char  *q = QRYseq->sequence();
    char  *g = GENseq->sequence() + GENlo;

    if (GENhi > GENseq->sequenceLength())
      GENhi = GENseq->sequenceLength();

    u32bit   qlen = seq->sequenceLength();
    u32bit   glen = GENhi - GENlo;

    bool     doForward =  theHits[h]._status & AHIT_DIRECTION_MASK;
    bool     doReverse = !doForward;

    if (doReverse) {
      reverseComplementSequence(q, qlen);
    }

    fprintf(stderr, "align QRYlen="u32bitFMT" GEN="u32bitFMT"-"u32bitFMT" GENlen="u32bitFMT"\n",
            qlen, GENlo, GENhi, glen);

    if ((qlen * 3 > glen) && ((qlen / 1024) * (glen / 1024) < 4 * 1024))
      matrix.dpAlign(q, qlen, g, glen, &match);

    if (doReverse) {
      reverseComplementSequence(q, qlen);

      u32bit x = match.begI;
      match.begI = qlen - match.endI;
      match.endI = qlen - x;
    }


    //  Build the proper match if it's even remotely good
    //
    if (match.matches > 0) {
      sim4polish      p;
      sim4polishExon  e;

      theHits[h]._status |= AHIT_VERIFIED;

      p.estID    = QRYseq->getIID();
      p.estLen   = QRYseq->sequenceLength();
      p.estPolyA = 0;
      p.estPolyT = 0;

      p.genID            = GENseq->getIID();
      p.genRegionOffset  = GENlo;
      p.genRegionLength  = GENhi - GENlo;

      p.numMatches   = match.matches;
      p.numMatchesN  = 0;
      p.numCovered   = match.endI - match.begI;
      
      p.percentIdentity  = 0;
      p.querySeqIdentity = 0;

      p.matchOrientation  = (doReverse) ? SIM4_MATCH_COMPLEMENT : SIM4_MATCH_FORWARD;
      p.strandOrientation = SIM4_STRAND_UNKNOWN;

      p.comment    = NULL;
      p.estDefLine = QRYseq->header();
      p.genDefLine = GENseq->header();

      p.numExons = 1;
      p.exons    = &e;

      e.estFrom           = match.begI + 1;
      e.estTo             = match.endI;
      e.genFrom           = match.begJ + GENlo + 1;
      e.genTo             = match.endJ + GENlo;
      e.numMatches        = match.matches;
      e.numMatchesN       = 0;
      e.percentIdentity   = 0;
      e.intronOrientation = SIM4_INTRON_NONE;
      e.estAlignment      = match.alignA;
      e.genAlignment      = match.alignB;

      s4p_updateAlignmentScores(&p);

      //  Save it if it is truely good.
      if ((p.percentIdentity  >= config._minMatchIdentity) &&
          (p.querySeqIdentity >= config._minMatchCoverage)) {
        char *pstr = s4p_polishToString(&p);

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

        //  Save the best scores
        //
        u32bit pi = p.percentIdentity;
        u32bit pc = p.querySeqIdentity;

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
      }
    }

    state->polished++;
  }  //  over all hits

  state->polishTime += getTime() - startTime;

  return(out);
}

