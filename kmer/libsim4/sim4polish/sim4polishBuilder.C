#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "bio++.H"
#include "sim4polishBuilder.H"


sim4polishBuilder::sim4polishBuilder() {
  it = 0L;

  exPos = 0;
  exMax = 32;
  exAli = 0;
  ex    = new sim4polishExon * [exMax];

  for (u32bit i=0; i<exMax; i++)
    ex[i] = 0L;
}

sim4polishBuilder::~sim4polishBuilder() {
  delete it;

  for (u32bit i=0; i<exMax; i++)
    delete ex[i];

  delete [] ex;
}


void
sim4polishBuilder::create(u32bit estid, u32bit estlen,
                          u32bit genid, u32bit genlo, u32bit genhi) {

  //  Someone didn't call release()!!
  //
  if (it) {
    fprintf(stderr, "sim4polishBuilder::create()-- WARNING:  release() not called.  Polish killed.\n");
    delete it;
  }

  it = new sim4polish;

  it->_estID    = estid;
  it->_estLen   = estlen;
  it->_estPolyA = 0;
  it->_estPolyT = 0;

  it->_genID           = genid;
  it->_genRegionOffset = genlo;
  it->_genRegionLength = genhi - genlo;

  it->_numMatches        = 0;
  it->_numMatchesN       = 0;
  it->_numCovered        = 0;
  it->_percentIdentity   = 0;
  it->_querySeqIdentity  = 0;
  it->_matchOrientation  = SIM4_MATCH_ERROR;
  it->_strandOrientation = SIM4_STRAND_ERROR;

  it->_comment    = 0L;
  it->_estDefLine = 0L;
  it->_genDefLine = 0L;

  it->_numExons = 0;
  it->_exons    = 0L;
}

void
sim4polishBuilder::setPolyTails(u32bit pa, u32bit pt) {
  it->_estPolyA = pa;
  it->_estPolyT = pt;
}

void
sim4polishBuilder::setESTdefline(char *defline) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setESTdefline()-- no polish to build; create() not called\n");
    return;
  }
  delete [] it->_estDefLine;
  it->_estDefLine = new char [strlen(defline) + 1];
  memcpy(it->_estDefLine, defline, sizeof(char) * (strlen(defline) + 1));
}


void
sim4polishBuilder::setGENdefline(char *defline) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setGENdefline()-- no polish to build; create() not called\n");
    return;
  }
  delete [] it->_genDefLine;
  it->_genDefLine = new char [strlen(defline) + 1];
  memcpy(it->_genDefLine, defline, sizeof(char) * (strlen(defline) + 1));
}


void
sim4polishBuilder::setNumberOfMatches(u32bit nummatches, u32bit nummatchesN) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setNumberOfMatches()-- no polish to build; create() not called\n");
    return;
  }
  it->_numMatches  = nummatches;
  it->_numMatchesN = nummatchesN;
}


void
sim4polishBuilder::setPercentIdentity(u32bit id) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setPercentIdentitysetPercentIdentity()-- no polish to build; create() not called\n");
    return;
  }
  it->_percentIdentity  = id;
}


void
sim4polishBuilder::setMatchOrientation(char o) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setMatchOrientation()-- no polish to build; create() not called\n");
    return;
  }
  switch (o) {
    case SIM4_MATCH_ERROR:
    case SIM4_MATCH_FORWARD:
    case SIM4_MATCH_COMPLEMENT:
      it->_matchOrientation = o;
      break;
    default:
      fprintf(stderr, "sim4polishBuilder::setMatchOrientation()-- invalid match orientation\n");
      break;
  }
}

void
sim4polishBuilder::setStrandOrientation(char o) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setStrandOrientation()-- no polish to build; create() not called\n");
    return;
  }
  switch (o) {
    case SIM4_STRAND_ERROR:
    case SIM4_STRAND_POSITIVE:
    case SIM4_STRAND_NEGATIVE:
    case SIM4_STRAND_UNKNOWN:
    case SIM4_STRAND_INTRACTABLE:
    case SIM4_STRAND_FAILED:
      it->_strandOrientation = o;
      break;
    default:
      fprintf(stderr, "sim4polishBuilder::setStrandOrientation()-- invalid match orientation\n");
      break;
  }
}


void
sim4polishBuilder::addExon(u32bit estlo, u32bit esthi,
                           u32bit genlo, u32bit genhi,
                           u32bit nummatches, u32bit nummatchesN, u32bit percentid,
                           char intronorientation) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::addExon()-- no polish to build; create() not called\n");
    return;
  }

  //  If we need more space for exons, reallocate the list of pointers
  //
  if (exPos >= exMax) {
    exMax *= 2;
    sim4polishExon **t = new sim4polishExon* [exMax];
    memcpy(t, ex, exPos * sizeof(sim4polishExon *));
    delete [] ex;
    ex = t;
    for (u32bit i=exPos; i<exMax; i++)
      ex[i] = 0L;
  }

  if (ex[exPos] == 0L) {
    ex[exPos] = new sim4polishExon;
  } else {
    //  Just in case someone didn't clean up after themselves.
    delete [] ex[exPos]->_estAlignment;
    delete [] ex[exPos]->_genAlignment;
  }

  ex[exPos]->_estAlignment      = 0L;
  ex[exPos]->_genAlignment      = 0L;

  ex[exPos]->_estFrom           = estlo;
  ex[exPos]->_estTo             = esthi;
  ex[exPos]->_genFrom           = genlo + it->_genRegionOffset;
  ex[exPos]->_genTo             = genhi + it->_genRegionOffset;
  ex[exPos]->_numMatches        = nummatches;
  ex[exPos]->_numMatchesN       = nummatchesN;
  ex[exPos]->_percentIdentity   = percentid;
  ex[exPos]->_intronOrientation = intronorientation;

  ex[exPos]->_estAlignment      = 0L;
  ex[exPos]->_genAlignment      = 0L;

  exPos++;
}


void
sim4polishBuilder::addExonAlignment(char *estalign,
                                    char *genalign) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::addExonAlignment()-- no polish to build; create() not called\n");
    return;
  }

  if (exAli >= exPos) {
    fprintf(stderr, "sim4polishBuilder::addExonAlignment()-- tried to add alignment for exon %u which doesn't exist\n", exAli);
    exit(1);
  }

  ex[exAli]->_estAlignment = (char *)memdup(estalign, (strlen(estalign) + 1) * sizeof(char));
  ex[exAli]->_genAlignment = (char *)memdup(genalign, (strlen(genalign) + 1) * sizeof(char));

  exAli++;
}

sim4polish*
sim4polishBuilder::release(void) {
  sim4polish   *retval = it;

  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::release()-- no polish to build; create() not called\n");
    return(0L);
  }

  if (exPos == 0)
    return(0L);

  it->_numCovered = 0;
  it->_numExons   = exPos;
  it->_exons      = new sim4polishExon [exPos];

  for (u32bit i=0; i<exPos; i++) {
    memcpy(it->_exons + i, ex[i], sizeof(sim4polishExon));
    ex[i]->_estAlignment = 0L;  //  Owned by 'it' now
    ex[i]->_genAlignment = 0L;

    it->_numCovered += (ex[i]->_estTo - ex[i]->_estFrom + 1);
  }

  //  Last, compute the querySeqIdentity using other fields (like our
  //  just updated numCovered).
  //
  it->_querySeqIdentity = it->s4p_percentCoverageApprox();

  it    = 0L;

  exPos = 0;
  exAli = 0;

  return(retval);
}
