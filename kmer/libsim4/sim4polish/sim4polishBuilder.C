#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "bio++.H"
#include "sim4polish.h"
#include "sim4polishBuilder.H"


sim4polishBuilder::sim4polishBuilder() {
  it = 0L;

  exPos = 0;
  exMax = 32;
  exAli = 0;
  ex    = new sim4polishExon* [exMax];

  for (u32bit i=0; i<exMax; i++)
    ex[i] = 0L;
}

sim4polishBuilder::~sim4polishBuilder() {
  s4p_destroyPolish(it);

  for (u32bit i=0; i<exMax; i++) {
    if (ex[i]) {
      free(ex[i]->estAlignment);
      free(ex[i]->genAlignment);
    }
    delete ex[i];
  }

  delete [] ex;
}


void
sim4polishBuilder::create(u32bit estid, u32bit estlen,
                          u32bit genid, u32bit genlo, u32bit genhi) {

  //  Someone didn't call release()!!
  //
  if (it) {
    fprintf(stderr, "sim4polishBuilder::create()-- WARNING:  release() not called.  Polish killed.\n");
    s4p_destroyPolish(it);
  }

  errno = 0;
  it = (sim4polish *)malloc(sizeof(sim4polish));
  if (errno) {
    fprintf(stderr, "malloc() error in sim4polishBuilder::create -- can't allocate sim4polish\n%s\n", strerror(errno));
    abort();
    exit(1);
  }

  it->estID    = estid;
  it->estLen   = estlen;
  it->estPolyA = 0;
  it->estPolyT = 0;

  it->genID = genid;
  it->genLo = genlo;
  it->genHi = genhi;

  it->numMatches        = 0;
  it->numMatchesN       = 0;
  it->percentIdentity   = 0;
  it->querySeqIdentity  = 0;
  it->matchOrientation  = SIM4_MATCH_ERROR;
  it->strandOrientation = SIM4_STRAND_ERROR;

  it->comment    = 0L;
  it->estDefLine = 0L;
  it->genDefLine = 0L;

  it->numExons = 0;
  it->exons    = 0L;
}

void
sim4polishBuilder::setPolyTails(u32bit pa, u32bit pt) {
  it->estPolyA = pa;
  it->estPolyT = pt;
}

void
sim4polishBuilder::setESTdefline(char *defline) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setESTdefline()-- no polish to build; create() not called\n");
    return;
  }
  free(it->estDefLine);
  it->estDefLine = (char *)memdup(defline, (strlen(defline) + 1) * sizeof(char));
}


void
sim4polishBuilder::setGENdefline(char *defline) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setGENdefline()-- no polish to build; create() not called\n");
    return;
  }
  free(it->genDefLine);
  it->genDefLine = (char *)memdup(defline, (strlen(defline) + 1) * sizeof(char));
}


void
sim4polishBuilder::setNumberOfMatches(u32bit nummatches, u32bit nummatchesN) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setNumberOfMatches()-- no polish to build; create() not called\n");
    return;
  }
  it->numMatches  = nummatches;
  it->numMatchesN = nummatchesN;
}


void
sim4polishBuilder::setPercentIdentity(u32bit id) {
  if (it == 0L) {
    fprintf(stderr, "sim4polishBuilder::setPercentIdentitysetPercentIdentity()-- no polish to build; create() not called\n");
    return;
  }
  it->percentIdentity  = id;
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
      it->matchOrientation = o;
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
      it->strandOrientation = o;
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
    free(ex[exPos]->estAlignment);
    free(ex[exPos]->genAlignment);
  }

  ex[exPos]->estAlignment      = 0L;
  ex[exPos]->genAlignment      = 0L;

  ex[exPos]->estFrom           = estlo;
  ex[exPos]->estTo             = esthi;
  ex[exPos]->genFrom           = genlo;
  ex[exPos]->genTo             = genhi;
  ex[exPos]->numMatches        = nummatches;
  ex[exPos]->numMatchesN       = nummatchesN;
  ex[exPos]->percentIdentity   = percentid;
  ex[exPos]->intronOrientation = intronorientation;

  ex[exPos]->estAlignment      = 0L;
  ex[exPos]->genAlignment      = 0L;

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

  ex[exAli]->estAlignment = (char *)memdup(estalign, (strlen(estalign) + 1) * sizeof(char));
  ex[exAli]->genAlignment = (char *)memdup(genalign, (strlen(genalign) + 1) * sizeof(char));

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

  it->numExons = exPos;
  it->exons    = (sim4polishExon *)malloc(exPos * sizeof(sim4polishExon));

  for (u32bit i=0; i<exPos; i++) {
    memcpy(it->exons + i, ex[i], sizeof(sim4polishExon));
    ex[i]->estAlignment = 0L;  //  Owned by 'it' now
    ex[i]->genAlignment = 0L;
  }

  //  Last, compute the querySeqIdentity using other fields.
  //
  it->querySeqIdentity = it->querySeqIdentity = (int)floor(100 * (double)(it->numMatches) / (double)(it->estLen - it->estPolyA - it->estPolyT));

  it    = 0L;

  exPos = 0;
  exAli = 0;

  return(retval);
}
