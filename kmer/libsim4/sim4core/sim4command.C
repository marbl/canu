#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sim4command.H"



//  Run a single EST against a genomic range
//
//  XXX: We should pull out the EST and GEN from the FastAWrapper,
//  and store them as the "two char*" method.
//
sim4command::sim4command(u32bit        ESTid,
                         FastAWrapper *ESTs,
                         u32bit        GENid,
                         u32bit        GENlo,
                         u32bit        GENhi,
                         FastAWrapper *GENs,
                         bool          doForward,
                         bool          doReverse) {

  _estIdx = ESTid;

  _ESTs              = ESTs;
  _ESTloaded         = 0L;
  _ESTsequence       = 0L;
  _ESTsequenceLength = 0;

  _genIdx = GENid;
  _genLo  = GENlo;
  _genHi  = GENhi;

  _GENs              = GENs;
  _GENloaded         = 0L;
  _GENsequence       = 0L;
  _GENsequenceLength = 0;

  _doForward = doForward;
  _doReverse = doReverse;

  _strandIndicator = 0L;
}


sim4command::sim4command(FastASequenceInCore  *EST,
                         FastASequenceInCore  *GEN,
                         u32bit                GENlo,
                         u32bit                GENhi,
                         bool                  doForward,
                         bool                  doReverse) {

  _estIdx = EST->getIID();

  _ESTs              = 0L;
  _ESTloaded         = EST;
  _ESTsequence       = 0L;
  _ESTsequenceLength = 0;

  _genIdx = GEN->getIID();
  _genLo  = GENlo;
  _genHi  = GENhi;

  _GENs              = 0L;
  _GENloaded         = GEN;
  _GENsequence       = 0L;
  _GENsequenceLength = 0;

  _doForward = doForward;
  _doReverse = doReverse;

  _strandIndicator = 0L;
}


//  Use two char*'s for sequence sources
//
sim4command::sim4command(char             *EST,
                         u32bit            ESTlen,
                         char             *GEN,
                         u32bit            GENlen,
                         u32bit            GENlo,
                         u32bit            GENhi,
                         bool              doForward,
                         bool              doReverse) {
  _estIdx = 0;

  _ESTs              = 0L;
  _ESTloaded         = 0L;
  _ESTsequence       = EST;
  _ESTsequenceLength = ESTlen;

  _genIdx = 0;
  _genLo  = GENlo;
  _genHi  = GENhi;

  _GENs              = 0L;
  _GENloaded         = 0L;
  _GENsequence       = GEN;
  _GENsequenceLength = GENlen;

  _doForward = doForward;
  _doReverse = doReverse;

  _strandIndicator = 0L;
}


sim4command::~sim4command() {
  if (_ESTs)
    delete _ESTloaded;
  if (_GENs)
    delete _GENloaded;
  delete [] _strandIndicator;
}


//  Make absolutely sure that the genomic sequence start and end
//  positions are within the actual sequence.  Ideally, this should
//  be checked by whatever generates the input, but it probably
//  isn't.
//
//  If the end position is too big, make it the same as the sequence
//  length.
//
//  If the start position is bigger than the (corrected) end
//  position, make it 100K less than the end position.
//
//  This has the side-effect of loading the genomic sequence.
//
void
sim4command::finalize(void) {

  if (_genHi > getGENlength())
    _genHi = getGENlength();

  if (_genLo > _genHi)
    if (_genHi > 100000)
      _genLo = _genHi - 100000;
    else
      _genLo = 0;
}



//  get() routines have multple cases
//
//  if no fastawrapper, they can quickly return
//  otherwise
//  if nothing loaded or the thing loaded isn't right:
//    delete the current
//    load the correct
//

void
sim4command::loadEST(void) {
  if ((_ESTloaded == 0L) ||
      (_ESTloaded->getIID() != _estIdx)) {
    delete _ESTloaded;
    if (_ESTs->find(_estIdx) == false) {
      fprintf(stderr, "ERROR: Can't find IID %lu in the set of ESTs\n", _estIdx);
      exit(1);
    }
    _ESTloaded = _ESTs->getSequence();
  }
}


u32bit
sim4command::getESTidx(void) {
  if (_ESTsequence)
    return(0);
  return(_estIdx);
}

char*
sim4command::getESTheader(void) {
  if (_ESTsequence)
    return("anonymous cDNA sequence");
  loadEST();
  return(_ESTloaded->header());
}

char*
sim4command::getESTsequence(void) {
  if (_ESTsequence)
    return(_ESTsequence);
  loadEST();
  return(_ESTloaded->sequence());
}

u32bit
sim4command::getESTlength(void) {
  if (_ESTsequence)
    return(_ESTsequenceLength);
  loadEST();
  return(_ESTloaded->sequenceLength());
}





void
sim4command::loadGEN(void) {
  if ((_GENloaded == 0L) ||
      (_GENloaded->getIID() != _genIdx)) {
    delete _GENloaded;
    if (_GENs->find(_genIdx) == false) {
      fprintf(stderr, "ERROR: Can't find IID %lu in the set of genomic sequences\n", _genIdx);
      exit(1);
    }
    _GENloaded = _GENs->getSequence();
  }
}

char*
sim4command::getGENheader(void) {
  if (_GENsequence)
    return("anonymous genomic sequence");
  loadGEN();
  return(_GENloaded->header());
}

char*
sim4command::getGENsequence(void) {
  if (_GENsequence)
    return(_GENsequence);
  loadGEN();
  return(_GENloaded->sequence());
}

u32bit
sim4command::getGENlength(void) {
  if (_GENsequence)
    return(_GENsequenceLength);
  loadGEN();
  return(_GENloaded->sequenceLength());
}

