
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "meryl.H"



merylInput::merylInput(merylOperation *o) {
  _operation   = o;
  _stream      = NULL;
  _sequence    = NULL;
#ifdef CANU
  _store       = NULL;
#endif

  _isMultiSet  = false;  //  set in initialize().

  _value       = 0;
  _valid       = false;

#ifdef CANU
  _sqBgn       = 0;
  _sqEnd       = 0;

  _read        = NULL;
  _readData    = NULL;
  _readID      = 0;
  _readPos     = UINT32_MAX;
#endif

  memset(_name, 0, FILENAME_MAX+1);
  strncpy(_name, toString(_operation->getOperation()), FILENAME_MAX);
}



merylInput::merylInput(const char *n, kmerCountFileReader *s, uint32 threadFile) {
  _operation   = NULL;
  _stream      = s;
  _sequence    = NULL;
#ifdef CANU
  _store       = NULL;
#endif

  _isMultiSet  = false;  //  set in initialize().

  if (threadFile != UINT32_MAX)
    _stream->enableThreads(threadFile);

  _value       = 0;
  _valid       = false;

#ifdef CANU
  _sqBgn       = 0;
  _sqEnd       = 0;

  _read        = NULL;
  _readData    = NULL;
  _readID      = 0;
  _readPos     = UINT32_MAX;
#endif

  memset(_name, 0, FILENAME_MAX+1);
  strncpy(_name, n, FILENAME_MAX);
}



merylInput::merylInput(const char *n, dnaSeqFile *f) {
  _operation   = NULL;
  _stream      = NULL;
  _sequence    = f;
#ifdef CANU
  _store       = NULL;
#endif

  _isMultiSet  = false;

  _value       = 0;
  _valid       = true;    //  Trick nextMer into doing something without a valid mer.

#ifdef CANU
  _sqBgn       = 0;
  _sqEnd       = 0;

  _read        = NULL;
  _readData    = NULL;
  _readID      = 0;
  _readPos     = UINT32_MAX;
#endif

  memset(_name, 0, FILENAME_MAX+1);
  strncpy(_name, n, FILENAME_MAX);
}



#ifdef CANU
merylInput::merylInput(const char *n, sqStore *s, uint32 segment, uint32 segmentMax) {
  _operation   = NULL;
  _stream      = NULL;
  _sequence    = NULL;
  _store       = s;

  _isMultiSet  = false;

  _value       = 0;
  _valid       = true;    //  Trick nextMer into doing something without a valid mer.

  _sqBgn       = 1;                                   //  C-style, not the usual
  _sqEnd       = _store->sqStore_getNumReads() + 1;   //  sqStore semantics!

  if (segmentMax > 1) {
    uint64  nBases = 0;

    for (uint32 ss=1; ss <= _store->sqStore_getNumReads(); ss++)
      nBases += _store->sqStore_getRead(ss)->sqRead_sequenceLength();

    uint64  nBasesPerSeg = nBases / segmentMax;

    _sqBgn = 0;
    _sqEnd = 0;

    nBases = 0;

    for (uint32 ss=1; ss <= _store->sqStore_getNumReads(); ss++) {
      nBases += _store->sqStore_getRead(ss)->sqRead_sequenceLength();

      if ((_sqBgn == 0) && ((nBases / nBasesPerSeg) == segment - 1))
        _sqBgn = ss;

      if ((_sqEnd == 0) && ((nBases / nBasesPerSeg) == segment))
        _sqEnd = ss;
    }

    if (segment == segmentMax)                      //  Annoying special case; if the last segment,
      _sqEnd = _store->sqStore_getNumReads() + 1;   //  sqEnd is set to the last read, not N+1.

    fprintf(stderr, "merylInput-- segment %u/%u picked reads %u-%u out of %u\n",
            segment, segmentMax, _sqBgn, _sqEnd, _store->sqStore_getNumReads());
  }

  _read        = NULL;
  _readData    = new sqReadData;
  _readID      = _sqBgn - 1;       //  Incremented before loading the first read
  _readPos     = 0;

  memset(_name, 0, FILENAME_MAX+1);
  strncpy(_name, n, FILENAME_MAX);
}
#endif



merylInput::~merylInput() {

  delete _stream;
  delete _operation;
  delete _sequence;

#ifdef CANU
  delete _readData;

  _store->sqStore_close();
#endif
}



void
merylInput::initialize(void) {
  if (_operation) {
    _operation->initialize();
    _isMultiSet = _operation->isMultiSet();
  }

  if (_stream) {
    _isMultiSet = _stream->isMultiSet();
  }
}



void
merylInput::nextMer(void) {
  char kmerString[256];

  if (_stream) {
    //fprintf(stderr, "merylIn::nextMer(%s)-- (stream)\n", _name);

    _valid = _stream->nextMer();
    _kmer  = _stream->theFMer();
    _value = _stream->theValue();
  }

  if (_operation) {
    //fprintf(stderr, "merylIn::nextMer(%s)-- (operation)\n", _name);

    _valid = _operation->nextMer();
    _kmer  = _operation->theFMer();
    _value = _operation->theValue();
  }

  //fprintf(stderr, "merylIn::nextMer(%s)-- now have valid=" F_U32 " kmer %s count " F_U64 "\n",
  //        _name, _valid, _kmer.toString(kmerString), _value);
  //fprintf(stderr, "\n");
}



bool
merylInput::loadBases(char    *seq,
                      uint64   maxLength,
                      uint64  &seqLength,
                      bool    &endOfSequence) {

  if (_stream) {
    return(false);
  }

  if (_operation) {
    return(false);
  }

  if (_sequence) {
    return(_sequence->loadBases(seq, maxLength, seqLength, endOfSequence));
  }

#ifdef CANU
  if (_store) {

    //  If no read currently loaded, load one, or return that we're done.
    //  We need to loop so we can ignore the length zero reads in seqStore
    //  that exist after correction/trimming.

    while ((_read    == NULL) ||
        (_readPos >= _read->sqRead_sequenceLength())) {
      _readID++;

      if (_readID >= _sqEnd)  //  C-style iteration, not usual sqStore semantics.
        return(false);

      _read    = _store->sqStore_getRead(_readID);
      _readPos = 0;

      _store->sqStore_loadReadData(_read, _readData);
    }

    //  How much of the read is left to return?

    uint32  len = _read->sqRead_sequenceLength() - _readPos;

    assert(len > 0);

    //  If the output space is big enough to hold the rest of the read, copy it,
    //  flagging it as the end of a sequence, and setup to load the next read.

    if (len < maxLength) {
      memcpy(seq, _readData->sqReadData_getSequence() + _readPos, sizeof(char) * len);

      _read          = NULL;

      seqLength      = len;
      endOfSequence  = true;
    }

    //  Otherwise, only part of the data will fit in the output space.

    else {
      memcpy(seq, _readData->sqReadData_getSequence() + _readPos, sizeof(char) * maxLength);

      _readPos      += maxLength;

      seqLength      = maxLength;
      endOfSequence  = false;
    }

    return(true);
  }
#endif

  return(false);
}
