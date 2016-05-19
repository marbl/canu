
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
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "fastqStdin.H"
#include "dnaAlphabets.H"

fastqStdin::fastqStdin(const char *filename) {
  clear();

#ifdef DEBUG
  fprintf(stderr, "fastqStdin::fastqStdin()-- '%s'\n", (filename) ? filename : "NULLPOINTER");
#endif

  if (filename == 0L) {
    strcpy(_filename, "(stdin)");
    _rb = new readBuffer("-");

  } else {
    _pipe = popen(filename, "r");
    _rb = new readBuffer(_pipe);
  }
}


fastqStdin::fastqStdin() {
  clear();
}



fastqStdin::~fastqStdin() {
  delete    _rb;
  delete [] _header;
  delete [] _sequence;
}



seqFile *
fastqStdin::openFile(const char *filename) {

#ifdef DEBUG
  fprintf(stderr, "fastqStdin::openFile()-- '%s'\n", (filename) ? filename : "NULLPOINTER");
#endif

  if (((filename == 0L) && (isatty(fileno(stdin)) == 0)) ||
      ((filename != 0L) && (filename[0] == '-') && (filename[1] == 0)))
    return(new fastqStdin(0L));

  //  The stdin variants also handle compressed inputs (because we can't seek in these).

  if (filename == 0L)
    return(0L);

  uint32       fl   = strlen(filename);
  char         cmd[32 + fl];

  if      ((filename[fl-3] == '.') && (filename[fl-2] == 'g') && (filename[fl-1] == 'z'))
    sprintf(cmd, "gzip -dc %s", filename);

  else if ((filename[fl-4] == '.') && (filename[fl-3] == 'b') && (filename[fl-2] == 'z') && (filename[fl-1] == '2'))
    sprintf(cmd, "bzip2 -dc %s", filename);

  else if ((filename[fl-3] == '.') && (filename[fl-2] == 'x') && (filename[fl-1] == 'z'))
    sprintf(cmd, "xz -dc %s", filename);

  else
    return(0L);

  return(new fastqStdin(cmd));
}



uint32
fastqStdin::getNumberOfSequences(void) {
  if (_rb->peek() == 0)
    return(_nextIID);
  else
    return(_nextIID + 1);
}


uint32
fastqStdin::find(const char *sequencename) {
  fprintf(stderr, "fastqStdin::find()-- ERROR!  Used for random access on sequence '%s'.\n", sequencename);
  assert(0);
  return(~uint32ZERO);
}



uint32
fastqStdin::getSequenceLength(uint32 iid) {

  if (iid == _nextIID)
    if (loadNextSequence(_header, _headerLen, _headerMax, _sequence, _sequenceLen, _sequenceMax) == false)
      return(0);

  if (iid + 1 != _nextIID) {
    fprintf(stderr, "fastqStdin::getSequence()-- ERROR!  Used for random access.  Requested iid=%u, at iid=%u\n",
            iid, _nextIID);
    assert(0);
  }

  return(strlen(_sequence));
}



bool
fastqStdin::getSequence(uint32 iid,
                        char *&h, uint32 &hLen, uint32 &hMax,
                        char *&s, uint32 &sLen, uint32 &sMax) {
  bool  ret = true;

#ifdef DEBUG
  fprintf(stderr, "fastqStdin::getSequence(full)-- "F_U32"\n", iid);
#endif

  if (iid == _nextIID)
    if (loadNextSequence(_header, _headerLen, _headerMax, _sequence, _sequenceLen, _sequenceMax) == false)
      return(false);

  if (iid + 1 != _nextIID) {
    fprintf(stderr, "fastqStdin::getSequence(full)-- ERROR!  Used for random access.  Requested iid=%u, at iid=%u\n",
            iid, _nextIID);
    assert(0);
  }

  if (hLen < _headerMax) {
    delete [] h;
    hMax = _headerMax;
    h    = new char [hMax];
  }

  if (sLen < _sequenceMax) {
    delete [] s;
    sMax = _sequenceMax;
    s    = new char [sMax];
  }

  memcpy(h, _header, _headerLen + 1);
  hLen = _headerLen;

  memcpy(s, _sequence, _sequenceLen + 1);
  sLen = _sequenceLen;

  return(true);
}



bool
fastqStdin::getSequence(uint32 iid,
                        uint32 bgn, uint32 end, char *s) {

  fprintf(stderr, "fastqStdin::getSequence(part)-- ERROR!  Used for random access on iid "F_U32" from position "F_U32"-"F_U32".\n", iid, bgn, end);
  assert(0);
  return(false);
}



void
fastqStdin::clear(void) {
  memset(_filename, 0, FILENAME_MAX);
  memset(_typename, 0, FILENAME_MAX);

  _randomAccessSupported = false;

  strcpy(_typename, "FastQstream");

  _numberOfSequences = ~uint32ZERO;

  _rb          = 0L;
  _nextIID     = 0;
  _pipe        = 0L;

  _header      = 0L;
  _headerLen   = 0;
  _headerMax   = 0;

  _sequence    = 0L;
  _sequenceLen = 0;
  _sequenceMax = 0;
}



bool
fastqStdin::loadNextSequence(char *&h, uint32 &hLen, uint32 &hMax,
                             char *&s, uint32 &sLen, uint32 &sMax) {

  if (hMax == 0) {
    hMax = 2048;
    h    = new char [hMax];
  }

  if (sMax == 0) {
    sMax = 2048;
    s    = new char [sMax];
  }

  hLen = 0;
  sLen = 0;

  char   x   = _rb->read();

  //  Skip whitespace at the start of the sequence.
  while ((_rb->eof() == false) && (alphabet.isWhitespace(x) == true))
    x = _rb->read();

  //  We should be at a '@' character now.  Fail if not.
  if (_rb->eof() == true)
    return(false);
  if (x != '@')
    fprintf(stderr, "fastqStdin::loadNextSequence(part)-- ERROR: In %s, expected '@' at beginning of defline, got '%c' instead.\n",
            _filename, x), exit(1);

  //  Skip the '@' in the defline
  x = _rb->read();

  //  Skip whitespace between the '@' and the defline
  while ((_rb->eof() == false) && (alphabet.isWhitespace(x) == true) && (x != '\r') && (x != '\n'))
    x = _rb->read();

  //  Copy the defline, until the first newline.
  while ((_rb->eof() == false) && (x != '\r') && (x != '\n')) {
    h[hLen++] = x;
    if (hLen >= hMax) {
      //fprintf(stderr, "realloc header\n");
      hMax += 2048;
      char *H = new char [hMax];
      memcpy(H, h, hLen);
      delete [] h;
      h = H;
    }
    x = _rb->read();
  }
  h[hLen] = 0;

  //  Skip whitespace between the defline and the sequence.
  while ((_rb->eof() == false) && (alphabet.isWhitespace(x) == true))
    x = _rb->read();

  //  Copy the sequence, until EOF or the start of the QV bases.
  while ((_rb->eof() == false) && (_rb->peek() != '+')) {
    if (alphabet.isWhitespace(x) == false) {
      s[sLen++] = x;
      if (sLen >= sMax) {
        //fprintf(stderr, "realloc sequence\n");
        sMax *= 2;
        char *S = new char [sMax];
        memcpy(S, s, sLen);
        delete [] s;
        s = S;
      }
    }
    x = _rb->read();
  }
  s[sLen] = 0;

  //  Skip the rest of the QV id line and then the entire QV line.

  //x = _rb->read();
  assert((_rb->eof() == true) || (x == '+'));

  while ((_rb->eof() == false) && (x != '\r') && (x != '\n'))
    x = _rb->read();
  x = _rb->read();
  while ((_rb->eof() == false) && (x != '\r') && (x != '\n'))
    x = _rb->read();

  _nextIID++;

  return(true);
}
