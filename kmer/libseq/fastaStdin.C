#include "fastaStdin.H"
#include "alphabet.h"

fastaStdin::fastaStdin(const char *filename) {
  clear();

#ifdef DEBUG
  fprintf(stderr, "fastaStdin::fastaStdin()-- '%s'\n", (filename) ? filename : "NULLPOINTER");
#endif

  strcpy(_filename, "(stdin)");

  _rb = new readBuffer("-");
}



fastaStdin::fastaStdin() {
  clear();
}



fastaStdin::~fastaStdin() {
  delete _rb;
}



seqFile *
fastaStdin::openFile(const char *filename) {

#ifdef DEBUG
  fprintf(stderr, "fastaStdin::openFile()-- '%s'\n", (filename) ? filename : "NULLPOINTER");
#endif

  if (((filename == 0L) && (isatty(fileno(stdin)) == 0)) ||
      ((filename != 0L) && (filename[0] == '-') && (filename[1] == 0)))
    return(new fastaStdin(0L));

  return(0L);
}



u32bit
fastaStdin::find(const char *sequencename) {
  fprintf(stderr, "fastaStdin::find()-- WARNING!  Used for random access.\n");
  return(~u32bitZERO);
}



u32bit
fastaStdin::getSequenceLength(u32bit iid) {
  fprintf(stderr, "fastaStdin::getSequenceLength()-- WARNING!  Used for random access.\n");
  return(0);
}



bool
fastaStdin::getSequence(u32bit iid,
                             char *&h, u32bit &hLen, u32bit &hMax,
                             char *&s, u32bit &sLen, u32bit &sMax) {

#ifdef DEBUG
  fprintf(stderr, "fastaStdin::getSequence(full)-- "u32bitFMT"\n", iid);
#endif

  if (iid != _thisIID)
    fprintf(stderr, "fastaStdin::getSequence(full)-- WARNING!  Used for random access.\n"), exit(1);

  if (sMax == 0) {
    sMax = 2048;
    s    = new char [sMax];
  }

  if (hMax == 0) {
    hMax = 2048;
    h    = new char [hMax];
  }

  return(loadNextSequence(h, hLen, hMax, s, sLen, sMax));
}



bool
fastaStdin::getSequence(u32bit iid,
                             u32bit bgn, u32bit end, char *s) {

#ifdef DEBUG
  fprintf(stderr, "fastaStdin::getSequence(part)-- "u32bitFMT"\n", iid);
#endif
  assert(0);
  return(false);
}



void
fastaStdin::clear(void) {
  memset(_filename, 0, FILENAME_MAX);
  memset(_typename, 0, FILENAME_MAX);

  strcpy(_typename, "FastAstream");

  _numberOfSequences = ~u32bitZERO;

  _rb = 0L;
  _thisIID = 0;
}



bool
fastaStdin::loadNextSequence(char *&h, u32bit &hLen, u32bit &hMax,
                                  char *&s, u32bit &sLen, u32bit &sMax) {
  char   x   = _rb->read();

  //  Skip whitespace at the start of the sequence.
  while ((_rb->eof() == false) && (whitespaceSymbol[x] == true))
    x = _rb->read();

  //  We should be at a '>' character now.  Fail if not.
  if (_rb->eof() == true)
    return(false);
  if (x != '>')
    fprintf(stderr, "fastaStdin::loadNextSequence(part)-- ERROR: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
            _filename, x), exit(1);

  //  Skip the '>' in the defline
  x = _rb->read();

  //  Skip whitespace between the '>' and the defline
  while ((_rb->eof() == false) && (whitespaceSymbol[x] == true) && (x != '\r') && (x != '\n'))
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
  while ((_rb->eof() == false) && (whitespaceSymbol[x] == true))
    x = _rb->read();

  //  Copy the sequence, until EOF or the next '>'.
  while ((_rb->eof() == false) && (_rb->peek() != '>')) {
    if (whitespaceSymbol[x] == false) {
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

  _thisIID++;

  return(true);
}
