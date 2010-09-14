#include "sim4polish.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>


//  Utility for reading a whole line, safely, from a file.  Newlines at the end are trimmed out.
//
class _line {
public:
  _line() {
    _len = 0;
    _max = 1024;
    _str = new char [_max];
    _num = 0;
  };
  ~_line() {
    delete [] _str;
  };

  void   readLine(FILE *F) {
    _len = 0;

    if (fgets(_str, _max, F)) {
      _num++;
      _len = strlen(_str);

      while ((_len + 1 >= _max) &&
             (_str[_len - 1] != '\n')) {
        _max *= 2;
        char *nnn = new char [_max];
        memcpy(nnn, _str, sizeof(char) * _len);
        delete [] _str;
        _str = nnn;

        fgets(_str + _len, _max - _len, F);
        _len += strlen(_str + _len);
      }

      while ((_len > 0) &&
             ((_str[_len - 1] == '\r') ||
              (_str[_len - 1] == '\n')))
        _str[--_len] = 0;
    }
  };

public:
  u32bit   _len;  //  length of the line
  u32bit   _max;  //  maximum space
  char    *_str;  //  the string
  u32bit   _num;  //  line number
};



void
sim4polish::s4p_readPolish(FILE *F) {
  _line            l;

  //  Clear this polish.

  _numExons = 0;

  delete [] _comment;     _comment    = 0L;
  delete [] _estDefLine;  _estDefLine = 0L;
  delete [] _genDefLine;  _genDefLine = 0L;
  delete [] _exons;       _exons      = 0L;

  //  Decide the type of record we're reading.

  //  Read it.

  u32bit    numLines = 10240;
  u32bit    curLine  = 0;

  char    **lines   = new char * [numLines + 1];
  u32bit   *lengths = new u32bit [numLines + 1];

  memset(lines,   0, sizeof(char *) * numLines);
  memset(lengths, 0, sizeof(u32bit) * numLines);

  l.readLine(F);
  while(!feof(F) && strcmp(l._str, "sim4begin")) {
    fprintf(stderr, "sim4reader: Got '%s', expecting 'sim4begin' at line "u32bitFMT"\n",
            l._str, l._num);
    l.readLine(F);
  }

  //  Stash the 'sim4begin' line into the lines array.

  u32bit    lineNumber = l._num;

  lines[curLine]   = new char [l._len + 1];
  lengths[curLine] = l._len;

  memcpy(lines[curLine++], l._str, sizeof(char) * (l._len + 1));

  //  Until we hit 'sim4end' stash lines into lines.  Yes, we test the previous line, then read the
  //  next.  At the end of the loop, we'll read 'sim4end', stash it in lines[], then test.

  while(!feof(F) && strcmp(l._str, "sim4end")) {
    l.readLine(F);

    if (curLine >= numLines) {
#warning LAZY PROGRAMMER did not extend an array
      fprintf(stderr, "ERROR: too many lines, lazy programmer.\n");
      exit(1);
    }

    lines[curLine]   = new char [l._len + 1];
    lengths[curLine] = l._len;

    memcpy(lines[curLine++], l._str, sizeof(char) * (l._len + 1));
  }

  if (feof(F) == false)
    s4p_linesToPolish(lineNumber, numLines, lines, lengths);

  for (u32bit i=0; i<curLine; i++)
    delete [] lines[i];

  delete [] lines;
  delete [] lengths;
}
