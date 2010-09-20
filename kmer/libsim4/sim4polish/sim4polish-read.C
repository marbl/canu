#include "sim4polish.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>


void
sim4polish::s4p_readPolishS4DB(readBuffer *rb) {

  //  Clear this polish.

  _numExons = 0;

  delete [] _comment;     _comment    = 0L;
  delete [] _estDefLine;  _estDefLine = 0L;
  delete [] _genDefLine;  _genDefLine = 0L;
  delete [] _exons;       _exons      = 0L;

  //  Decide the type of record we're reading.

  //  Read it.

  u64bit    startPosition = rb->tell();

  u64bit    thisLineMax  = 1048576;
  u64bit    thisLineLen  = 0;
  char     *thisLine     = new char [thisLineMax];

  u32bit    numLines = 10240;
  u32bit    curLine  = 0;

  char    **lines   = new char * [numLines + 1];
  u32bit   *lengths = new u32bit [numLines + 1];

  memset(lines,   0, sizeof(char *) * numLines);
  memset(lengths, 0, sizeof(u32bit) * numLines);

  thisLineLen = rb->read(thisLine, thisLineMax, '\n');
  chompL(thisLine, thisLineLen);

  while (!rb->eof() && strcmp(thisLine, "sim4begin")) {
    fprintf(stderr, "sim4reader: Got '%s', expecting 'sim4begin' at byte "u64bitFMT"\n",
            thisLine, startPosition);
    thisLineLen = rb->read(thisLine, thisLineMax, '\n');
    chompL(thisLine, thisLineLen);
  }

  //  Stash the 'sim4begin' line into the lines array.
  lines[curLine]   = new char [thisLineLen + 1];
  lengths[curLine] = thisLineLen;
  memcpy(lines[curLine++], thisLine, sizeof(char) * (thisLineLen + 1));

  //  Until we hit 'sim4end' stash lines into lines.  Yes, we test the previous line, then read the
  //  next.  At the end of the loop, we'll read 'sim4end', stash it in lines[], then test.

  while (!rb->eof() && strcmp(thisLine, "sim4end")) {
    thisLineLen = rb->read(thisLine, thisLineMax, '\n');
    chompL(thisLine, thisLineLen);

    if (curLine >= numLines) {
#warning LAZY PROGRAMMER did not extend an array
      fprintf(stderr, "ERROR: too many lines, lazy programmer.\n");
      exit(1);
    }

    //  Stash the line in the lines array.
    lines[curLine]   = new char [thisLineLen + 1];
    lengths[curLine] = thisLineLen;
    memcpy(lines[curLine++], thisLine, sizeof(char) * (thisLineLen + 1));
  }

  delete [] thisLine;

  if (numLines > 0)
    s4p_linesToPolishS4DB(startPosition, numLines, lines, lengths);

  for (u32bit i=0; i<curLine; i++)
    delete [] lines[i];

  delete [] lines;
  delete [] lengths;
}



void
sim4polish::s4p_readPolishGFF3(readBuffer *rb) {
}



void
sim4polish::s4p_readPolishATAC(readBuffer *rb) {
}
