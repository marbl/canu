#include "sim4polishReader.H"

#include "util++.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "sim4polishWriter.H"

sim4polishReader::sim4polishReader(const char *name, sim4polishWriter *writer) {

  if (name)
    _rb = new readBuffer(name);
  else
    _rb = new readBuffer(writer->surrenderToReader());

  //  Attempt to decide on the style of the input, based on the first line.

  char          firstLine[1024];
  splitToWords  firstWords;

  _rb->read(firstLine, 1024, '\n');
  _rb->seek(0);

  //  This fixes a bug in split to words, that white space at the end isn't trimmed.
  chomp(firstLine);

  firstWords.split(firstLine);

  if        (strcmp(firstWords[0], "sim4begin") == 0) {
    _style = sim4polishS4DB;

  } else if (strcmp(firstWords[0], "##gff-version") == 0) {
    if (strcmp(firstWords[1], "3") == 0)
      _style = sim4polishGFF3;
    else
      fprintf(stderr, "sim4polishReader()-- GFF format version %s not supported; only version 3 is supported.\n",
              firstWords[1]), exit(1);

  } else if ((strcmp(firstWords[0], "!format") == 0) &&
             (strcmp(firstWords[1], "atac") == 0)) {
  if (strcmp(firstWords[2], "1.0") == 0)
      _style = sim4polishATAC;
    else
      fprintf(stderr, "sim4polishReader()-- ATAC format version %s not supported; only version 1.0 is supported.\n",
              firstWords[2]), exit(1);

  } else {
    fprintf(stderr, "sim4polishReader()-- Failed to open '%s' for reading: unknown format.\n",
            _rb->filename()), exit(1);
  }
}


sim4polishReader::~sim4polishReader() {
  delete _rb;
  _rb = 0L;
}


sim4polish *
sim4polishReader::nextAlignment(void) {
  sim4polish *p = 0L;

  if (_rb->eof())
    return(p);

  p = new sim4polish(_rb, _style);

  if (p->_numExons == 0) {
    delete p;
    p = 0L;
  }

  return(p);
}


bool
sim4polishReader::nextAlignment(sim4polish * &p) {

  delete p;

  if (_rb->eof())
    return(false);

  p = new sim4polish(_rb, _style);

  if (p->_numExons == 0) {
    delete p;
    return(false);
  }

  return(true);
}
