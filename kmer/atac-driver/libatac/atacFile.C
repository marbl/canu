// This file is part of A2Amapper.
// Copyright (c) 2005, 2006 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bio++.H"
#include "atac.H"


static
bool
isHeader(char *inLine) {
  return((inLine[0] == '!') ||
         (inLine[0] == '#') ||
         (inLine[0] == '/'));
}


atacFileStream::atacFileStream(char const *filename) {

  if (filename == 0L)
    return;

  _inFile = stdin;

  if ((filename != 0L) && (strcmp(filename, "-") != 0)) {
    errno = 0;
    _inFile = fopen(filename, "r");
    if (errno)
      fprintf(stderr, "atacFileStream::atacFileStream()-- failed to open %s: %s\n", filename, strerror(errno)), exit(1);
  }

  readHeader(_inLine, _inFile);
}

atacFileStream::~atacFileStream() {
};


atacMatch*
atacFileStream::nextMatch(char type) {

  fgets(_inLine, 1024, _inFile);

  if (feof(_inFile))
    return(0L);

  while (!feof(_inFile)) {
    if (_inLine[0] == 'M') {
      _theMatch.decode(_inLine);
      if (_theMatch.type[0] == type)
        return(&_theMatch);
    }

    fgets(_inLine, 1024, _inFile);
  }

  return(0L);
}


atacFeature*
atacFileStream::nextFeature(char type[4]) {

  fgets(_inLine, 1024, _inFile);

  if (feof(_inFile))
    return(0L);

  while (!feof(_inFile)) {
    if (_inLine[0] == 'F') {
      _theFeature.decode(_inLine);

      //  Return the feature if it is the correct type.  0 is the
      //  wildcard.
      //
      if (((_theFeature.type[0] == 0) || (type[0] == 0) || (_theFeature.type[0] == type[0])) &&
          ((_theFeature.type[1] == 0) || (type[1] == 0) || (_theFeature.type[1] == type[1])) &&
          ((_theFeature.type[2] == 0) || (type[2] == 0) || (_theFeature.type[2] == type[2])) &&
          ((_theFeature.type[3] == 0) || (type[3] == 0) || (_theFeature.type[3] == type[3])))
        return(&_theFeature);
    }

    fgets(_inLine, 1024, _inFile);
  }

  return(0L);
}




atacFile::atacFile(char const *filename) {

  if (filename == 0L)
    return;

  FILE    *inFile = stdin;
  char     inLine[1024];

  if ((filename != 0L) && (strcmp(filename, "-") != 0)) {
    errno = 0;
    inFile = fopen(filename, "r");
    if (errno)
      fprintf(stderr, "atacFile::atacFile()-- failed to load %s: %s\n", filename, strerror(errno)), exit(1);
  }

  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in fileA and fileB.
  //
  readHeader(inLine, inFile);


  while (!feof(inFile)) {
    switch(inLine[0]) {
      case 'M':
        {
          atacMatch m(inLine);

          if (m.sanity(fastaA(), fastaB(), inLine)) {
            if        (m.type[0] == 'u') {
              _matches.add(m);
            } else if (m.type[0] == 'r') {
              _runs.add(m);
            } else if (m.type[0] == 'c') {
              _clumps.add(m);
            } else {
            }
          }
        }
        break;

      case 'F':
        {
        }
        break;

      default:
        {
        }
        break;
    }

    fgets(inLine, 1024, inFile);
  }
}

atacFile::~atacFile() {
}








atacFileBase::atacFileBase() {
  _fileA[0] = 0;
  _fileB[0] = 0;

  _labelA[0] = 0;
  _labelB[0] = 0;

  _seqA = 0L;
  _seqB = 0L;
}


atacFileBase::~atacFileBase() {
  delete _seqA;
  delete _seqB;
}


void
atacFileBase::readHeader(char *inLine, FILE *in) {

  fgets(inLine, 1024, in);
  while (!feof(in) && isHeader(inLine)) {
    chomp(inLine);

    if (inLine[0] == '/') {
      char *key = inLine + 1;
      char *val = inLine + 1;

      while (isspace(*key))  key++;  //  Skip whitespace between "/" and the key
      while (*val != '=')    val++;  //  Move to the "="
      *val++ = 0;                    //  Terminate the key
      while (isspace(*val))  val++;  //  Skip whitespace between "=" and the val

      chomp(key);
      chomp(val);

      //fprintf(stderr, "key='%s' val='%s'\n", key, val);

#if 0
      string K = key;
      string V = val;
      _params[K] = V;
#endif

      //  Save ones we use

      if (strncmp(key, "assemblyFile1", 14) == 0)
        strcpy(_fileA, val);

      if (strncmp(key, "assemblyFile2", 14) == 0)
        strcpy(_fileB, val);

      if (strncmp(key, "assemblyId1", 12) == 0)
        strcpy(_labelA, val);

      if (strncmp(key, "assemblyId2", 12) == 0)
        strcpy(_labelB, val);
    }

    //  Otherwise, it's a comment or the header

    fgets(inLine, 1024, in);
  }

  //fprintf(stderr, "assemblyFile1 = '%s'\n", _fileA);
  //fprintf(stderr, "assemblyFile2 = '%s'\n", _fileB);
  //fprintf(stderr, "assemblyId1   = '%s'\n", _labelA);
  //fprintf(stderr, "assemblyId2   = '%s'\n", _labelB);

  //  Open some FastAWrappers for each of the files
  //
  if (_fileA && _fileA[0]) {
    _seqA = new FastAWrapper(_fileA);
    _seqA->openIndex();
  }
  if (_fileB && _fileB[0]) {
    _seqB = new FastAWrapper(_fileB);
    _seqB->openIndex();
  }
}
