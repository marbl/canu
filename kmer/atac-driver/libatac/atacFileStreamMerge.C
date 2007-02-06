// This file is part of A2Amapper.
// Copyright (c) 2007 J. Craig Venter Institute
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


class afsm {
public:
  atacFileStream  *_theFile;
  atacMatch       *_theMatch;
  atacFeature     *_theFeature;

  bool             _endOfFile;

  afsm() {
    _theFile    = 0L;
    _theMatch   = 0L;
    _theFeature = 0L;
    _endOfFile  = false;
  };

  ~afsm() {
  };
};



atacFileStreamMerge::atacFileStreamMerge(void) {
  _filesLen      = 0;
  _filesMax      = 4;
  _files         = new afsm [_filesMax];
  _theMatchIID   = 0;
  _theFeatureIID = 0;
}

atacFileStreamMerge::~atacFileStreamMerge(void) {
  for (u32bit i=0; i<_filesLen; i++)
    delete _files[i]._theFile;
  delete [] _files;

  //  But wait!  Unless we munge our copies of data, we'll free things
  //  again when atacFileBase destructs!
  //
  _seqA = 0L;
  _seqB = 0L;
}



void
atacFileStreamMerge::writeHeader(FILE *out) {
  if (_files[0]._theFile != 0L)
    _files[0]._theFile->writeHeader(out);
}


void
atacFileStreamMerge::addFile(char const *filename) {

  if (filename == 0L)
    return;

  if (_filesLen >= _filesMax) {
    _filesMax *= 2;
    afsm *F = new afsm [_filesMax];
    memcpy(F, _files, sizeof(afsm) * _filesLen);
    delete [] _files;
    _files = F;
  }

  _files[_filesLen]._theFile = new atacFileStream(filename);

  //  Duplicate a bunch of stuff to our file.
  //
  if (_filesLen == 0) {
    strcpy(_fileA, _files[_filesLen]._theFile->assemblyFileA());
    strcpy(_fileB, _files[_filesLen]._theFile->assemblyFileB());

    strcpy(_labelA, _files[_filesLen]._theFile->labelA());
    strcpy(_labelB, _files[_filesLen]._theFile->labelB());

    //_params = _files[_filesLen]._theFile->_params;

    _seqA = _files[_filesLen]._theFile->fastaA();
    _seqB = _files[_filesLen]._theFile->fastaB();
  }

  _filesLen++;
}



atacMatch*
atacFileStreamMerge::nextMatch(char type) {
  atacMatch  *theMatch;
  u32bit      theMatchIdx;

  //  Make sure everyone has a match
  //
  for (u32bit i=0; i<_filesLen; i++) {
    if (_files[i]._endOfFile == false) {
      if (_files[i]._theMatch == 0L)
        _files[i]._theMatch = _files[i]._theFile->nextMatch(type);
      if (_files[i]._theMatch == 0L)
        _files[i]._endOfFile = true;
    }
  }

  //  Pick the smallest.
  //
  theMatch    = _files[0]._theMatch;
  theMatchIdx = 0;


  //  need to set matchIID
  //  should probably also make a new match UID, or better, fix seatac to make UIDs


  for (u32bit i=1; i<_filesLen; i++) {
    if (_files[i]._theMatch) {
      if (theMatch == 0L) {
        theMatch    = _files[i]._theMatch;
        theMatchIdx = i;
      }

      if (theMatch) {
        if (_files[i]._theMatch->iid1 <  theMatch->iid1) {
          theMatch    = _files[i]._theMatch;
          theMatchIdx = i;
        }

        if ((_files[i]._theMatch->iid1 <= theMatch->iid1) &&
            (_files[i]._theMatch->iid2 <= theMatch->iid2)) {
          theMatch    = _files[i]._theMatch;
          theMatchIdx = i;
        }
      }
    }
  }

  //  Mark it as used
  //
  _files[theMatchIdx]._theMatch = 0L;

  return(theMatch);
}


atacFeature*
atacFileStreamMerge::nextFeature(char type[4]) {
  return(0L);
}
