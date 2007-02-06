// This file is part of A2Amapper.
// Copyright (c) 2005 J. Craig Venter Institute
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

#ifndef ATAC_COMMON_H
#define ATAC_COMMON_H

#include <string>
#include <map>
using namespace std;

#include "bio++.H"

#include "atacMatch.H"
#include "atacMatchList.H"
#include "atacMatchOrder.H"

#include "atacFeature.H"
#include "atacFeatureList.H"


class atacFileBase {
public:
  atacFileBase();
  virtual ~atacFileBase();

  char           *labelA(void) { return(_labelA); };
  char           *labelB(void) { return(_labelB); };

  char           *assemblyFileA(void) { return(_fileA); };
  char           *assemblyFileB(void) { return(_fileB); };

  FastAWrapper   *fastaA(void) { return(_seqA); };
  FastAWrapper   *fastaB(void) { return(_seqB); };

protected:
  void  readHeader(char *inLine, FILE *in);
public:
  void  writeHeader(FILE *out);

protected:
  char               _fileA[1024];  //  The name of our genome files
  char               _fileB[1024];

  char               _labelA[256];  //  The label of each of the sequences
  char               _labelB[256];

  map<string,string> _params;

  FastAWrapper      *_seqA;
  FastAWrapper      *_seqB;
};


class atacFileStream : public atacFileBase {
public:
  atacFileStream(char const *filename);
  ~atacFileStream();

  atacMatch   *nextMatch(char type);
  atacFeature *nextFeature(char type[4]);

private:
  FILE        *_inFile;
  char         _inLine[1024];

  u32bit       _theMatchIID;
  u32bit       _theFeatureIID;

  atacMatch    _theMatch;
  atacFeature  _theFeature;
};


class afsm;
class atacFileStreamMerge : public atacFileBase {
public:
  atacFileStreamMerge(void);
  ~atacFileStreamMerge();

  void         addFile(char const *filename);

  atacMatch   *nextMatch(char type);
  atacFeature *nextFeature(char type[4]);

  void         writeHeader(FILE *out);

private:
  u32bit       _filesLen;
  u32bit       _filesMax;
  afsm        *_files;

  u32bit       _theMatchIID;
  u32bit       _theFeatureIID;
};


class atacFile : public atacFileBase {
public:
  atacFile(char const *filename);
  ~atacFile();

  atacMatchList  *matches(void) { return(&_matches); };
  atacMatchList  *runs(void)    { return(&_runs);    };
  atacMatchList  *clumps(void)  { return(&_clumps);  };

private:
  atacMatchList _matches;
  atacMatchList _runs;
  atacMatchList _clumps;
};

#endif  //  ATAC_COMMON_H
