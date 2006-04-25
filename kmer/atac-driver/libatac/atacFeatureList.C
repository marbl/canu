// This file is part of A2Amapper.
// Copyright (c) 2006 J. Craig Venter Institute
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


atacFeatureList::atacFeatureList(char *filename, char *type, bool saveLine) {

  FILE *inFile = stdin;
  if ((filename != 0L) && (strcmp(filename, "-") != 0)) {
    errno = 0;
    inFile = fopen(filename, "r");
    if (errno)
      fprintf(stderr, "atacFeatureList::atacFeatureList()-- failed to load %s: %s\n", filename, strerror(errno)), exit(1);
  }

  //  Read the preamble, look for our data sources.  This leaves us with
  //  the first match in the inLine, and fills in file1 and file2.
  //
  char    inLine[1024];
  readHeader(inLine, inFile, _file, 0L, 0L);

  _label[0] = 0;

  //  Open some FastAWrappers for each of the files -- we use these
  //  only to get the length of the sequence.
  //
  _seq = 0L;
  if (_file && _file[0]) { 
    _seq = new FastAWrapper(_file);
    _seq->openIndex();
  }

  _featuresLen = 0;
  _featuresMax = 32 * 1048576;
  _features    = new atacFeature [_featuresMax];

  while (!feof(inFile)) {
    if (inLine[0] == 'F') {
      splitToWords  S(inLine);

      //  Save the name/label
      if (_label[0] == 0)
        decodeAtacName(S[4], _label);

      u32bit  iid=0, pos=0, len=0;

      decodeFeature(S, iid, pos, len);

      bool  featureOK = true;
      if (_seq) {
        if ((pos) > _seq->sequenceLength(iid) || (pos + len) > _seq->sequenceLength(iid)) {
          chomp(inLine);
          fprintf(stderr, "Feature longer than sequence (by "u32bitFMT"bp): seqLen="u32bitFMTW(8)" %s\n",
                  pos + len - _seq->sequenceLength(iid),
                  _seq->sequenceLength(iid), inLine);
          featureOK = false;
        }

        if (iid >= _seq->getNumberOfSequences()) {
          chomp(inLine);
          fprintf(stderr, "Feature references invalid sequence iid: %s\n", inLine);
          featureOK = false;
        }
      }

      if (featureOK) {

        //  Add it to our list of matches
        //
        if (_featuresLen > _featuresMax) {
          fprintf(stderr, "SORRY!  I don't feel like reallocating matches.  Increase\n");
          fprintf(stderr, "the preallocated size in %s\n", __FILE__);
          exit(1);
        }

        strncpy(_features[_featuresLen].featureuid,  S[2], 16);
        strncpy(_features[_featuresLen].parentuid, S[3], 16);

        _features[_featuresLen].featureuid[15]  = 0;
        _features[_featuresLen].parentuid[15] = 0;

        _features[_featuresLen].featureiid = _featuresLen;

        _features[_featuresLen].type[0] = S[1][0];
        _features[_featuresLen].type[1] = S[1][1];
        if (S[1][1])
          _features[_featuresLen].type[2] = S[1][2];
        _features[_featuresLen].type[3] = 0;

        _features[_featuresLen].iid = iid;
        _features[_featuresLen].pos = pos;
        _features[_featuresLen].len = len;

        _featuresLen++;
      }
    }

    fgets(inLine, 1024, inFile);
  }
}


static
int
sort_(const void *a, const void *b) {
  const atacFeature *A = (const atacFeature *)a;
  const atacFeature *B = (const atacFeature *)b;

  if (A->iid < B->iid)  return(-1);
  if (A->iid > B->iid)  return(1);
  if (A->pos < B->pos)  return(-1);
  if (A->pos > B->pos)  return(1);
  if (A->len > B->len)  return(-1);
  if (A->len < B->len)  return(1);
  return(0);
}

static
int
sortfeatureuid_(const void *a, const void *b) {
  const atacFeature *A = (const atacFeature *)a;
  const atacFeature *B = (const atacFeature *)b;

  int  r = strcmp(A->featureuid, B->featureuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);
  r = strcmp(A->parentuid, B->parentuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);

  return(0);
}

static
int
sortparentuid_(const void *a, const void *b) {
  const atacFeature *A = (const atacFeature *)a;
  const atacFeature *B = (const atacFeature *)b;

  int  r = strcmp(A->parentuid, B->parentuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);
  r = strcmp(A->featureuid, B->featureuid);
  if (r < 0)  return(-1);
  if (r > 0)  return(1);
  
  return(0);
}


void
atacFeatureList::sort(u32bit first, u32bit len) {
  if (len == 0) len = _featuresLen;
  qsort(_features + first, len, sizeof(atacFeature), sort_);
}

void
atacFeatureList::sortFeatureUID(u32bit first, u32bit len) {
  if (len == 0) len = _featuresLen;
  qsort(_features + first, len, sizeof(atacFeature), sortfeatureuid_);
}

void
atacFeatureList::sortParentUID(u32bit first, u32bit len) {
  if (len == 0) len = _featuresLen;
  qsort(_features + first, len, sizeof(atacFeature), sortparentuid_);
}
