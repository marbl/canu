
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received (LICENSE.txt) a copy of the GNU General Public 
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
#include "AS_MER_stream.h"
#include <stdlib.h>

merStream::merStream(cds_uint32 merSize, char *fragstore) {
  _theAlloc        = 1024;
  _theSeq          = new char [_theAlloc];
  _theQlt          = new char [_theAlloc];
  _theLen          = 0;
  _thePos          = 0;
  _endPos          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;
  _theMerMask      = CDS_UINT64_MASK(_merSize << 1);

  _theFMer         = 0;
  _theRMer         = 0;

  _theRMerShift    = (_merSize << 1) - 2;

  //  Open the fragstore
  //
  _fs = openFragStore(fragstore, "r+");
  if (_fs == NULLSTOREHANDLE) {
    fprintf(stderr, "ERROR:  couldn't open the fragStore '%s'\n", fragstore);
    exit(1);
  }

  _fsBufferSize = 1024 * 1024;
  _fsBuffer     = new char [_fsBufferSize];

  _fsh = openFragStream(_fs, _fsBuffer, _fsBufferSize);
  _rs  = new_ReadStruct();

  _skipNum = 1;
  loadMer(_merSize - 1);
}


merStream::merStream(cds_uint32 merSize, char *fragstore, int skipNum) {
  _theAlloc        = 1024;
  _theSeq          = new char [_theAlloc];
  _theQlt          = new char [_theAlloc];
  _theLen          = 0;
  _thePos          = 0;
  _endPos          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;
  _theMerMask      = CDS_UINT64_MASK(_merSize << 1);

  _theFMer         = 0;
  _theRMer         = 0;

  _theRMerShift    = (_merSize << 1) - 2;

  //  Open the fragstore
  //
  _fs = openFragStore(fragstore, "r+");
  if (_fs == NULLSTOREHANDLE) {
    fprintf(stderr, "ERROR:  couldn't open the fragStore '%s'\n", fragstore);
    exit(1);
  }

  _fsBufferSize = 1024 * 1024;
  _fsBuffer     = new char [_fsBufferSize];

  _fsh = openFragStream(_fs, _fsBuffer, _fsBufferSize);
  _rs  = new_ReadStruct();

  _skipNum=skipNum;
  loadMer(_merSize - 1);

}

merStream::~merStream() {
  delete_ReadStruct(_rs);

  closeFragStream(_fsh);
  delete [] _fsBuffer;

  closeFragStore(_fs);

  delete [] _theSeq;
  delete [] _theQlt;
}

#if 0
char const *
merStream::theFMerString(void) {
  for (cds_uint32 i=0; i<_merSize; i++)
    _theMerString[_merSize-i-1] = decompressSymbol[(_theFMer >> (2*i)) & 0x03];
  _theMerString[_merSize] = 0;
  return(_theMerString);
}

char const *
merStream::theRMerString(void) {
  for (cds_uint32 i=0; i<_merSize; i++)
    _theMerString[_merSize-i-1] = decompressSymbol[(_theRMer >> (2*i)) & 0x03];
  _theMerString[_merSize] = 0;
  return(_theMerString);
}
#endif
