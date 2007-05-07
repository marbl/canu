
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

#ifndef AS_MER_GKPSTORE_TO_FASTABASE
#define AS_MER_GKPSTORE_TO_FASTABASE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include "AS_PER_gkpStore.h"
}

#include "bio++.H"

//  Transform a wgs-assembler gkpStore object into a kmer FastABase
//  object.
//
//  About IIDs:
//    gkpStore IIDs start at 1
//    kmer IIDs start at 0
//
//  This interface adds in a zero-length sequence to be kmer IID 0,
//  and does not change the IID mapping.
//
//  The "name" in openFile() should refer to a gatekeeper store
//  directory.  It can optionally include ":bgn-end:clr" to return a
//  subset of the reads, and to return only the clear range.  For
//  example:
//
//    dros.gkpStore:500-1000:obtini
//
//  would return fragments 500 through 1000 inclusive, and return the
//  initial OBT clear range.
//
//  NOTE!  The IID range violates a lot of assumptions!  Use at risk!
//  (assumptions like: sequences start at iid=0, go to iid=M-1, and
//  there are exactly M sequences in the file)

class gkpStoreSequence : public seqFile {
protected:
  gkpStoreSequence(char const *filename, uint32 bgn, uint32 end, uint32 clr);

public:
  gkpStoreSequence();
  ~gkpStoreSequence();

  seqFile     *openFile(const char *name);

  void         openIndex(u32bit indextypetoload=FASTA_INDEX_ANY) {};

  bool         rewind(void)   { _iid = _bgn; };
  bool         eof(void)      { return(_eof); };

  seqInCore   *getSequenceInCore(void);
  seqOnDisk   *getSequenceOnDisk(void);

  u32bit       currentIID(void) { return(_iid); };

  bool         find(seqIID  iid);
  bool         find(char   *uid);

  //  This must return the length of the sequence returned by
  //  getSequenceInCore().
  u32bit       sequenceLength(seqIID iid) {
    return(_clrEnd[iid] - _clrBeg[iid]);
  };

  //  But overmerry needs to know the untrimmed length so it can
  //  reverse-complement positions.
  u32bit       sequenceLengthUntrimmed(seqIID iid) {
    return(_seqLen[iid]);
  };

  u32bit       nameLength(seqIID iid) {
    return(64);
  };

  u32bit       getNumberOfSequences(void) {
    return(getLastElemFragStore(_gkp) + 1);
  }

  char const  *getSourceName(void) {
    return(_gkp->storePath);
  };

  void         printDescription(FILE *out, char *name) {
    fprintf(stderr, "Not implemented.\n");
    assert(0);
  };

  bool         isIndexed(void)  { return(true); };
  bool         isSqueezed(void) { return(true); };

public:
  uint16       clrBeg(seqIID iid) { return(_clrBeg[iid]); };
  uint16       clrEnd(seqIID iid) { return(_clrEnd[iid]); };

private:
  bool         getSequence(uint32 &hLen, char *&h,
                           uint32 &sLen, char *&s);

  GateKeeperStore  *_gkp;
  fragRecord       *_frg;
  uint32            _bgn;
  uint32            _end;
  uint32            _clr;
  uint32            _iid;
  bool              _eof;

  uint16           *_seqLen;
  uint16           *_clrBeg;
  uint16           *_clrEnd;
};


#endif  //  AS_MER_GKPSTORE_TO_FASTABASE
