
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

#include "AS_PER_gkpStore.H"

#include "seqFile.H"
#include "seqFactory.H"

//  Transform a wgs-assembler gkpStore object into a kmer seqFile
//  object.

//  About IIDs:
//    gkpStore IIDs start at 1
//    kmer IIDs start at 0
//
//  This interface adds in a zero-length sequence to be kmer IID 0,
//  and does not change the IID mapping.

//  The "name" in openFile() should refer to a gatekeeper store
//  directory.  It can optionally include ":bgn-end:clr" to return a
//  subset of the reads, and to return only the clear range.  For
//  example:
//
//    dros.gkpStore:500-1000:obtini
//
//  would return fragments 500 through 1000 inclusive, and return the
//  initial OBT clear range.

class gkpStoreFile : public seqFile {
protected:
  gkpStoreFile(char const *filename, uint32 bgn, uint32 end, uint32 clr);
  gkpStoreFile();

public:
  ~gkpStoreFile();

  static void         registerFile(void) {
    seqFactory::instance()->registerFile(new gkpStoreFile());
  };

protected:
  seqFile            *openFile(const char *filename);

public:
  uint32              find(const char *sequencename) {
    fprintf(stderr, "gkpStoreFile::find()--  Not implemented.\n");
    exit(1);
    return(~uint32ZERO);
  };

  uint32              getSequenceLength(uint32 iid) {
    return(_clrEnd[iid] - _clrBeg[iid]);
  };

  bool                getSequence(uint32 iid,
                                  char *&h, uint32 &hLen, uint32 &hMax,
                                  char *&s, uint32 &sLen, uint32 &sMax);
  bool                getSequence(uint32 iid,
                                  uint32 bgn, uint32 end, char *s);

  uint16              clrBeg(uint32 iid) { return(_clrBeg[iid]); };
  uint16              clrEnd(uint32 iid) { return(_clrEnd[iid]); };

private:
  void                clear(void);

  gkStore            *_gkp;
  gkFragment          _frg;
  uint32              _bgn;
  uint32              _end;

  uint16             *_clrBeg;
  uint16             *_clrEnd;
};


#endif  //  AS_MER_GKPSTORE_TO_FASTABASE
