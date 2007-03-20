
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
//  TODO:
//
//  Deleted fragments should be returned as zero-length sequence.
//
//  Take a clear range label as a parameter, return only that sequence.
//
//  This is probably from not ignoring deleted reads somewhere:
//    chainedSequence::sequenceNumberOfPosition()-- WARNING! Position p=44307369 too big; only 44307361 positions.
//



//  openIndex() should cache the sequence lengths (and header lengths?)
//
//  About IIDs: gkpStore starts IIDs at 1, while kmer starts at 0.
//  The _iid in gkpStoreSequence is a gkpStore IID, but externally to
//  gkpStoreSequence, we need this to be a kmer IID.

class gkpStoreSequence : public FastABase {
public:
  gkpStoreSequence(char const *gkpName);
  ~gkpStoreSequence();

  FastASequenceInCore *getSequence(void);
  FastASequenceOnDisk *getSequenceOnDisk(void);

  u32bit       sequenceLength(IID_t iid);
  u32bit       headerLength(IID_t iid);
  u32bit       currentIID(void) { return(_iid-1); };

  bool         eof(void) { return(_eof); };

  bool         find(IID_t  iid);
  bool         find(char  *id);


  bool         isIndexValid(u32bit indextypetoload=FASTA_INDEX_ANY,
                            bool   beVerbose=false) {
    return(true);
  };

  void         openIndex(u32bit indextypetoload=FASTA_INDEX_ANY);

  void         optimizeRandomAccess(void) {
    //  NOP.
  };

  u32bit       getNumberOfSequences(void) {
    return(getLastElemFragStore(_gkp) + 1);
  };

  char const  *getSourceName(void) {
    return(_gkp->storePath);
  };

  void         printDescription(FILE *out, char *name) {
    fprintf(stderr, "Not implemented.\n");
    assert(0);
  };

  bool         isRandomAccess(void) { return(true); };
  bool         isSqueezed(void)     { return(true); };

private:

  GateKeeperStore  *_gkp;
  FragStream       *_stm;

  fragRecord       *_frg;

  uint32            _iid;
  bool              _eof;

  uint32           *_seqLengths;
};


#endif  //  AS_MER_GKPSTORE_TO_FASTABASE
