
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

#ifndef AS_MER_GKPSTORECHAIN
#define AS_MER_GKPSTORECHAIN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AS_PER_gkpStore.H"

#include "seqFile.H"
#include "seqFactory.H"

//  Transform a wgs-assembler gkpStore object into a kmer seqFile object, by placing all sequence
//  into one virtual sequence, with fragments separated by a single '.'.
//
//  Position information (what read, where) is lost.


class gkpStoreChain : public seqFile {
protected:
  gkpStoreChain(char const *filename, uint32 clr);
  gkpStoreChain();

public:
  ~gkpStoreChain();

  static void         registerFile(void) {
    seqFactory::instance()->registerFile(new gkpStoreChain());
  };

protected:
  seqFile            *openFile(const char *filename);

public:
  u32bit              find(const char *sequencename) {
    fprintf(stderr, "gkpStoreChain::find()--  Not implemented.\n");
    exit(1);
    return(~u32bitZERO);
  };

  u32bit              getSequenceLength(u32bit iid) {
    return(_chainLen[iid]);
  };

  bool                getSequence(u32bit iid,
                                  char *&h, u32bit &hLen, u32bit &hMax,
                                  char *&s, u32bit &sLen, u32bit &sMax);
  bool                getSequence(u32bit iid,
                                  u32bit bgn, u32bit end, char *s);

private:
  void                clear(void);

  gkStore            *_gkp;
  gkFragment          _frg;

  u32bit              _maxChains;
  u32bit              _maxChainLen;

  u32bit              _chainIID;
  u32bit             *_chainLen;
  u32bit             *_chainBgnFrg;

  u32bit             *_frgLengths;

  u32bit              _lastFrg;
  u32bit              _lastPos;
};


#endif  //  AS_MER_GKPSTORECHAIN
