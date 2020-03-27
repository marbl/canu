
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

/*
 * https://github.com/PacificBiosciences/FALCON/blob/master/src/c/falcon.c
 *
 * Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted (subject to the limitations in the
 * disclaimer below) provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above
 *  copyright notice, this list of conditions and the following
 *  disclaimer in the documentation and/or other materials provided
 *  with the distribution.
 *
 *  * Neither the name of Pacific Biosciences nor the names of its
 *  contributors may be used to endorse or promote products derived
 *  from this software without specific prior written permission.
 *
 * NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
 * GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
 * BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#include "runtime.H"

#include "tgStore.H"  //  for tgTig

#include <vector>
#include <string>
#include <cstdlib>
#include <stdint.h>

#include <vector>
#include <algorithm>

using namespace std;

#include "falconConsensus-alignTag.H"
#include "falconConsensus-msa.H"

#include "system.H"

#ifndef FALCONCONSENSUS_H
#define FALCONCONSENSUS_H

#undef BRI
#undef BRIOUT

class falconInput {
public:
  falconInput() {
    read       = NULL;
    readLength = 0;

    ident      = 0;

    placedBgn  = 0;
    placedEnd  = 0;

    alignedBgn = 0;
    alignedEnd = 0;
  };

  void addInput(uint32  ident_,
                char   *read_,
                uint32  readLen_,
                uint32  bgn_,
                uint32  end_) {

    read       = new char [readLen_ + 1];
    readLength = readLen_;

    memcpy(read, read_, readLen_);

    ident      = ident_;

    placedBgn  = bgn_;
    placedEnd  = end_;

    alignedBgn = 0;
    alignedEnd = 0;
  };

  ~falconInput() {
    delete [] read;
  };

  char            *read;
  int32            readLength;

  uint32           ident;

  int32            placedBgn;
  int32            placedEnd;

  int32            alignedBgn;
  int32            alignedEnd;
};



class falconData {
public:
  falconData() {
    seq = NULL;
    eqv = NULL;
    pos = NULL;
    len = 0;
  };

  falconData(int32 len_) {
    seq    = new char  [len_];
    eqv    = new int32 [len_];
    pos    = new int32 [len_];
    len    = 0;
  };

  ~falconData() {
    delete [] seq;
    delete [] eqv;
    delete [] pos;
  };

  char            *seq;
  int32           *eqv;
  int32           *pos;    //  For tracking original read positions in corrected read.
  int32            len;
};



class falconConsensus {
public:
  falconConsensus(uint32               minOutputCoverage_,
                  uint32               minOutputLength_,
                  double               minOlapIdentity_,
                  uint32               minOlapLength_,
                  bool                 restrictToOverlap_ = true) {
    minOutputCoverage   = minOutputCoverage_;
    minOutputLength     = minOutputLength_;
    minOlapIdentity     = minOlapIdentity_;
    minOlapLength       = minOlapLength_;
    restrictToOverlap   = restrictToOverlap_;
    minRSS              = 0;
    maxRSS              = 0;
  };

  ~falconConsensus() {
  };

private:
  falconData *getConsensus(uint32         tagsLen,
                           alignTagList **tags,
                           uint32         templateLen);

private:
  void        setRSS(void) {
    uint64 active    = getBytesAllocated();

    //fprintf(stderr, "setRSS()-- %lu\n", active);

    minRSS = active;
    maxRSS = active;
  };

  void        updateRSS(void) {
    uint64 active    = getBytesAllocated();

    if (active < minRSS)
      minRSS = active;

    if (maxRSS < active)
      maxRSS = active;

    //fprintf(stdout, "updateRSS()-- minRSS %lu maxRSS %lu (delta %lu)\n", minRSS, maxRSS, active - before);
  };

public:
  uint64      getRSS(void) {
    return(maxRSS - minRSS);
  };



public:
  falconData *generateConsensus(falconInput         *evidence,
                                uint32               evidenceLen);

  uint64      estimateMemoryUsage(uint32 evidenceLen,
                                  uint64 nBasesInOlaps,
                                  uint32 templateLen);

  uint64      analyzeLength(tgTig        *layout,
                            uint32       &correctedLength,
                            uint64       &memoryNeeded);

private:
  uint32               minOutputCoverage;
  uint32               minOutputLength;

  double               minOlapIdentity;
  uint32               minOlapLength;

  bool                 restrictToOverlap;

  msa_vector_t         msa;

  uint64               minRSS;
  uint64               maxRSS;
};


#endif  //  FALCONCONSENSUS_H
