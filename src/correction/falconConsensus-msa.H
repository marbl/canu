
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

#ifndef FALCONCONSENSUS_MSA_H
#define FALCONCONSENSUS_MSA_H

class align_tag_col_t {
public:
  align_tag_col_t() {
    size          =  8;
    p_t_pos       =  new int32   [size];
    p_delta       =  new uint16  [size];
    p_q_base      =  new char    [size];
    link_count    =  new uint16  [size];
    clean();
  };

  ~align_tag_col_t() {
    delete [] p_t_pos;
    delete [] p_delta;
    delete [] p_q_base;
    delete [] link_count;
  };

  void   clean(void) {
    n_link         =  0;
    count          =  0;
    best_p_t_pos   = -1;
    best_p_delta   = -1;
    best_p_q_base  = -1;
    score          =  DBL_MIN;
  };

  void  addEntry(alignTag *tag) {

    if (n_link >= size) {
      int32  ns;  //  Need to keep resetting this back to the real number of elements in each array

      ns = size;  resizeArray(p_t_pos,    n_link, ns, size + 16);
      ns = size;  resizeArray(p_delta,    n_link, ns, size + 16);
      ns = size;  resizeArray(p_q_base,   n_link, ns, size + 16);
      ns = size;  resizeArray(link_count, n_link, ns, size + 16);

      size = size + 16;
    }

    p_t_pos   [n_link]  = tag->p_t_pos;
    p_delta   [n_link]  = tag->p_delta;
    p_q_base  [n_link]  = tag->p_q_base;
    link_count[n_link]  = 1;

    n_link++;
  };

  double     score;

  int32     *p_t_pos;        // the tag position of the previous base
  uint16    *p_delta;        // the tag delta of the previous base
  char      *p_q_base;       // the previous base
  uint16    *link_count;

  int32      best_p_t_pos;

  uint16     best_p_delta;
  uint16     best_p_q_base;  // encoded base
  uint16     count;          //  Number of times we've encountered this base
  uint16     size;           //  Number of items allocated in the arrays
  uint16     n_link;         //  Number of items used in the arrays
};



class  msa_base_group_t {
public:
  msa_base_group_t() {
  };

  ~msa_base_group_t() {
  };

  void                clean(void) {
    base[0].clean();  //  'A'
    base[1].clean();  //  'C'
    base[2].clean();  //  'G'
    base[3].clean();  //  'T'
    base[4].clean();  //  '-' and everything else
  }

  align_tag_col_t  base[5];
};



class msa_delta_group_t {
public:
  msa_delta_group_t() {
    coverage     = 0;

    deltaAlloc   = 8;                                    //  Allocate pointers for 8 new items
    deltaLen     = 0;
    delta        = new msa_base_group_t * [deltaAlloc];

    deltaPtrsLen = 1;
    deltaPtrs[0] = new msa_base_group_t   [deltaAlloc];  //  Allocate a block of items

    for (uint32 ii=1; ii<128; ii++)                      //  Clear the other pointers.
      deltaPtrs[ii] = 0;

    for (uint32 ii=0; ii<deltaAlloc; ii++)               //  And set pointers to the items
      delta[ii] = deltaPtrs[0] + ii;
  };

  ~msa_delta_group_t() {
    for (uint32 ii=0; ii<128; ii++)
      delete [] deltaPtrs[ii];

    delete [] delta;
  };

  void       increaseDeltaGroup(uint16 newMax) {
    uint32  newLen = newMax + 1;

    if (newLen <= deltaLen)    //  Requested group is already used.
      return;

    if (newLen < deltaAlloc) { //  Requested group is already allocated.
      deltaLen = newLen;
      return;
    }

    //  Figure out how much we want to grow.  deltaAlloc is limited to 65,536, and we have 128
    //  pointers to play with.  The first has 8.  We'll do a bunch of size 8, then size 16, then 32,
    //  and finally give up and do 1024.

    uint32  incr = 0;

    if      (deltaPtrsLen < 16)    //  128 here
      incr = 8;
    else if (deltaPtrsLen < 32)   //   256 here
      incr = 16;
    else if (deltaPtrsLen < 64)   //  1024 here
      incr = 32;
    else                          // 65535 here
      incr = 1024;

    //  Make sure it's really big enough.

    while (deltaAlloc + incr < newLen)
      incr *= 2;

    //  Then reallocate pointers and grab new blocks.  We copy all deltaAlloc elements, not just
    //  deltaLen, because (a) deltaLen is the last used position, and (b) we've already allocated
    //  space up to deltaAlloc.

    resizeArray(delta, deltaAlloc, deltaAlloc, deltaAlloc+incr);

    deltaPtrs[deltaPtrsLen] = new msa_base_group_t [incr];

    //  Finally, assign pointers.  (deltaAlloc-incr because incr was already added to it)

    for (uint32 ii=0; ii<incr; ii++)
      delta[deltaAlloc+ii-incr] = deltaPtrs[deltaPtrsLen] + ii;

    deltaPtrsLen++;

    //deltaAlloc += incr;   //  Done by resizeArray().
    deltaLen = newLen;
  };


  void    clean(void) {
    for (uint32 j=0; j<deltaAlloc; j++)  //  Would deltaLen suffice?
      delta[j]->clean();

    coverage = 0;
    deltaLen = 0;
  }


  uint16             coverage;
  uint16             deltaAlloc;       //  Size of 'delta' array
  uint16             deltaLen;         //  Number of 'delta' positions actually used
  uint16             deltaPtrsLen;     //  Next available spot 'deltaPtrs';

  msa_base_group_t **delta;            //  Just pointers into deltaPtrs;
  msa_base_group_t  *deltaPtrs[128];   //  Holds the actual delta storage
};



class msa_vector_t {
public:
  msa_vector_t() {
    dgLen = 0;
    dgMax = 0;
    dg    = NULL;
  };

  ~msa_vector_t() {
    delete [] dg;
  };

  void    resize(uint32 templateLen) {
    dgLen = templateLen;

    if (dgMax < dgLen) {
      delete [] dg;

      dgMax = dgLen;
      dg    = new msa_delta_group_t [dgMax];
    }

    for (uint32 i=0; i<dgLen; i++)    //  Clean out old data
      dg[i].clean();
  };

  msa_delta_group_t  *operator[](int32 i) {
    assert(i < dgLen);
    return(dg + i);
  };

private:
  uint32              dgLen;    //  Last used.
  uint32              dgMax;    //  Space allocated.
  msa_delta_group_t  *dg;
};

#endif  //  FALCONCONSENSUS_MSA_H
