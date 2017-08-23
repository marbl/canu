
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/falcon_sense/libfalcon/falcon.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2017-AUG-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
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

#ifndef FALCONCONSENSUS_ALIGNTAG_H
#define FALCONCONSENSUS_ALIGNTAG_H


class falconInput;



class alignTag {       //  Completely initialized in get_align_tags()
public:
  alignTag() {
    q_id     = INT32_MAX;
    t_pos    = INT32_MAX;
    p_t_pos  = INT32_MAX;
    delta    = UINT16_MAX;
    p_delta  = UINT16_MAX;
    q_base   = '.';
    p_q_base = '.';
  };
  alignTag(uint32 id, int32  tp, int32 ptp, uint16 d, uint16 pd, char qb, char pqb) {
    q_id     = id;
    t_pos    = tp;
    p_t_pos  = ptp;
    delta    = d;
    p_delta  = pd;
    q_base   = qb;
    p_q_base = pqb;
  };


  uint32  q_id;

  int32   t_pos;
  int32   p_t_pos;         // the tag position of the previous base

  uint16  delta;
  uint16  p_delta;         // the tag delta of the previous base

  char    q_base;
  char    p_q_base;        // the previous base
};



//  A list of aligned tags for a single evidence read.
//
class alignTagList {
public:
  alignTagList(uint32 l) {
    tagsLen  = 0;
    tags     = (l == 0) ? NULL : new alignTag [l + 1];
  };

  ~alignTagList() {
    delete [] tags;
  };

  alignTag      *operator[](int32 i) { return(&tags[i]); };
  int32          numberOfTags(void)  { return(tagsLen);  };

  void           setTag(uint32 id, int32  tp, int32 ptp, uint16 d, uint16 pd, char qb, char pqb) {
    tags[tagsLen].q_id     = id;
    tags[tagsLen].t_pos    = tp;
    tags[tagsLen].p_t_pos  = ptp;
    tags[tagsLen].delta    = d;
    tags[tagsLen].p_delta  = pd;
    tags[tagsLen].q_base   = qb;
    tags[tagsLen].p_q_base = pqb;

    tagsLen++;
  };

  void           addTag(alignTag tag) {
    tags[tagsLen++] = tag;
  };

private:
  int32          tagsLen;
  alignTag      *tags;
};



alignTagList **
alignReadsToTemplate(falconInput    *evidence,
                     uint32          evidenceLen,
                     double          minIdentity);

#endif  //  FALCONCONSENSUS_ALIGNTAG_H
