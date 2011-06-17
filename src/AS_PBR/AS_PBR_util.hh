/*
Copyright (C) 2011, Battelle National Biodefense Institute (BNBI);
all rights reserved. Authored by: Sergey Koren

This Software was prepared for the Department of Homeland Security
(DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
part of contract HSHQDC-07-C-00020 to manage and operate the National
Biodefense Analysis and Countermeasures Center (NBACC), a Federally
Funded Research and Development Center.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of the Battelle National Biodefense Institute nor
  the names of its contributors may be used to endorse or promote
  products derived from this software without specific prior written
  permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

using namespace std;

#include "AS_global.h"
#include "AS_OVS_overlapStore.h"

uint32 olapLengthOVL(OVSoverlap ovl, uint32 alen, uint32 blen) {
  int32   ah = ovl.dat.ovl.a_hang;
  int32   bh = ovl.dat.ovl.b_hang;
  uint32  le = 0;

  if (ah < 0) {
    if (bh < 0)
      le = alen + bh;
    else
      le = blen + ah - bh;
  } else {
    if (bh < 0)
      le = alen + bh - ah;
    else
      le = alen - ah;
  }

  return(le);
}

uint32 olapLengthOBT(OVSoverlap ovl, uint32 alen, uint32 blen) {
   uint32 aovl = ovl.dat.obt.a_end - ovl.dat.obt.a_beg;
   uint32 bend = ovl.dat.obt.b_end_hi >> 9 | ovl.dat.obt.b_end_lo;
   uint32 bbgn = MIN(ovl.dat.obt.b_beg, bend);
   bend = MAX(ovl.dat.obt.b_beg, bend);

   return MIN(aovl, (bend-bbgn));
}

extern uint32
olapLength(OVSoverlap ovl, uint32 alen, uint32 blen) {
   if (ovl.dat.ovl.type == AS_OVS_TYPE_OVL) {
      return olapLengthOVL(ovl, alen, blen);
   } else if (ovl.dat.ovl.type == AS_OVS_TYPE_OBT) {
      return olapLengthOBT(ovl, alen, blen);
   }
   return 0;
}

bool isOlapBadOVL(OVSoverlap olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
    //  The overlap is ALWAYS bad if the original error rate is above what we initially required
    //  overlaps to be at.  We shouldn't have even seen this overlap.  This is a bug in the
    //  overlapper.
    //
    if (olap.dat.ovl.orig_erate > AS_OVS_encodeQuality(AS_CNS_ERROR_RATE))
       return true;

    //  The overlap is GOOD (false == not bad) if the corrected error rate is below the requested
    //  erate.
    //
    if (olap.dat.ovl.corr_erate <= AS_OVS_encodeQuality(AS_UTG_ERROR_RATE)) {
      return false;
    }

    //  If we didn't allow fixed-number-of-errors, the overlap is now bad.  Just a slight
    //  optimization.
    //
    if (AS_UTG_ERROR_LIMIT <= 0)
      return true;

    double olen = olapLength(olap, alen, blen);
    double nerr = olen * AS_OVS_decodeQuality(olap.dat.ovl.corr_erate);

    assert(nerr >= 0);

    if (nerr <= AS_UTG_ERROR_LIMIT) {
       return false;
    }

    return true;
}

bool isOlapBadOBT(OVSoverlap olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
   if (olap.dat.obt.erate > AS_OVS_encodeQuality(AS_CNS_ERROR_RATE))
      return true;

   if (olap.dat.obt.erate <= AS_OVS_encodeQuality(AS_UTG_ERROR_RATE))
      return false;

   if (AS_UTG_ERROR_LIMIT <= 0) 
      return true;

   double olen = olapLength(olap, alen, blen);
   double nerr = olen * AS_OVS_decodeQuality(olap.dat.obt.erate);
   if (nerr <= AS_UTG_ERROR_LIMIT)
      return false;

   return true;
}

extern bool isOlapBad(OVSoverlap olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
    if (olap.dat.ovl.type == AS_OVS_TYPE_OVL) {
       return isOlapBadOVL(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE);
    } else if (olap.dat.ovl.type == AS_OVS_TYPE_OBT) {
       return isOlapBadOBT(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE);
    }
    return true;
}

uint64  scoreOverlapOVL(const OVSoverlap& olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {

    //  BPW's newer new score.  For the most part, we use the length of the overlap, but we also
    //  want to break ties with the higher quality overlap.
    //
    //  The high 20 bits are the length of the overlap.
    //  The next 12 are the corrected error rate.
    //  The last 12 are the original error rate.
    //
    //  (Well, 12 == AS_OVS_ERRBITS)
    if (isOlapBad(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE)) {
       return 0;
    }

    uint64  leng = 0;
    uint64  corr = (AS_OVS_MAX_ERATE - olap.dat.ovl.corr_erate);
    uint64  orig = (AS_OVS_MAX_ERATE - olap.dat.ovl.orig_erate);

    //  Shift AFTER assigning to a 64-bit value to avoid overflows.
    corr <<= AS_OVS_ERRBITS;

    //  Containments - the length of the overlaps are all the same.  We return the quality.
    //
    if (((olap.dat.ovl.a_hang >= 0) && (olap.dat.ovl.b_hang <= 0)) ||
        ((olap.dat.ovl.a_hang <= 0) && (olap.dat.ovl.b_hang >= 0)))
      return(corr | orig);

    //  Dovetails - the length of the overlap is the score, but we bias towards lower error.
    //  (again, shift AFTER assigning to avoid overflows)
    //
    leng   = olapLength(olap, alen, blen);
    leng <<= (2 * AS_OVS_ERRBITS);

    return(leng | corr | orig);
}

uint64  scoreOverlapOBT(const OVSoverlap& olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
    if (isOlapBad(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE)) {
       return 0;
    }

    uint64  leng = 0;
    uint64  orig = (AS_OVS_MAX_ERATE - olap.dat.obt.erate);

    leng   = olapLength(olap, alen, blen);
    leng <<= AS_OVS_ERRBITS;

    return(leng | orig);
}

extern
  uint64  scoreOverlap(const OVSoverlap& olap, uint32 alen, uint32 blen, double AS_UTG_ERROR_RATE, double AS_UTG_ERROR_LIMIT, double AS_CNS_ERROR_RATE) {
   if (olap.dat.ovl.type == AS_OVS_TYPE_OVL) {
      return scoreOverlapOVL(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE);
   } else if (olap.dat.ovl.type == AS_OVS_TYPE_OBT) {
      return scoreOverlapOBT(olap, alen, blen, AS_UTG_ERROR_RATE, AS_UTG_ERROR_LIMIT, AS_CNS_ERROR_RATE);
   }
   return 0;
}
