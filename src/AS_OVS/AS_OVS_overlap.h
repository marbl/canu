
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

#ifndef AS_OVS_OVERLAP_H
#define AS_OVS_OVERLAP_H

#include "AS_global.h"
#include "AS_MSG_pmesg.h"  //  pretty heavy just to get OverlapMesg.

#define AS_OVS_HNGBITS   12
#define AS_OVS_POSBITS   11
#define AS_OVS_ERRBITS   12

//  Convert q between a condensed/encoded integer and a floating point
//  value.
//
//  Q should be a floating point value between 0.000 and 1.000, and is
//  the fraction error in this alignment.  We are able to encode error
//  up to 0.4000 (40%), with up to four significant figures.
//
//  Previous versions of the overlap store stored any error, with
//  three significant figures, but used 16 bits to do it.  You can get
//  the same effect by using 1000.0 instead of 10000.0 below.

#define AS_OVS_MAX_ERATE          ((1 << AS_OVS_ERRBITS) - 1)

#define AS_OVS_decodeQuality(E)   ((E) / 10000.0)
#define AS_OVS_encodeQuality(Q)   (((Q) < AS_OVS_decodeQuality(AS_OVS_MAX_ERATE)) ? (int)(10000.0 * (Q) + 0.5) : AS_OVS_MAX_ERATE)

#define AS_OVS_TYPE_OVL   0x00
#define AS_OVS_TYPE_OBT   0x01
#define AS_OVS_TYPE_MER   0x02
#define AS_OVS_TYPE_UNS   0x03
#define AS_OVS_TYPE_ANY   0xff

typedef union {
  uint64   dat;
  struct {
    uint64  datpad:5;
    uint64  seed_value:8;
    uint64  flipped:1;
    int64   a_hang:AS_OVS_HNGBITS;
    int64   b_hang:AS_OVS_HNGBITS;
    uint64  orig_erate:AS_OVS_ERRBITS;
    uint64  corr_erate:AS_OVS_ERRBITS;
    uint64  type:2;
  } ovl;
  struct {
    uint64  datpad:19;
    uint64  compression_length:3;
    uint64  fwd:1;
    uint64  palindrome:1;
    uint64  a_pos:AS_OVS_POSBITS;
    uint64  b_pos:AS_OVS_POSBITS;
    uint64  k_count:8;
    uint64  k_len:8;
    uint64  type:2;
  } mer;
  struct {
    uint64  datpad:5;
    uint64  fwd:1;
    uint64  a_beg:AS_OVS_POSBITS;
    uint64  a_end:AS_OVS_POSBITS;
    uint64  b_beg:AS_OVS_POSBITS;
    uint64  b_end:AS_OVS_POSBITS;
    uint64  erate:AS_OVS_ERRBITS;
    uint64  type:2;
  } obt;
} OVSoverlapDAT;

typedef struct {
  uint32         b_iid;
  OVSoverlapDAT  dat;
} OVSoverlapINT;

typedef struct {
  uint32         a_iid;
  uint32         b_iid;
  OVSoverlapDAT  dat;
} OVSoverlap;


void  AS_OVS_convertOverlapMesgToOVSoverlap(OverlapMesg *omesg, OVSoverlap *ovs);
int   AS_OVS_convertOVLdumpToOVSoverlap(char *line, OVSoverlap *olap);
int   AS_OVS_convertOBTdumpToOVSoverlap(char *line, OVSoverlap *olap);


#endif  //  AS_OVS_OVERLAP_H
