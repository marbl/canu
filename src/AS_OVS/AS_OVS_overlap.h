
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

#define AS_OVS_HNGBITS   12
#define AS_OVS_POSBITS   11
#define AS_OVS_ERRBITS   12

//  With 16 bits for storing the error, we can store up to 65% error,
//  with three decimal points of precision, 65.XXX.  If space become
//  tight later on, we could store less precise overlaps, 12 bits
//  would get us 40.XX.

#define AS_OVS_MAX_ERATE        ((1 << AS_OVS_ERRBITS) - 1)

//  Convert q from/to a condensed form / floating point equivalent
//
#define Expand_Quality(Q)   ((Q) / 100.0)
#define Shrink_Quality(Q)   (((Q) < Expand_Quality(AS_OVS_MAX_ERATE)) ? (int)(100.0 * (Q) + 0.5) : AS_OVS_MAX_ERATE)

#define AS_OVS_TYPE_OVL   0x00
#define AS_OVS_TYPE_OBT   0x01
#define AS_OVS_TYPE_MER   0x02
#define AS_OVS_TYPE_UNS   0x03

typedef union {
  uint64   dat;
  struct {
    uint64  datpad:13;
    uint64  flipped:1;
    int64   a_hang:AS_OVS_HNGBITS;
    int64   b_hang:AS_OVS_HNGBITS;
    uint64  orig_erate:AS_OVS_ERRBITS;
    uint64  corr_erate:AS_OVS_ERRBITS;
    uint64  type:2;
  } ovl;
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
  struct {
    uint64  datpad:20;
    uint64  fwd:1;
    uint64  palindrome:1;
    uint64  a_pos:AS_OVS_POSBITS;
    uint64  b_pos:AS_OVS_POSBITS;
    uint64  k_count:10;
    uint64  k_len:8;
    uint64  type:2;
  } mer;
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

#endif  //  AS_OVS_OVERLAP_H
