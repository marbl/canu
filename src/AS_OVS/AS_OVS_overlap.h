
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


#define AS_OVS_POSBITS   11
#define AS_OVS_ERRBITS   16

#define MAX_ERATE        ((1 << AS_OVS_ERRBITS) - 1)


//  Convert q from/to a condensed form / floating point equivalent
//
#define Expand_Quality(Q)   ((Q) / 1000.0)
#define Shrink_Quality(Q)   (((Q) < Expand_Quality(MAX_ERATE)) ? (int)(1000.0 * (Q) + 0.5) : MAX_ERATE)


typedef union {
  uint64   dat;
  struct {
    uint64  datpad:9;
    uint64  flipped:1;
    int64   Ahang:AS_OVS_POSBITS;
    int64   Bhang:AS_OVS_POSBITS;
    uint64  origE:AS_OVS_ERRBITS;
    uint64  corrE:AS_OVS_ERRBITS;
  } ovl;
  struct {
    uint64  datpad:3;
    uint64  ori:1;
    uint64  Abeg:AS_OVS_POSBITS;
    uint64  Aend:AS_OVS_POSBITS;
    uint64  Bbeg:AS_OVS_POSBITS;
    uint64  Bend:AS_OVS_POSBITS;
    uint64  erate:AS_OVS_ERRBITS;
  } obt;
  struct {
    uint64  datpad:22;
    uint64  ori:1;
    uint64  palindrome:1;
    uint64  Apos:AS_OVS_POSBITS;
    uint64  Bpos:AS_OVS_POSBITS;
    uint64  kcount:10;
    uint64  klen:8;
  } mer;
} OVSoverlapDAT;


typedef struct {
  uint32         bid;
  OVSoverlapDAT  dat;
} OVSoverlapINT;


typedef struct {
  uint32         aid;
  uint32         bid;
  OVSoverlapDAT  dat;
} OVSoverlap;


#endif  //  AS_OVS_OVERLAP_H


