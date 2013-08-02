
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2013, J. Craig Venter Institute.
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

#ifndef AS_OBT_OVERLAPS_H
#define AS_OBT_OVERLAPS_H

#include "AS_global.H"

#include "AS_OVS_overlapStore.H"

//  Class to load all ovelaps for some read, optionally adjust to
//  an updated OBT clear range, and optionally filter out suspected
//  subread containing reads.
//
//  Initial trim does not use overlaps.
//  Final trim uses just the basic OBT overlap.
//  Chimer/spur detection also needs the 'hang'.
//  Dedupe also uses just the basic OBT overlap.

uint32
loadOverlaps(uint32         iid,
             OVSoverlap   *&ovl,
             uint32        &ovlLen,
             uint32        &ovlMax,
             OverlapStore  *ovlPrimary,
             OverlapStore  *ovlSecondary);


#if 0
class overlapForOBT {
public:
  overlapForOBT();

  overlapForOBT(uint64 Aiid_, uint32 Alhang_, uint32 Abeg_, uint32 Aend_, uint32 Arhang_,
                uint64 Biid_, uint32 Blhang_, uint32 Bbeg_, uint32 Bend_, uint32 Brhang_, char ori_);

  void   print(FILE *out) const;


public:
  uint64       flipped  : 1;
  uint64       deleted  : 1;
  uint64       style    : 4;  //  For chimer/spur detection

  uint64       Alhang   : AS_READ_MAX_NORMAL_LEN_BITS;  //  For chimer/spur detection
  uint64       Abeg     : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64       Aend     : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64       Arhang   : AS_READ_MAX_NORMAL_LEN_BITS;  //  For chimer/spur detection

  uint64       Blhang   : AS_READ_MAX_NORMAL_LEN_BITS;  //  For chimer/spur detection
  uint64       Bbeg     : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64       Bend     : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64       Brhang   : AS_READ_MAX_NORMAL_LEN_BITS;  //  For chimer/spur detection

  uint64       erate    : AS_OVS_ERR_BITS;

  uint32       Aiid;
  uint32       Biid;
};
#endif



#endif  //  AS_OBT_OVERLAPS_H
