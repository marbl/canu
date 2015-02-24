
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

#include "ovStore.H"

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
             ovsOverlap   *&ovl,
             uint32        &ovlLen,
             uint32        &ovlMax,
             ovStore       *ovsPrimary,
             ovStore       *ovsSecondary);

#endif  //  AS_OBT_OVERLAPS_H
