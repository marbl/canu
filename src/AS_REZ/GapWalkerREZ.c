
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

static char fileID[] = "$Id: GapWalkerREZ.c,v 1.14 2008-06-27 06:29:19 brianwalenz Exp $";

#include "AS_global.h"
#include "InputDataTypes_CGW.h"
#include "ScaffoldGraph_CGW.h"

LengthT
FindGapLength(ChunkInstanceT * lchunk,
              ChunkInstanceT * rchunk,
              int verbose) {
  LengthT gapSize;
  float lchunkMaxOffset, lchunkMaxVariance;
  float rchunkMinOffset, rchunkMinVariance;

  if (lchunk->offsetAEnd.mean < lchunk->offsetBEnd.mean) {
    lchunkMaxOffset = lchunk->offsetBEnd.mean;
    lchunkMaxVariance = lchunk->offsetBEnd.variance;
  } else {
    lchunkMaxOffset = lchunk->offsetAEnd.mean;
    lchunkMaxVariance = lchunk->offsetAEnd.variance;
  }

  if (rchunk->offsetAEnd.mean < rchunk->offsetBEnd.mean) {
    rchunkMinOffset = rchunk->offsetAEnd.mean;
    rchunkMinVariance = rchunk->offsetAEnd.variance;
  } else {
    rchunkMinOffset = rchunk->offsetBEnd.mean;
    rchunkMinVariance = rchunk->offsetBEnd.variance;
  }

  gapSize.mean = rchunkMinOffset - lchunkMaxOffset;
  gapSize.variance = rchunkMinVariance - lchunkMaxVariance;

  return gapSize;
}
