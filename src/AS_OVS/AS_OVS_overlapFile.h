
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

#ifndef AS_OVS_OVERLAPFILE_H
#define AS_OVS_OVERLAPFILE_H

#include <stdio.h>

#include "AS_global.h"
#include "AS_OVS_overlap.h"

typedef struct {
  int           bufferLen;  //  length of valid data in the buffer
  int           bufferPos;  //  position the read is at in the buffer
  int           bufferMax;  //  allocated size of the buffer
  uint32       *buffer;
  int           isOutput;   //  if true, we can AS_OVS_writeOverlap()
  int           isSeekable; //  if true, we can AS_OVS_seekOverlap()
  int           isPopened;  //  if true, we need to pclose()
  int           isInternal; //  if true, 3 words per overlap, else 4
  FILE         *file;
} BinaryOverlapFile;

BinaryOverlapFile *AS_OVS_openBinaryOverlapFile(const char *name, int isInternal);
BinaryOverlapFile *AS_OVS_createBinaryOverlapFile(const char *name, int isInternal);

void               AS_OVS_flushBinaryOverlapFile(BinaryOverlapFile *bof);
void               AS_OVS_closeBinaryOverlapFile(BinaryOverlapFile *bof);

void               AS_OVS_writeOverlap(BinaryOverlapFile *bof, OVSoverlap *overlap);
int                AS_OVS_readOverlap(BinaryOverlapFile *bof, OVSoverlap *overlap);

void               AS_OVS_seekOverlap(BinaryOverlapFile *bof, uint32 overlap);

#endif  //  AS_OVS_OVERLAPFILE_H
