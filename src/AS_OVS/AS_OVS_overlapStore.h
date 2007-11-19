
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

#ifndef AS_OVS_OVERLAPSTORE_H
#define AS_OVS_OVERLAPSTORE_H

#include <stdio.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"

typedef struct {
  uint64    ovsMagic;
  uint64    ovsVersion;
  uint64    numOverlapsPerFile;  //  on create, how big to make each file
  uint64    smallestIID;         //  smallest frag iid in the store
  uint64    largestIID;          //  largest frag iid in the store
  uint64    numOverlapsTotal;    //  number of overlaps in the store
  uint64    highestFileIndex;
} OverlapStoreInfo;

typedef struct {
  uint64   nSamples;
  uint64   histogram[65536];     //  Assumes input values are 16-bit unsigned ints!

  uint16   median;

  double   mean;
  double   stddev;

  uint16   mode;
  double   mad;
} OverlapStoreHistogram;

typedef struct {
  uint64                   numOVL;
  uint64                   numOBT;
  uint64                   numMER;

  //  Specific to OVL type
  OverlapStoreHistogram    orig_erate[16];
  OverlapStoreHistogram    corr_erate[16];
  OverlapStoreHistogram    length[16];

  //  Specific to OBT type
  OverlapStoreHistogram    obtAbeg;
  OverlapStoreHistogram    obtAend;
  OverlapStoreHistogram    obtAlength;
  OverlapStoreHistogram    obtBbeg;
  OverlapStoreHistogram    obtBend;
  OverlapStoreHistogram    obtBlength;
  OverlapStoreHistogram    obtErate;

  //  Specific to MER type
  OverlapStoreHistogram    merApos;
  OverlapStoreHistogram    merBpos;
  OverlapStoreHistogram    merKcount;
} OverlapStoreStats;

typedef struct {
  uint32    a_iid;
  uint32    fileno;    //  the file that contains this a_iid
  uint32    offset;    //  offset to the first overlap for this iid
  uint32    numOlaps;  //  number of overlaps for this iid
} OverlapStoreOffsetRecord;

typedef struct {
  char                        storePath[FILENAME_MAX];
  int                         isOutput;
  char                        useBackup;
  int                         saveSpace;

  OverlapStoreInfo            ovs;

  OverlapStoreStats           stats;
  int                         statsUpdated;

  FILE                       *offsetFile;
  OverlapStoreOffsetRecord    offset;
  OverlapStoreOffsetRecord    missing;

  uint32                      firstIIDrequested;
  uint32                      lastIIDrequested;

  int                         overlapsThisFile;
  int                         currentFileIndex;
  BinaryOverlapFile          *bof;

  GateKeeperStore            *gkp;

  uint16                     *fragClearBegin;
  uint16                     *fragClearEnd;
  uint16                     *fragClearLength;
} OverlapStore;

OverlapStore      *AS_OVS_openOverlapStorePrivate(const char *name, int useBackup, int saveSpace);
void               AS_OVS_closeOverlapStore(OverlapStore *ovs);
void               AS_OVS_restoreBackup(OverlapStore *ovs);

#define            AS_OVS_openOverlapStore(N)  AS_OVS_openOverlapStorePrivate((N), FALSE, FALSE)

int                AS_OVS_readOverlapFromStore(OverlapStore *ovs, OVSoverlap *overlap, uint32 type);

void               AS_OVS_setRangeOverlapStore(OverlapStore *ovs, uint32 low, uint32 high);
void               AS_OVS_resetRangeOverlapStore(OverlapStore *ovs);

uint64             AS_OVS_numOverlapsInRange(OverlapStore *ovs);

static
uint32             AS_OVS_lastFragInStore(OverlapStore *ovs) {
  return(ovs->ovs.largestIID);
}

//  The mostly private interface for dealing with stats

void   AS_OVS_histogramAdd(OverlapStoreHistogram *h, uint16 val);
void   AS_OVS_histogramCompute(OverlapStoreHistogram *h);
void   AS_OVS_histogramShow(char *label, char *type, OverlapStoreHistogram *h);

void   AS_OVS_accumulateStats(OverlapStore *ovs, OVSoverlap *ovl);



//  The mostly private interface for creating an overlap store.

OverlapStore      *AS_OVS_createOverlapStore(const char *name, int failOnExist);
void               AS_OVS_writeOverlapToStore(OverlapStore *ovs, OVSoverlap *olap);

#endif  //  AS_OVS_OVERLAPSTORE_H
