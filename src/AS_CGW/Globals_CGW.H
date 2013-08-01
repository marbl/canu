
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

#ifndef GLOBALS_CGW_H
#define GLOBALS_CGW_H

static const char *rcsid_GLOBALS_CGW_H = "$Id$";

#include "AS_global.H"

//  These are the global data structures for the CGW Module

class Globals_CGW {
public:
  Globals_CGW();
  ~Globals_CGW();

  //  Returns last ckp number for this runCA directory.
  int    setPrefix(char *runCAroot);

  int    verbose;
  int    debugLevel;

  int    repeatRezLevel;
  int    stoneLevel;

  int32  minWeightToMerge;

  uint32 shatterLevel;
  double mergeScaffoldMissingMates;

  int    minSamplesForOverride;

  int    outputOverlapOnlyContigEdges;
  int    demoteSingletonScaffolds;
  int    checkRepeatBranchPattern;

  int    ignoreChaffUnitigs;
  int    performCleanupScaffolds;

  int    cgbUniqueCutoff;
  int    cgbDefinitelyUniqueCutoff;
  int    cgbApplyMicrohetCutoff;
  double cgbMicrohetProb;

  int    allowDemoteMarkedUnitigs;

  int    closurePlacement;
  
  int    removeNonOverlapingContigsFromScaffold;
  int    doUnjiggleWhenMerging;

  int    mergeFilterLevel;

  char   outputPrefix[FILENAME_MAX];

  char   gkpStoreName[FILENAME_MAX];
  char   ovlStoreName[FILENAME_MAX];
  char   tigStoreName[FILENAME_MAX];

  char   unitigOverlaps[FILENAME_MAX];
};

extern Globals_CGW  *GlobalData;

#endif
