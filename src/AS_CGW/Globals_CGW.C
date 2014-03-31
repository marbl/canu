
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

static char *rcsid = "$Id$";

#include "AS_global.H"
#include "Globals_CGW.H"

#include "AS_UTL_fileIO.H"

Globals_CGW *GlobalData = NULL;

Globals_CGW::Globals_CGW() {
  verbose                                 = 0;
  debugLevel                              = 0;

  repeatRezLevel                          = repeatRezLevel;
  stoneLevel                              = 0;

  minWeightToMerge                        = 2;

  shatterLevel                            = 0;
  mergeScaffoldMissingMates               = 0;

  minSamplesForOverride                   = 100;

  outputOverlapOnlyContigEdges            = 0;
  demoteSingletonScaffolds                = TRUE;
  checkRepeatBranchPattern                = FALSE;

  ignoreChaffUnitigs                      = 0;
  performCleanupScaffolds                 = 1;

  cgbUniqueCutoff                         = CGB_UNIQUE_CUTOFF;
  cgbDefinitelyUniqueCutoff               = CGB_UNIQUE_CUTOFF;

  allowDemoteMarkedUnitigs                = TRUE;       // allow toggled unitigs to be demoted to be repeat if they were marked unique

  closurePlacement                        = 0;

  removeNonOverlapingContigsFromScaffold  = 0;
  doUnjiggleWhenMerging                   = 0;

  //  Generally, higher values are more strict.
  mergeFilterLevel                        = 1;

  memset(outputPrefix, 0, FILENAME_MAX);

  memset(gkpStoreName, 0, FILENAME_MAX);
  memset(ovlStoreName, 0, FILENAME_MAX);
  memset(tigStoreName, 0, FILENAME_MAX);

  memset(unitigOverlaps, 0, FILENAME_MAX);
}


Globals_CGW::~Globals_CGW() {
}


int
Globals_CGW::setPrefix(char *runCAroot) {
  char  testname[FILENAME_MAX];
  int   ckp = -1;

  sprintf(GlobalData->outputPrefix, "7-CGW/%s",    runCAroot);
  sprintf(GlobalData->gkpStoreName, "%s.gkpStore", runCAroot);
  sprintf(GlobalData->ovlStoreName, "%s.ovlStore", runCAroot);
  sprintf(GlobalData->tigStoreName, "%s.tigStore", runCAroot);

  //  Find the checkpoint number by testing what files open.  We assume checkpoints are numbered
  //  contiguously, and stop after the first non-contiguous block -- e.g., "4, 5, 6" would return 6.

  for (int32 i=0; i<1024; i++) {
    sprintf(testname, "%s.ckp.%d", GlobalData->outputPrefix, i);

    if (AS_UTL_fileExists(testname, FALSE, FALSE))
      ckp = i;
    else
      if (ckp != -1)
        break;
  }

  if (ckp == -1)
    fprintf(stderr, "Globals_CGW::setPrefix()-- I couldn't find any checkpoints in '7-CGW/%s/'.\n",
            runCAroot), exit(1);

  return(ckp);
}
