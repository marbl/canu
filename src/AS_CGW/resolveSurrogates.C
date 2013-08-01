
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

const char *mainid = "$Id: resolveSurrogates.c,v 1.23 2009-09-12 22:35:58 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.H"
#include "AS_CGW_dataTypes.H"
#include "ScaffoldGraph_CGW.H"
#include "ScaffoldGraphIterator_CGW.H"
#include "fragmentPlacement.H"


int
main(int argc, char **argv) {

  // if 1, aggressively place fragments in surrogates that are only
  // used once in the assembly; "aggressively" means place all the
  // fragments in the unitig, regardless of mate status, alignment
  // quality etc
  //
  int     placeAllFragsInSinglePlacedSurros = 0;
  double  cutoffToInferSingleCopyStatus     = 1.0;
  int     ckptNum                           = 0;

  GlobalData = new Globals_CGW();

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->outputPrefix, argv[++arg]);
    } else if (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->gkpStoreName, argv[++arg]);
    } else if (strcmp(argv[arg], "-n") == 0) {
      ckptNum = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-S") == 0) {
      cutoffToInferSingleCopyStatus=atof(argv[++arg]);
    } else if (strcmp(argv[arg], "-1") == 0) {
      placeAllFragsInSinglePlacedSurros = 1;
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err) ||
      (GlobalData->outputPrefix[0] == 0) ||
      (GlobalData->gkpStoreName[0] == 0) ||
      (ckptNum == 0)) {
    fprintf(stderr, "usage: %s -g <gkp> -c <ckp> -n <num> opts\n",argv[0]);
    fprintf(stderr, "  -S x   place all frags in singly-placed surrogates if\n");
    fprintf(stderr, "         at least fraction x can be placed.\n");
    fprintf(stderr, "  -1     place all frags in singly-placed surrogates\n");
    fprintf(stderr, "         aggressively; equivalent to -S 0.0\n");
    exit(1);
  }

  LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix, ckptNum, TRUE);

  resolveSurrogates(placeAllFragsInSinglePlacedSurros, cutoffToInferSingleCopyStatus);

  CheckpointScaffoldGraph("resolveSurrogates", "after resolveSurrogates");

  delete GlobalData;

  exit(0);
}
