
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
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

const char *mainid = "$Id: getNumScaffolds.c,v 1.8 2008-10-08 22:02:57 brianwalenz Exp $";

//  Not part of Overlap Based Trimming proper, but used by the script
//  that runs OBT.  The build system cannot (apparently) handle
//  buildnig executables in subdirectories (e.g., AS_RUN/runCA-OBT/)
//  which is where this should go.

#include <stdio.h>

#include "AS_global.h"

#include "Globals_CGW.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"

int
main(int argc , char **argv) {

  argc = AS_configure(argc, argv);

  if (argc != 4) {
    fprintf(stderr, "usage: %s gkpStore cgw-checkpoint-prefix checkpoint-number\n", argv[0]);
    exit(1);
  }

  GlobalData = CreateGlobal_CGW();
  GlobalData->stderrc  = stderr;
  GlobalData->timefp   = stderr;
  GlobalData->aligner  = DP_Compare;

  strcpy(GlobalData->Gatekeeper_Store_Name, argv[1]);
  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(argv[2], atoi(argv[3]), FALSE);
  fprintf(stdout, "%d\n", GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));

  exit(0);
}

