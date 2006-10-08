
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
  Global_CGW *data;

  GlobalData = data = CreateGlobal_CGW();
  data->stderrc  = stderr;
  data->timefp   = stderr;
  data->aligner  = DP_Compare;

  if (argc != 3) {
    fprintf(stderr, "usage: %s cgw-checkpoint-prefix checkpoint-number\n", argv[0]);
    exit(1);
  }

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(argv[1], atoi(argv[2]), FALSE);
  fprintf(stdout, "%d\n", GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));

  exit(0);
}

