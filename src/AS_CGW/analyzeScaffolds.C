
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

const char *mainid = "$Id: analyzeScaffolds.C,v 1.2 2012-07-16 09:13:17 brianwalenz Exp $";

#include "AS_global.H"

#include "Globals_CGW.H"
#include "AS_CGW_dataTypes.H"
#include "ScaffoldGraph_CGW.H"
#include "ScaffoldGraphIterator_CGW.H"

#include "CIScaffoldT_Analysis.H"


int main (int argc, char *argv[]) {
  int32       checkpointVers           = 0;
  int32       tigStoreVers             = 0;

  GlobalData = new Globals_CGW();

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->gkpStoreName, argv[++arg]);

    } else if (strcmp(argv[arg], "-t") == 0) {
      strcpy(GlobalData->tigStoreName, argv[++arg]);
      tigStoreVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->outputPrefix, argv[++arg]);
      checkpointVers = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((GlobalData->gkpStoreName[0] == 0) ||
      (GlobalData->tigStoreName[0] == 0) ||
      (err)) {
    fprintf(stderr, "usage: %s -g gkpStore [-o prefix] [-s firstUID] [-n namespace] [-E server] [-h]\n", argv[0]);
    fprintf(stderr, "  -g gkpStore             mandatory path to the gkpStore\n");
    fprintf(stderr, "  -t tigStore version     mandatory path to the tigStore and version\n");
    fprintf(stderr, "  -c checkpoint version   optional path to a checkpoint and version\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  LoadScaffoldGraphFromCheckpoint(GlobalData->outputPrefix, checkpointVers, FALSE);


  vector<instrumentLIB>   lib;

  for (int32 i=0; i<GetNumDistTs(ScaffoldGraph->Dists); i++) {
    DistT *dptr = GetDistT(ScaffoldGraph->Dists, i);

    lib.push_back(instrumentLIB(i, dptr->mu, dptr->sigma, true));
  }

  GraphNodeIterator   scaffolds;
  CIScaffoldT        *scaffold;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    if(scaffold->type != REAL_SCAFFOLD)
      continue;

    //if (scaffold->id != 14)
    //  continue;

    instrumentSCF  scf(scaffold);

    scf.analyze(lib);
    scf.report();
  }


  DestroyScaffoldGraph(ScaffoldGraph);

  return(0);
}
