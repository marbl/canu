
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
static char CM_ID[] = "$Id: Checkpoint.c,v 1.4 2005-03-22 19:48:35 jason_miller Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"

int main(int argc, char *argv[]){
  char *checkPointName = argv[1];
  int32 checkPointNumber = atoi(argv[2]);
  ScaffoldGraphT *scaffoldGraph = NULL;


  if(argc != 3){
    fprintf(stderr,"*usage: checkpoint <checkPointName> <checkPointNumber>\n");    
    fprintf(stderr,"*       e.g.:  checkpoint dros 3, loads checkpoint file dros.ckp.3\n");
  }

  GlobalData = CreateGlobal_CGW();
  strcpy( GlobalData->Frag_Store_Name, "");
  strcpy( GlobalData->Gatekeeper_Store_Name, "");


  fprintf(stderr,"* name %s number %d\n", checkPointName, checkPointNumber);
  scaffoldGraph = LoadScaffoldGraphFromCheckpoint(checkPointName,checkPointNumber, FALSE);

  ReportMemorySize(ScaffoldGraph, stderr);

  fprintf(stderr,"* Deleting Scaffold Graph *\n");
  DestroyScaffoldGraph(scaffoldGraph);
  return 0;
}
