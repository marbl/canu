
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
static char CM_ID[] = "$Id: testSplitChunks.c,v 1.3 2005-03-22 19:04:22 jason_miller Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "Globals_CGW.h"
#include "SplitChunks_CGW.h"

int main(int argc, char ** argv)
{
  int32  checkPointNumber = 0;
  ScaffoldGraphT * graph = NULL;
  
  if( argc != 5 )
  {
    fprintf(stderr,
            "Usage: %s <checkPointName> <checkPointNumber> "
            "<gkpStoreName> <fragStoreName>\n",
            argv[0] );
    fprintf(stderr,"\ttests unitig splitting code.\n" );
    return 1;
  }

  GlobalData = CreateGlobal_CGW();
  GlobalData->stderrc = stderr;
  strcpy( GlobalData->File_Name_Prefix, argv[1] );
  checkPointNumber = atoi( argv[2] );
  strcpy( GlobalData->Gatekeeper_Store_Name, argv[3] );
  strcpy( GlobalData->Frag_Store_Name, argv[4] );

  graph =
    LoadScaffoldGraphFromCheckpoint( GlobalData->File_Name_Prefix,
                                     checkPointNumber,
                                     FALSE );
  if( graph == NULL )
  {
    fprintf( stderr, "Failed to load scaffold graph.\n" );
    return 1;
  }

  SplitInputUnitigs(graph);
  DeleteGlobal_CGW(GlobalData);

  return 0;
}
