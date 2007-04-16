
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
static char CM_ID[] = "$Id: examineFragStore.c,v 1.10 2007-04-16 17:36:32 brianwalenz Exp $";


/*********************************************************************
 * Module:  based on AS_CGW_LoadCheckpoint
 * Description:
 *    For use with debugger to query values in a checkpoint
 * 
 *    Reference: 
 *
 *    Command Line Interface:
 *        $ loadcgw checkpointPath
 *
 *       CGBInputFiles: The file with new IUM,OUM, etc records to process. 
 *
 *       Checkpoint File: File named <outputName>.ckp.n
 *
 * 
 *********************************************************************/
//#define DEBUG 1
//#define DEBUG_BUCIS 1
//#define DEBUG_MERGE_SCAF 1
#define FRAG_POSS_SIZE 200000

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "FbacREZ.h"
#include "PublicAPI_CNS.h"

int printClearRanges( int fragIid );

int debug = 0;

static ReadStructp fsread = NULL;
GateKeeperStore *gkpStore;

int examineClearRanges( int fragIid )
{
  uint32 clr_bgn_orig, clr_bgn_latest, clr_end_orig, clr_end_latest;
  int setStatus = 0;
  
  if( fsread == NULL ) fsread = new_ReadStruct();

  if ( fragIid != -1)
    {
      getFrag( gkpStore, fragIid, fsread, FRAG_S_INF);	

      getClearRegion_ReadStruct( fsread, &clr_bgn_orig, &clr_end_orig, READSTRUCT_ORIGINAL);
      getClearRegion_ReadStruct( fsread, &clr_bgn_latest, &clr_end_latest, READSTRUCT_LATEST);

      if (( clr_bgn_orig != clr_bgn_latest ) || ( clr_end_orig != clr_end_latest ))
        printClearRanges( fragIid );	
    }
  return ( setStatus ); 
}

int printClearRanges( CDS_CID_t fragIid )
{
  uint32 clr_bgn, clr_end;
  int setStatus = 0;
  
  if( fsread == NULL ) fsread = new_ReadStruct();

  if ( fragIid != -1)
    {
      getFrag( gkpStore, fragIid, fsread, FRAG_S_ALL);

      getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_ORIGINAL);
      fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_ORG clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
               fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

      getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_OVL);
      fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_OVL clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
               fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

      getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CNS);
      fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_CNS clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
               fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

      getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_CGW);
      fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_CGW clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
               fragIid, clr_bgn, clr_end, clr_end - clr_bgn);

      getClearRegion_ReadStruct( fsread, &clr_bgn, &clr_end, READSTRUCT_LATEST);
      fprintf( stderr, "printClearRanges, frag %8" F_CIDP " READSTRUCT_LAT clr_bgn: %5" F_U32P ", clr_end %5" F_U32P ", len: %4" F_U32P "\n",
               fragIid, clr_bgn, clr_end, clr_end - clr_bgn);
    }

  return ( setStatus ); 
}

int main( int argc, char *argv[])
{
  Global_CGW *data;
  int ch, errflg=0;
  time_t t1;
  int64 i;
  
  if (argc != 2) {
    fprintf(stderr, "usage: %s gkpStore\n", argv[0]);
    exit(1);
  }

  strcpy( Gatekeeper_Store_Name, argv[1]);

  t1 = time(0);
  fprintf( stderr, "====> Starting at %s\n", ctime(&t1));

  gkpStore = openGateKeeperStore( Gatekeeper_Store_Name, FALSE);
  if ( gkpStore == NULLSTOREHANDLE )
    assert(0);

  fprintf( stderr, "examining frags from " F_S64 " to " F_S64 "\n",
           getFirstElemFragStore(gkpStore), getLastElemFragStore(gkpStore));
  
  for ( i = getFirstElemFragStore(gkpStore); i <= getLastElemFragStore(gkpStore); i++)
    {
      examineClearRanges( i );
      if ( (i + 1) % 10000 == 0)
        fprintf( stderr, "examined frag " F_S64 "\n", i);
    }
  
  closeGateKeeperStore( gkpStore );
  return 0;
}
