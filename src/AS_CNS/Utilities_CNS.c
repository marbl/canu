
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
/**************************************************************
 * AS_CNS/Utilities_CNS.c
 *
 * Utility functions for the CNS (consensus) subsystem 
 * of the Celera WGS assembler.
 *
 **************************************************************/

/*********************************************************************
 $Id: Utilities_CNS.c,v 1.1.1.1 2004-04-14 13:51:22 catmandew Exp $
 *********************************************************************/
// Operating System includes:
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

/***********************************

// Celera assembler includes:

// Consensus includes:
#include "Globals_CNS.h"
#include "Utilities_CNS.h"


//------------------------------------------------------------
static int statusToReportOnTermination = 1;



//-------------------------------------------------
// Set the exit status to report should consensus crash.
// Jason, 7/01.
//-------------------------------------------------
// This function is public, visible outside CNS.
void CNS_setExitStatus (int exitStatus) {
    // This function is protected, visible only within CNS.
    statusToReportOnTermination = exitStatus;
}

//-------------------------------------------------
// Set parameters to their default values.
// Jason, 7/01.
//-------------------------------------------------
void CNS_setConsensusParametersToDefault() {
     CNS_setConsensusParametersIndividually(
					0, // use sdb
					0,  // use part sdb
					0.02,  // seq err rate
					0.0003,  // snp rate
					1,  // haplotypes
					0,  // use public
					0,  // favor public
					1,  // debugging
					0); // partitioned
}

//-------------------------------------------------
// Set parameters used by the consensus algorithm.
// See setConsensusParametersToDefault().
// Jason, 7/01.
//-------------------------------------------------
void CNS_setConsensusParametersIndividually(
    int useSDB,
    int usePartSDB,
    float priorSequencingErrorRate,
    float priorSNPRate,
    int numHaplotypes,
    int basecallingUsesPublicData,
    int basecallingFavorsPublicData,
    int outputDebugging,
    int isPartitioned) {

    USE_SDB=useSDB;

    USE_SDB_PART=usePartSDB;

    // Used to calculate '-' probability
    CNS_SEQUENCING_ERROR_EST = priorSequencingErrorRate;

    // Used to calculate BIAS
    CNS_SNP_RATE   = priorSNPRate;

    // Used to calculate BIAS
    CNS_HAPLOTYPES = numHaplotypes;

    // Used to direct basecalling to include public data
    CNS_USE_PUBLIC = basecallingUsesPublicData;

    // Used to direct basecalling to favor public data
    CNS_CALL_PUBLIC =basecallingFavorsPublicData;

    debug_out = outputDebugging;
    partitioned=isPartitioned;    
}



//----------------------------------------------------------
// Consensus utility function:
// Terminate execution cleanly after some error.
// See setStatusToReportOnTermination().
//----------------------------------------------------------
// Moved here from Globals_CNS.h.
void CleanExit(char *mesg, int lineno, int rc) {
  char command[100+FILENAME_MAX];
  int i;
  fprintf(stderr,"%s at line: %d, rc: %d\n",mesg,lineno,rc);
  if( cnsout != NULL && ! std_output ){
    fclose(cnsout);
    if ( debug_out ) {
      sprintf(command,"mv -f %s %s.dbg",OutputFileName,OutputFileName);
    } else {
      sprintf(command,"rm -f %s",OutputFileName);
    }
    fprintf(stderr,"%s\n",command);
    system(command);
    sprintf(command,"touch %s",OutputFileName);
    system(command);
  }
  if( cnslog != NULL && ! std_error_log) {
    fclose(cnslog);
    if ( debug_out ) {
      sprintf(command,"mv -f %s %s.dbg",LogFileName,LogFileName);
    } else {
      sprintf(command,"rm -f %s",LogFileName);
    }
    fprintf(stderr,"%s\n",command);
    system(command);
  }
  exit(statusToReportOnTermination);
  //exit(rc);
}



**************************************/
