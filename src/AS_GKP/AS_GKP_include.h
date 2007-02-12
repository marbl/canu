
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

/* $Id: AS_GKP_include.h,v 1.13 2007-02-12 22:16:57 brianwalenz Exp $ */

#ifndef AS_GKP_INCLUDE_H
#define AS_GKP_INCLUDE_H

#define GATEKEEPER_SUCCESS 0
#define GATEKEEPER_WARNING 1
#define GATEKEEPER_FAILURE 2

#include <stdio.h>
#include <errno.h>
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_PHash.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_SequenceBucket.h"

#define GATEKEEPER_SCREENER_MIN_LENGTH  40
#define GATEKEEPER_MAX_ERROR_RATE        0.025 
#define GATEKEEPER_QV_WINDOW_WIDTH      50
#define GATEKEEPER_QV_WINDOW_THRESH      0.03

#define AS_ASSEMBLER_GRANDE  ((int)'A')
#define AS_ASSEMBLER_OBT     ((int)'T')

extern GateKeeperStore  *gkpStore;

int
Check_BatchMesg(BatchMesg           *bat_mesg,
                int                 *currentBatchID,
                time_t               currentTime,
                int                  verbose);

int
Check_LibraryMesg(DistanceMesg      *dst_mesg,
                  CDS_CID_t          currentBatchID,
                  int                verbose);

int
Check_FragMesg(FragMesg            *frg_mesg,  
               int                   check_qvs,
               int32                 batchID,
               time_t                currentTime,
               int                   assembler,
               int                   verbose);

int
Check_LinkMesg(LinkMesg             *lkg_mesg,
               CDS_CID_t             batchID,
               time_t                currentTime,
               int                   verbose);


//
//  Gatekeeper Errors
//

typedef enum GKPErrorType_tag{
  GKPError_FirstMessageBAT = 0,
  GKPError_BadUniqueBAT,
  GKPError_BadUniqueFRG,
  GKPError_BadUniqueLIB,
  GKPError_MissingFRG,
  GKPError_MissingLIB,
  GKPError_DeleteFRG,
  GKPError_DeleteLIB,
  GKPError_DeleteLNK,
  GKPError_Time,
  GKPError_Action,
  GKPError_Scalar,
  GKPError_FRGSequence,
  GKPError_FRGQuality,
  GKPError_FRGLength,
  GKPError_FRGClrRange,
  GKPError_FRGLocalPos,
  GKPError_FRGQualityWindow,
  GKPError_FRGQualityGlobal,
  GKPError_FRGQualityTail,
  GKPError_LNKFragLibMismatch,
  GKPError_LNKOneLink,
  GKPError_DSTValues,
  GKPError_MAX
} GKPErrorType;

#define MAX_GKPERROR ((int)(GKPError_MAX )-1)

void printGKPError(FILE *fout, GKPErrorType type);
void printAllGKPErrors(FILE *fout);

int32 CheckNmerProbabilities(FILE *, double threshhold);

#endif
