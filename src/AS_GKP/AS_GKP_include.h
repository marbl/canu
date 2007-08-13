
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

/* $Id: AS_GKP_include.h,v 1.29 2007-08-13 05:46:41 brianwalenz Exp $ */

#ifndef AS_GKP_INCLUDE_H
#define AS_GKP_INCLUDE_H

#include <stdio.h>
#include <errno.h>
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Var.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_SequenceBucket.h"

#define GATEKEEPER_MAX_ERROR_RATE        0.025 
#define GATEKEEPER_QV_WINDOW_WIDTH      50
#define GATEKEEPER_QV_WINDOW_THRESH      0.03

#define AS_ASSEMBLER_GRANDE  ((int)'A')
#define AS_ASSEMBLER_OBT     ((int)'T')

#define GATEKEEPER_SUCCESS 0
#define GATEKEEPER_FAILURE 1

extern GateKeeperStore  *gkpStore;
extern FILE             *errorFP;

int
Check_BatchMesg(BatchMesg           *bat_mesg);

int
Check_DistanceMesg(DistanceMesg     *dst_mesg);

int
Check_LibraryMesg(LibraryMesg       *dst_mesg);

int
Check_FragMesg(FragMesg            *frg_mesg,  
               int                   check_qvs,
               int                   assembler);

int
Check_LinkMesg(LinkMesg             *lkg_mesg);



void
dumpGateKeeperInfo(char       *gkpStoreName);

void
dumpGateKeeperBatches(char       *gkpStoreName,
                      CDS_IID_t   begIID,
                      CDS_IID_t   endIID,
                      char       *iidToDump, 
                      int         asTable);

void
dumpGateKeeperLibraries(char       *gkpStoreName,
                        CDS_IID_t   begIID,
                        CDS_IID_t   endIID,
                        char       *iidToDump, 
                        int         asTable);

void
dumpGateKeeperFragments(char       *gkpStoreName,
                        CDS_IID_t   begIID,
                        CDS_IID_t   endIID,
                        char       *iidToDump, 
                        int         dumpWithSequence,
                        int         dumpClear,
                        int         asTable);

void
dumpGateKeeperAsFasta(char       *gkpStoreName,
                      CDS_IID_t   begIID,
                      CDS_IID_t   endIID,
                      char       *iidToDump, 
                      int         dumpFastaAllReads,
                      int         dumpFastaClear,
                      int         dumpFastaQuality);

void
dumpGateKeeperAsFRG(char       *gkpStoreName,
                    int         dumpFormat,
                    CDS_IID_t   begIID,
                    CDS_IID_t   endIID,
                    char       *iidToDump,
                    int         doNotFixMates,
                    int         dumpFRGClear);



int
Build_Partition(char      *gatekeeperName,
                char      *partitionFile,
                int32      flags);
             
int
rebuildMap(char *hashFileName,
           char *gkpStoreName);

void
rearrangeStore(char *uidFile,
               char *gkpStore,
               char *newStore);

void
updateVectorClear(char *vectorClearFile, char *gkpStoreName);

void
editStore(char *editsFileName, char *gkpStoreName, int update);

#endif
