
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

/* $Id: AS_GKP_include.h,v 1.40 2008-06-12 03:41:29 brianwalenz Exp $ */

#ifndef AS_GKP_INCLUDE_H
#define AS_GKP_INCLUDE_H

#include <stdio.h>
#include <errno.h>
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Var.h"
#include "AS_MSG_pmesg.h"

#define GATEKEEPER_MAX_ERROR_RATE        0.025 
#define GATEKEEPER_QV_WINDOW_WIDTH      50
#define GATEKEEPER_QV_WINDOW_THRESH      0.03

#define AS_ASSEMBLER_GRANDE  ((int)'A')
#define AS_ASSEMBLER_OBT     ((int)'T')

extern GateKeeperStore  *gkpStore;
extern FILE             *errorFP;

int
Check_BatchMesg(BatchMesg           *bat_mesg);

int
Check_DistanceMesg(DistanceMesg     *dst_mesg,
                   int                believeInputStdDev);

int
Check_LibraryMesg(LibraryMesg       *dst_mesg,
                  int                believeInputStdDev);

int
Check_FragMesg(FragMesg            *frg_mesg,  
               int                   assembler);

int
Check_LinkMesg(LinkMesg             *lkg_mesg);


int
Load_SFF(FILE *sff, int searchForLinker);


void
dumpGateKeeperInfo(char       *gkpStoreName);

void
dumpGateKeeperBatches(char       *gkpStoreName,
                      AS_IID      begIID,
                      AS_IID      endIID,
                      char       *iidToDump, 
                      int         asTable);

void
dumpGateKeeperLibraries(char       *gkpStoreName,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump, 
                        int         asTable);

void
dumpGateKeeperFragments(char       *gkpStoreName,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump, 
                        int         dumpWithSequence,
                        int         dumpClear,
                        int         asTable);

void
dumpGateKeeperAsFasta(char       *gkpStoreName,
                      AS_IID      begIID,
                      AS_IID      endIID,
                      char       *iidToDump, 
                      int         dumpFastaAllReads,
                      int         dumpFastaClear,
                      int         dumpFastaQuality);

void
dumpGateKeeperAsFRG(char       *gkpStoreName,
                    int         dumpFormat,
                    AS_IID      begIID,
                    AS_IID      endIID,
                    char       *iidToDump,
                    int         doNotFixMates,
                    int         dumpFRGClear);

void
dumpGateKeeperAsNewbler(char       *gkpStoreName,
                        char       *prefix,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         doNotFixMates,
                        int         dumpFRGClear);


int
Build_Partition(char      *gatekeeperName,
                char      *partitionFile,
                int32      flags);
             
int
rebuildMap(char *gkpStoreName);

void
rearrangeStore(char *uidFile,
               char *gkpStore,
               char *newStore);

void
updateVectorClear(char *vectorClearFile, char *gkpStoreName);

void
editStore(char *editsFileName, char *gkpStoreName, int update);


//  Error handling

void
AS_GKP_reportError(int error, ...);

int
AS_GKP_summarizeErrors(void);


#define AS_GKP_BAT_ZERO_UID              1
#define AS_GKP_BAT_EXISTS                2

#define AS_GKP_FRG_INVALID_CHAR_SEQ      3
#define AS_GKP_FRG_INVALID_CHAR_QLT      4
#define AS_GKP_FRG_INVALID_LENGTH        5
#define AS_GKP_FRG_ZERO_UID              6
#define AS_GKP_FRG_EXISTS                7
#define AS_GKP_FRG_SEQ_TOO_LONG          8
#define AS_GKP_FRG_SEQ_TOO_SHORT         9
#define AS_GKP_FRG_CLR_BGN              10
#define AS_GKP_FRG_CLR_END              11
#define AS_GKP_FRG_CLR_TOO_SHORT        12
#define AS_GKP_FRG_UNKNOWN_LIB          13
#define AS_GKP_FRG_LOADED_DELETED       14
#define AS_GKP_FRG_DOESNT_EXIST         15
#define AS_GKP_FRG_HAS_MATE             16
#define AS_GKP_FRG_UNKNOWN_ACTION       17

#define AS_GKP_LIB_ILLEGAL_MEAN_STDDEV  18
#define AS_GKP_LIB_INVALID_MEAN         19
#define AS_GKP_LIB_INVALID_STDDEV       20
#define AS_GKP_LIB_STDDEV_TOO_BIG       21
#define AS_GKP_LIB_STDDEV_TOO_SMALL     22
#define AS_GKP_LIB_EXISTS               23
#define AS_GKP_LIB_ZERO_UID             24
#define AS_GKP_LIB_DOESNT_EXIST_UPDATE  25
#define AS_GKP_LIB_UNKNOWN_ACTION       28

#define AS_GKP_LKG_SELF_LINK            29
#define AS_GKP_LKG_UNSUPPORTED_TYPE     30
#define AS_GKP_LKG_FRG_DOESNT_EXIST     31
#define AS_GKP_LKG_FRG_DELETED          32
#define AS_GKP_LKG_ALREADY_MATED        33
#define AS_GKP_LKG_LIB_DOESNT_EXIST     36
#define AS_GKP_LKG_DIFFERENT_LIB        37
#define AS_GKP_LKG_UNKNOWN_ACTION       38

#define AS_GKP_SFF_ALREADY_EXISTS       39
#define AS_GKP_SFF_UID_ERROR            40
#define AS_GKP_SFF_TOO_SHORT            41
#define AS_GKP_SFF_TOO_LONG             42
#define AS_GKP_SFF_N                    43

#define AS_GKP_UNKNOWN_MESSAGE          50

#endif
