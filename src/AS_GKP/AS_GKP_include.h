
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

/* $Id: AS_GKP_include.h,v 1.61 2011-02-11 05:48:13 brianwalenz Exp $ */

#ifndef AS_GKP_INCLUDE_H
#define AS_GKP_INCLUDE_H

static const char *rcsid_AS_GKP_INCLUDE_H = "$Id: AS_GKP_include.h,v 1.61 2011-02-11 05:48:13 brianwalenz Exp $";

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

extern gkStore     *gkpStore;
extern gkFragment  *gkFrag1;
extern gkFragment  *gkFrag2;
extern FILE        *errorFP;
extern char         illuminaUIDmapName[FILENAME_MAX];
extern FILE        *illuminaUIDmap;

int
Check_DistanceMesg(DistanceMesg     *dst_mesg,
                   int                believeInputStdDev);

void
checkLibraryForIlluminaPointers(LibraryMesg *lib_mesg);

int
Check_LibraryMesg(LibraryMesg       *dst_mesg,
                  int                believeInputStdDev);

int
Check_FragMesg(FragMesg            *frg_mesg,
               int                   assembler);

int
Check_LinkMesg(LinkMesg             *lkg_mesg);

int
Check_PlacementMesg(PlacementMesg   *plc_mesg);

void
dumpGateKeeperInfo(char       *gkpStoreName,
                   int         asTable,
                   int         withoutUIDs);

void
dumpGateKeeperLibraries(char       *gkpStoreName,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         asTable,
                        int         withoutUIDs);

void
dumpGateKeeperFragments(char       *gkpStoreName,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         dumpWithSequence,
                        int         dumpClear,
                        int         asTable,
                        int         withoutUIDs);

void
dumpGateKeeperAsFasta(char       *gkpStoreName,
                      char       *prefix,
                      AS_IID      begIID,
                      AS_IID      endIID,
                      char       *iidToDump,
                      int         doNotFixMates,
                      int         dumpAllReads,
                      int         dumpAllBases,
                      int         dumpClear);

void
dumpGateKeeperAsFRG(char       *gkpStoreName,
                    int         dumpFormat,
                    AS_IID      begIID,
                    AS_IID      endIID,
                    char       *iidToDump,
                    int         doNotFixMates,
                    int         dumpAllReads,
                    int         dumpFRGClear);

void
dumpGateKeeperAsNewbler(char       *gkpStoreName,
                        char       *prefix,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         doNotFixMates,
                        int         dumpAllReads,
                        int         dumpAllBases,
                        int         dumpClear);

void
dumpGateKeeperAsFastQ(char       *gkpStoreName,
                      char       *prefix,
                      AS_IID      begIID,
                      AS_IID      endIID,
                      char       *iidToDump,
                      int         doNotFixMates,
                      int         dumpAllReads,
                      int         dumpAllBases,
                      int         dumpClear);

int
dumpGateKeeperIsFeatureSet(char       *gkpStoreName,
                           AS_IID     featureLibIID,
                           char      *featureName);

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
revertClearRange(char *clearRegionName, char *gkpStoreName);

void
editStore(char *editsFileName, char *gkpStoreName, int update);



//  Error handling

void
AS_GKP_reportError(int error, ...);

int
AS_GKP_summarizeErrors(void);


#define AS_GKP_BAT_ZERO_UID              1
#define AS_GKP_BAT_EXISTS                2

#define AS_GKP_FRG_INVALID_CHAR_SEQ     10
#define AS_GKP_FRG_INVALID_CHAR_QLT     11
#define AS_GKP_FRG_INVALID_LENGTH       12
#define AS_GKP_FRG_ZERO_UID             13
#define AS_GKP_FRG_EXISTS               14
#define AS_GKP_FRG_SEQ_TOO_LONG         15
#define AS_GKP_FRG_SEQ_TOO_SHORT        16
#define AS_GKP_FRG_CLR_BGN              17
#define AS_GKP_FRG_CLR_END              18
#define AS_GKP_FRG_CLR_TOO_SHORT        19
#define AS_GKP_FRG_CLR_INVALID          20
#define AS_GKP_FRG_UNKNOWN_LIB          21
#define AS_GKP_FRG_LOADED_DELETED       22
#define AS_GKP_FRG_DOESNT_EXIST         23
#define AS_GKP_FRG_HAS_MATE             24
#define AS_GKP_FRG_UNKNOWN_ACTION       25

#define AS_GKP_ILL_NOT_SEQ_START_LINE   30
#define AS_GKP_ILL_NOT_QLT_START_LINE   31
#define AS_GKP_ILL_SEQ_QLT_NAME_DIFFER  32
#define AS_GKP_ILL_SEQ_QLT_LEN_DIFFER   33
#define AS_GKP_ILL_SEQ_TOO_LONG         34
#define AS_GKP_ILL_CANT_OPEN_INPUT      35
#define AS_GKP_ILL_BAD_QV               36

#define AS_GKP_LIB_ILLEGAL_MEAN_STDDEV  50
#define AS_GKP_LIB_INVALID_MEAN         51
#define AS_GKP_LIB_INVALID_STDDEV       52
#define AS_GKP_LIB_STDDEV_TOO_BIG       53
#define AS_GKP_LIB_STDDEV_TOO_SMALL     54
#define AS_GKP_LIB_EXISTS               55
#define AS_GKP_LIB_ZERO_UID             56
#define AS_GKP_LIB_DOESNT_EXIST_UPDATE  57
#define AS_GKP_LIB_UNKNOWN_ACTION       58

#define AS_GKP_LKG_SELF_LINK            70
#define AS_GKP_LKG_UNSUPPORTED_TYPE     71
#define AS_GKP_LKG_FRG_DOESNT_EXIST     72
#define AS_GKP_LKG_FRG_DELETED          73
#define AS_GKP_LKG_ALREADY_MATED        74
#define AS_GKP_LKG_LIB_DOESNT_EXIST     75
#define AS_GKP_LKG_DIFFERENT_LIB        76
#define AS_GKP_LKG_UNMATED_LIB          77
#define AS_GKP_LKG_DIFFERENT_ORIENT     78
#define AS_GKP_LKG_UNKNOWN_ACTION       79

#define AS_GKP_SFF_ALREADY_EXISTS       90
#define AS_GKP_SFF_UID_ERROR            91
#define AS_GKP_SFF_TOO_SHORT            92
#define AS_GKP_SFF_TOO_LONG             93
#define AS_GKP_SFF_N                    94

#define AS_GKP_PLC_SAME_CONSTRAINT      100
#define AS_GKP_PLC_SELF_CONSTRAINT      101
#define AS_GKP_PLC_FRG_DOESNT_EXIST     102
#define AS_GKP_PLC_FRG_DELETED          103
#define AS_GKP_PLC_ALREADY_CONSTRAINED  104
#define AS_GKP_PLC_UNKNOWN_ACTION       105

#define AS_GKP_UNKNOWN_MESSAGE          120

#endif
