
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

static char CM_ID[] = "$Id: dumpSingletons.c,v 1.16 2007-02-20 21:57:59 brianwalenz Exp $";

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
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"

#include "SYS_UIDcommon.h"
#include "SYS_UIDclient.h"



CDS_UID_t
getFragmentClear(int    iid,
                 int    reversecomplement,
                 char  *toprint) {

  static fragRecord  *fs = NULL;
  unsigned int  clr_bgn, clr_end;

  if (fs == NULL)
    fs = new_fragRecord();

  if (getFrag(ScaffoldGraph->gkpStore, iid, fs, FRAG_S_SEQ) != 0) {
    fprintf(stderr,"Couldn't get fragment from frgStore for iid %d\n", iid);
    assert(0);
  }
 
  clr_bgn = getFragRecordClearRegionBegin(fs, AS_READ_CLEAR_LATEST);
  clr_end = getFragRecordClearRegionEnd  (fs, AS_READ_CLEAR_LATEST);

  strcpy(toprint, getFragRecordSequence(fs) + clr_bgn);
  toprint[clr_end - clr_bgn] = 0;

  if (reversecomplement)
    Complement_Seq(toprint);

  return(getFragRecordUID(fs));
}



CDS_UID_t
getUID(int realUID) {
  static uint64        blockSize = 0;
  static int           UIDstart  = 1230000;
  static CDS_UID_t     interval_UID[4];

  CDS_UID_t            uid       = 0;;

  if (blockSize == 0) {
    blockSize = 300;
    set_start_uid(UIDstart); /* used if realUID == FALSE */
    get_uids(blockSize, interval_UID, realUID);
  }

  if (get_next_uid(&uid, realUID) != UID_CODE_OK) {
    get_uids(blockSize, interval_UID, realUID);
    if (get_next_uid(&uid, realUID) != UID_CODE_OK) {
      fprintf(stderr, "Could not get UID!\n");
      assert(0);
    }
  }

  return(uid);
}



int
main( int argc, char **argv) {
  int ckptNum           = NULLINDEX;
  int realUID           = 0;
  int makeMiniScaffolds = 1;

  GlobalData          = CreateGlobal_CGW();
  GlobalData->stderrc = stderr;
  GlobalData->timefp  = stderr;

  GlobalData->File_Name_Prefix[0] = 0;
  GlobalData->Gatekeeper_Store_Name[0] = 0;

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-p") == 0) {
      ckptNum = SetFileNamePrefix_CGW(GlobalData, argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->File_Name_Prefix, argv[++arg]);
    } else if (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->Gatekeeper_Store_Name, argv[++arg]);
    } else if (strcmp(argv[arg], "-n") == 0) {
      ckptNum = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-U") == 0) {
      realUID = 1;
    } else if (strcmp(argv[arg], "-S") == 0) {
      makeMiniScaffolds = 0;
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err = 1;
    }
    arg++;
  }

  if ((GlobalData->File_Name_Prefix[0]      == 0) ||
      (GlobalData->Gatekeeper_Store_Name[0] == 0)) {
    fprintf(stderr, "usage: %s [[-p prefix] | [-c name -g gkpstore -n ckptNum]] [-U] [-S]\n", argv[0]);
    fprintf(stderr, "  -p      Attempt to locate the last checkpoint in directory 7-CGW.\n");
    fprintf(stderr, "  -c      Look for checkpoints in 'name'\n");
    fprintf(stderr, "  -g      Path to gkpStore\n");
    fprintf(stderr, "  -n      Checkpoint number to load\n");
    fprintf(stderr, "  -U      Use real UIDs for miniscaffolds, otherwise, UIDs start at 1230000\n");
    fprintf(stderr, "  -S      Do NOT make mini scaffolds.\n");
    exit(1);
  }

  char *toprint1   = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LEN);
  char *toprint2   = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LEN);

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(GlobalData->File_Name_Prefix, ckptNum, FALSE);

  int ifrag;
  for (ifrag=0; ifrag < GetNumVA_CIFragT(ScaffoldGraph->CIFrags); ifrag++) {
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, ifrag);
    CIFragT *mate = NULL;

    assert(frag->cid != NULLINDEX);
    assert((frag->flags.bits.hasMate == 0) || (frag->mateOf != NULLINDEX));

    //  Fix for missing mates -- OBT used to not delete mate links, leaving
    //  dangling mates.  Somebody else seems to be doing this too.
    //
    if (frag->flags.bits.hasMate) {
      mate = GetCIFragT(ScaffoldGraph->CIFrags, frag->mateOf);
      if (mate == NULL)
        frag->flags.bits.hasMate = 0;
    }

    //  If this fragment is not chaff, we have nothing to do here.
    //
    if (GetGraphNode(ScaffoldGraph->CIGraph,frag->cid)->flags.bits.isChaff == 0)
      continue;

    //  Print a singleton if there is no mate, the mate isn't chaff,
    //  or we were told to not make miniscaffolds.
    //
    if ((mate == NULL) ||
        (mate->flags.bits.isChaff == 0) ||
        (makeMiniScaffolds == 0)) {
      CDS_UID_t  fUID = getFragmentClear(frag->iid, 0, toprint1);

      fprintf(stdout, ">"F_S64" /type=singleton\n%s\n",
             fUID, toprint1);
    } else if ((mate != NULL) &&
               (mate->flags.bits.isChaff == 1) &&
               (makeMiniScaffolds == 1) &&
               (frag->iid < mate->iid)) {
      CDS_UID_t  fUID = getFragmentClear(frag->iid, 0, toprint1);
      CDS_UID_t  mUID = getFragmentClear(mate->iid, 1, toprint2);

      //  make sure the following chain of Ns is divisible by three;
      //  the exact length is arbitrary but Doug Rusch points out that
      //  by making it divisible by 3, we can get lucky and maintain
      //  the phase of a protein ...  which helps in the
      //  auto-annotation of environmental samples

      fprintf(stdout, ">"F_S64" /type=mini_scaffold /frgs=("F_S64","F_S64")\n%sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN%s\n",
              getUID(realUID),
              fUID, mUID, toprint1, toprint2);
    }
  }

  exit(0);
}
