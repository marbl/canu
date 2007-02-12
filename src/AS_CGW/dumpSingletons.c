
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

static char CM_ID[] = "$Id: dumpSingletons.c,v 1.13 2007-02-12 22:16:56 brianwalenz Exp $";

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
getFragmentClear(int iid,
                 int reversecomplement,
                 int  *alloclen,
                 char **seq,
                 char **qul,
                 char **toprint) {

  static ReadStructp               fs = NULL;
  static GateKeeperFragmentRecord  gs;

#warning someone please rewrite me!

  //  we don't need to get both getFrag() and the gatekeeper frag; they're the same.

  if (fs == NULL)
    fs = new_ReadStruct();

  if (getFrag(ScaffoldGraph->gkpStore,
              iid,
              fs, FRAG_S_ALL) != 0) {
    fprintf(stderr,"Couldn't get fragment from frgStore for iid %d\n", iid);
    assert(0);
  }

  if (getGateKeeperFragmentStore(ScaffoldGraph->gkpStore->frg,
                                 iid,
                                 &gs) != 0) {
    fprintf(stderr,"Couldn't get fragment from gkpStore for iid %d\n", iid);
    assert(0);
  }

  unsigned int  clr_bgn, clr_end;
  
  getClearRegion_ReadStruct(fs, &clr_bgn, &clr_end, READSTRUCT_LATEST);
  while (getSequence_ReadStruct(fs, *seq, *qul, *alloclen) != 0) {
    *alloclen *= 2;
    *seq       = (char*)safe_realloc(*seq,     *alloclen * sizeof(char));
    *qul       = (char*)safe_realloc(*qul,     *alloclen * sizeof(char));
    *toprint   = (char*)safe_realloc(*toprint, *alloclen * sizeof(char));
  }
  
  strcpy(*toprint, *seq + clr_bgn);
  (*toprint)[clr_end - clr_bgn] = 0;

  if (reversecomplement)
    Complement_Seq(*toprint);

  return(gs.UID);
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

  int   alloclen1  = 2048;
  char *seq1       = (char *)safe_malloc(sizeof(char) * alloclen1);
  char *qul1       = (char *)safe_malloc(sizeof(char) * alloclen1);
  char *toprint1   = (char *)safe_malloc(sizeof(char) * alloclen1);

  int   alloclen2  = 2048;
  char *seq2       = (char *)safe_malloc(sizeof(char) * alloclen2);
  char *qul2       = (char *)safe_malloc(sizeof(char) * alloclen2);
  char *toprint2   = (char *)safe_malloc(sizeof(char) * alloclen2);


  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(GlobalData->File_Name_Prefix, ckptNum, FALSE);

  int ifrag;
  for (ifrag=0; ifrag < GetNumVA_CIFragT(ScaffoldGraph->CIFrags); ifrag++) {
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, ifrag);
    CIFragT *mate = NULL;

    assert(frag->cid != NULLINDEX);
    assert((frag->numLinks == 0) || (frag->numLinks == 1));

    //  Fix for missing mates -- OBT used to not delete mate links, leaving
    //  dangling mates.  Somebody else seems to be doing this too.
    //
    if (frag->numLinks > 0) {
      mate = GetCIFragT(ScaffoldGraph->CIFrags, frag->mateOf);
      if (mate == NULL)
        frag->numLinks = 0;
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
      CDS_UID_t  fUID = getFragmentClear(frag->iid, 0, &alloclen1, &seq1, &qul1, &toprint1);

      fprintf(stdout, ">"F_S64" /type=singleton\n%s\n",
             fUID, toprint1);
    } else if ((mate != NULL) &&
               (mate->flags.bits.isChaff == 1) &&
               (makeMiniScaffolds == 1) &&
               (frag->iid < mate->iid)) {
      CDS_UID_t  fUID = getFragmentClear(frag->iid, 0, &alloclen1, &seq1, &qul1, &toprint1);
      CDS_UID_t  mUID = getFragmentClear(mate->iid, 1, &alloclen2, &seq2, &qul2, &toprint2);

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
