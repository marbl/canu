
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

const char *mainid = "$Id: ExamineSDB.c,v 1.2 2009-09-04 20:25:53 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"
#include "AS_MSG_pmesg.h"
#include "AS_SDB_SequenceDB.h"


int
main(int argc, char **argv) {
  char  *storeName         = NULL;
  int    storeVersion      = -1;

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-seqstore") == 0) {
      storeName = argv[++arg];

    } else if (strcmp(argv[arg], "-version") == 0) {
      storeVersion = atoi(argv[++arg]);

    } else {
      err++;
    }
    arg++;
  }
  if ((storeName    == NULL) ||
      (storeVersion == -1) ||
      (err)) {
    fprintf(stderr, "usage: %s -seqstore s.SeqStore -version n\n", argv[0]);
    exit(1);
  }

  fprintf(stderr, "sizeof(IntMultiPos)   = "F_S64"\n", sizeof(IntMultiPos));
  fprintf(stderr, "sizeof(IntUnitigPos)  = "F_S64"\n", sizeof(IntUnitigPos));
  exit(1);

  {
    char cmd[256];
    sprintf(cmd, "ls -l %s", storeName);
    system(cmd);
  }


  for (int32 v=0; v<storeVersion; v++) {
    tSequenceDB  *sequenceDB = openSequenceDB(storeName, FALSE, v);
    int32         nutg       = numUnitigsInSequenceDB(sequenceDB);
    int32         nctg       = numContigsInSequenceDB(sequenceDB);

    fprintf(stderr, "version: %3d  nutg: %6d  nctg: %6d\n", v, nutg, nctg);

    for (int32 i=0; i<nutg; i++) {
      tMARecord *tma = GettMARecord(sequenceDB->Unitigs, i);

      if ((tma->storeID > 0) && (tma->storeID == v))
        fprintf(stderr, "TMAU: i=%d storeID=%d  multiAlignID=%d  isDeleted=%d  offset="F_S64"\n",
                i, tma->storeID, tma->multiAlignID, tma->isDeleted, tma->offset);

      loadMultiAlignTFromSequenceDB(sequenceDB, i, TRUE);
    }

    for (int32 i=0; i<nctg; i++) {
      tMARecord *tma = GettMARecord(sequenceDB->Contigs, i);

      if ((tma->storeID > 0) && (tma->storeID == v))
        fprintf(stderr, "TMAC: i=%d storeID=%d  multiAlignID=%d  isDeleted=%d  offset="F_S64"\n",
                i, tma->storeID, tma->multiAlignID, tma->isDeleted, tma->offset);

      loadMultiAlignTFromSequenceDB(sequenceDB, i, FALSE);
    }

    deleteSequenceDB(sequenceDB);
  }
}

