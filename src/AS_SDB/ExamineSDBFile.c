
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

const char *mainid = "$Id: ExamineSDBFile.c,v 1.1 2009-09-04 20:25:53 brianwalenz Exp $";

//  Violates the SDB interface and examines a single SDB version directly.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"
#include "AS_MSG_pmesg.h"
#include "AS_SDB_SequenceDB.h"







tSequenceDB *
openSequenceDBFile(char *path, int readWrite, int revision){
  FILE *F;
  char N[FILENAME_MAX + 30];
  FILE *datafp;
  tSequenceDB *db = (tSequenceDB *)safe_calloc(1, sizeof(tSequenceDB));
  int i;

  db->path              = strdup(path);
  db->currentRevision   = revision;
  db->Unitigs           = NULL;
  db->Contigs           = NULL;
  db->UnitigStore       = CreateMultiAlignStoreT();
  db->ContigStore       = CreateMultiAlignStoreT();
  db->dataFileLen       = 0;
  db->dataFileMax       = 1024;
  db->dataFile          = (FILE **)safe_calloc(db->dataFileMax, sizeof(FILE *));

  sprintf(N,"%s/seqDB.v%03d.utg", db->path, revision);
  errno = 0;
  F = fopen(N,"r");
  if (errno)
    fprintf(stderr, "OpenSequenceDB()-- Failed to open '%s' for reading: %s\n",
            N, strerror(errno)), exit(1);
  db->Unitigs = CreateFromFileVA_tMARecord(F);
  fclose(F);

  sprintf(N,"%s/seqDB.v%03d.ctg", db->path, revision);
  errno = 0;
  F = fopen(N,"r");
  if (errno)
    fprintf(stderr, "OpenSequenceDB()-- Failed to open '%s' for reading: %s\n",
            N, strerror(errno)), exit(1);
  db->Contigs = CreateFromFileVA_tMARecord(F);
  fclose(F);

  sprintf(N, "%s/seqDB.v%03d.dat", db->path, revision);
  errno = 0;
  F = fopen(N, "r");
  if (errno)
    fprintf(stderr, "OpenSequenceDB()-- Failed to open '%s' for reading: %s\n",
            N, strerror(errno)), exit(1);
  db->dataFile[db->dataFileLen++] = F;

  return db;
}








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


  fprintf(stderr, "sizeof(IntMultiPos)  = %d\n", sizeof(IntMultiPos));
  fprintf(stderr, "sizeof(IntUnitigPos)  = %d\n", sizeof(IntUnitigPos));


  tSequenceDB  *sequenceDB = openSequenceDBFile(storeName, FALSE, storeVersion);
  int32         nutg       = numUnitigsInSequenceDB(sequenceDB);
  int32         nctg       = numContigsInSequenceDB(sequenceDB);


  fprintf(stderr, "version: %3d  nutg: %6d  nctg: %6d\n", storeVersion, nutg, nctg);


  for (int32 i=0; i<nutg; i++) {
    tMARecord *tma = GettMARecord(sequenceDB->Unitigs, i);

    if ((i >= 1086000) && (i <= 1086010))
    fprintf(stderr, "TMAU: i=%d storeID=%d  multiAlignID=%d  isDeleted=%d  offset="F_S64"\n",
            i, tma->storeID, tma->multiAlignID, tma->isDeleted, tma->offset);

    if (i != 1086001)
      loadMultiAlignTFromSequenceDB(sequenceDB, i, TRUE);
  }


  for (int32 i=0; i<nctg; i++) {
    tMARecord *tma = GettMARecord(sequenceDB->Contigs, i);

    if ((i >= 1086000) && (i <= 1086010))
    fprintf(stderr, "TMAC: i=%d storeID=%d  multiAlignID=%d  isDeleted=%d  offset="F_S64"\n",
            i, tma->storeID, tma->multiAlignID, tma->isDeleted, tma->offset);

    loadMultiAlignTFromSequenceDB(sequenceDB, i, FALSE);
  }

  
  deleteSequenceDB(sequenceDB);
}
