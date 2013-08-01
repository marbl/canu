
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 20008, J. Craig Venter Institute.
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

const char *mainid = "$Id$";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "AS_global.H"
#include "AS_UTL_Var.H"
#include "AS_UTL_fileIO.H"
#include "AS_UTL_splitToWords.H"
#include "AS_MSG_pmesg.H"
#include "AS_PER_gkpStore.H"

gkStore   *gkpStore = 0L;

#include "refAlignment.H"

fragmentData    *frg    = 0L;
uint32           frgLen = 0;

mappingData     *ali    = 0L;
uint32           aliLen = 0;
uint32           aliMax = 0;

uint32           minOverlap = 30;

//  Read all mappings
//  Analyze mappings to find correct trim (smallest, median, largest, best overlap)
//  Discard mappings smaller than trim, adjust larger
//  Output overlaps


static
void
readMapping(char *filename) {
  char  L[1024];
  FILE *F;

  errno = 0;
  F = fopen(filename, "r");
  if (errno)
    fprintf(stderr, "Failed to open input mapping '%s': %s\n", filename, strerror(errno)), exit(1);

  fgets(L, 1024, F);

  while (!feof(F)) {
    if (aliLen >= aliMax) {
      if (aliMax == 0)
        aliMax = 1048576 / 2;
      aliMax *= 2;
      mappingData *A = new mappingData [aliMax];
      memcpy(A, ali, sizeof(mappingData) * aliLen);
      delete [] ali;
      ali = A;
    }

    mappingDataParseLine(ali + aliLen++, L);

    fgets(L, 1024, F);
  }

  fclose(F);

  if (ali == 0L)
    fprintf(stderr, "Failed to read alignments from '%s'.\n", filename);

  qsort(ali, aliLen, sizeof(mappingData), mappingDataSortByRefPosition);

  //  Add sentinel to the end of ali, so that the main processing loop
  //  will find the last contig naturally.

  ali[aliLen].frgIID  = 0;
  ali[aliLen].frgUID  = AS_UID_undefined();
  ali[aliLen].frgBgn  = 0;
  ali[aliLen].frgEnd  = 0;
  ali[aliLen].fwd     = 0;
  ali[aliLen].refUID  = AS_UID_undefined();
  ali[aliLen].refBgn  = 0;
  ali[aliLen].refEnd  = 0;

  aliLen++;

  fprintf(stderr, "Read "F_U32" alignments.\n", aliLen);
}


static
void
loadFragments(void) {
  gkFragment  fr;
  uint32      err=0, alive=0, dead=0, unmapped=0;

  frgLen = getNumGateKeeperFragments(gkpStore) + 1;
  frg    = (fragmentData *)safe_calloc(frgLen, sizeof(fragmentData));

  memset(frg, 0, sizeof(fragmentData) * frgLen);

  //  Load the fragment data.  Set the mappingDataIndex to an invalid
  //  value, we'll set the real value in the loop after this.

  for (uint32 i=1; i<frgLen; i++) {
    getFrag(gkpStore, i, &fr, FRAG_S_INF);

    if (getFragRecordIsDeleted(&fr)) {
      dead++;
      frg[i].fragUID    = AS_UID_undefined();
      frg[i].mateIID    = 0;
      frg[i].libraryIID = 0;
    } else {
      alive++;
      frg[i].fragUID    = getFragRecordUID(&fr);
      frg[i].mateIID    = getFragRecordMateIID(&fr);
      frg[i].libraryIID = getFragRecordLibraryIID(&fr);
    }

    frg[i].mappingDataIndex = aliLen;
  }

  //  For each mapping, get the read IID, and update the
  //  mappingDataIndex pointer from the fragment data to the mapping
  //  data.
  //
  //  aliLen-1 since the last thing is a sentinel end of list marker,
  //  not a real alignment.

  unmapped = alive;

  for (uint32 i=0; i<aliLen-1; i++) {
    uint32  iid = ali[i].frgIID;

    if (AS_UID_compare(frg[iid].fragUID, AS_UID_undefined()) == 0)
      fprintf(stderr, "ERROR:  Fragment %s,%d is deleted in gkpStore, but used in the mapping.\n",
              AS_UID_toString(frg[iid].fragUID), iid), err++;

    if (frg[iid].mappingDataIndex != aliLen)
      fprintf(stderr, "ERROR:  Fragment %s,%d appears more than once in the mapping.\n",
              AS_UID_toString(frg[iid].fragUID), iid), err++;

    unmapped--;

    frg[iid].mappingDataIndex = i;
  }

  if (err)
    fprintf(stderr, "There were errors in the mapping.  Fail.\n"), exit(1);

  fprintf(stderr, "Found "F_U32" alive, "F_U32" dead, and "F_U32" unmapped fragments.\n", alive, dead, unmapped);
}







int
main(int argc, char **argv) {
  char    *mappingFileName   = 0L;
  char    *gkpStoreName      = 0L;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      mappingFileName = argv[++arg];

    } else if (strcmp(argv[arg], "-minoverlap") == 0) {
      minOverlap = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "Invalid option: '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpStoreName == 0L) || (mappingFileName == 0L)) {
    fprintf(stderr, "usage: %s [-U | -S] -g gkpStore -m mapping\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -g gkpStore\n");
    fprintf(stderr, "  -m mapping\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -minoverlap    fragments must overlap by at least this much\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  gkpStore = gkStore_open(gkpStoreName, false);
  if (gkpStore == 0L)
    fprintf(stderr, "Failed to open gkpStore '%s'.\n", gkpStoreName);

  readMapping(mappingFileName);
  loadFragments();

  outputOverlaps();

  exit(0);
}
