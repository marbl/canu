
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2010, J. Craig Venter Institute
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

const char *mainid = "$Id: merTrimApply.C,v 1.2 2011-01-06 19:52:08 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

#include "AS_UTL_splitToWords.H"

int
main(int argc, char **argv) {
  char     *gkpStoreName = NULL;
  gkStore  *gkpStore     = NULL;

  char     *listName     = NULL;
  FILE     *listFile     = NULL;

  char     *logName      = NULL;
  FILE     *logFile      = NULL;

  char      resultFileName[FILENAME_MAX];
  FILE     *resultFile   = NULL;
  char      resultLine[1024];

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-L") == 0) {
      listName = argv[++arg];

    } else if (strcmp(argv[arg], "-l") == 0) {
      logName = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if (err) {
    fprintf(stderr, "usage: %s -g gkpStore -L merTrimOutputList -l output.log\n", argv[0]);
    exit(1);
  }

  errno = 0;
  listFile = fopen(listName, "r");
  if (errno)
    fprintf(stderr, "Failed to open -L '%s' for reading: %s\n", listName, strerror(errno)), exit(1);

  errno = 0;
  logFile = fopen(logName, "w");
  if (errno)
    fprintf(stderr, "Failed to open -l '%s' for writing: %s\n", logName, strerror(errno)), exit(1);

  gkpStore = new gkStore(gkpStoreName, FALSE, TRUE);

  gkpStore->gkStore_metadataCaching(true);
  gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_OBTINITIAL);

  //  Over every result file, read each line, change the clear range.

  splitToWords  S;

  gkFragment    gkf;

  fgets(resultFileName, FILENAME_MAX, listFile);
  while (!feof(listFile)) {
    chomp(resultFileName);

    errno = 0;
    FILE *resultFile = fopen(resultFileName, "r");
    if (errno)
      fprintf(stderr, "Failed to open result file '%s' for reading: %s\n", resultFileName, strerror(errno)), exit(1);

    fgets(resultLine, 1024, resultFile);
    while (!feof(resultFile)) {
      chomp(resultLine);

      //  process the result

      //  Change the comma between the UID and the IID to a space.
      for (int i=0; resultLine[i]; i++)
        if (resultLine[i] == ',') {
          resultLine[i] = ' ';
          break;
        }

      S.split(resultLine);

      uint32  iid = S(1);
      uint32  bgn = S(2);
      uint32  end = S(3);

      gkpStore->gkStore_getFragment(iid, &gkf, GKFRAGMENT_INF);

      //  Update the clear range, or just delete the fragment.

      if (S.numWords() == 4) {
        gkf.gkFragment_setClearRegion(bgn, end, AS_READ_CLEAR_OBTINITIAL);
        gkpStore->gkStore_setFragment(&gkf);

        fprintf(logFile, "%s,"F_U32"\t"F_U32"\t"F_U32"\n", S[0], iid, bgn, end);

      } else {
        gkpStore->gkStore_delFragment(iid);

        fprintf(logFile, "%s,"F_U32"\t"F_U32"\t"F_U32"\t(deleted)\n", S[0], iid, bgn, end);
      }

      fgets(resultLine, 1024, resultFile);
    }

    fclose(resultFile);
    fgets(resultFileName, FILENAME_MAX, listFile);
  }

  fclose(listFile);
  fclose(logFile);

  delete gkpStore;

  return(0);
}

