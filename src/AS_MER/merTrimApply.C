
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

const char *mainid = "$Id: merTrimApply.C,v 1.5 2011-12-29 09:26:03 brianwalenz Exp $";

#include "AS_global.H"
#include "AS_PER_gkpStore.H"
#include "AS_UTL_splitToWords.H"

#include "merTrimResult.H"


void
merTrimApply(char *gkpStoreName,
             char *listName,
             char *logName) {
  char  resultFileName[FILENAME_MAX] = {0};
  char  resultLine[1024] = {0};

  errno = 0;
  FILE *listFile = fopen(listName, "r");
  if (errno)
    fprintf(stderr, "Failed to open -L '%s' for reading: %s\n", listName, strerror(errno)), exit(1);

  errno = 0;
  FILE *logFile = (logName != NULL) ? fopen(logName, "w") : NULL;
  if (errno)
    fprintf(stderr, "Failed to open -l '%s' for writing: %s\n", logName, strerror(errno)), exit(1);

  gkStore *gkpStore = new gkStore(gkpStoreName, FALSE, TRUE);

  gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_OBTINITIAL);
  gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_TNT);
  gkpStore->gkStore_metadataCaching(true);

  bool            tntEnabled = false;

  //  Over every result file, read each line, change the clear range.

  gkFragment      gkf;
  mertrimResult   res;

  fgets(resultFileName, FILENAME_MAX, listFile);
  while (!feof(listFile)) {
    chomp(resultFileName);

    errno = 0;
    FILE *resultFile = fopen(resultFileName, "r");
    if (errno)
      fprintf(stderr, "Failed to open result file '%s' for reading: %s\n", resultFileName, strerror(errno)), exit(1);

    while (res.readResult(resultFile)) {
      gkpStore->gkStore_getFragment(res.readIID, &gkf, GKFRAGMENT_INF);

      if (res.chimer)
        gkf.gkFragment_setClearRegion(res.chmBgn, res.chmEnd, AS_READ_CLEAR_TNT);

      gkf.gkFragment_setClearRegion(res.clrBgn, res.clrEnd, AS_READ_CLEAR_OBTINITIAL);

      if (res.deleted)
        gkpStore->gkStore_delFragment(res.readIID);
      else
        gkpStore->gkStore_setFragment(&gkf);

      res.print(logFile);
    }

    fclose(resultFile);
    fgets(resultFileName, FILENAME_MAX, listFile);
  }

  fclose(listFile);
  fclose(logFile);

  delete gkpStore;
}




void
merTrimShow(char *gkpStoreName,
            char *resultName) {

  errno = 0;
  FILE *resultFile = fopen(resultName, "r");
  if (errno)
    fprintf(stderr, "Failed to open result file '%s' for reading: %s\n", resultName, strerror(errno)), exit(1);

  mertrimResult   res;

  while (res.readResult(resultFile))
    res.print(stdout);
}





int
main(int argc, char **argv) {
  char     *gkpStoreName = NULL;
  char     *listName     = NULL;
  char     *logName      = NULL;
  char     *resultName   = NULL;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-L") == 0) {
      listName = argv[++arg];

    } else if (strcmp(argv[arg], "-l") == 0) {
      logName = argv[++arg];

    } else if (strcmp(argv[arg], "-d") == 0) {
      resultName = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if (((gkpStoreName == NULL) && (listName != NULL)) ||
      ((gkpStoreName != NULL) && (listName == NULL)))
    err++;
  if ((listName == NULL) && (resultName == NULL))
    err++;
  if ((listName != NULL) && (resultName != NULL))
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -L merTrimOutputList -g gkpStore [-l output.log]\n", argv[0]);
    fprintf(stderr, "       %s -d merTrimOutput\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  The first form will read a list of merTrim output names from\n");
    fprintf(stderr, "  merTrimOuptutList, and apply the results to gkpStore.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  The second form will read a single merTrimOutput file and decode\n");
    fprintf(stderr, "  the results to stdout.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    if (((gkpStoreName == NULL) && (listName != NULL)) ||
        ((gkpStoreName != NULL) && (listName == NULL)))
      fprintf(stderr, "ERROR:  First form needs both -L and -g.\n");
    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR:  No gatekeeper store supplied with -g.\n");
    if ((listName == NULL) && (resultName == NULL))
      fprintf(stderr, "ERROR:  Exactly one of -L and -d must be specified.\n");
    if ((listName != NULL) && (resultName != NULL))
      fprintf(stderr, "ERROR:  Only one of -L and -d can be specified.\n");
    fprintf(stderr, "\n");
    exit(1);
  }


  if (listName)
    merTrimApply(gkpStoreName, listName, logName);

  if (resultName)
    merTrimShow(gkpStoreName, resultName);

  return(0);
}

