
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2012, J. Craig Venter Institute.
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

#include "AS_global.H"

#include "gkStore.H"

#include "splitToWords.H"
#include "findKeyAndValue.H"

#include "AS_UTL_fileIO.H"



void
getLine(char *inLine, uint32 inLineLen, FILE *inFile) {

  do {
    fgets(inLine, inLineLen, inFile);
    chomp(inLine);
  } while ((inLine[0] == '#') || (inLine[0] == 0));
}



int
main(int argc, char **argv) {
  char            *gkpStoreName      = NULL;
  char            *outPrefix         = NULL;

  bool             ignoreClear       = true;

  uint32           firstFileArg      = 0;

  char             errorLogName[FILENAME_MAX];
  char             uidMapName[FILENAME_MAX];


  //argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-useclear") == 0) {
      ignoreClear = false;


    } else if (strcmp(argv[arg], "--") == 0) {
      firstFileArg = arg++;
      break;

    } else if (argv[arg][0] == '-') {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;

    } else {
      firstFileArg = arg;
      break;
    }
    arg++;
  }

  if (gkpStoreName == NULL)
    err++;
  if (firstFileArg == 0)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [...] -o gkpStore\n", argv[0]);
    fprintf(stderr, "  -o gkpStore         create this gkpStore\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  -useclear           enforce clear ranges on the name lines\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  \n");

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-g) supplied.\n");
    if (firstFileArg == 0)
      fprintf(stderr, "ERROR: no input files supplied.\n");

    exit(1);
  }


  gkStore     *gkpStore = new gkStore(gkpStoreName, gkStore_extend);
  gkRead      *gkpRead;
  gkReadData   gkpReadData;
  gkLibrary   *gkpLibrary;

  uint32       inLineLen = 1024;
  char         inLine[1024];

  errno = 0;

  sprintf(errorLogName, "%s.errorLog",    gkpStoreName);
  FILE    *errorLog = fopen(errorLogName, "w");
  if (errno)
    fprintf(stderr, "ERROR: cannot open error file '%s': %s\n", errorLogName, strerror(errno)), exit(1);

  sprintf(uidMapName,   "%s.fastqUIDmap", gkpStoreName);
  FILE    *uidMap   = fopen(uidMapName,   "w");
  if (errno)
    fprintf(stderr, "ERROR: cannot open uid map file '%s': %s\n", uidMapName, strerror(errno)), exit(1);


  for (; firstFileArg < argc; firstFileArg++) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Starting file '%s'.\n", argv[firstFileArg]);

    compressedFileReader *inFile = new compressedFileReader(argv[firstFileArg]);
    char                 *line   = new char [10240];
    KeyAndValue           keyval;

    fgets(line, 10240, inFile->file());
    chomp(line);

    while (!feof(inFile->file())) {
      fprintf(stderr, "LINE '%s'\n", line);

      keyval.find(line);

      if        (keyval.key() == NULL) {
        //  Just a damn comment line to ignore.

      } else if (strcasecmp(keyval.key(), "name") == 0) {
        gkpLibrary = gkpStore->gkStore_addEmptyLibrary(keyval.value());

      } else if (strcasecmp(keyval.key(), "preset") == 0) {

      } else if (strcasecmp(keyval.key(), "qv") == 0) {

      } else if (strcasecmp(keyval.key(), "isNotRandom") == 0) {

      } else if (strcasecmp(keyval.key(), "doNotTrustHomopolymerRuns") == 0) {

      } else if (strcasecmp(keyval.key(), "initialTrim") == 0) {

      } else if (strcasecmp(keyval.key(), "removeDuplicateReads") == 0) {

      } else if (strcasecmp(keyval.key(), "finalTrim") == 0) {

      } else if (strcasecmp(keyval.key(), "removeSpurReads") == 0) {

      } else if (strcasecmp(keyval.key(), "removeChimericReads") == 0) {

      } else if (strcasecmp(keyval.key(), "removeSubReads") == 0) {

      } else if (AS_UTL_fileExists(keyval.key(), false, false)) {

      } else {
        fprintf(stderr, "line '%s' invalid.\n", line);
      }

      fgets(line, 10240, inFile->file());
      chomp(line);
    }
  }

 noMoreLines:

  delete gkpStore;

  fclose(uidMap);
  fclose(errorLog);

  //return(AS_GKP_summarizeErrors());
  exit(0);
}
