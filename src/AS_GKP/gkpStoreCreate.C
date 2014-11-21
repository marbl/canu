
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
#include "AS_PER_gkpStore.H"

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

  AS_IID           libToDump         = 0;
  uint32           clrToDump         = AS_READ_CLEAR_LATEST;

  AS_IID           bgnIID            = 1;
  AS_IID           endIID            = AS_IID_MAX;

  bool             dumpAllBases      = true;
  bool             dumpAllReads      = false;
  bool             ignoreClear       = true;

  uint32           packedLength      = 150;

  uint32           firstFileArg      = 0;

  argc = AS_configure(argc, argv);

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

    if (gkpStoreName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-g) supplied.\n");
    if (firstFileArg == 0)
      fprintf(stderr, "ERROR: no input files supplied.\n");

    exit(1);
  }


  gkStore     *gkpStore = new gkStore(gkpStoreName, packedLength);
  gkFragment   gkpFrag;
  gkFragment   gkpMate;
  gkLibrary    gkpLibrary;

  uint32       inLineLen = 1024;
  char         inLine[1024];

  //  This is a special case for gatekeeper; we never call gkStore_getFragment() and so we never set
  //  up the gkFragment.  This also enables some of the set() methods that are allowed only when
  //  creating a gkpStore.

  gkpFrag.gkFragment_enableGatekeeperMode(gkpStore);
  gkpMate.gkFragment_enableGatekeeperMode(gkpStore);

  char     errorLogName[FILENAME_MAX];
  char     uidMapName[FILENAME_MAX];

  sprintf(errorLogName, "%s.errorLog",    gkpStoreName);
  sprintf(uidMapName,   "%s.fastqUIDmap", gkpStoreName);

  errno = 0;

  FILE    *errorLog = fopen(errorLogName, "w");
  if (errno)
    fprintf(stderr, "ERROR: cannot open error file '%s': %s\n", errorLogName, strerror(errno)), exit(1);

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
      keyval.find(line);

      if        (strcasecmp(keyval.key(), "libraryname") == 0) {
        //gkpStore->gkStore_addLibrary(lib_mesg->eaccession, &gkpl);

#if 0
    gkpl.libraryUID                 = lib_mesg->eaccession;

    strcpy(gkpl.libraryName, AS_UID_toString(gkpl.libraryUID));

    gkpl.mean                       = lib_mesg->mean;
    gkpl.stddev                     = lib_mesg->stddev;

    gkpl.orientation                = AS_READ_ORIENT_UNKNOWN;

    if (lib_mesg->link_orient.isInnie())
      gkpl.orientation = AS_READ_ORIENT_INNIE;

    if (lib_mesg->link_orient.isOuttie())
      gkpl.orientation = AS_READ_ORIENT_OUTTIE;

    if (lib_mesg->link_orient.isNormal())
      gkpl.orientation = AS_READ_ORIENT_NORMAL;

    if (lib_mesg->link_orient.isAnti())
      gkpl.orientation = AS_READ_ORIENT_ANTINORMAL;

    gkpl.gkLibrary_decodeFeatures(lib_mesg);

    gkpStore->gkStore_addLibrary(lib_mesg->eaccession, &gkpl);

    //  If this library specifies fastq reads, load them now.

    iid = gkpStore->gkStore_getUIDtoIID(lib_mesg->eaccession, NULL);
    assert(iid > 0);

    checkLibraryForFastQPointers(lib_mesg,
                                 iid,
                                 packedLength);
#endif



      } else if (strcasecmp(keyval.key(), "librarypreset") == 0) {

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


    //  Get the library name


#if 0
    getLine(inLine, inLineLen, inFile);
    {
      char    *nmPtr = gkpLibrary.libraryName;
      uint32   inPos = 0;
      uint32   nmLen = 0;

      while (inLine[inPos]) {
        if        (inLine[inPos] == '/') {
          nmPtr[nmLen++] = '_';
          nmPtr[nmLen++] = '-';
          nmPtr[nmLen++] = '_';
        } else if (isspace(inLine[inPos]) == 0) {
          nmPtr[nmLen++] = inLine[inPos];
        } else {
          nmPtr[nmLen++] = '_';
        }

        if (nmLen >= LIBRARY_NAME_SIZE) {
          nmPtr[LIBRARY_NAME_SIZE-1] = 0;
          break;
        }

        inPos++;
      }
    }

    //  Insert size and orientation

    bool reverseThese = false;

    gkpLibrary.orientation = AS_READ_ORIENT_INNIE;

    getLine(inLine, inLineLen, inFile);
    {
      if (strcasecmp(inLine, "fragment") == 0) {
      } else {
        splitToWords  split(inLine);

        if (split.numWords() != 4) {
          fprintf(stderr, "ERROR: bad insert line for library '%s': '%s'\n",
                  gkpLibrary.libraryName, inLine);
          exit(1);
        }

        gkpLibrary.mean   = atof(split[0]);
        gkpLibrary.stddev = atof(split[2]);

        switch (split[3][0]) {
          case 'I':
          case 'i':
            reverseThese = false;
            break;
          case 'O':
          case 'o':
            reverseThese = true;
            break;
          default:
            fprintf(stderr, "ERROR: unsupported library orientation '%s'\n", split[3]);
            exit(1);
            break;
        }
      }
    }

    //  Library features and inputs

  anotherLine:
    getLine(inLine, inLineLen, inFile);

    if (feof(inFile))
      goto noMoreLines;

    char *hzd = inLine;  //  Iterating through the line
    char *fea = inLine;  //  Pointer to feature string
    char *val = inLine;  //  Pointer to value string

    //  Bump fea up to the first non space
    while (isspace(*hzd) == true)
      hzd++;

    fea = hzd;

    if (fea[0] == '#')
      //  Just a comment line.
      goto anotherLine;

    //  Skip up until the first space or '='
    while ((*hzd != '=') &&
           (isspace(*hzd) == false))
      *hzd++;

    //  Skip spaces and '='
    while ((*hzd == '=') ||
           (isspace(*hzd) == true))
      *hzd++ = 0;

    val = hzd;

    //  Terminate the value.

    while ((*hzd) &&
           (isspace(*hzd) == false))
      hzd++;

    *hzd = 0;

    fprintf(stderr, "FEATURE '%s'  VALUE '%s'\n",
            fea, val);

    if ((strcasecmp(fea, "fastqQualityValues") == 0) ||
        (strcasecmp(fea, "fastqOrientation") == 0) ||
        (strcasecmp(fea, "fastqMates") == 0) ||
        (strcasecmp(fea, "fastqReads") == 0))
      fprintf(stderr, "FASTQ - '%s'\n", fea);

    else
  
    goto anotherLine;
#endif


  }

 noMoreLines:

  delete gkpStore;

  fclose(uidMap);
  fclose(errorLog);

  //return(AS_GKP_summarizeErrors());
  exit(0);
}
