
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
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

//  Compute a quality trimming based on quality scores.
//  Intersect the quality trim with any vector trim.
//  Output a binary clear range file.

#include "trim.H"
#include "clearRangeFile.H"

int
main(int argc, char **argv) {
  char   *gkpName             = NULL;
  char   *vecClrName          = NULL;
  char   *iniClrName          = NULL;

  char   *logFileName         = NULL;
  FILE   *logFile             = NULL;

  bool    doUpdate            = true;

  double  maxErate            = 0.06;
  uint32  minReadLength       = 64;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-log", 2) == 0) {
      logFileName = argv[++arg];

    } else if (strncmp(argv[arg], "-G", 2) == 0) {
      gkpName = argv[++arg];

    } else if (strncmp(argv[arg], "-V", 2) == 0) {
      vecClrName = argv[++arg];

    } else if (strncmp(argv[arg], "-I", 2) == 0) {
      iniClrName = argv[++arg];

    } else if (strncmp(argv[arg], "-n", 2) == 0) {
      doUpdate = false;

    } else if (strncmp(argv[arg], "-e", 2) == 0) {
      maxErate = atof(argv[++arg]);
    } else if (strncmp(argv[arg], "-l", 2) == 0) {
      minReadLength = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "Invalid option: '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err) || (!gkpName)) {
    fprintf(stderr, "usage: %s [-q quality] [-update] [-replace] [-log logfile] -frg some.gkpStore\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G   gkpStore   Operate on this gkpStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -V   vecClr     Load vector clear ranges from this (binary) file\n");
    fprintf(stderr, "  -I   vecClr     Save initial trimming clear ranges to this (binary) file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -q quality    Find quality trim points using 'quality' as the base.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -log X        Report the id, original trim and new quality trim\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  A report of the trimming is printed to stdout:\n");
    fprintf(stderr, "    id originalBegin originalEnd newBegin newEnd\n");
    fprintf(stderr, "    id origBegin origEnd qualBegin qualEnd vecBeg vecEnd newBegin newEnd\n");
    exit(1);
  }

  gkStore  *gkpStore = new gkStore(gkpName);

  clearRangeFile  *vecClr = (vecClrName == NULL) ? NULL : new clearRangeFile(vecClrName, gkpStore->gkStore_getNumReads());
  clearRangeFile  *iniClr = (vecClrName == NULL) ? NULL : new clearRangeFile(iniClrName, gkpStore->gkStore_getNumReads());


  if ((logFileName != NULL) && (strcmp(logFileName, "-") == 0)) {
    logFile = stdout;
  }

  if ((logFileName != NULL) && (strcmp(logFileName, "-") == 0)) {
    errno=0;
    logFile = fopen(logFileName, "w");
    if (errno)
      fprintf(stderr, "Failed to open %s for writing the log: %s\n",
              logFileName, strerror(errno)), exit(1);
  }

  //  One per read
  uint32       stat_alreadyDeleted = 0;
  uint32       stat_noTrimming     = 0;
  uint32       stat_merBased       = 0;
  uint32       stat_flowBased      = 0;
  uint32       stat_qualityBased   = 0;

  //  One per read
  uint32       stat_noVecClr       = 0;
  uint32       stat_noHQnonVec     = 0;
  uint32       stat_HQtrim5        = 0;
  uint32       stat_HQtrim3        = 0;
  uint32       stat_LQtrim5        = 0;
  uint32       stat_LQtrim3        = 0;
  uint32       stat_tooShort       = 0;

  //  Used to be a library parameter.
  double       minQuality = qual.lookupNumber(12);

  gkReadData   readData;  //  Can be reused.

  if (logFile)
    fprintf(logFile, "id\torigL\torigR\tqltL\tqltR\tfinalL\tfinalR\tvecL\tvecR\tdeleted?\n");

  for (uint32 id=1; id<=gkpStore->gkStore_getNumReads(); id++) {
    gkRead     *read = gkpStore->gkStore_getRead(id);
    gkLibrary  *libr = gkpStore->gkStore_getLibrary(read->gkRead_libraryID());

    assert(libr != NULL);

    //  A quality based trimming.  Computed here.
    uint32  qltL     = 0;
    uint32  qltR     = read->gkRead_sequenceLength();

    //  A vector based trimming.  Loaded.
    uint32  vecL     = 0;
    uint32  vecR     = UINT32_MAX;

    //  The final initial trimming.  Intersection of quality and vector based trimming.
    //  Initially, unset.
    uint32  finL     = iniClr->bgn(id);
    uint32  finR     = iniClr->end(id);
    bool    updateIt = false;
    bool    deleteIt = false;

    //  Decide on trimming.

    if ((vecClr) && (vecClr->isDeleted(id) == true)) {
      stat_alreadyDeleted++;   //  Vector clear exists, and it already deleted the read.
      deleteIt = true;
      goto update;
    }

    if (libr->gkLibrary_initialTrim() == INITIALTRIM_NONE) {
      stat_noTrimming++;
      goto update;
    }

    if (libr->gkLibrary_initialTrim() == INITIALTRIM_MER_BASED) {
      stat_merBased++;
      goto update;
    }

    //if (libr->gkLibrary_initialTrim() == INITIALTRIM_FLOW_BASED) {
    //  stat_flowBased++;
    //  goto update;
    //}

    if (libr->gkLibrary_initialTrim() == INITIALTRIM_QUALITY_BASED) {
      stat_qualityBased++;

      gkpStore->gkStore_loadReadData(read, &readData);

      doTrim(read, &readData, minQuality, qltL, qltR);
      updateIt = true;
    }

    //  If no vector clear, we're done.

    if ((vecClr == NULL) || (vecClr->isUndefined(id) == true)) {
      stat_noVecClr++;

      finL = qltL;
      finR = qltR;
      goto update;
    }

    //  Otherwise, intersect with the vector clear.

    vecL = vecClr->bgn(id);
    vecR = vecClr->end(id);

    if ((vecL > finR) || (vecR < finL)) {
      //  No intersection; trust nobody.  Delete the read.
      stat_noHQnonVec++;

      deleteIt = true;
      goto update;
    }

    //  They intersect.  Pick the largest begin and the smallest end

    if (finL < vecL) {
      stat_HQtrim5++;
      finL = vecL;
    } else {
      stat_LQtrim5++;
    }

    if (vecR < finR) {
      stat_HQtrim3++;
      finR = vecR;
    } else {
      stat_LQtrim3++;
    }
 


    //  Update the clear ranges
  update:

    if ((finL + minReadLength) > finR)
      stat_tooShort++;

    if (doUpdate) {
      iniClr->setbgn(id) = finL;
      iniClr->setend(id) = finR;
    }

    if (logFile)
      fprintf(logFile, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\n",
              id, 0, read->gkRead_sequenceLength(), qltL, qltR, vecL, vecR, finL, finR);
  }

  delete gkpStore;

  fprintf(stdout, "Reads trimmed using:\n");
  fprintf(stdout, "  alreadyDeleted   "F_U32"\n", stat_alreadyDeleted);
  fprintf(stdout, "  noTrimming       "F_U32"\n", stat_noTrimming);
  fprintf(stdout, "  merBased         "F_U32"\n", stat_merBased);
  fprintf(stdout, "  flowBased        "F_U32"\n", stat_flowBased);
  fprintf(stdout, "  qualityBased     "F_U32"\n", stat_qualityBased);
  
  fprintf(stdout, "Trimming result:\n");
  fprintf(stdout, " no vector clear range known:  "F_U32" (trimed to quality clear)\n", stat_noVecClr);
  fprintf(stdout, " no HQ non-vector sequence:    "F_U32" (deleted)\n", stat_noHQnonVec);
  fprintf(stdout, " HQ vector trimmed:            "F_U32" (trimmed to intersection)\n", stat_HQtrim5 + stat_HQtrim3);
  fprintf(stdout, "     5' end:                   "F_U32"\n", stat_HQtrim5);
  fprintf(stdout, "     3' end:                   "F_U32"\n", stat_HQtrim3);
  fprintf(stdout, " LQ vector trimmed:            "F_U32" (trimmed to intersection)\n", stat_LQtrim5 + stat_LQtrim3);
  fprintf(stdout, "     5' end:                   "F_U32"\n", stat_LQtrim5);
  fprintf(stdout, "     3' end:                   "F_U32"\n", stat_LQtrim3);
  fprintf(stdout, " final clear range too short:  "F_U32" (deleted)\n", stat_tooShort);

  return(0);
}
