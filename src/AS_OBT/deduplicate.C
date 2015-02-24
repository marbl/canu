
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

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "clearRangeFile.H"


#define FRAG_HANG_SLOP  0
#define DEFAULT_ERATE   2.0 / 100.0


int
main(int argc, char **argv) {
  uint32             errorLimit   = AS_OVS_encodeQuality(DEFAULT_ERATE);

  char              *gkpName      = NULL;
  char              *ovsName      = NULL;

  char              *iniClrName   = NULL;
  char              *outClrName   = NULL;

  char              *summaryName  = NULL;
  char              *logName      = NULL;

  bool               doUpdate     = true;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovsName = argv[++arg];


    } else if (strcmp(argv[arg], "-I") == 0) {
      iniClrName = argv[++arg];

    } else if (strcmp(argv[arg], "-D") == 0) {
      outClrName = argv[++arg];


    } else if (strcmp(argv[arg], "-erate") == 0) {
      double erate = atof(argv[++arg]);
      if (erate >= AS_MAX_ERROR_RATE)
        fprintf(stderr, "Error rate %s too large; must be 'fraction error' and below %f\n", argv[arg], AS_MAX_ERROR_RATE), exit(1);
      errorLimit = AS_OVS_encodeQuality(erate);


    } else if (strcmp(argv[arg], "-summary") == 0) {
      summaryName = argv[++arg];

    } else if (strcmp(argv[arg], "-log") == 0) {
      summaryName = argv[++arg];

    } else if (strcmp(argv[arg], "-n") == 0) {
      doUpdate = false;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((gkpName == 0L) || (ovsName == 0L) || (err)) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "  -G gkpStore     \n");
    fprintf(stderr, "  -O obtStore     \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -I iniClr       initial clear range input\n");
    fprintf(stderr, "  -D gkpStore     deduplicated clear ranges output\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -erate E        filter overlaps above this fraction error; default 0.015 (== 1.5%% error)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -summary S      write a summary of the fixes to S\n");
    fprintf(stderr, "  -report R       write a detailed report of the fixes to R\n");
    exit(1);
  }


  gkStore   *gkp = new gkStore(gkpName);
  ovStore   *ovs = new ovStore(ovsName);

  clearRangeFile  *iniClr = (iniClrName == NULL) ? NULL : new clearRangeFile(iniClrName);
  clearRangeFile  *outClr = (outClrName == NULL) ? NULL : new clearRangeFile(outClrName);

  FILE  *summaryFile = NULL;
  FILE  *logFile     = NULL;

  if (summaryName) {
    errno = 0;
    summaryFile = fopen(summaryName, "w");
    if (errno)
      fprintf(stderr, "Failed to open summary file '%s' for writing: %s\n", summaryFile, strerror(errno)), exit(1);
  }

  if (logName) {
    errno = 0;
    logFile = fopen(logName, "w");
    if (errno)
      fprintf(stderr, "Failed to open log file '%s' for writing: %s\n", logFile, strerror(errno)), exit(1);
  }





  bool            nothingToDo = true;

  for (uint32 i=1; i<=gkp->gkStore_getNumLibraries(); i++) {
    gkLibrary  *gkl = gkp->gkStore_getLibrary(i);

    if (gkl->gkLibrary_removeDuplicateReads() == true) {
      if (summaryFile)
        fprintf(summaryFile, "Checking library %s for duplicates.\n", gkl->gkLibrary_libraryName());
      nothingToDo = false;
    } else {
      if (summaryFile)
        fprintf(summaryFile, "Ignoring library %s.\n", gkl->gkLibrary_libraryName());
    }
  }

#if 0
  //  We never set a range...
  uint64  numOlaps += ovs->numOverlapsInRange();

  if (numOlaps == 0)
    nothingToDo = true;
#endif

  //  Stats on overlaps.

  uint64  nFwd  = 0;
  uint64  nRev  = 0;

  uint64  nLQ   = 0;
  uint64  nDelA = 0;
  uint64  nDelB = 0;
  uint64  nLib  = 0;
  uint64  nNoD  = 0;

  //  Stats on duplicates.

  uint32  nDups = 0;


  //
  //  Need to skip overlaps for libraries we aren't dedup'ing
  //



  uint32        ovlLen    = 0;
  uint32        ovlMax    = 4 * 1024 * 1024;
  ovsOverlap   *ovlBuffer = new ovsOverlap [ovlMax];

  ovlLen = ovs->readOverlaps(ovlBuffer, ovlMax, false);

  while (ovlLen > 0) {

    for (uint32 oo=0; oo<ovlLen; oo++) {
      ovsOverlap *ovl = ovlBuffer + oo;

      if (ovl->flipped() == true) {
        //  Dups must be forward...
        nRev++;
        continue;
      }

      nFwd++;

      if (ovl->evalue() > errorLimit) {
        //  ...and of good quality
        nLQ++;
        continue;
      }

      if ((iniClr) && (iniClr->isDeleted(ovl->a_iid))) {
        //  ...and not deleted already
        nDelA++;
        continue;
      }

      if ((iniClr) && (iniClr->isDeleted(ovl->b_iid))) {
        //  ...and not deleted already
        nDelB++;
        continue;
      }

      uint32  aLibID = gkp->gkStore_getRead(ovl->a_iid)->gkRead_libraryID();
      uint32  bLibID = gkp->gkStore_getRead(ovl->b_iid)->gkRead_libraryID();

      if (aLibID != bLibID) {
        //  ...and in the same library
        nLib++;
        continue;
      }

      if (gkp->gkStore_getLibrary(aLibID)->gkLibrary_removeDuplicateReads() == false) {
        //  ...and allowed to be duplicates
        nNoD++;
        continue;
      }

      //  Some of these variables might be leftover from the mate-based duplicate detection.

      int32  ab        = ovl->a_bgn();
      int32  ae        = ovl->a_end(gkp);
      int32  bb        = ovl->b_bgn();
      int32  be        = ovl->b_end(gkp);

      int32  aclrbgn   = (iniClr) ? (iniClr->bgn(ovl->a_iid)) : 0;
      int32  bclrbgn   = (iniClr) ? (iniClr->bgn(ovl->b_iid)) : 0;

      int32  aclrend   = (iniClr) ? (iniClr->end(ovl->a_iid)) : (gkp->gkStore_getRead(ovl->a_iid)->gkRead_sequenceLength());
      int32  bclrend   = (iniClr) ? (iniClr->end(ovl->b_iid)) : (gkp->gkStore_getRead(ovl->b_iid)->gkRead_sequenceLength());

      int32  abeg     = ab + aclrbgn;
      int32  bbeg     = bb + bclrbgn;
      int32  aend     = ae + aclrbgn;
      int32  bend     = be + bclrbgn;

      int32  ahang    = bbeg - abeg;
      int32  bhang    = bend - aend;

      int32  abegdiff = ab;
      int32  bbegdiff = bb;
      int32  aenddiff = aclrend - aclrbgn - ae;
      int32  benddiff = bclrend - bclrbgn - be;

      double error    = ovl->erate();

      //  For unmated reads, delete if it is a near perfect prefix of something else.
      //
      //  Since these are partial overlaps, we need to check both that the overlap covers about the
      //  same piece of each fragment, and that it extends to the start of each fragment.
      //
      //  To pick the longer fragment, we then want to make sure the overlap extends to the end of
      //  this fragment, and that this fragment is contained in the other.

      if ((ahang >= -FRAG_HANG_SLOP) && (ahang <= FRAG_HANG_SLOP) &&
          (abegdiff <= FRAG_HANG_SLOP) &&
          (bbegdiff <= FRAG_HANG_SLOP) &&
          (aenddiff <= FRAG_HANG_SLOP) && (bhang >= 0) &&
          (error    <= 0.025)) {
        fprintf(logFile, "Delete %u DUPof %u  a %d,%d  b %d,%d  hang %d,%d  diff %d,%d  error %f\n",
                ovl->a_iid,
                ovl->b_iid,
                abeg, aend,
                bbeg, bend,
                ahang, bhang,
                abegdiff, bbegdiff,
                error);

        nDups++;

        if ((outClr) && (doUpdate))
          outClr->setDeleted(ovl->a_iid);
      }
    }

    ovlLen = ovs->readOverlaps(ovlBuffer, ovlMax, false);
  }

  delete    iniClr;
  delete    outClr;

  delete [] ovlBuffer;

  delete    gkp;
  delete    ovs;


  if (logFile)
    fclose(logFile);


  if (summaryFile) {
    fprintf(summaryFile, "READ STATS\n");
    fprintf(summaryFile, "duplicateReads:    "F_U32"\n", nDups);
    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "OVERLAP STATS\n");
    fprintf(summaryFile, "forward orient      "F_U64" (can indicate duplicates)\n", nFwd);
    fprintf(summaryFile, "reverse orient      "F_U64" (can't indicate duplicates)\n", nRev);
    fprintf(summaryFile, "low quality         "F_U64"\n", nLQ);
    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "A read deleted (target)         "F_U64"\n", nDelA);
    fprintf(summaryFile, "B read deleted (evidence)       "F_U64"\n", nDelB);
    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "different libraries "F_U64"\n", nLib);
    fprintf(summaryFile, "dedup not allowed   "F_U64"\n", nNoD);
    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "\n");

    fclose(summaryFile);
  }

  exit(0);
}
