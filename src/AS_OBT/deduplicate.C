
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
#include "adjustOverlaps.H"


int
main(int argc, char **argv) {
  uint32  evalueLimit  = AS_OVS_encodeQuality(0.02);
  uint32  diff5        = 0;
  uint32  diff3        = 0;

  char   *gkpName      = NULL;
  char   *ovsName      = NULL;

  char   *iniClrName   = NULL;
  char   *outClrName   = NULL;

  char   *outputPrefix = NULL;
  char    logName[FILENAME_MAX] = {0};
  char    sumName[FILENAME_MAX] = {0};
  FILE   *logFile = 0L;
  FILE   *sumFile = 0L;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovsName = argv[++arg];


    } else if (strcmp(argv[arg], "-Ci") == 0) {
      iniClrName = argv[++arg];

    } else if (strcmp(argv[arg], "-Co") == 0) {
      outClrName = argv[++arg];


    } else if (strcmp(argv[arg], "-e") == 0) {
      double erate = atof(argv[++arg]);
      evalueLimit = AS_OVS_encodeQuality(erate);

    } else if (strcmp(argv[arg], "-d5") == 0) {
      diff5 = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-d3") == 0) {
      diff3 = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

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
    fprintf(stderr, "  -Ci iniClr      initial clear range input (usually not supplied)\n");
    fprintf(stderr, "  -Co outClr      deduplicated clear ranges output\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e E            filter overlaps above this fraction error\n");
    fprintf(stderr, "                    The default is 0.02 (or.0%% error)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -d5 diff        overlap position must differ by more than 'diff' bp at the 5' end of the reads\n");
    fprintf(stderr, "  -d3 diff        overlap position must differ by more than 'diff' bp at the 3' end of the reads\n");
    fprintf(stderr, "                    e.g., an overlap from bases 3-10 on the A read and from 4-8 on the B read\n");
    fprintf(stderr, "                          would have a difference of 1 at the 5' end and 2 on the 3' end\n");
    fprintf(stderr, "                    The default for both is 0\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o prefix       write logging to prefix.log and prefix.summary\n");
    exit(1);
  }


  gkStore         *gkp = new gkStore(gkpName);
  ovStore         *ovs = new ovStore(ovsName);

  //  Output logging outputs.

  if (outputPrefix) {
    sprintf(logName, "%s.log",     outputPrefix);
    sprintf(sumName, "%s.summary", outputPrefix);

    errno = 0;
    logFile = fopen(logName, "w");
    if (errno)
      fprintf(stderr, "Failed to open log file '%s' for writing: %s\n", logName, strerror(errno)), exit(1);

    fprintf(logFile, "evidence <= deleted\n");

    sumFile = fopen(sumName, "w");
    if (errno)
      fprintf(stderr, "Failed to open summary file '%s' for writing: %s\n", sumName, strerror(errno)), exit(1);
  }

  //  Open inputs and outputs.  The outClr must always exist, but the code looks so much better this way.

  clearRangeFile  *iniClr = (iniClrName == NULL) ? NULL : new clearRangeFile(iniClrName, gkp);
  clearRangeFile  *outClr = (outClrName == NULL) ? NULL : new clearRangeFile(outClrName, gkp);

  //  Copy clear ranges to the output.  The dedup process will only delete reads.

  outClr->copy(iniClr);


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
  uint64  nSkipped = 0;
  uint64  nNotSame = 0;
  uint64  nNotEtoE = 0;

  //  Stats on duplicates.

  uint32  nDups = 0;


  //
  //  Need to skip overlaps for libraries we aren't dedup'ing
  //

#if 0
  bool            nothingToDo = true;

  for (uint32 i=1; i<=gkp->gkStore_getNumLibraries(); i++) {
    gkLibrary  *gkl = gkp->gkStore_getLibrary(i);

    if (gkl->gkLibrary_removeDuplicateReads() == true) {
      if (sumFile)
        fprintf(sumFile, "Checking library '%s' for duplicates.\n", gkl->gkLibrary_libraryName());
      nothingToDo = false;

    } else {
      if (sumFile)
        fprintf(sumFile, "Ignoring library '%s'.\n", gkl->gkLibrary_libraryName());
    }
  }

  if (nothingToDo) {
    if (sumFile)
      fprintf(sumFile, "\nAll libraies skipped.  Bye.\n");

    if (logFile)
      fclose(logFile);
    if (sumFile)
      fclose(sumFile);

    exit(0);
  }
#endif


  uint32        ovlLen    = 0;
  uint32        ovlMax    = 4 * 1024 * 1024;
  ovsOverlap   *ovlBuffer = new ovsOverlap [ovlMax];

  ovlLen = ovs->readOverlaps(ovlBuffer, ovlMax, false);

  //  We don't need all overlaps for a single read.  All overlaps are independent, and can be processed in any ol' order.
  //  That's the 'false' flag above, just load 4m overlaps.

  while (ovlLen > 0) {
    for (uint32 oo=0; oo<ovlLen; oo++) {
      ovsOverlap *ovl = ovlBuffer + oo;

      //  Dups must be forward...

      if (ovl->flipped() == true) {
        nRev++;
        continue;
      }

      nFwd++;

      //  ...and of good quality

      if (ovl->evalue() > evalueLimit) {
        nLQ++;
        continue;
      }

      //  ...and not deleted already

      if ((iniClr) && (iniClr->isDeleted(ovl->a_iid))) {
        nDelA++;
        continue;
      }

      //  ...and not deleted already (also)

      if ((iniClr) && (iniClr->isDeleted(ovl->b_iid))) {
        nDelB++;
        continue;
      }

      uint32  aLibID = gkp->gkStore_getRead(ovl->a_iid)->gkRead_libraryID();
      uint32  bLibID = gkp->gkStore_getRead(ovl->b_iid)->gkRead_libraryID();

      //  ...and in the same library

      if (aLibID != bLibID) {
        nLib++;
        continue;
      }

      //  ...and allowed to be duplicates

      if (gkp->gkStore_getLibrary(aLibID)->gkLibrary_removeDuplicateReads() == false) {
        nNoD++;
        continue;
      }

      //  ...and not marked to ignore.  In normal usage, this will catch nothing, as the testa above caught the problem ones.

      if (ovl->forDUP() == false) {
        nSkipped++;
        continue;
      }

      //  All overlaps are forward here.  Adjust the overlap by any trimming already done.  If there is no iniClr range
      //  supplied, there was no trimming, and the overlap isn't adjusted.

      uint32  aovlbgn = ovl->a_bgn(gkp);
      uint32  aovlend = ovl->a_end(gkp);
      uint32  bovlbgn = ovl->b_bgn(gkp);
      uint32  bovlend = ovl->b_end(gkp);

      uint32  aclrbgn = 0;
      uint32  aclrend = gkp->gkStore_getRead(ovl->a_iid)->gkRead_sequenceLength();
      uint32  bclrbgn = 0;
      uint32  bclrend = gkp->gkStore_getRead(ovl->b_iid)->gkRead_sequenceLength();

      if ((iniClr) && (adjustNormal(iniClr, ovl,
                                    aovlbgn, aovlend, bovlbgn, bovlend,
                                    aclrbgn, aclrend, bclrbgn, bclrend) == false)) {
        //nTrimmed++;
        continue;
      }

      //  A read is a duplicate if it is contained in the other one, and the other one is contained in it.

      if ((aclrbgn + diff5 <= aovlbgn) ||
          (bclrbgn + diff5 <= bovlbgn) ||
          (aovlend + diff3 <= aclrend) ||
          (bovlend + diff3 <= bclrend)) {
        nNotSame++;
        continue;
      }


#if 0
      int32  bgnDiff   = (abgn < bbgn) ? (bbgn - abgn) : (abgn - bbgn);
      int32  endDiff   = (aend < bend) ? (bend - aend) : (aend - bend);

      assert(bgnDiff >= 0);
      assert(endDiff >= 0);

      //  ...and the overlap starts at the same point in each read
      //  (which must happen for the next rule to pass)

      if ((bgnDiff > diff5) ||
          (endDiff > diff3)) {
        nNotSame++;
        continue;
      }

      uint32  alen = aend - abgn;  //gkp->gkStore_getRead(ovl->a_iid)->gkRead_sequenceLength();
      uint32  blen = bend - bbgn;  //gkp->gkStore_getRead(ovl->b_iid)->gkRead_sequenceLength();

      //  ...and the overlap is end-to-end

      if ((abgn >= bgnDiff) || (aend + bgnDiff >= al) ||
          (bbgn >= bgnDiff) || (bend + bgnDiff >= bl)) {
        nNotEtoE++;
        continue;
      }
#endif


      //  Otherwise, looks like a duck, quacks like a duck...a witch!
      //  But don't overcount duplicates!

      fprintf(logFile, "%u <= %u\n", ovl->b_iid, ovl->a_iid);

      if (outClr->isDeleted(ovl->a_iid) == false) {
        nDups++;
        outClr->setDeleted(ovl->a_iid);
      }
    }

    ovlLen = ovs->readOverlaps(ovlBuffer, ovlMax, false);
  }

  if (logFile)
    fclose(logFile);

  delete [] ovlBuffer;

  delete    gkp;
  delete    ovs;


  uint32  nDeld = 0;
  uint32  nActv = 0;

  for (uint32 ii=1; ii <= gkp->gkStore_getNumReads(); ii++) {
    if (outClr->isDeleted(ii))
      nDeld++;
    else
      nActv++;
  }

  delete    iniClr;
  delete    outClr;

  if (sumFile) {
    fprintf(sumFile, "READ STATS\n");
    fprintf(sumFile, "\n");
    fprintf(sumFile, "duplicateReads:     %16"F_U32P" (reads deleted here)\n", nDups);
    fprintf(sumFile, "\n");
    fprintf(sumFile, "deletedReads        %16"F_U32P" (final state)\n", nDeld);
    fprintf(sumFile, "activeReads         %16"F_U32P" (final state)\n", nActv);
    fprintf(sumFile, "\n");
    fprintf(sumFile, "OVERLAP STATS\n");
    fprintf(sumFile, "\n");
    fprintf(sumFile, "forward orient      %16"F_U64P" (can indicate duplicates)\n", nFwd);
    fprintf(sumFile, "reverse orient      %16"F_U64P" (can't indicate duplicates)\n", nRev);
    fprintf(sumFile, "low quality         %16"F_U64P" (error rate too high)\n", nLQ);
    fprintf(sumFile, "\n");
    fprintf(sumFile, "A read deleted      %16"F_U64P" (this is the read being tested)\n", nDelA);
    fprintf(sumFile, "B read deleted      %16"F_U64P" (this is the evidence read)\n", nDelB);
    fprintf(sumFile, "\n");
    fprintf(sumFile, "different overlap   %16"F_U64P" (overlap is to different regions on the two reads)\n", nNotSame);
    fprintf(sumFile, "not end-to-end      %16"F_U64P" (overlap doesn't cover both reads fully)\n", nNotEtoE);
    fprintf(sumFile, "\n");
    fprintf(sumFile, "different libraries %16"F_U64P" (overlap is between reads in different libraries)\n", nLib);
    fprintf(sumFile, "dedup not allowed   %16"F_U64P" (we did a bad job filtering overlaps, and tried to dedupe something we shouldn't have)\n", nNoD);
    fprintf(sumFile, "\n");
    fprintf(sumFile, "skipped             %16"F_U64P" (marked not useful for dedupe)\n", nSkipped);
    fclose(sumFile);
  }

  exit(0);
}
