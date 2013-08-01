
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

static char const *rcsid = "$Id: AS_GKP_dump.c,v 1.75 2012-09-19 20:06:51 skoren Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_GKP_include.H"
#include "AS_UTL_fasta.H"
#include "AS_UTL_reverseComplement.H"
#include "AS_PER_gkpStore.H"


void
dumpGateKeeperInfo(char       *gkpStoreName,
                   int         asTable) {

  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE, TRUE);

  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  gkStoreStats  *stats = new gkStoreStats(gkp);

  if (asTable == 0) {
    fprintf(stdout, "LOAD STATS\n");
    fprintf(stdout, "\n");
    fprintf(stdout, F_U32"\tlibInput\n",    gkp->inf.libInput);
    fprintf(stdout, F_U32"\tlibLoaded\n",   gkp->inf.libLoaded);
    fprintf(stdout, F_U32"\tlibErrors\n",   gkp->inf.libErrors);
    fprintf(stdout, F_U32"\tlibWarnings\n", gkp->inf.libWarnings);
    fprintf(stdout, "\n");
    fprintf(stdout, F_U32"\tfrgInput\n",    gkp->inf.frgInput);
    fprintf(stdout, F_U32"\tfrgLoaded\n",   gkp->inf.frgLoaded);
    fprintf(stdout, F_U32"\tfrgErrors\n",   gkp->inf.frgErrors);
    fprintf(stdout, F_U32"\tfrgWarnings\n", gkp->inf.frgWarnings);
    fprintf(stdout, "\n");
    fprintf(stdout, F_U32"\tlkgInput\n",    gkp->inf.lkgInput);
    fprintf(stdout, F_U32"\tlkgLoaded\n",   gkp->inf.lkgLoaded);
    fprintf(stdout, F_U32"\tlkgErrors\n",   gkp->inf.lkgErrors);
    fprintf(stdout, F_U32"\tlkgWarnings\n", gkp->inf.lkgWarnings);
    fprintf(stdout, "\n");
    fprintf(stdout, F_U32"\tsffInput\n",    gkp->inf.sffInput);
    fprintf(stdout, F_U32"\tsffLoaded\n",   gkp->inf.sffLoaded);
    fprintf(stdout, F_U32"\tsffErrors\n",   gkp->inf.sffErrors);
    fprintf(stdout, F_U32"\tsffWarnings\n", gkp->inf.sffWarnings);
    fprintf(stdout, "\n");
    fprintf(stdout, F_U32"\tsffLibCreated\n", gkp->inf.sffLibCreated);
    fprintf(stdout, "\n");
    fprintf(stdout, F_U32"\tplcInput\n",    gkp->inf.plcInput);
    fprintf(stdout, F_U32"\tplcLoaded\n",   gkp->inf.plcLoaded);
    fprintf(stdout, F_U32"\tplcErrors\n",   gkp->inf.plcErrors);
    fprintf(stdout, F_U32"\tplcWarnings\n", gkp->inf.plcWarnings);
    fprintf(stdout, "\n");
    fprintf(stdout, F_U32"\tnumRandom\n",   gkp->inf.numRandom);
    fprintf(stdout, F_U32"\tnumPacked\n",   gkp->inf.numPacked);
    fprintf(stdout, F_U32"\tnumNormal\n",   gkp->inf.numNormal);
    fprintf(stdout, F_U32"\tnumStrobe\n",   gkp->inf.numStrobe);
    fprintf(stdout, "\n");
    fprintf(stdout, "GLOBAL STATS\n");
    fprintf(stdout, "\n");
    fprintf(stdout, F_S32"\tFRG\n", (uint32)gkp->gkStore_getNumFragments());
    fprintf(stdout, F_S32"\tLIB\n", (uint32)gkp->gkStore_getNumLibraries());
    fprintf(stdout, "\n");
  }

  //  Header

  fprintf(stdout, "libIID\tbgnIID\tendIID\tactive\tdeleted\tmated\ttotLen\tclrLen\tlibName\n");

  //  Global

  fprintf(stdout, "0\t%d\t%d\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U64"\t"F_U64"\tGLOBAL\n",
          1, stats->numActiveFrag + stats->numDeletedFrag, stats->numActiveFrag, stats->numDeletedFrag, stats->numMatedFrag, stats->readLength, stats->clearLength);

  //  Per Library

  for (uint32 j=0; j<gkp->gkStore_getNumLibraries() + 1; j++)
    fprintf(stdout, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U64"\t"F_U64"\t%s\n",
            j,
            stats->lowestIID[j],
            stats->highestIID[j],
            stats->numActivePerLib[j],
            stats->numDeletedPerLib[j],
            stats->numMatedPerLib[j],
            stats->readLengthPerLib[j],
            stats->clearLengthPerLib[j],
            (j == 0 ? "LegacyUnmatedReads" : gkp->gkStore_getLibrary(j)->libraryName));

  delete stats;
  delete gkp;
}


void
dumpGateKeeperLibraries(char       *gkpStoreName,
                        AS_IID      bgnIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         asTable,
                        int         withoutUIDs) {
  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE, withoutUIDs);

  int       i;

  if (bgnIID < 1)
    bgnIID = 1;
  if (gkp->gkStore_getNumLibraries() < endIID)
    endIID = gkp->gkStore_getNumLibraries();

  if (asTable)
    fprintf(stdout, "UID\tIID\tOrientation\tMean\tStdDev\tNumFeatures\n");

  for (i=bgnIID; i<=endIID; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      gkLibrary      *gkpl = gkp->gkStore_getLibrary(i);
      LibraryMesg     lmesg;

      gkpl->gkLibrary_encodeFeatures(&lmesg);

      if (asTable) {
        fprintf(stdout, "%s\t"F_IID"\t%s\t%.3f\t%.3f\t%d\n",
                (withoutUIDs == TRUE ? "NA" : AS_UID_toString(gkpl->libraryUID)), i,
                AS_READ_ORIENT_NAMES[gkpl->orientation],
                gkpl->mean,
                gkpl->stddev,
                lmesg.num_features);
      } else {
        uint32 f;

        fprintf(stdout, "libraryIdent         = %s,"F_IID"\n", AS_UID_toString(gkpl->libraryUID), i);
        fprintf(stdout, "libraryOrientation   = %s\n", AS_READ_ORIENT_NAMES[gkpl->orientation]);
        fprintf(stdout, "libraryMean          = %.3f\n", gkpl->mean);
        fprintf(stdout, "libraryStdDev        = %.3f\n", gkpl->stddev);
        fprintf(stdout, "libraryNumFeatures   = %d\n", lmesg.num_features);

        for (f=0; f<lmesg.num_features; f++)
          fprintf(stdout, "libraryFeature[%d]    = %s=%s\n", f, lmesg.features[f], lmesg.values[f]);
      }

      gkpl->gkLibrary_encodeFeaturesCleanup(&lmesg);
    }
  }

  delete gkp;
}



void
dumpGateKeeperFragments(char       *gkpStoreName,
                        AS_IID      bgnIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         dumpWithSequence,
                        int         dumpClear,
                        int         asTable,
                        int         withoutUIDs) {
  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE, withoutUIDs);

  if (bgnIID < 1)
    bgnIID = 1;
  if (gkp->gkStore_getNumFragments() < endIID)
    endIID = gkp->gkStore_getNumFragments();

  gkFragment  fr;
  gkStream   *fs = new gkStream(gkp, bgnIID, endIID, (!dumpWithSequence || asTable) ? GKFRAGMENT_INF : GKFRAGMENT_QLT);

  //int           i;

  if (asTable)
    fprintf(stdout, "UID\tIID\tmateUID\tmateIID\tlibUID\tlibIID\tisDeleted\tisNonRandom\tOrient\tLength\tclrBegin%s\tclrEnd%s\n",
            AS_READ_CLEAR_NAMES[dumpClear], AS_READ_CLEAR_NAMES[dumpClear]);

  while (fs->next(&fr)) {
    if ((iidToDump == NULL) || (iidToDump[fr.gkFragment_getReadIID()])) {
      AS_IID     mateiid = fr.gkFragment_getMateIID();
      AS_UID     mateuid = AS_UID_undefined();

      if (mateiid > 0)
        mateuid = gkp->gkStore_getIIDtoUID(mateiid, AS_IID_FRG);

      AS_IID     libiid = fr.gkFragment_getLibraryIID();
      AS_UID     libuid = AS_UID_undefined();

      uint32     clrBgn = fr.gkFragment_getClearRegionBegin(dumpClear);
      uint32     clrEnd = fr.gkFragment_getClearRegionEnd  (dumpClear);

      if (libiid > 0)
        libuid = gkp->gkStore_getLibrary(libiid)->libraryUID;

      if (asTable) {
        fprintf(stdout, "%s\t"F_IID"\t%s\t"F_IID"\t%s\t"F_IID"\t%d\t%d\t%s\t%d\t%d\t%d\n",
                (withoutUIDs == TRUE ? "NA" : AS_UID_toString(fr.gkFragment_getReadUID())), fr.gkFragment_getReadIID(),
                (withoutUIDs == TRUE ? "NA" : AS_UID_toString(mateuid)), mateiid,
                (withoutUIDs == TRUE ? "NA" : AS_UID_toString(libuid)), libiid,
                fr.gkFragment_getIsDeleted(),
                fr.gkFragment_getIsNonRandom(),
                AS_READ_ORIENT_NAMES[fr.gkFragment_getOrientation()],
                fr.gkFragment_getSequenceLength(),
                clrBgn, clrEnd);
      } else {
        fprintf(stdout, "fragmentIdent           = %s,"F_IID"\n", AS_UID_toString(fr.gkFragment_getReadUID()), fr.gkFragment_getReadIID());
        fprintf(stdout, "fragmentMate            = %s,"F_IID"\n", AS_UID_toString(mateuid), mateiid);
        fprintf(stdout, "fragmentLibrary         = %s,"F_IID"\n", AS_UID_toString(libuid), libiid);

        fprintf(stdout, "fragmentIsDeleted       = %d\n", fr.gkFragment_getIsDeleted());
        fprintf(stdout, "fragmentIsNonRandom     = %d\n", fr.gkFragment_getIsNonRandom());
        fprintf(stdout, "fragmentOrientation     = %s\n", AS_READ_ORIENT_NAMES[fr.gkFragment_getOrientation()]);

        fprintf(stdout, "fragmentSeqLen          = %d\n", fr.gkFragment_getSequenceLength());

        fprintf(stdout, "fragmentClear           = %d,%d\n",
                fr.gkFragment_getClearRegionBegin(),
                fr.gkFragment_getClearRegionEnd());

        for (uint32 i=0; i<AS_READ_CLEAR_NUM; i++) {
          uint32 bgn, end;
          fr.gkFragment_getClearRegion(bgn, end, i);
          if (bgn < end)
            fprintf(stdout, "fragmentClear           = %s,%d,%d\n",
                    AS_READ_CLEAR_NAMES[i], bgn, end);
        }

        if (dumpWithSequence) {
          char          *seq = fr.gkFragment_getSequence();
          char          *qlt = fr.gkFragment_getQuality();

          for (uint32 i=0; i<clrBgn; i++)
            seq[i] = tolower(seq[i]);
          for (uint32 i=clrEnd; seq[i]; i++)
            seq[i] = tolower(seq[i]);

          fprintf(stdout, "fragmentSequence      = %s\n", seq);
          fprintf(stdout, "fragmentQuality       = %s\n", qlt);
        }

        fprintf(stdout, "fragmentSeqOffset       = "F_U64"\n", fr.gkFragment_getSequenceOffset());
        fprintf(stdout, "fragmentQltOffset       = "F_U64"\n", fr.gkFragment_getQualityOffset());
      }
    }
  }

  delete fs;
  delete gkp;
}



static
void
adjustBeginEndAddMates(gkStore *gkp,
                       AS_IID  &bgnIID,
                       AS_IID  &endIID,
                       int32  *&libToDump,
                       char   *&iidToDump,
                       AS_UID *&frgUID,
                       int      doNotFixMates,
                       int      dumpAllReads,
                       int      withoutUIDs) {

  if (bgnIID < 1)
    bgnIID = 1;

  if (gkp->gkStore_getNumFragments() < endIID)
    endIID = gkp->gkStore_getNumFragments();

  bool  iidToDumpIsPrecious = true;

  if (iidToDump == NULL) {
    iidToDump = (char *)safe_calloc(endIID + 1, sizeof(char));

    for (AS_IID i=bgnIID; i<=endIID; i++)
      iidToDump[i] = 1;

    iidToDumpIsPrecious = false;
  }

  libToDump = (int32  *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(int));

  if (withoutUIDs == 0)
    frgUID    = (AS_UID *)safe_calloc(gkp->gkStore_getNumFragments() + 1, sizeof(AS_UID));

  gkStream   *fs = new gkStream(gkp, bgnIID, endIID, GKFRAGMENT_INF);
  gkFragment  fr;

  int32       mateAdded = 0;

  fprintf(stderr, "Scanning store to find libraries used and reads to dump.\n");

  while (fs->next(&fr)) {
    AS_IID  iid = fr.gkFragment_getReadIID();
    AS_IID  mid = fr.gkFragment_getMateIID();

    if (withoutUIDs == 0)
      frgUID[iid] = fr.gkFragment_getReadUID();

    if ((iidToDumpIsPrecious == false) &&
        (dumpAllReads == false) &&
        (fr.gkFragment_getIsDeleted() == true))
      iidToDump[iid] = 0;
        
    if (iidToDump[iid]) {
      libToDump[fr.gkFragment_getLibraryIID()]++;

      if ((mid > 0) && (iidToDump[mid] == 0)) {
        mateAdded++;
        if (doNotFixMates == 0)
          iidToDump[fr.gkFragment_getMateIID()] = 1;
      }
    }
  }

  delete fs;

  fprintf(stderr, "%sdded %d reads to maintain mate relationships.\n",
          (doNotFixMates) ? "Would have a" : "A", mateAdded);

  fprintf(stderr, "Dumping %d fragments from unknown library (version 1 has these)\n",
          libToDump[0]);

  for (int32 i=1; i<=gkp->gkStore_getNumLibraries(); i++)
    fprintf(stderr, "Dumping %d fragments from library IID %d\n",
            libToDump[i], i);
}




void
dumpGateKeeperAsFasta(char       *gkpStoreName,
                      char       *prefix,
                      AS_IID      bgnIID,
                      AS_IID      endIID,
                      char       *iidToDump,
                      int         doNotFixMates,
                      int         dumpAllReads,
                      int         dumpAllBases,
                      int         dumpClear,
                      int         withoutUIDs) {

  gkStore   *gkp       = new gkStore(gkpStoreName, FALSE, FALSE);
  int32     *libToDump = NULL;
  AS_UID    *frgUID    = NULL;

  adjustBeginEndAddMates(gkp, bgnIID, endIID, libToDump, iidToDump, frgUID, doNotFixMates, dumpAllReads, withoutUIDs);

  char  fname[FILENAME_MAX];  sprintf(fname, "%s.fasta",      prefix);
  char  qname[FILENAME_MAX];  sprintf(qname, "%s.fasta.qv",   prefix);
  char  Qname[FILENAME_MAX];  sprintf(Qname, "%s.fasta.qual", prefix);

  errno = 0;
  FILE *f = fopen(fname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", fname, strerror(errno)), exit(1);

  FILE *q = fopen(qname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", qname, strerror(errno)), exit(1);

  FILE *Q = fopen(Qname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", Qname, strerror(errno)), exit(1);

  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.

  while ((bgnIID < endIID) && (iidToDump[bgnIID] == 0))
    bgnIID++;

  while ((bgnIID < endIID) && (iidToDump[endIID] == 0))
    endIID--;

  gkStream     *fs = new gkStream(gkp, bgnIID, endIID, GKFRAGMENT_QLT);
  gkFragment    fr;

  char          dl[1024];

  while (fs->next(&fr)) {
    int32   lclr = fr.gkFragment_getClearRegionBegin(dumpClear) + 1;
    int32   rclr = fr.gkFragment_getClearRegionEnd  (dumpClear);

    AS_IID  id1  = fr.gkFragment_getReadIID();
    AS_IID  id2  = fr.gkFragment_getMateIID();

    AS_IID  libIID = fr.gkFragment_getLibraryIID();
    AS_UID  libUID = (libIID > 0) ? gkp->gkStore_getLibrary(libIID)->libraryUID : AS_UID_undefined();

    if (iidToDump[id1] == 0) {
      //  Fragment isn't marked for dumping, don't dump.  Skip ahead to the next dumpable fragment.

      id2 = id1 + 1;
      while ((id2 < endIID) && (iidToDump[id2] == 0))
        id2++;

      if (id2 > id1 + 1) {
        //fprintf(stderr, "skip from "F_IID" to "F_IID".\n", id1, id2);
        fs->reset(id2, endIID);
      }

      continue;
    }

    if ((lclr >= rclr) && (dumpAllBases == false))
      //  Fragment has null or invalid clear range, don't dump.
      continue;

    AS_UID  fUID = AS_UID_undefined();
    AS_UID  mUID = AS_UID_undefined();

    if (withoutUIDs == 0) {
      fUID = fr.gkFragment_getReadUID();

      if (id2 > 0)
        mUID = gkp->gkStore_getIIDtoUID(id2, AS_IID_FRG);
    }

    char *seq = fr.gkFragment_getSequence() + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(dumpClear) : 0);
    char *qlt = fr.gkFragment_getQuality()  + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(dumpClear) : 0);

    int32 len = (dumpAllBases == false) ? fr.gkFragment_getClearRegionLength(dumpClear) : fr.gkFragment_getSequenceLength();

    seq[len] = 0;
    qlt[len] = 0;

    //  Lowercase sequence not in the clear range.

    if (fr.gkFragment_getIsDeleted()) {
        if (dumpAllBases == true) {
            for (int i=0; seq[i]; i++)
                seq[i] = tolower(seq[i]);
        } else {
            continue;
        }
    }

    //  If allBases == true, seq is the whole read.  If false, seq is limited to just the clear
    //  bases.

    if (dumpAllBases == true) {
      for (int i=0; i<=lclr; i++)
        seq[i] = tolower(seq[i]); 
      for (int i=rclr; seq[i]; i++)
        seq[i] = tolower(seq[i]);
    }

    AS_UTL_writeFastA(f, seq, len, 0,
                      ">%s,"F_IID" mate=%s,"F_IID" lib=%s,"F_IID" clr=%s,%d,%d deleted=%d\n",
                      AS_UID_toString(fUID), id1,
                      AS_UID_toString(mUID), id2,
                      AS_UID_toString(libUID), libIID,
                      AS_READ_CLEAR_NAMES[dumpClear], lclr, rclr,
                      fr.gkFragment_getIsDeleted());

    AS_UTL_writeFastA(q, qlt, len, 0,
                      ">%s,"F_IID" mate=%s,"F_IID" lib=%s,"F_IID" clr=%s,%d,%d deleted=%d\n",
                      AS_UID_toString(fUID), id1,
                      AS_UID_toString(mUID), id2,
                      AS_UID_toString(libUID), libIID,
                      AS_READ_CLEAR_NAMES[dumpClear], lclr, rclr,
                      fr.gkFragment_getIsDeleted());

    AS_UTL_writeQVFastA(Q, qlt, len, 0,
                        ">%s,"F_IID" mate=%s,"F_IID" lib=%s,"F_IID" clr=%s,%d,%d deleted=%d\n",
                        AS_UID_toString(fUID), id1,
                        AS_UID_toString(mUID), id2,
                        AS_UID_toString(libUID), libIID,
                        AS_READ_CLEAR_NAMES[dumpClear], lclr, rclr,
                        fr.gkFragment_getIsDeleted());
  }

  delete fs;
  delete gkp;
}



void
dumpGateKeeperAsFRG(char       *gkpStoreName,
                    int         dumpFormat,
                    AS_IID      bgnIID,
                    AS_IID      endIID,
                    char       *iidToDump,
                    int         doNotFixMates,
                    int         dumpAllReads,
                    int         dumpFRGClear,
                    int         withoutUIDs) {
  gkFragment        fr;
  gkStream       *fs = NULL;

  unsigned int      firstElem = 0;
  unsigned int      lastElem = 0;

  GenericMesg       pmesg;

  int               i;

  int              *libToDump;
  AS_UID           *libUID;
  AS_UID           *frgUID;
  int               mateAdded = 0;

  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE);

  //  Someone could adjust the iidToDump and frgUID arrays to be
  //  relative to bgnIID.

  if (bgnIID < 1)
    bgnIID = 1;
  if (gkp->gkStore_getNumFragments() < endIID)
    endIID = gkp->gkStore_getNumFragments();

  if (iidToDump == NULL) {
    iidToDump = (char *)safe_calloc(gkp->gkStore_getNumFragments() + 1, sizeof(char));
    for (i=bgnIID; i<=endIID; i++)
      iidToDump[i] = 1;
  }

  //  Pass 1: Scan the fragments, build a list of libraries to dump,
  //  also add in any mates that were omitted from the input.

  fprintf(stderr, "Scanning store to find libraries used.\n");

  libToDump = (int    *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(int));
  libUID    = (AS_UID *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(AS_UID));
  frgUID    = (AS_UID *)safe_calloc(gkp->gkStore_getNumFragments() + 1, sizeof(AS_UID));

  fs = new gkStream(gkp, bgnIID, endIID, GKFRAGMENT_INF);

  while (fs->next(&fr)) {
    if ((iidToDump[fr.gkFragment_getReadIID()]) &&
        (dumpAllReads || !fr.gkFragment_getIsDeleted())) {
      AS_IID frgIID  = fr.gkFragment_getReadIID();
      frgUID[frgIID] = fr.gkFragment_getReadUID();

      libToDump[fr.gkFragment_getLibraryIID()]++;

      if ((fr.gkFragment_getMateIID() > 0) && (iidToDump[fr.gkFragment_getMateIID()] == 0)) {
        mateAdded++;
        if (doNotFixMates == 0) {
          // pull the mate fragment so we can get its UID and store it in our lookup table (needed for LKG record)
          gkp->gkStore_getFragment(fr.gkFragment_getMateIID(), &fr, GKFRAGMENT_INF);
          iidToDump[fr.gkFragment_getReadIID()] = 1;
          frgUID[fr.gkFragment_getReadIID()] = fr.gkFragment_getReadUID();
        }
      }
      if (gkp->gkStore_getFRGtoPLC(frgIID) != 0) {
         // pull the placement fragment information so we can get their UIDs and store in the lookup table (needed for PLC record)
         gkPlacement *gkpl = gkp->gkStore_getReadPlacement(frgIID);
         
         gkp->gkStore_getFragment(gkpl->bound1, &fr, GKFRAGMENT_INF);
         iidToDump[fr.gkFragment_getReadIID()] = 1;
         frgUID[fr.gkFragment_getReadIID()] = fr.gkFragment_getReadUID();
         
         gkp->gkStore_getFragment(gkpl->bound2, &fr, GKFRAGMENT_INF);
         iidToDump[fr.gkFragment_getReadIID()] = 1;
         frgUID[fr.gkFragment_getReadIID()] = fr.gkFragment_getReadUID();
      }
    }
  }

  fprintf(stderr, "%sdded %d reads to maintain mate relationships.\n",
          (doNotFixMates) ? "Would have a" : "A", mateAdded);

  fprintf(stderr, "Dumping %d fragments from unknown library (version 1 has these)\n", libToDump[0]);

  for (i=1; i<=gkp->gkStore_getNumLibraries(); i++)
    fprintf(stderr, "Dumping %d fragments from library IID %d\n", libToDump[i], i);

  //  Dump the format message
  //
  if (dumpFormat > 1) {
    AS_MSG_setFormatVersion(dumpFormat);

    VersionMesg  vmesg;
    vmesg.version = dumpFormat;

    pmesg.m = &vmesg;
    pmesg.t = MESG_VER;

    WriteProtoMesg_AS(stdout, &pmesg);
  }


  //  Dump libraries.
  //
  for (i=1; i<=gkp->gkStore_getNumLibraries(); i++) {
    if (libToDump[i]) {
      DistanceMesg              dmesg;
      LibraryMesg               lmesg;
      gkLibrary  *gkpl = gkp->gkStore_getLibrary(i);

      //  We don't really need to cache UIDs anymore;
      //  gkStore_getLibrary should be doing this for us.  It's
      //  already implemented, and slightly more efficient, so BPW
      //  left it in.

      libUID[i] = gkpl->libraryUID;

      if (dumpFormat == 1) {
        pmesg.m = &dmesg;
        pmesg.t = MESG_DST;

        dmesg.action     = AS_ADD;
        dmesg.eaccession = gkpl->libraryUID;
        dmesg.mean       = gkpl->mean;
        dmesg.stddev     = gkpl->stddev;
      } else {
        pmesg.m = &lmesg;
        pmesg.t = MESG_LIB;

        lmesg.action       = AS_ADD;
        lmesg.eaccession   = gkpl->libraryUID;
        lmesg.mean         = gkpl->mean;
        lmesg.stddev       = gkpl->stddev;
        lmesg.source       = NULL;

        lmesg.link_orient.setIsUnknown();
        switch(gkpl->orientation) {
          case AS_READ_ORIENT_INNIE:
            lmesg.link_orient.setIsInnie();
            break;
          case AS_READ_ORIENT_OUTTIE:
            lmesg.link_orient.setIsOuttie();
            break;
          case AS_READ_ORIENT_NORMAL:
            lmesg.link_orient.setIsNormal();
            break;
          case AS_READ_ORIENT_ANTINORMAL:
            lmesg.link_orient.setIsAnti();
            break;
          case AS_READ_ORIENT_UNKNOWN:
            lmesg.link_orient.setIsUnknown();
            break;
          default:
            //  Cannot happen, unless someone adds a new orientation to gkFragment.
            assert(0);
            break;
        }

        gkpl->gkLibrary_encodeFeatures(&lmesg);
      }

      WriteProtoMesg_AS(stdout, &pmesg);

      if (dumpFormat != 1)
        gkpl->gkLibrary_encodeFeaturesCleanup(&lmesg);
    }
  }
  delete fs;

  //  Load clear regions, if they exist.
  //
  gkClearRange  *clrVec = new gkClearRange(gkp, AS_READ_CLEAR_VEC);
  gkClearRange  *clrMax = new gkClearRange(gkp, AS_READ_CLEAR_MAX);

  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.
  //
  //  Also, dump any constraints the fragments have
  //
  fs = new gkStream(gkp, bgnIID, endIID, GKFRAGMENT_QLT);

  while (fs->next(&fr)) {
    FragMesg      fmesg;
    LinkMesg      lmesg;
    PlacementMesg cmesg;

    if ((iidToDump[fr.gkFragment_getReadIID()]) &&
        (dumpAllReads || !fr.gkFragment_getIsDeleted())) {
      pmesg.m = &fmesg;
      pmesg.t = MESG_FRG;

      //  This code used in AS_GKP_dump.c (dumpFRG), and in AS_FGB_io.c
      fmesg.action              = fr.gkFragment_getIsDeleted() ? AS_DELETE : AS_ADD;
      fmesg.eaccession          = fr.gkFragment_getReadUID();
      fmesg.library_uid         = libUID[fr.gkFragment_getLibraryIID()];
      fmesg.library_iid         = fr.gkFragment_getLibraryIID();
      fmesg.plate_uid           = AS_UID_undefined();  // fr.gkFragment_getPlateUID();
      fmesg.plate_location      = 0;                   // fr.gkFragment_getPlateLocation();
      fmesg.type                = AS_READ;
      fmesg.is_random           = (fr.gkFragment_getIsNonRandom()) ? 0 : 1;
      fmesg.status_code         = 'X';
      fmesg.clear_rng.bgn       = fr.gkFragment_getClearRegionBegin();
      fmesg.clear_rng.end       = fr.gkFragment_getClearRegionEnd  ();
      fmesg.clear_vec.bgn       = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC);
      fmesg.clear_vec.end       = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC);
      fmesg.clear_max.bgn       = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX);
      fmesg.clear_max.end       = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX);
      fmesg.contamination.bgn   = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT);
      fmesg.contamination.end   = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT);
      fmesg.source              = NULL; //fr.gkFragment_getSource();
      fmesg.sequence            = fr.gkFragment_getSequence();
      fmesg.quality             = fr.gkFragment_getQuality();
      fmesg.hps                 = NULL; //fr.gkFragment_getHPS();
      fmesg.iaccession          = firstElem;

      WriteProtoMesg_AS(stdout, &pmesg);

      if ((fr.gkFragment_getMateIID() > 0) &&
          (fr.gkFragment_getMateIID() < fr.gkFragment_getReadIID())) {
        pmesg.m = &lmesg;
        pmesg.t = MESG_LKG;

        lmesg.action      = AS_ADD;
        lmesg.type.setIsMatePair();
        lmesg.link_orient.setIsUnknown();
        lmesg.frag1       = frgUID[fr.gkFragment_getMateIID()];
        lmesg.frag2       = fr.gkFragment_getReadUID();
        lmesg.distance    = libUID[fr.gkFragment_getLibraryIID()];

        lmesg.link_orient.setIsUnknown();

        switch(fr.gkFragment_getOrientation()) {
          case AS_READ_ORIENT_INNIE:
            lmesg.link_orient.setIsInnie();
            break;
          case AS_READ_ORIENT_OUTTIE:
            lmesg.link_orient.setIsOuttie();
            break;
          case AS_READ_ORIENT_NORMAL:
            lmesg.link_orient.setIsNormal();
            break;
          case AS_READ_ORIENT_ANTINORMAL:
            lmesg.link_orient.setIsAnti();
            break;
          case AS_READ_ORIENT_UNKNOWN:
            lmesg.link_orient.setIsUnknown();
            break;
          default:
            //  Cannot happen, unless someone adds a new orientation to gkFragment.
            assert(0);
            break;
        }

        WriteProtoMesg_AS(stdout, &pmesg);
      }
      
      gkPlacement *gkpl = gkp->gkStore_getReadPlacement(fr.gkFragment_getReadIID());
      if (gkpl != NULL && dumpFormat != 1) {
         pmesg.m = &cmesg;
         pmesg.t = MESG_PLC;
         
         cmesg.action  = AS_ADD;
         cmesg.frag    = frgUID[gkpl->frag];
         cmesg.bound1  = frgUID[gkpl->bound1];
         cmesg.bound2  = frgUID[gkpl->bound2];
         
         WriteProtoMesg_AS(stdout, &pmesg);
      }
    }
  }

  delete fs;
  delete gkp;
}



void
dumpGateKeeperAsNewbler(char       *gkpStoreName,
                        char       *prefix,
                        AS_IID      bgnIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         doNotFixMates,
                        int         dumpAllReads,
                        int         dumpAllBases,
                        int         dumpClear,
                        int         withoutUIDs) {

  gkStore   *gkp       = new gkStore(gkpStoreName, FALSE, FALSE);
  int32     *libToDump = NULL;
  AS_UID    *frgUID    = NULL;

  adjustBeginEndAddMates(gkp, bgnIID, endIID, libToDump, iidToDump, frgUID, doNotFixMates, dumpAllReads, withoutUIDs);

  char  fname[FILENAME_MAX];  sprintf(fname, "%s.fna",      prefix);
  char  qname[FILENAME_MAX];  sprintf(qname, "%s.fna.qual", prefix);

  errno = 0;
  FILE *f = fopen(fname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", fname, strerror(errno)), exit(1);

  FILE *q = fopen(qname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", qname, strerror(errno)), exit(1);

  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.

  gkStream     *fs = new gkStream(gkp, bgnIID, endIID, GKFRAGMENT_QLT);
  gkFragment    fr;

  char          dl[1024];

  while (fs->next(&fr)) {
    int32   lclr = fr.gkFragment_getClearRegionBegin(dumpClear) + 1;
    int32   rclr = fr.gkFragment_getClearRegionEnd  (dumpClear);

    AS_IID  id1  = fr.gkFragment_getReadIID();
    AS_IID  id2  = fr.gkFragment_getMateIID();

    AS_IID  libIID = fr.gkFragment_getLibraryIID();
    AS_UID  libUID = (libIID > 0) ? gkp->gkStore_getLibrary(libIID)->libraryUID : AS_UID_undefined();

    if (iidToDump[id1] == 0) {
      //  Fragment isn't marked for dumping, don't dump.
      continue;
    }

    if ((lclr >= rclr) && (dumpAllBases == false))
      //  Fragment has null or invalid clear range, don't dump.
      continue;

    if (id2 > 0)
      sprintf(dl, ">%s template=%d+%d dir=%c library=%s trim=%d-%d\n",
              //  ID
              AS_UID_toString(fr.gkFragment_getReadUID()),
              //  template
              (id1 < id2) ? id1 : id2,
              (id1 < id2) ? id2 : id1,
              //  dir
              (id1 < id2) ? 'F' : 'R',
              //  library
              AS_UID_toString(libUID),
              //  trim
              lclr,
              rclr);

    else
      sprintf(dl, ">%s trim=%d-%d\n",
              //  ID
              AS_UID_toString(fr.gkFragment_getReadUID()),
              //  trim
              lclr,
              rclr);

    AS_UTL_writeFastA(f,
                      fr.gkFragment_getSequence(),
                      fr.gkFragment_getSequenceLength(), 0,
                      dl, NULL);

    AS_UTL_writeQVFastA(q,
                        fr.gkFragment_getQuality(),
                        fr.gkFragment_getSequenceLength(), 0,
                        dl, NULL);
  }

  delete fs;
  delete gkp;
}



void
dumpGateKeeperAsFastQ(char       *gkpStoreName,
                      char       *prefix,
                      AS_IID      bgnIID,
                      AS_IID      endIID,
                      char       *iidToDump,
                      int         doNotFixMates,
                      int         dumpAllReads,
                      int         dumpAllBases,
                      int         dumpClear,
                      int         withoutUIDs) {

  gkStore   *gkp       = new gkStore(gkpStoreName, FALSE, FALSE);
  int32     *libToDump = NULL;
  AS_UID    *frgUID    = NULL;

  adjustBeginEndAddMates(gkp, bgnIID, endIID, libToDump, iidToDump, frgUID, doNotFixMates, dumpAllReads, withoutUIDs);

  char  aname[FILENAME_MAX];  sprintf(aname, "%s.1.fastq",       prefix);
  char  bname[FILENAME_MAX];  sprintf(bname, "%s.2.fastq",       prefix);
  char  pname[FILENAME_MAX];  sprintf(pname, "%s.paired.fastq",  prefix);
  char  uname[FILENAME_MAX];  sprintf(uname, "%s.unmated.fastq", prefix);

  errno = 0;
  FILE *a = fopen(aname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", aname, strerror(errno)), exit(1);

  errno = 0;
  FILE *b = fopen(bname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", bname, strerror(errno)), exit(1);

  errno = 0;
  FILE *p = fopen(pname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", pname, strerror(errno)), exit(1);

  FILE *u = fopen(uname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", uname, strerror(errno)), exit(1);

  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.

  while ((bgnIID < endIID) && (iidToDump[bgnIID] == 0))
    bgnIID++;

  while ((bgnIID < endIID) && (iidToDump[endIID] == 0))
    endIID--;

  gkStream   *fs = new gkStream(gkp, bgnIID, endIID, GKFRAGMENT_QLT);
  gkFragment  fr;

  while (fs->next(&fr)) {
    int32   lclr   = fr.gkFragment_getClearRegionBegin(dumpClear) + 1;
    int32   rclr   = fr.gkFragment_getClearRegionEnd  (dumpClear);

    AS_IID  id1    = fr.gkFragment_getReadIID();
    AS_IID  id2    = fr.gkFragment_getMateIID();

    AS_IID  libIID = fr.gkFragment_getLibraryIID();
    AS_UID  libUID = (libIID > 0) ? gkp->gkStore_getLibrary(libIID)->libraryUID : AS_UID_undefined();

    bool    flipReads = gkp->gkStore_getLibrary(libIID)->readsAreReversed;

    if (iidToDump[id1] == 0) {
      //  Fragment isn't marked for dumping, don't dump.  Skip ahead to the next dumpable fragment.

      id2 = id1 + 1;
      while ((id2 < endIID) && (iidToDump[id2] == 0))
        id2++;

      if (id2 > id1 + 1) {
        //fprintf(stderr, "skip from "F_IID" to "F_IID".\n", id1, id2);
        fs->reset(id2, endIID);
      }

      continue;
    }

    if ((lclr >= rclr) && (dumpAllBases == false))
      //  Fragment has null or invalid clear range, don't dump.
      continue;

    if ((id2 != 0) && (id2 < id1))
      //  Mated, and the mate is the first frag.  We've already reported this one.
      continue;

    char *seq = fr.gkFragment_getSequence() + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(dumpClear) : 0);
    char *qlt = fr.gkFragment_getQuality()  + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(dumpClear) : 0);

    int32 len = (dumpAllBases == false) ? fr.gkFragment_getClearRegionLength(dumpClear) : fr.gkFragment_getSequenceLength();

    seq[len] = 0;
    qlt[len] = 0;

    uint32  clrBgn = fr.gkFragment_getClearRegionBegin(dumpClear);
    uint32  clrEnd = fr.gkFragment_getClearRegionEnd  (dumpClear);

    uint32  vecBgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC);
    uint32  vecEnd = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC);

    uint32  maxBgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX);
    uint32  maxEnd = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX);

    uint32  tntBgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT);
    uint32  tntEnd = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT);

    //  If dumping only the clear range, reset clear ranges.
    //
    if (dumpAllBases == false) {
      if (vecBgn < vecEnd) {
        vecBgn = (vecBgn < clrBgn) ?   0 : (vecBgn - clrBgn);
        vecEnd = (vecEnd < clrBgn) ?   0 : (vecEnd - clrBgn);
        vecEnd = (len    < vecEnd) ? len : (vecEnd);
      }

      if (maxBgn < maxEnd) {
        maxBgn = (maxBgn < clrBgn) ?   0 : (maxBgn - clrBgn);
        maxEnd = (maxEnd < clrBgn) ?   0 : (maxEnd - clrBgn);
        maxEnd = (len    < maxEnd) ? len : (maxEnd);
      }

      if (tntBgn < tntEnd) {
        tntBgn = (tntBgn < clrBgn) ?   0 : (tntBgn - clrBgn);
        tntEnd = (tntEnd < clrBgn) ?   0 : (tntEnd - clrBgn);
        tntEnd = (len    < tntEnd) ? len : (tntEnd);
      }

      clrBgn = 0;
      clrEnd = len;

      //  The easiest fix is to just ignore clear ranges.
      //vecBgn = maxBgn = tntBgn = 1;
      //vecEnd = maxEnd = tntEnd = 0;
    }

    //  If flipping, do the flip and reset clear ranges.
    //
    if (flipReads) {
      uint32 bgn;

      reverseComplement(seq, qlt, len);

      if (clrBgn < clrEnd) {
        bgn    = clrBgn;
        clrBgn = len - clrEnd;
        clrEnd = len - bgn;
      }

      if (vecBgn < vecEnd) {
        bgn    = vecBgn;
        vecBgn = len - vecEnd;
        vecEnd = len - bgn;
      }

      if (maxBgn < maxEnd) {
        bgn    = maxBgn;
        maxBgn = len - maxEnd;
        maxEnd = len - bgn;
      }

      if (tntBgn < tntEnd) {
        bgn    = tntBgn;
        tntBgn = len - tntEnd;
        tntEnd = len - bgn;
      }
    }

    if (id2 == 0) {
      //  Unmated read, dump to the unmated reads file.
      //
      AS_UTL_writeFastQ(u, seq, len, qlt, len,
                        "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                        AS_UID_toString(fr.gkFragment_getReadUID()),
                        clrBgn, clrEnd, vecBgn, vecEnd, maxBgn, maxEnd, tntBgn, tntEnd,
                        fr.gkFragment_getIsNonRandom() ? 'f' : 't');
      continue;
    }


    //  We must now have a valid mate pair.  Note that a after we dump,
    //  the iidToDump[] array is updated so we never see this pair again.


    //  Write the first fragment (twice).
    AS_UTL_writeFastQ(a, seq, len, qlt, len,
                      "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      AS_UID_toString(fr.gkFragment_getReadUID()),
                      clrBgn, clrEnd, vecBgn, vecEnd, maxBgn, maxEnd, tntBgn, tntEnd,
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    AS_UTL_writeFastQ(p, seq, len, qlt, len,
                      "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      AS_UID_toString(fr.gkFragment_getReadUID()),
                      clrBgn, clrEnd, vecBgn, vecEnd, maxBgn, maxEnd, tntBgn, tntEnd,
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    //  Grab the second fragment.

    gkp->gkStore_getFragment(id2, &fr, GKFRAGMENT_QLT);

    lclr = fr.gkFragment_getClearRegionBegin(dumpClear) + 1;
    rclr = fr.gkFragment_getClearRegionEnd  (dumpClear);

    seq = fr.gkFragment_getSequence() + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(dumpClear) : 0);
    qlt = fr.gkFragment_getQuality()  + ((dumpAllBases == false) ? fr.gkFragment_getClearRegionBegin(dumpClear) : 0);
    len = (dumpAllBases == false) ? fr.gkFragment_getClearRegionLength(dumpClear) : fr.gkFragment_getSequenceLength();

    seq[len] = 0;
    qlt[len] = 0;

    clrBgn = fr.gkFragment_getClearRegionBegin(dumpClear);
    clrEnd = fr.gkFragment_getClearRegionEnd  (dumpClear);

    vecBgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC);
    vecEnd = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC);

    maxBgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX);
    maxEnd = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX);

    tntBgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT);
    tntEnd = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT);

    //  If dumping only the clear range, reset clear ranges.
    //
    if (dumpAllBases == false) {
      if (vecBgn < vecEnd) {
        vecBgn = (vecBgn < clrBgn) ?   0 : (vecBgn - clrBgn);
        vecEnd = (vecEnd < clrBgn) ?   0 : (vecEnd - clrBgn);
        vecEnd = (len    < vecEnd) ? len : (vecEnd);
      }

      if (maxBgn < maxEnd) {
        maxBgn = (maxBgn < clrBgn) ?   0 : (maxBgn - clrBgn);
        maxEnd = (maxEnd < clrBgn) ?   0 : (maxEnd - clrBgn);
        maxEnd = (len    < maxEnd) ? len : (maxEnd);
      }

      if (tntBgn < tntEnd) {
        tntBgn = (tntBgn < clrBgn) ?   0 : (tntBgn - clrBgn);
        tntEnd = (tntEnd < clrBgn) ?   0 : (tntEnd - clrBgn);
        tntEnd = (len    < tntEnd) ? len : (tntEnd);
      }

      clrBgn = 0;
      clrEnd = len;

      //  The easiest fix is to just ignore clear ranges.
      //vecBgn = maxBgn = tntBgn = 1;
      //vecEnd = maxEnd = tntEnd = 0;
    }

    //  If flipping, do the flip and reset clear ranges.
    //
    if (flipReads) {
      uint32 bgn;

      reverseComplement(seq, qlt, len);

      if (clrBgn < clrEnd) {
        bgn    = clrBgn;
        clrBgn = len - clrEnd;
        clrEnd = len - bgn;
      }

      if (vecBgn < vecEnd) {
        bgn    = vecBgn;
        vecBgn = len - vecEnd;
        vecEnd = len - bgn;
      }

      if (maxBgn < maxEnd) {
        bgn    = maxBgn;
        maxBgn = len - maxEnd;
        maxEnd = len - bgn;
      }

      if (tntBgn < tntEnd) {
        bgn    = tntBgn;
        tntBgn = len - tntEnd;
        tntEnd = len - bgn;
      }
    }

    //  Write the second fragment (twice).
    AS_UTL_writeFastQ(b, seq, len, qlt, len,
                      "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      AS_UID_toString(fr.gkFragment_getReadUID()),
                      clrBgn, clrEnd, vecBgn, vecEnd, maxBgn, maxEnd, tntBgn, tntEnd,
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    AS_UTL_writeFastQ(p, seq, len, qlt, len,
                      "@%s clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      AS_UID_toString(fr.gkFragment_getReadUID()),
                      clrBgn, clrEnd, vecBgn, vecEnd, maxBgn, maxEnd, tntBgn, tntEnd,
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    //  Mark the pair as dumped.  This is a cheap way around testing if we've already dumped a pair.
    //
    iidToDump[id1] = 0;
    iidToDump[id2] = 0;
  }

  delete fs;
  delete gkp;
}



int
dumpGateKeeperIsFeatureSet(char      *gkpStoreName,
                           AS_IID     libIID,
                           char      *featureName)
{
   assert(featureName != NULL);
  
   gkStore *gkp = new gkStore(gkpStoreName, FALSE, FALSE, TRUE);
   int isSet = 0;
   uint32 i;
   
   uint32 bgnIID, endIID;
   bgnIID = endIID = libIID;
   if (libIID <= 0 || libIID > gkp->gkStore_getNumLibraries()) {   
      bgnIID = 1;
      endIID = gkp->gkStore_getNumLibraries();
   }

   for (i=bgnIID; i<=endIID; i++) {
      gkLibrary      *gkpl = gkp->gkStore_getLibrary(i);
      LibraryMesg     lmesg;
      uint32          f;
      gkpl->gkLibrary_encodeFeatures(&lmesg);

      for (f=0; f<lmesg.num_features; f++)
         if (strcasecmp(lmesg.features[f], featureName) == 0)
            isSet |= atoi(lmesg.values[f]);

      gkpl->gkLibrary_encodeFeaturesCleanup(&lmesg);
   }
   
   delete gkp;
   return isSet;
}
