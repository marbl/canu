
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

static char const *rcsid = "$Id: AS_GKP_dump.c,v 1.48 2009-06-10 18:05:13 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_GKP_include.h"
#include "AS_UTL_fasta.h"
#include "AS_PER_gkpStore.h"

void
dumpGateKeeperInfo(char       *gkpStoreName,
                   int         asTable) {

  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

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
    fprintf(stdout, F_U32"\tnumRandom\n",   gkp->inf.numRandom);
    fprintf(stdout, F_U32"\tnumShort\n",    gkp->inf.numShort);
    fprintf(stdout, F_U32"\tnumMedium\n",   gkp->inf.numMedium);
    fprintf(stdout, F_U32"\tnumLong\n",     gkp->inf.numLong);
  }

  gkFragment    fr;
  gkStream     *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  int           i, j;

  uint32  numActiveFrag     = 0;
  uint32  numDeletedFrag    = 0;
  uint32  numMatedFrag      = 0;
  uint64  readLength        = 0;
  uint64  clearLength       = 0;

  uint32  *numActivePerLib   = (uint32 *)safe_calloc(sizeof(uint32), (gkp->gkStore_getNumLibraries() + 1));
  uint32  *numDeletedPerLib  = (uint32 *)safe_calloc(sizeof(uint32), (gkp->gkStore_getNumLibraries() + 1));
  uint32  *numMatedPerLib    = (uint32 *)safe_calloc(sizeof(uint32), (gkp->gkStore_getNumLibraries() + 1));
  uint64  *readLengthPerLib  = (uint64 *)safe_calloc(sizeof(uint64), (gkp->gkStore_getNumLibraries() + 1));
  uint64  *clearLengthPerLib = (uint64 *)safe_calloc(sizeof(uint64), (gkp->gkStore_getNumLibraries() + 1));

  while (fs->next(&fr)) {
    AS_IID     lib = fr.gkFragment_getLibraryIID();

    if (fr.gkFragment_getIsDeleted()) {
      numDeletedFrag++;
      numDeletedPerLib[lib]++;
    } else {
      numActiveFrag++;
      numActivePerLib[lib]++;

      if (fr.gkFragment_getMateIID() > 0) {
        numMatedFrag++;
        numMatedPerLib[lib]++;
      }

      readLength             += fr.gkFragment_getSequenceLength();
      readLengthPerLib[lib]  += fr.gkFragment_getSequenceLength();

      clearLength            += fr.gkFragment_getClearRegionLength();
      clearLengthPerLib[lib] += fr.gkFragment_getClearRegionLength();
    }
  }

  if (asTable == 0) {
    fprintf(stdout, "\n");
    fprintf(stdout, "GLOBAL STATS\n");
    fprintf(stdout, "\n");
    fprintf(stdout, F_S32"\tFRG\n", (uint32)gkp->gkStore_getNumFragments());
    fprintf(stdout, F_S32"\tLIB\n", (uint32)gkp->gkStore_getNumLibraries());
    fprintf(stdout, "\n");
  }

  //  Header

  fprintf(stdout, "LibraryName\tnumActiveFRG\tnumDeletedFRG\tnumMatedFRG\treadLength\tclearLength");

  //  Global

  fprintf(stdout, "GLOBAL\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U64"\t"F_U64"\n",
          numActiveFrag, numDeletedFrag, numMatedFrag, readLength, clearLength);

  //  Per Library

  for (j=0; j<gkp->gkStore_getNumLibraries() + 1; j++) {
    fprintf(stdout, "%s\t"F_U32"\t"F_U32"\t"F_U32,
            (j == 0) ? "LegacyUnmatedReads" : AS_UID_toString(gkp->gkStore_getLibrary(j)->libraryUID),
            numActivePerLib[j], numDeletedPerLib[j], numMatedPerLib[j], readLengthPerLib[j], clearLengthPerLib[j]);
  }

  safe_free(numActivePerLib);
  safe_free(numDeletedPerLib);
  safe_free(numMatedPerLib);
  safe_free(readLengthPerLib);
  safe_free(clearLengthPerLib);

  delete fs;
  delete gkp;
}


void
dumpGateKeeperLibraries(char       *gkpStoreName,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         asTable) {
  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE);

  int       i;

  if (begIID < 1)
    begIID = 1;
  if (gkp->gkStore_getNumLibraries() < endIID)
    endIID = gkp->gkStore_getNumLibraries();

  if (asTable)
    fprintf(stdout, "UID\tIID\tOrientation\tMean\tStdDev\tNumFeatures\n");

  for (i=begIID; i<=endIID; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      gkLibrary      *gkpl = gkp->gkStore_getLibrary(i);
      LibraryMesg     lmesg;
      int              nf;

      gkpl->gkLibrary_encodeFeatures(&lmesg);

      if (asTable) {
        fprintf(stdout, "%s\t"F_IID"\t%s\t%.3f\t%.3f\t%d\n",
                AS_UID_toString(gkpl->libraryUID), i,
                AS_READ_ORIENT_NAMES[gkpl->orientation],
                gkpl->mean,
                gkpl->stddev,
                nf);
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
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         dumpWithSequence,
                        int         dumpClear,
                        int         asTable) {
  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE);

  if (begIID < 1)
    begIID = 1;
  if (gkp->gkStore_getNumFragments() < endIID)
    endIID = gkp->gkStore_getNumFragments();

  gkFragment  fr;
  gkStream   *fs = new gkStream(gkp, begIID, endIID, (!dumpWithSequence || asTable) ? GKFRAGMENT_INF : GKFRAGMENT_QLT);

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
                AS_UID_toString(fr.gkFragment_getReadUID()), fr.gkFragment_getReadIID(),
                AS_UID_toString(mateuid), mateiid,
                AS_UID_toString(libuid), libiid,
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


void
dumpGateKeeperAsFasta(char       *gkpStoreName,
                      AS_IID      begIID,
                      AS_IID      endIID,
                      char       *iidToDump,
                      int         dumpAllReads,
                      int         dumpClear,
                      int         dumpQuality) {
  gkStore      *gkp = new gkStore(gkpStoreName, FALSE, FALSE);

  int           i;

  if (begIID < 1)
    begIID = 1;
  if (gkp->gkStore_getNumFragments() < endIID)
    endIID = gkp->gkStore_getNumFragments();

  gkFragment    fr;
  gkStream     *fs  = new gkStream(gkp, begIID, endIID, (dumpQuality) ? GKFRAGMENT_QLT : GKFRAGMENT_SEQ);

  while (fs->next(&fr)) {
    if ((iidToDump == NULL) || (iidToDump[fr.gkFragment_getReadIID()])) {
      if (dumpAllReads || !fr.gkFragment_getIsDeleted()) {
        AS_IID     mateiid = fr.gkFragment_getMateIID();
        AS_UID     mateuid = {0};

        if (mateiid > 0)
          mateuid = gkp->gkStore_getIIDtoUID(mateiid, AS_IID_FRG);

        AS_IID     libiid = fr.gkFragment_getLibraryIID();
        AS_UID     libuid = {0};

        if (libiid > 0)
          libuid = gkp->gkStore_getLibrary(libiid)->libraryUID;

        unsigned int   clrBeg   = fr.gkFragment_getClearRegionBegin();
        unsigned int   clrEnd   = fr.gkFragment_getClearRegionEnd  ();
        char          *seqStart = fr.gkFragment_getSequence();

        if (dumpQuality)
          seqStart = fr.gkFragment_getQuality();

        char          *seq      = seqStart+clrBeg;
        seq[clrEnd] = 0;

        if (dumpQuality >=2) {
          AS_UTL_writeQVFastA(stdout, seq, clrEnd-clrBeg,
                              ">%s,"F_IID" mate=%s,"F_IID" lib=%s,"F_IID" clr=%s,%d,%d deleted=%d\n",
                              AS_UID_toString(fr.gkFragment_getReadUID()), fr.gkFragment_getReadIID(),
                              AS_UID_toString(mateuid), mateiid,
                              AS_UID_toString(libuid), libiid,
                              AS_READ_CLEAR_NAMES[dumpClear], clrBeg, clrEnd,
                              fr.gkFragment_getIsDeleted());

        } else {
          AS_UTL_writeFastA(stdout, seq, clrEnd-clrBeg,
                            ">%s,"F_IID" mate=%s,"F_IID" lib=%s,"F_IID" clr=%s,%d,%d deleted=%d\n",
                            AS_UID_toString(fr.gkFragment_getReadUID()), fr.gkFragment_getReadIID(),
                            AS_UID_toString(mateuid), mateiid,
                            AS_UID_toString(libuid), libiid,
                            AS_READ_CLEAR_NAMES[dumpClear], clrBeg, clrEnd,
                            fr.gkFragment_getIsDeleted());
	}
      }
    }
  }

  delete fs;
  delete gkp;
}


void
dumpGateKeeperAsFRG(char       *gkpStoreName,
                    int         dumpFormat,
                    AS_IID      begIID,
                    AS_IID      endIID,
                    char       *iidToDump,
                    int         doNotFixMates,
                    int         dumpAllReads,
                    int         dumpFRGClear) {
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
  //  relative to begIID.

  if (begIID < 1)
    begIID = 1;
  if (gkp->gkStore_getNumFragments() < endIID)
    endIID = gkp->gkStore_getNumFragments();

  if (iidToDump == NULL) {
    iidToDump = (char *)safe_calloc(endIID + 1, sizeof(char));
    for (i=begIID; i<=endIID; i++)
      iidToDump[i] = 1;
  }

  //  Pass 1: Scan the fragments, build a list of libraries to dump,
  //  also add in any mates that were omitted from the input.

  fprintf(stderr, "Scanning store to find libraries used.\n");

  libToDump = (int    *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(int));
  libUID    = (AS_UID *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(AS_UID));
  frgUID    = (AS_UID *)safe_calloc(endIID+1,                         sizeof(AS_UID));

  fs = new gkStream(gkp, begIID, endIID, GKFRAGMENT_INF);

  while (fs->next(&fr)) {
    if ((iidToDump[fr.gkFragment_getReadIID()]) &&
        (dumpAllReads || !fr.gkFragment_getIsDeleted())) {

      frgUID[fr.gkFragment_getReadIID()] = fr.gkFragment_getReadUID();

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
#warning unsafe conversion of orient
        lmesg.link_orient = (OrientType)AS_READ_ORIENT_NAMES[fr.gkFragment_getOrientation()][0];

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
  fs = new gkStream(gkp, begIID, endIID, GKFRAGMENT_QLT);

  while (fs->next(&fr)) {
    FragMesg  fmesg;
    LinkMesg  lmesg;

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
        lmesg.type        = AS_MATE;
#warning unsafe conversion of orient
        lmesg.link_orient = (OrientType)AS_READ_ORIENT_NAMES[fr.gkFragment_getOrientation()][0];
        lmesg.frag1       = frgUID[fr.gkFragment_getMateIID()];
        lmesg.frag2       = fr.gkFragment_getReadUID();
        lmesg.distance    = libUID[fr.gkFragment_getLibraryIID()];

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
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         doNotFixMates,
                        int         dumpFRGClear) {
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
  //  relative to begIID.

  if (begIID < 1)
    begIID = 1;
  if (gkp->gkStore_getNumFragments() < endIID)
    endIID = gkp->gkStore_getNumFragments();

  if (iidToDump == NULL) {
    iidToDump = (char *)safe_calloc(endIID + 1, sizeof(char));
    for (i=begIID; i<=endIID; i++)
      iidToDump[i] = 1;
  }

  //  Pass 1: Scan the fragments, build a list of libraries to dump,
  //  also add in any mates that were omitted from the input.

  fprintf(stderr, "Scanning store to find libraries used.\n");

  libToDump = (int    *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(int));
  libUID    = (AS_UID *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(AS_UID));
  frgUID    = (AS_UID *)safe_calloc(endIID+1,                         sizeof(AS_UID));


  fs = new gkStream(gkp, begIID, endIID, GKFRAGMENT_INF);

  while (fs->next(&fr)) {
    frgUID[fr.gkFragment_getReadIID()] = fr.gkFragment_getReadUID();

    if (iidToDump[fr.gkFragment_getReadIID()]) {
      libToDump[fr.gkFragment_getLibraryIID()]++;

      if ((fr.gkFragment_getMateIID() > 0) && (iidToDump[fr.gkFragment_getMateIID()] == 0)) {
        mateAdded++;
        if (doNotFixMates == 0)
          iidToDump[fr.gkFragment_getMateIID()] = 1;
      }
    }
  }

  delete fs;

  fprintf(stderr, "%sdded %d reads to maintain mate relationships.\n",
          (doNotFixMates) ? "Would have a" : "A", mateAdded);

  fprintf(stderr, "Dumping %d fragments from unknown library (version 1 has these)\n", libToDump[0]);

  for (i=1; i<=gkp->gkStore_getNumLibraries(); i++)
    fprintf(stderr, "Dumping %d fragments from library IID %d\n", libToDump[i], i);

  char  fname[FILENAME_MAX];
  char  qname[FILENAME_MAX];

  sprintf(fname, "%s.fna",      prefix);
  sprintf(qname, "%s.fna.qual", prefix);

  errno = 0;
  FILE *f = fopen(fname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", fname, strerror(errno)), exit(1);

  FILE *q = fopen(qname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", qname, strerror(errno)), exit(1);


  for (i=1; i<=gkp->gkStore_getNumLibraries(); i++) {
    if (libToDump[i]) {
      gkLibrary  *gkpl = gkp->gkStore_getLibrary(i);

      //  We don't really need to cache UIDs anymore;
      //  gkStore_getLibrary should be doing this for us.  It's
      //  already implemented, and slightly more efficient, so BPW
      //  left it in.

      libUID[i] = gkpl->libraryUID;
    }
  }


  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.
  //
  fs = new gkStream(gkp, begIID, endIID, GKFRAGMENT_QLT);

  while (fs->next(&fr)) {
    FragMesg  fmesg;
    LinkMesg  lmesg;

    //  Newbler is not happy at all with clear ranges < 1 base.  We do
    //  not dump those.

    int  lclr = fr.gkFragment_getClearRegionBegin(dumpFRGClear) + 1;
    int  rclr = fr.gkFragment_getClearRegionEnd  (dumpFRGClear);

    if ((iidToDump[fr.gkFragment_getReadIID()]) && (lclr < rclr)) {
      char    defline[1024];

      if (fr.gkFragment_getMateIID()) {
        AS_IID  id1 = fr.gkFragment_getReadIID();
        AS_IID  id2 = fr.gkFragment_getMateIID();

        sprintf(defline, ">%s template=%d+%d dir=%c library=%s trim=%d-%d\n",
                //  ID
                AS_UID_toString(fr.gkFragment_getReadUID()),
                //  template
                (id1 < id2) ? id1 : id2,
                (id1 < id2) ? id2 : id1,
                //  dir
                (id1 < id2) ? 'F' : 'R',
                //  library
                AS_UID_toString(libUID[fr.gkFragment_getLibraryIID()]),
                //  trim
                lclr,
                rclr);
      } else {
        sprintf(defline, ">%s trim=%d-%d\n",
                //  ID
                AS_UID_toString(fr.gkFragment_getReadUID()),
                //  trim
                lclr,
                rclr);
      }

      AS_UTL_writeFastA(f,
                        fr.gkFragment_getSequence(),
                        fr.gkFragment_getSequenceLength(),
                        defline, NULL);

      AS_UTL_writeQVFastA(q,
			  fr.gkFragment_getQuality(),
                          fr.gkFragment_getSequenceLength(),
                          defline, NULL);
    }
  }

  delete fs;
  delete gkp;
}

void
dumpGateKeeperAsVelvet(char       *gkpStoreName,
                        char       *prefix,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         doNotFixMates,
                        int         dumpFRGClear) {
  gkFragment        fr;
  gkFragment        mate;
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
  //  relative to begIID.

  if (begIID < 1)
    begIID = 1;
  if (gkp->gkStore_getNumFragments() < endIID)
    endIID = gkp->gkStore_getNumFragments();

  if (iidToDump == NULL) {
    iidToDump = (char *)safe_calloc(endIID + 1, sizeof(char));
    for (i=begIID; i<=endIID; i++)
      iidToDump[i] = 1;
  }

  //  Pass 1: Scan the fragments, build a list of libraries to dump,
  //  also add in any mates that were omitted from the input.

  fprintf(stderr, "Scanning store to find libraries used.\n");

  libToDump = (int    *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(int));
  libUID    = (AS_UID *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(AS_UID));
  frgUID    = (AS_UID *)safe_calloc(endIID+1,                           sizeof(AS_UID));

  fs = new gkStream(gkp, begIID, endIID, GKFRAGMENT_INF);

  while (fs->next(&fr)) {
    frgUID[fr.gkFragment_getReadIID()] = fr.gkFragment_getReadUID();

    if (iidToDump[fr.gkFragment_getReadIID()]) {
      libToDump[fr.gkFragment_getLibraryIID()]++;

      if ((fr.gkFragment_getMateIID() > 0) && (iidToDump[fr.gkFragment_getMateIID()] == 0)) {
        mateAdded++;
        if (doNotFixMates == 0)
          iidToDump[fr.gkFragment_getMateIID()] = 1;
      }
    }
  }

  fprintf(stderr, "%sdded %d reads to maintain mate relationships.\n",
          (doNotFixMates) ? "Would have a" : "A", mateAdded);

  fprintf(stderr, "Dumping %d fragments from unknown library (version 1 has these)\n", libToDump[0]);

  for (i=1; i<=gkp->gkStore_getNumLibraries(); i++)
    fprintf(stderr, "Dumping %d fragments from library IID %d\n", libToDump[i], i);

  char  fname[FILENAME_MAX];
  char  uname[FILENAME_MAX];

  sprintf(fname, "%s.fastq",              prefix);
  sprintf(uname, "%s.unmated.fastq",      prefix);

  errno = 0;
  FILE *f = fopen(fname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", fname, strerror(errno)), exit(1);

  FILE *u = fopen(uname, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", uname, strerror(errno)), exit(1);

  for (i=1; i<=gkp->gkStore_getNumLibraries(); i++) {
    if (libToDump[i]) {
      gkLibrary  *gkpl = gkp->gkStore_getLibrary(i);

      //  We don't really need to cache UIDs anymore;
      //  gkStore_getLibrary should be doing this for us.  It's
      //  already implemented, and slightly more efficient, so BPW
      //  left it in.

      libUID[i] = gkpl->libraryUID;
    }
  }

  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.
  //
  fs = new gkStream(gkp, begIID, endIID, GKFRAGMENT_QLT);

  while (fs->next(&fr)) {
    FragMesg  fmesg;
    LinkMesg  lmesg;

    //  Newbler is not happy at all with clear ranges < 1 base.  We do
    //  not dump those.

    int  lclr = fr.gkFragment_getClearRegionBegin(dumpFRGClear) + 1;
    int  rclr = fr.gkFragment_getClearRegionEnd  (dumpFRGClear);

    if (!fr.gkFragment_getIsDeleted() && (iidToDump[fr.gkFragment_getReadIID()]) && (lclr < rclr)) {
      char    defline[1024];

      if (fr.gkFragment_getMateIID()) {
        AS_IID  id1 = fr.gkFragment_getReadIID();
        AS_IID  id2 = fr.gkFragment_getMateIID();

        sprintf(defline, "@%s template=%d+%d dir=%c library=%s trim=%d-%d\n",
                //  ID
                AS_UID_toString(fr.gkFragment_getReadUID()),
                //  template
                (id1 < id2) ? id1 : id2,
                (id1 < id2) ? id2 : id1,
                //  dir
                (id1 < id2) ? 'F' : 'R',
                //  library
                AS_UID_toString(libUID[fr.gkFragment_getLibraryIID()]),
                //  trim
                lclr,
                rclr);

         AS_UTL_writeFastAWithBreaks(f,
                           fr.gkFragment_getSequence(),
                           fr.gkFragment_getSequenceLength(),
                           0,
                           defline, NULL);
   
         // for fastq file we need different header character for quality
         defline[0] = '+';
         AS_UTL_writeQVFastQWithBreaks(f,
              fr.gkFragment_getQuality(),
                             fr.gkFragment_getSequenceLength(),
                             0,
                             defline, NULL);

         // now dump the mate
         gkp->gkStore_getFragment(id2, &mate, GKFRAGMENT_QLT);
         lclr = mate.gkFragment_getClearRegionBegin(dumpFRGClear) + 1;
         rclr = mate.gkFragment_getClearRegionEnd  (dumpFRGClear);
         sprintf(defline, "@%s template=%d+%d dir=%c library=%s trim=%d-%d\n",
                //  ID
                 AS_UID_toString(mate.gkFragment_getReadUID()),
                //  template
                (id1 < id2) ? id1 : id2,
                (id1 < id2) ? id2 : id1,
                //  dir
                (id1 < id2) ? 'F' : 'R',
                //  library
                AS_UID_toString(libUID[mate.gkFragment_getLibraryIID()]),
                //  trim
                lclr,
                rclr);
                
         AS_UTL_writeFastAWithBreaks(f,
                                     mate.gkFragment_getSequence(),
                                     mate.gkFragment_getSequenceLength(),
                                     0,
                                     defline, NULL);
   
         // for fastq file we need different header character for quality
         defline[0] = '+';
         AS_UTL_writeQVFastQWithBreaks(f,
                                       mate.gkFragment_getQuality(),
                                       mate.gkFragment_getSequenceLength(),
                                       0,
                                       defline, NULL);
         iidToDump[id2] = 0; // we dumped the mate we don't need it anymore
      } else {
        // unmated need to be in a separate file
        sprintf(defline, "@%s trim=%d-%d\n",
                //  ID
                AS_UID_toString(fr.gkFragment_getReadUID()),
                //  trim
                lclr,
                rclr);
         AS_UTL_writeFastAWithBreaks(u,
                                     fr.gkFragment_getSequence(),
                                     fr.gkFragment_getSequenceLength(),
                                     0,
                                     defline, NULL);
   
         // for fastq file we need different header character for quality
         defline[0] = '+';
         AS_UTL_writeQVFastQWithBreaks(u,
                                       fr.gkFragment_getQuality(),
                                       fr.gkFragment_getSequenceLength(),
                                       0,
                                       defline, NULL);

      }
    }
  }

  delete fs;
  delete gkp;
}
