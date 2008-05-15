
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

static char const *rcsid = "$Id: AS_GKP_dump.c,v 1.35 2008-05-15 21:42:10 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_GKP_include.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_UTL_fasta.h"


void
dumpGateKeeperInfo(char       *gkpStoreName) {

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  fprintf(stdout, "LOAD STATS\n");
  fprintf(stdout, "\n");
  fprintf(stdout, F_U32"\tbatInput\n",    gkp->gkp.batInput);
  fprintf(stdout, F_U32"\tbatLoaded\n",   gkp->gkp.batLoaded);
  fprintf(stdout, F_U32"\tbatErrors\n",   gkp->gkp.batErrors);
  fprintf(stdout, F_U32"\tbatWarnings\n", gkp->gkp.batWarnings);
  fprintf(stdout, "\n");
  fprintf(stdout, F_U32"\tlibInput\n",    gkp->gkp.libInput);
  fprintf(stdout, F_U32"\tlibLoaded\n",   gkp->gkp.libLoaded);
  fprintf(stdout, F_U32"\tlibErrors\n",   gkp->gkp.libErrors);
  fprintf(stdout, F_U32"\tlibWarnings\n", gkp->gkp.libWarnings);
  fprintf(stdout, "\n");
  fprintf(stdout, F_U32"\tfrgInput\n",    gkp->gkp.frgInput);
  fprintf(stdout, F_U32"\tfrgLoaded\n",   gkp->gkp.frgLoaded);
  fprintf(stdout, F_U32"\tfrgErrors\n",   gkp->gkp.frgErrors);
  fprintf(stdout, F_U32"\tfrgWarnings\n", gkp->gkp.frgWarnings);
  fprintf(stdout, "\n");
  fprintf(stdout, F_U32"\tlkgInput\n",    gkp->gkp.lkgInput);
  fprintf(stdout, F_U32"\tlkgLoaded\n",   gkp->gkp.lkgLoaded);
  fprintf(stdout, F_U32"\tlkgErrors\n",   gkp->gkp.lkgErrors);
  fprintf(stdout, F_U32"\tlkgWarnings\n", gkp->gkp.lkgWarnings);
  fprintf(stdout, "\n");
  fprintf(stdout, F_U32"\tsffInput\n",    gkp->gkp.sffInput);
  fprintf(stdout, F_U32"\tsffLoaded\n",   gkp->gkp.sffLoaded);
  fprintf(stdout, F_U32"\tsffErrors\n",   gkp->gkp.sffErrors);
  fprintf(stdout, F_U32"\tsffWarnings\n", gkp->gkp.sffWarnings);
  fprintf(stdout, "\n");
  fprintf(stdout, F_U32"\tsffLibCreated\n", gkp->gkp.sffLibCreated);

  fragRecord    fr;
  FragStream   *fs = openFragStream(gkp, FRAG_S_INF);

  int           i, j;

  uint32  numActiveFrag   = 0;
  uint32  numDeletedFrag  = 0;
  uint32  numMatedFrag    = 0;
  uint64 *clearLength     = (uint64 *)safe_calloc(sizeof(uint64), (AS_READ_CLEAR_UNTRIM + 1));

  uint32  *numActivePerLib  = (uint32  *)safe_calloc(sizeof(uint32),   (getNumGateKeeperLibraries(gkp) + 1));
  uint32  *numDeletedPerLib = (uint32  *)safe_calloc(sizeof(uint32),   (getNumGateKeeperLibraries(gkp) + 1));
  uint32  *numMatedPerLib   = (uint32  *)safe_calloc(sizeof(uint32),   (getNumGateKeeperLibraries(gkp) + 1));
  uint64 **clearLengthByLib = (uint64 **)safe_calloc(sizeof(uint64 *), (getNumGateKeeperLibraries(gkp) + 1));

  for (j=0; j<getNumGateKeeperLibraries(gkp) + 1; j++)
    clearLengthByLib[j] = (uint64 *)safe_calloc(sizeof(uint64), (AS_READ_CLEAR_UNTRIM + 1));

  while (nextFragStream(fs, &fr)) {
    AS_IID     lib = getFragRecordLibraryIID(&fr);

    if (getFragRecordIsDeleted(&fr)) {
      numDeletedFrag++;
      numDeletedPerLib[lib]++;
    } else {
      numActiveFrag++;
      numActivePerLib[lib]++;

      if (getFragRecordMateIID(&fr) > 0) {
        numMatedFrag++;
        numMatedPerLib[lib]++;
      }

      for (i=0; i<AS_READ_CLEAR_UNTRIM + 1; i++) {
        uint32 clrlen = getFragRecordClearRegionEnd(&fr, i) - getFragRecordClearRegionBegin(&fr, i);

        //  What a ucky special case.
        if ((i == AS_READ_CLEAR_QLT) && (fr.gkfr.hasQualityClear == 0))
          continue;
        if ((i == AS_READ_CLEAR_VEC) && (fr.gkfr.hasVectorClear == 0))
          continue;

        clearLength[i]           += clrlen;
        clearLengthByLib[lib][i] += clrlen;
      }
    }
  }

  fprintf(stdout, "\n");
  fprintf(stdout, "GLOBAL STATS\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "numberFRG="F_S32"\n", (uint32)getNumGateKeeperFragments(gkp));
  fprintf(stdout, "numberLIB="F_S32"\n", (uint32)getNumGateKeeperLibraries(gkp));
  fprintf(stdout, "\n");
  fprintf(stdout, "activeFRG="F_U32" deletedFRG="F_U32" matedFRG="F_U32" ", numActiveFrag, numDeletedFrag, numMatedFrag);
  for (i=0; i<AS_READ_CLEAR_UNTRIM + 1; i++)
    fprintf(stdout, "%s="F_U64" ", AS_READ_CLEAR_NAMES[i], clearLength[i]);
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "PER LIBRARY STATS\n");
  fprintf(stdout, "\n");
  for (j=0; j<getNumGateKeeperLibraries(gkp) + 1; j++) {
    fprintf(stdout, "libraryIID="F_IID" activeFRG="F_U32" deletedFRG="F_U32" matedFRG="F_U32" ", j, numActivePerLib[j], numDeletedPerLib[j], numMatedPerLib[j]);
    for (i=0; i<AS_READ_CLEAR_UNTRIM + 1; i++)
      fprintf(stdout, "%s="F_U64" ", AS_READ_CLEAR_NAMES[i], clearLengthByLib[j][i]);
    fprintf(stdout, "\n");
  }

  //  yeah, a memory leak.

  closeFragStream(fs);
  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperBatches(char       *gkpStoreName,
                      AS_IID      begIID,
                      AS_IID      endIID,
                      char       *iidToDump, 
                      int         asTable) {
  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  int       i;

  if (begIID < getFirstElemStore(gkp->bat))
    begIID = getFirstElemStore(gkp->bat);
  if (getLastElemStore(gkp->bat) < endIID)
    endIID = getLastElemStore(gkp->bat);

  if (asTable)
    fprintf(stdout, "UID\tIID\tName\tnumFrags\tnumLibs\tnumShadowLibs\n");

  for (i=begIID; i<=endIID; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      GateKeeperBatchRecord gkpb;

      getGateKeeperBatch(gkp, i, &gkpb);

      if (asTable) {
        fprintf(stdout, "%s\t"F_IID"\t%s\n",
                AS_UID_toString(gkpb.batchUID), i,
                (gkpb.name[0]) ? gkpb.name : ".");
      } else {
        fprintf(stdout, "batchIdent   = %s,"F_IID"\n", AS_UID_toString(gkpb.batchUID), i);
        fprintf(stdout, "batchName    = %s\n", gkpb.name);
        chomp(gkpb.comment);
        fprintf(stdout, "batchComment\n");
        if (gkpb.comment[0] != 0)
          fprintf(stdout, "%s\n", gkpb.comment);
        fprintf(stdout, "batchCommentEnd\n");
      }
    }
  }

  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperLibraries(char       *gkpStoreName,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump, 
                        int         asTable) {
  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  int       i;

  if (begIID < getFirstElemStore(gkp->lib))
    begIID = getFirstElemStore(gkp->lib);
  if (getLastElemStore(gkp->lib) < endIID)
    endIID = getLastElemStore(gkp->lib);

  if (asTable)
    fprintf(stdout, "UID\tIID\tOrientation\tMean\tStdDev\tNumFeatures\n");

  for (i=begIID; i<=endIID; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      GateKeeperLibraryRecord *gkpl = getGateKeeperLibrary(gkp, i);
      LibraryMesg              lmesg;
      int                      nf;

      AS_PER_encodeLibraryFeatures(gkpl, &lmesg);

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

        chomp(gkpl->comment);
        fprintf(stdout, "libraryComment\n");
        if (gkpl->comment[0] != 0)
          fprintf(stdout, "%s\n", gkpl->comment);
        fprintf(stdout, "libraryCommentEnd\n");
      }

      AS_PER_encodeLibraryFeaturesCleanup(&lmesg);
    }
  }       

  closeGateKeeperStore(gkp);
}



void
dumpGateKeeperFragments(char       *gkpStoreName,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump, 
                        int         dumpWithSequence,
                        int         dumpClear,
                        int         asTable) {
  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  fragRecord    fr = {0};
  FragStream   *fs = openFragStream(gkp, (!dumpWithSequence || asTable) ? FRAG_S_INF : FRAG_S_ALL);

  int           i;

  if (begIID < getFirstElemStore(gkp->frg))
    begIID = getFirstElemStore(gkp->frg);
  if (getLastElemStore(gkp->frg) < endIID)
    endIID = getLastElemStore(gkp->frg);

  resetFragStream(fs, begIID, endIID);

  if (asTable)
    fprintf(stdout, "UID\tIID\tmateUID\tmateIID\tlibUID\tlibIID\tisDeleted\tisNonRandom\tStatus\tOrient\tLength\tclrBegin%s\tclrEnd%s\n",
            AS_READ_CLEAR_NAMES[dumpClear], AS_READ_CLEAR_NAMES[dumpClear]);

  while (nextFragStream(fs, &fr)) {
    if ((iidToDump == NULL) || (iidToDump[getFragRecordIID(&fr)])) {
      AS_IID     mateiid = getFragRecordMateIID(&fr);
      AS_UID     mateuid = {0};

      if (mateiid > 0)
        mateuid = getGatekeeperIIDtoUID(gkp, mateiid, AS_IID_FRG);

      AS_IID     libiid = getFragRecordLibraryIID(&fr);
      AS_UID     libuid = {0};

      if (libiid > 0)
        libuid = getGateKeeperLibrary(gkp, libiid)->libraryUID;

      if (asTable) {
        fprintf(stdout, "%s\t"F_IID"\t%s\t"F_IID"\t%s\t"F_IID"\t%d\t%d\t%s\t%s\t%d\t%d\t%d\n",
                AS_UID_toString1(getFragRecordUID(&fr)), getFragRecordIID(&fr),
                AS_UID_toString2(mateuid), mateiid,
                AS_UID_toString3(libuid), libiid,
                getFragRecordIsDeleted(&fr),
                getFragRecordIsNonRandom(&fr),
                AS_READ_STATUS_NAMES[fr.gkfr.status],
                AS_READ_ORIENT_NAMES[fr.gkfr.orientation],
                getFragRecordSequenceLength(&fr),
                getFragRecordClearRegionBegin(&fr, dumpClear),
                getFragRecordClearRegionEnd  (&fr, dumpClear));
      } else {
        fprintf(stdout, "fragmentIdent           = %s,"F_IID"\n", AS_UID_toString(getFragRecordUID(&fr)), getFragRecordIID(&fr));
        fprintf(stdout, "fragmentMate            = %s,"F_IID"\n", AS_UID_toString(mateuid), mateiid);
        fprintf(stdout, "fragmentLibrary         = %s,"F_IID"\n", AS_UID_toString(libuid), libiid);

        fprintf(stdout, "fragmentIsDeleted       = %d\n", getFragRecordIsDeleted(&fr));
        fprintf(stdout, "fragmentIsNonRandom     = %d\n", getFragRecordIsNonRandom(&fr));
        fprintf(stdout, "fragmentStatus          = %s\n", AS_READ_STATUS_NAMES[fr.gkfr.status]);
        fprintf(stdout, "fragmentOrientation     = %s\n", AS_READ_ORIENT_NAMES[fr.gkfr.orientation]);

        fprintf(stdout, "fragmentHasVectorClear  = %d\n", fr.gkfr.hasVectorClear);
        fprintf(stdout, "fragmentHasQualityClear = %d\n", fr.gkfr.hasQualityClear);

        fprintf(stdout, "fragmentPlate           = %s\n", AS_UID_toString(fr.gkfr.plateUID));
        fprintf(stdout, "fragmentPlateLocation   = %d\n", fr.gkfr.plateLocation);

        fprintf(stdout, "fragmentSeqLen          = %d\n", getFragRecordSequenceLength(&fr));
        fprintf(stdout, "fragmentHPSLen          = %d\n", getFragRecordHPSLength(&fr));
        fprintf(stdout, "fragmentSrcLen          = %d\n", getFragRecordSourceLength(&fr));

        if (dumpWithSequence) {
          unsigned int   clrBeg = getFragRecordClearRegionBegin(&fr, dumpClear);
          unsigned int   clrEnd = getFragRecordClearRegionEnd  (&fr, dumpClear);
          char          *seq = getFragRecordSequence(&fr);
          char          *qlt = getFragRecordQuality(&fr);
          char          *hps = getFragRecordHPS(&fr);
          char          *src = getFragRecordSource(&fr);

          chomp(src);

          seq[clrEnd] = 0;
          qlt[clrEnd] = 0;

          fprintf(stdout, "fragmentSequence%-6s  = %s\n", AS_READ_CLEAR_NAMES[dumpClear], seq + clrBeg);
          fprintf(stdout, "fragmentQuality%-6s   = %s\n", AS_READ_CLEAR_NAMES[dumpClear], qlt + clrBeg);
          fprintf(stdout, "fragmentHPS             = %s\n", hps);
          fprintf(stdout, "fragmentSource          = %s\n", src);
        }

        for (i=0; i <= AS_READ_CLEAR_LATEST; i++) {
          fprintf(stdout, "fragmentClear%-6s     = %d,%d\n", 
                  AS_READ_CLEAR_NAMES[i],
                  fr.gkfr.clearBeg[i], fr.gkfr.clearEnd[i]);
        }

        fprintf(stdout, "fragmentSeqOffset       = "F_U64"\n", fr.gkfr.seqOffset);
        fprintf(stdout, "fragmentQltOffset       = "F_U64"\n", fr.gkfr.qltOffset);
        fprintf(stdout, "fragmentHpsOffset       = "F_U64"\n", fr.gkfr.hpsOffset);
        fprintf(stdout, "fragmentSrcOffset       = "F_U64"\n", fr.gkfr.srcOffset);
      }
    }
  }

  closeFragStream(fs);
  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperAsFasta(char       *gkpStoreName,
                      AS_IID      begIID,
                      AS_IID      endIID,
                      char       *iidToDump, 
                      int         dumpAllReads,
                      int         dumpClear,
                      int         dumpQuality) {
  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  fragRecord    fr;
  FragStream   *fs = openFragStream(gkp, (dumpQuality) ? FRAG_S_QLT : FRAG_S_SEQ);

  int           i;

  if (begIID < getFirstElemStore(gkp->frg))
    begIID = getFirstElemStore(gkp->frg);
  if (getLastElemStore(gkp->frg) < endIID)
    endIID = getLastElemStore(gkp->frg);

  resetFragStream(fs, begIID, endIID);

  while (nextFragStream(fs, &fr)) {
    if ((iidToDump == NULL) || (iidToDump[getFragRecordIID(&fr)])) {
      if (dumpAllReads || !getFragRecordIsDeleted(&fr)) {
        AS_IID     mateiid = getFragRecordMateIID(&fr);
        AS_UID     mateuid = {0};

        if (mateiid > 0)
          mateuid = getGatekeeperIIDtoUID(gkp, mateiid, AS_IID_FRG);

        AS_IID     libiid = getFragRecordLibraryIID(&fr);
        AS_UID     libuid = {0};

        if (libiid > 0)
          libuid = getGateKeeperLibrary(gkp, libiid)->libraryUID;

        unsigned int   clrBeg   = getFragRecordClearRegionBegin(&fr, dumpClear);
        unsigned int   clrEnd   = getFragRecordClearRegionEnd  (&fr, dumpClear);
        char          *seqStart = getFragRecordSequence(&fr);

        if (dumpQuality) 
          seqStart = getFragRecordQuality(&fr);

        char          *seq      = seqStart+clrBeg;
        seq[clrEnd] = 0;

        if (dumpQuality >=2) {
           fprintf(stdout, ">%s,"F_IID" mate=%s,"F_IID" lib=%s,"F_IID" clr=%s,%d,%d deleted=%d\n",
               	AS_UID_toString1(getFragRecordUID(&fr)), getFragRecordIID(&fr),
                AS_UID_toString2(mateuid), mateiid,
                AS_UID_toString3(libuid), libiid,
                AS_READ_CLEAR_NAMES[dumpClear], clrBeg, clrEnd,
                getFragRecordIsDeleted(&fr));

           int i = 0;
           int chars = 0;
           for (i = clrBeg; i < clrEnd; i++) {
              fprintf(stdout, "%#2d ", (((int)seq[i])-'0'));
              chars += 3;
              if (chars % 60 == 0) { fprintf(stdout, "\n");}
           }
           fprintf(stdout, "\n");
        } else {
	   AS_UTL_writeFastA(stdout, seq, clrEnd-clrBeg,
		">%s,"F_IID" mate=%s,"F_IID" lib=%s,"F_IID" clr=%s,%d,%d deleted=%d\n",
                	AS_UID_toString1(getFragRecordUID(&fr)), getFragRecordIID(&fr),
                	AS_UID_toString2(mateuid), mateiid,
                	AS_UID_toString3(libuid), libiid,
                	AS_READ_CLEAR_NAMES[dumpClear], clrBeg, clrEnd,
                	getFragRecordIsDeleted(&fr));
	}
      }
    }
  }

  closeFragStream(fs);
  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperAsFRG(char       *gkpStoreName,
                    int         dumpFormat,
                    AS_IID      begIID,
                    AS_IID      endIID,
                    char       *iidToDump,
                    int         doNotFixMates,
                    int         dumpFRGClear) {
  fragRecord        fr;
  FragStream       *fs = NULL;

  unsigned int      firstElem = 0;
  unsigned int      lastElem = 0;

  GenericMesg       pmesg;

  LibraryMesg       libMesg;
  DistanceMesg      dstMesg;
  FragMesg          frgMesg;
  LinkMesg          lnkMesg;

  int               i;

  int              *libToDump;
  AS_UID           *libUID;
  AS_UID           *frgUID;
  int               mateAdded = 0;

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  //  Someone could adjust the iidToDump and frgUID arrays to be
  //  relative to begIID.

  if (begIID < getFirstElemStore(gkp->frg))
    begIID = getFirstElemStore(gkp->frg);
  if (getLastElemStore(gkp->frg) < endIID)
    endIID = getLastElemStore(gkp->frg);

  if (iidToDump == NULL) {
    iidToDump = (char *)safe_calloc(endIID + 1, sizeof(char));
    for (i=begIID; i<=endIID; i++)
      iidToDump[i] = 1;
  }

  //  Pass 1: Scan the fragments, build a list of libraries to dump,
  //  also add in any mates that were omitted from the input.

  fprintf(stderr, "Scanning store to find libraries used.\n");

  libToDump = (int    *)safe_calloc(getLastElemStore(gkp->lib)+1, sizeof(int));
  libUID    = (AS_UID *)safe_calloc(getLastElemStore(gkp->lib)+1, sizeof(AS_UID));
  frgUID    = (AS_UID *)safe_calloc(endIID+1,                     sizeof(AS_UID));

  fs = openFragStream(gkp, FRAG_S_INF);
  resetFragStream(fs, begIID, endIID);

  while (nextFragStream(fs, &fr)) {
    frgUID[getFragRecordIID(&fr)] = getFragRecordUID(&fr);

    if (iidToDump[getFragRecordIID(&fr)]) {
      libToDump[getFragRecordLibraryIID(&fr)]++;

      if ((getFragRecordMateIID(&fr) > 0) && (iidToDump[getFragRecordMateIID(&fr)] == 0)) {
        mateAdded++;
        if (doNotFixMates == 0)
          iidToDump[getFragRecordMateIID(&fr)] = 1;
      }
    }
  }

  fprintf(stderr, "%sdded %d reads to maintain mate relationships.\n",
          (doNotFixMates) ? "Would have a" : "A", mateAdded);

  fprintf(stderr, "Dumping %d fragments from unknown library (version 1 has these)\n", libToDump[0]);

  for (i=1; i<=getLastElemStore(gkp->lib); i++)
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
  for (i=1; i<=getLastElemStore(gkp->lib); i++) {
    if (libToDump[i]) {
      DistanceMesg              dmesg;
      LibraryMesg               lmesg;
      GateKeeperLibraryRecord  *gkpl = getGateKeeperLibrary(gkp, i);

      //  We don't really need to cache UIDs anymore;
      //  getGateKeeperLibrary should be doing this for us.  It's
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
        lmesg.source       = gkpl->comment;
        lmesg.link_orient  = AS_READ_ORIENT_NAMES[gkpl->orientation][0];

        AS_PER_encodeLibraryFeatures(gkpl, &lmesg);
      }

      WriteProtoMesg_AS(stdout, &pmesg);

      if (dumpFormat != 1)
        AS_PER_encodeLibraryFeaturesCleanup(&lmesg);
    }
  }
  closeFragStream(fs);


  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.
  //
  fs = openFragStream(gkp, FRAG_S_ALL);
  resetFragStream(fs, begIID, endIID);

  while (nextFragStream(fs, &fr)) {
    FragMesg  fmesg;
    LinkMesg  lmesg;

    if (iidToDump[getFragRecordIID(&fr)]) {
      pmesg.m = &fmesg;
      pmesg.t = MESG_FRG;

      //  This code used in AS_GKP_dump.c (dumpFRG), and in AS_FGB_io.c
      fmesg.action          = getFragRecordIsDeleted(&fr) ? AS_DELETE : AS_ADD;
      fmesg.eaccession      = getFragRecordUID(&fr);
      fmesg.library_uid     = libUID[getFragRecordLibraryIID(&fr)];
      fmesg.library_iid     = getFragRecordLibraryIID(&fr);
      fmesg.plate_uid       = fr.gkfr.plateUID;
      fmesg.plate_location  = fr.gkfr.plateLocation;
      fmesg.type            = AS_READ;
      fmesg.is_random       = (getFragRecordIsNonRandom(&fr)) ? 0 : 1;
      fmesg.status_code     = AS_READ_STATUS_NAMES[fr.gkfr.status][0];
      fmesg.clear_rng.bgn   = getFragRecordClearRegionBegin(&fr, dumpFRGClear);
      fmesg.clear_rng.end   = getFragRecordClearRegionEnd  (&fr, dumpFRGClear);
      fmesg.clear_vec.bgn   = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_VEC);
      fmesg.clear_vec.end   = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_VEC);
      fmesg.clear_qlt.bgn   = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_QLT);
      fmesg.clear_qlt.end   = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_QLT);
      fmesg.source          = getFragRecordSource(&fr);
      fmesg.sequence        = getFragRecordSequence(&fr);
      fmesg.quality         = getFragRecordQuality(&fr);
      fmesg.hps             = getFragRecordHPS(&fr);
      fmesg.iaccession      = firstElem;

      if (fr.gkfr.hasVectorClear == 0) {
        fmesg.clear_vec.bgn   = 1;
        fmesg.clear_vec.end   = 0;
      }

      if (fr.gkfr.hasQualityClear == 0) {
        fmesg.clear_qlt.bgn   = 1;
        fmesg.clear_qlt.end   = 0;
      }

      WriteProtoMesg_AS(stdout, &pmesg);

      if ((getFragRecordMateIID(&fr) > 0) &&
          (getFragRecordMateIID(&fr) < getFragRecordIID(&fr))) {
        pmesg.m = &lmesg;
        pmesg.t = MESG_LKG;

        lmesg.action      = AS_ADD;
        lmesg.type        = AS_MATE;
        lmesg.link_orient = AS_READ_ORIENT_NAMES[fr.gkfr.orientation][0];
        lmesg.frag1       = frgUID[getFragRecordMateIID(&fr)];
        lmesg.frag2       = getFragRecordUID(&fr);
        lmesg.distance    = libUID[getFragRecordLibraryIID(&fr)];

        WriteProtoMesg_AS(stdout, &pmesg);
      }
    }
  }

  closeFragStream(fs);
  closeGateKeeperStore(gkp);
}




void
dumpGateKeeperAsNewbler(char       *gkpStoreName,
                        char       *prefix,
                        AS_IID      begIID,
                        AS_IID      endIID,
                        char       *iidToDump,
                        int         doNotFixMates,
                        int         dumpFRGClear) {
  fragRecord        fr;
  FragStream       *fs = NULL;

  unsigned int      firstElem = 0;
  unsigned int      lastElem = 0;

  GenericMesg       pmesg;

  LibraryMesg       libMesg;
  DistanceMesg      dstMesg;
  FragMesg          frgMesg;
  LinkMesg          lnkMesg;

  int               i;

  int              *libToDump;
  AS_UID           *libUID;
  AS_UID           *frgUID;
  int               mateAdded = 0;

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  //  Someone could adjust the iidToDump and frgUID arrays to be
  //  relative to begIID.

  if (begIID < getFirstElemStore(gkp->frg))
    begIID = getFirstElemStore(gkp->frg);
  if (getLastElemStore(gkp->frg) < endIID)
    endIID = getLastElemStore(gkp->frg);

  if (iidToDump == NULL) {
    iidToDump = (char *)safe_calloc(endIID + 1, sizeof(char));
    for (i=begIID; i<=endIID; i++)
      iidToDump[i] = 1;
  }

  //  Pass 1: Scan the fragments, build a list of libraries to dump,
  //  also add in any mates that were omitted from the input.

  fprintf(stderr, "Scanning store to find libraries used.\n");

  libToDump = (int    *)safe_calloc(getLastElemStore(gkp->lib)+1, sizeof(int));
  libUID    = (AS_UID *)safe_calloc(getLastElemStore(gkp->lib)+1, sizeof(AS_UID));
  frgUID    = (AS_UID *)safe_calloc(endIID+1,                     sizeof(AS_UID));


  fs = openFragStream(gkp, FRAG_S_INF);
  resetFragStream(fs, begIID, endIID);

  while (nextFragStream(fs, &fr)) {
    frgUID[getFragRecordIID(&fr)] = getFragRecordUID(&fr);

    if (iidToDump[getFragRecordIID(&fr)]) {
      libToDump[getFragRecordLibraryIID(&fr)]++;

      if ((getFragRecordMateIID(&fr) > 0) && (iidToDump[getFragRecordMateIID(&fr)] == 0)) {
        mateAdded++;
        if (doNotFixMates == 0)
          iidToDump[getFragRecordMateIID(&fr)] = 1;
      }
    }
  }

  fprintf(stderr, "%sdded %d reads to maintain mate relationships.\n",
          (doNotFixMates) ? "Would have a" : "A", mateAdded);

  fprintf(stderr, "Dumping %d fragments from unknown library (version 1 has these)\n", libToDump[0]);

  for (i=1; i<=getLastElemStore(gkp->lib); i++)
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


  for (i=1; i<=getLastElemStore(gkp->lib); i++) {
    if (libToDump[i]) {
      GateKeeperLibraryRecord  *gkpl = getGateKeeperLibrary(gkp, i);

      //  We don't really need to cache UIDs anymore;
      //  getGateKeeperLibrary should be doing this for us.  It's
      //  already implemented, and slightly more efficient, so BPW
      //  left it in.

      libUID[i] = gkpl->libraryUID;
    }
  }

  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.
  //
  fs = openFragStream(gkp, FRAG_S_ALL);
  resetFragStream(fs, begIID, endIID);

  while (nextFragStream(fs, &fr)) {
    FragMesg  fmesg;
    LinkMesg  lmesg;

    if (iidToDump[getFragRecordIID(&fr)]) {
      char    defline[1024];

      if (getFragRecordMateIID(&fr)) {
        AS_IID  id1 = getFragRecordIID(&fr);
        AS_IID  id2 = getFragRecordMateIID(&fr);

        sprintf(defline, ">%s template=%d+%d dir=%c library=%s trim=%d-%d\n",
                //  ID
                AS_UID_toString1(getFragRecordUID(&fr)),
                //  template
                (id1 < id2) ? id1 : id2,
                (id1 < id2) ? id2 : id1,
                //  dir
                (id1 < id2) ? 'F' : 'R',
                //  library
                AS_UID_toString2(libUID[getFragRecordLibraryIID(&fr)]),
                //  trim
                getFragRecordClearRegionBegin(&fr, dumpFRGClear) + 1,
                getFragRecordClearRegionEnd  (&fr, dumpFRGClear));
      } else {
        sprintf(defline, ">%s trim=%d-%d\n",
                //  ID
                AS_UID_toString1(getFragRecordUID(&fr)),
                //  trim
                getFragRecordClearRegionBegin(&fr, dumpFRGClear) + 1,
                getFragRecordClearRegionEnd  (&fr, dumpFRGClear));
      }

      AS_UTL_writeFastA(f,
                        getFragRecordSequence(&fr),
                        getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_UNTRIM),
                        defline, NULL);

      AS_UTL_writeQVFastA(q,
                          getFragRecordSequence(&fr),
                          getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_UNTRIM),
                          defline, NULL);
    }
  }

  closeFragStream(fs);
  closeGateKeeperStore(gkp);
}
