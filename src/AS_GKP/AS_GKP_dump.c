
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

/* $Id: AS_GKP_dump.c,v 1.12 2007-04-20 08:57:35 brianwalenz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_GKP_include.h"


//  perl's chomp is pretty nice
#define chomp(S) { char *t=S; while (*t) t++; t--; while (isspace(*t)) *t--=0; }

static const char *ctimec(uint64 createdtime) {
  static char  timestring[256];
  time_t       createdtimetime = createdtime;
  strcpy(timestring, ctime(&createdtimetime));
  chomp(timestring);
  return(timestring);
}


void
dumpGateKeeperInfo(char       *gkpStoreName) {

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  fprintf(stdout, "num fragments        = "F_S32"\n", getNumGateKeeperFragments(gkp->frg));
  fprintf(stdout, "num libraries        = "F_S32"\n", getNumGateKeeperLibrarys(gkp->lib));
  fprintf(stdout, "num shadow libraries = "F_S32"\n", getNumGateKeeperLibrarys(gkp->lis));

  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperBatches(char       *gkpStoreName,
                      CDS_IID_t   begIID,
                      CDS_IID_t   endIID,
                      char       *iidToDump, 
                      int         asTable) {
  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  StoreStat stat;
  int       i;

  statsStore(gkp->bat, &stat);

  if (begIID < stat.firstElem)
    begIID = stat.firstElem;
  if (stat.lastElem < endIID)
    endIID = stat.lastElem;

  if (asTable)
    fprintf(stdout, "UID\tIID\tName\tCreated\tCreated\tnumFrags\tnumLibs\tnumShadowLibs\n");

  for (i=begIID; i<=endIID; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      GateKeeperBatchRecord gkpb;

      getGateKeeperBatchStore(gkp->bat, i, &gkpb);

      if (asTable) {
        fprintf(stdout, F_UID"\t"F_IID"\t%s\t"F_U64"\t%s\t"F_S32"\t"F_S32"\t"F_S32"\n",
                gkpb.batchUID, i,
                (gkpb.name[0]) ? gkpb.name : ".",
                gkpb.created, ctimec(gkpb.created),
                gkpb.numFragments,
                gkpb.numLibraries,
                gkpb.numLibraries_s);
      } else {
        fprintf(stdout, "batchIdent   = "F_UID","F_IID"\n", gkpb.batchUID, i);
        fprintf(stdout, "batchName    = %s\n", gkpb.name);
        fprintf(stdout, "batchCreated = "F_U64" (%s)\n", gkpb.created, ctimec(gkpb.created));
        fprintf(stdout, "batchNFrags  = "F_S32"\n", gkpb.numFragments);
        fprintf(stdout, "batchNLibs   = "F_S32"\n", gkpb.numLibraries);
        fprintf(stdout, "batchNLibsS  = "F_S32"\n", gkpb.numLibraries_s);
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
                        CDS_IID_t   begIID,
                        CDS_IID_t   endIID,
                        char       *iidToDump, 
                        int         asTable) {
  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  StoreStat stat;
  int       i;

  statsStore(gkp->lib, &stat);

  if (begIID < stat.firstElem)
    begIID = stat.firstElem;
  if (stat.lastElem < endIID)
    endIID = stat.lastElem;

  if (asTable)
    fprintf(stdout, "UID\tIID\tName\tCreated\tCreated\tisDeleted\tisRedefined\tOrientation\tMean\tStdDev\n");

  for (i=begIID; i<=endIID; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      GateKeeperLibraryRecord gkpl;

      getGateKeeperLibraryStore(gkp->lib, i, &gkpl);

      if (asTable) {
        fprintf(stdout, F_UID"\t"F_IID"\t%s\t"F_U64"\t%s\t%d\t%d\t%s\t%f\t%f\n",
                gkpl.libraryUID, i,
                (gkpl.name[0]) ? gkpl.name : ".",
                gkpl.created, ctimec(gkpl.created),
                gkpl.deleted,
                gkpl.redefined,
                AS_READ_ORIENT_NAMES[gkpl.orientation],
                gkpl.mean,
                gkpl.stddev);
      } else {
        fprintf(stdout, "libraryIdent         = "F_UID","F_IID"\n", gkpl.libraryUID, i);
        fprintf(stdout, "libraryName          = %s\n", gkpl.name);
        fprintf(stdout, "libraryCreated       = "F_U64" (%s)\n", gkpl.created, ctimec(gkpl.created));
        fprintf(stdout, "libraryDeleted       = %d\n", gkpl.deleted);
        fprintf(stdout, "libraryRedefined     = %d\n", gkpl.redefined);
        fprintf(stdout, "libraryOrientation   = %s\n", AS_READ_ORIENT_NAMES[gkpl.orientation]);
        fprintf(stdout, "libraryMean          = %f\n", gkpl.mean);
        fprintf(stdout, "libraryStdDev        = %f\n", gkpl.stddev);
        fprintf(stdout, "libraryNumFeatures   = %d\n", gkpl.numFeatures);
        fprintf(stdout, "libraryPrevInstance  = "F_IID"\n", gkpl.prevInstanceID);
        fprintf(stdout, "libraryPrevID        = "F_IID"\n", gkpl.prevID);
        fprintf(stdout, "libraryBirthBatch    = %d\n", gkpl.birthBatch);
        fprintf(stdout, "libraryDeathBatch    = %d\n", gkpl.deathBatch);
        chomp(gkpl.comment);
        fprintf(stdout, "libraryComment\n");
        if (gkpl.comment[0] != 0)
          fprintf(stdout, "%s\n", gkpl.comment);
        fprintf(stdout, "libraryCommentEnd\n");
      }
    }
  }       

  statsStore(gkp->lis, &stat);

  for (i=stat.firstElem; i<=stat.lastElem; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      GateKeeperLibraryRecord gkpl;

      getGateKeeperLibraryStore(gkp->lis, i, &gkpl);

      //  This is a verbatim copy of the above block, except we change
      //  the label to indiacte it is a shadow library.
      //
      if (asTable) {
        fprintf(stdout, F_UID"\t"F_IID"\t%s\t"F_U64"\t%s\t%d\t%d\t%s\t%f\t%f\tshadow\n",
                gkpl.libraryUID, i,
                (gkpl.name[0]) ? gkpl.name : ".",
                gkpl.created, ctimec(gkpl.created),
                gkpl.deleted,
                gkpl.redefined,
                AS_READ_ORIENT_NAMES[gkpl.orientation],
                gkpl.mean,
                gkpl.stddev);
      } else {
        fprintf(stdout, "shadowLibraryIdent         = "F_UID","F_IID"\n", gkpl.libraryUID, i);
        fprintf(stdout, "shadowLibraryName          = %s\n", gkpl.name);
        fprintf(stdout, "shadowLibraryCreated       = "F_U64" (%s)\n", gkpl.created, ctimec(gkpl.created));
        fprintf(stdout, "shadowLibraryDeleted       = %d\n", gkpl.deleted);
        fprintf(stdout, "shadowLibraryRedefined     = %d\n", gkpl.redefined);
        fprintf(stdout, "shadowLibraryOrientation   = %s\n", AS_READ_ORIENT_NAMES[gkpl.orientation]);
        fprintf(stdout, "shadowLibraryMean          = %f\n", gkpl.mean);
        fprintf(stdout, "shadowLibraryStdDev        = %f\n", gkpl.stddev);
        fprintf(stdout, "shadowLibraryNumFeatures   = %d\n", gkpl.numFeatures);
        fprintf(stdout, "shadowLibraryPrevInstance  = "F_IID"\n", gkpl.prevInstanceID);
        fprintf(stdout, "shadowLibraryPrevID        = "F_IID"\n", gkpl.prevID);
        fprintf(stdout, "shadowLibraryBirthBatch    = %d\n", gkpl.birthBatch);
        fprintf(stdout, "shadowLibraryDeathBatch    = %d\n", gkpl.deathBatch);
        chomp(gkpl.comment);
        fprintf(stdout, "shadowLibraryComment\n");
        if (gkpl.comment[0] != 0)
          fprintf(stdout, "%s\n", gkpl.comment);
        fprintf(stdout, "shadowLibraryCommentEnd\n");
      }
    }
  }

  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperFragments(char       *gkpStoreName,
                        CDS_IID_t   begIID,
                        CDS_IID_t   endIID,
                        char       *iidToDump, 
                        int         dumpWithSequence,
                        int         dumpClear,
                        int         asTable) {
  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  fragRecord   *fr = new_fragRecord();
  FragStream   *fs = openFragStream(gkp, (!dumpWithSequence || asTable) ? FRAG_S_INF : FRAG_S_ALL);
  StoreStat     stat;

  int           i;

  statsStore(gkp->frg, &stat);

  if (begIID < stat.firstElem)
    begIID = stat.firstElem;
  if (stat.lastElem < endIID)
    endIID = stat.lastElem;

  resetFragStream(fs, begIID, endIID);

  if (asTable)
    fprintf(stdout, "UID\tIID\tmateUID\tmateIID\tlibUID\tlibIID\tisDeleted\tisNonRandom\tStatus\tOrient\tLength\tclrBegin\tclrEnd\n");

  while (nextFragStream(fs, fr)) {
    if ((iidToDump == NULL) || (iidToDump[getFragRecordIID(fr)])) {
      if (asTable) {
        fprintf(stdout, F_UID"\t"F_IID"\t"F_UID"\t"F_IID"\t"F_UID"\t"F_IID"\t%d\t%d\t%s\t%s\t%d\t%d\t%d\n",
                getFragRecordUID(fr), getFragRecordIID(fr),
                (CDS_UID_t)0, getFragRecordMateIID(fr),
                (CDS_UID_t)0, getFragRecordLibraryIID(fr),
                getFragRecordIsDeleted(fr),
                getFragRecordIsNonRandom(fr),
                AS_READ_STATUS_NAMES[fr->gkfr.status],
                AS_READ_ORIENT_NAMES[fr->gkfr.orientation],
                getFragRecordSequenceLength(fr),
                getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_LATEST),
                getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_LATEST));
      } else {
        fprintf(stdout, "fragmentIdent           = "F_UID","F_IID"\n", getFragRecordUID(fr), getFragRecordIID(fr));
        fprintf(stdout, "fragmentMate            = "F_UID","F_IID"\n", (CDS_UID_t)0, getFragRecordMateIID(fr));
        fprintf(stdout, "fragmentLibrary         = "F_UID","F_IID"\n", (CDS_UID_t)0, getFragRecordLibraryIID(fr));

        fprintf(stdout, "fragmentIsDeleted       = %d\n", getFragRecordIsDeleted(fr));
        fprintf(stdout, "fragmentIsNonRandom     = %d\n", getFragRecordIsNonRandom(fr));
        fprintf(stdout, "fragmentStatus          = %s\n", AS_READ_STATUS_NAMES[fr->gkfr.status]);
        fprintf(stdout, "fragmentOrientation     = %s\n", AS_READ_ORIENT_NAMES[fr->gkfr.orientation]);
        fprintf(stdout, "fragmentPlate           = "F_UID"\n", fr->gkfr.plateUID);
        fprintf(stdout, "fragmentPlateLocation   = %d\n", fr->gkfr.plateLocation);

        fprintf(stdout, "fragmentSeqLen          = %d\n", getFragRecordSequenceLength(fr));
        fprintf(stdout, "fragmentHPSLen          = %d\n", getFragRecordHPSLength(fr));
        fprintf(stdout, "fragmentSrcLen          = %d\n", getFragRecordSourceLength(fr));

        fprintf(stdout, "fragmentBirthBatch      = %d\n", fr->gkfr.birthBatch);
        fprintf(stdout, "fragmentDeathBatch      = %d\n", fr->gkfr.birthBatch);

        if (dumpWithSequence) {
          unsigned int   clrBeg = getFragRecordClearRegionBegin(fr, dumpClear);
          unsigned int   clrEnd = getFragRecordClearRegionEnd  (fr, dumpClear);
          char          *seq = getFragRecordSequence(fr);
          char          *qlt = getFragRecordQuality(fr);
          char          *hps = getFragRecordHPS(fr);
          char          *src = getFragRecordSource(fr);

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
                  fr->gkfr.clearBeg[i], fr->gkfr.clearEnd[i]);
        }

        fprintf(stdout, "fragmentSeqOffset       = "F_U64"\n", fr->gkfr.seqOffset);
        fprintf(stdout, "fragmentQltOffset       = "F_U64"\n", fr->gkfr.qltOffset);
        fprintf(stdout, "fragmentHpsOffset       = "F_U64"\n", fr->gkfr.hpsOffset);
        fprintf(stdout, "fragmentSrcOffset       = "F_U64"\n", fr->gkfr.srcOffset);
      }
    }
  }

  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperAsFasta(char       *gkpStoreName,
                      CDS_IID_t   begIID,
                      CDS_IID_t   endIID,
                      char       *iidToDump, 
                      int         dumpAllReads,
                      int         dumpClear,
                      int         dumpQuality) {
  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  fragRecord   *fr = new_fragRecord();
  FragStream   *fs = openFragStream(gkp, (dumpQuality) ? FRAG_S_QLT : FRAG_S_SEQ);
  StoreStat     stat;

  int           i;

  statsStore(gkp->frg, &stat);

  if (begIID < stat.firstElem)
    begIID = stat.firstElem;
  if (stat.lastElem < endIID)
    endIID = stat.lastElem;

  resetFragStream(fs, begIID, endIID);

  while (nextFragStream(fs, fr)) {
    if ((iidToDump == NULL) || (iidToDump[getFragRecordIID(fr)])) {
      if (dumpAllReads || !getFragRecordIsDeleted(fr)) {
        unsigned int   clrBeg = getFragRecordClearRegionBegin(fr, dumpClear);
        unsigned int   clrEnd = getFragRecordClearRegionEnd  (fr, dumpClear);
        char          *seq = getFragRecordSequence(fr);

        if (dumpQuality)
          seq = getFragRecordQuality(fr);

        seq[clrEnd] = 0;

        fprintf(stdout, ">"F_UID","F_IID" mate="F_UID","F_IID" lib="F_UID","F_IID" clr=%s,%d,%d deleted=%d\n%s\n",
                getFragRecordUID(fr), getFragRecordIID(fr),
                (uint64)0, getFragRecordMateIID(fr),
                (uint64)0, getFragRecordLibraryIID(fr),
                AS_READ_CLEAR_NAMES[dumpClear], clrBeg, clrEnd,
                getFragRecordIsDeleted(fr),
                seq + clrBeg);
      }
    }
  }

  del_fragRecord(fr);
  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperAsOFG(char       *gkpStoreName) {
  fragRecord       *fr = new_fragRecord();
  unsigned int      firstElem = 0;
  unsigned int      lastElem = 0;

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  firstElem = getFirstElemFragStore(gkp);
  lastElem  = getLastElemFragStore(gkp) + 1;

  for (; firstElem < lastElem; firstElem++) {
    getFrag(gkp, firstElem, fr, FRAG_S_INF);

    if (getFragRecordIsDeleted(fr)) {
      printf("{OFG\n"
             "act:D\n"
             "acc:("F_UID","F_IID")\n"
             "}\n",
             getFragRecordUID(fr),
             getFragRecordIID(fr));
    } else {
      printf("{OFG\n"
             "act:A\n"
             "acc:("F_UID","F_IID")\n"
             "typ:%c\n"
             "src:\n"
             "(omit)\n.\n"
             "etm:"F_TIME_T"\n"
             "clr:"F_U32","F_U32"\n"
             "}\n",
             getFragRecordUID(fr),
             getFragRecordIID(fr),
             AS_READ,
             0,
             getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_OBT),
             getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_OBT));
    }
  }

  del_fragRecord(fr);
  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperAsFRG(char       *gkpStoreName) {
  fragRecord       *fr = new_fragRecord();
  unsigned int      firstElem = 0;
  unsigned int      lastElem = 0;

  GenericMesg       pmesg;
  FragMesg          fmesg;

  fprintf(stderr, "WARNING:  dumpGateKeeperAsFRG isn't fully implemented.\n");

  pmesg.m = &fmesg;
  pmesg.t = MESG_FRG;

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  firstElem = getFirstElemFragStore(gkp);
  lastElem  = getLastElemFragStore(gkp) + 1;

  for (; firstElem < lastElem; firstElem++) {
    getFrag(gkp, firstElem, fr, FRAG_S_ALL);

    fmesg.action        = AS_ADD;
    fmesg.eaccession    = getFragRecordUID(fr);
    fmesg.type          = AS_READ;
    fmesg.entry_time    = 0;
    fmesg.clear_rng.bgn = getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_ORIG);
    fmesg.clear_rng.end = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_ORIG);
    fmesg.source        = getFragRecordSource(fr);
    fmesg.sequence      = getFragRecordSequence(fr);
    fmesg.quality       = getFragRecordQuality(fr);
    fmesg.iaccession    = firstElem;

    WriteProtoMesg_AS(stdout, &pmesg);
  }

  del_fragRecord(fr);
  closeGateKeeperStore(gkp);
}
