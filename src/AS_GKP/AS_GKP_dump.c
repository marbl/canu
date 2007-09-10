
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

static char const *rcsid = "$Id: AS_GKP_dump.c,v 1.23 2007-09-10 19:44:41 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_GKP_include.h"


//  perl's chomp is pretty nice
#define chomp(S) { char *t=S; while (*t) t++; t--; while (isspace(*t)) *t--=0; }


void
dumpGateKeeperInfo(char       *gkpStoreName) {

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  fprintf(stdout, "num fragments        = "F_S32"\n", getNumGateKeeperFragments(gkp));
  fprintf(stdout, "num libraries        = "F_S32"\n", getNumGateKeeperLibraries(gkp));

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
    fprintf(stdout, "UID\tIID\tName\tnumFrags\tnumLibs\tnumShadowLibs\n");

  for (i=begIID; i<=endIID; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      GateKeeperBatchRecord gkpb;

      getGateKeeperBatch(gkp, i, &gkpb);

      if (asTable) {
        fprintf(stdout, F_UID"\t"F_IID"\t%s\n",
                gkpb.batchUID, i,
                (gkpb.name[0]) ? gkpb.name : ".");
      } else {
        fprintf(stdout, "batchIdent   = "F_UID","F_IID"\n", gkpb.batchUID, i);
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
    fprintf(stdout, "UID\tIID\tOrientation\tMean\tStdDev\tNumFeatures\n");

  for (i=begIID; i<=endIID; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      GateKeeperLibraryRecord *gkpl = getGateKeeperLibrary(gkp, i);
      LibraryMesg              lmesg;
      int                      nf;

      AS_PER_encodeLibraryFeatures(gkpl, &lmesg);

      if (asTable) {
        fprintf(stdout, F_UID"\t"F_IID"\t%s\t%.3f\t%.3f\t%d\n",
                gkpl->libraryUID, i,
                AS_READ_ORIENT_NAMES[gkpl->orientation],
                gkpl->mean,
                gkpl->stddev,
                nf);
      } else {
        uint32 f;

        fprintf(stdout, "libraryIdent         = "F_UID","F_IID"\n", gkpl->libraryUID, i);
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
    fprintf(stdout, "UID\tIID\tmateUID\tmateIID\tlibUID\tlibIID\tisDeleted\tisNonRandom\tStatus\tOrient\tLength\tclrBegin%s\tclrEnd%s\n",
            AS_READ_CLEAR_NAMES[dumpClear], AS_READ_CLEAR_NAMES[dumpClear]);

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
                getFragRecordClearRegionBegin(fr, dumpClear),
                getFragRecordClearRegionEnd  (fr, dumpClear));
      } else {
        fprintf(stdout, "fragmentIdent           = "F_UID","F_IID"\n", getFragRecordUID(fr), getFragRecordIID(fr));
        fprintf(stdout, "fragmentMate            = "F_UID","F_IID"\n", (CDS_UID_t)0, getFragRecordMateIID(fr));
        fprintf(stdout, "fragmentLibrary         = "F_UID","F_IID"\n", (CDS_UID_t)0, getFragRecordLibraryIID(fr));

        fprintf(stdout, "fragmentIsDeleted       = %d\n", getFragRecordIsDeleted(fr));
        fprintf(stdout, "fragmentIsNonRandom     = %d\n", getFragRecordIsNonRandom(fr));
        fprintf(stdout, "fragmentStatus          = %s\n", AS_READ_STATUS_NAMES[fr->gkfr.status]);
        fprintf(stdout, "fragmentOrientation     = %s\n", AS_READ_ORIENT_NAMES[fr->gkfr.orientation]);

        fprintf(stdout, "fragmentHasVectorClear  = %d\n", fr->gkfr.hasVectorClear);
        fprintf(stdout, "fragmentHasQualityClear = %d\n", fr->gkfr.hasQualityClear);

        fprintf(stdout, "fragmentPlate           = "F_UID"\n", fr->gkfr.plateUID);
        fprintf(stdout, "fragmentPlateLocation   = %d\n", fr->gkfr.plateLocation);

        fprintf(stdout, "fragmentSeqLen          = %d\n", getFragRecordSequenceLength(fr));
        fprintf(stdout, "fragmentHPSLen          = %d\n", getFragRecordHPSLength(fr));
        fprintf(stdout, "fragmentSrcLen          = %d\n", getFragRecordSourceLength(fr));

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

  closeFragStream(fs);
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
  closeFragStream(fs);
  closeGateKeeperStore(gkp);
}


void
dumpGateKeeperAsFRG(char       *gkpStoreName,
                    int         dumpFormat,
                    CDS_IID_t   begIID,
                    CDS_IID_t   endIID,
                    char       *iidToDump,
                    int         doNotFixMates,
                    int         dumpFRGClear) {
  fragRecord       *fr = new_fragRecord();
  FragStream       *fs = NULL;
  StoreStat         stat;

  unsigned int      firstElem = 0;
  unsigned int      lastElem = 0;

  GenericMesg       pmesg;

  LibraryMesg       libMesg;
  DistanceMesg      dstMesg;
  FragMesg          frgMesg;
  LinkMesg          lnkMesg;

  int               i;

  int              *libToDump;
  CDS_UID_t        *libUID;
  CDS_UID_t        *frgUID;
  int               mateAdded = 0;

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  statsStore(gkp->frg, &stat);

  if (begIID < stat.firstElem)
    begIID = stat.firstElem;
  if (stat.lastElem < endIID)
    endIID = stat.lastElem;

  if (iidToDump == NULL) {
    iidToDump = (char *)safe_calloc(stat.lastElem+1, sizeof(char));
    for (i=begIID; i<=endIID; i++)
      iidToDump[i] = 1;
  }

  //  Pass 1: Scan the fragments, build a list of libraries to dump,
  //  also add in any mates that were omitted from the input.

  fprintf(stderr, "Scanning store to find libraries used.\n");

  statsStore(gkp->lib, &stat);

  libToDump = (int       *)safe_calloc(stat.lastElem+1, sizeof(int));
  libUID    = (CDS_UID_t *)safe_calloc(stat.lastElem+1, sizeof(CDS_UID_t));
  frgUID    = (CDS_UID_t *)safe_calloc(endIID+1, sizeof(CDS_UID_t));

  fs = openFragStream(gkp, FRAG_S_INF);
  resetFragStream(fs, begIID, endIID);

  while (nextFragStream(fs, fr)) {
    frgUID[getFragRecordIID(fr)] = getFragRecordUID(fr);

    if (iidToDump[getFragRecordIID(fr)]) {
      libToDump[getFragRecordLibraryIID(fr)]++;

      if ((getFragRecordMateIID(fr) > 0) && (iidToDump[getFragRecordMateIID(fr)] == 0)) {
        mateAdded++;
        if (doNotFixMates == 0)
          iidToDump[getFragRecordMateIID(fr)] = 1;
      }
    }
  }

  fprintf(stderr, "%sdded %d reads to maintain mate relationships.\n",
          (doNotFixMates) ? "Would have a" : "A", mateAdded);

  fprintf(stderr, "Dumping %d fragments from unknown library (version 1 has these)\n", libToDump[0]);

  for (i=1; i<=stat.lastElem; i++)
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
  for (i=1; i<=stat.lastElem; i++) {
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

  while (nextFragStream(fs, fr)) {
    FragMesg  fmesg;
    LinkMesg  lmesg;

    if (iidToDump[getFragRecordIID(fr)]) {
      pmesg.m = &fmesg;
      pmesg.t = MESG_FRG;

      //  This code used in AS_GKP_dump.c (dumpFRG), and in AS_FGB_io.c
      fmesg.action          = getFragRecordIsDeleted(fr) ? AS_DELETE : AS_ADD;
      fmesg.eaccession      = getFragRecordUID(fr);
      fmesg.library_uid     = libUID[getFragRecordLibraryIID(fr)];
      fmesg.library_iid     = getFragRecordLibraryIID(fr);
      fmesg.plate_uid       = fr->gkfr.plateUID;
      fmesg.plate_location  = fr->gkfr.plateLocation;
      fmesg.type            = AS_READ;
      fmesg.is_random       = (getFragRecordIsNonRandom(fr)) ? 0 : 1;
      fmesg.status_code     = AS_READ_STATUS_NAMES[fr->gkfr.status][0];
      fmesg.clear_rng.bgn   = getFragRecordClearRegionBegin(fr, dumpFRGClear);
      fmesg.clear_rng.end   = getFragRecordClearRegionEnd  (fr, dumpFRGClear);
      fmesg.clear_vec.bgn   = getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_VEC);
      fmesg.clear_vec.end   = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_VEC);
      fmesg.clear_qlt.bgn   = getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_QLT);
      fmesg.clear_qlt.end   = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_QLT);
      fmesg.source          = getFragRecordSource(fr);
      fmesg.sequence        = getFragRecordSequence(fr);
      fmesg.quality         = getFragRecordQuality(fr);
      fmesg.hps             = getFragRecordHPS(fr);
      fmesg.iaccession      = firstElem;

      if (fr->gkfr.hasVectorClear == 0) {
        fmesg.clear_vec.bgn   = 1;
        fmesg.clear_vec.end   = 0;
      }

      if (fr->gkfr.hasQualityClear == 0) {
        fmesg.clear_qlt.bgn   = 1;
        fmesg.clear_qlt.end   = 0;
      }

      WriteProtoMesg_AS(stdout, &pmesg);

      if ((getFragRecordMateIID(fr) > 0) &&
          (getFragRecordMateIID(fr) < getFragRecordIID(fr))) {
        pmesg.m = &lmesg;
        pmesg.t = MESG_LKG;

        lmesg.action      = AS_ADD;
        lmesg.type        = AS_MATE;
        lmesg.link_orient = AS_READ_ORIENT_NAMES[fr->gkfr.orientation][0];
        lmesg.frag1       = frgUID[getFragRecordMateIID(fr)];
        lmesg.frag2       = getFragRecordUID(fr);
        lmesg.distance    = libUID[getFragRecordLibraryIID(fr)];

        WriteProtoMesg_AS(stdout, &pmesg);
      }
    }
  }

  del_fragRecord(fr);
  closeFragStream(fs);
  closeGateKeeperStore(gkp);
}
