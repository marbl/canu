
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

/* $Id: AS_GKP_dump.c,v 1.6 2007-03-05 23:42:11 brianwalenz Exp $ */

#include "AS_GKP_include.h"

void
dumpGateKeeperAsCompatible(char       *gkpStoreName,
                           CDS_IID_t   begIID,
                           CDS_IID_t   endIID) {
  fprintf(stderr, "Not yet.  Sorry.\n");
}



void
dumpGateKeeperAsFasta(char       *gkpStoreName,
                      CDS_IID_t   begIID,
                      CDS_IID_t   endIID,
                      int         dumpFastaAllReads,
                      int         dumpFastaAllBases,
                      int         dumpFastaClear) {
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

  if (firstElem < begIID)
    firstElem = begIID;
  if (endIID < lastElem)
    lastElem = endIID;

  for (; firstElem < lastElem; firstElem++) {
    unsigned int      deleted = 0;

    getFrag(gkp, firstElem, fr, FRAG_S_SEQ | FRAG_S_QLT);

    deleted = getFragRecordIsDeleted(fr);

    if (dumpFastaAllReads || !deleted) {
      CDS_UID_t      readUID = getFragRecordUID(fr);
      CDS_UID_t      mateUID = 0;
      CDS_UID_t      libUID  = 0;

      CDS_IID_t      readIID = getFragRecordIID(fr);
      CDS_IID_t      mateIID = getFragRecordMateIID(fr);
      CDS_IID_t      libIID  = getFragRecordLibraryIID(fr);

      unsigned int   clrBeg = getFragRecordClearRegionBegin(fr, dumpFastaClear);
      unsigned int   clrEnd = getFragRecordClearRegionEnd  (fr, dumpFastaClear);

      char          *seq = getFragRecordSequence(fr);
      char          *qlt = getFragRecordQuality(fr);

      if (dumpFastaAllBases) {
        fprintf(stdout, ">"F_UID","F_IID" mate="F_UID","F_IID" lib="F_UID","F_IID" deleted=%d clrBeg="F_U32" clrEnd="F_U32"\n%s\n",
                readUID, readIID,
                mateUID, mateIID,
                libUID,  libIID,
                deleted,
                clrBeg, clrEnd,
                seq);
        fprintf(stdout, ">"F_UID","F_IID" mate="F_UID","F_IID" lib="F_UID","F_IID" deleted=%d clrBeg="F_U32" clrEnd="F_U32"\n%s\n",
                readUID, readIID,
                mateUID, mateIID,
                libUID,  libIID,
                deleted,
                clrBeg, clrEnd,
                qlt);
      } else {
        seq[clrEnd] = 0;
        fprintf(stdout, ">"F_UID","F_IID" mate="F_UID","F_IID" lib="F_UID","F_IID" deleted=%d\n%s\n",
                readUID, readIID,
                mateUID, mateIID,
                libUID,  libIID,
                deleted,
                seq + clrBeg);
        fprintf(stdout, ">"F_UID","F_IID" mate="F_UID","F_IID" lib="F_UID","F_IID" deleted=%d\n%s\n",
                readUID, readIID,
                mateUID, mateIID,
                libUID,  libIID,
                deleted,
                qlt + clrBeg);
      }
    }
  }

  del_fragRecord(fr);
  closeGateKeeperStore(gkp);
}





void
dumpGateKeeperAsXML(char       *gkpStoreName,
                    CDS_IID_t   begIID,
                    CDS_IID_t   endIID) {

  int dumpBatchInfo    = 1;
  int dumpLibraryInfo  = 1;
  int dumpFragmentInfo = 1;

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  if (dumpBatchInfo) {
    StoreStat stat;
    int       i;

    statsStore(gkp->bat, &stat);
     
    for (i=stat.firstElem; i<=stat.lastElem; i++) {
      GateKeeperBatchRecord gkpb;
      time_t                createdtime;

      getGateKeeperBatchStore(gkp->bat, i, &gkpb);

      createdtime = (time_t)gkpb.created;

      fprintf(stdout, "\n\n");
      fprintf(stdout, "<batch>\n");
      fprintf(stdout, "<ident>"F_UID","F_IID"</ident>\n", gkpb.UID, i);
      fprintf(stdout, "<name>%s</name>\n", gkpb.name);
      fprintf(stdout, "<created>("F_U64") %s</created>\n", gkpb.created, ctime(&createdtime));
      if (gkpb.comment[0])
        fprintf(stdout, "<comment>\n%s\n</comment>\n", gkpb.comment);
      fprintf(stdout, "<nFrags>"F_S32"</nFrags>\n", gkpb.numFragments);
      fprintf(stdout, "<nLibs>"F_S32"</nLibs)\n", gkpb.numLibraries);
      fprintf(stdout, "<nLibsS>"F_S32"</nLibsS>\n", gkpb.numLibraries_s);
      fprintf(stdout, "</batch>\n");
    }

    fprintf(stdout, "\n\n");
    fprintf(stdout, "<total>\n");
    fprintf(stdout, "<nFragsTotal>"F_S32"</nFragsTotal>\n", getNumGateKeeperFragments(gkp->frg));
    fprintf(stdout, "<nLibsTotal>"F_S32"</nLibsTotal>\n", getNumGateKeeperLibrarys(gkp->lib));
    fprintf(stdout, "<nLibsSTotal>"F_S32"</nLibsSTotal>\n", getNumGateKeeperLibrarys(gkp->lis));
    fprintf(stdout, "</total>\n");
  }


  if (dumpLibraryInfo) {
    StoreStat stat;
    int       i;

    statsStore(gkp->lib, &stat);

    for (i=stat.firstElem; i<=stat.lastElem; i++) {
      GateKeeperLibraryRecord gkpl;
      time_t                  createdtime;

      getGateKeeperLibraryStore(gkp->lib, i, &gkpl);

      createdtime = (time_t)gkpl.created;

      fprintf(stdout, "\n\n");
      fprintf(stdout, "<library>\n");
      fprintf(stdout, "<ident>"F_UID","F_IID"</ident>\n", gkpl.UID, i);
      fprintf(stdout, "<name>%s</name>\n", gkpl.name);
      fprintf(stdout, "<created>("F_U64") %s</created>\n", gkpl.created, ctime(&createdtime));
      if (gkpl.comment[0])
        fprintf(stdout, "<comment>\n%s\n</comment>\n", gkpl.comment);
      fprintf(stdout, "<deleted>%d</deleted>\n", gkpl.deleted);
      fprintf(stdout, "<redefined>%d</redefined>\n", gkpl.redefined);
      fprintf(stdout, "<orientation>%d</orientation>\n", gkpl.orientation);
      fprintf(stdout, "<mean>%f</mean>\n", gkpl.mean);
      fprintf(stdout, "<stddev>%f</stddev>\n", gkpl.stddev);
      fprintf(stdout, "<numFeatures>%d</numFeatures>\n", gkpl.numFeatures);
      fprintf(stdout, "<prevInstanceID>"F_IID"</prevInstanceID>\n", gkpl.prevInstanceID);
      fprintf(stdout, "<prevID>"F_IID"</prevID>\n", gkpl.prevID);
      fprintf(stdout, "<birthBatch>%d</birthBatch>\n", gkpl.birthBatch);
      fprintf(stdout, "<deathBatch>%d</deathBatch>\n", gkpl.deathBatch);
      fprintf(stdout, "</library>\n");
    }
       
    statsStore(gkp->lis, &stat);

    for (i=stat.firstElem; i<=stat.lastElem; i++) {
      GateKeeperLibraryRecord gkpl;
      time_t                  createdtime;

      getGateKeeperLibraryStore(gkp->lis, i, &gkpl);

      createdtime = (time_t)gkpl.created;

      fprintf(stdout, "\n\n");
      fprintf(stdout, "<shadow_library>\n");
      fprintf(stdout, "<ident>"F_UID","F_IID"</ident>\n", gkpl.UID, i);
      fprintf(stdout, "<name>%s</name>\n", gkpl.name);
      fprintf(stdout, "<created>("F_U64") %s</created>\n", gkpl.created, ctime(&createdtime));
      if (gkpl.comment[0])
        fprintf(stdout, "<comment>\n%s\n</comment>\n", gkpl.comment);
      fprintf(stdout, "<deleted>%d</deleted>\n", gkpl.deleted);
      fprintf(stdout, "<redefined>%d</redefined>\n", gkpl.redefined);
      fprintf(stdout, "<orientation>%d</orientation>\n", gkpl.orientation);
      fprintf(stdout, "<mean>%f</mean>\n", gkpl.mean);
      fprintf(stdout, "<stddev>%f</stddev>\n", gkpl.stddev);
      fprintf(stdout, "<numFeatures>%d</numFeatures>\n", gkpl.numFeatures);
      fprintf(stdout, "<prevInstanceID>"F_IID"</prevInstanceID>\n", gkpl.prevInstanceID);
      fprintf(stdout, "<prevID>"F_IID"</prevID>\n", gkpl.prevID);
      fprintf(stdout, "<birthBatch>%d</birthBatch>\n", gkpl.birthBatch);
      fprintf(stdout, "<deathBatch>%d</deathBatch>\n", gkpl.deathBatch);
      fprintf(stdout, "</shadow_library>\n");
    }
  }


  if (dumpFragmentInfo) {
    StreamHandle             frags = openStream(gkp->frg, NULL, 0);
    GateKeeperFragmentRecord gkpf;
    StoreStat                stat;
    int                      i;

    statsStore(gkp->frg, &stat);

    if (begIID < stat.firstElem)
      begIID = stat.firstElem;
    if (endIID > stat.lastElem)
      endIID = stat.lastElem;

    resetStream(frags, begIID, endIID);
    i = begIID;

    while (nextStream(frags, &gkpf)) {

      fprintf(stdout, "\n\n");
      fprintf(stdout, "<fragment>\n");
      fprintf(stdout, "<ident>"F_UID","F_IID"</ident>\n", gkpf.readUID, gkpf.readIID);
      fprintf(stdout, "<mate>"F_UID","F_IID"</mate>\n", (CDS_UID_t)0, gkpf.mateIID);
      fprintf(stdout, "<plate>"F_UID"</plate>\n", gkpf.plateUID);
      fprintf(stdout, "<library>"F_IID"</library>\n", gkpf.libraryIID);

      //  XXX  flags

      fprintf(stdout, "<sequence>\n%s\n</sequence>\n", "(notyet)");
      fprintf(stdout, "<quality>\n%s\n</quality>\n", "(notyet)");

      for (i=0; i<AS_READ_CLEAR_NUM; i++) {
        fprintf(stdout, "<clear_%s>%d,%d</clear_%s>\n", 
                AS_READ_CLEAR_NAMES[i],
                gkpf.clearBeg[i], gkpf.clearEnd[i],
                AS_READ_CLEAR_NAMES[i]);
      }

      fprintf(stdout, "<seqOffset>"F_U64"</seqOffset)\n", gkpf.seqOffset);
      fprintf(stdout, "<qltOffset>"F_U64"</qltOffset)\n", gkpf.qltOffset);
      fprintf(stdout, "<hpsOffset>"F_U64"</hpsOffset)\n", gkpf.hpsOffset);
      fprintf(stdout, "<srcOffset>"F_U64"</srcOffset)\n", gkpf.srcOffset);

      fprintf(stdout, "</fragment>\n");

      i++;
    }
  }
}





void
dumpGateKeeperAsOFG(char *gkpStoreName) {
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
dumpGateKeeperAsFRG(char *gkpStoreName) {
  fragRecord       *fr = new_fragRecord();
  unsigned int      firstElem = 0;
  unsigned int      lastElem = 0;

  GenericMesg       pmesg;
  FragMesg          fmesg;

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
