
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

static char const *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_GKP_include.H"
#include "AS_UTL_fasta.H"
#include "AS_UTL_reverseComplement.H"
#include "AS_PER_gkpStore.H"

#include <vector>
using namespace std;


//  Container for all files of a specific type (fasta, qual, fastq, etc).
//  Ask for the file for libIID, type.
//  Special case for 'monolithic' where libIID is ignored.

class dumpFileFile {
public:
  dumpFileFile() {
    memset(_N, 0, sizeof(char)   * 5 * FILENAME_MAX);
    memset(_F, 0, sizeof(FILE *) * 5);
  };

  char    _N[5][FILENAME_MAX];
  FILE   *_F[5];
};


class dumpFile {
public:
  dumpFile(gkStore *gkp, char *prefix, char *extension) {

    sprintf(_mono._N[0], "%s.%s",         prefix, extension);
    sprintf(_mono._N[1], "%s.1.%s",       prefix, extension);
    sprintf(_mono._N[2], "%s.2.%s",       prefix, extension);
    sprintf(_mono._N[3], "%s.paired.%s",  prefix, extension);
    sprintf(_mono._N[4], "%s.unmated.%s", prefix, extension);

    _lib.resize(gkp->gkStore_getNumLibraries() + 1);

    for (uint32 ll=0; ll<=gkp->gkStore_getNumLibraries(); ll++) {
      sprintf(_lib[ll]._N[0], "%s.%03u.%s.%s",         prefix, ll, gkp->gkStore_getLibrary(ll)->libraryName, extension);
      sprintf(_lib[ll]._N[1], "%s.%03u.%s.1.%s",       prefix, ll, gkp->gkStore_getLibrary(ll)->libraryName, extension);
      sprintf(_lib[ll]._N[2], "%s.%03u.%s.2.%s",       prefix, ll, gkp->gkStore_getLibrary(ll)->libraryName, extension);
      sprintf(_lib[ll]._N[3], "%s.%03u.%s.paired.%s",  prefix, ll, gkp->gkStore_getLibrary(ll)->libraryName, extension);
      sprintf(_lib[ll]._N[4], "%s.%03u.%s.unmated.%s", prefix, ll, gkp->gkStore_getLibrary(ll)->libraryName, extension);
    }
  };

  ~dumpFile() {
    for (uint32 ii=0; ii<5; ii++)
      if (_mono._F[ii])  fclose(_mono._F[ii]);

    for (uint32 ll=0; ll<_lib.size(); ll++)
      for (uint32 ii=0; ii<5; ii++)
        if (_lib[ll]._F[ii])  fclose(_lib[ll]._F[ii]);
  };

  FILE   *getFile(bool withLibrary, uint32 libIID, char type) {
    uint32   tt = 0;

    switch (type) {
      case 'f': tt=0; break;
      case '1': tt=1; break;
      case '2': tt=2; break;
      case 'p': tt=3; break;
      case 'u': tt=4; break;
      default:
        assert(0);
        break;
    }

    dumpFileFile  *FF = (withLibrary == false) ?  &_mono : &_lib[libIID];
    
    if (FF->_F[tt] == NULL) {
      errno = 0;
      FF->_F[tt] = fopen(FF->_N[tt], "w");

      if (errno)
        fprintf(stderr, "Failed to open '%s': %s\n", FF->_N[tt], strerror(errno)), exit(1);
    }

    return(FF->_F[tt]);
  }

public:
  dumpFileFile            _mono;
  vector<dumpFileFile>    _lib;
};



void
dumpGateKeeperInfo(char       *gkpStoreName,
                   int         asTable) {

  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE);

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
                        int         asTable) {
  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE);

  if (bgnIID < 1)
    bgnIID = 1;
  if (gkp->gkStore_getNumLibraries() < endIID)
    endIID = gkp->gkStore_getNumLibraries();

  if (asTable)
    fprintf(stdout, "IID\tOrientation\tMean\tStdDev\tNumFeatures\n");

  for (int32 i=bgnIID; i<=endIID; i++) {
    if ((iidToDump == NULL) || (iidToDump[i])) {
      gkLibrary      *gkpl = gkp->gkStore_getLibrary(i);
      LibraryMesg     lmesg;

      gkpl->gkLibrary_encodeFeatures(&lmesg);

      if (asTable) {
        fprintf(stdout, "%s\t"F_IID"\t%s\t%.3f\t%.3f\t%d\n",
                gkpl->libraryName, i,
                AS_READ_ORIENT_NAMES[gkpl->orientation],
                gkpl->mean,
                gkpl->stddev,
                lmesg.num_features);

      } else {
        fprintf(stdout, "libraryIdent         = %s,"F_IID"\n", gkpl->libraryName, i);
        fprintf(stdout, "libraryOrientation   = %s\n", AS_READ_ORIENT_NAMES[gkpl->orientation]);
        fprintf(stdout, "libraryMean          = %.3f\n", gkpl->mean);
        fprintf(stdout, "libraryStdDev        = %.3f\n", gkpl->stddev);
        fprintf(stdout, "libraryNumFeatures   = %d\n", lmesg.num_features);

        for (uint32 f=0; f<lmesg.num_features; f++)
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
                        int         asTable) {
  gkStore   *gkp = new gkStore(gkpStoreName, FALSE, FALSE);

  if (bgnIID < 1)
    bgnIID = 1;
  if (gkp->gkStore_getNumFragments() < endIID)
    endIID = gkp->gkStore_getNumFragments();

  gkFragment  fr;
  gkStream   *fs = new gkStream(gkp, bgnIID, endIID, (!dumpWithSequence || asTable) ? GKFRAGMENT_INF : GKFRAGMENT_QLT);

  //int           i;

  if (asTable)
    fprintf(stdout, "IID\tmateIID\tlibIID\tisDeleted\tisNonRandom\tOrient\tLength\tclrBegin%s\tclrEnd%s\n",
            AS_READ_CLEAR_NAMES[dumpClear], AS_READ_CLEAR_NAMES[dumpClear]);

  while (fs->next(&fr)) {
    if ((iidToDump == NULL) || (iidToDump[fr.gkFragment_getReadIID()])) {
      AS_IID     mateiid = fr.gkFragment_getMateIID();

      AS_IID     libiid = fr.gkFragment_getLibraryIID();

      uint32     clrBgn = fr.gkFragment_getClearRegionBegin(dumpClear);
      uint32     clrEnd = fr.gkFragment_getClearRegionEnd  (dumpClear);

      if (asTable) {
        fprintf(stdout, F_IID"\t"F_IID"\t"F_IID"\t%d\t%d\t%s\t%d\t%d\t%d\n",
                fr.gkFragment_getReadIID(), mateiid, libiid,
                fr.gkFragment_getIsDeleted(),
                fr.gkFragment_getIsNonRandom(),
                AS_READ_ORIENT_NAMES[fr.gkFragment_getOrientation()],
                fr.gkFragment_getSequenceLength(),
                clrBgn, clrEnd);
      } else {
        fprintf(stdout, "fragmentIdent           = "F_IID"\n", fr.gkFragment_getReadIID());
        fprintf(stdout, "fragmentMate            = "F_IID"\n", mateiid);
        fprintf(stdout, "fragmentLibrary         = "F_IID"\n", libiid);

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
                       int      doNotFixMates,
                       int      dumpAllReads) {

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

  gkStream   *fs = new gkStream(gkp, bgnIID, endIID, GKFRAGMENT_INF);
  gkFragment  fr;

  int32       mateAdded = 0;

  fprintf(stderr, "Scanning store to find libraries used and reads to dump.\n");

  while (fs->next(&fr)) {
    AS_IID  iid = fr.gkFragment_getReadIID();
    AS_IID  mid = fr.gkFragment_getMateIID();

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
                      int         withLibName,
                      AS_IID      bgnIID,
                      AS_IID      endIID,
                      char       *iidToDump,
                      int         doNotFixMates,
                      int         dumpAllReads,
                      int         dumpAllBases,
                      int         dumpClear) {

  gkStore   *gkp       = new gkStore(gkpStoreName, FALSE, FALSE);
  int32     *libToDump = NULL;

  adjustBeginEndAddMates(gkp, bgnIID, endIID, libToDump, iidToDump, doNotFixMates, dumpAllReads);

  dumpFile *fasta = new dumpFile(gkp, prefix, "fasta");
  dumpFile *qvals = new dumpFile(gkp, prefix, "fasta.qv");
  dumpFile *quals = new dumpFile(gkp, prefix, "fasta.qual");

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

    AS_UTL_writeFastA(fasta->getFile(withLibName, libIID, 'f'), seq, len, 0,
                      ">"F_IID" mate="F_IID" lib="F_IID" clr=%s,%d,%d deleted=%d\n",
                      id1, id2, libIID,
                      AS_READ_CLEAR_NAMES[dumpClear], lclr, rclr,
                      fr.gkFragment_getIsDeleted());

    AS_UTL_writeFastA(qvals->getFile(withLibName, libIID, 'f'), qlt, len, 0,
                      ">"F_IID" mate="F_IID" lib="F_IID" clr=%s,%d,%d deleted=%d\n",
                      id1, id2, libIID,
                      AS_READ_CLEAR_NAMES[dumpClear], lclr, rclr,
                      fr.gkFragment_getIsDeleted());

    AS_UTL_writeQVFastA(quals->getFile(withLibName, libIID, 'f'), qlt, len, 0,
                        ">"F_IID" mate="F_IID" lib="F_IID" clr=%s,%d,%d deleted=%d\n",
                        id1, id2, libIID,
                        AS_READ_CLEAR_NAMES[dumpClear], lclr, rclr,
                        fr.gkFragment_getIsDeleted());
  }

  safe_free(iidToDump);
  safe_free(libToDump);

  delete fasta;
  delete qvals;
  delete quals;

  delete fs;
  delete gkp;
}





void
dumpGateKeeperAsFastQ(char       *gkpStoreName,
                      char       *prefix,
                      int         withLibName,
                      AS_IID      bgnIID,
                      AS_IID      endIID,
                      char       *iidToDump,
                      int         doNotFixMates,
                      int         dumpAllReads,
                      int         dumpAllBases,
                      int         dumpClear) {

  gkStore   *gkp       = new gkStore(gkpStoreName, FALSE, FALSE);
  int32     *libToDump = NULL;

  adjustBeginEndAddMates(gkp, bgnIID, endIID, libToDump, iidToDump, doNotFixMates, dumpAllReads);

  dumpFile *fastq = new dumpFile(gkp, prefix, "fastq");

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
      AS_UTL_writeFastQ(fastq->getFile(withLibName, libIID, 'u'), seq, len, qlt, len,
                        "@%u clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                        fr.gkFragment_getReadIID(),
                        clrBgn, clrEnd, vecBgn, vecEnd, maxBgn, maxEnd, tntBgn, tntEnd,
                        fr.gkFragment_getIsNonRandom() ? 'f' : 't');
      continue;
    }


    //  We must now have a valid mate pair.  Note that a after we dump,
    //  the iidToDump[] array is updated so we never see this pair again.


    //  Write the first fragment (twice).
    AS_UTL_writeFastQ(fastq->getFile(withLibName, libIID, '1'), seq, len, qlt, len,
                      "@%u clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      fr.gkFragment_getReadIID(),
                      clrBgn, clrEnd, vecBgn, vecEnd, maxBgn, maxEnd, tntBgn, tntEnd,
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    AS_UTL_writeFastQ(fastq->getFile(withLibName, libIID, 'p'), seq, len, qlt, len,
                      "@%u clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      fr.gkFragment_getReadIID(),
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
    AS_UTL_writeFastQ(fastq->getFile(withLibName, libIID, '2'), seq, len, qlt, len,
                      "@%u clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      fr.gkFragment_getReadIID(),
                      clrBgn, clrEnd, vecBgn, vecEnd, maxBgn, maxEnd, tntBgn, tntEnd,
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    AS_UTL_writeFastQ(fastq->getFile(withLibName, libIID, 'p'), seq, len, qlt, len,
                      "@%u clr="F_U32","F_U32" clv="F_U32","F_U32" max="F_U32","F_U32" tnt="F_U32","F_U32" rnd=%c\n",
                      fr.gkFragment_getReadIID(),
                      clrBgn, clrEnd, vecBgn, vecEnd, maxBgn, maxEnd, tntBgn, tntEnd,
                      fr.gkFragment_getIsNonRandom() ? 'f' : 't');

    //  Mark the pair as dumped.  This is a cheap way around testing if we've already dumped a pair.
    //
    iidToDump[id1] = 0;
    iidToDump[id2] = 0;
  }

  safe_free(iidToDump);
  safe_free(libToDump);

  delete fastq;

  delete fs;
  delete gkp;
}



int
dumpGateKeeperIsFeatureSet(char      *gkpStoreName,
                           AS_IID     libIID,
                           char      *featureName)
{
   assert(featureName != NULL);
  
   gkStore *gkp = new gkStore(gkpStoreName, FALSE, FALSE);
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
