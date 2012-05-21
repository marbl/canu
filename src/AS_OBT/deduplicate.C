
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

const char *mainid = "$Id: deduplicate.C,v 1.19 2012-05-21 04:51:47 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "readOverlap.H"

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"

#define F_U32W(X)  "%" #X F_U32P
#define F_U64W(X)  "%" #X F_U64P

#define FRAG_HANG_SLOP  0
#define MATE_HANG_SLOP  0
#define DEFAULT_ERATE   2.0 / 100.0

FILE    *summaryFile = stderr;
FILE    *reportFile  = stdout;

uint32   duplicateFrags  = 0;
uint32   duplicateMates  = 0;

uint32   mateOvlTypes[4][4] = { {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} };


//  Stores the overlap to some fragment.  The overlap is to the biid fragment.  If the overlap
//  starts at the start of both fragments, 'a' is set.  If the overlap ends at the end of both
//  fragments, 'b' is set.
//
//  'a' : -------->            'b' : -------->
//        -------------->               ----->
//
class olapT {
public:
  AS_IID   biid;
  uint32   a:1;
  uint32   b:1;
};

class fragT {
public:

  void addOlap(uint32 biid, uint32 ahang, uint32 bhang) {
    if (ovlmax == 0) {
      ovlmax = 4;
      ovl = new olapT [ovlmax];
    }

    if (ovllen >= ovlmax) {
      ovlmax *= 2;
      olapT *O = new olapT [ovlmax];
      memcpy(O, ovl, sizeof(olapT) * ovllen);
      delete [] ovl;
      ovl = O;
    }

    ovl[ovllen].biid = biid;
    ovl[ovllen].a    = ahang ? 1 : 0;
    ovl[ovllen].b    = bhang ? 1 : 0;

    ovllen++;
  };

  AS_IID   mateIID;

  //  This could be a 32-bit value, for BITS=11 (leaving 8 bits for library IID), but
  //  any larger value of BITS reduces the library IID too much.  By BITS=15, there is
  //  no space at all.  So we just eat the extra memory and leave it 64-bit wide.

  uint64   libraryIID      : 64 - 2 * AS_READ_MAX_NORMAL_LEN_BITS - 2;
  uint64   matePatternLeft : 1;
  uint64   isDeleted       : 1;

  uint64   clrbeg          : AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   clrlen          : AS_READ_MAX_NORMAL_LEN_BITS;

  uint32   ovllen;
  uint32   ovlmax;
  olapT   *ovl;
};



fragT *
loadFragments(gkStore *gkp) {
  gkStream       *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);
  gkFragment      fr;
  fragT          *frag = new fragT [gkp->gkStore_getNumFragments() + 1];

  fprintf(stderr, "loading fragment data\n");

  while (fs->next(&fr)) {
    AS_IID  iid  = fr.gkFragment_getReadIID();

    frag[iid].matePatternLeft  = 0;
    frag[iid].isDeleted        = fr.gkFragment_getIsDeleted() ? 1 : 0;

    if (fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTINITIAL) < fr.gkFragment_getClearRegionEnd(AS_READ_CLEAR_OBTINITIAL)) {
      frag[iid].clrbeg = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTINITIAL);
      frag[iid].clrlen = fr.gkFragment_getClearRegionLength(AS_READ_CLEAR_OBTINITIAL);
    } else {
      frag[iid].clrbeg = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_CLR);
      frag[iid].clrlen = fr.gkFragment_getClearRegionLength(AS_READ_CLEAR_CLR);
    }

    frag[iid].mateIID    = fr.gkFragment_getMateIID();
    frag[iid].libraryIID = fr.gkFragment_getLibraryIID();

    frag[iid].ovllen = 0;
    frag[iid].ovlmax = 0;
    frag[iid].ovl    = 0L;
  }

  delete fs;

  return(frag);
}






uint32
loadOverlaps(OverlapStore *store,
             uint32        ovlLen,
             uint32        ovlMax,
             OVSoverlap   *ovlBuffer,
             AS_IID       &last) {
  uint32  nr = 0;

  //  If there is space left, load overlaps.  If any were loaded, remember the ID if the first one
  //  loaded.  We limit mate processing to be strictly less than this ID.

  if (ovlLen < ovlMax) {
    nr = AS_OVS_readOverlapsFromStore(store, ovlBuffer + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_OBT, false);

    if (nr > 0)
      last = ovlBuffer[ovlLen].a_iid;
    else
      last = AS_IID_MAX;

    ovlLen += nr;
  }

  //fprintf(stderr, "loadOverlaps()-- read %8"F_U32P" overlaps from store 0x%016p (now %8"F_U32P" overlaps)\n",
  //        nr, store, ovlLen);

  return(ovlLen);
}



//  Process all the mated overlaps.
//
//  Delete if one read has a near zero ahang, and the other read has a near
//  zero bhang.
//
//  i   ----------->      <----------  j
//  iod ----->         <-------------  jod
//
//  If imd (the mate of iod) is the same as jod (and implicitly jmd == iod), then
//  we can delete one of these mates.  Which?
//
void
processMatedFragment(gkStore *gkp, fragT *frag, AS_IID iid, AS_IID mid) {

  for (uint32 i=0; i<frag[iid].ovllen; i++) {
    uint32  iod = frag[iid].ovl[i].biid;   //  Read IID of my overlapping fragment
    uint32  imd = frag[iod].mateIID;       //  Mate IID of that overlapping framgnet

    if ((frag[iod].isDeleted == 1) ||
        (imd == -1) ||
        (frag[imd].isDeleted == 1))
      //  Overlapping frag already deleted, or no mate, or mate already deleted
      continue;

    for (uint32 j=0; j<frag[mid].ovllen; j++) {
      uint32  jod = frag[mid].ovl[j].biid;
      uint32  jmd = frag[jod].mateIID;

      if (imd != jod)
        //  Mate of my overlapping fragment is different than the overlapping fragment of my mate
        continue;

      {
        uint32 a = frag[iid].ovl[i].a + frag[iid].ovl[i].b * 2;
        uint32 b = frag[mid].ovl[j].a + frag[mid].ovl[j].b * 2;

        mateOvlTypes[a][b]++;
      }

      //  If the proper overlap pattern is found, delete me.
      //
      if (frag[iid].ovl[i].a && frag[mid].ovl[j].a) {
        fprintf(reportFile, "Delete %d <-> %d DUPof %d <-> %d %d%d%d%d\n",
                iid,
                mid,
                iod,
                jod,
                frag[iid].ovl[i].a,
                frag[iid].ovl[i].b,
                frag[mid].ovl[j].a,
                frag[mid].ovl[j].b);
        duplicateMates++;
        frag[iid].isDeleted = 1;
        frag[mid].isDeleted = 1;
        i = frag[iid].ovllen;
        j = frag[mid].ovllen;
      }
    }
  }

  delete [] frag[iid].ovl;  frag[iid].ovl = NULL;
  delete [] frag[mid].ovl;  frag[mid].ovl = NULL;
}



AS_IID
processMatedFragmentsInline(gkStore *gkp, fragT *frag, AS_IID fid, AS_IID current) {

  for (; fid < current; fid++) {
    AS_IID mid = frag[fid].mateIID;

    if (mid == 0)
      //  fid is not mated.
      continue;

    if (fid < mid)
      //  we haven't seen the mate overlaps yet
      continue;

    processMatedFragment(gkp, frag, fid, mid);
  }

  return(fid);
}


AS_IID
processOverlaps(gkStore      *gkp,
                gkLibrary   **libs,
                OVSoverlap   *ovlBuffer,
                uint32        ovlLen,
                uint32        errorLimit,
                fragT        *frag,
                AS_IID        currM,
                AS_IID        lastP,
                AS_IID        lastS) {
  uint64  nRev = 0;
  uint64  nLQ  = 0;
  uint64  nDel = 0;
  uint64  nLib = 0;
  uint64  nNoD = 0;

  uint64  nMate = 0;
  uint64  nMaFr = 0;
  uint64  nFrag = 0;

  fprintf(stderr, "processing "F_U32" overlaps currM="F_IID" lastP="F_IID" lastS="F_IID" ",
          ovlLen, currM, lastP, lastS);

  for (uint32 oo=0; oo<ovlLen; oo++) {
    OVSoverlap *ovl = ovlBuffer + oo;

    if (ovl->dat.obt.fwd == 0) {
      //  Dups must be forward
      nRev++;
      continue;
    }

    if (ovl->dat.obt.erate > errorLimit) {
      //  And of good quality
      nLQ++;
      continue;
    }

    if ((frag[ovl->a_iid].isDeleted) || (frag[ovl->b_iid].isDeleted)) {
      //  And not deleted already
      nDel++;
      continue;
    }

    if ((frag[ovl->a_iid].libraryIID != frag[ovl->b_iid].libraryIID)) {
      //  And in the same library
      nLib++;
      continue;
    }

    if (libs[frag[ovl->a_iid].libraryIID]->doRemoveDuplicateReads == 0) {
      //  And marked for deduplication (lib 0 is init to not dedup)
      nNoD++;
      continue;
    }

    int32 ab = ovl->dat.obt.a_beg;
    int32 ae = ovl->dat.obt.a_end;
    int32 bb = ovl->dat.obt.b_beg;
    int32 be = (ovl->dat.obt.b_end_hi << 9) | (ovl->dat.obt.b_end_lo);

    int32  abeg     = ab + frag[ovl->a_iid].clrbeg;
    int32  bbeg     = bb + frag[ovl->b_iid].clrbeg;
    int32  ahang    = bbeg - abeg;
    int32  abegdiff = ab;
    int32  bbegdiff = bb;

    int32  aend     = ae + frag[ovl->a_iid].clrbeg;
    int32  bend     = be + frag[ovl->b_iid].clrbeg;
    int32  bhang    = bend - aend;
    int32  aenddiff = frag[ovl->a_iid].clrlen - ae;
    int32  benddiff = frag[ovl->b_iid].clrlen - be;

    double error    = AS_OVS_decodeQuality(ovl->dat.obt.erate);

    bool   isMateA  = (frag[ovl->a_iid].mateIID != 0);
    bool   isMateB  = (frag[ovl->b_iid].mateIID != 0);


    //  No duplicates between mated and unmated fragments.
    //
    //  There is one (rare?) case that cannot be detected, when a circle is duplicated, but one
    //  duplicate forms an unmated fragment:
    //
    //  mated      ------------->     <-----
    //  fragment   ----------------->    <-- (this frag too short)
    //  fragment   ->    <------------------
    //
    if (isMateA != isMateB) {
      nMaFr++;
    }

    //  If both reads in the overlap are mated, save the overlap for later processing.  We can't
    //  process mated reads until we have all the overlaps for both this read and it's mate.
    //
    else if ((isMateA == true) && (isMateB == true)) {
      nMate++;

      frag[ovl->a_iid].addOlap(ovl->b_iid,
                               (ahang >= -MATE_HANG_SLOP) && (ahang <= MATE_HANG_SLOP) && (abegdiff <= MATE_HANG_SLOP) && (bbegdiff <= MATE_HANG_SLOP),
                               (bhang >= -MATE_HANG_SLOP) && (bhang <= MATE_HANG_SLOP) && (aenddiff <= MATE_HANG_SLOP) && (benddiff <= MATE_HANG_SLOP));
    }

    //  For unmated reads, delete if it is a near perfect prefix of something else.
    //
    //  Since these are partial overlaps, we need to check both that the overlap covers about the
    //  same piece of each fragment, and that it extends to the start of each fragment.
    //
    //  To pick the longer fragment, we then want to make sure the overlap extends to the end of
    //  this fragment, and that this fragment is contained in the other.
    //
    else {
      nFrag++;

      assert(isMateA == false);
      assert(isMateB == false);

      if ((ahang >= -FRAG_HANG_SLOP) && (ahang <= FRAG_HANG_SLOP) &&
          (abegdiff <= FRAG_HANG_SLOP) &&
          (bbegdiff <= FRAG_HANG_SLOP) &&
          (aenddiff <= FRAG_HANG_SLOP) && (bhang >= 0) &&
          (error    <= 0.025)) {
        fprintf(reportFile, "Delete %u DUPof %u  a %d,%d  b %d,%d  hang %d,%d  diff %d,%d  error %f\n",
                ovl->a_iid,
                ovl->b_iid,
                abeg, aend,
                bbeg, bend,
                ahang, bhang,
                abegdiff, bbegdiff,
                error);
        duplicateFrags++;
        frag[ovl->a_iid].isDeleted = 1;
      }
    }
  }

  fprintf(stderr, "-- nMate="F_U64" nMaFR="F_U64" nFrag="F_U64" -- nRev="F_U64" nLQ="F_U64" nDel="F_U64" nLib="F_U64" nNoD="F_U64"\n",
          nMate, nMaFr, nFrag, nRev, nLQ, nDel, nLib, nNoD);

  currM = processMatedFragmentsInline(gkp, frag, currM, MIN(lastP, lastS));

  return(currM);
}






//  Update the gkpStore with any deletions.  This is a special case -- if a frag is deleted, its
//  mate (yes, if mated) is always deleted.  We need to update our local deleted cache though.
//
void
deleteFragments(gkStore *gkp, fragT *frag) {
  for (uint32 iid=0; iid<=gkp->gkStore_getNumFragments(); iid++) {
    if (frag[iid].isDeleted) {
      uint32 mid = frag[iid].mateIID;

      gkp->gkStore_delFragment(iid, true);
      frag[iid].isDeleted = 0;
      frag[mid].isDeleted = 0;
    }
  }
}



int
main(int argc, char **argv) {
  uint32             errorLimit   = AS_OVS_encodeQuality(DEFAULT_ERATE);
  gkStore           *gkp          = 0L;
  OverlapStore      *ovsprimary   = 0L;
  OverlapStore      *ovssecondary = 0L;

  bool               doUpdate     = true;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-gkp", 2) == 0) {
      gkp = new gkStore(argv[++arg], FALSE, doUpdate);
      gkp->gkStore_metadataCaching(true);
      fprintf(stderr, "gkpStore opened.\n");

    } else if (strncmp(argv[arg], "-ovs", 2) == 0) {
      if (ovsprimary == NULL) {
        ovsprimary = AS_OVS_openOverlapStore(argv[++arg]);
        fprintf(stderr, "primary opened.\n");
      } else if (ovssecondary == NULL) {
        ovssecondary = AS_OVS_openOverlapStore(argv[++arg]);
        fprintf(stderr, "secondary opened.\n");
      } else {
        fprintf(stderr, "Only two obtStores allowed.\n");
        err++;
      }

    } else if (strncmp(argv[arg], "-erate", 2) == 0) {
      double erate = atof(argv[++arg]);
      if (erate >= AS_MAX_ERROR_RATE)
        fprintf(stderr, "Error rate %s too large; must be 'fraction error' and below %f\n", argv[arg], AS_MAX_ERROR_RATE), exit(1);
      errorLimit = AS_OVS_encodeQuality(erate);

    } else if (strncmp(argv[arg], "-summary", 2) == 0) {
      ++arg;
      errno = 0;
      summaryFile = fopen(argv[arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for writing: %s\n", argv[arg], strerror(errno)), exit(1);

    } else if (strncmp(argv[arg], "-report", 2) == 0) {
      ++arg;
      errno = 0;
      reportFile  = fopen(argv[arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for writing: %s\n", argv[arg], strerror(errno)), exit(1);

    } else if (strncmp(argv[arg], "-n", 2) == 0) {
      doUpdate = false;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkp == 0L) || (ovsprimary == 0L) || (err)) {
    fprintf(stderr, "usage: %s [-1] -gkp <gkpStore> -ovs <ovsStore> [opts]\n", argv[0]);
    fprintf(stderr, "  -erate E        filter overlaps above this fraction error; default 0.015 (== 1.5%% error)\n");
    fprintf(stderr, "  -summary S      write a summary of the fixes to S\n");
    fprintf(stderr, "  -report R       write a detailed report of the fixes to R\n");
    exit(1);
  }

  bool            nothingToDo = true;

  for (uint32 i=1; i<=gkp->gkStore_getNumLibraries(); i++) {
    gkLibrary  *gkl = gkp->gkStore_getLibrary(i);

    if (gkl->doRemoveDuplicateReads == true) {
      if (summaryFile)
        fprintf(summaryFile, "Checking library %s for duplicates.\n", gkl->libraryName);
      nothingToDo = false;
    } else {
      if (summaryFile)
        fprintf(summaryFile, "Ignoring library %s.\n", gkl->libraryName);
    }
  }

  if (nothingToDo == false) {
    fragT        *frag = loadFragments(gkp);
    gkLibrary   **libs = new gkLibrary * [gkp->gkStore_getNumLibraries() + 1];

    for (uint32 i=1; i<=gkp->gkStore_getNumLibraries(); i++)
      libs[i] = gkp->gkStore_getLibrary(i);

    AS_IID        currMate      = 0;
    AS_IID        lastPrimary   = 0;
    AS_IID        lastSecondary = 0;

    uint32        ovlLen    = 0;
    uint32        ovlMax    = 4 * 1024 * 1024;
    OVSoverlap   *ovlBuffer = new OVSoverlap [ovlMax];

    ovlLen = 0;
    ovlLen = loadOverlaps(ovsprimary,   ovlLen, ovlMax/2, ovlBuffer, lastPrimary);
    ovlLen = loadOverlaps(ovssecondary, ovlLen, ovlMax,   ovlBuffer, lastSecondary);

    while (ovlLen > 0) {
      currMate = processOverlaps(gkp, libs, ovlBuffer, ovlLen, errorLimit, frag, currMate, lastPrimary, lastSecondary);

      if (lastPrimary < lastSecondary) {
        ovlLen = 0;
        ovlLen = loadOverlaps(ovsprimary,   ovlLen, ovlMax, ovlBuffer, lastPrimary);
        ovlLen = loadOverlaps(ovssecondary, ovlLen, ovlMax, ovlBuffer, lastSecondary);
      } else {
        ovlLen = 0;
        ovlLen = loadOverlaps(ovssecondary, ovlLen, ovlMax, ovlBuffer, lastSecondary);
        ovlLen = loadOverlaps(ovsprimary,   ovlLen, ovlMax, ovlBuffer, lastPrimary);
      }
    }

    currMate = processMatedFragmentsInline(gkp, frag, currMate, gkp->gkStore_getNumFragments() + 2);

    if (doUpdate)
      deleteFragments(gkp, frag);

    delete [] ovlBuffer;
    delete [] frag;
  }

  delete    gkp;

  if (summaryFile) {
    fprintf(summaryFile, "duplicateFrags:    "F_U32"\n", duplicateFrags);
    fprintf(summaryFile, "duplicateMates:    "F_U32"\n", duplicateMates);
  }

  if (summaryFile) {
    char *label[4] = { "~a~b", "a~b", "~ab", "ab" };

    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "\n");
    fprintf(summaryFile, "Duplicate Mate Overlap Types\n");
    fprintf(summaryFile, "\t~a~b\ta~b\t~ab\tab\n");

    for (uint32 i=0; i<4; i++) {
      fprintf(summaryFile, "%s", label[i]);

      for (uint32 j=0; j<4; j++)
        fprintf(summaryFile, "\t"F_U32, mateOvlTypes[i][j]);

      fprintf(summaryFile, "\n");
    }
  }

  exit(0);
}
