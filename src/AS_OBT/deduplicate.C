
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

const char *mainid = "$Id: deduplicate.C,v 1.15 2011-12-29 09:26:03 brianwalenz Exp $";

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
  AS_IID   libraryIID;

  uint64   matePatternLeft:1;
  uint64   isDeleted:1;

  uint64   clrbeg:AS_READ_MAX_NORMAL_LEN_BITS;
  uint64   clrlen:AS_READ_MAX_NORMAL_LEN_BITS;

  uint32   ovllen;
  uint32   ovlmax;
  olapT   *ovl;
};



fragT *
loadFragments(gkStore *gkp) {
  gkStream       *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);
  gkFragment      fr;
  fragT          *frag = new fragT [gkp->gkStore_getNumFragments() + 1];

  while (fs->next(&fr)) {
    AS_IID  iid  = fr.gkFragment_getReadIID();

    frag[iid].matePatternLeft  = 0;
    frag[iid].isDeleted        = fr.gkFragment_getIsDeleted() ? 1 : 0;
    frag[iid].clrbeg           = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_OBTINITIAL);
    frag[iid].clrlen           = fr.gkFragment_getClearRegionLength(AS_READ_CLEAR_OBTINITIAL);
    frag[iid].mateIID          = fr.gkFragment_getMateIID();
    frag[iid].libraryIID       = fr.gkFragment_getLibraryIID();

    frag[iid].ovllen           = 0;
    frag[iid].ovlmax           = 0;
    frag[iid].ovl              = 0L;
  }

  delete fs;

  return(frag);
}



void
readOverlapsAndProcessFragments(gkStore      *gkp,
                                OverlapStore *ovsprimary,
                                OverlapStore *ovssecondary,
                                uint32        errorLimit,
                                fragT        *frag) {

  gkLibrary  *libs = new gkLibrary [gkp->gkStore_getNumLibraries() + 1];

  for (uint32 i=1; i<=gkp->gkStore_getNumLibraries(); i++)
    gkp->gkStore_getLibrary(i, libs + i);


  for (OVSoverlap *ovl = readOverlap(ovsprimary, ovssecondary);
       ovl;
       ovl = readOverlap(ovsprimary, ovssecondary)) {

    if (ovl->dat.obt.fwd == 0)
      //  Dups must be forward
      continue;

    if (ovl->dat.obt.erate > errorLimit)
      //  And of good quality
      continue;

    if ((frag[ovl->a_iid].isDeleted) || (frag[ovl->b_iid].isDeleted))
      //  And not deleted already
      continue;

    if ((frag[ovl->a_iid].libraryIID != frag[ovl->b_iid].libraryIID))
      //  And in the same library
      continue;

    if (libs[frag[ovl->a_iid].libraryIID].doRemoveDuplicateReads == 0)
      //  And marked for deduplication (lib 0 is init to not dedup)
      continue;


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

    //  We can't process mated reads until we have all the overlaps for both this read and it's
    //  mate.  We save the overlaps for later consumption.
    //
    //  There is one (rare?) case that cannot be detected, when a circle is duplicated, but one
    //  duplicate forms an unmated fragment:
    //
    //  mated      ------------->     <-----
    //  fragment   ----------------->    <-- (this frag too short)
    //  fragment   ->    <------------------
    //
    if ((frag[ovl->a_iid].mateIID != 0) && (frag[ovl->b_iid].mateIID != 0))
      frag[ovl->a_iid].addOlap(ovl->b_iid,
                               (ahang >= -MATE_HANG_SLOP) && (ahang <= MATE_HANG_SLOP) && (abegdiff <= MATE_HANG_SLOP) && (bbegdiff <= MATE_HANG_SLOP),
                               (bhang >= -MATE_HANG_SLOP) && (bhang <= MATE_HANG_SLOP) && (aenddiff <= MATE_HANG_SLOP) && (benddiff <= MATE_HANG_SLOP));


    if ((frag[ovl->a_iid].mateIID != 0) || (frag[ovl->b_iid].mateIID != 0))
      //  Finally, for fragment duplicates, we require that both fragments be unmated.  We
      //  could let the b_iid fragment be mated, hoping to catch another (rare?) case, where
      //  this fragment read came from a mate read where the b_iid mate is bad:
      //
      //  mated      ------------->     <-----
      //  fragment   --------->    <--          (this frag too short, or bad quality, etc)
      //
      continue;

    //  For unmated reads, delete if it is a near perfect prefix of something else.
    //
    //  Since these are partial overlaps, we need to check both that the overlap covers about the
    //  same piece of each fragment, and that it extends to the start of each fragment.
    //
    //  To pick the longer fragment, we then want to make sure the overlap extends to the end of
    //  this fragment, and that this fragment is contained in the other.
    //
    if (frag[ovl->a_iid].mateIID == 0) {
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

  delete [] libs;
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
processMatedFragments(gkStore *gkp, fragT *frag) {
  for (uint32 iid=0; iid<=gkp->gkStore_getNumFragments(); iid++) {
    uint32 mid = frag[iid].mateIID;

    if ((frag[iid].mateIID   == 0) ||
        (frag[iid].isDeleted == 1) || (frag[iid].ovllen == 0) ||
        (frag[mid].isDeleted == 1) || (frag[mid].ovllen == 0))
      //  Not mated, or already deleted, or no overlaps.
      continue;

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

        //  If the proper overlap pattern is found, delete me.
        if ((frag[iid].ovl[i].a && frag[mid].ovl[j].b) ||
            (frag[iid].ovl[i].b && frag[mid].ovl[j].a)) {
          fprintf(reportFile, "Delete %d <-> %d DUPof %d <-> %d\n",
                  iid,
                  mid,
                  iod,
                  jod);
          duplicateMates++;
          frag[iid].isDeleted = 1;
          frag[mid].isDeleted = 1;
          i = frag[iid].ovllen;
          j = frag[mid].ovllen;
        }
      }
    }
  }
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

      //  The cache is not enabled, as we don't expect many changes to the store.
      gkp->gkStore_metadataCaching(true);

    } else if (strncmp(argv[arg], "-ovs", 2) == 0) {
      if (ovsprimary == NULL)
        ovsprimary = AS_OVS_openOverlapStore(argv[++arg]);
      else if (ovssecondary == NULL)
        ovssecondary = AS_OVS_openOverlapStore(argv[++arg]);
      else {
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
    fragT  *frag = loadFragments(gkp);

    readOverlapsAndProcessFragments(gkp, ovsprimary, ovssecondary, errorLimit, frag);
    processMatedFragments(gkp, frag);

    if (doUpdate)
      deleteFragments(gkp, frag);

    delete [] frag;
  }

  delete    gkp;

  if (summaryFile) {
    fprintf(summaryFile, "duplicateFrags:    "F_U32"\n", duplicateFrags);
    fprintf(summaryFile, "duplicateMates:    "F_U32"\n", duplicateMates);
  }

  exit(0);
}
