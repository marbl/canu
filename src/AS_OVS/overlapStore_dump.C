
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>

#include "AS_global.H"
#include "AS_UTL_fileIO.H"
#include "AS_OVS_overlap.H"
#include "AS_OVS_overlapFile.H"
#include "AS_OVS_overlapStore.H"

#include "AS_PER_gkpStore.H"

#include "overlapStore.H"



void
dumpStore(char *ovlName, uint32 dumpBinary, double dumpERate, uint32 dumpType, uint32 dumpLength, uint32 bgnIID, uint32 endIID, uint32 qryIID) {
  OverlapStore  *ovlStore = AS_OVS_openOverlapStore(ovlName);
  OVSoverlap     overlap;
  uint64         erate     = AS_OVS_encodeQuality(dumpERate / 100.0);
  char           ovlString[1024];

  AS_OVS_setRangeOverlapStore(ovlStore, bgnIID, endIID);

  //  Length filtering is expensive to compute, need to load both reads to get their length.
  //
  //if ((dumpLength > 0) && (dumpLength < overlapLength(overlap)))
  //  continue;
  
  while (AS_OVS_readOverlapFromStore(ovlStore, &overlap, AS_OVS_TYPE_ANY) == TRUE) {
    if ((qryIID != 0) && (qryIID != overlap.b_iid))
      continue;

    switch (overlap.dat.ovl.type) {
      case AS_OVS_TYPE_OVL:
        if (overlap.dat.ovl.corr_erate > erate)
          continue;
        if (((dumpType & DUMP_5p) == 0) &&
            (overlap.dat.ovl.a_hang < 0) && (overlap.dat.ovl.b_hang < 0))
          continue;
        if (((dumpType & DUMP_3p) == 0) &&
            (overlap.dat.ovl.a_hang > 0) && (overlap.dat.ovl.b_hang > 0))
          continue;
        if (((dumpType & DUMP_CONTAINS) == 0) &&
            (overlap.dat.ovl.a_hang >= 0) && (overlap.dat.ovl.b_hang <= 0))
          continue;
        if (((dumpType & DUMP_CONTAINED) == 0) &&
            (overlap.dat.ovl.a_hang <= 0) && (overlap.dat.ovl.b_hang >= 0))
          continue;
        break;
      case AS_OVS_TYPE_OBT:
        if (overlap.dat.obt.erate > erate)
          continue;
        break;
      case AS_OVS_TYPE_MER:
        break;
      default:
        assert(0);
        break;
    }

    //  All the slow for this dump is in sprintf() within AS_OVS_toString().
    //    Without both the puts() and AS_OVS_toString(), a dump ran in 3 seconds.
    //    With both, 138 seconds.
    //    Without the puts(), 127 seconds.

    if (dumpBinary)
      AS_UTL_safeWrite(stdout, &overlap, "dumpStore", sizeof(OVSoverlap), 1);
    else
      //fprintf(stdout, "%s\n", AS_OVS_toString(ovlString, overlap));
      puts(AS_OVS_toString(ovlString, overlap));
  }

  AS_OVS_closeOverlapStore(ovlStore);
}




int
sortOBT(const void *a, const void *b) {
  OVSoverlap const *A = (OVSoverlap const *)a;
  OVSoverlap const *B = (OVSoverlap const *)b;

  if (A->dat.obt.a_beg < B->dat.obt.a_beg)  return(-1);
  if (A->dat.obt.a_beg > B->dat.obt.a_beg)  return(1);

  if (A->dat.obt.a_end < B->dat.obt.a_end)  return(-1);
  if (A->dat.obt.a_end > B->dat.obt.a_end)  return(1);

  return(0);
}



void
dumpPicture(OVSoverlap *overlaps, uint64 novl, gkStore *gkpStore, uint32 clearRegion, uint32 qryIID) {
  char           ovl[256] = {0};
  gkFragment     A;
  gkFragment     B;

  uint32         clrBgnA, clrEndA;

  uint32         MHS = 6;  //  Max Hang Size, amount of padding for "+### "

#if (AS_READ_MAX_NORMAL_LEN_BITS > 13)
  MHS = 7;
#endif

  gkpStore->gkStore_getFragment(qryIID, &A, GKFRAGMENT_INF);
  A.gkFragment_getClearRegion(clrBgnA, clrEndA, clearRegion);

  uint32  frgLenA = clrEndA - clrBgnA;

  for (int32 i=0; i<256; i++)
    ovl[i] = ' ';

  for (int32 i=0; i<100; i++)
    ovl[i + MHS] = '-';
  ovl[ 99 + MHS] = '>';
  ovl[100 + MHS] = 0;

#if (AS_READ_MAX_NORMAL_LEN_BITS > 13)
  fprintf(stdout, "%8d  A: %5d %5d                                           %s\n",
          qryIID,
          clrBgnA, clrEndA,
          ovl);
#else
  fprintf(stdout, "%8d  A: %4d %4d                                       %s\n",
          qryIID,
          clrBgnA, clrEndA,
          ovl);
#endif

  qsort(overlaps, novl, sizeof(OVSoverlap), sortOBT);


  for (uint32 o=0; o<novl; o++) {
    uint32 clrBgnB, clrEndB;

    gkpStore->gkStore_getFragment(overlaps[o].b_iid, &B, GKFRAGMENT_INF);
    B.gkFragment_getClearRegion(clrBgnB, clrEndB, clearRegion);

    uint32 frgLenB = clrEndB - clrBgnB;

    uint32 ovlBgnA = (overlaps[o].dat.obt.a_beg);
    uint32 ovlEndA = (overlaps[o].dat.obt.a_end);

    uint32 ovlBgnB = (overlaps[o].dat.obt.b_beg);
    uint32 ovlEndB = (overlaps[o].dat.obt.b_end_hi << 9) | (overlaps[o].dat.obt.b_end_lo);

    uint32 ovlStrBgn = (ovlBgnA < ovlEndA) ? ovlBgnA : ovlEndA;
    uint32 ovlStrEnd = (ovlBgnA < ovlEndA) ? ovlEndA : ovlBgnA;

    ovlStrBgn = ovlStrBgn * 100 / frgLenA + MHS;
    ovlStrEnd = ovlStrEnd * 100 / frgLenA + MHS;

    for (int32 i=0; i<256; i++)
      ovl[i] = ' ';

    for (uint32 i=ovlStrBgn; i<ovlStrEnd; i++)
      ovl[i] = '-';

    if (ovlBgnB < ovlEndB)
      ovl[ovlStrEnd-1] = '>';
    else
      ovl[ovlStrBgn] = '<';

    ovl[ovlStrEnd] = 0;

    uint32 ovlBgnHang = 0;
    uint32 ovlEndHang = 0;

    if (ovlBgnB < ovlEndB) {
      ovlBgnHang = ovlBgnB;
      ovlEndHang = frgLenB - ovlEndB;
    } else {
      ovlBgnHang = frgLenB - ovlBgnB;
      ovlEndHang = ovlEndB;
    }

    if (ovlBgnHang > 0) {
      char  str[256];
      int32 len;

      sprintf(str, "+%d", ovlBgnHang);
      len = strlen(str);

      for (int32 i=0; i<len; i++)
        ovl[ovlStrBgn - len - 1 + i] = str[i];
    }

    if (ovlEndHang > 0) {
      sprintf(ovl + ovlStrEnd, " +%d", ovlEndHang);
    }

#if (AS_READ_MAX_NORMAL_LEN_BITS > 13)
    fprintf(stdout, "%8d  A: %5d %5d (%5d)  B: %5d %5d (%5d)  %5.2f%%   %s\n",
            overlaps[o].b_iid,
            ovlBgnA, ovlEndA, frgLenA,
            ovlBgnB, ovlEndB, frgLenB,
            AS_OVS_decodeQuality(overlaps[o].dat.obt.erate) * 100.0,
            ovl);
#else
    fprintf(stdout, "%8d  A: %4d %4d (%4d)  B: %4d %4d (%4d)  %5.2f%%   %s\n",
            overlaps[o].b_iid,
            ovlBgnA, ovlEndA, frgLenA,
            ovlBgnB, ovlEndB, frgLenB,
            AS_OVS_decodeQuality(overlaps[o].dat.obt.erate) * 100.0,
            ovl);
#endif
  }
}




void
dumpPicture(char *ovlName, char *gkpName, uint32 clearRegion, double dumpERate, uint32 dumpLength, uint32 dumpType, uint32 qryIID) {
  fprintf(stderr, "DUMPING PICTURE for ID "F_IID" in store %s (gkp %s clear %s)\n",
          qryIID, ovlName, gkpName, AS_READ_CLEAR_NAMES[clearRegion]);

  OverlapStore  *ovlStore = AS_OVS_openOverlapStore(ovlName);
  gkStore       *gkpStore = new gkStore(gkpName, FALSE, FALSE);

  gkFragment     A;
  gkFragment     B;

  uint32  clrBgnA, clrEndA;

  gkpStore->gkStore_getFragment(qryIID, &A, GKFRAGMENT_INF);
  A.gkFragment_getClearRegion(clrBgnA, clrEndA, clearRegion);

  uint32  frgLenA = clrEndA - clrBgnA;

  AS_OVS_setRangeOverlapStore(ovlStore, qryIID, qryIID);

  uint64         novl     = 0;
  OVSoverlap     overlap;
  OVSoverlap    *overlaps = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * AS_OVS_numOverlapsInRange(ovlStore));
  uint64         erate    = AS_OVS_encodeQuality(dumpERate / 100.0);

  //  Load all the overlaps so we can sort by the A begin position.

  while (AS_OVS_readOverlapFromStore(ovlStore, &overlap, AS_OVS_TYPE_ANY) == TRUE) {

    //  For OBT, the only filter is erate (and length, at the bottom)

    if (overlap.dat.obt.type == AS_OVS_TYPE_OBT) {
      if (overlap.dat.obt.erate > erate)
        continue;
    }

    //  For OVL, we can also filter on overlap type.

    if (overlap.dat.obt.type == AS_OVS_TYPE_OVL) {
      if (overlap.dat.ovl.corr_erate > erate)
        continue;

      if (((dumpType & DUMP_5p) == 0) &&
          (overlap.dat.ovl.a_hang < 0) && (overlap.dat.ovl.b_hang < 0))
        continue;

      if (((dumpType & DUMP_3p) == 0) &&
          (overlap.dat.ovl.a_hang > 0) && (overlap.dat.ovl.b_hang > 0))
        continue;

      if (((dumpType & DUMP_CONTAINS) == 0) &&
          (overlap.dat.ovl.a_hang >= 0) && (overlap.dat.ovl.b_hang <= 0))
        continue;

      if (((dumpType & DUMP_CONTAINED) == 0) &&
          (overlap.dat.ovl.a_hang <= 0) && (overlap.dat.ovl.b_hang >= 0))
        continue;
    }

    //  Convert the OVL to an OBT overlap, so we can check length.  We need to do it anyway for
    //  dumpPicture(), so it isn't a horrible waste.

    if (overlap.dat.obt.type == AS_OVS_TYPE_OVL) {
      uint32  clrBgnB, clrEndB;

      gkpStore->gkStore_getFragment(overlap.b_iid, &B, GKFRAGMENT_INF);
      B.gkFragment_getClearRegion(clrBgnB, clrEndB, clearRegion);

      uint32  frgLenB = clrEndB - clrBgnB;

      int32   ahang   = overlap.dat.ovl.a_hang;
      int32   bhang   = overlap.dat.ovl.b_hang;

      uint32  abgn    = (ahang < 0) ? (0)               : (ahang);
      uint32  aend    = (bhang < 0) ? (frgLenA + bhang) : (frgLenA);
      uint32  bbgn    = (ahang < 0) ? (-ahang)          : (0);
      uint32  bend    = (bhang < 0) ? (frgLenB)         : (frgLenB - bhang);

      uint64  erate   = overlap.dat.ovl.corr_erate;

      if (overlap.dat.ovl.flipped) {
        bbgn = frgLenB - bbgn;
        bend = frgLenB - bend;
      }

      overlap.dat.obt.type     = AS_OVS_TYPE_OBT;

      overlap.dat.obt.a_beg    = abgn;
      overlap.dat.obt.a_end    = aend;
      overlap.dat.obt.b_beg    = bbgn;
      overlap.dat.obt.b_end_hi = bend >> 9;
      overlap.dat.obt.b_end_lo = bend & 0x1ff;

      overlap.dat.obt.erate    = erate;
    }

    //  Finally, check the length.

    if (overlap.dat.obt.type == AS_OVS_TYPE_OBT) {
      if ((overlap.dat.obt.b_end_hi << 9 | overlap.dat.obt.b_end_lo) - overlap.dat.obt.b_beg < dumpLength)
        continue;
      if (overlap.dat.obt.a_end - overlap.dat.obt.a_beg < dumpLength)
        continue;
    }

    overlaps[novl++] = overlap;
  }


#if 0
  //  dumpPictureOBT() (removed 21 Aug 2014) used to have this chunk of code, also ifdef'd out.  It
  //  seems to work without, but I don't know.

  for (uint32 o=0; o<novl; o++) {
    uint32 abgn = (overlaps[o].dat.obt.a_beg);
    uint32 aend = (overlaps[o].dat.obt.a_end);

    uint32 bbgn = (overlaps[o].dat.obt.b_beg);
    uint32 bend = (overlaps[o].dat.obt.b_end_hi << 9) | (overlaps[o].dat.obt.b_end_lo);

    if (overlaps[o].dat.obt.fwd) {
      //  Does nothing, just replaces with the same stuff.
      overlaps[o].dat.obt.a_beg    = abgn;
      overlaps[o].dat.obt.a_end    = aend;
      overlaps[o].dat.obt.b_beg    = bbgn;
      overlaps[o].dat.obt.b_end_hi = bend >> 9;
      overlaps[o].dat.obt.b_end_lo = bend & 0x1ff;
    } else {
      //  Reverse the coords for B to indicate a reversed overlap.
      overlaps[o].dat.obt.a_beg    = abgn;
      overlaps[o].dat.obt.a_end    = aend;
      overlaps[o].dat.obt.b_beg    = bend;
      overlaps[o].dat.obt.b_end_hi = bbgn >> 9;
      overlaps[o].dat.obt.b_end_lo = bbgn & 0x1ff;
    }
  }
#endif


  if      (novl == 0)
    fprintf(stderr, "no overlaps to show.\n");

  else if (overlaps[0].dat.ovl.type == AS_OVS_TYPE_MER)
    fprintf(stderr, "mer overlaps cannot be visualized.\n");

  else if (overlaps[0].dat.obt.type == AS_OVS_TYPE_OBT)
    dumpPicture(overlaps, novl, gkpStore, clearRegion, qryIID);

  else if (overlaps[0].dat.ovl.type == AS_OVS_TYPE_OVL)
    dumpPicture(overlaps, novl, gkpStore, clearRegion, qryIID);
}
