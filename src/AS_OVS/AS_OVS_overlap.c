
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

static const char *rcsid = "$Id: AS_OVS_overlap.c,v 1.15 2011-06-27 17:13:42 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_OVS_overlap.h"



static
int
stringSplit(char *string, char **ptrs, int ptrsLen) {
  int  ptrsFound = 0;
  int  isFirst   = 1;  //  true if this is the first letter in a word

  while (*string) {
    if ((*string != ' ') &&
        (*string != '\t')) {

      if (isFirst) {
        ptrs[ptrsFound++] = string;
        isFirst           = 0;
      }
    } else {
      *string = 0;
      isFirst = 1;
    }
    string++;
  }

  return(ptrsFound);
}



void
AS_OVS_convertOverlapMesgToOVSoverlap(OverlapMesg *omesg, OVSoverlap *ovs) {

  //  The asserts below check for encoding errors -- the dat
  //  structure only saves a small number of bits for each field,
  //  and we want to be sure we stored all the bits.

  ovs->a_iid = omesg->aifrag;
  ovs->b_iid = omesg->bifrag;
  ovs->dat.ovl.orig_erate = AS_OVS_encodeQuality(omesg->quality);
  ovs->dat.ovl.corr_erate = ovs->dat.ovl.orig_erate;
  ovs->dat.ovl.type = AS_OVS_TYPE_OVL;

  if (omesg->orientation.isNormal()) {
    ovs->dat.ovl.a_hang   = omesg->ahg;
    ovs->dat.ovl.b_hang   = omesg->bhg;
    ovs->dat.ovl.flipped  = FALSE;

    assert(ovs->dat.ovl.a_hang  == omesg->ahg);
    assert(ovs->dat.ovl.b_hang  == omesg->bhg);

  } else if (omesg->orientation.isInnie()) {
    ovs->dat.ovl.a_hang   = omesg->ahg;
    ovs->dat.ovl.b_hang   = omesg->bhg;
    ovs->dat.ovl.flipped  = TRUE;

    assert(ovs->dat.ovl.a_hang  == omesg->ahg);
    assert(ovs->dat.ovl.b_hang  == omesg->bhg);

  } else if (omesg->orientation.isOuttie()) {
    ovs->dat.ovl.a_hang   = -omesg->bhg;
    ovs->dat.ovl.b_hang   = -omesg->ahg;
    ovs->dat.ovl.flipped  = TRUE;

    assert(ovs->dat.ovl.a_hang  == -omesg->bhg);
    assert(ovs->dat.ovl.b_hang  == -omesg->ahg);

  } else if (omesg->orientation.isAnti()) {
    ovs->dat.ovl.a_hang   = -omesg->bhg;
    ovs->dat.ovl.b_hang   = -omesg->ahg;
    ovs->dat.ovl.flipped  = FALSE;

    assert(ovs->dat.ovl.a_hang  == -omesg->bhg);
    assert(ovs->dat.ovl.b_hang  == -omesg->ahg);

  } else {
    fprintf(stderr, "YIKES:  Bad overlap orientation = %d for a = %d  b = %d\n",
            omesg->orientation.toLetter(), omesg->aifrag, omesg->bifrag);
    assert(0);
  }
}



int
AS_OVS_convertOVLdumpToOVSoverlap(char *line, OVSoverlap *olap) {
  char *ptrs[16] = {0};
  int   items    = stringSplit(line, ptrs, 16);

  if ((items == 7) ||
      (items == 8)) {
    olap->a_iid              = atoi(ptrs[0]);
    olap->b_iid              = atoi(ptrs[1]);
    olap->dat.ovl.flipped    = (ptrs[2][0] == 'i') || (ptrs[2][0] == 'I');
    olap->dat.ovl.a_hang     = atoi(ptrs[3]);
    olap->dat.ovl.b_hang     = atoi(ptrs[4]);
    olap->dat.ovl.orig_erate = AS_OVS_encodeQuality(atof(ptrs[5]) / 100.0);
    olap->dat.ovl.corr_erate = AS_OVS_encodeQuality(atof(ptrs[6]) / 100.0);
    olap->dat.ovl.seed_value = (ptrs[7] == NULL) ? 0 : atoi(ptrs[7]);
    olap->dat.ovl.type       = AS_OVS_TYPE_OVL;

  } else {
    fprintf(stderr, "AS_OVS_convertOVLdumpToOVSoverlap()-- invalid line (%d items):", items);
    for (uint32 i=0; i<items; i++)
      fprintf(stderr, " %s", ptrs[i]);
    fprintf(stderr, "\n");

    return(false);
  }

  assert(olap->dat.ovl.a_hang     == atoi(ptrs[3]));
  assert(olap->dat.ovl.b_hang     == atoi(ptrs[4]));
  assert(olap->dat.ovl.orig_erate == AS_OVS_encodeQuality(atof(ptrs[5]) / 100.0));
  assert(olap->dat.ovl.corr_erate == AS_OVS_encodeQuality(atof(ptrs[6]) / 100.0));

  return(true);
}



int
AS_OVS_convertOBTdumpToOVSoverlap(char *line, OVSoverlap *olap) {
  char *ptrs[16] = {0};
  int   items    = stringSplit(line, ptrs, 16);

  if (items == 8) {
    olap->a_iid  = atoi(ptrs[0]);
    olap->b_iid  = atoi(ptrs[1]);

    olap->dat.obt.fwd      = (ptrs[2][0] == 'f');
    olap->dat.obt.a_beg    = atoi(ptrs[3]);
    olap->dat.obt.a_end    = atoi(ptrs[4]);
    olap->dat.obt.b_beg    = atoi(ptrs[5]);
    olap->dat.obt.b_end_hi = atoi(ptrs[6]) >> 9;
    olap->dat.obt.b_end_lo = atoi(ptrs[6]) & 0x1ff;
    olap->dat.obt.erate    = AS_OVS_encodeQuality(atof(ptrs[7]) / 100.0);
    olap->dat.ovl.type     = AS_OVS_TYPE_OBT;

  } else {
    fprintf(stderr, "AS_OVS_convertOBTdumpToOVSoverlap()-- invalid line (%d items):", items);
    for (uint32 i=0; i<items; i++)
      fprintf(stderr, " %s", ptrs[i]);
    fprintf(stderr, "\n");

    return(false);
  }

  assert(olap->dat.obt.a_beg == atoi(ptrs[3]));
  assert(olap->dat.obt.a_end == atoi(ptrs[4]));
  assert(olap->dat.obt.b_beg == atoi(ptrs[5]));
  assert(((olap->dat.obt.b_end_hi << 9) | (olap->dat.obt.b_end_lo)) == atoi(ptrs[6]));

  return(true);
}






char *
AS_OVS_toString(char *outstr, OVSoverlap &olap) {

  //  Even though the b_end_hi | b_end_lo is uint64 in the struct, the result
  //  of combining them doesn't appear to be 64-bit.  The cast is necessary.

  switch (olap.dat.ovl.type) {
    case AS_OVS_TYPE_OVL:
      sprintf(outstr, "%8d %8d  %c  %5"F_S64P" %5"F_S64P"  %4.2f  %4.2f",
              olap.a_iid,
              olap.b_iid,
              olap.dat.ovl.flipped ? 'I' : 'N',
              olap.dat.ovl.a_hang,
              olap.dat.ovl.b_hang,
              AS_OVS_decodeQuality(olap.dat.ovl.orig_erate) * 100.0,
              AS_OVS_decodeQuality(olap.dat.ovl.corr_erate) * 100.0);
      break;
    case AS_OVS_TYPE_OBT:
      sprintf(outstr, "%7d %7d  %c  %4"F_U64P" %4"F_U64P"  %4"F_U64P" %4"F_U64P"  %5.2f",
              olap.a_iid, olap.b_iid,
              olap.dat.obt.fwd ? 'f' : 'r',
              olap.dat.obt.a_beg,
              olap.dat.obt.a_end,
              olap.dat.obt.b_beg,
              (uint64)((olap.dat.obt.b_end_hi << 9) | (olap.dat.obt.b_end_lo)),
              AS_OVS_decodeQuality(olap.dat.obt.erate) * 100.0);
      break;
    case AS_OVS_TYPE_MER:
      sprintf(outstr, "%7d %7d  %c  "F_U64"  %4"F_U64P" %4"F_U64P"  %4"F_U64P" %4"F_U64P,
              olap.a_iid, olap.b_iid,
              olap.dat.mer.palindrome ? 'p' : (olap.dat.mer.fwd ? 'f' : 'r'),
              olap.dat.mer.compression_length,
              olap.dat.mer.a_pos,
              olap.dat.mer.b_pos,
              olap.dat.mer.k_count,
              olap.dat.mer.k_len);
      break;
    case AS_OVS_TYPE_UNS:
    default:
      break;
  }

  return(outstr);
}
