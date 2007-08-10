
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

static char CM_ID[] = "$Id: AS_OVS_overlap.c,v 1.4 2007-08-10 06:47:14 brianwalenz Exp $";

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

  switch (omesg->orientation) {
    case  AS_NORMAL:
      ovs->dat.ovl.a_hang   = omesg->ahg;
      ovs->dat.ovl.b_hang   = omesg->bhg;
      ovs->dat.ovl.flipped  = FALSE;
          
      assert(ovs->dat.ovl.a_hang  == omesg->ahg);
      assert(ovs->dat.ovl.b_hang  == omesg->bhg);
      break;

    case  AS_INNIE:
      ovs->dat.ovl.a_hang   = omesg->ahg;
      ovs->dat.ovl.b_hang   = omesg->bhg;
      ovs->dat.ovl.flipped  = TRUE;

      assert(ovs->dat.ovl.a_hang  == omesg->ahg);
      assert(ovs->dat.ovl.b_hang  == omesg->bhg);
      break;

    case  AS_OUTTIE:
      ovs->dat.ovl.a_hang   = -omesg->bhg;
      ovs->dat.ovl.b_hang   = -omesg->ahg;
      ovs->dat.ovl.flipped  = TRUE;

      assert(ovs->dat.ovl.a_hang  == -omesg->bhg);
      assert(ovs->dat.ovl.b_hang  == -omesg->ahg);
      break;

    case  AS_ANTI:
      ovs->dat.ovl.a_hang   = -omesg->bhg;
      ovs->dat.ovl.b_hang   = -omesg->ahg;
      ovs->dat.ovl.flipped  = FALSE;

      assert(ovs->dat.ovl.a_hang  == -omesg->bhg);
      assert(ovs->dat.ovl.b_hang  == -omesg->ahg);
      break;

    default:
      fprintf(stderr, "YIKES:  Bad overlap orientation = %d for a = %d  b = %d\n",
              omesg->orientation, omesg->aifrag, omesg->bifrag);
      break;
  }
}



int
AS_OVS_convertOVLdumpToOVSoverlap(char *line, OVSoverlap *olap) {
  char *ptrs[16] = {0};
  int   items    = stringSplit(line, ptrs, 16);

  if (items == 7) {
    olap->a_iid              = atoi(ptrs[0]);
    olap->b_iid              = atoi(ptrs[1]);
    olap->dat.ovl.flipped    = (ptrs[2][0] == 'i') || (ptrs[2][0] == 'I');
    olap->dat.ovl.a_hang     = atoi(ptrs[3]);
    olap->dat.ovl.b_hang     = atoi(ptrs[4]);
    olap->dat.ovl.orig_erate = AS_OVS_encodeQuality(atof(ptrs[5]) / 100.0);
    olap->dat.ovl.corr_erate = AS_OVS_encodeQuality(atof(ptrs[6]) / 100.0);
    olap->dat.ovl.type       = AS_OVS_TYPE_OVL;

#if 0
    fprintf(stderr, "%u %u %f\n",
            (uint32)olap->dat.ovl.orig_erate,
            (uint32)AS_OVS_encodeQuality(atof(ptrs[5]) / 100.0),
            atof(ptrs[5]) / 100.0);
#endif

    assert(olap->dat.ovl.a_hang     == atoi(ptrs[3]));
    assert(olap->dat.ovl.b_hang     == atoi(ptrs[4]));
    assert(olap->dat.ovl.orig_erate == AS_OVS_encodeQuality(atof(ptrs[5]) / 100.0));
    assert(olap->dat.ovl.corr_erate == AS_OVS_encodeQuality(atof(ptrs[6]) / 100.0));
  } else {
    //  Should report the line, but we munged it.
  }

  return(items == 7);
}



int
AS_OVS_convertOBTdumpToOVSoverlap(char *line, OVSoverlap *olap) {
  char *ptrs[16] = {0};
  int   items    = stringSplit(line, ptrs, 16);

  if (items == 10) {
    olap->a_iid  = atoi(ptrs[0]);
    olap->b_iid  = atoi(ptrs[1]);

    olap->dat.obt.fwd    = (ptrs[2][0] == 'f');
    olap->dat.obt.a_beg  = atoi(ptrs[3]);
    olap->dat.obt.a_end  = atoi(ptrs[4]);
    olap->dat.obt.b_beg  = atoi(ptrs[5]);
    olap->dat.obt.b_end  = atoi(ptrs[6]);
    olap->dat.obt.erate  = AS_OVS_encodeQuality(atof(ptrs[7]) / 100.0);
    olap->dat.ovl.type   = AS_OVS_TYPE_OBT;

    assert(olap->dat.obt.a_beg == atoi(ptrs[3]));
    assert(olap->dat.obt.a_end == atoi(ptrs[4]));
    assert(olap->dat.obt.b_beg == atoi(ptrs[5]));
    assert(olap->dat.obt.b_end == atoi(ptrs[6]));
  } else {
    //  Should report the line, but we munged it.
  }

  return(items == 10);
}

