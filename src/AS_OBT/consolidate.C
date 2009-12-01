
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

const char *mainid = "$Id: consolidate.C,v 1.20 2009-12-01 22:27:17 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "util++.H"

#include "AS_global.h"
#include "AS_OVS_overlapStore.h"

#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)


typedef struct {
  uint32 iid;     //  A iid for these overlaps
  uint32 min5;    //  minimum value we've ever seen
  uint32 minm5;   //  minimum value we've ever seen more than once
  uint32 minm5c;  //  number of times we've seen minm5
  uint32 mode5;   //  mode
  uint32 mode5c;  //  number of time we've seen the mode
  uint32 max3;
  uint32 maxm3;
  uint32 maxm3c;
  uint32 mode3;
  uint32 mode3c;
} obtConsolidate_t;



//  sort the position values on the 5' end -- this sorts increasingly
int
position_compare5(const void *a, const void *b) {
  OVSoverlap  *A = (OVSoverlap *)a;
  OVSoverlap  *B = (OVSoverlap *)b;

  if (A->dat.obt.a_beg < B->dat.obt.a_beg)  return(-1);
  if (A->dat.obt.a_beg > B->dat.obt.a_beg)  return(1);
  return(0);
}


//  sort the position values on the 3' end -- this sorts decreasingly
int
position_compare3(const void *a, const void *b) {
  OVSoverlap  *A = (OVSoverlap *)a;
  OVSoverlap  *B = (OVSoverlap *)b;

  if (A->dat.obt.a_end < B->dat.obt.a_end)  return(1);
  if (A->dat.obt.a_end > B->dat.obt.a_end)  return(-1);
  return(0);
}


void
consolidate(OverlapStore *ovsprimary, OverlapStore *ovssecondary, obtConsolidate_t &cd) {
}


void
consolidate(OVSoverlap *ovl, int32 ovlLen, obtConsolidate_t &cd) {

  cd.iid    = 0xffffffff;
  cd.min5   = 0xffffffff;
  cd.minm5  = 0xffffffff;
  cd.minm5c = 0;
  cd.mode5  = 0xffffffff;
  cd.mode5c = 1;
  cd.max3   = 0xffffffff;
  cd.maxm3  = 0xffffffff;
  cd.maxm3c = 0;
  cd.mode3  = 0xffffffff;
  cd.mode3c = 1;

  if ((ovl == NULL) || (ovlLen == 0))
    return;

  cd.iid    = ovl[0].a_iid;

  int32  mtmp;  //  Mode we are computing, value
  int32  mcnt;  //  Mode we are computing, number of times we've seen it

  qsort(ovl, ovlLen, sizeof(OVSoverlap), position_compare5);

  cd.min5   = ovl[0].dat.obt.a_beg;
  cd.mode5  = mtmp = ovl[0].dat.obt.a_beg;
  cd.mode5c = mcnt = 1;

  for (uint32 i=1; i<ovlLen; i++) {

    //  5' end.  Scan the list, remembering the best mode we've seen so far.  When a better one
    //  arrives, we copy it to the saved one -- and keep copying it as it gets better.

    if (mtmp == ovl[i].dat.obt.a_beg) {  //  Same mode?  Count.
      mcnt++;
    } else {
      mtmp = ovl[i].dat.obt.a_beg;  //  Different mode, restart.
      mcnt = 1;
    }

    if (mcnt > cd.mode5c) {  //  Bigger mode?  Save it.
      cd.mode5  = mtmp;
      cd.mode5c = mcnt;
    }

    //  If our mode is more than one and we've not seen a multiple hit before
    //  save this position.
    //
    if ((cd.mode5c > 1) && (cd.minm5 == 0xffffffff))
      cd.minm5  = cd.mode5;

    if (cd.minm5 == cd.mode5)
      cd.minm5c = cd.mode5c;
  }

  //  Do it all again for the 3' -- remember that we've sorted this decreasingly.
  //
  qsort(ovl, ovlLen, sizeof(OVSoverlap), position_compare3);

  cd.max3   = ovl[0].dat.obt.a_end;
  cd.mode3  = mtmp = ovl[0].dat.obt.a_end;
  cd.mode3c = mcnt = 1;

  for (uint32 i=1; i<ovlLen; i++) {
    if (mtmp == ovl[i].dat.obt.a_end) {
      mcnt++;
    } else {
      mtmp = ovl[i].dat.obt.a_end;
      mcnt = 1;
    }

    if (mcnt > cd.mode3c) {
      cd.mode3  = mtmp;
      cd.mode3c = mcnt;
    }

    if ((cd.mode3c > 1) && (cd.maxm3 == 0xffffffff))
      cd.maxm3  = cd.mode3;

    if (cd.maxm3 == cd.mode3)
      cd.maxm3c = cd.mode3c;
  }
}






int
main(int argc, char **argv) {
  OverlapStore  *ovsprimary   = 0L;
  OverlapStore  *ovssecondary = 0L;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-ovs") == 0) {
      if (ovsprimary == NULL)
        ovsprimary = AS_OVS_openOverlapStore(argv[++arg]);
      else if (ovssecondary == NULL)
        ovssecondary = AS_OVS_openOverlapStore(argv[++arg]);
      else {
        fprintf(stderr, "Only two obtStores allowed.\n");
        err++;
      }
    } else {
      fprintf(stderr, "%s: unknown arg '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((ovsprimary == NULL) || err)
    fprintf(stderr, "usage: %s -ovs obtStore > asm.ovl.consolidated\n", argv[0]), exit(1);


  uint32      ovlMax = MAX_OVERLAPS_PER_FRAG;
  uint32      ovlLen = 0;
  OVSoverlap *ovl    = (OVSoverlap *)safe_malloc(sizeof(OVSoverlap) * ovlMax);

  ovlLen += AS_OVS_readOverlapsFromStore(ovsprimary,   ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
  ovlLen += AS_OVS_readOverlapsFromStore(ovssecondary, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);

  //  NOTE!  We DO get multiple overlaps for the same pair of fragments in the partial overlap
  //  output.  We used to pick one of the overlaps (the first seen) and ignore the rest.  We do not
  //  do that anymore.

  while (ovlLen > 0) {
    obtConsolidate_t  cd;

    consolidate(ovl, ovlLen, cd);

    fprintf(stdout, F_U32"  "F_U32" "F_U32" "F_U32" "F_U32" "F_U32"  "F_U32" "F_U32" "F_U32" "F_U32" "F_U32"\n",
            cd.iid,
            cd.min5, cd.minm5, cd.minm5c, cd.mode5, cd.mode5c,
            cd.max3, cd.maxm3, cd.maxm3c, cd.mode3, cd.mode3c);

    ovlLen  = 0;
    ovlLen += AS_OVS_readOverlapsFromStore(ovsprimary,   ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
    ovlLen += AS_OVS_readOverlapsFromStore(ovssecondary, ovl + ovlLen, ovlMax - ovlLen, AS_OVS_TYPE_ANY);
  }

  exit(0);
}
