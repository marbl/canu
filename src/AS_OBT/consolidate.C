
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "util++.H"

extern "C" {
#include "AS_global.h"
#include "AS_OVS_overlapStore.h"
}

#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)


//  sort the position values on the 5' end -- this sorts increasingly
int
position_compare5(const void *a, const void *b) {
  uint32  A = (uint32)*((uint32 *)a);
  uint32  B = (uint32)*((uint32 *)b);

  if (A < B)  return(-1);
  if (A > B)  return(1);
  return(0);
}


//  sort the position values on the 3' end -- this sorts decreasingly
int
position_compare3(const void *a, const void *b) {
  uint32  A = (uint32)*((uint32 *)a);
  uint32  B = (uint32)*((uint32 *)b);

  if (A < B)  return(1);
  if (A > B)  return(-1);
  return(0);
}


void
sortAndOutput(uint32 fid, uint32 numOvl, uint32 *left, uint32 *right) {

  qsort(left,  numOvl, sizeof(uint32), position_compare5);
  qsort(right, numOvl, sizeof(uint32), position_compare3);

  //  XXX:  Print the 5' and 3' stuff
  //
  //  We might as well find the
  //    min/max
  //    min/max with more than one hit
  //    mode
  //  since we have everything here
  //
  //  minN    -- minimum value we've ever seen
  //  minmN   -- minimum value we've seen more than once
  //  minmNc  -- number of times we've seen minm
  //  modeN   -- mode
  //  modeNc  -- number of times we've seen the mode
  //  modeNt  -- temp copy of the mode
  //  modeNtc -- temp copy of the mode, number of times
  //
  uint32  min5, minm5, minm5c,  mode5, mode5c,  mode5t, mode5tc;
  uint32  max3, maxm3, maxm3c,  mode3, mode3c,  mode3t, mode3tc;

  min5 = left[0];
  max3 = right[0];

  minm5 = 9999;       minm5c  = 0;
  maxm3 = 9999;       maxm3c  = 0;

  mode5 = left[0];    mode5c  = 1;
  mode3 = right[0];   mode3c  = 1;

  mode5t = left[0];   mode5tc = 1;
  mode3t = right[0];  mode3tc = 1;

  for (uint32 i=1; i<numOvl; i++) {

    //  5' end.  We scan the list, remembering the best mode
    //  we've seen so far.  When a better one arrives, we copy
    //  it to the saved one -- and keep copying it as it gets
    //  better.
    //
    if (mode5t == left[i]) {  //  Same mode?  Count.
      mode5tc++;
    } else {
      mode5t  = left[i];  //  Different mode, restart.
      mode5tc = 1;
    }
    if (mode5tc > mode5c) {  //  Bigger mode?  Save it.
      mode5  = mode5t;
      mode5c = mode5tc;
    }

    //  If our mode is more than one and we've not seen a multiple hit before
    //  save this position.
    //
    if ((mode5c > 1) && (minm5 == 9999))
      minm5  = mode5;
    if (minm5 == mode5)
      minm5c = mode5c;


    //  Do it all again for the 3' -- remember that we've
    //  sorted this decreasingly.


    if (mode3t == right[i]) {  //  Same mode?  Count.
      mode3tc++;
    } else {
      mode3t  = right[i];  //  Different mode, restart.
      mode3tc = 1;
    }
    if (mode3tc > mode3c) {  //  Bigger mode?  Save it.
      mode3  = mode3t;
      mode3c = mode3tc;
    }
    if ((mode3c > 1) && (maxm3 == 9999))
      maxm3  = mode3;
    if (maxm3 == mode3)
      maxm3c = mode3c;
  }

  //  Output!
  //
  fprintf(stdout, F_U32"  "F_U32" "F_U32" "F_U32" "F_U32" "F_U32"  "F_U32" "F_U32" "F_U32" "F_U32" "F_U32"",
          fid,
          min5, minm5, minm5c, mode5, mode5c,
          max3, maxm3, maxm3c, mode3, mode3c);

  //  Save all the overlaps too
  //
#if 0
  fprintf(stdout, "  "F_U32"", numOvl);
  for (uint32 i=0; i<numOvl; i++)
    fprintf(stdout, "  "F_U32" "F_U32"",
            left[i], right[i]);
#endif

  fprintf(stdout, "\n");
}


int
main(int argc, char **argv) {
  bool           beVerbose = false;
  OverlapStore  *ovs = NULL;
  OVSoverlap     ovl;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;
    } else if (strcmp(argv[arg], "-ovs") == 0) {
      ovs = AS_OVS_openOverlapStore(argv[++arg]);
    } else {
      fprintf(stderr, "%s: unknown arg '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((ovs == NULL) || err)
    fprintf(stderr, "usage: %s [-v] -ovs obtStore > asm.ovl.consolidated\n", argv[0]), exit(1);

  uint32   idAlast     = 0;
  uint32   numOverlaps = 0;
  uint32  *left        = new uint32 [MAX_OVERLAPS_PER_FRAG];
  uint32  *right       = new uint32 [MAX_OVERLAPS_PER_FRAG];

  while (AS_OVS_readOverlapFromStore(ovs, &ovl, AS_OVS_TYPE_OBT)) {
    uint64   wa, wb;

    //  If we see a different idA than we had last time, process
    //  the previous read.
    //
    if ((idAlast != ovl.a_iid) &&
        (numOverlaps > 0)) {
      assert(ovl.a_iid > idAlast);
      sortAndOutput(idAlast, numOverlaps, left, right);
      numOverlaps = 0;
    }

    idAlast            = ovl.a_iid;
    left[numOverlaps]  = ovl.dat.obt.a_beg;
    right[numOverlaps] = ovl.dat.obt.a_end;
    numOverlaps++;
  }

  //  Don't forget to do the last batch!
  if (numOverlaps > 0) {
    sortAndOutput(idAlast, numOverlaps, left, right);
    numOverlaps = 0;
  }

  delete [] left;
  delete [] right;
}
