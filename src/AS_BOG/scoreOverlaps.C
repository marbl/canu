
/**************************************************************************
 * Copyright (C) 2005, J Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: scoreOverlaps.C,v 1.1 2009-06-18 14:55:50 brianwalenz Exp $";

//  Scores all overlaps for a given frag IID on the command line

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_UnitigGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"
#include "AS_BOG_MateChecker.hh"

char         *gkpName = NULL;
gkStore      *gkpStore;

char         *ovlName = NULL;
OverlapStore *ovlStore;

uint32        consensusCutoff;
uint32        mismatchCutoff;

#define END_5  '5'
#define END_3  '3'
#define END_C  'C'

uint32 AEnd(const OVSoverlap& olap) {
  if (olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 0)
    return END_5;
  if (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > 0)
    return END_3;
  return END_C;
}

uint32 BEnd(const OVSoverlap& olap) {
  if (olap.dat.ovl.a_hang < 0 && olap.dat.ovl.b_hang < 0)
    return((olap.dat.ovl.flipped) ? END_5 : END_3);
  if (olap.dat.ovl.a_hang > 0 && olap.dat.ovl.b_hang > 0)
    return((olap.dat.ovl.flipped) ? END_3 : END_5);
  return(END_C);
}


uint32
getFragmentLength(uint32 iid) {
  gkFragment fr;

  gkpStore->gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);
  return(fr.gkFragment_getClearRegionLength());
}

//  True if read IID is contained in something.
uint32
isContained(uint32 iid) {
  OverlapStore *ovs = AS_OVS_openOverlapStore(ovlName);
  OVSoverlap    ovl;
  uint32        isc = 0;

  AS_OVS_setRangeOverlapStore(ovs, iid, iid);
  while (AS_OVS_readOverlapFromStore(ovs, &ovl, AS_OVS_TYPE_OVL)) {
    if (ovl.dat.ovl.orig_erate > consensusCutoff)
      continue;
    if (ovl.dat.ovl.corr_erate > mismatchCutoff)
      continue;

    int a_hang = ovl.dat.ovl.a_hang;
    int b_hang = ovl.dat.ovl.b_hang;

    if (((a_hang >= 0) && (b_hang <= 0)) ||
        ((a_hang <= 0) && (b_hang >= 0)))
      isc = 1;
  }

  AS_OVS_closeOverlapStore(ovs);

  return(isc);
}



//  Unfortunately, this is copied from AS_BOG_BestOverlapGraph.hh,
//  with slight modifications since we don't have the fragInfo
//  structure loaded.
//
uint64
scoreOverlap(const OVSoverlap& olap) {

  if (olap.dat.ovl.orig_erate > consensusCutoff)
    return 0;
  if (olap.dat.ovl.corr_erate > mismatchCutoff)
    return 0;

  uint64  leng = 0;
  uint64  corr = (AS_OVS_MAX_ERATE - olap.dat.ovl.corr_erate);
  uint64  orig = (AS_OVS_MAX_ERATE - olap.dat.ovl.orig_erate);

  //  Shift AFTER assigning to a 64-bit value to avoid overflows.
  corr <<= AS_OVS_ERRBITS;

  int a_hang = olap.dat.ovl.a_hang;
  int b_hang = olap.dat.ovl.b_hang;

  //  Containments - the length of the overlaps are all the same.
  //  We return the quality.  We possibly do not need this test, as
  //  the leng is only set for dovetail overlaps, but lets play
  //  safe.
  //
  if (((a_hang >= 0) && (b_hang <= 0)) ||
      ((a_hang <= 0) && (b_hang >= 0)))
    return(corr | orig);

  //  Dovetails - the length of the overlap is the score, but we
  //  bias towards lower error.
  //
  if ((a_hang < 0) && (b_hang < 0))
    leng = getFragmentLength(olap.a_iid) + b_hang;

  if ((a_hang > 0) && (b_hang > 0))
    leng = getFragmentLength(olap.a_iid) - a_hang;

  //  Shift AFTER assigning to a 64-bit value to avoid overflows.
  leng <<= (2 * AS_OVS_ERRBITS);

  return(leng | corr | orig);
};




int
main(int argc, char **argv) {
  OVSoverlap    ovl;
  uint32        iid = 0;
  double        AS_UTG_ERROR_RATE = 0.06;

  argc = AS_configure(argc, argv);

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];
    } else if (strcmp(argv[arg], "-o") == 0) {
      ovlName = argv[++arg];
    } else if (strcmp(argv[arg], "-i") == 0) {
      iid = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      AS_UTG_ERROR_RATE = atof(argv[++arg]);
    } else {
      fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  ovlStore = AS_OVS_openOverlapStore(ovlName);

  mismatchCutoff  = AS_OVS_encodeQuality(AS_UTG_ERROR_RATE);
  consensusCutoff = AS_OVS_encodeQuality(AS_CNS_ERROR_RATE);

  AS_OVS_setRangeOverlapStore(ovlStore, iid, iid);
  while (AS_OVS_readOverlapFromStore(ovlStore, &ovl, AS_OVS_TYPE_OVL)) {
    uint64  score = scoreOverlap(ovl);

    fprintf(stdout, "%10d %c%c %10d %c%c  %c  %5d %5d %.2f %.2f -- len %12llu score %12llu\n",
            ovl.a_iid, isContained(ovl.a_iid) ? 'C' : 'D', AEnd(ovl),
            ovl.b_iid, isContained(ovl.b_iid) ? 'C' : 'D', BEnd(ovl),
            ovl.dat.ovl.flipped ? 'I' : 'N',
            ovl.dat.ovl.a_hang,
            ovl.dat.ovl.b_hang,
            AS_OVS_decodeQuality(ovl.dat.ovl.orig_erate) * 100.0,
            AS_OVS_decodeQuality(ovl.dat.ovl.corr_erate) * 100.0,
            score >> (2 * AS_OVS_ERRBITS),
            score,
            isContained(ovl.b_iid) ? 'C' : 'D');
  }

  exit(0);
}
