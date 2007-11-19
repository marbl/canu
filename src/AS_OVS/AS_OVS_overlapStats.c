
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

//  $Id: AS_OVS_overlapStats.c,v 1.1 2007-11-19 13:18:29 brianwalenz Exp $

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"

#include "overlapStore.h"



void
AS_OVS_histogramAdd(OverlapStoreHistogram *h, uint16 val) {
  h->nSamples++;
  h->histogram[val]++;
}


void
AS_OVS_histogramCompute(OverlapStoreHistogram *h) {
  uint32   i = 0;
  uint64   n = 0;
  uint16   t[65536] = {0};

  h->median  = 0;
  h->mean    = 0.0;
  h->stddev  = 0.0;
  h->mode    = 0;
  h->mad     = 0.0;

  //  median - count up to nSamples/2
  //
  for (i=0, n=0; n < h->nSamples/2; i++)
    n += h->histogram[i];
  h->median = i - 1;

  //  mean - sum, divide by n
  //
  for (i=0; i<65536; i++)
    h->mean += (double)h->histogram[i] * (double)i;
  h->mean /= h->nSamples;

  //  stddev -- we have h[i] values, and the value is i.
  //
  for (i=0; i<65536; i++)
    h->stddev += h->histogram[i] * (h->mean - i) * (h->mean - i);
  h->stddev  = sqrt(h->stddev / h->nSamples);

  //  mode - just find the max
  //
  for (i=0, n=0; i < 65536; i++)
    if (h->histogram[n] < h->histogram[i])
      n = i;
  h->mode = n;

  //  mad - 1.4826*median(abs(median(v)-v))
  //
  for (i=0; i<65536; i++)
    t[(h->median < i) ? (i - h->median) : (h->median - i)] += h->histogram[i];
  for (i=0, n=0; n < h->nSamples/2; i++)
    n += t[i];
  h->mad = 1.4826 * (i - 1);
}


void
AS_OVS_histogramShow(char *label, char *type, OverlapStoreHistogram *h) {
  fprintf(stdout, "\n");
  fprintf(stdout, "%s -- %s\n", type, label);
  fprintf(stdout, "nSamples      "F_U64"\n", h->nSamples);
  fprintf(stdout, "median        %hu\n", h->median);
  fprintf(stdout, "mean/stddev   %f +- %f\n", h->mean, h->stddev);
  fprintf(stdout, "mode/mad      %f +- %f\n", h->mode, h->mad);

  uint16  max = 0;
  uint16  min = 0;
  uint16  rng = 0;
  uint16  xxx = 0;

  for (min=0;     (min < 65535) && (h->histogram[min] == 0); min++)
    ;
  for (max=65535; (max > 0)     && (h->histogram[max] == 0); max--)
    ;

  rng = (max - min) / 10 + 1;

  while (min < max) {
    uint64  sum = 0;
    uint16  lim = min + rng;

    while (min < lim)
      sum += h->histogram[min++];

    fprintf(stdout, "%hu - %hu : "F_U64"\n", lim - rng, lim, sum);
  }
}



static
void
AS_OVS_loadClearLengths(OverlapStore *ovs, OVSoverlap *ovl) {

  //  Load the fragment lengths if we don't already have them.  This
  //  needed to be delayed until now, because the choice of clear
  //  region depends on the type of overlap.

  if (ovs->fragClearLength != NULL)
    return;

  FragStream      *frgStream = openFragStream(ovs->gkp, FRAG_S_INF);
  fragRecord       fr        = {0};
  uint64           maxIID    = getLastElemFragStore(ovs->gkp) + 1;

  ovs->fragClearBegin  = (uint16 *)safe_calloc(maxIID, sizeof(uint16));
  ovs->fragClearEnd    = (uint16 *)safe_calloc(maxIID, sizeof(uint16));
  ovs->fragClearLength = (uint16 *)safe_calloc(maxIID, sizeof(uint16));

  int  typ = AS_READ_CLEAR_OBT;

  switch (ovl->dat.ovl.type) {
    case AS_OVS_TYPE_OVL:
      typ = AS_READ_CLEAR_OBT;
      break;
    case AS_OVS_TYPE_OBT:
      typ = AS_READ_CLEAR_OBTINI;
      break;
    case AS_OVS_TYPE_MER:
      typ = AS_READ_CLEAR_UNTRIM;
      break;
    default:
      fprintf(stderr, "Unknown type %d in overlap.\n", ovl->dat.ovl.type);
      exit(1);
      break;
  }

  while (nextFragStream(frgStream, &fr)) {
    ovs->fragClearBegin [getFragRecordIID(&fr)] = getFragRecordClearRegionBegin(&fr, typ);
    ovs->fragClearEnd   [getFragRecordIID(&fr)] = getFragRecordClearRegionEnd  (&fr, typ);
    ovs->fragClearLength[getFragRecordIID(&fr)] = (getFragRecordClearRegionEnd  (&fr, typ) -
                                                   getFragRecordClearRegionBegin(&fr, typ));
  }

  closeFragStream(frgStream);
}



void
AS_OVS_accumulateStats(OverlapStore *ovs, OVSoverlap *ovl) {

  //  No stats unless we know gkpStore.  We could generate histograms
  //  of erates without knowing fragment length, but it's better to
  //  just force the user to always have a gkpStore.

  assert(ovs->gkp != NULL);

  if (ovs->fragClearLength == NULL)
    AS_OVS_loadClearLengths(ovs, ovl);

  ovs->statsUpdated = 1;

  if (ovl->dat.ovl.type == AS_OVS_TYPE_OVL) {
    ovs->stats.numOVL++;

    int32  ah = ovl->dat.ovl.a_hang;
    int32  bh = ovl->dat.ovl.b_hang;
    int32  tp = 0;

    if (ah == 0)
      tp |= 0x00000004;
    else if (ah < 0)
      tp |= 0x00000008;
    else
      tp |= 0x0000000c;

    if (bh == 0)
      tp |= 0x00000001;
    else if (bh < 0)
      tp |= 0x00000002;
    else
      tp |= 0x00000003;

    //  a_hang   b_hang     type    label
    //  (compared to zero)          (describes A)
    //
    //  any      any        0x00    all overlaps
    //                      %0000
    //
    //  =        =          0x05    degenerate     ----------
    //                      %0101                  ----------
    //
    //  =        <          0x06    5' containee   ----------
    //                      %0110                  ---
    //
    //  =        >          0x07    5' contained   ---
    //                      %0111                  ----------
    //
    //  <        =          0x09    3' contained          ---
    //                      %1001                  ----------
    //
    //  <        <          0x0a    5' dovetal        ----------
    //                      %1010                  ----------
    //
    //  <        >          0x0b    contained         ----
    //                      %1011                  ----------
    //
    //  >        =          0x0d    3' containee   ----------
    //                      %1101                         ---
    //
    //  >        <          0x0e    containee      ----------
    //                      %1110                     ----
    //
    //  >        >          0x0f    3' dovetail    ----------
    //                      %1111                     ----------
    //
    //  unused types -- 0x01, 0x02, 0x03, 0x04, 0x08, 0x0c
    //

    //  Swiped from AS_BOG/AS_BOG_BestOverlapGraph.cc::olapLength()
    int32 length = 0;
    if (ah < 0) {
      if (bh < 0)
        length = ovs->fragClearLength[ovl->a_iid] + bh;
      else
        length = ovs->fragClearLength[ovl->b_iid] + ah - bh;
    } else {
      if (bh < 0)
        length = ovs->fragClearLength[ovl->a_iid] + bh - ah;
      else
        length = ovs->fragClearLength[ovl->a_iid] - ah;
    }


    AS_OVS_histogramAdd(&ovs->stats.orig_erate[0], ovl->dat.ovl.orig_erate);
    AS_OVS_histogramAdd(&ovs->stats.corr_erate[0], ovl->dat.ovl.corr_erate);
    AS_OVS_histogramAdd(&ovs->stats.length[0], length);

    AS_OVS_histogramAdd(&ovs->stats.orig_erate[tp], ovl->dat.ovl.orig_erate);
    AS_OVS_histogramAdd(&ovs->stats.corr_erate[tp], ovl->dat.ovl.corr_erate);
    AS_OVS_histogramAdd(&ovs->stats.length[tp], length);
  }


  if (ovl->dat.ovl.type == AS_OVS_TYPE_OBT) {
    ovs->stats.numOBT++;

    AS_OVS_histogramAdd(&ovs->stats.obtAbeg,    ovl->dat.obt.a_beg);
    AS_OVS_histogramAdd(&ovs->stats.obtAend,    ovl->dat.obt.a_end);
    AS_OVS_histogramAdd(&ovs->stats.obtAlength, ovl->dat.obt.a_end - ovl->dat.obt.a_beg);

    if (ovl->dat.obt.b_beg < ovl->dat.obt.b_end) {
      AS_OVS_histogramAdd(&ovs->stats.obtBbeg,    ovl->dat.obt.b_beg);
      AS_OVS_histogramAdd(&ovs->stats.obtBend,    ovl->dat.obt.b_end);
      AS_OVS_histogramAdd(&ovs->stats.obtBlength, ovl->dat.obt.b_end - ovl->dat.obt.b_beg);
    } else {
      AS_OVS_histogramAdd(&ovs->stats.obtBbeg,    ovl->dat.obt.b_end);
      AS_OVS_histogramAdd(&ovs->stats.obtBend,    ovl->dat.obt.b_beg);
      AS_OVS_histogramAdd(&ovs->stats.obtBlength, ovl->dat.obt.b_beg - ovl->dat.obt.b_end);
    }

    AS_OVS_histogramAdd(&ovs->stats.obtErate,   ovl->dat.obt.erate);
  }


  if (ovl->dat.ovl.type == AS_OVS_TYPE_MER) {
    ovs->stats.numMER++;

    AS_OVS_histogramAdd(&ovs->stats.merApos,   ovl->dat.mer.a_pos);
    AS_OVS_histogramAdd(&ovs->stats.merBpos,   ovl->dat.mer.b_pos);
    AS_OVS_histogramAdd(&ovs->stats.merKcount, ovl->dat.mer.k_count);
  }

  if (ovl->dat.ovl.type == AS_OVS_TYPE_UNS) {
    assert(0);
  }
}
