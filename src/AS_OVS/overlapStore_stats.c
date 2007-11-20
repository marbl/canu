
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

//  $Id: overlapStore_stats.c,v 1.6 2007-11-20 07:08:28 brianwalenz Exp $

#include "AS_OVS_overlapStore.h"
#include "overlapStore.h"


void
dumpStats(char *storeName) {
  OverlapStore  *ovs = AS_OVS_openOverlapStore(storeName);
  int            i;
  char          *labels[16] = { "all overlaps",
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                "degenerate",
                                "5' containee",
                                "5' contained",
                                NULL,
                                "3' contained",
                                "5' dovetail",
                                "contained",
                                NULL,
                                "3'containee",
                                "containee",
                                "3' dovetail" };

  fprintf(stdout, "numOVL: "F_U64"\n", ovs->stats.numOVL);
  fprintf(stdout, "numOBT: "F_U64"\n", ovs->stats.numOBT);
  fprintf(stdout, "numMER: "F_U64"\n", ovs->stats.numMER);

  if (ovs->stats.numOVL > 0) {
    for (i=0; i<16; i++) {
      if (labels[i]) {
        AS_OVS_histogramShow(labels[i], "orig_erate", &ovs->stats.orig_erate[i]);
        AS_OVS_histogramShow(labels[i], "corr_erate", &ovs->stats.corr_erate[i]);
        AS_OVS_histogramShow(labels[i], "length",     &ovs->stats.length[i]);
      }
    }
  }

  if (ovs->stats.numOBT > 0) {
    AS_OVS_histogramCompute(&ovs->stats.obtAbeg);
    AS_OVS_histogramCompute(&ovs->stats.obtAend);
    AS_OVS_histogramCompute(&ovs->stats.obtAlength);

    AS_OVS_histogramCompute(&ovs->stats.obtBbeg);
    AS_OVS_histogramCompute(&ovs->stats.obtBend);
    AS_OVS_histogramCompute(&ovs->stats.obtBlength);

    AS_OVS_histogramCompute(&ovs->stats.obtErate);


    AS_OVS_histogramShow("A beg",      "", &ovs->stats.obtAbeg);
    AS_OVS_histogramShow("A end",      "", &ovs->stats.obtAend);
    AS_OVS_histogramShow("A length",   "", &ovs->stats.obtAlength);

    AS_OVS_histogramShow("B beg",      "", &ovs->stats.obtBbeg);
    AS_OVS_histogramShow("B end",      "", &ovs->stats.obtBend);
    AS_OVS_histogramShow("B length",   "", &ovs->stats.obtBlength);

    AS_OVS_histogramShow("error rate", "", &ovs->stats.obtErate);
  }

  if (ovs->stats.numMER > 0) {
    AS_OVS_histogramCompute(&ovs->stats.merApos);
    AS_OVS_histogramCompute(&ovs->stats.merBpos);
    AS_OVS_histogramCompute(&ovs->stats.merKcount);

    AS_OVS_histogramShow("A position", "", &ovs->stats.merApos);
    AS_OVS_histogramShow("B position", "", &ovs->stats.merBpos);
    AS_OVS_histogramShow("K count",    "", &ovs->stats.merKcount);
  }

  AS_OVS_closeOverlapStore(ovs);
}



void
rebuildStats(char *storeName, char *gkpName) {
  OverlapStore  *ovs = AS_OVS_openOverlapStore(storeName);
  OVSoverlap     ovl;
  int            i;

  ovs->gkp = openGateKeeperStore(gkpName, FALSE);

  //  Clear the stat struct
  memset(&ovs->stats, 0, sizeof(OverlapStoreStats));

  //  Read all the overlaps, compute stats on them
  while (AS_OVS_readOverlapFromStore(ovs, &ovl, AS_OVS_TYPE_ANY) == TRUE)
    AS_OVS_accumulateStats(ovs, &ovl);

  //  Finish the stats
  for (i=0; i<16; i++) {
    AS_OVS_histogramCompute(&ovs->stats.orig_erate[i]);
    AS_OVS_histogramCompute(&ovs->stats.corr_erate[i]);
    AS_OVS_histogramCompute(&ovs->stats.length[i]);
  }

  //  Close the store, updating the stats
  AS_OVS_closeOverlapStore(ovs);
}
