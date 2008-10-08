
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

const char *mainid = "$Id: filterOverlap.c,v 1.5 2008-10-08 22:02:58 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"

#define  FORMAT_NONE      -1
#define  FORMAT_OVL       AS_READ_CLEAR_OBT
#define  FORMAT_OBT       AS_READ_CLEAR_OBTINI

void   filterOVL(void);
void   filterOBT(void);

uint32  minLength      = 0;
uint32  maxError       = 0;
uint32  noDovetail     = 0;
uint32  noContainment  = 0;
char   *gkpStoreName   = NULL;

typedef struct {
  uint64   len:12;
  uint64   beg:12;
  uint64   end:12;
} fragInfo;

uint32    numReads       = 0;
fragInfo *readLength     = NULL;

int
main(int argc, char **argv) {
  int    toBinary      = 0;
  int    toASCII       = 0;
  int    format        = FORMAT_NONE;

  argc = AS_configure(argc, argv);

  maxError = AS_OVS_encodeQuality(1.0);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-ovl") == 0) {
      format = FORMAT_OVL;
    } else if (strcmp(argv[arg], "-obt") == 0) {
      format = FORMAT_OBT;

    } else if (strncmp(argv[arg], "-minlength", 5) == 0) {
      minLength = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-maxerror", 5) == 0) {
      double e = atof(argv[++arg]);
      maxError = AS_OVS_encodeQuality(e);
    } else if (strncmp(argv[arg], "-nocontainment", 4) == 0) {
      //  aka, only dovetail
      noContainment = 1;
    } else if (strncmp(argv[arg], "-nodovetail", 4) == 0) {
      //  aka, only containment
      noDovetail = 1;

    } else if (strcmp(argv[arg], "-gkp") == 0) {
      gkpStoreName = argv[++arg];

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((format == FORMAT_NONE) ||
      (err)) {
    fprintf(stderr, "usage: %s [-ovl | -obt] < input > output\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Filters overlaps (raw binary files, not the store) based on\n");
    fprintf(stderr, "length, error, dovetail or containment.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -ovl             -- overlaps are OVL\n");
    fprintf(stderr, "  -obt             -- overlaps are OBT\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -minlength l     -- throw out overlaps shorter than l\n");
    fprintf(stderr, "  -maxerror e      -- throw out overlaps with more than fraction e error\n");
    fprintf(stderr, "  -nocontainment   -- throw out containment overlaps\n");
    fprintf(stderr, "  -nodovetail      -- throw out dovetail overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -gkp gkpstore     Needed for -ovl or -nocontainment or -nodovetail\n");
    exit(1);
  }

  if (noDovetail || noContainment || (format == FORMAT_OVL)) {
    GateKeeperStore *gkp = openGateKeeperStore(gkpStoreName, FALSE);
    FragStream      *fs  = openFragStream(gkp, FRAG_S_INF);
    fragRecord       fr;

    numReads   = getLastElemFragStore(gkp);
    readLength = (fragInfo *)safe_malloc(sizeof(fragInfo) * numReads);

    fprintf(stderr, "Reading gkpStore to get clear ranges for "F_U32" reads.\n", numReads);

    while (nextFragStream(fs, &fr)) {
      readLength[getFragRecordIID(&fr)].len = getFragRecordSequenceLength(&fr);
      readLength[getFragRecordIID(&fr)].beg = getFragRecordClearRegionBegin(&fr, format);
      readLength[getFragRecordIID(&fr)].end = getFragRecordClearRegionEnd  (&fr, format);
    }

    closeFragStream(fs);
    closeGateKeeperStore(gkp);
  }

  fprintf(stderr, "WARNING:\n");
  fprintf(stderr, "WARNING:  This has not been fully tested.  Only -obt -minlength\n");
  fprintf(stderr, "WARNING:  is guaranteed (unless it doesn't work).  Though,\n");
  fprintf(stderr, "WARNING:  -maxerror is pretty trivial.\n");
  fprintf(stderr, "WARNING:\n");

  switch (format) {
    case FORMAT_OVL:
      filterOVL();
      break;
    case FORMAT_OBT:
      filterOBT();
      break;
    default:
      fprintf(stderr, "%s: unknown format (%d)?!?!\n", argv[0], format);
      return(1);
      break;
  }

  return(0);
}


void
filterOVL(void) {
  OVSoverlap    olap;
  uint64        discard = 0;
  uint64        total   = 0;
  OverlapMesg   omesg;
  GenericMesg   pmesg;

  omesg.alignment_trace = NULL;

  pmesg.m = &omesg;
  pmesg.t = MESG_OVL;

  BinaryOverlapFile  *input  = AS_OVS_openBinaryOverlapFile(NULL, FALSE);
  BinaryOverlapFile  *output = AS_OVS_createBinaryOverlapFile(NULL, FALSE);

  while (AS_OVS_readOverlap(input, &olap)) {
    int   doOutput = 1;

    int32   ahng = olap.dat.ovl.a_hang;
    int32   bhng = olap.dat.ovl.b_hang;


    if (((ahng > 0) && (bhng > 0)) ||
        ((ahng < 0) && (bhng < 0))) {
      //  is Dovetail
      if (noDovetail)
        doOutput = 0;
    } else {
      //  not Dovetail
      if (noContainment)
        doOutput = 0;
    }

    if (olap.dat.ovl.corr_erate > maxError)
      doOutput = 0;

    if (ahng < 0)
      ahng = -ahng;
    if (bhng < 0)
      bhng = -bhng;

    if (((readLength[olap.a_iid].end - readLength[olap.a_iid].beg - ahng) +
         (readLength[olap.b_iid].end - readLength[olap.b_iid].beg - bhng)) < minLength + minLength)
      doOutput = 0;

    if (doOutput)
      AS_OVS_writeOverlap(output, &olap);
    else
      discard++;

    total++;

    if ((total % 1000000) == 0) {
      fprintf(stderr, F_U64"/"F_U64" = %f%% trash.\n", discard, total, 100.0 * discard / total);
    }
  }

  AS_OVS_closeBinaryOverlapFile(input);
  AS_OVS_closeBinaryOverlapFile(output);
}


void
filterOBT(void) {
  OVSoverlap    olap;
  uint64        discard = 0;
  uint64        total   = 0;

  BinaryOverlapFile  *input  = AS_OVS_openBinaryOverlapFile(NULL, FALSE);
  BinaryOverlapFile  *output = AS_OVS_createBinaryOverlapFile(NULL, FALSE);

  while (AS_OVS_readOverlap(input, &olap)) {
    int   doOutput = 1;

    uint32   abeg       = olap.dat.obt.a_beg;
    uint32   aend       = olap.dat.obt.a_end;
    uint32   bbeg       = olap.dat.obt.b_beg;
    uint32   bend       = olap.dat.obt.b_end;

    if (olap.dat.obt.fwd == 0) {
      bbeg = olap.dat.obt.b_end;
      bend = olap.dat.obt.b_beg;
    }

#if 0
        fprintf(stdout, "a: %d %d %d -%c- b: %d %d %d\n",
                (int)abeg, (int)aend, (int)readLength[olap.a_iid].end,
                olap.dat.obt.fwd ? 'f' : 'r',
                (int)bbeg, (int)bend, (int)readLength[olap.b_iid].end);
#endif

    //  -------------
    //          |   |
    //          -------------
    //
    if (readLength) {
      if (((bbeg < 10) && (readLength[olap.a_iid].end < (aend + 10))) ||
          ((abeg < 10) && (readLength[olap.b_iid].end < (bend + 10)))) {
        //  is Dovetail
        if (noDovetail)
          doOutput = 0;
      } else {
        //  not Dovetail
        if (noContainment)
          doOutput = 0;
      }
    }

    if (olap.dat.obt.erate > maxError)
      doOutput = 0;

    if ((aend - abeg < minLength) ||
        (bend - bbeg < minLength))
      doOutput = 0;

    if (doOutput)
      AS_OVS_writeOverlap(output, &olap);
    else
      discard++;

    total++;

    if ((total % 1000000) == 0) {
      fprintf(stderr, F_U64"/"F_U64" = %f%% trash.\n", discard, total, 100.0 * discard / total);
    }
  }

  AS_OVS_closeBinaryOverlapFile(input);
  AS_OVS_closeBinaryOverlapFile(output);
}
