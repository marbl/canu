
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

//  $Id: convertOverlap.c,v 1.17 2008-10-07 15:07:49 brianwalenz Exp $

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"

//  Converts overlaps to/from binary/ascii.
//
//  Reads several input ascii formats:
//    OVL messages
//    OVL dump
//    OBT (partial overlaps from overlapper)
//    MER (overlap seeds)
//
//  Writes binary messages specific to each input format.  You cannot
//  "convert" an OVL message into an OBT message.


#define  FORMAT_NONE      0
#define  FORMAT_OVL       1
#define  FORMAT_OVLDUMP   2
#define  FORMAT_OBT       3
#define  FORMAT_MER       4


void   convertOVLtoBinary(void);
void   convertOVLDUMPtoBinary(void);
void   convertOBTtoBinary(void);
void   convertMERtoBinary(void);


void   convertOVLtoASCII(void);
void   convertOVLDUMPtoASCII(void);
void   convertOBTtoASCII(void);
void   convertMERtoASCII(void);


int
main(int argc, char **argv) {

  int    toBinary     = 0;
  int    toASCII      = 0;
  int    format       = FORMAT_NONE;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      toASCII  = 1;
      toBinary = 0;
    } else if (strcmp(argv[arg], "-b") == 0) {
      toASCII  = 0;
      toBinary = 1;
    } else if (strcmp(argv[arg], "-ovl") == 0) {
      format = FORMAT_OVL;
    } else if (strcmp(argv[arg], "-ovldump") == 0) {
      format = FORMAT_OVLDUMP;
    } else if (strcmp(argv[arg], "-obt") == 0) {
      format = FORMAT_OBT;
    } else if (strcmp(argv[arg], "-mer") == 0) {
      format = FORMAT_MER;
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((format == FORMAT_NONE) ||
      (err)) {
    fprintf(stderr, "usage: %s [-a | -b] [-ovl | -obt | -mer] < input > output\n", argv[0]);
    fprintf(stderr, "  -a    convert to ASCII, from BINARY.\n");
    fprintf(stderr, "  -b    convert to BINARY, from ASCII.\n");
    fprintf(stderr, "  -ovl\n");
    fprintf(stderr, "  -obt\n");
    fprintf(stderr, "  -mer\n");
    exit(1);
  }

  switch (format) {
    case FORMAT_OVL:
      if (toBinary)
        convertOVLtoBinary();
      if (toASCII)
        convertOVLtoASCII();
      break;
    case FORMAT_OVLDUMP:
      if (toBinary)
        convertOVLDUMPtoBinary();
      if (toASCII)
        convertOVLDUMPtoASCII();
      break;
    case FORMAT_OBT:
      if (toBinary)
        convertOBTtoBinary();
      if (toASCII)
        convertOBTtoASCII();
      break;
    case FORMAT_MER:
      if (toBinary)
        convertMERtoBinary();
      if (toASCII)
        convertMERtoASCII();
      break;
    default:
      fprintf(stderr, "%s: unknown format (%d)?!?!\n", argv[0], format);
      return(1);
      break;
  }

  return(0);
}





void
convertOVLtoBinary(void) {
  GenericMesg  *pmesg;
  OVSoverlap    olap;

  BinaryOverlapFile  *output = AS_OVS_createBinaryOverlapFile(NULL, FALSE);

  while (EOF != ReadProtoMesg_AS(stdin, &pmesg)) {
    if (pmesg->t == MESG_OVL) {
      AS_OVS_convertOverlapMesgToOVSoverlap((OverlapMesg *)pmesg->m, &olap);
      AS_OVS_writeOverlap(output, &olap);
    }
  }

  AS_OVS_closeBinaryOverlapFile(output);
}


void   convertOVLDUMPtoBinary(void) {
  char          line[1024];
  OVSoverlap    olap;

  BinaryOverlapFile  *output = AS_OVS_createBinaryOverlapFile(NULL, FALSE);

  fgets(line, 1024, stdin);
  while (!feof(stdin)) {
    if (AS_OVS_convertOVLdumpToOVSoverlap(line, &olap))
      AS_OVS_writeOverlap(output, &olap);
    fgets(line, 1024, stdin);
  }

  AS_OVS_closeBinaryOverlapFile(output);
}


void   convertOBTtoBinary(void) {
  char         *ptrs[16];
  char          line[1024];
  OVSoverlap    olap;

  BinaryOverlapFile  *output = AS_OVS_createBinaryOverlapFile(NULL, FALSE);

  fgets(line, 1024, stdin);
  while (!feof(stdin)) {
    if (AS_OVS_convertOBTdumpToOVSoverlap(line, &olap))
      AS_OVS_writeOverlap(output, &olap);
    fgets(line, 1024, stdin);
  }

  AS_OVS_closeBinaryOverlapFile(output);
}


void   convertMERtoBinary(void) {
  fprintf(stderr, "not implemented.\n");
  exit(1);
}







void   convertOVLtoASCII(void) {
  OVSoverlap    olap;
  OverlapMesg   omesg;
  GenericMesg   pmesg;

  omesg.alignment_trace = NULL;

  pmesg.m = &omesg;
  pmesg.t = MESG_OVL;

  BinaryOverlapFile  *input = AS_OVS_openBinaryOverlapFile(NULL, FALSE);

  while (AS_OVS_readOverlap(input, &olap)) {
    //  orient:   AS_NORMAL, AS_OUTTIE, AS_INNIE, AS_ANTI
    //  otype:    AS_CONTAINMENT, AS_DOVETAIL
    omesg.aifrag       = olap.a_iid;
    omesg.bifrag       = olap.b_iid;
    omesg.ahg          = olap.dat.ovl.a_hang;
    omesg.bhg          = olap.dat.ovl.b_hang;
    omesg.orientation  = (olap.dat.ovl.flipped) ? AS_INNIE : AS_NORMAL;
    omesg.overlap_type = AS_DOVETAIL;
    omesg.quality      = AS_OVS_decodeQuality(olap.dat.ovl.orig_erate);
    omesg.min_offset   = 0;
    omesg.max_offset   = 0;
    omesg.polymorph_ct = 0;
    WriteProtoMesg_AS(stdout, &pmesg);
  }

  AS_OVS_closeBinaryOverlapFile(input);
}


void   convertOVLDUMPtoASCII(void) {
  OVSoverlap    olap;

  BinaryOverlapFile  *input = AS_OVS_openBinaryOverlapFile(NULL, FALSE);

  while (AS_OVS_readOverlap(input, &olap)) {
    fprintf(stdout, "    %8d %8d %c %5d %5d %4.1f %4.1f\n",
            olap.a_iid,
            olap.b_iid,
            olap.dat.ovl.flipped ? 'I' : 'N',
            olap.dat.ovl.a_hang,
            olap.dat.ovl.b_hang,
            AS_OVS_decodeQuality(olap.dat.ovl.orig_erate) * 100.0,
            AS_OVS_decodeQuality(olap.dat.ovl.corr_erate) * 100.0);
  }

  AS_OVS_closeBinaryOverlapFile(input);
}


void   convertOBTtoASCII(void) {
  OVSoverlap    olap;

  BinaryOverlapFile  *input = AS_OVS_openBinaryOverlapFile(NULL, FALSE);

  while (AS_OVS_readOverlap(input, &olap)) {
    fprintf(stdout, "%7d %7d  %c %4d %4d  %4d %4d  %5.2f\n",
            olap.a_iid, olap.b_iid,
            olap.dat.obt.fwd ? 'f' : 'r',
            olap.dat.obt.a_beg,
            olap.dat.obt.a_end,
            olap.dat.obt.b_beg,
            olap.dat.obt.b_end,
            AS_OVS_decodeQuality(olap.dat.obt.erate) * 100.0);
  }

  AS_OVS_closeBinaryOverlapFile(input);
}


void   convertMERtoASCII(void) {
  OVSoverlap    olap;

  BinaryOverlapFile  *input = AS_OVS_openBinaryOverlapFile(NULL, FALSE);

  while (AS_OVS_readOverlap(input, &olap)) {
    fprintf(stdout, "%7d %7d  %c %c  %d  %4d %4d  %4d %4d\n",
            olap.a_iid, olap.b_iid,
            olap.dat.mer.fwd ? 'f' : 'r',
            olap.dat.mer.palindrome ? 'p' : 'n',
            olap.dat.mer.compression_length,
            olap.dat.mer.a_pos,
            olap.dat.mer.b_pos,
            olap.dat.mer.k_count,
            olap.dat.mer.k_len);
  }

  AS_OVS_closeBinaryOverlapFile(input);
}

