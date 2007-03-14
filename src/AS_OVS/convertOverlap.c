
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

static char CM_ID[] = "$Id: convertOverlap.c,v 1.10 2007-03-14 19:07:29 brianwalenz Exp $";

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
convertOVLtoBinary(void) {
  GenericMesg  *pmesg;
  OverlapMesg  *omesg;
  OVSoverlap    olap;

  BinaryOverlapFile  *output = AS_OVS_createBinaryOverlapFile(NULL, FALSE);

  while (EOF != ReadProtoMesg_AS(stdin, &pmesg)) {
    if (pmesg->t == MESG_OVL) {
      omesg = pmesg->m;

      //  The asserts below check for encoding errors -- the dat
      //  structure only saves a small number of bits for each field,
      //  and we want to be sure we stored all the bits.

      olap.a_iid = omesg->aifrag;
      olap.b_iid = omesg->bifrag;
      olap.dat.ovl.orig_erate = AS_OVS_encodeQuality(omesg->quality);
      olap.dat.ovl.corr_erate = olap.dat.ovl.orig_erate;
      olap.dat.ovl.type = AS_OVS_TYPE_OVL;

      switch (omesg->orientation) {
        case  AS_NORMAL:
          olap.dat.ovl.a_hang   = omesg->ahg;
          olap.dat.ovl.b_hang   = omesg->bhg;
          olap.dat.ovl.flipped  = FALSE;
          
          assert(olap.dat.ovl.a_hang  == omesg->ahg);
          assert(olap.dat.ovl.b_hang  == omesg->bhg);
          break;

        case  AS_INNIE:
          olap.dat.ovl.a_hang   = omesg->ahg;
          olap.dat.ovl.b_hang   = omesg->bhg;
          olap.dat.ovl.flipped  = TRUE;

          assert(olap.dat.ovl.a_hang  == omesg->ahg);
          assert(olap.dat.ovl.b_hang  == omesg->bhg);
          break;

        case  AS_OUTTIE:
          olap.dat.ovl.a_hang   = -omesg->bhg;
          olap.dat.ovl.b_hang   = -omesg->ahg;
          olap.dat.ovl.flipped  = TRUE;

          assert(olap.dat.ovl.a_hang  == -omesg->bhg);
          assert(olap.dat.ovl.b_hang  == -omesg->ahg);
          break;

        case  AS_ANTI:
          olap.dat.ovl.a_hang   = -omesg->bhg;
          olap.dat.ovl.b_hang   = -omesg->ahg;
          olap.dat.ovl.flipped  = FALSE;

          assert(olap.dat.ovl.a_hang  == -omesg->bhg);
          assert(olap.dat.ovl.b_hang  == -omesg->ahg);
          break;

        default:
          fprintf(stderr, "YIKES:  Bad overlap orientation = %d for a = %d  b = %d\n",
                  omesg->orientation, omesg->aifrag, omesg->bifrag);
          continue;
      }

      AS_OVS_writeOverlap(output, &olap);
    }
  }

  AS_OVS_closeBinaryOverlapFile(output);
}


void   convertOVLDUMPtoBinary(void) {
  char         *ptrs[16];
  char          line[1024];
  OVSoverlap    olap;

  BinaryOverlapFile  *output = AS_OVS_createBinaryOverlapFile(NULL, FALSE);

  //  Comment from GrowOlapStoreOVL:
  //  Open file in  Dump_Format_Input_Path  and read its overlaps
  //  and add them to the store files.  The data in the file is
  //  assumed to have no duplicates and each overlap should appear
  //  only once (i.e., A overlapping B should not also be reported
  //  as B overlapping A.

  fgets(line, 1024, stdin);
  while (!feof(stdin)) {
    int  items = stringSplit(line, ptrs, 16);

    if (items == 7) {
      olap.a_iid              = atoi(ptrs[0]);
      olap.b_iid              = atoi(ptrs[1]);
      olap.dat.ovl.flipped    = (ptrs[2][0] == 'n') || (ptrs[2][0] == 'N');
      olap.dat.ovl.a_hang     = atoi(ptrs[3]);
      olap.dat.ovl.b_hang     = atoi(ptrs[4]);
      olap.dat.ovl.orig_erate = atoi(ptrs[4]);
      olap.dat.ovl.corr_erate = atoi(ptrs[4]);
      olap.dat.ovl.type       = AS_OVS_TYPE_OVL;

      assert(olap.dat.ovl.a_hang     == atoi(ptrs[3]));
      assert(olap.dat.ovl.b_hang     == atoi(ptrs[4]));
      assert(olap.dat.ovl.orig_erate == atoi(ptrs[4]));
      assert(olap.dat.ovl.corr_erate == atoi(ptrs[4]));

      AS_OVS_writeOverlap(output, &olap);
    } else {
      //  Should report the line, but we munged it.
    }

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
    int  items = stringSplit(line, ptrs, 16);

    if (items == 10) {
      olap.a_iid  = atoi(ptrs[0]);
      olap.b_iid  = atoi(ptrs[1]);

      olap.dat.obt.fwd    = (ptrs[2][0] == 'f');
      olap.dat.obt.a_beg  = atoi(ptrs[3]);
      olap.dat.obt.a_end  = atoi(ptrs[4]);
      olap.dat.obt.b_beg  = atoi(ptrs[6]);
      olap.dat.obt.b_end  = atoi(ptrs[7]);
      olap.dat.obt.erate  = AS_OVS_encodeQuality(atof(ptrs[9]) / 100.0);
      olap.dat.ovl.type   = AS_OVS_TYPE_OBT;

      assert(olap.dat.obt.a_beg == atoi(ptrs[3]));
      assert(olap.dat.obt.a_end == atoi(ptrs[4]));
      assert(olap.dat.obt.b_beg == atoi(ptrs[6]));
      assert(olap.dat.obt.b_end == atoi(ptrs[7]));

      AS_OVS_writeOverlap(output, &olap);
    } else {
      //  Should report the line, but we munged it.
    }

    fgets(line, 1024, stdin);
  }

  AS_OVS_closeBinaryOverlapFile(output);
}


void   convertMERtoBinary(void) {
}













void   convertOVLtoASCII(void) {
  OVSoverlap    olap;
  OverlapMesg   omesg;
  GenericMesg   pmesg;
  signed char   deltas[16] = {0};

  pmesg.m = &omesg; 
  pmesg.t = MESG_OVL;

  omesg.aifrag       = 0;
  omesg.bifrag       = 0;
  omesg.ahg          = 0;
  omesg.bhg          = 0;
  omesg.orientation  = 0;
  omesg.overlap_type = 0;
  omesg.quality      = 0;
  omesg.min_offset   = 0;
  omesg.max_offset   = 0;
  omesg.polymorph_ct = 0;
  omesg.delta        = deltas;

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
    omesg.delta        = deltas;

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
    fprintf(stdout, "%7d %7d  %c %4d %4d %4d  %4d %4d %4d  %5.2f\n",
            olap.a_iid, olap.b_iid,
            olap.dat.obt.fwd ? 'f' : 'r',
            olap.dat.obt.a_beg,
            olap.dat.obt.a_end,
            666,
            olap.dat.obt.b_beg,
            olap.dat.obt.b_end,
            666,
            AS_OVS_decodeQuality(olap.dat.obt.erate) * 100.0);
  }

  AS_OVS_closeBinaryOverlapFile(input);
}


void   convertMERtoASCII(void) {
}

