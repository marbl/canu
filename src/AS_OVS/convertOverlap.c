
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

static char CM_ID[] = "$Id: convertOverlap.c,v 1.1 2007-02-28 14:07:33 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AS_global.h"
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
    fprintf(stderr, "usage: %s .....\n", argv[0]);
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
  OverlapMesg  *omesg;

  OVSoverlap    fwdolap;
  OVSoverlap    revolap;

  BinaryOverlapFile  *output = AS_OVS_createBinaryOverlapFile(NULL, TRUE);

  while (EOF != ReadProtoMesg_AS(stdin, &pmesg)) {
    if (pmesg->t == MESG_OVL) {
      omesg = pmesg->m;

      fwdolap.aid = omesg->aifrag;
      fwdolap.bid = omesg->bifrag;
      fwdolap.dat.ovl.origE = fwdolap.dat.ovl.corrE = Shrink_Quality(omesg->quality);

      revolap.aid = omesg->bifrag;
      revolap.bid = omesg->aifrag;
      revolap.dat.ovl.origE = revolap.dat.ovl.corrE = Shrink_Quality(omesg->quality);

      switch (omesg->orientation) {
        case  AS_NORMAL:
          fwdolap.dat.ovl.Ahang   = omesg->ahg;
          fwdolap.dat.ovl.Bhang   = omesg->bhg;
          fwdolap.dat.ovl.flipped  = FALSE;

          revolap.dat.ovl.Ahang   = -omesg->ahg;
          revolap.dat.ovl.Bhang   = -omesg->bhg;
          revolap.dat.ovl.flipped  = FALSE;
          break;

        case  AS_INNIE:
          fwdolap.dat.ovl.Ahang   = omesg->ahg;
          fwdolap.dat.ovl.Bhang   = omesg->bhg;
          fwdolap.dat.ovl.flipped  = TRUE;

          revolap.dat.ovl.Ahang   = omesg->bhg;
          revolap.dat.ovl.Bhang   = omesg->ahg;
          revolap.dat.ovl.flipped  = TRUE;
          break;

        case  AS_OUTTIE:
          fwdolap.dat.ovl.Ahang   = -omesg->ahg;
          fwdolap.dat.ovl.Bhang   = -omesg->bhg;
          fwdolap.dat.ovl.flipped  = TRUE;

          revolap.dat.ovl.Ahang   = -omesg->bhg;
          revolap.dat.ovl.Bhang   = -omesg->ahg;
          revolap.dat.ovl.flipped  = TRUE;
          break;

        case  AS_ANTI:
          fwdolap.dat.ovl.Ahang   = -omesg->bhg;
          fwdolap.dat.ovl.Bhang   = -omesg->ahg;
          fwdolap.dat.ovl.flipped  = FALSE;

          revolap.dat.ovl.Ahang   = omesg->bhg;
          revolap.dat.ovl.Bhang   = omesg->ahg;
          revolap.dat.ovl.flipped  = FALSE;
          break;

        default:
          fprintf(stderr, "YIKES:  Bad overlap orientation = %d for a = %d  b = %d\n",
                  omesg->orientation, omesg->aifrag, omesg->bifrag);
          continue;
      }

      AS_OVS_writeOverlap(output, &fwdolap);
      AS_OVS_writeOverlap(output, &revolap);
    }
  }

  AS_OVS_closeBinaryOverlapFile(output);
}
  

void   convertOVLDUMPtoBinary(void) {
}


void   convertOBTtoBinary(void) {
}


void   convertMERtoBinary(void) {
}







void   convertOVLtoASCII(void) {
}


void   convertOVLDUMPtoASCII(void) {
}


void   convertOBTtoASCII(void) {
}


void   convertMERtoASCII(void) {
}

