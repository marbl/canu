
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

const char *mainid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "AS_global.H"
#include "AS_MSG_pmesg.H"
#include "AS_OVS_overlap.H"
#include "AS_OVS_overlapFile.H"

#define  FORMAT_NONE      0
#define  FORMAT_OVL       2
#define  FORMAT_OBT       3
#define  FORMAT_MER       4



void
toBINARY(char *out, int format) {
  char         *ptrs[16];
  char          line[1024];
  OVSoverlap    olap;

  BinaryOverlapFile  *output = AS_OVS_createBinaryOverlapFile(out, FALSE);

  fgets(line, 1024, stdin);
  while (!feof(stdin)) {
    switch (format) {
      case FORMAT_OVL:
        if (AS_OVS_convertOVLdumpToOVSoverlap(line, &olap))
          AS_OVS_writeOverlap(output, &olap);
        break;
      case FORMAT_OBT:
        if (AS_OVS_convertOBTdumpToOVSoverlap(line, &olap))
          AS_OVS_writeOverlap(output, &olap);
        break;
      case FORMAT_MER:
        break;
    }

    fgets(line, 1024, stdin);
  }


  AS_OVS_closeBinaryOverlapFile(output);
}



void
toASCII(char *out) {
  BinaryOverlapFile  *input = AS_OVS_openBinaryOverlapFile(out, FALSE);
  OVSoverlap          olap;
  char                olapstring[256];

  while (AS_OVS_readOverlap(input, &olap))
    fprintf(stdout, "%s\n", AS_OVS_toString(olapstring, olap));

  AS_OVS_closeBinaryOverlapFile(input);
}



int
main(int argc, char **argv) {
  int    toB = 0;
  int    toA = 0;
  int    fmt = FORMAT_NONE;
  char  *out = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-a") == 0) {
      toA++;

    } else if (strcmp(argv[arg], "-ovl") == 0) {
      toB++;
      fmt = FORMAT_OVL;

    } else if (strcmp(argv[arg], "-obt") == 0) {
      toB++;
      fmt = FORMAT_OBT;

    } else if (strcmp(argv[arg], "-mer") == 0) {
      toB++;
      fmt = FORMAT_MER;

    } else if (strcmp(argv[arg], "-out") == 0) {
      out = argv[++arg];

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }

  if ((err) ||
      (toA + toB != 1) ||
      out == NULL) {
    fprintf(stderr, "usage: %s [-a | -ovl | -obt | -mer] < input -out output\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "MANDATORY:  specify what to convert\n");
    fprintf(stderr, "  -a           convert to ASCII, from a BINARY overlap file.\n");
    fprintf(stderr, "  -ovl         convert to BINARY, from an ASCII overlap file.\n");
    fprintf(stderr, "  -obt         convert to BINARY, from an ASCII partial overlap file.\n");
    fprintf(stderr, "  -mer         convert to BINARY, from an ASCII mer overlap file.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "ASCII formats are:\n");
    fprintf(stderr, "  OVL:   aIID bIID [I|N] aHang bHang error error_corrected\n");
    fprintf(stderr, "  OBT:   aIID bIID [f|r] aBgn aEnd bBgn bEnd error\n");
    fprintf(stderr, "  MER:   aIID bIID [p|f|r] compression_length aPos bPos kCount kLen\n");
    fprintf(stderr, "\n");
    if (toA + toB == 0)
      fprintf(stderr, "ERROR:  what to do?  Supply exactly one of -a, -ovl, -obt and -mer.\n");
    if (toA + toB > 1)
      fprintf(stderr, "ERROR:  conflicting options.  Supply exactly one of -a, -ovl, -obt and -mer.\n");
    if (out == NULL)
      fprintf(stderr, "ERROR:  no output file specified. Supply one with -out <filename>\n");
    exit(1);
  }

  if (toA)
    toASCII(out);
  else
    toBINARY(out, fmt);

  return(0);
}

