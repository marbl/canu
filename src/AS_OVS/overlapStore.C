
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

#include "overlapStore.H"
#include "AS_OVS_overlap.H"   //  Just to know the sizes of structs
#include "AS_PER_gkpStore.H"  //  Just to know clear region labels

int
main(int argc, char **argv) {
  uint32          operation   = OP_NONE;
  char           *storeName   = NULL;
  char           *gkpName     = NULL;
  uint32          clearRegion = AS_READ_CLEAR_ERROR;

  uint32          dumpBinary  = FALSE;
  double          dumpERate   = 100.0;
  uint32          dumpLength  = 0;
  uint32          dumpType    = 0;

  char           *erateFile   = NULL;

  uint32          bgnIID      = 0;
  uint32          endIID      = UINT32_MAX;
  uint32          qryIID      = 0;

  bool            beVerbose   = false;

  //  Genome size parameters
  uint32          gs_ovlLimit = 100;
  uint32          gs_winSize  = 100;
  uint32          gs_minOvl   = 40;
  uint32          gs_into     = 5;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {

    if        (strcmp(argv[arg], "-c") == 0) {
      fprintf(stderr, "ERROR: building (-c) option no longer supported here.  Use overlapStoreBuild or the parallel version instead.\n"), err++;

    } else if (strcmp(argv[arg], "-m") == 0) {
      fprintf(stderr, "ERROR: merging (-m) option no longer supported here.  There is currently no replacement.\n"), err++;

    } else if (strcmp(argv[arg], "-d") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_DUMP;

    } else if (strcmp(argv[arg], "-p") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      qryIID      = atoi(argv[++arg]);
      storeName   = argv[++arg];
      gkpName     = argv[++arg];
      clearRegion = gkStore_decodeClearRegionLabel(argv[++arg]);
      operation   = OP_DUMP_PICTURE;

    } else if (strcmp(argv[arg], "-G") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      storeName    = argv[++arg];
      gkpName      = argv[++arg];
      gs_ovlLimit  = atoi(argv[++arg]);
      operation    = OP_GENOME_LENGTH;

    } else if (strcmp(argv[arg], "-u") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S, -G or -u may be supplied.\n"), err++;
      storeName   = argv[++arg];
      operation   = OP_UPDATE_ERATES;

    } else if (strcmp(argv[arg], "-E") == 0) {
      dumpERate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-L") == 0) {
      dumpLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-d5") == 0) {
      dumpType |= DUMP_5p;

    } else if (strcmp(argv[arg], "-d3") == 0) {
      dumpType |= DUMP_3p;

    } else if (strcmp(argv[arg], "-dC") == 0) {
      dumpType |= DUMP_CONTAINS;

    } else if (strcmp(argv[arg], "-dc") == 0) {
      dumpType |= DUMP_CONTAINED;

    } else if (strcmp(argv[arg], "-B") == 0) {
      dumpBinary = TRUE;

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endIID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;

    } else if (strcmp(argv[arg], "-q") == 0) {
      if (storeName)
        fprintf(stderr, "ERROR: only one of -c, -m, -d, -q, -s, -S or -u may be supplied.\n"), err++;
      bgnIID    = atoi(argv[++arg]);
      endIID    = bgnIID;
      qryIID    = atoi(argv[++arg]);
      storeName = argv[++arg];
      operation = OP_DUMP;

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if ((argv[arg][0] != '-') && (argv[arg][1] != 0) && (erateFile == NULL)) {
      erateFile = argv[arg];

    } else {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((operation == OP_NONE) || (storeName == NULL) || (err)) {
    fprintf(stderr, "usage: %s -d storeName [-B] [-E erate] [-b beginIID] [-e endIID]\n", argv[0]);
    fprintf(stderr, "       %s -q aiid biid storeName\n", argv[0]);
    fprintf(stderr, "       %s -p iid storeName gkpStore clr\n", argv[0]);
    fprintf(stderr, "       %s -G storeName gkpStore ovlLimit\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "There are six modes of operation, selected by the first option:\n");
    fprintf(stderr, "  -d  dump a store\n");
    fprintf(stderr, "  -q  report the a,b overlap, if it exists.\n");
    fprintf(stderr, "  -p  dump a picture of overlaps to fragment 'iid', using clear region 'clr'.\n");
    fprintf(stderr, "  -G  estimate genome length\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DUMPING - report overlaps in the store\n");
    fprintf(stderr, "  -B                Dump the store as binary, suitable for input to create a new store.\n");
    fprintf(stderr, "  -E erate          Dump only overlaps <= erate error.\n");
    fprintf(stderr, "  -L length         Dump only overlaps that are larger than L bases (only for -p picture mode).\n");
    fprintf(stderr, "  -d5               Dump only overlaps off the 5' end of the A frag.\n");
    fprintf(stderr, "  -d3               Dump only overlaps off the 3' end of the A frag.\n");
    fprintf(stderr, "  -dC               Dump only overlaps that are contained in the A frag (B contained in A).\n");
    fprintf(stderr, "  -dc               Dump only overlaps that are containing the A frag (A contained in B).\n");
    fprintf(stderr, "  -b beginIID       Start dumping at 'beginIID'.\n");
    fprintf(stderr, "  -e endIID         Stop dumping after 'endIID'.\n");
    fprintf(stderr, "  -v                Report statistics (to stderr) on some dumps (-d).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "QUERYING - quickly ask if an overlap exists\n");
    fprintf(stderr, "  -q aiid biid storeName\n");
    fprintf(stderr, "                    If an overlap between fragments 'aiid' and 'biid' exists, it is printed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DUMPING PICTURES - draw a multi-alignment-like picture for a single fragment and its overlaps\n");
    fprintf(stderr, "  -p iid storeName gkpStore clr\n");
    fprintf(stderr, "                    clr is usually OBTINITIAL for obtStore.\n");
    fprintf(stderr, "                    clr is usually OBTCHIMERA for ovlStore when OBT is used.\n");
    fprintf(stderr, "                    clr is usually CLR        for ovlStore when OBT is not used.\n");
    fprintf(stderr, "\n");
    exit(1);
  }
  if (dumpType == 0)
    dumpType = DUMP_5p | DUMP_3p | DUMP_CONTAINED | DUMP_CONTAINS;

  switch (operation) {
    case OP_DUMP:
      dumpStore(storeName, dumpBinary, dumpERate, dumpLength, dumpType, bgnIID, endIID, qryIID, beVerbose);
      break;
    case OP_DUMP_PICTURE:
      dumpPicture(storeName, gkpName, clearRegion, dumpERate, dumpLength, dumpType, qryIID);
      break;
    case OP_GENOME_LENGTH:
      estimateGenomeLength(storeName,
                           gkpName,
                           gs_ovlLimit,
                           bgnIID,
                           endIID,
                           gs_into,
                           gs_winSize,
                           gs_minOvl);
      break;
    case OP_UPDATE_ERATES:
      updateErates(storeName, erateFile);
      break;
    default:
      break;
  }

  exit(0);
}
