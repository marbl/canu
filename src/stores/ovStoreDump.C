
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

#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"


enum dumpOp {
  OP_NONE           = 1,
  OP_DUMP           = 2,
  OP_DUMP_PICTURE   = 3
};


enum dumpFlags {
  DUMP_5p         = 1,
  DUMP_3p         = 2,
  DUMP_CONTAINED  = 4,
  DUMP_CONTAINS   = 8
};


//
//  Also accept a single ovStoreFile (output from overlapper) and dump.
//

//
//  Then need some way of loading ascii overlaps into a store, or converting ascii overlaps to
//  binary and use the normal store build.  The normal store build also needs to take sorted
//  overlaps and just rewrite as a store.
//

void
dumpStore(char   *ovlName,
          char   *gkpName,
          uint32  asBinary,
          double  dumpERate,
          uint32  dumpType,
          uint32  dumpLength,
          uint32  bgnID,
          uint32  endID,
          uint32  qryID,
          bool    asCoords,
          bool    beVerbose) {
  ovStore       *ovlStore = new ovStore(ovlName);
  gkStore       *gkpStore = new gkStore(gkpName);

  ovsOverlap     overlap;
  uint64         evalue = AS_OVS_encodeQuality(dumpERate / 100.0);
  char           ovlString[1024];

  uint32   ovlTooHighError = 0;
  uint32   ovlNot5p        = 0;
  uint32   ovlNot3p        = 0;
  uint32   ovlNotContainer = 0;
  uint32   ovlNotContainee = 0;
  uint32   ovlDumped       = 0;
  uint32   obtTooHighError = 0;
  uint32   obtDumped       = 0;
  uint32   merDumped       = 0;

  ovlStore->setRange(bgnID, endID);

  //  Length filtering is expensive to compute, need to load both reads to get their length.
  //
  //if ((dumpLength > 0) && (dumpLength < overlapLength(overlap)))
  //  continue;
  
  while (ovlStore->readOverlap(&overlap) == TRUE) {
    if ((qryID != 0) && (qryID != overlap.b_iid))
      continue;

    if (overlap.evalue() > evalue) {
      ovlTooHighError++;
      continue;
    }

    int32 ahang = overlap.a_hang();
    int32 bhang = overlap.b_hang();

    if (((dumpType & DUMP_5p) == 1) && (ahang < 0) && (bhang < 0)) {
      ovlNot5p++;
      continue;
    }
    if (((dumpType & DUMP_3p) == 1) && (ahang > 0) && (bhang > 0)) {
      ovlNot3p++;
      continue;
    }
    if (((dumpType & DUMP_CONTAINS) == 1) && (ahang >= 0) && (bhang <= 0)) {
      ovlNotContainer++;
      continue;
    }
    if (((dumpType & DUMP_CONTAINED) == 1) && (ahang <= 0) && (bhang >= 0)) {
      ovlNotContainee++;
      continue;
    }

    ovlDumped++;

    //  All the slow for this dump is in sprintf() within AS_OVS_toString().
    //    Without both the puts() and AS_OVS_toString(), a dump ran in 3 seconds.
    //    With both, 138 seconds.
    //    Without the puts(), 127 seconds.

    if (asBinary)
      AS_UTL_safeWrite(stdout, &overlap, "dumpStore", sizeof(ovsOverlap), 1);
    else
      fputs(overlap.toString(ovlString, gkpStore, asCoords), stdout);
  }

  delete ovlStore;
  delete gkpStore;

  if (beVerbose) {
    fprintf(stderr, "ovlTooHighError %u\n",  ovlTooHighError);
    fprintf(stderr, "ovlNot5p        %u\n",  ovlNot5p);
    fprintf(stderr, "ovlNot3p        %u\n",  ovlNot3p);
    fprintf(stderr, "ovlNotContainer %u\n",  ovlNotContainer);
    fprintf(stderr, "ovlNotContainee %u\n",  ovlNotContainee);
    fprintf(stderr, "ovlDumped       %u\n",  ovlDumped);
    fprintf(stderr, "obtTooHighError %u\n",  obtTooHighError);
    fprintf(stderr, "obtDumped       %u\n",  obtDumped);
    fprintf(stderr, "merDumped       %u\n",  merDumped);
  }
}



int
sortOBT(const void *a, const void *b) {
  ovsOverlap const *A = (ovsOverlap const *)a;
  ovsOverlap const *B = (ovsOverlap const *)b;

  if (A->a_bgn() < B->a_bgn())  return(-1);
  if (A->a_bgn() > B->a_bgn())  return(1);

  //  Don't know gkpStore here, dang.
  //if (A->a_end(gkpStore) < B->a_end(gkpStore))  return(-1);
  //if (A->a_end(gkpStore) > B->a_end(gkpStore))  return(1);

  return(0);
}




void
dumpPicture(ovsOverlap *overlaps,
            uint64      novl,
            gkStore    *gkpStore,
            uint32      qryID) {
  char     ovl[256] = {0};

  uint32   MHS = 7;  //  Max Hang Size, amount of padding for "+### "

  gkRead   *A      = gkpStore->gkStore_getRead(qryID);
  uint32   frgLenA = A->gkRead_sequenceLength();

  for (int32 i=0; i<256; i++)
    ovl[i] = ' ';

  for (int32 i=0; i<100; i++)
    ovl[i + MHS] = '-';
  ovl[ 99 + MHS] = '>';
  ovl[100 + MHS] = 0;

  fprintf(stdout, "%8d  A: %5d %5d                                           %s\n",
          qryID,
          0, frgLenA,
          ovl);

  qsort(overlaps, novl, sizeof(ovsOverlap), sortOBT);

  //  Build ascii representations for each overlapping read.

  for (uint32 o=0; o<novl; o++) {
    gkRead   *B      = gkpStore->gkStore_getRead(overlaps[o].b_iid);
    uint32   frgLenB = B->gkRead_sequenceLength();

    //  Find bgn/end points on each read.  If the overlap is reverse complement,
    //  the B coords are flipped so that bgn > end.

    uint32   ovlBgnA = overlaps[o].a_bgn(gkpStore);
    uint32   ovlEndA = overlaps[o].a_end(gkpStore);

    uint32   ovlBgnB = overlaps[o].b_bgn(gkpStore);
    uint32   ovlEndB = overlaps[o].b_end(gkpStore);

    assert(ovlBgnA < ovlEndA);  //  The A coordiantes are always forward

    if (overlaps[o].flipped() == false)
      assert(ovlBgnB < ovlEndB);  //  Forward overlaps are forward
    else
      assert(ovlEndB < ovlBgnB);  //  Flipped overlaps are reversed

    //  For the A read, find the points in our string representation where the overlap ends.

    uint32 ovlStrBgn = ovlBgnA * 100 / frgLenA + MHS;
    uint32 ovlStrEnd = ovlEndA * 100 / frgLenA + MHS;

    //  Fill the string representation with spaces, then fill the string with dashes where the read
    //  is, add an arrow, and terminate the string.

    for (int32 i=0; i<256; i++)
      ovl[i] = ' ';

    for (uint32 i=ovlStrBgn; i<ovlStrEnd; i++)
      ovl[i] = '-';

    if (overlaps[o].flipped() == true)
      ovl[ovlStrBgn] = '<';
    else
      ovl[ovlStrEnd-1] = '>';

    ovl[ovlStrEnd] = 0;

    //  For the B read, find how much is unaliged on each end.  Though the store directly keeps this information,
    //  we can't get to it, and have to reverse the compuitation.

    uint32  ovlBgnHang = 0;
    uint32  ovlEndHang = 0;

    if (overlaps[o].flipped() == false) {
      ovlBgnHang = ovlBgnB;
      ovlEndHang = frgLenB - ovlEndB;
    } else {
      ovlBgnHang = frgLenB - ovlBgnB;
      ovlEndHang = ovlEndB;
    }

    if (ovlBgnHang > 0) {
      char  str[256];
      int32 len;

      sprintf(str, "+%d", ovlBgnHang);
      len = strlen(str);

      for (int32 i=0; i<len; i++)
        ovl[ovlStrBgn - len - 1 + i] = str[i];
    }

    if (ovlEndHang > 0) {
      sprintf(ovl + ovlStrEnd, " +%d", ovlEndHang);
    }


    fprintf(stdout, "%8d  A: %5d %5d (%5d)  B: %5d %5d (%5d)  %5.2f%%   %s\n",
            overlaps[o].b_iid,
            ovlBgnA, ovlEndA, frgLenA,
            ovlBgnB, ovlEndB, frgLenB,
            overlaps[o].erate(),
            ovl);
  }
}




void
dumpPicture(char   *ovlName,
            char   *gkpName,
            double  dumpERate,
            uint32  dumpLength,
            uint32  dumpType,
            uint32  qryID) {

  fprintf(stderr, "DUMPING PICTURE for ID "F_U32" in store %s (gkp %s)\n",
          qryID, ovlName, gkpName);

  ovStore  *ovlStore = new ovStore(ovlName);
  gkStore  *gkpStore = new gkStore(gkpName);

  gkRead   *A      = gkpStore->gkStore_getRead(qryID);
  uint32   frgLenA = A->gkRead_sequenceLength();

  ovlStore->setRange(qryID, qryID);

  uint64         novl     = 0;
  ovsOverlap     overlap;
  ovsOverlap    *overlaps = (ovsOverlap *)safe_malloc(sizeof(ovsOverlap) * ovlStore->numOverlapsInRange());
  uint64         erate    = AS_OVS_encodeQuality(dumpERate / 100.0);

  //  Load all the overlaps so we can sort by the A begin position.

  while (ovlStore->readOverlap(&overlap) == TRUE) {

    if (overlap.erate() > erate)
      continue;

    if (((dumpType & DUMP_5p) == 0) &&
        (overlap.a_hang() < 0) && (overlap.b_hang() < 0))
      continue;

    if (((dumpType & DUMP_3p) == 0) &&
        (overlap.a_hang() > 0) && (overlap.b_hang() > 0))
      continue;

    if (((dumpType & DUMP_CONTAINS) == 0) &&
        (overlap.a_hang() >= 0) && (overlap.b_hang() <= 0))
      continue;

    if (((dumpType & DUMP_CONTAINED) == 0) &&
        (overlap.a_hang() <= 0) && (overlap.b_hang() >= 0))
      continue;

    if (overlap.b_end(gkpStore) - overlap.b_bgn(gkpStore) < dumpLength)
      continue;

    if (overlap.a_end(gkpStore) - overlap.a_bgn(gkpStore) < dumpLength)
      continue;

    overlaps[novl++] = overlap;
  }


  if (novl == 0)
    fprintf(stderr, "no overlaps to show.\n");
  else
    dumpPicture(overlaps, novl, gkpStore, qryID);


  delete ovlStore;
  delete gkpStore;
}




int
main(int argc, char **argv) {
  uint32          operation   = OP_NONE;

  char           *gkpName     = NULL;
  char           *ovlName     = NULL;

  uint32          asBinary    = false;
  double          dumpERate   = 100.0;
  uint32          dumpLength  = 0;
  uint32          dumpType    = 0;

  char           *erateFile   = NULL;

  uint32          bgnID       = 0;
  uint32          endID       = UINT32_MAX;
  uint32          qryID       = 0;

  bool            beVerbose   = false;

  bool            asCoords    = true;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {

    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-b") == 0) {
      bgnID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      endID = atoi(argv[++arg]);



      //  Dump as a picture, the next ID
      //  Should be easy to extend to using -b -e range
    } else if (strcmp(argv[arg], "-p") == 0) {
      operation  = OP_DUMP_PICTURE;
      bgnID      = atoi(argv[++arg]);
      endID      = bgnID;
      qryID      = bgnID;

      //  Query if the overlap for the next two integers exists
    } else if (strcmp(argv[arg], "-q") == 0) {
      operation  = OP_DUMP;
      bgnID      = atoi(argv[++arg]);
      endID      = bgnID;
      qryID      = atoi(argv[++arg]);


      //  Standard bulk dump of overlaps
    } else if (strcmp(argv[arg], "-d") == 0) {
      operation  = OP_DUMP;

    } else if (strcmp(argv[arg], "-coords") == 0) {
      asCoords = true;

    } else if (strcmp(argv[arg], "-hangs") == 0) {
      asCoords = false;

    } else if (strcmp(argv[arg], "-binary") == 0) {
      asBinary = true;

      //  standard bulk dump options
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

    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;


    } else {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((operation == OP_NONE) || (gkpName == NULL) || (ovlName == NULL) || (err)) {
    fprintf(stderr, "usage: %s -d storeName [-B] [-E erate] [-b beginID] [-e endID]\n", argv[0]);
    fprintf(stderr, "       %s -q aiid biid storeName\n", argv[0]);
    fprintf(stderr, "       %s -p iid storeName gkpStore clr\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   NOT UP TO DATE\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "There are three modes of operation, selected by the first option:\n");
    fprintf(stderr, "  -d  dump a store\n");
    fprintf(stderr, "  -q  report the a,b overlap, if it exists.\n");
    fprintf(stderr, "  -p  dump a picture of overlaps to fragment 'iid', using clear region 'clr'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DUMPING - report overlaps in the store, or in a binary overlap file\n");
    fprintf(stderr, "  -B                Dump the store as binary, suitable for input to create a new store.\n");
    fprintf(stderr, "  -E erate          Dump only overlaps <= erate error.\n");
    fprintf(stderr, "  -L length         Dump only overlaps that are larger than L bases (only for -p picture mode).\n");
    fprintf(stderr, "  -d5               Dump only overlaps off the 5' end of the A frag.\n");
    fprintf(stderr, "  -d3               Dump only overlaps off the 3' end of the A frag.\n");
    fprintf(stderr, "  -dC               Dump only overlaps that are contained in the A frag (B contained in A).\n");
    fprintf(stderr, "  -dc               Dump only overlaps that are containing the A frag (A contained in B).\n");
    fprintf(stderr, "  -b beginID        Start dumping at 'beginID'.\n");
    fprintf(stderr, "  -e endID          Stop dumping after 'endID'.\n");
    fprintf(stderr, "  -v                Report statistics (to stderr) on some dumps (-d).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -coords           default\n");
    fprintf(stderr, "  -hangs            \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "QUERYING - quickly ask if an overlap exists\n");
    fprintf(stderr, "  -q aiid biid storeName\n");
    fprintf(stderr, "                    If an overlap between fragments 'aiid' and 'biid' exists, it is printed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "PICTURES - draw a multi-alignment-like picture for a single fragment and its overlaps\n");
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
      dumpStore(ovlName, gkpName, asBinary, dumpERate, dumpLength, dumpType, bgnID, endID, qryID, asCoords, beVerbose);
      break;
    case OP_DUMP_PICTURE:
      for (qryID=bgnID; qryID <= endID; qryID++)
        dumpPicture(ovlName, gkpName, dumpERate, dumpLength, dumpType, qryID);
      break;
    default:
      break;
  }

  exit(0);
}
