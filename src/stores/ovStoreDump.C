
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_OVS/overlapStore.C
 *    src/AS_OVS/overlapStore.c
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2007-MAR-08 to 2013-AUG-01
 *      are Copyright 2007-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-DEC-16 to 2009-AUG-14
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Gregory Sims from 2012-FEB-01 to 2012-MAR-14
 *      are Copyright 2012 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-AUG-22 to 2015-JUN-25
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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
dumpStore(ovStore *ovlStore,
          gkStore *gkpStore,
          bool     asBinary,
          bool     asCounts,
          double   dumpERate,
          uint32   dumpType,
          uint32   dumpLength,
          uint32   bgnID,
          uint32   endID,
          uint32   qryID,
          ovOverlapDisplayType    type,
          bool     beVerbose) {

  ovOverlap     overlap(gkpStore);
  uint64         evalue = AS_OVS_encodeEvalue(dumpERate);
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

  uint32  *counts          = (asCounts) ? new uint32 [endID - bgnID + 1] : NULL;

  if (asCounts)
    for (uint32 ii=bgnID; ii<=endID; ii++)
      counts[ii - bgnID] = 0;

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

    //  The toString() method is quite slow, all from sprintf().
    //    Without both the puts() and AtoString(), a dump ran in 3 seconds.
    //    With both, 138 seconds.
    //    Without the puts(), 127 seconds.

    if (asCounts)
      counts[overlap.a_iid - bgnID]++;

    else if (asBinary)
      AS_UTL_safeWrite(stdout, &overlap, "dumpStore", sizeof(ovOverlap), 1);

    else
      fputs(overlap.toString(ovlString, type, true), stdout);
  }

  if (asCounts)
    for (uint32 ii=bgnID; ii<=endID; ii++)
      fprintf(stdout, "%u\t%u\n", ii + bgnID, counts[ii]);

  delete [] counts;

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
  ovOverlap const *A = (ovOverlap const *)a;
  ovOverlap const *B = (ovOverlap const *)b;

  if (A->a_bgn() < B->a_bgn())  return(-1);
  if (A->a_bgn() > B->a_bgn())  return(1);

  //  For overlaps off the 5' end, put the thinnest ones first.  Sadly, we don't know
  //  gkpStore here, and can't actually get the end coordinate.
  //
  //if (A->b_bgn() > B->b_bgn())  return(-1);
  //if (A->b_bgn() < B->b_bgn())  return(1);

  if (A->a_end() < B->a_end())  return(-1);
  if (A->a_end() > B->a_end())  return(1);

  return(0);
}




void
dumpPicture(ovOverlap *overlaps,
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

  qsort(overlaps, novl, sizeof(ovOverlap), sortOBT);

  //  Build ascii representations for each overlapping read.

  for (uint32 o=0; o<novl; o++) {
    gkRead   *B      = gkpStore->gkStore_getRead(overlaps[o].b_iid);
    uint32   frgLenB = B->gkRead_sequenceLength();

    //  Find bgn/end points on each read.  If the overlap is reverse complement,
    //  the B coords are flipped so that bgn > end.

    uint32   ovlBgnA = overlaps[o].a_bgn();
    uint32   ovlEndA = overlaps[o].a_end();

    uint32   ovlBgnB = overlaps[o].b_bgn();
    uint32   ovlEndB = overlaps[o].b_end();

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
            overlaps[o].erate() * 100.0,
            ovl);
  }
}




void
dumpPicture(ovStore  *ovlStore,
            gkStore  *gkpStore,
            double    dumpERate,
            uint32    dumpLength,
            uint32    dumpType,
            uint32    qryID) {

  //fprintf(stderr, "DUMPING PICTURE for ID "F_U32" in store %s (gkp %s)\n",
  //        qryID, ovlName, gkpName);

  gkRead   *A      = gkpStore->gkStore_getRead(qryID);
  uint32   frgLenA = A->gkRead_sequenceLength();

  ovlStore->setRange(qryID, qryID);

  uint64         novl     = 0;
  ovOverlap     overlap(gkpStore);
  ovOverlap    *overlaps = ovOverlap::allocateOverlaps(gkpStore, ovlStore->numOverlapsInRange());
  uint64         evalue   = AS_OVS_encodeEvalue(dumpERate);

  //  Load all the overlaps so we can sort by the A begin position.

  while (ovlStore->readOverlap(&overlap) == TRUE) {

    if (overlap.evalue() > evalue)
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

    if (overlap.b_end() - overlap.b_bgn() < dumpLength)
      continue;

    if (overlap.a_end() - overlap.a_bgn() < dumpLength)
      continue;

    overlaps[novl++] = overlap;
  }


  if (novl == 0)
    fprintf(stderr, "no overlaps to show.\n");
  else
    dumpPicture(overlaps, novl, gkpStore, qryID);

  delete [] overlaps;
}






int
main(int argc, char **argv) {
  uint32          operation   = OP_NONE;

  char           *gkpName     = NULL;
  char           *ovlName     = NULL;

  bool            asBinary    = false;
  bool            asCounts    = false;

  double          dumpERate   = 1.0;
  uint32          dumpLength  = 0;
  uint32          dumpType    = 0;

  char           *erateFile   = NULL;

  uint32          bgnID       = 0;
  uint32          endID       = UINT32_MAX;
  uint32          qryID       = 0;

  bool            beVerbose   = false;

  ovOverlapDisplayType  type = ovOverlapAsCoords;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {

    if      (strcmp(argv[arg], "-G") == 0)
      gkpName = argv[++arg];

    else if (strcmp(argv[arg], "-O") == 0)
      ovlName = argv[++arg];

    else if (strcmp(argv[arg], "-b") == 0)
      bgnID = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-e") == 0)
      endID = atoi(argv[++arg]);


    //  Standard bulk dump of overlaps
    else if (strcmp(argv[arg], "-d") == 0)
      operation  = OP_DUMP;

    //  Dump as a picture, the next ID
    //  Should be easy to extend to using -b -e range
    else if (strcmp(argv[arg], "-p") == 0) {
      operation  = OP_DUMP_PICTURE;
      bgnID      = atoi(argv[++arg]);
      endID      = bgnID;
      qryID      = bgnID;
    }

    //  Query if the overlap for the next two integers exists
    else if (strcmp(argv[arg], "-q") == 0) {
      operation  = OP_DUMP;
      bgnID      = atoi(argv[++arg]);
      endID      = bgnID;
      qryID      = atoi(argv[++arg]);
    }


    //  Format of the dump
    else if (strcmp(argv[arg], "-coords") == 0)
      type = ovOverlapAsCoords;

    else if (strcmp(argv[arg], "-hangs") == 0)
      type = ovOverlapAsHangs;

    else if (strcmp(argv[arg], "-raw") == 0)
      type = ovOverlapAsRaw;

    else if (strcmp(argv[arg], "-binary") == 0)
      asBinary = true;

    else if (strcmp(argv[arg], "-counts") == 0)
      asCounts = true;


    //  standard bulk dump options
    else if (strcmp(argv[arg], "-E") == 0)
      dumpERate = atof(argv[++arg]);

    else if (strcmp(argv[arg], "-L") == 0)
      dumpLength = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-d5") == 0)
      dumpType |= DUMP_5p;

    else if (strcmp(argv[arg], "-d3") == 0)
      dumpType |= DUMP_3p;

    else if (strcmp(argv[arg], "-dC") == 0)
      dumpType |= DUMP_CONTAINS;

    else if (strcmp(argv[arg], "-dc") == 0)
      dumpType |= DUMP_CONTAINED;

    else if (strcmp(argv[arg], "-v") == 0)
      beVerbose = true;


    else {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }

  if (operation == OP_NONE)
    err++;
  if (gkpName == NULL)
    err++;
  if (ovlName == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore [-b bgnID] [-e endID] ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "There are three modes of operation:\n");
    fprintf(stderr, "  -d         dump a store (range selected with -b and -e)\n");
    fprintf(stderr, "  -q a b     report the a,b overlap, if it exists.\n");
    fprintf(stderr, "  -p a       dump a picture of overlaps to fragment 'a'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  FORMAT (for -d)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -coords    dump overlap showing coordinates in the reads (default)\n");
    fprintf(stderr, "  -hangs     dump overlap showing dovetail hangs unaligned\n");
    fprintf(stderr, "  -raw       dump overlap showing its raw native format (four hangs)\n");
    fprintf(stderr, "  -binary    dump overlap as raw binary data\n");
    fprintf(stderr, "  -counts    dump the number of overlaps per read\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  MODIFIERS (for -d and -p)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -E erate          Dump only overlaps <= erate fraction error.\n");
    fprintf(stderr, "  -L length         Dump only overlaps that are larger than L bases (only for -p picture mode).\n");
    fprintf(stderr, "  -d5               Dump only overlaps off the 5' end of the A frag.\n");
    fprintf(stderr, "  -d3               Dump only overlaps off the 3' end of the A frag.\n");
    fprintf(stderr, "  -dC               Dump only overlaps that are contained in the A frag (B contained in A).\n");
    fprintf(stderr, "  -dc               Dump only overlaps that are containing the A frag (A contained in B).\n");
    fprintf(stderr, "  -v                Report statistics (to stderr) on some dumps (-d).\n");
    fprintf(stderr, "\n");

    if (operation == OP_NONE)
      fprintf(stderr, "ERROR: no operation (-d, -q or -p) supplied.\n");
    if (gkpName == NULL)
      fprintf(stderr, "ERROR: no input gkpStore (-G) supplied.\n");
    if (ovlName == NULL)
      fprintf(stderr, "ERROR: no input ovlStore (-O) supplied.\n");

    exit(1);
  }

  if (dumpType == 0)
    dumpType = DUMP_5p | DUMP_3p | DUMP_CONTAINED | DUMP_CONTAINS;

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);
  ovStore  *ovlStore = new ovStore(ovlName, gkpStore);

  if (endID > gkpStore->gkStore_getNumReads())
    endID = gkpStore->gkStore_getNumReads();

  if (endID < bgnID)
    fprintf(stderr, "ERROR: invalid bgn/end range bgn=%u end=%u; only %u reads in the store\n", bgnID, endID, gkpStore->gkStore_getNumReads()), exit(1);

  switch (operation) {
    case OP_DUMP:
      dumpStore(ovlStore, gkpStore, asBinary, asCounts, dumpERate, dumpLength, dumpType, bgnID, endID, qryID, type, beVerbose);
      break;
    case OP_DUMP_PICTURE:
      for (qryID=bgnID; qryID <= endID; qryID++)
        dumpPicture(ovlStore, gkpStore, dumpERate, dumpLength, dumpType, qryID);
      break;
    default:
      break;
  }

  delete ovlStore;

  gkpStore->gkStore_close();

  exit(0);
}
