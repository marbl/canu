
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

static char *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

#include "AS_global.H"
#include "AS_PER_genericStore.H"
#include "AS_PER_gkpStore.H"
#include "AS_PER_encodeSequenceQuality.H"
#include "AS_UTL_fileIO.H"

//  Special one-off to upgrade a v6 store to a v7 store.
//
//  The normal store routined are completely bypassed.  The upgrade reads in all the fpk data, and
//  splits it into fpk and qpk.  This is all done special case -- the normal data structures are all
//  private, and we cannot access them.
//
//  It MUST be run from within the gkpStore directory.


class gkStoreInfoOriginal {
public:
  uint64    gkMagic;
  uint64    gkVersion;

  uint32    gkLibrarySize;
  uint32    gkPackedFragmentSize;
  uint32    gkNormalFragmentSize;
  uint32    gkStrobeFragmentSize;
  uint32    gkPlacementSize;

  uint32    readMaxLenBits;

  //  Statistics on our load

  uint32    libInput;
  uint32    libLoaded;
  uint32    libErrors;
  uint32    libWarnings;

  uint32    frgInput;
  uint32    frgLoaded;
  uint32    frgErrors;
  uint32    frgWarnings;

  uint32    lkgInput;
  uint32    lkgLoaded;
  uint32    lkgErrors;
  uint32    lkgWarnings;
   
  uint32    sffInput;
  uint32    sffLoaded;
  uint32    sffErrors;
  uint32    sffWarnings;

  uint32    sffLibCreated;

  uint32    plcInput;
  uint32    plcLoaded;
  uint32    plcErrors;
  uint32    plcWarnings;

  //  Counts of types of things we have loaded

  uint32    numRandom;

  //  The IID space is broken into three classes.  See comments in AS_PER_gkStore_IID.C

  uint32    numPacked;
  uint32    numNormal;
  uint32    numStrobe;
};


class gkPackedFragmentOriginal {
public:
  AS_UID           readUID;

  AS_IID           readIID;
  AS_IID           mateIID;
  AS_IID           libraryIID;

  uint32           pad         : 4;
  uint32           deleted     : 1;
  uint32           nonrandom   : 1;
  uint32           orientation : 2;
  uint32           seqLen      : 8;  //  AS_READ_MAX_PACKED_LEN_BITS;
  uint32           clearBeg    : 8;  //  AS_READ_MAX_PACKED_LEN_BITS;
  uint32           clearEnd    : 8;  //  AS_READ_MAX_PACKED_LEN_BITS;

  char             enc[104];  //  AS_READ_MAX_PACKED_LEN
};




//  Because the real classes in AS_PER have all their data PRIVATE, we need to duplicate
//  them here.  This IS one-off code, remember.


class  gkPackedSequenceNew {
public:
  char             enc[AS_READ_MAX_PACKED_LEN];
};

class gkPackedFragmentNew {
public:
  AS_UID           readUID;

  AS_IID           readIID;
  AS_IID           mateIID;
  AS_IID           libraryIID;

  uint32           pad         : 4;
  uint32           deleted     : 1;
  uint32           nonrandom   : 1;
  uint32           orientation : 2;
  uint32           seqLen      : AS_READ_MAX_PACKED_LEN_BITS;
  uint32           clearBeg    : AS_READ_MAX_PACKED_LEN_BITS;
  uint32           clearEnd    : AS_READ_MAX_PACKED_LEN_BITS;
};





int
main(int argc, char **argv) {

  //  Fail if the backups (or the new file) are already there

  if (AS_UTL_fileExists("fpk.original", FALSE, FALSE))
    fprintf(stderr, "fpk backup file exists, cannot proceed (store already converted?)\n"), exit(1);

  if (AS_UTL_fileExists("qpk", FALSE, FALSE))
    fprintf(stderr, "qpk file exists, cannot proceed (store already converted?)\n"), exit(1);

  //  Make backups of the originals

  errno = 0;
  rename("fpk", "fpk.original");
  if (errno)
    fprintf(stderr, "Failed to rename 'fpk' to 'fpk.original': %s\n", strerror(errno)), exit(1);

  errno = 0;
  rename("inf", "inf.original");
  if (errno)
    fprintf(stderr, "Failed to rename 'inf' to 'inf.original': %s\n", strerror(errno)), exit(1);

  //  Make damn sure the originals aren't there.

  if (AS_UTL_fileExists("fpk", FALSE, FALSE))
    fprintf(stderr, "fpk file exists, cannot proceed (wanted to overwrite data!)\n"), exit(1);

  if (AS_UTL_fileExists("inf", FALSE, FALSE))
    fprintf(stderr, "inf file exists, cannot proceed (wanted to overwrite data!)\n"), exit(1);

  //  Update the inf block -- argh, we added one int, and need to copy every field.

  errno = 0;
  FILE *ioO = fopen("inf.original", "r");
  if (errno)
    fprintf(stderr, "failed to open 'inf.original': %s\n", strerror(errno)), exit(1);

  FILE *ioN = fopen("inf", "w");
  if (errno)
    fprintf(stderr, "failed to open 'inf': %s\n", strerror(errno)), exit(1);

  gkStoreInfoOriginal io;
  gkStoreInfo         in;

  if (1 != AS_UTL_safeRead(ioO, &io, "ioO", sizeof(gkStoreInfoOriginal), 1))
    fprintf(stderr, "failed to read 'inf.original': %s\n", strerror(errno)), exit(1);

  assert(io.gkVersion == 6);

  in.gkMagic              = io.gkMagic;
  in.gkVersion            = 7;  //  private to AS_PER_gkStore.C!!  AS_GKP_CURRENT_VERSION;
  in.gkLibrarySize        = io.gkLibrarySize;
  in.gkPackedSequenceSize = sizeof(gkPackedSequence);
  in.gkPackedFragmentSize = sizeof(gkPackedFragment);
  in.gkNormalFragmentSize = io.gkNormalFragmentSize;
  in.gkStrobeFragmentSize = io.gkStrobeFragmentSize;
  in.gkPlacementSize      = io.gkPlacementSize;
  in.readMaxLenBits       = io.readMaxLenBits;
  in.libInput             = io.libInput;
  in.libLoaded            = io.libLoaded;
  in.libErrors            = io.libErrors;
  in.libWarnings          = io.libWarnings;
  in.frgInput             = io.frgInput;
  in.frgLoaded            = io.frgLoaded;
  in.frgErrors            = io.frgErrors;
  in.frgWarnings          = io.frgWarnings;
  in.lkgInput             = io.lkgInput;
  in.lkgLoaded            = io.lkgLoaded;
  in.lkgErrors            = io.lkgErrors;
  in.lkgWarnings          = io.lkgWarnings;
  in.sffInput             = io.sffInput;
  in.sffLoaded            = io.sffLoaded;
  in.sffErrors            = io.sffErrors;
  in.sffWarnings          = io.sffWarnings;
  in.sffLibCreated        = io.sffLibCreated;
  in.plcInput             = io.plcInput;
  in.plcLoaded            = io.plcLoaded;
  in.plcErrors            = io.plcErrors;
  in.plcWarnings          = io.plcWarnings;
  in.numRandom            = io.numRandom;
  in.numPacked            = io.numPacked;
  in.numNormal            = io.numNormal;
  in.numStrobe            = io.numStrobe;

  AS_UTL_safeWrite(ioN, &in, "ioF", sizeof(gkStoreInfo), 1);

  //  Dump the rest of the data (this is more or less copied from AS_PER_gkStore.C).

  if (!feof(ioO)) {
    uint32 nr = io.numPacked + io.numNormal + io.numStrobe + 1;
    uint32 na = 0;
    uint32 nb = 0;

    uint8   *IIDtoTYPE = (uint8  *)safe_malloc(sizeof(uint8)  * nr);
    uint32  *IIDtoTIID = (uint32 *)safe_malloc(sizeof(uint32) * nr);

    na = AS_UTL_safeRead(ioO, IIDtoTYPE, "gkStore_open:header", sizeof(uint8), nr);
    nb = AS_UTL_safeRead(ioO, IIDtoTIID, "gkStore_open:header", sizeof(uint32), nr);

    //  If EOF was hit, and nothing was read, there is no index saved.  Otherwise, something was
    //  read, and we fail if either was too short.

    if ((feof(ioO)) && (na == 0) && (nb == 0)) {
      safe_free(IIDtoTYPE);
      safe_free(IIDtoTIID);
    } else if ((na != nr) || (nb != nr)) {
      fprintf(stderr, "couldn't read the IID maps: %s\n", strerror(errno)), exit(1);
    }

    AS_UTL_safeWrite(ioN, IIDtoTYPE, "ioF", sizeof(uint8),  na);
    AS_UTL_safeWrite(ioN, IIDtoTIID, "ioF", sizeof(uint32), nb);
  }

  fclose(ioO);
  fclose(ioN);


  ////////////////////////////////////////
  //
  //  Open the old store.

  StoreStruct *fpkoriginal = openStore("fpk.original", "r");

  //  Create the new store.

  StoreStruct *fpk = createIndexStore("fpk", "fpk", sizeof(gkPackedFragmentNew), 1);
  StoreStruct *qpk = createIndexStore("qpk", "fpk", sizeof(gkPackedSequenceNew), 1);

  //  Loop over the whole old store, loading data, copying to the new struct

  gkPackedFragmentOriginal  fo;

  gkPackedFragmentNew       fn;
  gkPackedSequenceNew       sn;

  int32 f = getFirstElemStore(fpkoriginal);
  int32 e = getLastElemStore(fpkoriginal);

  for (int32 i=f; i<=e; i++) {
    getIndexStore(fpkoriginal, i, &fo);

    //  Convert

    fn.readUID     = fo.readUID;
    fn.readIID     = fo.readIID;
    fn.mateIID     = fo.mateIID;
    fn.libraryIID  = fo.libraryIID;
    fn.pad         = fo.pad;
    fn.deleted     = fo.deleted;
    fn.nonrandom   = fo.nonrandom;
    fn.orientation = fo.orientation;
    fn.seqLen      = fo.seqLen;
    fn.clearBeg    = fo.clearBeg;
    fn.clearEnd    = fo.clearEnd;

    memset(sn.enc, 0,      sizeof(char) * AS_READ_MAX_PACKED_LEN);
    memcpy(sn.enc, fo.enc, sizeof(char) * 104);

    //  Write

    appendIndexStore(fpk, &fn);
    appendIndexStore(qpk, &sn);
  }

  //  Clean up

  closeStore(fpkoriginal);
  closeStore(fpk);
  closeStore(qpk);

  exit(0);
}
