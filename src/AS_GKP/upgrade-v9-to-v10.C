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

static const char* rcsid = "$Id:  $";

#include "AS_global.H"
#include "AS_PER_genericStore.H"
#include "AS_PER_gkpStore.H"
#include "AS_UTL_fileIO.H"
#include "AS_UTL_Hash.H"

#include <map>

using namespace std;


class gkLibraryOriginal {
public:
  AS_UID         libraryUID;
  char           libraryName[LIBRARY_NAME_SIZE];

  double         mean;
  double         stddev;

  uint64         spareUN3;
  uint64         spareUN2;
  uint64         spareUN1;

  uint64         spareUTG:62;
  uint64         forceBOGunitigger:1;
  uint64         isNotRandom:1;

  uint64         spareALN:63;
  uint64         doNotTrustHomopolymerRuns:1;

  uint64         spareOBT:53;

  uint64         doTrim_initialNone:1;
  uint64         doTrim_initialMerBased:1;
  uint64         doTrim_initialFlowBased:1;
  uint64         doTrim_initialQualityBased:1;

  uint64         doRemoveDuplicateReads:1;

  uint64         doTrim_finalLargestCovered:1;
  uint64         doTrim_finalEvidenceBased:1;
  uint64         doTrim_finalBestEdge:1;

  uint64         doRemoveSpurReads:1;
  uint64         doRemoveChimericReads:1;

  uint64         doConsensusCorrection:1;

  uint64         spareGKP:63;
  uint64         forceShortReadFormat:1;

  uint64         spareLIB:59;
  uint64         readsAreReversed:1;
  uint64         constantInsertSize:1;
  uint64         orientation:3;
};



int main(int argc, char** argv) {

  //  Fail if the backups are already there

  if (AS_UTL_fileExists("inf.original", FALSE, FALSE))
    fprintf(stderr, "inf backup file exists, cannot proceed (store already converted?)\n"), exit(1);

  if (AS_UTL_fileExists("lib.original", FALSE, FALSE))
    fprintf(stderr, "lib backup file exists, cannot proceed (store already converted?)\n"), exit(1);

  //  Make backups of the originals

  errno = 0;
  rename("inf", "inf.original");
  if (errno)
    fprintf(stderr, "Failed to rename 'inf' to 'inf.original': %s\n", strerror(errno)), exit(1);

  errno = 0;
  rename("lib", "lib.original");
  if (errno)
    fprintf(stderr, "Failed to rename 'lib' to 'lib.original': %s\n", strerror(errno)), exit(1);

  //  Make damn sure the originals aren't there.

  if (AS_UTL_fileExists("inf", FALSE, FALSE))
    fprintf(stderr, "inf file exists, cannot proceed (wanted to overwrite data!)\n"), exit(1);

  if (AS_UTL_fileExists("lib", FALSE, FALSE))
    fprintf(stderr, "lib file exists, cannot proceed (wanted to overwrite data!)\n"), exit(1);

  //  Update the inf block -- argh, just to change the version (and struct sizes, maybe)

  errno = 0;
  FILE *ioO = fopen("inf.original", "r");
  if (errno)
    fprintf(stderr, "failed to open 'inf.original': %s\n", strerror(errno)), exit(1);

  FILE *ioN = fopen("inf", "w");
  if (errno)
    fprintf(stderr, "failed to open 'inf': %s\n", strerror(errno)), exit(1);

  gkStoreInfo         io;  //  Original
  gkStoreInfo         in;  //  Updated

  if (1 != AS_UTL_safeRead(ioO, &io, "ioO", sizeof(gkStoreInfo), 1))
    fprintf(stderr, "failed to read 'inf.original': %s\n", strerror(errno)), exit(1);

  assert(io.gkVersion == 9);

  in               = io;

  in.gkVersion     = 10;

  AS_UTL_safeWrite(ioN, &in, "inf", sizeof(gkStoreInfo), 1);

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

  //
  //  Do the conversion
  //

  StoreStruct *oldLib = convertStoreToMemoryStore(openStore("lib.original", "r"));
  StoreStruct *newLib = createIndexStore("lib", "lib", sizeof(gkLibrary), 1);

  for (AS_IID iid=1; iid <= oldLib->lastElem - oldLib->firstElem + 1; iid++) {
    gkLibraryOriginal lo;
    gkLibrary         ln;
  
    getIndexStore(oldLib, iid, &lo);
  
    ln.libraryUID                 = lo.libraryUID;
    ln.mean                       = lo.mean;
    ln.stddev                     = lo.stddev;
    ln.spareUN3                   = lo.spareUN3;
    ln.spareUN2                   = lo.spareUN2;
    ln.spareUN1                   = lo.spareUN1;
    ln.spareUTG                   = lo.spareUTG;
    ln.forceBOGunitigger          = lo.forceBOGunitigger;
    ln.isNotRandom                = lo.isNotRandom;
    ln.spareALN                   = lo.spareALN;
    ln.doNotTrustHomopolymerRuns  = lo.doNotTrustHomopolymerRuns;
    ln.spareOBT                   = lo.spareOBT;
    ln.doTrim_initialNone         = lo.doTrim_initialNone;
    ln.doTrim_initialMerBased     = lo.doTrim_initialMerBased;
    ln.doTrim_initialFlowBased    = lo.doTrim_initialFlowBased;
    ln.doTrim_initialQualityBased = lo.doTrim_initialQualityBased;
    ln.doRemoveDuplicateReads     = lo.doRemoveDuplicateReads;
    ln.doTrim_finalLargestCovered = lo.doTrim_finalLargestCovered;
    ln.doTrim_finalEvidenceBased  = lo.doTrim_finalEvidenceBased;
    ln.doTrim_finalBestEdge       = lo.doTrim_finalBestEdge;
    ln.doRemoveSpurReads          = lo.doRemoveSpurReads;
    ln.doRemoveChimericReads      = lo.doRemoveChimericReads;
    ln.doCheckForSubReads         = 0;
    ln.doConsensusCorrection      = lo.doConsensusCorrection;
    ln.spareGKP                   = lo.spareGKP;
    ln.forceShortReadFormat       = lo.forceShortReadFormat;
    ln.spareLIB                   = lo.spareLIB;
    ln.readsAreReversed           = lo.readsAreReversed;
    ln.constantInsertSize         = lo.constantInsertSize;
    ln.orientation                = lo.orientation;

    strcpy(ln.libraryName, lo.libraryName);
  
    appendIndexStore(newLib, &ln);
  }
 
  closeStore(oldLib);
  closeStore(newLib);

  fprintf(stderr, "Success!\n");

  exit(0);
}
