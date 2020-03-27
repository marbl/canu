
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_TigVector.H"

#include "AS_BAT_ReadInfo.H"


void
promoteToSingleton(TigVector &tigs) {
  uint32  nPromoted = 0;

  for (uint32 fi=1; fi<=RI->numReads(); fi++) {
    if (tigs.inUnitig(fi) != 0)  // Placed.
      continue;

    if (RI->readLength(fi) == 0)  //  Deleted.
      continue;

    nPromoted++;

    Unitig *utg = tigs.newUnitig(false);
    ufNode  read;

    read.ident             = fi;
    read.contained         = 0;
    read.parent            = 0;
    read.ahang             = 0;
    read.bhang             = 0;
    read.position.bgn      = 0;
    read.position.end      = RI->readLength(fi);

    utg->addRead(read, 0, false);

    utg->_isUnassembled = true;
  }

  writeStatus("promoteToSingleton()-- Moved " F_U32 " unplaced read%s to singleton tigs.\n",
              nPromoted, (nPromoted == 1) ? "" : "s");
}
