
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

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_PlaceReadUsingOverlaps.H"

#include "tgStore.H"



void
writeTigsToStore(TigVector     &tigs,
                 char const    *filePrefix,
                 char const    *storeName,
                 bool           isFinal) {
  char        filename[FILENAME_MAX] = {0};

  snprintf(filename, FILENAME_MAX, "%s.%sStore", filePrefix, storeName);
  tgStore     *tigStore = new tgStore(filename);
  tgTig       *tig      = new tgTig;

  for (uint32 ti=0; ti<tigs.size(); ti++) {
    Unitig  *utg = tigs[ti];

    if ((utg == NULL) || (utg->getNumReads() == 0))
      continue;

    assert(utg->getLength() > 0);

    //  Initialize the output tig.

    tig->clear();

    tig->_tigID           = utg->id();

    //  Set the class and some flags.

    tig->_class           = (utg->_isUnassembled == true) ? tgTig_unassembled : tgTig_contig;
    tig->_suggestRepeat   = utg->_isRepeat;
    tig->_suggestCircular = utg->_isCircular;
    tig->_suggestBubble   = utg->_isBubble;
    tig->_circularLength  = utg->_circularLength;

    tig->_layoutLen       = utg->getLength();

    tig->_trimBgn         = 0;
    tig->_trimEnd         = utg->getLength();

    //  Transfer reads from the bogart tig to the output tig.

    tig->allocateChildren(utg->ufpath.size());

    for (uint32 ti=0; ti<utg->ufpath.size(); ti++) {
      ufNode        *frg   = &utg->ufpath[ti];

      tig->addChild()->set(frg->ident,
                           frg->parent, frg->ahang, frg->bhang,
                           frg->position.bgn, frg->position.end);
    }

    //  And write to the store

    tigStore->insertTig(tig, false);
  }

  delete    tig;
  delete    tigStore;
}
