
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
 *    src/AS_BAT/AS_BAT_Outputs.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2014-MAR-31
 *      are Copyright 2010-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUN-05
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-04
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-MAR-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_OverlapCache.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"

#include "AS_BAT_PlaceReadUsingOverlaps.H"

#include "tgStore.H"
#include <string>
#include <sstream>



void
writeTigsToStore(TigVector     &tigs,
                 char          *filePrefix,
                 char          *storeName,
                 bool           isFinal) {
  char        filename[FILENAME_MAX] = {0};

  snprintf(filename, FILENAME_MAX, "%s.%sStore", filePrefix, storeName);
  tgStore     *tigStore = new tgStore(filename);
  tgTig       *tig      = new tgTig;

  snprintf(filename, FILENAME_MAX, "%s.%s.paths.gfa", filePrefix, storeName);
  FILE *gfa_paths = fopen(filename, "w");
  for (uint32 ti = 0; ti < tigs.size(); ti++) {
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

    tig->_layoutLen       = utg->getLength();

    //  Transfer reads from the bogart tig to the output tig.

    resizeArray(tig->_children, tig->_childrenLen, tig->_childrenMax, utg->ufpath.size(), resizeArray_doNothing);

    std::string delim = "";
    std::stringstream segment_names;
    std::stringstream overlap_descs;

    for (uint32 ti=0; ti<utg->ufpath.size(); ti++) {
      ufNode        *frg   = &utg->ufpath[ti];

      tig->addChild()->set(frg->ident,
                           frg->parent, frg->ahang, frg->bhang,
                           frg->position.bgn, frg->position.end);
      segment_names << delim << frg->ident << (frg->position.bgn <= frg->position.end ? '+' : '-');
      overlap_descs << delim << '*';
    }

    fprintf(gfa_paths, "P\t%d\t%s\t%s\n", tig->_tigID, segment_names.str().c_str(), overlap_descs.str().c_str());

    //  And write to the store

    tigStore->insertTig(tig, false);
  }

  fclose(gfa_paths);

  delete    tig;
  delete    tigStore;
}
