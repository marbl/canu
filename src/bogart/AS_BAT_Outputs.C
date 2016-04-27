
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Outputs.H"



void
unitigToTig(tgTig       *tig,
            uint32       tigid,
            Unitig      *utg) {

  //  Initialize the output tig.

  tig->clear();

  tig->_tigID           = tigid;
  utg->_tigID           = tigid;

  tig->_coverageStat    = 1.0;  //  Default to just barely unique
  tig->_microhetProb    = 1.0;  //  Default to 100% probability of unique

  //  Set the class.

  if      (utg->_isUnassembled == true)
    tig->_class = tgTig_unassembled;

  else if (utg->_isBubble == true)
    tig->_class = tgTig_bubble;

  else
    tig->_class = tgTig_contig;

  tig->_suggestRepeat   = (utg->_isRepeat   == true);
  tig->_suggestCircular = (utg->_isCircular == true);

  tig->_layoutLen       = utg->getLength();

  //  Transfer reads from the bogart tig to the output tig.

  resizeArray(tig->_children, tig->_childrenLen, tig->_childrenMax, utg->ufpath.size(), resizeArray_doNothing);

  for (uint32 ti=0; ti<utg->ufpath.size(); ti++) {
    ufNode        *frg   = &utg->ufpath[ti];

    tig->addChild()->set(frg->ident,
                         frg->parent, frg->ahang, frg->bhang,
                         frg->position.bgn, frg->position.end);
  }
}



void
writeUnitigsToStore(UnitigVector  &unitigs,
                    char          *fileprefix,
                    char          *tigStorePath,
                    uint32         frg_count_target,
                    bool           isFinal) {
  uint32      utg_count              = 0;
  uint32      frg_count              = 0;
  uint32      prt_count              = 1;
  char        filename[FILENAME_MAX] = {0};

  // Open up the initial output file

  sprintf(filename, "%s.iidmap", fileprefix);
  FILE *iidm = fopen(filename, "w");
  assert(NULL != iidm);

  sprintf(filename, "%s.partitioning", fileprefix);
  FILE *part = fopen(filename, "w");
  assert(NULL != part);

  sprintf(filename, "%s.partitioningInfo", fileprefix);
  FILE *pari = fopen(filename, "w");
  assert(NULL != pari);

  //  Step through all the unitigs once to build the partition mapping and IID mapping.

  tgStore     *tigStore = new tgStore(tigStorePath);
  tgTig       *tig      = new tgTig;

  for (uint32 tigID=0, ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if ((utg == NULL) || (utg->getNumFrags() == 0))
      continue;

    assert(utg->getLength() > 0);

    //  Convert the bogart tig to a tgTig and save to the store.

    unitigToTig(tig, (isFinal) ? tigID : ti, utg);
    tigID++;

    tigStore->insertTig(tig, false);

    //  Increment the partition if the current one is too large.

    if ((frg_count + utg->getNumFrags() >= frg_count_target) &&
        (frg_count                      >  0)) {
      fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",
              prt_count, utg_count, frg_count);

      prt_count++;
      utg_count = 0;
      frg_count = 0;
    }

    //  Note that the tig is included in this partition.

    utg_count += 1;
    frg_count += utg->getNumFrags();

    //  Map the tig to a partition, and log both the tig-to-partition map and the partition-to-read map.

    fprintf(iidm, "bogart "F_U32" -> tig "F_U32" (in partition "F_U32" with "F_U32" frags)\n",
            utg->id(),
            utg->tigID(),
            prt_count,
            utg->getNumFrags());

    for (uint32 fragIdx=0; fragIdx<utg->getNumFrags(); fragIdx++)
      fprintf(part, "%d\t%d\n", prt_count, utg->ufpath[fragIdx].ident);
  }

  fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",   //  Don't forget to log the last partition!
          prt_count, utg_count, frg_count);

  fclose(pari);
  fclose(part);
  fclose(iidm);

  delete    tig;
  delete    tigStore;
}
