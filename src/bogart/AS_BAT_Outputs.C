
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static const char *rcsid = "$Id$";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Instrumentation.H"

#include "AS_BAT_Outputs.H"

//#include "AS_CGB_histo.H"
//#include "tgStore.H"




//  Massage the Unitig into a MultiAlignT (also used in SplitChunks_CGW.c)
void
unitigToTig(tgTig       *tig,
            uint32       tigid,
            Unitig      *utg) {

  tig->clear();

  tig->_tigID           = tigid;
  tig->_coverageStat    = 1.0;  //  Default to just barely unique
  tig->_microhetProb    = 1.0;  //  Default to 100% probability of unique

  tig->_suggestRepeat   = false;
  tig->_suggestUnique   = false;
  tig->_suggestCircular = false;
  tig->_suggestHaploid  = false;

  tig->_layoutLen       = utg->getLength();

  resizeArray(tig->_children, tig->_childrenLen, tig->_childrenMax, utg->ufpath.size(), resizeArray_doNothing);

  for (uint32 fi=0; fi<utg->ufpath.size(); fi++) {
    ufNode        *frg = &utg->ufpath[fi];
    tgPosition    *pos =  tig->addChild();

    pos->set(frg->ident, frg->parent, frg->ahang, frg->bhang, frg->position.bgn, frg->position.end);
  }

  //fprintf(stderr, "unitigToTig()--  tig %u has %u children\n", tig->_tigID, tig->_childrenLen);
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
  uint32     *partmap                = new uint32 [unitigs.size()];

  //  This code closely follows that in AS_CGB_unitigger.c::output_the_chunks()

  if (isFinal)
    checkUnitigMembership(unitigs);

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

  memset(partmap, 0xff, sizeof(uint32) * unitigs.size());

  for (uint32 iumiid=0, ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];
    uint32   nf  = (utg) ? utg->getNumFrags() : 0;

    if ((utg == NULL) || (nf == 0))
      continue;

    assert(utg->getLength() > 0);
    assert(nf == utg->ufpath.size());

    if ((frg_count + nf >= frg_count_target) &&
        (frg_count      >  0)) {
      fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",
              prt_count, utg_count, frg_count);

      prt_count++;
      utg_count = 0;
      frg_count = 0;
    }

    uint32 tigid = (isFinal) ? iumiid : ti;

    assert(tigid < unitigs.size());
    partmap[tigid] = prt_count;

    fprintf(iidm, "Unitig "F_U32" == IUM "F_U32" (in partition "F_U32" with "F_U32" frags)\n",
            utg->id(),
            (tigid),
            partmap[(tigid)],
            nf);

    for (uint32 fragIdx=0; fragIdx<nf; fragIdx++) {
      ufNode  *f = &utg->ufpath[fragIdx];

      fprintf(part, "%d\t%d\n", prt_count, f->ident);
    }

    utg_count += 1;
    frg_count += nf;

    iumiid++;
  }

  fprintf(pari, "Partition %d has %d unitigs and %d fragments.\n",
          prt_count, utg_count, frg_count);

  fclose(pari);
  fclose(part);
  fclose(iidm);

  //  Step through all the unitigs once to build the partition mapping and IID mapping.

  tgStore     *tigStore = new tgStore(tigStorePath);
  tgTig       *tig      = new tgTig;

  for (uint32 iumiid=0, ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];
    uint32   nf  = (utg) ? utg->getNumFrags() : 0;

    if ((utg == NULL) || (nf == 0))
      continue;

    unitigToTig(tig, (isFinal) ? iumiid : ti, utg);

    tigStore->insertTig(tig, false);

    iumiid++;
  }

  delete    tig;
  delete    tigStore;
  delete [] partmap;
}


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
//
//  Wow, this is ancient.
//
void
writeOverlapsUsed(UnitigVector &unitigs,
                  char         *fileprefix) {
  char         filename[FILENAME_MAX] = {0};
#if 0
  GenericMesg  pmesg;
  OverlapMesg  omesg;
#endif

  sprintf(filename, "%s.unused.ovl", fileprefix);
  FILE *file = fopen(filename, "w");
  assert(file != NULL);

#if 0
  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->ufpath.size(); fi++) {
      ufNode  *frg = &utg->ufpath[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = OG->getBestEdgeOverlap(frg->ident, false);
      BestEdgeOverlap *bestedge3 = OG->getBestEdgeOverlap(frg->ident, true);

      int              bestident5 = 0;
      int              bestident3 = 0;

      if (bestedge5) {
        bestident5 = bestedge5->fragId();

        if ((bestident5 > 0) && (utg->fragIn(bestident5) != utg->id())) {
          omesg.aifrag          = frg->ident;
          omesg.bifrag          = bestident5;
          omesg.ahg             = bestedge5->ahang();
          omesg.bhg             = bestedge5->bhang();
          omesg.orientation.setIsUnknown();
          omesg.overlap_type    = AS_DOVETAIL;
          omesg.quality         = 0.0;
          omesg.min_offset      = 0;
          omesg.max_offset      = 0;
          omesg.polymorph_ct    = 0;
          omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
          omesg.alignment_delta = NULL;
#endif

          //  This overlap is off of the 5' end of this fragment.
          if (bestedge5->frag3p() == false)
            omesg.orientation.setIsOuttie();
          if (bestedge5->frag3p() == true)
            omesg.orientation.setIsAnti();

          pmesg.t = MESG_OVL;
          pmesg.m = &omesg;

          WriteProtoMesg_AS(file, &pmesg);
        }
      }

      if (bestedge3) {
        bestident3 = bestedge3->fragId();

        if ((bestident3 > 0) && (utg->fragIn(bestident3) != utg->id())) {
          omesg.aifrag          = frg->ident;
          omesg.bifrag          = bestident3;
          omesg.ahg             = bestedge3->ahang();
          omesg.bhg             = bestedge3->bhang();
          omesg.orientation.setIsUnknown();
          omesg.overlap_type    = AS_DOVETAIL;
          omesg.quality         = 0.0;
          omesg.min_offset      = 0;
          omesg.max_offset      = 0;
          omesg.polymorph_ct    = 0;
          omesg.alignment_trace = NULL;
#ifdef AS_MSG_USE_OVL_DELTA
          omesg.alignment_delta = NULL;
#endif

          //  This overlap is off of the 3' end of this fragment.
          if (bestedge3->frag3p() == false)
            omesg.orientation.setIsNormal();
          if (bestedge3->frag3p() == true)
            omesg.orientation.setIsInnie();

          pmesg.t = MESG_OVL;
          pmesg.m = &omesg;

          WriteProtoMesg_AS(file, &pmesg);
        }
      }
    }
  }
#endif

  fclose(file);
}
