
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
 *    src/AS_BAT/AS_BAT_Unitig.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"

static std::map<uint32,int>* containPartialOrder;

uint32* Unitig::_inUnitig     = NULL;
uint32* Unitig::_pathPosition = NULL;


void
Unitig::reverseComplement(bool doSort) {

  //  If there are contained fragments, we need to sort by position to place them correctly after
  //  their containers.  If there are no contained fragments, sorting can break the initial unitig
  //  building.  When two frags start at position zero, we'll exchange the order.  Initial unitig
  //  building depends on having the first fragment added become the last fragment in the unitig
  //  after reversing.

  for (uint32 fi=0; fi<ufpath.size(); fi++) {
    ufNode  *frg = &ufpath[fi];

    frg->position.bgn = getLength() - frg->position.bgn;
    frg->position.end = getLength() - frg->position.end;

    //if (frg->contained != 0)
    //  doSort = true;

    assert(frg->position.bgn >= 0);
    assert(frg->position.end >= 0);
  }

  //  We've updated the positions of everything.  Now, sort or reverse the list, and rebuild the
  //  pathPosition map.

  if (doSort) {
    sort();
  } else {
    std::reverse(ufpath.begin(), ufpath.end());

    for (uint32 fi=0; fi<ufpath.size(); fi++)
      _pathPosition[ufpath[fi].ident] = fi;
  }
}




class olapDat {
public:
  olapDat(uint32 b, uint32 e, double er) {
    bgn    = b;
    end    = e;
    erate  = er;
  };

  bool operator<(const olapDat &that)     const { return(bgn < that.bgn); };

  uint32  bgn;
  uint32  end;
  double  erate;
};


void
Unitig::computeErrorProfile(void) {
  //uint32          ti = id();


  //fprintf(stderr, "Find error profile for tig "F_U32".\n", id());

  vector<olapDat>  olaps;
  vector<uint32>   coords;

  for (uint32 fi=0; fi<ufpath.size(); fi++) {
    ufNode     *rdA    = &ufpath[fi];
    bool        rdAfwd = (rdA->position.bgn < rdA->position.end);
    int32       rdAlo  = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
    int32       rdAhi  = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

    uint32      ovlLen =  0;
    BAToverlap *ovl    =  OC->getOverlaps(rdA->ident, AS_MAX_ERATE, ovlLen);

    for (uint32 oi=0; oi<ovlLen; oi++) {
      if (id() != Unitig::fragIn(ovl[oi].b_iid))
        //  Reads in different tigs.  Don't care.
        continue;

      ufNode  *rdB    = &ufpath[ Unitig::pathPosition(ovl[oi].b_iid) ];
      bool     rdBfwd = (rdB->position.bgn < rdB->position.end);
      int32    rdBlo  = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
      int32    rdBhi  = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

      if ((rdAhi < rdBlo) || (rdBhi < rdAlo))
        continue;

      //  Now figure out what region is covered by the overlap.

      int32    tiglo = 0;
      int32    tighi = FI->fragmentLength(rdA->ident);

      if (ovl[oi].a_hang > 0)
        tiglo += ovl[oi].a_hang;  //  Postiive hang!

      if (ovl[oi].b_hang < 0)
        tighi += ovl[oi].b_hang;  //  Negative hang!

      assert(0     <= tiglo);
      assert(0     <= tighi);
      assert(tiglo <= tighi);
      assert(tiglo <= FI->fragmentLength(rdA->ident));
      assert(tighi <= FI->fragmentLength(rdA->ident));

      //  Offset and adjust to tig coordinates

      //  Beacuse the read is placed with a lot of fudging in the positions, we need
      //  to scale the coordinates we compute here.
      double      sc = (rdAhi - rdAlo) / (double)FI->fragmentLength(rdA->ident);

      olapDat     olap((uint32)floor(rdAlo + sc * tiglo),
                       (uint32)floor(rdAlo + sc * tighi),
                       ovl[oi].erate);

      coords.push_back(olap.bgn);
      coords.push_back(olap.end);

      olaps.push_back(olap);
    }
  }

  std::sort(coords.begin(), coords.end());
  std::sort(olaps.begin(),  olaps.end());

  //fprintf(stderr, "Generated "F_SIZE_T" coords and "F_SIZE_T" olaps.\n", coords.size(), olaps.size());

  //  Convert coordinates into intervals.  Conceptually, squish out the duplicate numbers, then
  //  create an interval for every adjacent pair.

  for (uint32 bb=0, ii=1; ii<coords.size(); ii++) {
    if (coords[bb] == coords[ii])
      continue;

    errorProfile.push_back(epValue(coords[bb], coords[ii]));

    bb = ii;
  }

  //fprintf(stderr, "Generated "F_SIZE_T" profile regions.\n", profile.size());

  //  Transfer overlap error rates to the regions.  Find the first region we intersect (by
  //  construction, the first interval will have olap.bgn == profile.bgn), then add the overlap
  //  error rate to each interval it covers (again, by construction, the last interval we add to
  //  will have olap.end == profile.end).

  uint64   nPieces = 0;

  for (uint32 bb=0, ii=0; ii<olaps.size(); ii++) {

    //  Find first region for this overlap.  All later overlaps will start at or after this region.
    for (; (bb < errorProfile.size()) && (errorProfile[bb].bgn < olaps[ii].bgn); bb++)
      ;

    //if (ii > 1771000)
    //  fprintf(stderr, "olap ii=%u %u-%u bb=%u profile %u-%u\n", ii, olaps[ii].bgn, olaps[ii].end, bb, profile[bb].bgn, profile[bb].end);

    //  Until the overlap stops overlapping, add its error rate to the regions.
    for (uint32 ee=bb; (ee < errorProfile.size()) && (errorProfile[ee].end <= olaps[ii].end); ee++, nPieces++)
      errorProfile[ee].dev.update(olaps[ii].erate);
  }

  //fprintf(stderr, "Generated "F_SIZE_T" profile regions with "F_U64" overlap pieces.\n", profile.size(), nPieces);
}



void
Unitig::reportErrorProfile(char *prefix) {
  char  N[FILENAME_MAX];

  sprintf(N, "%s.%08u.profile", prefix, id());

  FILE *F = fopen(N, "w");

  if (F == NULL)
    return;

  for (uint32 ii=0; ii<errorProfile.size(); ii++)
    fprintf(F, "%u %u %f +- %f (%u overlaps)\n",
            errorProfile[ii].bgn,        errorProfile[ii].end,
            errorProfile[ii].dev.mean(), errorProfile[ii].dev.stddev(),
            errorProfile[ii].dev.size());

  fclose(F);
}
