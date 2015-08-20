
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
 *    src/AS_BAT/AS_BAT_MateLocation.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz beginning on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_MATE_LOCATION
#define INCLUDE_AS_BAT_MATE_LOCATION

static const char *rcsid_INCLUDE_AS_BAT_MATELOCATION = "$Id$";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_InsertSizes.H"


class MateLocationEntry {
public:
  SeqInterval mlePos1;
  SeqInterval mlePos2;
  uint32      mleFrgID1;
  uint32      mleFrgID2;
  uint32      mleUtgID1;
  uint32      mleUtgID2; // in the future the table might be across unitigs
  bool        isGrumpy;

  //bool operator==(MateLocation &that) {
  //  return((mlePos1 == that.mlePos1) && (mlePos2 == that.mlePos2));
  //};

  bool operator<(MateLocationEntry const &that) const {

#if 0
    if (mlePos1  < that.mlePos1)                           return true;
    if (mlePos1 == that.mlePos1 && mlePos2 < that.mlePos2) return true;
    else                                                   return false;
#endif

    if (mleFrgID1 == 0)
      return(true);
    if (that.mleFrgID1 == 0)
      return(false);

    int32  m1 = MIN(     mlePos1.bgn,      mlePos1.end);
    int32  t1 = MIN(that.mlePos1.bgn, that.mlePos1.end);
    int32  m2 = MIN(     mlePos2.bgn,      mlePos2.end);
    int32  t2 = MIN(that.mlePos2.bgn, that.mlePos2.end);

    return((m1 < t1) || ((m1 == t1) && (m2 < t2)));
  };
};



//  The MateLocation table builds a table of positions of mated reads.
//    o  Unmated reads are NOT in the table.
//    o  Mates in other unitigs are not in the table.  The fragment
//       in this unitig is present, but the mate is NULL.
//    o  Mates in the same unitig are in the table.
//
class MateLocation {
public:
  MateLocation(UnitigVector &unitigs, Unitig *utg);
  ~MateLocation();

  MateLocationEntry getById(uint32 fragId) {
    map<uint32,uint32>::const_iterator  e = _iidToTableEntry.find(fragId);

    if (e == _iidToTableEntry.end())
      return(_table[0]);
    else
      return(_table[e->second]);
  };

  uint32             numMates(void) {
    return(_numMates);
  };

  void               dumpHappiness(const char *prefix, const char *name);
  void               dumpStats(void);

  int32                     nunmated;

  int32  *good,             ngood[3];
  int32  *badFwd,           nbadFwd[3];
  int32  *badRev,           nbadRev[3];

  int32  *goodExternalFwd,  ngoodExternalFwd[3];
  int32  *goodExternalRev,  ngoodExternalRev[3];

  int32  *badExternalFwd,   nbadExternalFwd[3];
  int32  *badExternalRev,   nbadExternalRev[3];

  int32  *badCompressed,    nbadCompressed[3];
  int32  *badStretched,     nbadStretched[3];
  int32  *badNormal,        nbadNormal[3];
  int32  *badAnti,          nbadAnti[3];
  int32  *badOuttie,        nbadOuttie[3];

private:
  void buildTable(void);
  void buildHappinessGraphs(UnitigVector &unitigs);

  void incrRange(int32 *graph, int32 val, int32 n, int32 m) {
    n = MAX(n, 0);
    m = MIN(m, _tigLen);

    assert(n <= _tigLen);
    assert(0 <= m);

    //  Earlier versions asserted n<m (and even earlier versions used i<=m in the loop below, which
    //  made this far more complicated than necessary).  Now, we don't care.  We'll adjust n and m
    //  to the min/max possible, and ignore out of bounds cases.  Those happen when, for example,
    //  fragments are the same orientation.  If one of those is the last fragment in the unitig,
    //  we'll call incrRange with n=(the higher coord)=(_tigLen), and m=(the lower coord + max
    //  insert size).  We threshold m to _tigLen, and correctly do nothing in the loop.

    for (int32 i=n; i<m; i++)
      graph[i] += val;
    for (int32 i=m; i<n; i++)
      graph[i] += val;
  };

  Unitig                    *_tig;
  int32                      _tigLen;

  uint32                     _numMates;

  vector<MateLocationEntry>  _table;
  map<uint32,uint32>         _iidToTableEntry;
};


#endif // INCLUDE_AS_BAT_MATE_LOCATION
