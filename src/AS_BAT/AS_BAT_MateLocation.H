
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

#ifndef INCLUDE_AS_BAT_MATE_LOCATION
#define INCLUDE_AS_BAT_MATE_LOCATION

static const char *rcsid_INCLUDE_AS_BAT_MATELOCATION = "$Id: AS_BAT_MateLocation.H,v 1.2 2010-12-06 08:03:48 brianwalenz Exp $";

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
