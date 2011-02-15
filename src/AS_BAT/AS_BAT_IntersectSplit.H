
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

#ifndef INCLUDE_AS_BAT_INTERSECTSPLIT
#define INCLUDE_AS_BAT_INTERSECTSPLIT

static const char *rcsid_INCLUDE_AS_BAT_INTERSECTSPLIT = "$Id: AS_BAT_IntersectSplit.H,v 1.2 2011-02-15 08:10:11 brianwalenz Exp $";




class intersectionEvidence {
public:
  intersectionEvidence() {
    frag5tig = 0;
    frag3tig = 0;

    frag5confirmed = false;
    frag3confirmed = false;

    frag5self = false;
    frag3self = false;
  };
  ~intersectionEvidence() {
  };


  BestEdgeOverlap   edge5;    //  fragID of the frag on our 5' best overlap
  BestEdgeOverlap   edge3;    //

  uint32   frag5tig;          //  tigID the frag on our 5' best overlap is in
  uint32   frag3tig;          //

  uint32   frag5confirmed:1;  //  true if our 5' end is confirmed by a best overlap in the same unitig
  uint32   frag3confirmed:1;  //

  uint32   frag5self:1;       //  true if our 5' end is intersecting the same unitig
  uint32   frag3self:1;
};




class intersectionPoint {
public:
  intersectionPoint() {
    isectFrg  = 0;
    isect3p   = false;

    invadFrg  = 0;
    invad3p   = false;

    isSelf    = false;
  };
  intersectionPoint(BestEdgeOverlap edge, uint32 sId, bool s3p, bool self) {
    isectFrg  = edge.fragId();
    isect3p   = edge.frag3p();

    invadFrg = sId;
    invad3p  = s3p;

    isSelf    = self;

    //fprintf(stderr, "intersectionPoint()-- %d/%c' from %d/%c' self=%d\n",
    //        isectFrg,  isect3p  ? '3' : '5',
    //        invadFrg, invad3p ? '3' : '5',
    //        isSelf);
  };
  ~intersectionPoint() {
  };

  bool operator<(const intersectionPoint that) const {
    return(isectFrg < that.isectFrg);
  };

  uint32   isectFrg;  //  Fragment that is intersected into, we split on this.
  bool     isect3p;   //  True if we intersected onto the 3' end of the fragment.

  uint32   invadFrg;
  bool     invad3p;

  bool     isSelf;
};




class intersectionList {
public:
  intersectionList();
  intersectionList(UnitigVector &unitigs);
  ~intersectionList();

  void               logIntersections(void);


  uint32             numIntersections(uint32  fid) {
    if (isectsNum.find(fid) == isectsNum.end())
      return(0);
    return(isectsNum[fid]);
  };

  intersectionPoint *getIntersection(uint32  fid, uint32 i) {

    if (isectsMap.find(fid) == isectsMap.end())
      return(NULL);

    return(&isects[isectsMap[fid] + i]);
  };

private:
  vector<intersectionPoint>   isects;
  map<uint32,uint32>          isectsNum;  //  Number of intersections for this fragment
  map<uint32,uint32>          isectsMap;  //  Start index in isects for this fragment
};

void breakUnitigs(UnitigVector &unitigs, char *output_prefix, bool enableIntersectionBreaking);

#endif  //  INCLUDE_AS_BAT_INTERSECTSPLIT
