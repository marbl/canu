
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
 *    src/AS_BAT/AS_BAT_IntersectSplit.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-DEC-06 to 2013-AUG-01
 *      are Copyright 2010-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_INTERSECTSPLIT
#define INCLUDE_AS_BAT_INTERSECTSPLIT




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
