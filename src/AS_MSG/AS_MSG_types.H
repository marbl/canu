
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
/* $Id: AS_MSG_types.h,v 1.6 2012-08-01 03:20:19 brianwalenz Exp $   */

#ifndef AS_MSG_PMESG_TYPES_H
#define AS_MSG_PMESG_TYPES_H

static const char *rcsid_AS_MSG_PMESG_TYPES_H = "$Id: AS_MSG_types.h,v 1.6 2012-08-01 03:20:19 brianwalenz Exp $";

#include <cstdio>
#include <ctime>

using namespace std;

#include "AS_global.H"
#include "AS_UTL_IID.H"

typedef AS_IID      IntLibrary_ID;
typedef AS_IID      IntDist_ID;
typedef AS_IID      IntFragment_ID;
typedef AS_IID      IntChunk_ID;
typedef AS_IID      IntUnitig_ID;
typedef AS_IID      IntContig_ID;
typedef AS_IID      IntScaffold_ID;




class LinkType {
public:
  LinkType()       { linkType = 'U'; };
  ~LinkType()      {};

  //  This operator is dangerous.  In a printf, where we'd usually want to have an implicit cast,
  //  the compiler warns that it doesn't cast, and will crash at runtime.  Use the toLetter() method
  //  instead.
  //
  //operator char()  { return(toLetter()); };

  char  toLetter(void) const {
    return(linkType);
  };

  bool  isMatePair(void) const    { assert(linkType != 0);  return(linkType == 'M'); };
  bool  isOverlap(void) const     { assert(linkType != 0);  return(linkType == 'X'); };  //  Really should be 'O', but that changes the asm

  void  setIsMatePair(void)  { linkType = 'M'; };
  void  setIsOverlap(void)   { linkType = 'X'; };

  bool  isValid(void) const {
    return(linkType != 0);
  };

private:
  char  linkType;
  char  pad[3];
};



//  PairOrient -- the orientation between two sequence objects.
//
//  Former names:
//    OrientType, ChunkOrientationType
//
//  Former values:
//    AS_UNKNOWN == XX_XX == 'U'
//    AS_INNIE   == AB_BA == 'I'
//    AS_OUTTIE  == BA_AB == 'O'
//    AS_NORMAL  == AB_AB == 'N'
//    AS_ANTI    == BA_BA == 'A'

class PairOrient {
public:
  PairOrient()    { orient = 'U'; pad[0] = pad[1] = pad[2] = 0; };
  ~PairOrient()   {};

  char  toLetter(void) const  { return(orient); };

  bool  isInnie(void) const   { assert(orient != 0);  return(orient == 'I'); };
  bool  isOuttie(void) const  { assert(orient != 0);  return(orient == 'O'); };
  bool  isNormal(void) const  { assert(orient != 0);  return(orient == 'N'); };
  bool  isAnti(void) const    { assert(orient != 0);  return(orient == 'A'); };
  bool  isUnknown(void) const { assert(orient != 0);  return(orient == 'U'); };

  void  setIsInnie(void)      { orient = 'I'; };
  void  setIsOuttie(void)     { orient = 'O'; };
  void  setIsNormal(void)     { orient = 'N'; };
  void  setIsAnti(void)       { orient = 'A'; };
  void  setIsUnknown(void)    { orient = 'U'; };

  bool  isAB_BA(void) const   { assert(orient != 0);  return(orient == 'I'); };
  bool  isBA_AB(void) const   { assert(orient != 0);  return(orient == 'O'); };
  bool  isAB_AB(void) const   { assert(orient != 0);  return(orient == 'N'); };
  bool  isBA_BA(void) const   { assert(orient != 0);  return(orient == 'A'); };

  void  setIsAB_BA(void)      { orient = 'I'; };
  void  setIsBA_AB(void)      { orient = 'O'; };
  void  setIsAB_AB(void)      { orient = 'N'; };
  void  setIsBA_BA(void)      { orient = 'A'; };

  void  invert(void) {
    assert(orient != 0);  

    switch (orient) {
      case 'I': orient = 'O'; break;
      case 'O': orient = 'I'; break;
      case 'N': orient = 'A'; break;
      case 'A': orient = 'N'; break;
    }
  };
  void  flip(void) {
    switch(orient){
      case 'N': orient = 'A'; break;
      case 'A': orient = 'N'; break;
    }
  };
  void  swap(void) {
    // Given orientation relating X and Y as ori(X)_ori(Y), return ori(Y)_ori(X)
    switch(orient){
      case 'I': orient = 'O'; break;
      case 'O': orient = 'I'; break;
    }
  };

  bool  operator==(const PairOrient &that) const {
    assert(orient      != 0);
    assert(that.orient != 0);
    return(orient == that.orient);
  };
  bool  operator!=(const PairOrient &that) const {
    assert(orient      != 0);
    assert(that.orient != 0);
    return(orient != that.orient);
  };

  bool  isValid(void) const {
    return(orient != 0);
  };

private:
  char  orient;
  char  pad[3];
};

//  SequenceOrient -- the orientation of a single sequence object.
//
//  Former names:
//    DirectionType, FragOrient, ChunkOrient, NodeOrient and CIOrient,
//
//  Former values:
//    AS_FORWARD == A_B == 'F'
//    AS_REVERSE == B_A == 'R'

class SequenceOrient {
public:
  SequenceOrient()    { orient = 'U'; };
  ~SequenceOrient()   {};

  char  toLetter(void) const      { return(orient); };

  bool  isForward(void) const     { assert(orient != 0);  return(orient == 'F'); };
  bool  isReverse(void) const     { assert(orient != 0);  return(orient == 'R'); };
  bool  isUnknown(void) const     { assert(orient != 0);  return(orient == 'U'); };

  void  setIsForward(bool f=true)   { orient = (f) ? 'F' : 'R'; };
  void  setIsReverse(bool r=true)   { orient = (r) ? 'R' : 'F'; };
  void  setIsUnknown(void)          { orient =       'U';       };

  void  flip(void) {
    assert(orient != 0);  

    switch (orient) {
      case 'F': orient = 'R'; break;
      case 'R': orient = 'F'; break;
    }
  };

  bool  operator==(const SequenceOrient &that) const {
    assert(orient      != 0);
    assert(that.orient != 0);
    return(orient == that.orient);
  };
  bool  operator!=(const SequenceOrient &that) const {
    assert(orient      != 0);
    assert(that.orient != 0);
    return(orient != that.orient);
  };

  bool  isValid(void) const {
    return(orient != 0);
  };

private:
  char  orient;
};






//
//  XXX When more types are added, or when the library starts telling
//  what type the reads are, it will be VERY useful to remove ALL of
//  these types to catch the numerous places where the type of the
//  read is assumed to be AS_READ.
//
//  After wasting a couple hours, BPW figured out why AS_UNITIG and
//  AS_CONTIG are FragTypes.  Consensus is using the IMP to store both
//  fragments and unitigs.
//

typedef enum {
  AS_READ    = (int)'R',  //  Celera Read
  AS_EXTR    = (int)'X',  //  External WGS read
  AS_TRNR    = (int)'T',  //  Transposon library read
  AS_UNITIG  = (int)'U',  //  Assembled unitig
  AS_CONTIG  = (int)'C'   //  Assembled contig
} FragType;

#define AS_FA_READ(type) 		((type == AS_READ) || (type == AS_EXTR))
#define AS_FA_RANDOM(type) 		((type == AS_READ) || (type == AS_EXTR))
#define AS_FA_SHREDDED(type) 		(0)
#define AS_FA_CAN_HAVE_MATE(type) 	((type == AS_READ) || (type == AS_EXTR) || (type == AS_TRNR))
#define AS_FA_GUIDE(type)         	(0)
#define AS_FA_TRUSTED_FOR_CNS(type) 	((type == AS_READ) || (type == AS_TRNR))
#define AS_FA_EXTERNAL_READ(type)	((type == AS_EXTR))

typedef enum {
  AS_UNIQUE_UNITIG   = (int)'U',  // U-Unique
  AS_ROCK_UNITIG     = (int)'R',  // Rock
  AS_STONE_UNITIG    = (int)'S',  // Stone
  AS_PEBBLE_UNITIG   = (int)'P',  // Pebble
  AS_SINGLE_UNITIG   = (int)'s',  // Singleton Unitig in Unplaced Contig
  AS_OTHER_UNITIG    = (int)'X'   // Unspecified surrogate unitig
} UnitigType;

/*OVL message*/

#define AS_LONGEST_DELTA     (126)
#define AS_LONG_DELTA_CODE  (-127)
#define AS_POLY_DELTA_CODE  (-128)
#define AS_ENDOF_DELTA_CODE    (0)

typedef enum {
  AS_DOVETAIL    = (int)'D',
  AS_CONTAINMENT = (int)'C',

//  These are dead.  BPW 2009-06

  AS_SUPERREPEAT = (int)'S',

  AS_DOVETAIL_TRANS = (int) 'X',
  // A dovetail overlap that is transitively inferrable by a path in
  // the fragment overlap graph, for example using a path of dovetail
  // overlaps followed by a path of containment overlaps.

  /* The following are created during de-chording the fragment overlap
     graph. They are not by default present in a path compressed
     fragment overlap graph but instead removed since they are easily
     inferred. */
  AS_DOVETAIL_CHORD  = (int)'d',
  // This dovetail overlap is known to be a chord in the dovetail
  // overlap sub-graph.
  AS_CONTAINMENT_CHORD = (int)'c'
  // This containment overlap is known to be a chord in the
  // containment overlap sub-graph.

} OverlapType;

typedef enum {
  /* The following are CGW overlap output classifications: */
  AS_NO_OVERLAP	    = (int) 'N', // Not used by the unitigger.
  AS_OVERLAP        = (int) 'O', // A dovetail overlap between unitigs.
  AS_TANDEM_OVERLAP = (int) 'T', // Not used by the unitigger.

  /* The following are unitigger overlap output classifications: */
  AS_1_CONTAINS_2_OVERLAP = (int) 'C',
  AS_2_CONTAINS_1_OVERLAP = (int) 'I',
  // Two types of containment overlap types.

  AS_TOUCHES_CONTAINED_OVERLAP   = (int) 'M',
  // A dovetail overlap to a singleton unitig composed of a contained
  // fragment from a non-contained unitig. This overlaps touch
  // multiply contained fragments and orphan contained fragments. This
  // is another kind of overlap necessary when processing containment
  // overlaps.

  AS_TRANSCHUNK_OVERLAP   = (int) 'X',
  // A dovetail overlap between unitigs that was marked for removal by
  // because it was transitively inferrable in the fragment overlap
  // graph but was not a chord in the dovetail overlap sub-graph.

  AS_DOVETAIL_CHORD_OVERLAP  = (int)'d',
  // This dovetail overlap is known to be a chord in the dovetail
  // overlap sub-graph.  This clasification is created during
  // de-chording the fragment overlap graph. They are by default not
  // present in a path compressed fragment overlap graph but instead
  // removed since they are easily inferred.

  AS_1_CONTAINS_2_CHORD_OVERLAP = (int)'c',
  // This containment overlap is known to be a chord in the
  // containment overlap sub-graph.  This clasification is created
  // during de-chording the fragment overlap graph. They are by
  // default not present in a path compressed fragment overlap graph
  // but instead removed since they are easily inferred.

  AS_2_CONTAINS_1_CHORD_OVERLAP = (int)'i',
  // This containment overlap is known to be a chord in the
  // containment overlap sub-graph.  This clasification is created
  // during de-chording the fragment overlap graph. They are by
  // default not present in a path compressed fragment overlap graph
  // but instead removed since they are easily inferred.

  AS_BETWEEN_CONTAINED_OVERLAP   = (int) 'Y',
  // A dovetail overlap between unitigs each spanned by a contained
  // fragment. Thus, this overlap is between multiply contained
  // fragments and/or orphan contained fragments.

  AS_1_CONTAINS_2_STACK_OVERLAP   = (int) 'Z'
  // A containment overlap between two globally contained fragments
  // that is not a chord in the containment overlap graph.

} UnitigOverlapType;


typedef enum {
  AS_UNIQUE =     (int)'U',
  AS_NOTREZ =     (int)'N',
  AS_SEP =        (int)'S',
  AS_UNASSIGNED = (int)'X'
} UnitigStatus;

typedef enum {
  AS_FORCED_NONE    = (int)'X',
  AS_FORCED_UNIQUE  = (int)'U',
  AS_FORCED_REPEAT  = (int)'R'
} UnitigFUR;

/* IMP message */

typedef enum {
  AS_PLACED	= (int)'P',
  AS_UNPLACED   = (int)'U'
} ContigStatus;

/* ICM */


typedef enum {
  AS_IN_ASSEMBLY	= (int) 'A',
  AS_POLYMORPHIRSM	= (int) 'P',
  AS_BAD		= (int) 'B',
  AS_UNKNOWN_IN_ASSEMBLY= (int) 'U'
} PlacementStatusType;


/* IAF message */

/* Grangers new list, roughly in order of precedence */
typedef enum {
  UNASSIGNED_MATE    = 'Z',
  GOOD_MATE          = 'G',
  BAD_SHORT_MATE     = 'C',
  BAD_LONG_MATE      = 'L',
  SAME_ORIENT_MATE   = 'S',
  OUTTIE_ORIENT_MATE = 'O',
  NO_MATE            = 'N',
  BOTH_CHAFF_MATE    = 'H',
  CHAFF_MATE         = 'A',
  BOTH_DEGEN_MATE    = 'D',
  DEGEN_MATE         = 'E',
  BOTH_SURR_MATE     = 'U',
  SURR_MATE          = 'R',
  DIFF_SCAFF_MATE    = 'F'
} MateStatType;


#endif
