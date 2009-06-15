
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

#ifndef INCLUDE_AS_BOG_DATATYPES
#define INCLUDE_AS_BOG_DATATYPES

static const char *rcsid_INCLUDE_AS_BOG_DATATYPES = "$Id: AS_BOG_Datatypes.hh,v 1.35 2009-06-15 05:52:49 brianwalenz Exp $";

#include <map>
#include <set>
#include <list>
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"

//  Assign values to the enum to show this bug
#warning there is a comparison assuming fragment_end_type FIVE_PRIME < THREE_PRIME
enum fragment_end_type {
  FIVE_PRIME,
  THREE_PRIME
};

typedef std::list<SeqInterval> IntervalList;

class FragmentEnd {
public:
  FragmentEnd(uint32 id=0, fragment_end_type end=FIVE_PRIME) {
    _id  = id;
    _end = end;
  };

  uint32             fragId(void)  const { return(_id); };
  fragment_end_type  fragEnd(void) const { return(_end); };

  bool operator==(FragmentEnd const that) const {
    return((fragId() == that.fragId()) && (fragEnd() == that.fragEnd()));
  };

  bool operator!=(FragmentEnd const that) const {
    return((fragId() != that.fragId()) || (fragEnd() != that.fragEnd()));
  };

  bool operator<(FragmentEnd const that) const {
    if (fragId() != that.fragId())
      return fragId() < that.fragId();
    else
      return fragEnd() < that.fragEnd();
  };

private:
  uint32            _id;
  fragment_end_type _end;
};


class BestEdgeOverlap{
public:
  uint32            frag_b_id;

  float             olap_score;
  short             olap_length;

  fragment_end_type bend;

  short             ahang;
  short             bhang;

  void              print(FILE *f) {
    fprintf(f, "BestEdgeOverlap()-- id=%d score=%f length=%hd bend=%c ahang=%hd bhang=%hd\n",
            frag_b_id, olap_score, olap_length, (bend == FIVE_PRIME) ? '5' : '3', ahang, bhang);
  };
};


// Contains information on what a known fragment overlaps.
// It is assumed that an index into an array of BestOverlap
// will tell us what fragment has this best overlap
class BestFragmentOverlap{
public:
  BestEdgeOverlap five_prime;
  BestEdgeOverlap three_prime;
};


// Contains what kind of containment relationship exists between
// fragment a and fragment b

class BestContainment{
public:
  BestContainment() {
  };
  ~BestContainment() {
    delete [] olaps;
  };

  uint32  container;

  float   contain_score;

  short   a_hang;
  short   b_hang;

  bool    sameOrientation;
  bool    isContained;
  bool    isPlaced;

  bool    olapsSorted;
  short   olapsLen;
  short   olapsMax;
  uint32 *olaps;
};


class FragmentInfo {
public:
  FragmentInfo(gkStore *gkpStore) {
    gkStream         *fs = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
    gkFragment        fr;

    _numLibraries = gkpStore->gkStore_getNumLibraries();
    _numFragments = gkpStore->gkStore_getNumFragments();

    _fragLength    = new int  [_numFragments + 1];
    _mateIID       = new uint32 [_numFragments + 1];
    _libIID        = new uint32 [_numFragments + 1];

    _mean          = new double [_numLibraries + 1];
    _stddev        = new double [_numLibraries + 1];

    _numFragsInLib = new int [_numLibraries + 1];
    _numMatesInLib = new int [_numLibraries + 1];

    for (int i=0; i<_numFragments + 1; i++) {
      _fragLength[i] = 0;
      _mateIID[i] = 0;
      _libIID[i] = 0;
    }

    for (int i=0; i<_numLibraries + 1; i++) {
      _mean[i]          = 0.0;
      _stddev[i]        = 0.0;
      _numFragsInLib[i] = 0;
      _numMatesInLib[i] = 0;
    }

    for (int i=1; i<_numLibraries + 1; i++) {
      _mean[i]          = gkpStore->gkStore_getLibrary(i)->mean;
      _stddev[i]        = gkpStore->gkStore_getLibrary(i)->stddev;
      _numFragsInLib[i] = 0;
      _numMatesInLib[i] = 0;
    }

    int numDeleted = 0;
    int numLoaded  = 0;

    while(fs->next(&fr)) {
      uint32 iid = fr.gkFragment_getReadIID();
      uint32 lib = fr.gkFragment_getLibraryIID();

      if (fr.gkFragment_getIsDeleted()) {
        numDeleted++;
      } else {
        _fragLength[iid] = fr.gkFragment_getClearRegionLength();
        _mateIID[iid]    = fr.gkFragment_getMateIID();;
        _libIID[iid]     = fr.gkFragment_getLibraryIID();

        _numFragsInLib[lib]++;

        if (_mateIID[iid])
          _numMatesInLib[lib]++;

        numLoaded++;
      }
    }

    for (int i=0; i<_numLibraries + 1; i++) {
      _numMatesInLib[i] /= 2;
    }

    fprintf(stderr, "Loaded %d alive fragments, skipped %d dead fragments.\n", numLoaded, numDeleted);

    delete fs;
  };
  ~FragmentInfo() {
    delete [] _fragLength;
    delete [] _mateIID;
    delete [] _libIID;
  };

  int     numFragments(void) { return(_numFragments); };
  int     numLibraries(void) { return(_numLibraries); };

  int     fragmentLength(uint32 iid) { return(_fragLength[iid]); };
  uint32  mateIID(uint32 iid)        { return(_mateIID[iid]); };
  uint32  libraryIID(uint32 iid)     { return(_libIID[iid]);  };

  double  mean(uint32 iid)   { return(_mean[iid]); };
  double  stddev(uint32 iid) { return(_stddev[iid]); };

  int     numMatesInLib(uint32 iid) { return(_numMatesInLib[iid]); };

private:
  int      _numFragments;
  int      _numLibraries;

  int     *_fragLength;
  uint32  *_mateIID;
  uint32  *_libIID;

  double  *_mean;
  double  *_stddev;

  int     *_numFragsInLib;
  int     *_numMatesInLib;
};

#endif

