
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

#include <map>
#include <list>
#include <vector>

extern "C" {
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"
}

enum fragment_end_type {
  FIVE_PRIME, 	// 5' End of fragment
  THREE_PRIME 	// 3' End of Fragment
};

typedef AS_IID    iuid;
const iuid NULL_FRAG_ID=0;

typedef std::list<SeqInterval> IntervalList;

class FragmentEnd {
public:
  FragmentEnd(iuid id=0, fragment_end_type end=FIVE_PRIME) {
    _id  = id;
    _end = end;
  };

  iuid               fragId(void)  const { return(_id); };
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
  iuid              _id;
  fragment_end_type _end;
};


class BestEdgeOverlap{
public:
  iuid              frag_b_id;

  float             olap_score;
  short             olap_length;

  fragment_end_type bend;                

  short             ahang;
  short             bhang;
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

  iuid    container;

  float   contain_score;

  short   a_hang;
  short   b_hang;

  bool    sameOrientation;
  bool    isContained;
  bool    isPlaced;

  bool    olapsSorted;
  short   olapsLen;
  short   olapsMax;
  iuid   *olaps;
};


class FragmentInfo {
public:
  FragmentInfo(GateKeeperStore *gkpStore) {
    FragStream       *fs = openFragStream(gkpStore, FRAG_S_INF);
    fragRecord        fr;

    _numLibraries = getNumGateKeeperLibraries(gkpStore);
    _numFragments = getNumGateKeeperFragments(gkpStore);

    _fragLength    = new int  [_numFragments + 1];
    _mateIID       = new iuid [_numFragments + 1];
    _libIID        = new iuid [_numFragments + 1];

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
      _mean[i]          = getGateKeeperLibrary(gkpStore, i)->mean;
      _stddev[i]        = getGateKeeperLibrary(gkpStore, i)->stddev;
      _numFragsInLib[i] = 0;
      _numMatesInLib[i] = 0;
    }

    int numDeleted = 0;
    int numLoaded  = 0;

    while(nextFragStream(fs, &fr)) {
      iuid  iid = getFragRecordIID(&fr);
      iuid  lib = getFragRecordLibraryIID(&fr);

      if (getFragRecordIsDeleted(&fr)) {
        numDeleted++;
      } else {
        _fragLength[iid] = (getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_OBT) -
                            getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_OBT));
        _mateIID[iid]    =  getFragRecordMateIID(&fr);;
        _libIID[iid]     =  getFragRecordLibraryIID(&fr);

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

    closeFragStream(fs); 
  };
  ~FragmentInfo() {
    delete [] _fragLength;
    delete [] _mateIID;
    delete [] _libIID;
  };

  int     numFragments(void) { return(_numFragments); };
  int     numLibraries(void) { return(_numLibraries); };

  int     fragmentLength(iuid iid) { return(_fragLength[iid]); };
  iuid    mateIID(iuid iid)        { return(_mateIID[iid]); };
  iuid    libraryIID(iuid iid)     { return(_libIID[iid]);  };

  double  mean(iuid iid)   { return(_mean[iid]); };
  double  stddev(iuid iid) { return(_stddev[iid]); };

  int     numMatesInLib(iuid iid) { return(_numMatesInLib[iid]); };

private:
  int      _numFragments;
  int      _numLibraries;

  int     *_fragLength;
  iuid    *_mateIID;
  iuid    *_libIID;

  double  *_mean;
  double  *_stddev;

  int     *_numFragsInLib;
  int     *_numMatesInLib;
};

#endif

