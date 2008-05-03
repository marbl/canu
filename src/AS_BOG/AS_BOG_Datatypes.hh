
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


struct BogOptions {
  static int badMateBreakThreshold;
  static bool unitigIntersectBreaking;
  static bool ejectUnhappyContained;
  static bool useGkpStoreLibStats;
};

#endif
