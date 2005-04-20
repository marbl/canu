// This file is part of A2Amapper.
// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef OVERLAP_STATS_H
#define OVERLAP_STATS_H

#include "overlap.H"

//  Statistics and Histograms
//
//  Index 1 is the assembly, 2 is the mapping.  Stats count the
//  number of bases covered, histograms are of the block sizes.
//
class histogram {
public:
  histogram(u32bit max=1048576) {
    histMax = max;
    hist    = new u32bit [histMax];
    sum     = 0;

    for (u32bit i=0; i<histMax; i++)
      hist[i] = 0;
  };
  ~histogram() {
    delete [] hist;
  };

  void  update(u32bit length) {
    sum += length;
    hist[ ((length) > histMax) ? 0 : length ]++;
  };

  void operator+=(u32bit length) {
    update(length);
  };

  u32bit  *hist;
  u32bit   histMax;
  u32bit   sum;
};

class overlapStats {
public:
  histogram  unmapped;
  histogram  unique[2];
  histogram  different;
  histogram  wilddiff;
  histogram  same;
  histogram  inconsistent;
};


#endif  //  OVERLAP_STATS_H
