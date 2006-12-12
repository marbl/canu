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

#ifndef ATAC_FEATURELIST_H
#define ATAC_FEATURELIST_H

#include <string>
#include <map>
using namespace std;

#include "bio++.H"


class atacFeatureList {
public:
  atacFeatureList();
  ~atacFeatureList();

  u32bit         numberOfFeatures(void) { return(_featuresLen);  };
  atacFeature   *getFeature(u32bit i)   { return(_features + i); };
  atacFeature   *operator[](u32bit i)   { return(_features + i); };

  void           add(atacFeature &m);
  //void           delete();

  //  Sort by the A assembly or B assembly location
  void           sort(u32bit first=0, u32bit len=0);
  void           sortFeatureUID(u32bit first=0, u32bit len=0);
  void           sortParentUID(u32bit first=0, u32bit len=0);

private:
  u32bit        _featuresLen;
  u32bit        _featuresMax;
  atacFeature  *_features;
};

#endif  //  ATAC_FEATURELIST_H
