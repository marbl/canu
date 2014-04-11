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

  uint32         numberOfFeatures(void) { return(_featuresLen);  };
  atacFeature   *getFeature(uint32 i)   { return(_features + i); };
  atacFeature   *operator[](uint32 i)   { return(_features + i); };

  void           add(atacFeature &m);
  //void           delete();

  //  Sort by the A assembly or B assembly location
  void           sort(uint32 first=0, uint32 len=0);
  void           sortFeatureUID(uint32 first=0, uint32 len=0);
  void           sortParentUID(uint32 first=0, uint32 len=0);

private:
  uint32        _featuresLen;
  uint32        _featuresMax;
  atacFeature  *_features;
};

#endif  //  ATAC_FEATURELIST_H
