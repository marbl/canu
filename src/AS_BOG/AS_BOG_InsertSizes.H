
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

#ifndef INCLUDE_AS_INSERT_SIZES
#define INCLUDE_AS_INSERT_SIZES

static const char *rcsid_INCLUDE_AS_BOG_INSERTSIZES = "$Id$";

#include "AS_BOG_Datatypes.H"

class Unitig;

class
InsertSizes {
public:
  InsertSizes();
  ~InsertSizes();

  int32      mean(uint32 libIID)    { assert(libIID > 0); assert(libIID <= _numLibs); return(_mean[libIID]);          };
  int32      stddev(uint32 libIID)  { assert(libIID > 0); assert(libIID <= _numLibs); return(_stddev[libIID]);        };
  bool       valid(uint32 libIID)   { assert(libIID > 0); assert(libIID <= _numLibs); return(_samples[libIID] >= 10); };

private:
  void       accumulateLibraryStats(Unitig *utg);

  uint32    _numLibs;

  int32   **_dist;
  int32    *_distLen;
  int32    *_distMax;

  int32    *_mean;
  int32    *_stddev;
  int32    *_samples;
};

#endif // INCLUDE_AS_INSERT_SIZES
