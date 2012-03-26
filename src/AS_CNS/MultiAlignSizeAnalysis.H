
/**************************************************************************
 * Copyright (C) 2012, J Craig Venter Institute. All rights reserved.
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

#ifndef MULTIALIGNSIZEANALYSIS
#define MULTIALIGNSIZEANALYSIS

static const char *rcsid_MULTIALIGNSIZEANALYSIS = "$Id: MultiAlignSizeAnalysis.H,v 1.1 2012-03-26 07:58:27 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"

#include <vector>
#include <algorithm>
using namespace std;

class sizeAnalysis {
public:
  sizeAnalysis(uint64 genomeSize);
  ~sizeAnalysis();

  void         evaluateTig(MultiAlignT *ma, bool isUnitig);
  void         finalize(void);

  void         printSummary(FILE *out, char *description, vector<uint32> &data);
  void         printSummary(FILE *out);

private:
  uint64          genomeSize;

  vector<uint32>  utgLenUnassigned;
  vector<uint32>  utgLenUnique;
  vector<uint32>  utgLenSep;
  vector<uint32>  utgLenNotRez;

  vector<uint32>  ctgLenPlaced;
  vector<uint32>  ctgLenUnplaced;
};


#endif  //  MULTIALIGNSIZEANALYSIS
