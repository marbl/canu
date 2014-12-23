
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

#ifndef TGTIGSIZEANALYSIS
#define TGTIGSIZEANALYSIS

static const char *rcsid_TGTIGSIZEANALYSIS = "$Id$";

#include "tgTig.H"

#include <vector>

using namespace std;


class tgTigSizeAnalysis {
public:
  tgTigSizeAnalysis(uint64 genomeSize);
  ~tgTigSizeAnalysis();

  void         evaluateTig(tgTig *tig);
  void         finalize(void);

  void         printSummary(FILE *out, char *description, vector<uint32> &data);
  void         printSummary(FILE *out);

private:
  uint64          genomeSize;

  vector<uint32>  lenSuggestRepeat;
  vector<uint32>  lenSuggestUnique;
  vector<uint32>  lenSuggestCircular;
  vector<uint32>  lenSuggestHaploid;

  vector<uint32>  lenSingleton;  //  Has exactly one child
  vector<uint32>  lenAssembled;  //  Has more than one child
};


#endif  //  TGTIGSIZEANALYSIS
