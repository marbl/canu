
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
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

#ifndef MULTIALIGNMATEPAIRANALYSIS
#define MULTIALIGNMATEPAIRANALYSIS

static const char *rcsid_MULTIALIGNMATEPAIRANALYSIS = "$Id: MultiAlignMatePairAnalysis.H,v 1.2 2011-12-18 07:04:08 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"

class mpaLibrary;


class matePairAnalysis {
public:
  matePairAnalysis(char *gkpStore);
  ~matePairAnalysis();

  void         evaluateTig(MultiAlignT *ma);
  void         finalize(void);

  void         printSummary(FILE *out);
  void         drawPlots(char *prefix);

  double       mean(uint32 libid);
  double       stddev(uint32 libid);

private:
  gkStore     *gkpStore;
  mpaLibrary  *libdata;

  AS_IID      *pairing;
  AS_IID      *library;
};


#endif  //  MULTIALIGNMATEPAIRANALYSIS
