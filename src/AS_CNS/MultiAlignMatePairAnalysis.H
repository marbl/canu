
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

static const char *rcsid_MULTIALIGNMATEPAIRANALYSIS = "$Id: MultiAlignMatePairAnalysis.H,v 1.4 2012-06-28 20:01:18 skoren Exp $";

#include "AS_global.H"
#include "AS_PER_gkpStore.H"
#include "MultiAlign.H"
#include "MultiAlignStore.H"
#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"

#include <pthread.h>

class mpaLibrary;


class matePairAnalysis {
public:
  matePairAnalysis(char *gkpStore);
  ~matePairAnalysis();

  void         evaluateTig(MultiAlignT *ma);
  void         evaluateTig(MultiAlignT *ma, pthread_mutex_t *mux);
  void         finalize(void);

  void         printSummary(FILE *out);
  void         writeUpdate(char *prefix);
  void         drawPlots(char *prefix);

  double       mean(uint32 libid);
  double       stddev(uint32 libid);
  uint64       numSamples(uint32 libid);

private:
  gkStore     *gkpStore;
  mpaLibrary  *libdata;

  AS_IID      *pairing;
  AS_IID      *library;
};


#endif  //  MULTIALIGNMATEPAIRANALYSIS
