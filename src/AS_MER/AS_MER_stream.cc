
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
#include "AS_MER_stream.h"
#include <stdlib.h>

merStream::merStream(cds_uint32 merSize, char *gkpstore) {
  _theSeq          = 0L;
  _theLen          = 0;
  _thePos          = 0;
  _endPos          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;
  _theMerMask      = CDS_UINT64_MASK(_merSize << 1);

  _theFMer         = 0;
  _theRMer         = 0;

  _theRMerShift    = (_merSize << 1) - 2;

  _fs = openGateKeeperStore(gkpstore, FALSE);
  if (_fs == NULL) {
    fprintf(stderr, "ERROR:  couldn't open the gatekeeper store '%s'\n", gkpstore);
    exit(1);
  }
  _fr  = new_fragRecord();

  _iid = 1;
  _max = getLastElemFragStore(_fs);

  _skipNum = 1;
  loadMer(_merSize - 1);
}


merStream::merStream(cds_uint32 merSize, char *gkpstore, int skipNum) {
  _theSeq          = 0L;
  _theLen          = 0;
  _thePos          = 0;
  _endPos          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;
  _theMerMask      = CDS_UINT64_MASK(_merSize << 1);

  _theFMer         = 0;
  _theRMer         = 0;

  _theRMerShift    = (_merSize << 1) - 2;

  _fs = openGateKeeperStore(gkpstore, FALSE);
  if (_fs == NULL) {
    fprintf(stderr, "ERROR:  couldn't open the gatekeeper store '%s'\n", gkpstore);
    exit(1);
  }
  _fr  = new_fragRecord();

  _iid = 1;
  _max = getLastElemFragStore(_fs);

  _skipNum=skipNum;
  loadMer(_merSize - 1);
}

merStream::~merStream() {
  del_fragRecord(_fr);
  closeGateKeeperStore(_fs);
}
