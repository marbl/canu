/**************************************************************************
 * This file is part of A2Amapper.
 * Copyright (c) 2004 Applera Corporation
 * Author: Clark Mobarry
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
**************************************************************************/


#ifndef GF_ALN_LOCALEXTEND_H
#define GF_ALN_LOCALEXTEND_H

#include "GF_ALN_local.h"

// Extend a local alignment out from seed coordinates x and y,
// return the trace of the alignment
// and, as side effects,
// 1) put the aligned substrings in Asub and Bsub
// 2) return the segment coordinates in seg
// 3) set spnt to the ahang of the Asub to Bsub alignment
int *extend_from_seed(const string &A,
		      const string &B,
		      int x,
		      int y,
		      float MaxDiffs,
		      int MaxiGap,
		      string &Asub,
		      string &Bsub,
		      Local_Segment *&seg,
		      int & spnt);


#endif
