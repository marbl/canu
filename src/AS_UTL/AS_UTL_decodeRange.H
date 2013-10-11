
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
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

#ifndef AS_UTL_DECODERANGE_H
#define AS_UTL_DECODERANGE_H

#include "AS_global.H"

#include <set>

using namespace std;

void  AS_UTL_decodeRange(char *range, set<uint64> &ranges);
void  AS_UTL_decodeRange(char *range, set<uint32> &ranges);

void  AS_UTL_decodeRange(char *range, uint64 &lo, uint64 &hi);
void  AS_UTL_decodeRange(char *range,  int64 &lo,  int64 &hi);
void  AS_UTL_decodeRange(char *range, uint32 &lo, uint32 &hi);
void  AS_UTL_decodeRange(char *range,  int32 &lo,  int32 &hi);

#endif  //  AS_UTL_DECODERANGE_H
