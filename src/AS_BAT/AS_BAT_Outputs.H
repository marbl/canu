
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

#ifndef INCLUDE_AS_BAT_OUTPUTS
#define INCLUDE_AS_BAT_OUTPUTS

static const char *rcsid_INCLUDE_AS_BAT_OUTPUTS = "$Id: AS_BAT_Outputs.H,v 1.3 2012-01-15 23:49:34 brianwalenz Exp $";

#include "MultiAlignStore.h"

void unitigToMA(MultiAlignT *ma, uint32 iumiid, Unitig *utg);

void writeIUMtoFile(UnitigVector &unitigs, char *fileprefix, char *tigStorePath, uint32 fragment_count_target, bool isFinal=true);
void writeOVLtoFile(UnitigVector &unitigs, char *fileprefix);

#endif  //  INCLUDE_AS_BAT_OUTPUTS
