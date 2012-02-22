
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

#ifndef INCLUDE_AS_BAT_MERGESPLITJOIN
#define INCLUDE_AS_BAT_MERGESPLITJOIN

static const char *rcsid_INCLUDE_AS_BAT_MERGESPLITJOIN = "$Id: AS_BAT_MergeSplitJoin.H,v 1.4 2012-02-22 19:16:36 brianwalenz Exp $";

void
mergeSplitJoin(UnitigVector &unitigs, const char *prefix, bool shatterRepeats);

void
extendByMates(UnitigVector &unitigs,
              double        erateGraph,
              double        elimitGraph);

void
reconstructRepeats(UnitigVector &unitigs,
                   double        erateGraph,
                   double        elimitGraph);

void
promoteToSingleton(UnitigVector &unitigs, bool enablePromoteToSingleton);


#endif  //  INCLUDE_AS_BAT_MERGESPLITJOIN
