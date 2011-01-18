
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

#ifndef INCLUDE_AS_BAT_PLACEZOMBIES
#define INCLUDE_AS_BAT_PLACEZOMBIES

static const char *rcsid_INCLUDE_AS_BAT_PLACEZOMBIES = "$Id: AS_BAT_PlaceZombies.H,v 1.2 2011-01-18 22:37:22 brianwalenz Exp $";

void
placeZombies(UnitigVector &unitigs,
             OverlapStore *ovlStoreUniq, OverlapStore *ovlStoreRept, double erate, double elimit);

#endif  //  INCLUDE_AS_BAT_PLACEZOMBIES
