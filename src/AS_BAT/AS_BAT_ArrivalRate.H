
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

#ifndef INCLUDE_AS_BAT_ARRIVALRATE
#define INCLUDE_AS_BAT_ARRIVALRATE

static const char *rcsid_INCLUDE_AS_BAT_ARRIVALRATE = "$Id: AS_BAT_ArrivalRate.H,v 1.2 2011-12-29 09:26:03 brianwalenz Exp $";

double getGlobalArrivalRate(UnitigVector &unitigs, long total_random_frags_in_genome=0, long genome_size=0);

#endif  //  INCLUDE_AS_BAT_ARRIVALRATE
