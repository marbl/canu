
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
/* 	$Id: Stats_CGW.h,v 1.6 2008-06-27 06:29:14 brianwalenz Exp $	 */
#ifndef STATS_CGW_H
#define STATS_CGW_H

/* Statistics -- see Stats_CGW.c */
void GenerateCIGraphStats(void);
void GenerateCIGraph_U_Stats(void);
void GeneratePlacedContigGraphStats(char *phase, int iteration);
void GenerateScaffoldGraphStats(char *phase, int interation);
void GenerateLinkStats(GraphCGW_T *graph, char *phase, int iteration);
void GenerateSurrogateStats(char *label);
void GenerateContigAlignmentStats(char *phase);

void ComputeFragmentMembershipStats(void);

#endif
