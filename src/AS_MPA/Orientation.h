
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
/* $Id: Orientation.h,v 1.1.1.1 2004-04-14 13:52:04 catmandew Exp $ */
#ifndef ORIENTATION_H
#define ORIENTATION_H

typedef enum
{
  SINGLE_A_B,
  SINGLE_B_A,
  SINGLE_NUM_ORIENTATIONS
} SingleOrientation_e;

typedef enum
{
  PAIR_INNIE,
  PAIR_OUTTIE,
  PAIR_NORMAL,
  PAIR_ANTINORMAL,
  PAIR_NUM_ORIENTATIONS
} PairOrientation_e;

typedef enum
{
  MPI_STRETCHED = 0,
  MPI_COMPRESSED,
  MPI_OUTTIE,
  MPI_NORMAL,
  MPI_ANTINORMAL,
  MPI_INVERSION,
  MPI_TRANSPOSITION,
  MPI_SATISFIED,
  MPI_INTERCHROMOSOME,
  MPI_UNKNOWN,
  MPI_NUM_INDICES
} MatePairIndex_e;

#endif
