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

#ifndef CHECKPOINTS_CGW_H
#define CHECKPOINTS_CGW_H

//  Defines the logical checkpoints

#define CHECKPOINT_AFTER_READING_INPUT 0
#define CHECKPOINT_AFTER_UNITIG_SPLITTING 1
#define CHECKPOINT_AFTER_BUILDING_EDGES 2
#define CHECKPOINT_AFTER_BUILDING_CIGRAPH 3
#define CHECKPOINT_AFTER_BUILDING_SCAFFOLDS 4
#define CHECKPOINT_AFTER_BUILDING_AND_CLEANING_SCAFFOLDS 5
#define CHECKPOINT_BEFORE_CONSERVATIVE_WALKING 6  //  NOP
#define CHECKPOINT_BEFORE_1ST_SCAFF_MERGE 7
#define CHECKPOINT_BEFORE_STONES 8
#define CHECKPOINT_BEFORE_AGGRESSIVE_WALKING 9  //  NOP
#define CHECKPOINT_BEFORE_2ND_SCAFF_MERGE 10
#define CHECKPOINT_BEFORE_FINAL_ROCKS 11
#define CHECKPOINT_BEFORE_PARTIAL_STONES 12
#define CHECKPOINT_BEFORE_FINAL_CONTAINED_STONES 13
#define CHECKPOINT_BEFORE_INTER_SCAFFOLD_WALKING 14
#define CHECKPOINT_BEFORE_FINAL_CLEANUP 15
#define CHECKPOINT_BEFORE_RESOLVE_SURROGATES 16

#endif
