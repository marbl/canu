
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
/* 	$Id: RepeatRez.h,v 1.3 2005-03-22 19:07:58 jason_miller Exp $	 */

/************************************************************************
 *  RepeatRez.h
 *  
 *  Saul A. Kravitz 10/99
 *
 *  This is the interface from CGW to REZ
 *
 ************************************************************************/
#ifndef REPEAT_REZ_H
#define REPEAT_REZ_H

#include "Globals_CGW.h"

int  Fill_Gaps
    (Global_CGW * data, char *, int, int redo_index);

int  Hurl_Contained_Rocks
    (char * prefix, int level, int redo_index);

#define AGGRESSIVE_WALKING_STD_DEVS 5.0
#define CONSERVATIVE_WALKING_STD_DEVS 3.0

int  Walk_Gaps
    (Global_CGW * data, char *, int, int startWalkFrom, double gapSizeStdDevs);

int  Throw_Stones
    (char *, int, int);

int  Toss_Contained_Stones
    (char * prefix, int level, int redo_index);

int Inter_Scaffold_Walking(void);

#endif
