
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
/*********************************************************************
 * Celera Whole Genome Shotgun Assembler.
 * AS_MSG - Directory for "Proto-I/O", the file-based data pipeline.
 * Utility - Utility functions for the I/O. 
 * Jason Miller, 7/2001.
 *********************************************************************/

/* $Id: AS_MSG_Utility.h,v 1.1.1.1 2004-04-14 13:52:14 catmandew Exp $ */

#ifndef AS_MSG_UTILITY_INCLUDE
#define AS_MSG_UTILITY_INCLUDE

#include "AS_MSG_pmesg.h"

FragType AS_MSG_SafeConvert_charToFragType (const char inch, bool strict) ;

#endif /* AS_MSG_UTILITY_INCLUDE */
