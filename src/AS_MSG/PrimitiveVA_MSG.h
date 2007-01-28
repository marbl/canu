
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
/* 	$Id: PrimitiveVA_MSG.h,v 1.6 2007-01-28 21:52:24 brianwalenz Exp $	 */
/***************************************************************************
 *  PrimitiveVA.h
 *  
 *  Saul A. Kravitz 9/99
 *
 *  Declarations of Variable Arrays (VA) for primitive types declared in MSG
 *  
 **************************************************************************/
#ifndef PRIMITIVE_VA_MSG
#define PRIMITIVE_VA_MSG

#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"

VA_DEF(IntMultiVar)
VA_DEF(IntMultiPos)
VA_DEF(IntElementPos)
VA_DEF(IntUnitigPos)
#endif
