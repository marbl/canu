
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
/* 	$Id: PrimitiveVA.h,v 1.2 2004-09-23 20:33:59 mcschatz Exp $	 */
/**************************************************************************
 *  PrimitiveVA.h
 *  
 *  Saul A. Kravitz 9/99
 *
 *  Declarations of Variable Arrays (VA) for primitive types
 *  
 *************************************************************************/
#ifndef PRIMITIVE_VA
#define PRIMITIVE_VA

#include "AS_UTL_Var.h"

typedef void *PtrT;

VA_DEF(char)
VA_DEF(short)
VA_DEF(int)
VA_DEF(uint)
VA_DEF(int16)
VA_DEF(cds_int16)
VA_DEF(uint16)
VA_DEF(cds_uint16)
VA_DEF(int32)
VA_DEF(cds_int32)
VA_DEF(uint32)
VA_DEF(cds_uint32)
VA_DEF(uint64)
VA_DEF(cds_uint64)
VA_DEF(PtrT)

VA_DEF(CDS_UID_t)
VA_DEF(CDS_IID_t)
VA_DEF(CDS_CID_t)
VA_DEF(CDS_COORD_t)
#endif
