
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

#ifndef UID_CLIENT_H
#define UID_CLIENT_H

//  The simple UID client interface, from AS_TER

int32    SYS_UIDgetLastUIDInterval(uint64* interval);
int32    SYS_UIDgetNewUIDInterval(uint64* interval);
int32    SYS_UIDgetMaxUIDSize(uint64* size);
void         SYS_UIDsetUIDSize(uint64 block_size);
int32    SYS_UIDgetNextUID(uint64* uid);
int32    SYS_UIDgetLastUID(uint64* uid);
void         SYS_UIDset_euid_server(const char * servers);
void         SYS_UIDset_euid_namespace(const char * namespaceName);


// Allocates blockSize many UIDs from the UID server if real is
// TRUE. Otherwise it allocates some dummy numbers.
//
void get_uids(uint64 blockSize, uint64 *interval, int32 real);


// Returns the next available uid. (A real if real==TRUE, o.w. a fake
// UID
//
int32 get_next_uid(uint64 *uid, int32 real);


//  Sets the initial UID to return, if get_uids() is returning dummy
//  numbers.
//
void    set_start_uid(uint64 s);
uint64  get_start_uid(void);

#endif






