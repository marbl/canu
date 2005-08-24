
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

/**********************************************************************
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/SYS_UIDclient.h,v $
$Revision: 1.4 $
$Date: 2005-08-24 10:57:43 $
$Name: not supported by cvs2svn $
$Author: brianwalenz $
$Log: not supported by cvs2svn $
Revision 1.3  2005/03/22 19:49:28  jason_miller
The TIGR tip as of March 22 2005. Commit by Jason Miller at TIGR.

Revision 1.3  2004/09/10 12:31:43  mschatz
Add standard copyright notice

Revision 1.2  2004/09/09 22:38:51  mschatz
USE_SOAP_UID support

Revision 1.1  2004/06/24 12:51:06  mpop
Added AS_UID

Revision 1.2  2003/05/09 21:04:01  mpop
Dos2unixed all files.
Modified c_make.as to set SEP_PATH relative to LOCAL_WORK

Revision 1.1.1.1  2003/05/08 18:40:11  aaronhalpern
versions from TIGR

Revision 1.2  2001/09/25 23:03:20  mpop
Dos2Unixed

Revision 1.1.1.1  2001/09/25 20:21:05  mpop
Celera Assembler

Revision 1.4  1999/07/14 17:24:16  stine
update_cds script was executed against these files.
Only one manual modification - to SYS_UIDcommon.h - I
put in a #include <cds.h> so that it would find the
newfangled cds_* typedefs. Previously, it must have been
using those defined elsewhere in the system.

Revision 1.3  1999/03/04 22:08:04  sdmurphy
added maxSize func

Revision 1.2  1999/01/13 14:30:37  sdmurphy
version 0 prelim

Revision 1.1  1998/12/30 20:03:56  sdmurphy
Renamed uid_client.h SYS_UIDclient.h

Revision 1.3  1998/12/21 19:02:32  sdmurphy
support for uid incrementer

Revision 1.2  1998/12/18 18:05:20  sdmurphy
changed return type of set_UID_size to void

Revision 1.1  1998/12/17 21:29:27  sdmurphy
include file for uid client API

**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#ifndef UID_CLIENT_H
#define UID_CLIENT_H

cds_int32    SYS_UIDgetLastUIDInterval(cds_uint64* interval);
cds_int32    SYS_UIDgetNewUIDInterval(cds_uint64* interval);
cds_int32    SYS_UIDgetMaxUIDSize(cds_uint64* size);
void         SYS_UIDsetUIDSize(cds_uint64 block_size);
cds_int32    SYS_UIDgetNextUID(cds_uint64* uid);
cds_int32    SYS_UIDgetLastUID(cds_uint64* uid);
void         SYS_UIDset_euid_server(const char * servers);

//  The simple UID client interface, from AS_TER
//


// Allocates blockSize many UIDs from the UID server if real is
// TRUE. Otherwise it allocates some dummy numbers.
//
int32 get_uids(uint64 blockSize, uint64 *interval, int32 real);

// Returns the next available uid. (A real if real==TRUE, o.w. a fake
// UID
//
int32 get_next_uid(uint64 *uid, int32 real);

// Checks whether the SYS_UID* variables are set.  If not it exits
// proposing a standard value
//
void check_environment(void);

//  Sets the initial UID to return, if get_uids() is returning dummy
//  numbers.
//
void    set_start_uid(uint64 s);
uint64  get_start_uid(void);


#endif






