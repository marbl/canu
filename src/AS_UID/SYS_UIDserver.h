
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDserver.h,v $
$Revision: 1.3 $
$Date: 2005-03-22 19:49:28 $
$Name: not supported by cvs2svn $
$Author: jason_miller $
$Log: not supported by cvs2svn $
Revision 1.2  2004/09/10 12:31:43  mschatz
Add standard copyright notice

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

Revision 1.4  1999/10/15 15:00:51  sdmurphy
added timeout info

Revision 1.3  1999/07/14 17:24:16  stine
update_cds script was executed against these files.
Only one manual modification - to SYS_UIDcommon.h - I
put in a #include <cds.h> so that it would find the
newfangled cds_* typedefs. Previously, it must have been
using those defined elsewhere in the system.

Revision 1.2  1999/01/13 14:31:46  sdmurphy
version 0 prelim

Revision 1.1  1998/12/30 20:07:15  sdmurphy
Renamed uid_server.h to SYS_UIDserver.h

Revision 1.2  1998/12/21 16:43:01  sdmurphy
added daemon logging support

Revision 1.1  1998/12/17 21:29:56  sdmurphy
include file for uid server

**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#ifndef UID_SERVER_H
#define UID_SERVER_H

cds_int32  SYS_UIDserverInitialize(cds_int32 argc, char** argv);
cds_int32  SYS_UIDserverStart(void);
void       SYS_UIDparseOptions(int argc, char** argv);

#endif




