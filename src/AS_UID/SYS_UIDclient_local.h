
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDclient_local.h,v $
$Revision: 1.3 $
$Date: 2005-03-22 19:49:28 $
$Name: not supported by cvs2svn $
$Author: jason_miller $
$Log: not supported by cvs2svn $
Revision 1.3  2004/09/10 12:31:43  mschatz
Add standard copyright notice

Revision 1.2  2004/06/25 04:05:04  ahalpern
Fixes to allow access to JTC uid server from linux

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

Revision 1.3  1999/03/04 22:14:52  sdmurphy
added code to client-server comm

Revision 1.2  1999/01/13 14:31:03  sdmurphy
version 0 prelim

Revision 1.1  1999/01/07 09:47:15  sdmurphy
local include for SYS_UIDclient


**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#ifndef UID_CLIENT_LOCAL_H
#define UID_CLIENT_LOCAL_H

const char JTC_GUID_URL[] = "http://guid.jtc.jcvsf.org:8080/guid/GuidClientServer?Request=GET&Size=";

// local functions /////////

static cds_int32                QueryServer(cds_int32 code, cds_uint64* interval);
size_t                          JTC_GUIDWriteMemoryCallback(void *ptr, size_t size, size_t nmemb, void *data);
CDS_UID_t                       getGUIDBlock(int guidRequestSize);
int                             findGuidStartFromHttpString(char* httpString);
int                             findGuidEndFromHttpString(char* httpString);
static cds_int32                CreateConnection(void);
static cds_int32                ConfigureConnection(void);
static void                     ReceiveServerMessage(cds_int32 code, cds_uint64* interval);
static void                     CloseConnection(void);
static cds_int32                GetServerHostInfo(void);

static cds_int32                CreateFailsafeConnection(void);
static cds_int32                ConfigureFailsafeConnection(void);
static void                 ReceiveFailsafeServerMessage(cds_int32 code, cds_uint64* interval);
static void                 CloseFailsafeConnection(void);
static cds_int32                GetFailsafeServerHostInfo(void);

static void                 SetUIDInterval(cds_uint64 a, cds_uint64 a_size, 
                                                cds_uint64 b, cds_uint64 b_size);
static void                 Initialize(void);

// local vars /////////

static cds_uint64               interval_UID[4];
static cds_int32                status                            = UID_CODE_START;
static cds_int32                server_request_mode               = UID_CODE_OK;

static cds_int32                server_connection_id              = 0;
static struct sockaddr_in   server_connection_data;
static cds_int32                server_connection_data_size       = 0;
static struct hostent*      server_host_info                  = NULL;
static cds_int32                server_tcp_port                   = 0;
char                        server_host_name[300];

static cds_int32                failsafe_server_connection_id     = 0;
static struct sockaddr_in   failsafe_server_connection_data;
static cds_int32                failsafe_server_connection_data_size = 0;
static struct hostent*      failsafe_server_host_info         = NULL;
static cds_int32                failsafe_server_tcp_port          = 0;
char                        failsafe_server_host_name[300];

static cds_uint64               size_UID                          = 1L;
static cds_uint64               increment_offset_UID              = 0L;
static cds_uint64               increment_max_UID                 = 0L;

#endif
