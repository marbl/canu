
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDserver_local.h,v $
$Revision: 1.1 $
$Date: 2004-09-23 20:32:58 $
$Name: not supported by cvs2svn $
$Author: mcschatz $
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

Revision 1.5  1999/10/15 15:00:51  sdmurphy
added timeout info

Revision 1.4  1999/07/14 17:24:17  stine
update_cds script was executed against these files.
Only one manual modification - to SYS_UIDcommon.h - I
put in a #include <cds.h> so that it would find the
newfangled cds_* typedefs. Previously, it must have been
using those defined elsewhere in the system.

Revision 1.3  1999/03/04 22:17:15  sdmurphy
new req and mail info

Revision 1.2  1999/01/13 14:32:02  sdmurphy
version 0 prelim

Revision 1.1  1999/01/07 09:47:40  sdmurphy
local include for SYS_UIDserver


**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#ifndef UID_SERVER_LOCAL_H
#define UID_SERVER_LOCAL_H


static char*   positionfile_name          = NULL;         
static cds_uint64  index_size                 = 1;                
static cds_uint64  index_update_increment     = 1;    
static cds_uint64  max_block_size             = 10000;            
static cds_int32   status_code                = UID_CODE_START;           
static cds_int32   client_connection_id;      
static cds_int32   connection_id;             
static cds_uint64  current_position_UID       = 0L;     
static cds_uint64  index_start_UID            = 0L;           
static cds_uint64  current_UID                = 0L;               
static cds_uint64  interval_UID[4];           
static char    unrecoverable_error_flag   = UID_OK;
static cds_int32   tcp_port;
static cds_int64   ma_alert_time;
static cds_int64   me_alert_time;
static cds_int32   ma_alert_interval          = 0;
static cds_int32   me_alert_interval          = 0;
static char*   ma_address_list            = NULL;
static char*   me_address_list            = NULL;
static char*   mf_address_list            = NULL;
static char    kill_option_invoked_flag   = 0;

static void    Usage(cds_int32 argc, char** argv);
static void    ParseOptions(cds_int32 argc, char** argv);
static cds_int32   UpdateFromPositionFile(void);
static cds_int32   PositionWrite(cds_uint64 position);
static cds_int32   PositionRead(void);
static cds_int32   IncrementUpdatePositionFile(void);
static cds_int32   UIDIsValid(cds_uint64 block_size);
static cds_int32   CleanConnectionPath(void);
static cds_int32   CreateConnection(void);
static cds_int32   RegisterConnection(void);
static cds_int32   ActivateConnection(void);
static cds_int32   AcceptClientConnection(void);
static cds_int32   SendClientMessage(void);
static void    ReadClientRequest(cds_int32* client_status, cds_uint64* request_size);
static void    WriteClientRequest(void);
static void    DebugClientMessage(void);
static void    GetNextUID(cds_uint64 size);
static void    CloseConnection(void);
static void    SetUIDInterval(cds_uint64 a, cds_uint64 a_size, cds_uint64 b, cds_uint64 b_size);
static void    SetStatus(cds_int32 status, char flag);
static void    InitializeAlerts(void);
static void    TriggerMaAlert(void);
static void    TriggerMeAlert(void);
static void    TriggerMfAlert(void);
static void    SendMailList(const char* subject, const char* emessage, const char* address_list);
static cds_int32   IssueKillSignal(void);
static void    ShutdownServer(void);
static int AttemptPortInitialization(void);

static void    FillAddressList(const char*    option, 
                               char**         list, 
                               cds_int32          argc, 
                               char**         argv,
                               cds_int32*         interval);

#endif
