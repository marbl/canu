
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
