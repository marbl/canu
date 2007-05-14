
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
static uint64  index_size                 = 1;                
static uint64  index_update_increment     = 1;    
static uint64  max_block_size             = 10000;            
static int32   status_code                = UID_CODE_START;           
static int32   client_connection_id;      
static int32   connection_id;             
static uint64  current_position_UID       = 0L;     
static uint64  index_start_UID            = 0L;           
static uint64  current_UID                = 0L;               
static uint64  interval_UID[4];           
static char    unrecoverable_error_flag   = UID_OK;
static int32   tcp_port;
static int64   ma_alert_time;
static int64   me_alert_time;
static int32   ma_alert_interval          = 0;
static int32   me_alert_interval          = 0;
static char*   ma_address_list            = NULL;
static char*   me_address_list            = NULL;
static char*   mf_address_list            = NULL;
static char    kill_option_invoked_flag   = 0;

static void    Usage(int32 argc, char** argv);
static void    ParseOptions(int32 argc, char** argv);
static int32   UpdateFromPositionFile(void);
static int32   PositionWrite(uint64 position);
static int32   PositionRead(void);
static int32   IncrementUpdatePositionFile(void);
static int32   UIDIsValid(uint64 block_size);
static int32   CleanConnectionPath(void);
static int32   CreateConnection(void);
static int32   RegisterConnection(void);
static int32   ActivateConnection(void);
static int32   AcceptClientConnection(void);
static int32   SendClientMessage(void);
static void    ReadClientRequest(int32* client_status, uint64* request_size);
static void    WriteClientRequest(void);
static void    DebugClientMessage(void);
static void    GetNextUID(uint64 size);
static void    CloseConnection(void);
static void    SetUIDInterval(uint64 a, uint64 a_size, uint64 b, uint64 b_size);
static void    SetStatus(int32 status, char flag);
static void    InitializeAlerts(void);
static void    TriggerMaAlert(void);
static void    TriggerMeAlert(void);
static void    TriggerMfAlert(void);
static void    SendMailList(const char* subject, const char* emessage, const char* address_list);
static int32   IssueKillSignal(void);
static void    ShutdownServer(void);
static int AttemptPortInitialization(void);

static void    FillAddressList(const char*    option, 
                               char**         list, 
                               int32          argc, 
                               char**         argv,
                               int32*         interval);

#endif
