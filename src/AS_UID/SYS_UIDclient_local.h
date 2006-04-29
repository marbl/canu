
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

#ifndef UID_CLIENT_LOCAL_H
#define UID_CLIENT_LOCAL_H

// local functions /////////

static cds_int32                QueryServer(cds_int32 code, cds_uint64* interval);
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
static void                     ReceiveFailsafeServerMessage(cds_int32 code, cds_uint64* interval);
static void                     CloseFailsafeConnection(void);
static cds_int32                GetFailsafeServerHostInfo(void);

static void                     SetUIDInterval(cds_uint64 a, cds_uint64 a_size, 
                                               cds_uint64 b, cds_uint64 b_size);
static void                     Initialize(void);

// local vars /////////

static cds_uint64               interval_UID[4];
static cds_int32                status                            = UID_CODE_START;
static cds_int32                server_request_mode               = UID_CODE_OK;

static cds_int32                server_connection_id              = 0;
static struct sockaddr_in       server_connection_data;
static cds_int32                server_connection_data_size       = 0;
static struct hostent*          server_host_info                  = NULL;
static cds_int32                server_tcp_port                   = 0;
char                            server_host_name[300];

static cds_int32                failsafe_server_connection_id     = 0;
static struct sockaddr_in       failsafe_server_connection_data;
static cds_int32                failsafe_server_connection_data_size = 0;
static struct hostent*          failsafe_server_host_info         = NULL;
static cds_int32                failsafe_server_tcp_port          = 0;
char                            failsafe_server_host_name[300];

static cds_uint64               size_UID                          = 1L;
static cds_uint64               increment_offset_UID              = 0L;
static cds_uint64               increment_max_UID                 = 0L;

#endif
