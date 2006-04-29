
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

#ifndef UID_COMMON_H
#define UID_COMMON_H

/*  New UID server/client disabled/modified the original
 *  code.  This enables those changes.
 *
 *  This also disables the need for these environment variables:
 *    SYS_UID_SERVER_HOSTNAME
 *    SYS_UID_SERVER_PORT
 *    SYS_UID_FAILSAFE_SERVER_HOSTNAME
 *    SYS_UID_FAILSAFE_SERVER_PORT
 */
#define NOT_IMPLEMENTED_JTC

#define UID_OK                          1
#define UID_FAILS                       0

#define UID_DEFAULT_SERVER_TCP_PORT              5001
#define UID_DEFAULT_SERVER_HOST_NAME             "localhost"

#define UID_DEFAULT_FAILSAFE_SERVER_TCP_PORT     5002
#define UID_DEFAULT_FAILSAFE_SERVER_HOST_NAME    "localhost"

#define UID_ERR_STR_SIZE             20000
#define UID_MESSAGE_SIZE                36

#define UID_CLIENT_SEND_TIMEOUT   20
#define UID_CLIENT_RECV_TIMEOUT   20
#define UID_SERVER_SEND_TIMEOUT   20
#define UID_SERVER_RECV_TIMEOUT   20

#define UID_NO_TYPE                     0
#define UID_SERVER_TYPE                 1
#define UID_CLIENT_TYPE                 2

/* both client and server - note - 0 must be */
/* an error code to detect NULL transmits    */
#define UID_CODE_OK                     101
#define UID_CODE_UNKNOWN_ERROR          102
#define UID_CODE_RAN_OUT_OF_SPACE       103
#define UID_CODE_POS_BOUNDS_ERROR       104
#define UID_CODE_POS_CONFIG_ERROR       105
#define UID_CODE_BLOCK_TOO_LARGE        106
#define UID_CODE_REQUEST_ERROR          107
#define UID_CODE_NEED_SIZE_INFO         108

/* server side - never seen by client */
#define UID_CODE_REGISTER_CONN_FAILED   201
#define UID_CODE_ACTIVATE_CONN_FAILED   202
#define UID_CODE_ACCEPT_CONN_FAILED     203
#define UID_CODE_SEND_ERROR             204
#define UID_CODE_SERVER_KILL            205

/* client side only */
#define UID_CODE_NULL_TRANSMISSION        0
#define UID_CODE_START                  301
#define UID_CODE_CREATE_CONN_FAILED     302
#define UID_CODE_CONFIGURE_CONN_FAILED  303
#define UID_CODE_CANT_CONNECT           304
#define UID_CODE_CANT_READ              305
#define UID_CODE_NULL_INTERVAL_PTR      306
#define UID_CODE_INCREMENT_OVERFLOW     307

/* New UID client does not implement SYS_UIDgetMaxUIDSize()
 * as a server query, instead, just returns this constant.
 */
#define UID_MAX_REQUEST_SIZE            1000000

#ifdef __cplusplus
extern "C" {
#endif

// MCS: Put cds.h before arpa/inet.h for C++ builds on alpha
#include <cds.h>

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <sys/un.h>
#include <errno.h>
#include <unistd.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <time.h>
#include <sys/time.h>
#ifdef _OSF_SOURCE
#include <sys/timers.h>
#endif
#include <signal.h>
#include <limits.h>


/* structs *******************************************************/

/* variables ********************************************************/
extern char   SYS_UIDerr_str[UID_ERR_STR_SIZE];
extern char   SYS_UIDmessage_array[UID_MESSAGE_SIZE];
extern char   SYS_UIDtype;
extern char   SYS_UIDdebug_flag;

/* functions */
void       SYS_UIDperr(const char* message);
cds_int32  SYS_UIDreadn(cds_int32 fd, char* ptr, cds_int32 nbytes);
cds_int32  SYS_UIDwriten(cds_int32 fd, char* ptr, cds_int32 nbytes);
cds_int32  SYS_UIDpackUIDMessageXdr(cds_uint64* uid, cds_int32 status);
cds_int32  SYS_UIDunpackUIDMessageXdr(cds_uint64* uid, cds_int32* status);
cds_int32  SYS_UIDpackUIDRequestXdr(char* writebuffer, cds_int32 status, cds_uint64 request_size);
cds_int32  SYS_UIDunpackUIDRequestXdr(char* readbuffer, cds_int32* status, cds_uint64* request_size);
void       SYS_UIDlogMessage(const char* message);

#ifdef __cplusplus
}
#endif

#endif






