
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDcommon.h,v $
$Revision: 1.1 $
$Date: 2004-09-23 20:32:58 $
$Name: not supported by cvs2svn $
$Author: mcschatz $
$Log: not supported by cvs2svn $
Revision 1.4  2004/09/10 12:31:43  mschatz
Add standard copyright notice

Revision 1.3  2004/09/09 22:37:40  mschatz
USE_SOAP_UID support

Revision 1.2  2004/07/19 14:54:46  mpop
Added a fix for CURL on linux

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

Revision 1.8  2000/07/25 14:32:32  cmobarry
Linting for PC/Linux

Revision 1.7  2000/03/02 17:49:51  sdmurphy
extended err msg buffer

Revision 1.6  1999/10/15 15:00:50  sdmurphy
added timeout info

Revision 1.5  1999/07/14 17:24:16  stine
update_cds script was executed against these files.
Only one manual modification - to SYS_UIDcommon.h - I
put in a #include <cds.h> so that it would find the
newfangled cds_* typedefs. Previously, it must have been
using those defined elsewhere in the system.

Revision 1.4  1999/03/05 21:01:01  jlscott
Changed C++ style comments to C style comments (mostly so that client programs
can compile with the -std1 option).

Revision 1.3  1999/03/04 22:15:48  sdmurphy
new size info code

Revision 1.2  1999/01/13 14:31:14  sdmurphy
version 0 prelim

Revision 1.1  1998/12/30 20:06:22  sdmurphy
Renamed uid_common.h to SYS_UIDcommon.h

Revision 1.4  1998/12/24 15:55:06  sdmurphy
externed certain data

Revision 1.3  1998/12/22 12:08:05  sdmurphy
added header (forgot previously)


**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#ifndef UID_COMMON_H
#define UID_COMMON_H


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

/* FOR JTC */
#define JTC_MAX_REQUEST_SIZE        1000000
#define JTC_GUID_REQUEST_URL_MAX_SIZE 2048
#define JTC_GUID_HTTP_RESPONSE_MAX_SIZE 4096
#define JTC_GUID_NUM_BUFFER_SIZE 100

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
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netdb.h>
#include <sys/un.h>
#include <errno.h>
#include <unistd.h>
#include <rpc/xdr.h>
#include <time.h>
#include <sys/time.h>
#ifdef _OSF_SOURCE
#include <sys/timers.h>
#endif
#include <signal.h>
#include <limits.h>

#ifndef USE_SOAP_UID
#include <curl/curl.h>
#include <curl/types.h>
#include <curl/easy.h>
#endif


/* structs *******************************************************/
struct JTC_GUIDMemoryStruct {
  char *memory;
  size_t size;
};

/* variables ********************************************************/
extern char   SYS_UIDerr_str[UID_ERR_STR_SIZE];
extern char   SYS_UIDmessage_array[UID_MESSAGE_SIZE];
extern char   SYS_UIDtype;
extern char   SYS_UIDdebug_flag;

/* functions */
void   SYS_UIDperr(const char* message);
cds_int32  SYS_UIDreadn(cds_int32 fd, char* ptr, cds_int32 nbytes);
cds_int32  SYS_UIDwriten(cds_int32 fd, char* ptr, cds_int32 nbytes);
cds_int32  SYS_UIDpackUIDMessageXdr(cds_uint64* uid, cds_int32 status);
cds_int32  SYS_UIDunpackUIDMessageXdr(cds_uint64* uid, cds_int32* status);
cds_int32  SYS_UIDpackUIDRequestXdr(char* writebuffer, cds_int32 status, cds_uint64 request_size);
cds_int32  SYS_UIDunpackUIDRequestXdr(char* readbuffer, cds_int32* status, cds_uint64* request_size);
void   SYS_UIDlogMessage(const char* message);

#ifndef i386
/* this is ifdef'ed out of the unistd.h for some reason and needs to be here */
//extern int usleep(useconds_t useconds);
#endif

#ifdef __cplusplus
}
#endif

#endif






