
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDclient_cpp.h,v $
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

Revision 1.4  1999/07/14 17:24:16  stine
update_cds script was executed against these files.
Only one manual modification - to SYS_UIDcommon.h - I
put in a #include <cds.h> so that it would find the
newfangled cds_* typedefs. Previously, it must have been
using those defined elsewhere in the system.

Revision 1.3  1999/03/04 22:09:56  sdmurphy
added code to client-server communication

Revision 1.2  1999/01/13 14:30:50  sdmurphy
version 0 prelim

Revision 1.1  1998/12/30 20:05:05  sdmurphy
Renamed uid_client_cpp.h to SYS_UIDclient_cpp.h

Revision 1.1  1998/12/23 19:15:14  sdmurphy
include for cpp uid client


**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#ifndef UID_CLIENT_CPP_H
#define UID_CLIENT_CPP_H

class SYS_UIDclient {
  
 public:

   SYS_UIDclient();
   SYS_UIDclient(cds_uint64 block_size);
   ~SYS_UIDclient();

   void          Initialize(void);
   cds_int32         GetMaxUIDSize(cds_uint64* size);
   cds_int32         GetLastUIDInterval(cds_uint64* last_UID);
   cds_int32         GetNewUIDInterval(cds_uint64* new_UID);
   void          SetUIDSize(cds_uint64 block_size);
   cds_int32         GetNextUID(cds_uint64* uid);
   cds_int32         GetLastUID(cds_uint64* uid);

 private:

  // functions /////////////////////////////////////////////
  void           GetNewUIDFromServer(void);
  void           HandleConnectionError(void);
  void           SetUIDInterval(cds_uint64 a, cds_uint64 a_size, 
				cds_uint64 b, cds_uint64 b_size);
  cds_int32          QueryServer(cds_int32 code, cds_uint64* interval);

  // regular types //////////////////////////////////////////
  cds_uint64         interval_UID[4];
  cds_int32          status;
  cds_uint64         size_UID;
  cds_uint64         increment_offset_UID;
  cds_uint64         increment_max_UID;

  // communication-related - custom per OS ///////////////////
  cds_int32          CreateConnection(void);
  cds_int32          ConfigureConnection(void);
  void           ReceiveServerMessage(cds_int32 code, cds_uint64* interval);
  void           CloseConnection(void);
  cds_int32          GetServerHostInfo(void);

  cds_int32               server_connection_id;
  struct sockaddr_in  server_connection_data;
  cds_int32               server_connection_data_size;
  struct hostent*     server_host_info;
  char                server_host_name[300];
  cds_int32               server_tcp_port;

  cds_int32          CreateFailsafeConnection(void);
  cds_int32          ConfigureFailsafeConnection(void);
  void           ReceiveFailsafeServerMessage(cds_int32 code, cds_uint64* interval);
  void           CloseFailsafeConnection(void);
  cds_int32          GetFailsafeServerHostInfo(void);

  cds_int32               failsafe_server_connection_id;
  struct sockaddr_in  failsafe_server_connection_data;
  cds_int32               failsafe_server_connection_data_size;
  struct hostent*     failsafe_server_host_info;
  char                failsafe_server_host_name[300];
  cds_int32               failsafe_server_tcp_port;

};

#endif






