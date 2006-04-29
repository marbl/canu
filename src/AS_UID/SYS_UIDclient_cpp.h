
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

#ifndef UID_CLIENT_CPP_H
#define UID_CLIENT_CPP_H

class SYS_UIDclient {
  
 public:

   SYS_UIDclient();
   SYS_UIDclient(cds_uint64 block_size);
   ~SYS_UIDclient();

   void              Initialize(void);
   cds_int32         GetMaxUIDSize(cds_uint64* size);
   cds_int32         GetLastUIDInterval(cds_uint64* last_UID);
   cds_int32         GetNewUIDInterval(cds_uint64* new_UID);
   void              SetUIDSize(cds_uint64 block_size);
   cds_int32         GetNextUID(cds_uint64* uid);
   cds_int32         GetLastUID(cds_uint64* uid);

 private:

  // functions /////////////////////////////////////////////
  void               GetNewUIDFromServer(void);
  void               HandleConnectionError(void);
  void               SetUIDInterval(cds_uint64 a, cds_uint64 a_size, 
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
  void               ReceiveServerMessage(cds_int32 code, cds_uint64* interval);
  void               CloseConnection(void);
  cds_int32          GetServerHostInfo(void);

  cds_int32           server_connection_id;
  struct sockaddr_in  server_connection_data;
  cds_int32           server_connection_data_size;
  struct hostent*     server_host_info;
  char                server_host_name[300];
  cds_int32           server_tcp_port;

  cds_int32          CreateFailsafeConnection(void);
  cds_int32          ConfigureFailsafeConnection(void);
  void               ReceiveFailsafeServerMessage(cds_int32 code, cds_uint64* interval);
  void               CloseFailsafeConnection(void);
  cds_int32          GetFailsafeServerHostInfo(void);

  cds_int32           failsafe_server_connection_id;
  struct sockaddr_in  failsafe_server_connection_data;
  cds_int32           failsafe_server_connection_data_size;
  struct hostent*     failsafe_server_host_info;
  char                failsafe_server_host_name[300];
  cds_int32           failsafe_server_tcp_port;

};

#endif






