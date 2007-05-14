
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
   SYS_UIDclient(uint64 block_size);
   ~SYS_UIDclient();

   void              Initialize(void);
   int32         GetMaxUIDSize(uint64* size);
   int32         GetLastUIDInterval(uint64* last_UID);
   int32         GetNewUIDInterval(uint64* new_UID);
   void              SetUIDSize(uint64 block_size);
   int32         GetNextUID(uint64* uid);
   int32         GetLastUID(uint64* uid);

 private:

  // functions /////////////////////////////////////////////
  void               GetNewUIDFromServer(void);
  void               HandleConnectionError(void);
  void               SetUIDInterval(uint64 a, uint64 a_size, 
                                    uint64 b, uint64 b_size);
  int32          QueryServer(int32 code, uint64* interval);

  // regular types //////////////////////////////////////////
  uint64         interval_UID[4];
  int32          status;
  uint64         size_UID;
  uint64         increment_offset_UID;
  uint64         increment_max_UID;

  // communication-related - custom per OS ///////////////////
  int32          CreateConnection(void);
  int32          ConfigureConnection(void);
  void               ReceiveServerMessage(int32 code, uint64* interval);
  void               CloseConnection(void);
  int32          GetServerHostInfo(void);

  int32           server_connection_id;
  struct sockaddr_in  server_connection_data;
  int32           server_connection_data_size;
  struct hostent*     server_host_info;
  char                server_host_name[300];
  int32           server_tcp_port;

  int32          CreateFailsafeConnection(void);
  int32          ConfigureFailsafeConnection(void);
  void               ReceiveFailsafeServerMessage(int32 code, uint64* interval);
  void               CloseFailsafeConnection(void);
  int32          GetFailsafeServerHostInfo(void);

  int32           failsafe_server_connection_id;
  struct sockaddr_in  failsafe_server_connection_data;
  int32           failsafe_server_connection_data_size;
  struct hostent*     failsafe_server_host_info;
  char                failsafe_server_host_name[300];
  int32           failsafe_server_tcp_port;

};

#endif






