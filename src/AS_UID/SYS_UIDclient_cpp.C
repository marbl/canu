
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDclient_cpp.C,v $
$Revision: 1.3 $
$Date: 2005-03-22 19:49:28 $
$Name: not supported by cvs2svn $
$Author: jason_miller $
$Log: not supported by cvs2svn $
Revision 1.2  2004/09/10 12:31:43  mschatz
Add standard copyright notice

Revision 1.1  2004/06/24 12:51:06  mpop
Added AS_UID

Revision 1.2  2003/05/09 21:04:03  mpop
Dos2unixed all files.
Modified c_make.as to set SEP_PATH relative to LOCAL_WORK

Revision 1.1.1.1  2003/05/08 18:40:11  aaronhalpern
versions from TIGR

Revision 1.2  2001/09/25 23:03:20  mpop
Dos2Unixed

Revision 1.1.1.1  2001/09/25 20:21:05  mpop
Celera Assembler

Revision 1.6  1999/07/14 17:24:33  stine
update_cds script was executed against these files.

Revision 1.5  1999/01/28 14:58:06  sdmurphy
added GetMaxUIDSize and QueryServer

Revision 1.4  1999/01/13 16:44:28  sdmurphy
took out unnecessary error msgs

Revision 1.3  1999/01/13 16:06:05  sdmurphy
fixed bug in Initialize()

Revision 1.2  1999/01/13 14:28:19  sdmurphy
version 0 prelim

Revision 1.1  1998/12/30 19:48:58  sdmurphy
Renamed uid_client_cpp.C SYS_UIDclient_cpp.C

Revision 1.1  1998/12/22 21:17:41  sdmurphy
starting c++ version



**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#include "SYS_UIDcommon.h"
#include "SYS_UIDerror.h"
#include <iostream.h>
#include "SYS_UIDclient_cpp.h"

/*******************************************************************************

Description: Default constructor


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
SYS_UIDclient::SYS_UIDclient()
{
   SetUIDSize(1L);
   Initialize();
}


/*******************************************************************************

Description: Constructor with size_UID argument


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
SYS_UIDclient::SYS_UIDclient(cds_uint64 request_size)
{
   SetUIDSize(request_size);
   Initialize();
}
/*******************************************************************************

Description: Constructor-independent initializer


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void SYS_UIDclient::Initialize(void)
{
   char* port_string;
   char* host_string;
   char* failsafe_port_string;
   char* failsafe_host_string;

   SYS_UIDtype = UID_CLIENT_TYPE;
   SetUIDInterval(0L, 0L, 0L, 0L);
   status = UID_CODE_START;

   ////// REGULAR SERVER /////////////////////////////////////////////
   if ( (port_string = getenv("SYS_UID_SERVER_PORT")) == NULL)
      server_tcp_port = UID_DEFAULT_SERVER_TCP_PORT;
   else
      server_tcp_port = atoi(port_string);

   if ( (host_string = getenv("SYS_UID_SERVER_HOST_NAME")) == NULL)
      strcpy(server_host_name, UID_DEFAULT_SERVER_HOST_NAME);
   else
   {
      if (strlen(host_string) < 300)
      {
         strcpy(server_host_name, host_string);
      }
      else
      {
         SYS_UIDerrorMsg("UID Client: host name too long, using default server host name\n");
         strcpy(server_host_name, UID_DEFAULT_SERVER_HOST_NAME);
      }   
   }

   ////// FAILSAFE SERVER /////////////////////////////////////////////
   if ( (failsafe_port_string = getenv("SYS_UID_FAILSAFE_SERVER_PORT")) == NULL)
      failsafe_server_tcp_port = UID_DEFAULT_FAILSAFE_SERVER_TCP_PORT;
   else
      failsafe_server_tcp_port = atoi(failsafe_port_string);

   if ( (failsafe_host_string = getenv("SYS_UID_FAILSAFE_SERVER_HOST_NAME")) == NULL)
      strcpy(failsafe_server_host_name, UID_DEFAULT_FAILSAFE_SERVER_HOST_NAME);
   else
   {
      if (strlen(failsafe_host_string) < 300)
      {
         strcpy(failsafe_server_host_name, failsafe_host_string);
      }
      else
      {
         SYS_UIDerrorMsg("UID Client: failsafe host name too long, using default failsafe\
 server host name\n");
         strcpy(failsafe_server_host_name, UID_DEFAULT_FAILSAFE_SERVER_HOST_NAME);
      }   
   }
}


/*******************************************************************************

Description: Destructor


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
SYS_UIDclient::~SYS_UIDclient()
{
}


/*******************************************************************************

Description: public method for setting request size


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void SYS_UIDclient::SetUIDSize(cds_uint64 size)
{
   size_UID = size;
}


/*******************************************************************************

Description: Gets maximum UID size


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  SYS_UIDclient::GetMaxUIDSize(cds_uint64* size)
{
   cds_uint64  interval_buffer[4];

   if (size == NULL)
      return UID_CODE_NULL_INTERVAL_PTR;

   return QueryServer(UID_CODE_NEED_SIZE_INFO, size);
}

/*******************************************************************************

Description: utility for setting interval values


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void SYS_UIDclient::SetUIDInterval(cds_uint64 a, cds_uint64 a_size,
                                 cds_uint64 b, cds_uint64 b_size)
{
   interval_UID[0] = a;
   interval_UID[1] = a_size;
   interval_UID[2] = b;
   interval_UID[3] = b_size;

   increment_offset_UID        = 0L;
   increment_max_UID           = a_size + b_size;
}


/*******************************************************************************

Description: 

   External function for getting the current increment UID. It returns the
   current value of the incrementer, and increments the value for next time.
   Note that even though it increments the value, it returns the
   pre-incremented value.


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32   SYS_UIDclient::GetNextUID(cds_uint64* uid)
{
   // check limit
   if (increment_offset_UID >= increment_max_UID)
   {
     *uid = 0L;
     return UID_CODE_INCREMENT_OVERFLOW;
   }
   
   // check which interval 
   if (increment_offset_UID < interval_UID[1])
     *uid = increment_offset_UID + interval_UID[0];
   else
     *uid = (increment_offset_UID - interval_UID[1]) + interval_UID[2];

   // increment for next call
   increment_offset_UID++;

   return status;
}

/*******************************************************************************

Description: External function for getting the last incremented UID


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32   SYS_UIDclient::GetLastUID(cds_uint64* uid)
{
   // check to see if incrementer already called
   if (increment_offset_UID == 0)
   {
     *uid = 0L;
     return UID_CODE_INCREMENT_OVERFLOW;
   }
   
   // check which interval 
   if (increment_offset_UID < interval_UID[1])
     *uid = increment_offset_UID + interval_UID[0];
   else
     *uid = (increment_offset_UID - interval_UID[1]) + interval_UID[2];

   return status;
}

/*******************************************************************************

Description: public method for getting the current (last) interval vals


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDclient::GetLastUIDInterval(cds_uint64* interval)
{
   if (interval == NULL)
      return UID_CODE_NULL_INTERVAL_PTR;
 
   if (status == UID_CODE_START)
      return UID_CODE_START;

   interval[0] = interval_UID[0];
   interval[1] = interval_UID[1];
   interval[2] = interval_UID[2];
   interval[3] = interval_UID[3];

   return status;
}

/*******************************************************************************

Description: public method for getting a new UID interval


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDclient::GetNewUIDInterval(cds_uint64* interval)
{
   cds_int32 result_status;

   if (interval == NULL)
      return UID_CODE_NULL_INTERVAL_PTR;

   SetUIDInterval(0L, 0L, 0L, 0L);

   result_status = QueryServer(UID_CODE_OK, interval_UID);

   increment_offset_UID = 0L;
   increment_max_UID = interval_UID[1] + interval_UID[3];

   interval[0] = interval_UID[0];
   interval[1] = interval_UID[1];
   interval[2] = interval_UID[2];
   interval[3] = interval_UID[3];

   return result_status;
}

/*******************************************************************************

Description: handles communication setup with server


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDclient::QueryServer(cds_int32 code, cds_uint64* interval)
{
   char failsafe_flag = 0;

   if (CreateConnection() == UID_FAILS)
   {
      failsafe_flag = 1;
   }
   else if (ConfigureConnection() == UID_FAILS)
   {
      CloseConnection();
      failsafe_flag = 1;
   }
   else 
   {
      status = UID_CODE_OK;
      ReceiveServerMessage(code, interval);
      CloseConnection();
      if (status != UID_CODE_OK)
         failsafe_flag = 1;
   }

   if (failsafe_flag == 1)
   {
      status = UID_CODE_OK;

      if (CreateFailsafeConnection() == UID_FAILS)
	{
        status = UID_CODE_CREATE_CONN_FAILED;
        return status;
        }
   
      if (ConfigureFailsafeConnection() == UID_FAILS)
      {
         CloseFailsafeConnection();
         status = UID_CODE_CONFIGURE_CONN_FAILED;
         return status;
      }
      ReceiveFailsafeServerMessage(code, interval);
      CloseFailsafeConnection(); 
   }

   return status;
}

/*******************************************************************************

Description: gets server host info


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDclient::GetServerHostInfo(void)
{
   // uses the member var 'server_host_info' of type struct hostent
   if ( (server_host_info = gethostbyname(server_host_name)) == NULL)
   {
      sprintf(SYS_UIDerr_str, "UID Client::get_server_host: could not find address \
for host %s \n", server_host_name);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      return UID_FAILS;
   }
   return UID_OK;
}
/*******************************************************************************

Description: gets failsafe server host info


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDclient::GetFailsafeServerHostInfo(void)
{
   // uses the member var 'failsafe_server_host_info' of type struct hostent
   if ( (failsafe_server_host_info = gethostbyname(failsafe_server_host_name)) == NULL)
   {
      sprintf(SYS_UIDerr_str, "UID Client: could not find address for \
failsafe host %s\n", failsafe_server_host_name);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      return UID_FAILS;
   }
   return UID_OK;
}

/*******************************************************************************

Description: creates server connection


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  SYS_UIDclient::CreateConnection(void)
{
   if ( (server_connection_id = socket(AF_INET, SOCK_STREAM, 0)) < 0)
   {
      sprintf(SYS_UIDerr_str, "UID Client::CreateConnection: could not open socket\n");
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      return UID_FAILS;
   }
   return UID_OK;
}

/*******************************************************************************

Description: creates failsafe server connection


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  SYS_UIDclient::CreateFailsafeConnection(void)
{
   if ( (failsafe_server_connection_id = socket(AF_INET, SOCK_STREAM, 0)) < 0)
   {
      sprintf(SYS_UIDerr_str, "UID Client: could not open socket for failsafe server\n");
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      return UID_FAILS;
   }
   return UID_OK;
}

/*******************************************************************************

Description: configures server connection


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  SYS_UIDclient::ConfigureConnection(void)
{
   if (GetServerHostInfo() == UID_FAILS)
      return UID_FAILS;
   bzero( (char*) &server_connection_data, sizeof(server_connection_data) );
   server_connection_data.sin_family = AF_INET;
   server_connection_data.sin_addr.s_addr = 
      inet_addr( inet_ntoa( *((struct in_addr*)(server_host_info->h_addr_list[0])) ) );
   server_connection_data.sin_port = htons(server_tcp_port);
   server_connection_data_size = sizeof(server_connection_data);
   return UID_OK;
}

/*******************************************************************************

Description: configures failsafe server connection 


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  SYS_UIDclient::ConfigureFailsafeConnection(void)
{
   if (GetFailsafeServerHostInfo() == UID_FAILS)
      return UID_FAILS;
   bzero( (char*) &failsafe_server_connection_data, sizeof(failsafe_server_connection_data) );
   failsafe_server_connection_data.sin_family = AF_INET;
   failsafe_server_connection_data.sin_addr.s_addr = 
      inet_addr( inet_ntoa( *((struct in_addr*)(failsafe_server_host_info->h_addr_list[0])) ) );
   failsafe_server_connection_data.sin_port = htons(failsafe_server_tcp_port);
   failsafe_server_connection_data_size = sizeof(failsafe_server_connection_data);
   return UID_OK;
}

/*******************************************************************************

Description: receives server message


Input:


Output:


Returns:


Globals:


Notes:

   No error messages because these are deferred to 
   ReceiveFailsafeServerMessage().


*******************************************************************************/
void  SYS_UIDclient::ReceiveServerMessage(cds_int32 code, cds_uint64* interval)
{
   char  writebuffer[12];

   if(connect(server_connection_id, (struct sockaddr*) &server_connection_data, 
      server_connection_data_size) < 0)
   {
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_CONNECT;
      return;
   }

   SYS_UIDpackUIDRequestXdr(writebuffer, code, size_UID);
 
   if (SYS_UIDwriten(server_connection_id, writebuffer, 12) == UID_FAILS)
   {
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_READ;
      return;
   }

   // read the message from the server
   if ( SYS_UIDreadn(server_connection_id, (char*)SYS_UIDmessage_array, 
		     UID_MESSAGE_SIZE) == UID_FAILS)
   {
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_READ;
      return;
   }

   // unpackage
   if (SYS_UIDunpackUIDMessageXdr(interval, &status) == UID_FAILS)
   {
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_READ;
      return;
   }
 
   return;
}

/*******************************************************************************

Description: receives failsafe server message


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void  SYS_UIDclient::ReceiveFailsafeServerMessage(cds_int32 code, cds_uint64* interval)
{
   char  writebuffer[12];

   if(connect(failsafe_server_connection_id, 
	      (struct sockaddr*) &failsafe_server_connection_data, 
      failsafe_server_connection_data_size) < 0)
   {
      SYS_UIDerrorMsg("UID Client: can't connect to failsafe server\n");
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_CONNECT;
      return;
   }

   SYS_UIDpackUIDRequestXdr(writebuffer, code, size_UID);
 
   if (SYS_UIDwriten(failsafe_server_connection_id, writebuffer, 12) == UID_FAILS)
   {
      SYS_UIDerrorMsg("UID Client: writen failed to failsafe server\n");
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_READ;
      return;
   }

   // read the message from the server
   if ( SYS_UIDreadn(failsafe_server_connection_id, (char*)SYS_UIDmessage_array, 
		     UID_MESSAGE_SIZE) == UID_FAILS)
   {
      SYS_UIDerrorMsg("UID Client: readn failed from failsafe server\n");
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_READ;
      return;
   }

   // unpackage
   if (SYS_UIDunpackUIDMessageXdr(interval, &status) == UID_FAILS)
   {
      SYS_UIDerrorMsg("UID Client: xdr unpack failed\n");
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_READ;
      return;
   }
 
   return;
}

/*******************************************************************************

Description: closes server connection


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void  SYS_UIDclient::CloseConnection(void)
{
   close(server_connection_id);
}

/*******************************************************************************

Description: closes failsafe server connection


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void  SYS_UIDclient::CloseFailsafeConnection(void)
{
   close(failsafe_server_connection_id);
}



