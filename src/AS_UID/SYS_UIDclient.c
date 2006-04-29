
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

#include "cds.h"
#include "SYS_UIDcommon.h"
#include "SYS_UIDerror.h"
#include "SYS_UIDclient.h"
#include "SYS_UIDclient_local.h"


/*******************************************************************************

Description: Utility func for setting interval message array.

*******************************************************************************/
void SetUIDInterval(cds_uint64 a, cds_uint64 a_size, cds_uint64 b, cds_uint64 b_size)
{
   interval_UID[0]        = a;
   interval_UID[1]        = a_size;
   interval_UID[2]        = b;
   interval_UID[3]        = b_size;
   increment_offset_UID   = 0L;
   increment_max_UID      = a_size + b_size;
}

/*******************************************************************************

Description: External function for setting current blocksize.

*******************************************************************************/
void SYS_UIDsetUIDSize(cds_uint64 block_size)
{
   size_UID      = block_size;
}


/*******************************************************************************

Description: Gets maximum UID size

*******************************************************************************/
cds_int32  SYS_UIDgetMaxUIDSize(cds_uint64* size)
{

   if (size == NULL)
      return UID_CODE_NULL_INTERVAL_PTR;
#ifndef NOT_IMPLEMENTED_JTC
   return QueryServer(UID_CODE_NEED_SIZE_INFO, size);
#else
   *size = UID_MAX_REQUEST_SIZE; 
   return UID_CODE_OK;
#endif
}


/*******************************************************************************

Description: 

   External function for getting the current increment UID. It returns the
   current value of the incrementer, and increments the value for next time.
   Note that even though it increments the value, it returns the
   pre-incremented value.

*******************************************************************************/
cds_int32   SYS_UIDgetNextUID(cds_uint64* uid)
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

*******************************************************************************/
cds_int32   SYS_UIDgetLastUID(cds_uint64* uid)
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

Description: External function for getting the current UID interval

*******************************************************************************/
cds_int32 SYS_UIDgetLastUIDInterval(cds_uint64* interval)
{
   if (interval == NULL)
      return UID_CODE_NULL_INTERVAL_PTR;

   if (status == UID_CODE_START)
      return status;

   interval[0] = interval_UID[0];
   interval[1] = interval_UID[1];
   interval[2] = interval_UID[2];
   interval[3] = interval_UID[3];

   return status;
}


/*******************************************************************************

Description: Initializes state of client on first NEW call

*******************************************************************************/
void Initialize(void)
{

#ifndef NOT_IMPLEMENTED_JTC
   char* port_string;
   char* host_string;
   char* failsafe_port_string;
   char* failsafe_host_string;

   SYS_UIDtype = UID_CLIENT_TYPE;
   
   // REGULAR SERVER **********************************************************
   if ( (port_string = getenv("SYS_UID_SERVER_PORT")) == NULL)
      server_tcp_port = UID_DEFAULT_SERVER_TCP_PORT;
   else
      server_tcp_port = atoi(port_string);

   if ( (host_string = getenv("SYS_UID_SERVER_HOST_NAME")) == NULL)
   {
      strcpy(server_host_name, UID_DEFAULT_SERVER_HOST_NAME);
   }
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

   // FAILSAFE SERVER **********************************************************
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
         SYS_UIDerrorMsg("UID Client: failsafe host name too long, using default \
failsafe server host name\n");
         strcpy(failsafe_server_host_name, UID_DEFAULT_FAILSAFE_SERVER_HOST_NAME);
      }   
   }
#endif
}

/*******************************************************************************

Description: External function for getting a new UID interval

*******************************************************************************/
cds_int32 SYS_UIDgetNewUIDInterval(cds_uint64* interval)
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

Description: Internal function for handling server communication

*******************************************************************************/
static cds_int32 QueryServer(cds_int32 code, cds_uint64* interval)
{
   static char FirstTime     = 1;
   char        failsafe_flag = 0;
   CDS_UID_t   newBlockStart = 0;

#ifndef NOT_IMPLEMENTED_JTC

   if (FirstTime)
   {
      Initialize();
      FirstTime = 0;
   }

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
         return UID_CODE_CREATE_CONN_FAILED;
   
      if (ConfigureFailsafeConnection() == UID_FAILS)
      {
         CloseFailsafeConnection();
         return UID_CODE_CONFIGURE_CONN_FAILED;
      }
      ReceiveFailsafeServerMessage(code, interval);
      CloseFailsafeConnection(); 
   }

#else
  
   newBlockStart = getGUIDBlock(size_UID);
   if (newBlockStart > 0) {
     code = UID_CODE_OK;
     interval[0] = newBlockStart;
     interval[1] = size_UID;
     interval[2] = 0L;
     interval[3] = 0L;
     status = UID_CODE_OK;
   } else {
     status = UID_CODE_CANT_CONNECT;
   }
   return status;

#endif
}

 

/*******************************************************************************

Description: Utility func for creating a socket connection to server

*******************************************************************************/
static cds_int32  CreateConnection(void)
{
  struct timeval send_time_out;
  struct timeval recv_time_out;
  char logmessage[200];

   /* obtain socket */
   if ( (server_connection_id = socket(AF_INET, SOCK_STREAM, 0)) < 0)
   {
      SYS_UIDperr("UID Client::create_connection: could not open socket\n");
      return UID_FAILS;
   }
   /* set socket timeout characteristics */
   send_time_out.tv_sec = UID_CLIENT_SEND_TIMEOUT;
   send_time_out.tv_usec = 0;
   if (setsockopt(server_connection_id,
                    SOL_SOCKET,
                    SO_SNDTIMEO,
                    (const void *)(&send_time_out),
                    sizeof(send_time_out)) < 0)
   {
      SYS_UIDperr("UID Client::create_connection: could not set socket timeout options\n");
      sprintf(logmessage,"client: could not set socket timeout options\n");
      SYS_UIDlogMessage(logmessage);
      return UID_FAILS;
   }
   recv_time_out.tv_sec = UID_CLIENT_RECV_TIMEOUT;
   recv_time_out.tv_usec = 0;
   if (setsockopt(server_connection_id,
                    SOL_SOCKET,
                    SO_RCVTIMEO,
                    (const void *)(&recv_time_out),
                    sizeof(recv_time_out)) < 0)
   {
      SYS_UIDperr("UID Client::create_connection: could not set socket timeout options\n");
      sprintf(logmessage,"client: could not set socket timeout options\n");
      SYS_UIDlogMessage(logmessage);
      return UID_FAILS;
   }
   return UID_OK;
}
/*******************************************************************************

Description: Utility func for creating a socket connection to failsafe server

*******************************************************************************/
static cds_int32  CreateFailsafeConnection(void)
{
   if ( (failsafe_server_connection_id = socket(AF_INET, SOCK_STREAM, 0)) < 0)
   {
      SYS_UIDperr("UID Client::CreateFailsafeConnection: could not open socket\n");
      return UID_FAILS;
   }
   return UID_OK;
}


/*******************************************************************************

Description: Utility func for setting up server host info

*******************************************************************************/
static cds_int32  GetServerHostInfo(void)
{
   // uses the member var 'server_host_info' of type struct hostent
   if ( (server_host_info = gethostbyname(server_host_name)) == NULL)
   {
      sprintf(SYS_UIDerr_str,"UID Client: could not find address for host %s \n",
	      server_host_name);
      SYS_UIDperr(SYS_UIDerr_str);
      return UID_FAILS;
   }
   return UID_OK;
}

/*******************************************************************************

Description: Utility func for setting up failsafe server host info

*******************************************************************************/
static cds_int32  GetFailsafeServerHostInfo(void)
{
   // uses the member var 'failsafe_server_host_info' of type struct hostent
   if ( (failsafe_server_host_info = gethostbyname(failsafe_server_host_name)) == NULL)
   {
      sprintf(SYS_UIDerr_str,"UID Client: could not find address for failsafe host %s \n",
	      failsafe_server_host_name);
      SYS_UIDperr(SYS_UIDerr_str);
      return UID_FAILS;
   }
   return UID_OK;
}


/*******************************************************************************

Description: Utility func for configuring socket to server

*******************************************************************************/
static cds_int32  ConfigureConnection(void)
{
   if (5000 > server_tcp_port || server_tcp_port > 65535)
      {
      SYS_UIDerrorMsg("UID server port not in valid range\n");
      return UID_FAILS;
      }
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

Description: Utility func for configuring socket to failsafe server

*******************************************************************************/
static cds_int32  ConfigureFailsafeConnection(void)
{
   if (5000 > failsafe_server_tcp_port || failsafe_server_tcp_port > 65535)
      {
      SYS_UIDerrorMsg("UID failsafe server port not in valid range\n");
      return UID_FAILS;
      }
   if (GetFailsafeServerHostInfo() == UID_FAILS)
      return UID_FAILS;
   bzero( (char*) &failsafe_server_connection_data, 
	  sizeof(failsafe_server_connection_data) );
   failsafe_server_connection_data.sin_family = AF_INET;
   failsafe_server_connection_data.sin_addr.s_addr = 
      inet_addr( inet_ntoa( *((struct in_addr*)
			      (failsafe_server_host_info->h_addr_list[0])) ) );
   failsafe_server_connection_data.sin_port = htons(failsafe_server_tcp_port);
   failsafe_server_connection_data_size = sizeof(failsafe_server_connection_data);
   return UID_OK;
}


/*******************************************************************************

Description: Function for retrieving UID message from server

*******************************************************************************/
static void  ReceiveServerMessage(cds_int32 code, cds_uint64* interval)
{
   char writebuffer[12];
   char logmessage[200];

   if(connect(server_connection_id, (struct sockaddr*) &server_connection_data, 
      server_connection_data_size) < 0)
      {
      sprintf(logmessage,"client: error %d on socket connect", errno);
      SYS_UIDlogMessage(logmessage);
      // no error message here b/c want to default to failsafe if not available
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_CONNECT;
      return;
      }

   SYS_UIDpackUIDRequestXdr(writebuffer, code, size_UID);
   
   if (SYS_UIDwriten(server_connection_id, (char*)writebuffer, 12) == UID_FAILS)
   {
      sprintf(logmessage,"client: error %d on socket writen", errno);
      SYS_UIDlogMessage(logmessage);
      SYS_UIDperr("UID Client::ReceiveServerMessage: writen failed\n");
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_READ;
      return;
   }

   if (SYS_UIDreadn(server_connection_id, (char*)SYS_UIDmessage_array, 
		    UID_MESSAGE_SIZE) == UID_FAILS)
   {
      sprintf(logmessage,"client: error %d on socket readn", errno);
      SYS_UIDlogMessage(logmessage);
      SYS_UIDperr("UID Client::ReceiveServerMessage: readn failed\n");
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
      SYS_UIDperr("UID Client::ReceiveServerMessage: xdr unpack failed\n");
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

Description: Function for retrieving UID message from failsafe server

*******************************************************************************/
static void  ReceiveFailsafeServerMessage(cds_int32 code, cds_uint64* interval)
{
   char writebuffer[12];
   char logmessage[200];

   if(connect(failsafe_server_connection_id, 
	      (struct sockaddr*) &failsafe_server_connection_data, 
      failsafe_server_connection_data_size) < 0)
      {
      sprintf(logmessage,"client: error %d on failsafe socket connect", errno);
      SYS_UIDlogMessage(logmessage);
      SYS_UIDperr("UID Client: can't connect to failsafe server\n");
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_CONNECT;
      return;
      }

   SYS_UIDpackUIDRequestXdr(writebuffer, code, size_UID);
   
   if (SYS_UIDwriten(failsafe_server_connection_id, 
		     (char*)writebuffer, 12) == UID_FAILS)
   {
      sprintf(logmessage,"client: error %d on failsafe socket writen", errno);
      SYS_UIDlogMessage(logmessage);
      SYS_UIDperr("UID Client:send_server_message: writen failed to failsafe server\n");
      interval[0] = 0L;
      interval[1] = 0L;
      interval[2] = 0L;
      interval[3] = 0L;
      status = UID_CODE_CANT_READ;
      return;
   }

   if (SYS_UIDreadn(failsafe_server_connection_id, (char*)SYS_UIDmessage_array, 
		    UID_MESSAGE_SIZE) == UID_FAILS)
   {
      sprintf(logmessage,"client: error %d on failsafe socket readn", errno);
      SYS_UIDlogMessage(logmessage);
      SYS_UIDperr("UID Client: readn failed from failsafe server\n");
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
      SYS_UIDperr("UID Client: xdr unpack failed\n");
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

Description: Utility func for closing connection

*******************************************************************************/
static void  CloseConnection(void)
{
   close(server_connection_id);
}

/*******************************************************************************

Description: Utility func for closing connection to failsafe server

*******************************************************************************/
static void  CloseFailsafeConnection(void)
{
   close(failsafe_server_connection_id);
}






