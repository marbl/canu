
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDclient.c,v $
$Revision: 1.1 $
$Date: 2004-09-23 20:32:58 $
$Name: not supported by cvs2svn $
$Author: mcschatz $
$Log: not supported by cvs2svn $
Revision 1.4  2004/09/10 12:31:43  mschatz
Add standard copyright notice

Revision 1.3  2004/09/09 22:39:03  mschatz
USE_SOAP_UID support

Revision 1.2  2004/06/25 04:05:04  ahalpern
Fixes to allow access to JTC uid server from linux

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

Revision 1.5  1999/10/13 19:02:56  sdmurphy
misc small changes for debugging

Revision 1.4  1999/07/14 17:24:33  stine
update_cds script was executed against these files.

Revision 1.3  1999/01/28 14:57:45  sdmurphy
added GetMaxUIDSize and QueryServer

Revision 1.2  1999/01/13 14:28:05  sdmurphy
version 0 prelim

Revision 1.1  1998/12/30 19:46:13  sdmurphy
Renamed uid_client.c to SYS_UIDclient.c

Revision 1.3  1998/12/21 19:04:04  sdmurphy
added support for uid incrementer

Revision 1.2  1998/12/18 18:02:50  sdmurphy
made set_UID_size return void type

Revision 1.1  1998/12/17 21:21:20  sdmurphy
Implements uid client API

**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#include "cds.h"
#include "SYS_UIDcommon.h"
#include "SYS_UIDerror.h"
#include "SYS_UIDclient.h"
#include "SYS_UIDclient_local.h"

#ifdef USE_SOAP_UID

#include "euidH.h"
#include "soapEUIDServerService.nsmap"
#include <assert.h>

// default server name
const char *TIGR_DefaultEuidServerNames = "tools.tigr.org:8190";
// pointer to runtime configuration
char * EuidServerNames = NULL;

#endif



/*******************************************************************************

Description: Utility func for setting interval message array.


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void SYS_UIDsetUIDSize(cds_uint64 block_size)
{
   size_UID      = block_size;
}


/*******************************************************************************

Description: Gets maximum UID size


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  SYS_UIDgetMaxUIDSize(cds_uint64* size)
{
   cds_uint64  interval_buffer[4];

   if (size == NULL)
      return UID_CODE_NULL_INTERVAL_PTR;
   /* Not implemented for JTC
   return QueryServer(UID_CODE_NEED_SIZE_INFO, size);
   */
   *size = JTC_MAX_REQUEST_SIZE; 
   return UID_CODE_OK;

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


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void Initialize(void)
{
  /* Not implemented for JTC 
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
  */
}

/*******************************************************************************

Description: External function for getting a new UID interval


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
static cds_int32 QueryServer(cds_int32 code, cds_uint64* interval)
{
   static char FirstTime     = 1;
   char        failsafe_flag = 0;
   CDS_UID_t   newBlockStart = 0;

   /* Not implemented for JTC

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

   */
  
   /* Begin new JTC section */
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
}

/***********************************************************************************/

size_t
JTC_GUIDWriteMemoryCallback(void *ptr, size_t size, size_t nmemb, void *data)
{
  register int realsize = size * nmemb;
  struct JTC_GUIDMemoryStruct *mem = (struct JTC_GUIDMemoryStruct *)data;
  
  mem->memory = (char *)realloc(mem->memory, mem->size + realsize + 1);
  if (mem->memory) {
    memcpy(&(mem->memory[mem->size]), ptr, realsize);
    mem->size += realsize;
    mem->memory[mem->size] = 0;
  }
  return realsize;
}


#ifdef USE_SOAP_UID
void SYS_UIDset_euid_server(const char * servers)
{
  EuidServerNames=strdup(servers);
  assert(EuidServerNames != NULL);
}
#endif

/***************************************************************************************/

CDS_UID_t getGUIDBlock(int guidRequestSize)
{
  CDS_UID_t guidStart = 0;

#ifndef USE_SOAP_UID
  CURL *curl_handle;
  char guidRequest[JTC_GUID_REQUEST_URL_MAX_SIZE];
  char httpResponse[JTC_GUID_HTTP_RESPONSE_MAX_SIZE];
  char guidNumResponse[JTC_GUID_NUM_BUFFER_SIZE]; /* needs to comfortably fit a long as a string */
  int guidNumLength = 0;
  struct JTC_GUIDMemoryStruct chunk;
  int i;
  int guidPositionStart = 0;
  int guidPositionEnd = 0;

  chunk.memory=NULL; /* we expect realloc(NULL, size) to work */
  chunk.size = 0;    /* no data at this point */
 
  curl_global_init(CURL_GLOBAL_ALL);
 
  /* init the curl session */
  curl_handle = curl_easy_init();
 
  /* specify URL to get */
  sprintf(guidRequest, "%s%d", JTC_GUID_URL, guidRequestSize);
  curl_easy_setopt(curl_handle, CURLOPT_URL, guidRequest);
 
  /* send all data to this function  */
  curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, JTC_GUIDWriteMemoryCallback);

  /* we pass our 'chunk' struct to the callback function */
  curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *)&chunk);
  
  /* get it! */
  curl_easy_perform(curl_handle);
 
  /* cleanup curl stuff */
  curl_easy_cleanup(curl_handle);

  if (chunk.size >= JTC_GUID_HTTP_RESPONSE_MAX_SIZE) {
    /* error */
    free(chunk.memory);
    return 0;
  }

  memcpy(httpResponse, chunk.memory, chunk.size);
  httpResponse[chunk.size] = '\0';

  /* HTTP response of this form is assumed:   */
  /* <html><title>SUCCESS</title><body>       */
  /* <h2>Guid Start:</h2>                     */
  /* 1089045040000                            */
  /* </body></html>                           */

  if (strncmp(httpResponse+13,"SUCCESS",7)==0) {
    guidPositionStart = findGuidStartFromHttpString(httpResponse);
    guidPositionEnd = findGuidEndFromHttpString(httpResponse);
    guidNumLength = guidPositionEnd - guidPositionStart;
    if (guidPositionStart == 0 || guidPositionEnd <= guidPositionStart) {
      /* error */
      free(chunk.memory);
      return 0;
    }
    memcpy(guidNumResponse, httpResponse + guidPositionStart, guidNumLength);
    for (i=guidNumLength;i<JTC_GUID_NUM_BUFFER_SIZE;i++) {
      guidNumResponse[i] = '\0';
    }
    guidStart = STR_TO_UID(guidNumResponse,NULL,10);
  } else {
    /* error */
    free(chunk.memory);
    return 0;
  }

  free(chunk.memory);
#else
  struct soap soap;
  struct impl__getEUIDBlockResponse euid;
  char dummy[256];
  char euidServerBuffer[2048];
  char *servers[10];
  int loop;
  int i;

  if (EuidServerNames == NULL)
  {
    // Have to copy because a static string lives in a read-only section on alpha
    strcpy(euidServerBuffer, TIGR_DefaultEuidServerNames);
  }
  else
  {
    strcpy(euidServerBuffer, EuidServerNames);
  }

  fprintf(stderr, "Parsing \"%s\"\n", euidServerBuffer);

  servers[0] = strtok(euidServerBuffer, ","); 
  if(servers[0] == NULL) 
  {
    SYS_UIDerrorMsg("EUIDService URL not specified\n");
    return UID_FAILS; // Error
  }

  fprintf(stderr, "servers[0]=\"%s\"\n", servers[0]);

  // parse rest of servers
  for(loop=1; loop < 10; loop++) 
  {
    servers[loop] = strtok(NULL, ",");
    if(servers[loop] == NULL) { break; }
    fprintf(stderr, "servers[%d]=\"%s\"\n", loop, servers[loop]);
  }

  fprintf(stderr, "Initializing soap... ");
  soap_init(&soap);
  fprintf(stderr, "ok.\n");

  for(i = 0; i < loop; i++) 
  {
    sprintf(dummy,"http://%s/axis/services/EUIDServer", servers[i]);
    fprintf(stderr, "Trying to contact %s\n", dummy);

    if (soap_call_impl__getEUIDBlock (&soap, dummy,	"", "TIGR_DEFAULT", 
                                      guidRequestSize, &euid ) == SOAP_OK) 
    {
      // got an euid
      guidStart = euid._getEUIDBlockReturn;
      fprintf(stderr, "EUID query suceeded -- starting EUID is "F_S64"\n", guidStart);
      break;
    }
  }
	  
  if(guidStart == 0)
  {
    // error condition
    soap_print_fault(&soap,stderr);
    SYS_UIDerrorMsg("EUIDService failure -- no servers available\n");
    return UID_FAILS;
  }
#endif
 
  return guidStart;
}

/***************************************************************************/

int findGuidStartFromHttpString(char* httpString) {
  /* the first position after the 2nd return character */
  int returnCount = 0;
  int i;
  for (i=0;i<JTC_GUID_NUM_BUFFER_SIZE;i++) {
    if (httpString[i] == '\n') {
      returnCount++;
    }
    if (returnCount == 2) {
      return i+1;
    }
  }
  /* error */
  return 0;
}

/***************************************************************************/

int findGuidEndFromHttpString(char* httpString) {
  /* the first position before the 3nd return character */
  int returnCount = 0;
  int i;
  for (i=0;i<JTC_GUID_NUM_BUFFER_SIZE;i++) {
    if (httpString[i] == '\n') {
      returnCount++;
    }
    if (returnCount == 3) {
      return i-1;
    }
  }
  /* error */
  return 0;
}

/*******************************************************************************

Description: Utility func for creating a socket connection to server


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


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


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
static void  CloseConnection(void)
{
   close(server_connection_id);
}

/*******************************************************************************

Description: Utility func for closing connection to failsafe server


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
static void  CloseFailsafeConnection(void)
{
   close(failsafe_server_connection_id);
}






