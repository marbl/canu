
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDserver.c,v $
$Revision: 1.5 $
$Date: 2005-06-16 20:28:05 $
$Name: not supported by cvs2svn $
$Author: brianwalenz $
$Log: not supported by cvs2svn $
Revision 1.4  2005/03/30 21:11:56  eliv
Fix a couple of gcc warnings.

Revision 1.3  2005/03/22 19:49:28  jason_miller
The TIGR tip as of March 22 2005. Commit by Jason Miller at TIGR.

Revision 1.2  2004/09/10 12:31:43  mschatz
Add standard copyright notice

Revision 1.1  2004/06/24 12:51:06  mpop
Added AS_UID

Revision 1.2  2003/05/09 21:04:03  mpop
Dos2unixed all files.
Modified c_make.as to set SEP_PATH relative to LOCAL_WORK

Revision 1.1.1.1  2003/05/08 18:40:12  aaronhalpern
versions from TIGR

Revision 1.2  2001/09/25 23:03:20  mpop
Dos2Unixed

Revision 1.1.1.1  2001/09/25 20:21:05  mpop
Celera Assembler

Revision 1.6  2000/03/02 17:48:48  sdmurphy
added full path warning

Revision 1.5  1999/10/13 19:02:57  sdmurphy
misc small changes for debugging

Revision 1.4  1999/07/14 17:24:34  stine
update_cds script was executed against these files.

Revision 1.3  1999/01/28 15:01:31  sdmurphy
added support for GetMaxUIDSize and beefed up messages

Revision 1.2  1999/01/13 14:29:45  sdmurphy
version 0 prelim

Revision 1.1  1998/12/30 19:36:50  sdmurphy
Renamed uid_server.c SYS_UIDserver.c

Revision 1.3  1998/12/21 16:42:39  sdmurphy
added daemon logging capability

Revision 1.2  1998/12/18 18:25:19  sdmurphy
removed vestigial var us_current_uid_block_size

Revision 1.1  1998/12/17 21:23:52  sdmurphy
Implements uid server

**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#include "cds.h"
#include "SYS_UIDcommon.h"
#include "SYS_UIDerror.h"
#include "SYS_UIDserver.h"
#include "SYS_UIDserver_local.h"

/*******************************************************************************

Description: prints usage if args don't match


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void Usage(cds_int32 argc, char** argv)
{

printf("\n                                                                 \n");
printf("Use -i to initialize a UID position file. Should theoretically only\n");
printf("be done ONCE when a machine is configured. Should NOT be re-run if \n");;
printf("the machine is turned off or server re-booted or may risk          \n");
printf("corrupting UID uniqueness.                                         \n");
printf("                                                                   \n");
printf("usage: %s -i <pos_filename> <initial_uid>                          \n", argv[0]);
printf("                                                                   \n");
printf("Use -r in the machine boot script to start the server daemon       \n");
printf("                                                                   \n");
printf("usage: %s -r                                                       \n", argv[0]);
printf("   <pos_filename> <start#> <size> <max_block> <update_freq> <port> \n");
/* Non implemented for JTC
printf("   [-ma <hour interval> <address list for email on all requests>]  \n");
printf("   [-me <hour interval> <address list for email on errors>]        \n");
printf("   [-mf <address list for email on fatal error>]                   \n");
*/
printf("   [-d debug]                                                      \n");
printf("                                                                   \n");
printf("Use -k to kill a uid_server currently running in the background    \n");
printf("                                                                   \n");
printf("usage: %s -k <port of server to kill>                            \n\n", argv[0]);

   exit(0);
}


/*******************************************************************************

Description: Parses command line options


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void SYS_UIDparseOptions(cds_int32 argc, char** argv)
{
   cds_int32 i, j, k;
   cds_int32 address_length;

   /* initial argc check */
   if (argc < 3) 
     Usage(argc, argv);

   // check for kill option
   kill_option_invoked_flag = 0;
   for (i=0;i<argc;i++)
      if (strcmp(argv[i],"-k") == 0)
      {
  	 // if -k invoked without a port number, give usage info
         if (i == (argc-1))
            Usage(argc, argv);

	 // use the 'tcp_port' field for the server to kill
         tcp_port = atoi(argv[i+1]);         
         
         // set flag and return regardless of other options
         kill_option_invoked_flag = 1;
         return;
      }

   // check for debug option
   SYS_UIDdebug_flag = 0;
   for (i=0;i<argc;i++)
      if (strcmp(argv[i],"-d") == 0)
      {
         SYS_UIDdebug_flag = 1;
         SYS_UIDerrorMsg("Debug flag ON\n");
      }
   /* Not implemented for JTC
   // ma alert option
   FillAddressList("-ma", &ma_address_list, argc, argv, &ma_alert_interval); 
   if (SYS_UIDdebug_flag == 1)
   {
       if (ma_address_list != NULL)
       {
          sprintf(SYS_UIDerr_str, "sending -ma alerts every %d hours to:%s\n", 
		  ma_alert_interval, ma_address_list);
          SYS_UIDerrorMsg(SYS_UIDerr_str);
       }
   }

   // me alert option
   FillAddressList("-me", &me_address_list, argc, argv, &me_alert_interval); 
   if (SYS_UIDdebug_flag == 1)
   {
       if (me_address_list != NULL)
       {
          sprintf(SYS_UIDerr_str, "sending -me alerts every %d hours to:%s\n", 
		  me_alert_interval, me_address_list);
          SYS_UIDerrorMsg(SYS_UIDerr_str);
       }
   }

   // mf alert option
   FillAddressList("-mf", &mf_address_list, argc, argv, NULL); 
   if (SYS_UIDdebug_flag == 1)
   {
       if (mf_address_list != NULL)
       {
          sprintf(SYS_UIDerr_str, "sending -mf alerts to:%s\n", mf_address_list);
          SYS_UIDerrorMsg(SYS_UIDerr_str);
       }
   }
   */

}

/*******************************************************************************

Description: Initializes alerts


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void FillAddressList(const char*    option, 
                     char**         list, 
                     cds_int32          argc, 
                     char**         argv,
                     cds_int32*         interval)
{
   cds_int32 i = 0;
   cds_int32 j = 0;
   cds_int32 k = 0;
   cds_int32 address_length = 0;
   cds_int32 interval_offset = 0;

   if (interval == NULL)
      interval_offset = 2;
   else
      interval_offset = 3;

   // start with an empty list
   *list = NULL;

   // check for option
   for (i=0;i<argc;i++)
   {
      if (strcmp(argv[i], option) == 0)
      {
   	 // check for up to two args, one numeric and one address
         if (i+interval_offset > argc)
            Usage(argc, argv);

	 // find the last argument to belong to the address list
         // ends with either - or runs out of arguments
         for(j=i+interval_offset;j<argc;j++)
            if (argv[j][0] == '-') 
               break;

         // get alert timer
         if (interval != NULL)
            *interval = atoi(argv[i+1]);

         // count space to insert the remaining args into ma address list
         address_length = 0;
         for(k=i+(interval_offset-1);k<j;k++)
            address_length += strlen(argv[k]);
         if (address_length == 0)
            Usage(argc, argv);
         *list = (char*)calloc((address_length + 
			(j - (i+(interval_offset-1))) + 2), sizeof(char));

         // insert addresses
         sprintf(*list,"");
         address_length = 0;
         for(k=i+(interval_offset-1);k<j;k++)
         {
            strcat(*list, argv[k]);
            strcat(*list, " ");
	 }
      }
   }
}

/*******************************************************************************

Description: Initializes alerts by obtaining the current system time


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void InitializeAlerts(void)
{
  ma_alert_time = 0;
  me_alert_time = 0;
}

/*******************************************************************************

Description: Activates MA alert if time requirement is met


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
/* Not implemented for JTC
void TriggerMaAlert(void)
{
   cds_int32 hours_elapsed;
   cds_int64 current_time;
   struct timespec tp;
   char emessage[800];

   getclock(TIMEOFDAY, &tp);
   current_time = tp.tv_sec;
   hours_elapsed = (current_time - ma_alert_time)/3600;

   if (hours_elapsed >= ma_alert_interval)
     {
     sprintf(emessage, "received UID request", tcp_port); 

     SendMailList("received UID request",
                  emessage,
		  ma_address_list);

     ma_alert_time = current_time;
     }
}
*/

/*******************************************************************************

Description: Activates ME alert if time requirement is met


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
/* Not implemented for JTC
void TriggerMeAlert(void)
{
   cds_int32 hours_elapsed;
   cds_int64 current_time;
   struct timespec tp;
   char emessage[800];

   getclock(TIMEOFDAY, &tp);
   current_time = tp.tv_sec;
   hours_elapsed = (current_time - me_alert_time)/3600;

   if (hours_elapsed >= me_alert_interval)
     {
     sprintf(emessage,"nonfatal UID server error", tcp_port);
     SendMailList("nonfatal UID server error", 
                  emessage,
		  ma_address_list);
     me_alert_time = current_time;
     }
}
*/

/*******************************************************************************

Description: Activates MF alert - always sends messages


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
/* Not implemented for JTC
void TriggerMfAlert(void)
{
   char emessage[800];

   sprintf(emessage,"fatal UID server error", tcp_port);
   SendMailList("fatal UID server error", emessage, mf_address_list);
}
*/

/*******************************************************************************

Description: Sends mail to address list with subject


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
/* Not implemented for JTC
void SendMailList(const char* subject, const char* emessage, const char* address_list)
{
   char sendbuffer[1000];

   if ( (strlen(subject) + +strlen(emessage) + strlen(address_list)) > 900)
      {
      SYS_UIDerrorMsg("UID Server: address list too long to email");
      return;
      }

   sprintf(sendbuffer, "echo \"%s\" | mailx -s \"%s\" %s ", 
	   emessage, subject, address_list);
   system(sendbuffer);
}
*/

/*******************************************************************************

Description: Sets the UID status and manages the unrecoverable_error_flag


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void SetStatus(cds_int32 status, char flag)
{
   status_code = status;
   if (unrecoverable_error_flag == UID_OK) // this is a one-way switch
      unrecoverable_error_flag = flag;
}

/*******************************************************************************

Description: initializes vars and parses arguments, branching appropriately


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDserverInitialize(cds_int32 argc, char** argv)
{
   // initialization of statics
   positionfile_name           = NULL;
   index_size                  = 0L;  
   index_update_increment      = 0L;  
   max_block_size              = 0L;  
   status_code                 = 0;   
   client_connection_id        = 0;   
   connection_id               = 0;   
   current_position_UID        = 0L;  
   index_start_UID             = 0L;  
   current_UID                 = 0L;  
   unrecoverable_error_flag    = UID_OK;

   // clear interval
   SetUIDInterval(0L, 0L, 0L, 0L);

   if (kill_option_invoked_flag == 1)
   {
      if (IssueKillSignal() == UID_OK)
      {
         exit(0);
      }
      else
      {
	 sprintf(SYS_UIDerr_str, 
		 "Could not connect to server at port %d to kill\n", tcp_port);
	 SYS_UIDerrorMsg(SYS_UIDerr_str);
         exit(1);
      }
   }

   //  check argument count
   if (argc < 4) {
     SYS_UIDerrorMsg("UID Server: incorrect argument usage...quiting\n");
     exit(0);
   }

   // Handle argument possibilities:

   // RUN SERVER ///////////////////////////////////////////////////
   if (argc >= 8 && argv[1][1] == 'r')
   {
      positionfile_name      =         argv[2];      
      index_start_UID        = strtoul(argv[3], (char**)NULL, 10);
      index_size             = strtoul(argv[4], (char**)NULL, 10);
      max_block_size         = strtoul(argv[5], (char**)NULL, 10);
      index_update_increment = strtoul(argv[6], (char**)NULL, 10);
      tcp_port               =    atoi(argv[7]);

      // very first thing, check the port for availability to avoid any other
      // server running on this port - don't want to overwrite this 
      // position file while it's running.
      sprintf(SYS_UIDerr_str,"Starting UID server on port %d...\n", tcp_port);
      SYS_UIDerrorMsg(SYS_UIDerr_str);

      if (AttemptPortInitialization() == UID_FAILS) {
        SYS_UIDerrorMsg("Port initialization failed...exiting\n");
        CloseConnection();
        exit(1);
      }

      // initialize the UID status
      SetStatus(UID_CODE_OK, UID_OK);
      InitializeAlerts();

      return UpdateFromPositionFile();
   }

   // INITIALIZE POSITION FILE //////////////////////////////////////
   else if (argc == 4 && argv[1][1] == 'i')
   {
      positionfile_name      = argv[2];
      index_start_UID        = strtoul(argv[3], (char**)NULL, 10);

      if (PositionWrite(index_start_UID) == UID_OK)
         exit(0);

      exit(1);
   }
   else {
     SYS_UIDerrorMsg("UID server: unable to parse arguments...quiting\n");
     exit(0);
   }
   
   return UID_OK;
}

/*******************************************************************************

Description: 

   Writes current_position_UID to position file.
 
Input:


Output:


Returns:


Globals:


Notes:

   Does not change the value: this value must be incremented before this 
   function is called if this is the desired effect.

*******************************************************************************/
cds_int32  PositionWrite(cds_uint64 position)
{
   FILE* write_fp;

   if ( (write_fp=fopen(positionfile_name, "wb")) == NULL)
   {
      sprintf(SYS_UIDerr_str,"UID Server: could not open position file %s to write - check for full path usage\n",
	      positionfile_name);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      SetStatus(UID_CODE_POS_CONFIG_ERROR, UID_FAILS);
      return UID_FAILS;
   }
   if (SYS_UIDdebug_flag == 1)
   {
      sprintf(SYS_UIDerr_str,"UID server: updating position file to %ld\n",
	      current_position_UID);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
   }    
   fwrite(&position, sizeof(cds_uint64), 1, write_fp);
   fclose(write_fp);
   return UID_OK;
}

/*******************************************************************************

Description: 

   Reads the UID_current_pos_value from the position file.
   Should only be called once per program invocation.

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  PositionRead(void)
{
   FILE* read_fp;

   if ( (read_fp=fopen(positionfile_name, "rb")) == NULL)
   {
      sprintf(SYS_UIDerr_str,"UID Server: could not open position file %s to read\n",
	      positionfile_name);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      SetStatus(UID_CODE_POS_CONFIG_ERROR, UID_FAILS);
      return UID_FAILS;
   }
   fread(&current_UID, sizeof(cds_uint64), 1, read_fp);
   if (SYS_UIDdebug_flag == 1)
   {
      sprintf(SYS_UIDerr_str,"UID Server: read position value %ld\n",
	      current_UID);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
   }
   if (feof(read_fp) != 0)
   {
      sprintf(SYS_UIDerr_str,"UID Server: position file %s not correct format\n",
	      positionfile_name);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      SetStatus(UID_CODE_POS_CONFIG_ERROR, UID_FAILS);
      return UID_FAILS;
   }
   fclose(read_fp);
   return UID_OK;
}

/*******************************************************************************

Description: 

   This function should only be called once, during the initialization
   of the server.

   This function must read from the position file the UID number that
   is guaranteed to be the highest one issued by the server, such that
   if we start with one more than this number, the continued uniqueness
   of the UID stream is guaranteed.

   If there are problems setting up UID_current, the uid_status_code is
   set appropriately.


Input:


Output:

   Returns UID_OK or UID_FAILS

Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  UpdateFromPositionFile(void)
{

   if (PositionRead() == UID_FAILS)
      return UID_FAILS;

   if (UIDIsValid(0L) == UID_FAILS)
   {
      // override standard UIDIsValid message to be a pos config err
      SetStatus(UID_CODE_POS_CONFIG_ERROR, UID_OK);
      return UID_FAILS;
   }

   // Now the UID in the position file has passed the error checks - and because we are
   // about to use it, it is now invalid. Therefore, we will update its value.
   // Note that this is happening before we actually use it, so that if the machine
   // is killed as this write takes place we are still OK in terms of guaranteeing
   // UID validity.

   return IncrementUpdatePositionFile();
}

/*******************************************************************************

Description: 

   Updates the position file. 

Input:


Output:


Returns:


Globals:


Notes:

   This keys off of the current_UID value to ensure that even if a
   request is well over the room left in the previous block, the
   new position always starts index_update_increment above the 
   first block that surpasses the previous block end. In other words,
   blocks are not contiguous, but have gaps proportional to the
   requests that trigger their release.

*******************************************************************************/
cds_int32 IncrementUpdatePositionFile()
{
   current_position_UID = current_UID + index_update_increment;
   if (PositionWrite(current_position_UID) == UID_FAILS)
      return UID_FAILS;
   return UID_OK;
}

/*******************************************************************************

Description: 

   This is the  service body and indefinite loop.

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDserverStart(void)
{
  /* This blocking function runs the server loop */

   while(1)
   {
      if (AcceptClientConnection() == UID_OK) 
      {
	/* Not implemented for JTC
         if (ma_address_list != NULL)
            TriggerMaAlert();
         if (SendClientMessage() == UID_FAILS)
	    {
            if (mf_address_list != NULL)
               TriggerMfAlert();
            SYS_UIDerrorMsg("UID Server: ERROR - must be killed and re-booted!\n"); 
            CloseConnection();
            return UID_FAILS;
            }
	*/
      }
      CloseConnection();
      /* Not implemented for JTC
     if (status_code != UID_CODE_OK && me_address_list != NULL)
        TriggerMeAlert();
      */
   }
   return UID_FAILS; // should never reach
}

/*******************************************************************************

Description: 

   Reads the request message from UID clients

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void  ReadClientRequest(cds_int32* client_status, cds_uint64* request_size)
{

   char  readbuffer[12];
   char logmessage[200];

   // read client request block size
   if (SYS_UIDreadn(client_connection_id, (char*)readbuffer, 12) == UID_FAILS)
   {
      sprintf(SYS_UIDerr_str, "UID Server: error %d on socket readn\n", errno);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      SetStatus(UID_CODE_REQUEST_ERROR, UID_OK);
      SetUIDInterval(0L, 0L, 0L, 0L);
      return;
   }

   // unpack request size
   if (SYS_UIDunpackUIDRequestXdr(readbuffer, client_status, request_size) == UID_FAILS)
      {
      SetStatus(UID_CODE_REQUEST_ERROR, UID_OK);
      SetUIDInterval(0L, 0L, 0L, 0L);
      return;
      }

   if (SYS_UIDdebug_flag == 1)
      {
      sprintf(SYS_UIDerr_str,"UID Server: received client request for code %d blocksize %ld\n",
	      *client_status, *request_size);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      }

   return;
}

/*******************************************************************************

Description: 

   Writes the request message to UID clients

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void  WriteClientRequest(void)
{
  char logmessage[200];

   // pack client message into the UID_message_array
   if (SYS_UIDpackUIDMessageXdr(interval_UID, status_code) == UID_FAILS)
   {
      SetStatus(UID_CODE_SEND_ERROR, UID_OK);
      SetUIDInterval(0L, 0L, 0L, 0L);
   }

   // send to client - attempt even if the pack failed to prevent client hang
   if (SYS_UIDwriten(client_connection_id, (char*)SYS_UIDmessage_array, 
		     UID_MESSAGE_SIZE) == UID_FAILS)
   {
      sprintf(SYS_UIDerr_str,"UID Server: error %d on socket writen\n", errno);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      SetStatus(UID_CODE_SEND_ERROR, UID_OK);
      SetUIDInterval(0L, 0L, 0L, 0L);
   }

   return;
}

/*******************************************************************************

Description: 

   Prints message debug information in hex rows and columns.

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void  DebugClientMessage(void)
{
   cds_int32     i;
   cds_int32     j;
   char      num_buffer[50];
   cds_int32     length = 0;

   // for each of the 5 lines it takes to print an outgoing message...
   for (i=0; i<5 ; i++)   
   {
      // reset the byte offset of the message line
      j = 0;
      // start the line by printing the first number
      sprintf(SYS_UIDerr_str, "UID server, debug client: wrote %2x ", 
         (unsigned char)SYS_UIDmessage_array[i*8 + j]);
      // the amount of space in SYS_UIDerr_str used so far for this line
      length = strlen(SYS_UIDerr_str);
      // for this line, print the remaining values written
      for (j=1 ; j<8 ; j++)
      {
  	 // if array position still within the message size
         if ((i*8 + j) < UID_MESSAGE_SIZE)
         {
            sprintf(num_buffer, "%2x ", 
		    (unsigned char)SYS_UIDmessage_array[i*8 + j]);
            length += strlen(num_buffer);
            if (length >= (UID_ERR_STR_SIZE - 10))
               return;
            strcat(SYS_UIDerr_str, num_buffer);
         }
      }
      strcat(SYS_UIDerr_str,"\n");
      SYS_UIDerrorMsg(SYS_UIDerr_str);
   }
}

/*******************************************************************************

Description: 

   Processes a UID block request from the client.

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  SendClientMessage(void)
{
   cds_uint64 request_size;
   cds_int32  client_status;

   // safeguard against any client communication in fatal error state
   if (unrecoverable_error_flag != UID_OK)
      return UID_FAILS;

   if (SYS_UIDdebug_flag == 1)
   {
      sprintf(SYS_UIDerr_str,"UID server: connected to client...\n");
      SYS_UIDerrorMsg(SYS_UIDerr_str);
   }

   ReadClientRequest(&client_status, &request_size);

   switch(client_status)
   {
   case UID_CODE_SERVER_KILL:
      ShutdownServer();
      break;
   case UID_CODE_NEED_SIZE_INFO:
      // set the status_code to OK and use
      // the first field of the interval to
      // return the value of max_block_size
      status_code = UID_CODE_OK;
      SetUIDInterval(max_block_size, 0L, 0L, 0L);
      break;
   default:
      if (status_code == UID_CODE_OK)
         GetNextUID(request_size);
      // if status_code is not OK, then do not
      // deliver another chunk of id space, but
      // send whatever the server status is as a 
      // default - do not keep the client
      // hanging
      break;
   }

   WriteClientRequest();

   /*   if (SYS_UIDdebug_flag == 1)
        DebugClientMessage();           UNCOMMENT TO SEE INFO IN LOG */

   // handle error possibility during execution of this function
   if (unrecoverable_error_flag != UID_OK)
      return UID_FAILS;

   return UID_OK;
}

/*******************************************************************************

Description: 

   closes connection with current client

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void CloseConnection(void)
{
   close(client_connection_id);
}

/*******************************************************************************

Description: 

   This function checks whether the desired blocksize would be out of
   bounds. If not, it leaves state vars alone such that future 
   requests that are within bounds will work.

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  UIDIsValid(cds_uint64 block_size)
{
   cds_uint64 end_of_block;
   cds_uint64 max_ul = CDS_UINT64_MAX;


   // do 64-bit bounds check
   if ((max_ul - current_UID) <= block_size)
   {
      sprintf(SYS_UIDerr_str,"UID Server: UID 64-bit overflow\n");
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      SetStatus(UID_CODE_POS_BOUNDS_ERROR, UID_OK);
      return UID_FAILS;     
   }

   end_of_block = current_UID + block_size;

   // check to see that current UID request is within allowable limits
   if (end_of_block > (index_start_UID + index_size) || 
       end_of_block < index_start_UID)
   {
      sprintf(SYS_UIDerr_str,"UID Server: %ld is not between %ld and %ld\n",
         end_of_block, index_start_UID, index_start_UID + index_size);
      SYS_UIDerrorMsg(SYS_UIDerr_str);
      SetStatus(UID_CODE_RAN_OUT_OF_SPACE, UID_OK);
      return UID_FAILS;
   }
   return UID_OK;
}   


/*******************************************************************************

Description: 

   opens server socket for business

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  CreateConnection(void)
{
  struct timeval send_time_out;
  struct timeval recv_time_out;
  char logmessage[200];

   if ( (connection_id = socket(AF_INET, SOCK_STREAM, 0)) < 0)
   {
      SYS_UIDerrorMsg("UID Server: could not open socket\n");
      SYS_UIDhandleCreateError(errno);
      return UID_FAILS;
   }
   /* set socket timeout characteristics */
   send_time_out.tv_sec = UID_SERVER_SEND_TIMEOUT;
   send_time_out.tv_usec = 0;
   if (setsockopt(connection_id,
                    SOL_SOCKET,
                    SO_SNDTIMEO,
                    (const void *)(&send_time_out),
                    sizeof(send_time_out)) < 0)
   {
     SYS_UIDerrorMsg("UID Server::create_connection: could not set socket timeout options\n");
     return UID_FAILS;
   }
   recv_time_out.tv_sec = UID_SERVER_RECV_TIMEOUT;
   recv_time_out.tv_usec = 0;
   if (setsockopt(connection_id,
                    SOL_SOCKET,
                    SO_RCVTIMEO,
                    (const void *)(&recv_time_out),
                    sizeof(recv_time_out)) < 0)
   {
      SYS_UIDerrorMsg("UID Server::create_connection: could not set socket timeout options\n");
      return UID_FAILS;
   }
   return UID_OK;
}

/*******************************************************************************

Description: 

   configures server socket

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  RegisterConnection(void)
{

   struct sockaddr_in   connection_data;       
   cds_int32                  connection_data_size;  

   if (5000 > tcp_port || tcp_port > 65535)
      return UID_FAILS;

   // configure the server socket address structure to be a local UNIX stream socket
   bzero((char*) &connection_data, sizeof(connection_data));
   connection_data.sin_family              = AF_INET;
   connection_data.sin_port                = htons(tcp_port);
   connection_data.sin_addr.s_addr         = htonl(INADDR_ANY);
   connection_data_size                    = sizeof(connection_data);

   // bind the socket to the OS
   if ( bind(connection_id, (struct sockaddr *) &connection_data, 
      connection_data_size) < 0)
   {
      SYS_UIDhandleRegisterError(errno);
      return UID_FAILS;
   }
   return UID_OK;
}

/*******************************************************************************

Description: 

   Activates the server socket as listening.

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  ActivateConnection(void)
{
   // listen - activate the socket - use 1/2 the maximum buffer size.
   // The max size for the Dec Alpha is 1024 as of revision 4.2.34.5 
   //  of socket.h from 10/97
   if (listen(connection_id, SOMAXCONN/2) < 0)
   {
      SYS_UIDhandleActivateError(errno);
      return UID_FAILS;
   }
   return UID_OK;
}

/*******************************************************************************

Description: 

   Blocks for a client to connect.


Input:


Output:


Returns:


Globals:


Notes:

   This is the one and only place where the status is reset.


*******************************************************************************/
cds_int32  AcceptClientConnection(void)
{
   struct sockaddr_in   client_connection_data;
   cds_int32                client_connection_data_size;
   static char          FirstTime = 1;

   if (FirstTime)
   {
      FirstTime = 0;
      client_connection_data_size = sizeof(client_connection_data);
   }
   // wait for client connection
   client_connection_id = accept(connection_id, 
                                 (struct sockaddr *) &client_connection_data, 
                                 &client_connection_data_size);

   // reset status - the one and only place!
   SetStatus(UID_CODE_OK, UID_OK);

   // do error checking
   if (client_connection_id < 0)
   {
      SetStatus(UID_CODE_ACCEPT_CONN_FAILED, UID_OK);
      SYS_UIDhandleAcceptError(errno);
      return UID_FAILS;
   }
   return UID_OK;
}


/*******************************************************************************

Description: 

   Utility func for setting interval message array


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void SetUIDInterval(cds_uint64 a, cds_uint64 a_size, 
                         cds_uint64 b, cds_uint64 b_size)
{
   interval_UID[0] = a;
   interval_UID[1] = a_size;
   interval_UID[2] = b;
   interval_UID[3] = b_size;
}

/*******************************************************************************

Description: 

   This is where the work is done to communicate with the client.

   This function must do several things:

   (1) serve the next UID
   (2) if necessary, update the position file
   (3) decide if there are any errors the client should know about

   This function updates, directly or indirectly, global vars:

     current_UID_current        : the beginning of the next UID block
     interval_UID               : the 4-number interval message buffer

   It is assumed that the calling function takes responsbility for
   validating error states before this function is called.

Input:


Output:


Returns:


Globals:


Notes:

   This function does not overwrite the status unless it itself 
   generates an error, leaving any current status error still valid.

*******************************************************************************/
void GetNextUID(cds_uint64 size)
{

   // First, check range validity of the size request
   if (size > max_block_size)
   {
      SetUIDInterval(0L, 0L, 0L, 0L);
      SetStatus(UID_CODE_BLOCK_TOO_LARGE, UID_OK);
      return;
   }

   // UID_current holds the starting value for the block we want to send
   if (UIDIsValid(size) == UID_FAILS)
   {
      SetUIDInterval(0L, 0L, 0L, 0L);
      SYS_UIDerrorMsg("UID Server: Server ran out of space\n");
      return;
   }

   SetUIDInterval(current_UID, size, 0L, 0L);
   current_UID += size;

   // Verified UID block within bounds. Will now check to see if position_file 
   // update is needed. If there is an error, must bail even though UID is
   // valid to ensure post-crash validity.
   if (current_UID >= current_position_UID) // need update?
      {
      if (IncrementUpdatePositionFile() != UID_OK)
	 {
         SetUIDInterval(0L, 0L, 0L, 0L);
         return;
         }
      }
   return;
}

/*******************************************************************************

Description: 

   This is a cleanup called between client connects - empty now
   but may be necessary in the future.

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 CleanConnectionPath(void)
{
    // INET version - automatically times out
    return UID_OK;
}

/*******************************************************************************

Description: 

   This function issues a UID_CODE_SERVER_KILL message to the server listening
   at the port specified of the local machine. If there is no server at that
   port, or a message cannot be established, an error is issued.

Input:


Output:


Returns:

   UID_OK if successful, otherwise UID_FAILS

Globals:


Notes:


*******************************************************************************/
cds_int32 IssueKillSignal(void)
{
   cds_int32                   server_connection_id;
   struct sockaddr_in      server_connection_data;
   cds_int32                   server_connection_data_size;
   struct hostent*         server_host_info;

   char writebuffer[12];

   // Create connection to server
   if ((server_connection_id = socket(AF_INET, SOCK_STREAM, 0)) < 0)
      return UID_FAILS;

   // Get server host info
   if ((server_host_info = gethostbyname("localhost")) == NULL)
      return UID_FAILS;

   // Fill in communications data
   bzero( (char*) &server_connection_data, sizeof(server_connection_data));
   server_connection_data.sin_family = AF_INET;
   server_connection_data.sin_addr.s_addr = 
      inet_addr( inet_ntoa( *((struct in_addr*) (server_host_info->h_addr_list[0]))));
   server_connection_data.sin_port = htons(tcp_port);
   server_connection_data_size = sizeof( server_connection_data );

   // Connect to server
   if (connect(server_connection_id, (struct sockaddr*) &server_connection_data,
	       server_connection_data_size) < 0)
      return UID_FAILS;
   
   // Package message
   SYS_UIDpackUIDRequestXdr(writebuffer, UID_CODE_SERVER_KILL, 0L);

   // Send message
   if (SYS_UIDwriten(server_connection_id, (char*)writebuffer, 12) == UID_FAILS)
      return UID_FAILS;

   // Close connection
   close( server_connection_id );

   return UID_OK;
}

/*******************************************************************************

Description: 

   Shuts down the server gracefully.

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void ShutdownServer(void)
{

   sprintf(SYS_UIDerr_str,"UID server on port %d: Received Shutdown instruction...\n",
	   tcp_port);
   SYS_UIDerrorMsg(SYS_UIDerr_str);

   // Close connections
   if (client_connection_id != 0)
      CloseConnection();
   if (connection_id != 0)
      close( connection_id );

   // exit without writing current_UID if a serious error has occurred
   if (unrecoverable_error_flag != UID_OK)
      exit(1);

   // if so far so good attempt to conserve UID space for next boot-up
   if (PositionWrite(current_UID) == UID_FAILS)
      exit(1);

   // if no errors
   exit(0);
}

/*******************************************************************************

Description: 

   Attempts to initialize the port

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
int AttemptPortInitialization(void)
{
      if (CleanConnectionPath() == UID_FAILS)
        return UID_FAILS;

      if (CreateConnection()    == UID_FAILS)
        return UID_FAILS;

      if (RegisterConnection()  == UID_FAILS)
        {
          CloseConnection();
          return UID_FAILS;
        }
      if (ActivateConnection()  == UID_FAILS)
        {
          CloseConnection();
          return UID_FAILS;
        }
      return UID_OK;
}

