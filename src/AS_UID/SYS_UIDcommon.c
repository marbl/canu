
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Attic/SYS_UIDcommon.c,v $
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

Revision 1.1.1.1  2003/05/08 18:40:12  aaronhalpern
versions from TIGR

Revision 1.2  2001/09/25 23:03:20  mpop
Dos2Unixed

Revision 1.1.1.1  2001/09/25 20:21:05  mpop
Celera Assembler

Revision 1.5  1999/10/13 19:02:57  sdmurphy
misc small changes for debugging

Revision 1.4  1999/07/14 17:24:33  stine
update_cds script was executed against these files.

Revision 1.3  1999/01/13 15:49:57  sdmurphy
added error.h

Revision 1.2  1999/01/13 14:29:16  sdmurphy
version 0 prelim

Revision 1.1  1998/12/30 19:56:07  sdmurphy
Renamed uid_common.c to SYS_UIDcommon.c

Revision 1.2  1998/12/22 21:20:11  sdmurphy
minor reorganization

Revision 1.1  1998/12/17 21:23:00  sdmurphy
common .c file for UID libs and progs

**********************************************************************/

/**********************************************************************
Module:

Description:

Assumptions:

**********************************************************************/

#include "SYS_UIDcommon.h"
#include "SYS_UIDerror.h"

char   SYS_UIDerr_str[UID_ERR_STR_SIZE];
char   SYS_UIDmessage_array[UID_MESSAGE_SIZE];
char   SYS_UIDtype = UID_NO_TYPE;
char   SYS_UIDdebug_flag = 0;

/*******************************************************************************

Description: utility for flushing error string to stderr


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
void SYS_UIDperr(const char* message)
{
   fprintf(stderr, message);
   fflush(stderr);
}

/*******************************************************************************

Description: reads a specified number of bytes from a buffered input stream


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  SYS_UIDreadn(cds_int32 fd, char* ptr, cds_int32 nbytes)
{
   cds_int32 nleft;
   cds_int32 nread;

   nleft = nbytes;
   while (nleft > 0) 
   {
      nread = recv(fd, ptr, nleft, 0);
      if (nread < 0)
         return UID_FAILS;
      else if (nread == 0)
         break;
      nleft -= nread;
      ptr += nread;
   }
   return UID_OK;    
}

/*******************************************************************************

Description: writes a specified number of bytes to a buffered output stream


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32  SYS_UIDwriten(cds_int32 fd, char* ptr, cds_int32 nbytes)
{
   cds_int32 nleft;
   cds_int32 nwritten;

   nleft = nbytes;
   while (nleft > 0) 
   {
      nwritten = send(fd, ptr, nleft, 0);
      if (nwritten <= 0)
         return UID_FAILS;
      nleft -= nwritten;
      ptr += nwritten;
   }
   return UID_OK;    
}

/*******************************************************************************

Description: unpacks the 8-byte cds_uint64 blocksize request


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDunpackUIDRequestXdr(char* readbuffer, cds_int32* status, cds_uint64* request_size)
{
   XDR xdr_stream;  
   char request_temp_array[12];
   cds_int32 i;

   memcpy(request_temp_array, readbuffer, 12);
   xdrmem_create(&xdr_stream, request_temp_array, 12, XDR_DECODE);
   if (xdr_int(&xdr_stream, status)==0)
   {
      SYS_UIDerrorMsg("unpack_UID_message_xdr: xdr_int() problem");
      return UID_FAILS;
   }
   if (xdr_u_hyper(&xdr_stream, request_size)==0)
   {
      SYS_UIDerrorMsg("unpack_UID_request_xdr: xdr_u_hyper() problem");
      return UID_FAILS;
   }
   xdr_destroy(&xdr_stream);
   return UID_OK;
}

/*******************************************************************************

Description: packs the 8-byte cds_uint64 blocksize request


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDpackUIDRequestXdr(char* writebuffer, cds_int32 status, cds_uint64 request_size)
{
   XDR xdr_stream;
   char request_temp_array[12];
   cds_int32 i;

   xdrmem_create(&xdr_stream, request_temp_array, 12, XDR_ENCODE);
   if (xdr_int(&xdr_stream, &status)==0)
   {
      SYS_UIDerrorMsg("unpack_UID_message_xdr: xdr_int() problem");
      return UID_FAILS;
   }
   if (xdr_u_hyper(&xdr_stream, &request_size)==0)
   {
      SYS_UIDerrorMsg("pack_UID_request_xdr: xdr_u_hyper() problem");
      return UID_FAILS;
   }
   memcpy(writebuffer, request_temp_array, 12);

   xdr_destroy(&xdr_stream);
   return UID_OK;
}

/*******************************************************************************

Description: unpacks the uid_interval array of 8-byte uint64s and also
             unpacks the cds_int32 status code

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDunpackUIDMessageXdr(cds_uint64* uid_interval, cds_int32* status)
{
   XDR xdr_stream;
   char uid_temp_array[UID_MESSAGE_SIZE];
   cds_int32 i;

   memcpy(uid_temp_array, SYS_UIDmessage_array, UID_MESSAGE_SIZE);
   xdrmem_create(&xdr_stream, uid_temp_array, UID_MESSAGE_SIZE, XDR_DECODE);
   if (xdr_u_hyper(&xdr_stream, uid_interval)==0)
   {
      SYS_UIDerrorMsg("unpack_UID_message_xdr: xdr_u_hyper() problem");
      return UID_FAILS;
   }
   if (xdr_u_hyper(&xdr_stream, uid_interval+1)==0)
   {
      SYS_UIDerrorMsg("unpack_UID_message_xdr: xdr_u_hyper() problem");
      return UID_FAILS;
   }
   if (xdr_u_hyper(&xdr_stream, uid_interval+2)==0)
   {
      SYS_UIDerrorMsg("unpack_UID_message_xdr: xdr_u_hyper() problem");
      return UID_FAILS;
   }
   if (xdr_u_hyper(&xdr_stream, uid_interval+3)==0)
   {
      SYS_UIDerrorMsg("unpack_UID_message_xdr: xdr_u_hyper() problem");
      return UID_FAILS;
   }
   if (xdr_int(&xdr_stream, status)==0)
   {
      SYS_UIDerrorMsg("unpack_UID_message_xdr: xdr_int() problem");
      return UID_FAILS;
   }
   xdr_destroy(&xdr_stream);
   return UID_OK;
}
   
/*******************************************************************************

Description:   packs the uid_interval array of 8-byte uint64s and also
               packs the cds_int32 status code

Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/
cds_int32 SYS_UIDpackUIDMessageXdr(cds_uint64* uid_interval, cds_int32 status)
{
   XDR xdr_stream;
   char uid_temp_array[UID_MESSAGE_SIZE];
   cds_int32 i;

   xdrmem_create(&xdr_stream, uid_temp_array, UID_MESSAGE_SIZE, XDR_ENCODE);
   if (xdr_u_hyper(&xdr_stream, &(uid_interval[0]))==0)
   {
      SYS_UIDerrorMsg("pack_UID_message_xdr: xdr_u_hyper() failed");
      return UID_FAILS;
   }
   if (xdr_u_hyper(&xdr_stream, &(uid_interval[1]))==0)
   {
      SYS_UIDerrorMsg("pack_UID_message_xdr: xdr_u_hyper() failed");
      return UID_FAILS;
   }
   if (xdr_u_hyper(&xdr_stream, &(uid_interval[2]))==0)
   {
      SYS_UIDerrorMsg("pack_UID_message_xdr: xdr_u_hyper() failed");
      return UID_FAILS;
   }
   if (xdr_u_hyper(&xdr_stream, &(uid_interval[3]))==0)
   {
      SYS_UIDerrorMsg("pack_UID_message_xdr: xdr_u_hyper() failed");
      return UID_FAILS;
   }
   if (xdr_int(&xdr_stream, &status)==0)
   {
      SYS_UIDerrorMsg("pack_UID_message_xdr: xdr_int() failed");
      return UID_FAILS;
   }
   memcpy(SYS_UIDmessage_array, uid_temp_array, UID_MESSAGE_SIZE);
   xdr_destroy(&xdr_stream);

   return UID_OK;
}
   

/*******************************************************************************

Description:   a log message facility for debugging certain timeout problems.
               writes to file /tmp/uidlog


Input:


Output:


Returns:


Globals:


Notes:


*******************************************************************************/

void   SYS_UIDlogMessage(const char* message)
{
  FILE* lfp;

  /* generate timestamp */
  system("date >> /tmp/uidlog");
  /* open file to append */
  if ((lfp = fopen("/tmp/uidlog", "a")) == NULL) {
    /* do nothing */
    return;
  } 
  fprintf(lfp, " : %s\n", message);
  fclose(lfp);
  return;
}




