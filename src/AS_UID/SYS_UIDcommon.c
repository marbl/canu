
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

#include "SYS_UIDcommon.h"
#include "SYS_UIDerror.h"

char   SYS_UIDerr_str[UID_ERR_STR_SIZE];
char   SYS_UIDmessage_array[UID_MESSAGE_SIZE];
char   SYS_UIDtype = UID_NO_TYPE;
char   SYS_UIDdebug_flag = 0;



//  Drop in replacement for missing xdr_u_hyper on FreeBSD4 and OS-X.
//  Set this in c_make.as if needed.
//
#ifdef NEEDXDRUHYPER

static
int
xdr_u_hyper(XDR *xdrs, uint64 *hp) {
  int res = 0;

  //  By definition of XDR, a 'long' is always 4 bytes, regardless of
  //  what the real long type is (according to OSF1's man page).

  unsigned long a, b;
  a = (*hp      ) & 0xffffffff;
  b = (*hp >> 32) & 0xffffffff;
  res += xdr_u_long(xdrs, &a);
  res += xdr_u_long(xdrs, &b);
  *hp   = b;
  *hp <<= 32;
  *hp  |= a;

  return(res);
}

#endif






/*******************************************************************************

Description: utility for flushing error string to stderr

*******************************************************************************/
void SYS_UIDperr(const char* message)
{
   fprintf(stderr, message);
   fflush(stderr);
}

/*******************************************************************************

Description: reads a specified number of bytes from a buffered input stream

*******************************************************************************/
int32  SYS_UIDreadn(int32 fd, char* ptr, int32 nbytes)
{
   int32 nleft;
   int32 nread;

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

*******************************************************************************/
int32  SYS_UIDwriten(int32 fd, char* ptr, int32 nbytes)
{
   int32 nleft;
   int32 nwritten;

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

Description: unpacks the 8-byte uint64 blocksize request

*******************************************************************************/
int32 SYS_UIDunpackUIDRequestXdr(char* readbuffer, int32* status, uint64* request_size)
{
   XDR xdr_stream;  
   char request_temp_array[12];
   int32 i;

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

Description: packs the 8-byte uint64 blocksize request

*******************************************************************************/
int32 SYS_UIDpackUIDRequestXdr(char* writebuffer, int32 status, uint64 request_size)
{
   XDR xdr_stream;
   char request_temp_array[12];
   int32 i;

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
             unpacks the int32 status code

*******************************************************************************/
int32 SYS_UIDunpackUIDMessageXdr(uint64* uid_interval, int32* status)
{
   XDR xdr_stream;
   char uid_temp_array[UID_MESSAGE_SIZE];
   int32 i;

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
               packs the int32 status code

*******************************************************************************/
int32 SYS_UIDpackUIDMessageXdr(uint64* uid_interval, int32 status)
{
   XDR xdr_stream;
   char uid_temp_array[UID_MESSAGE_SIZE];
   int32 i;

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
