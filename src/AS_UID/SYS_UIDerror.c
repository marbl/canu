
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

static const char *rcsid = "$Id: SYS_UIDerror.c,v 1.7 2008-10-08 22:03:00 brianwalenz Exp $";

#include <syslog.h>
#include "SYS_UIDcommon.h"
#include "SYS_UIDerror.h"

/*******************************************************************************

Description: error handling for "accept" socket step

*******************************************************************************/
void  SYS_UIDhandleAcceptError(int32 err_code)
{
  switch(err_code)
    {
    case EBADF:
      SYS_UIDerrorMsg("UID Server: server socket invalid\n");
      break;
    case ECONNABORTED:
      SYS_UIDerrorMsg("UID Server: client socket connection aborted\n");
      break;
    case EFAULT:
      SYS_UIDerrorMsg("UID Server: server socket address not accessible\n");
      break;
    case EINTR:
      SYS_UIDerrorMsg("UID Server: client socket accept() interrupted\n");
      break;
    case EINVAL:
      SYS_UIDerrorMsg("UID Server: server socket not accepting connections\n");
      break;
    case EMFILE:
      SYS_UIDerrorMsg("UID Server: too many open file descriptors\n");
      break;
    case ENFILE:
      SYS_UIDerrorMsg("UID Server: max file descriptors already open\n");
      break;
    case ENOMEM:
      SYS_UIDerrorMsg("UID Server: ran out of memory for process descriptor\n");
      break;
    case ENOTSOCK:
      SYS_UIDerrorMsg("UID Server: strange error: server socket marked as file\n");
      break;
    case EOPNOTSUPP:
      SYS_UIDerrorMsg("UID Server: server socket cannot accept any connections\n");
      break;
    case EWOULDBLOCK:
      SYS_UIDerrorMsg("UID Server: server socket is non-blocking\n");
      break;
    default:
      SYS_UIDerrorMsg("UID Server: unknown error accepting socket\n");
      break;
    }
}

/*******************************************************************************

Description: error handling for "register" socket step

*******************************************************************************/
void  SYS_UIDhandleRegisterError(int32 err_code)
{
  switch(err_code)
    {
    case EACCES:
      SYS_UIDerrorMsg("UID Server: address is protected - user does not have permission\n");
      break;
    case EADDRINUSE:
      SYS_UIDerrorMsg("UID Server: address is already in use\n");
      break;
    case EADDRNOTAVAIL:
      SYS_UIDerrorMsg("UID Server: address is not available on local machine\n");
      break;
    case EAFNOSUPPORT:
      SYS_UIDerrorMsg("UID Server: the socket is protected or shut down\n");
      break;
    case EBADF:
      SYS_UIDerrorMsg("UID Server: the socket parameter is not valid\n");
      break;
    case EDESTADDRREQ:
      SYS_UIDerrorMsg("UID Server: the address is NULL\n");
      break;
    case EFAULT:
      SYS_UIDerrorMsg("UID Server: the address is out of user address space\n");
      break;
    case EINVAL:
      SYS_UIDerrorMsg("UID Server: the socket is already bound to an address\n");
      break;
    case EIO:
      SYS_UIDerrorMsg("UID Server: an IO error occurred\n");
      break;
    case EISDIR:
      SYS_UIDerrorMsg("UID Server: the address argument is a NULL pointer\n");
      break;
    case ELOOP:
      SYS_UIDerrorMsg("UID Server: too many symbolic links in address\n");
      break;
    case ENAMETOOLONG:
      SYS_UIDerrorMsg("UID Server: pathname too long\n");
      break;
    case ENOENT:
      SYS_UIDerrorMsg("UID Server: pathname is empty or does not exist\n");
      break;
    case ENOTDIR:
      SYS_UIDerrorMsg("UID Server: a component of path prefix is not a directory\n");
      break;
    case ENOTSOCK:
      SYS_UIDerrorMsg("UID Server: the socket parameter refers to a file\n");
      break;
    case EROFS:
      SYS_UIDerrorMsg("UID Server: the name is on a read-only filesystem\n");
      break;
    default:
      SYS_UIDerrorMsg("UID Server: unknown error binding socket\n");
      break;
    }
}

/*******************************************************************************

Description: error handling for socket "create" step

*******************************************************************************/
void  SYS_UIDhandleCreateError(int32 err_code)
{
   switch(err_code)
     {
     case EACCES:
       SYS_UIDerrorMsg("UID Server:handle_create_error: process does not have priviledges\n");
       break;
     case EAFNOSUPPORT:
       SYS_UIDerrorMsg("UID Server:handle_create_error: required addresses not available to kernel\n");
       break;
     case EMFILE:
       SYS_UIDerrorMsg("UID Server:handle_create_error: per-process descriptor table full\n");
       break;
     case ENFILE:
       SYS_UIDerrorMsg("UID Server:handle_create_error: no more file descriptors available\n");
       break;
     case ENOBUFS:
       SYS_UIDerrorMsg("UID Server:handle_create_error: insufficient resources\n");
       break;
     case ENOMEM:
       SYS_UIDerrorMsg("UID Server:handle_create_error: no memory to increase descriptor table\n");
       break;
     case EPERM:
       SYS_UIDerrorMsg("UID Server:handle_create_error: superuser necessary for raw socket\n");
       break;
     case EPROTONOSUPPORT:
       SYS_UIDerrorMsg("UID Server:handle_create_error: this address family not supported\n");
       break;
     case EPROTOTYPE:
       SYS_UIDerrorMsg("UID Server:handle_create_error: socket type not supported by protocol\n");
       break;
     default:
       SYS_UIDerrorMsg("UID Server:handle_create_error: unknown error creating socket\n");
       break;
     }
}

/*******************************************************************************

Description: error handling for socket "activation" step

*******************************************************************************/
void  SYS_UIDhandleActivateError(int32 err_code)
{
   switch(err_code)
     {
     case EBADF:
       SYS_UIDerrorMsg("UID Server:handle_create_error: the socket parameter is not valid\n");
       break;
     case EDESTADDRREQ:
       SYS_UIDerrorMsg("UID Server:handle_create_error: socket not bound to local address\n");
       break;
     case EINVAL:
       SYS_UIDerrorMsg("UID Server:handle_create_error: socket is already connected\n");;
       break;
     case ENOTSOCK:
       SYS_UIDerrorMsg("UID Server:handle_create_error: the parameter refers to file, not socket\n");
       break;
     case EOPNOTSUPP:
       SYS_UIDerrorMsg("UID Server:handle_create_error: socket type does not support listening\n");
       break;
     default:
       SYS_UIDerrorMsg("UID Server:handle_activate_error: unknown error listening for client\n");
       break;
     }
}

/*******************************************************************************

Description: method for using "openlog" for writing daemon messages

*******************************************************************************/
void  SYS_UIDerrorMsg(const char* err_str)
{
  char buffer[LINE_MAX];

  /* protect overflow of syslog linebuffer */
  if (strlen(err_str) > (LINE_MAX - 10)) {
    strncpy(buffer, err_str, (LINE_MAX - 10));
    buffer[LINE_MAX - 10] = '\0';
  } else {
    strcpy(buffer, err_str);
  }

   switch(SYS_UIDtype)
   {
   case UID_NO_TYPE:
      SYS_UIDperr(buffer);
      break;
   case UID_CLIENT_TYPE:
      SYS_UIDperr(buffer);
      break;
   case UID_SERVER_TYPE:
     openlog("uid_server", LOG_CONS | LOG_PID, LOG_DAEMON);
     syslog(LOG_ERR, buffer);
     closelog();
     break;
   default:
      SYS_UIDperr(buffer);
      break;
   }
}





