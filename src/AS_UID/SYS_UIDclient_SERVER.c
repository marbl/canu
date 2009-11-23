
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

static const char *rcsid = "$Id: SYS_UIDclient_SERVER.c,v 1.1 2009-11-23 00:31:38 brianwalenz Exp $";

#include "uidserver_common.h"

void
SYS_UIDset_euid_server(const char * servers) {
}

void
SYS_UIDset_euid_namespace(const char * namespaceName) {
}

uint64
getGUIDBlock(int guidRequestSize) {
  UIDserverMessage        mesg     = {0};
  UIDserverMessage        resp     = {0};

  mesg.message = UIDserverMessage_REQUEST;
  mesg.bgnUID  = 0;
  mesg.endUID  = 0;
  mesg.numUID  = guidRequestSize;

  int32  server = connectToServer("blackdeath.home", "6969");
  sendMessage(server, &mesg);
  recvMessage(server, &resp);
  close(server);

  //fprintf(stderr, "Requested %d, got %d-%d\n", mesg.numUID, resp.bgnUID, resp.endUID);

  return(resp.bgnUID);
}
