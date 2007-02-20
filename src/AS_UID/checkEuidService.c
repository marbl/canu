
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
#include "SYS_UIDclient.h"
#include "SYS_UIDclient_local.h"
 
#define UID_CHECK_OK             0
#define UID_CHECK_ERROR          1

int main (int argc, char** argv) {
  char *option = NULL; 
  char *service = NULL; 
  cds_uint64 uid;
  int32 uidStatus;

  if(argc == 2) {
    option = argv[1];
    if(strcmp(option,"-h") == 0) {
       printf("USAGE: checkEuidService servicename\n");
       exit(1);
    }
    else {  
       service = option;
    }
  }
  else {
     printf("USAGE: checkEuidService servicename\n");
     exit(1);
  }
         
  printf("the service is %s\n", service);
  SYS_UIDset_euid_server(service);
  uidStatus = getGUIDBlock(1);
  
  if(uidStatus == UID_FAILS) 
  {
     printf("UID service %s is not available\n", service);
     return UID_CHECK_ERROR;
  }
  else
  {
     printf("UID service %s is available\n", service);
    return UID_CHECK_OK;
  } 
}
