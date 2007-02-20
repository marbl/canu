
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

#include "euidH.h"
#include "soapEUIDServerService.nsmap"
#include <assert.h>


#include "SYS_UIDcommon.h"
#include "SYS_UIDerror.h"

const char *TIGR_DefaultEuidServerNames = "http://tools2.tigr.org/servlet/axis/services/EUIDServer";
char * EuidServerNames = NULL;
const int DUMMY_BUFFER_LENGTH = 1000;


void SYS_UIDset_euid_server(const char * servers)
{
  EuidServerNames=strdup(servers);
  assert(EuidServerNames != NULL);
}


CDS_UID_t getGUIDBlock(int guidRequestSize)
{
  CDS_UID_t guidStart = 0;


  struct soap soap;
  struct impl__getEUIDBlockResponse euid;
  char dummy[DUMMY_BUFFER_LENGTH+1];
  char euidServerBuffer[2048];
  char *servers[10];
  int loop;
  int i,j;

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
    //// sprintf(dummy,"http://%s/axis/services/EUIDServer", servers[i]);

    // Jason Miller October 2005.
    // We now expect to be given the full URL.

    // Avoid write beyond buffer on unexpectedly long server names.
    // If URL is too long, truncate it. Let the SOAP request fail.
    strncpy (dummy, servers[i], DUMMY_BUFFER_LENGTH);

    fprintf(stderr, "Trying to contact %s, block size %ld\n", 
	    dummy, guidRequestSize);

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

  return guidStart;
}
