
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

//  A very simple, but not very good, uid "server".  Reads the next
//  available uid from a file, and then updates the file.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "AS_global.h"

void SYS_UIDset_euid_server(const char * servers)
{
	// do nothing for local server
}

void SYS_UIDset_euid_namespace(const char * namespaceName)
{
	// do nothing for local server
}

CDS_UID_t
getGUIDBlock(int guidRequestSize) {
  CDS_UID_t guidStart = 7180000;

  guidStart *= 1000;  //  Gets around integer overflow on 32-bit.
  guidStart *= 1000;

  errno = 0;
  FILE *F = fopen(".local-guid", "r+");
  if (errno == ENOENT) {
    //  Dang, doesn't exist!  Make it.
    F = fopen(".local-guid", "w+");
    fprintf(F, F_UID"\n", guidStart);
  } else if (errno) {
    fprintf(stderr, "getGUIDBlock()-- Can't open '.local-guid' for read/write, can't get last UID used: %s\n", strerror(errno));
    exit(1);
  } else {
    fscanf(F, F_UID"\n", &guidStart);
  }
  rewind(F);
  fprintf(F, F_UID"\n", guidStart + guidRequestSize);
  fclose(F);

  return(guidStart);
}
