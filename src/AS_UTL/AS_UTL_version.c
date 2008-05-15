
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

static const char CM_ID[] = "$Id: AS_UTL_version.c,v 1.14 2008-05-15 00:34:56 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_version.h"
#include "AS_UTL_fileIO.h"
#include "AS_MSG_pmesg.h"


int
VersionStampADT(AuditMesg *adt_mesg, int argc, char **argv) {
  time_t      t        = time(NULL);
  int         i        = 0;
  char       *identC   = (char *)safe_calloc(sizeof(char), FILENAME_MAX);
  char       *ident    = (char *)safe_calloc(sizeof(char), 1048576);
  char       *identP   = ident;
  AuditLine  *adt_line = (AuditLine *)safe_malloc(sizeof(AuditLine));

  strcpy(identP, "command:");
  crunch(identP);

  for (i=0; i<argc; i++) {
    *identP = ' ';
    identP++;
    *identP = 0;
    strcpy(identP, argv[i]);
    crunch(identP);
  }

  *identP = '\n';
  identP++;
  *identP = 0;

  sprintf(identP, "started: %s", ctime(&t));
  crunch(identP);
  sprintf(identP, "directory: %s\n", getcwd(identC, FILENAME_MAX));
  crunch(identP);

  AppendAuditLine_AS(adt_mesg, adt_line, t, argv[0], "(no version)", ident);

  safe_free(identC);
  safe_free(identP);
} 
