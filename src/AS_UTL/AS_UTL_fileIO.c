/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005, J. Craig Venter Institute. All rights reserved.
 * Author: Brian Walenz
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

//static char CM_ID[] = "$Id: AS_UTL_fileIO.c,v 1.4 2007-02-18 14:04:50 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#include "AS_global.h"

//  Provides a safe and reliable mechanism for reading / writing
//  binary data.
//
//  Split writes/reads into smaller pieces, check the result of each
//  piece.  Really needed by OSF1 (V5.1), useful on other platforms to
//  be a little more friendly (big writes are usually not
//  interruptable).

void
AS_UTL_safeWrite(FILE *file, const void *buffer, char *desc, size_t nbytes) {
  size_t  position = 0;
  size_t  length   = 32 * 1024 * 1024;
  size_t  towrite  = 0;
  size_t  written  = 0;
  int     filedes  = fileno(file);

  while (position < nbytes) {
    towrite = length;
    if (position + towrite > nbytes)
      towrite = nbytes - position;

    errno = 0;
    written = fwrite(((char *)buffer) + position, sizeof(char), towrite, file);

    if (errno) {
      fprintf(stderr, "safeWrite()-- Write failure on %s: %s\n", desc, strerror(errno));
      fprintf(stderr, "safeWrite()-- Wanted to write "F_SIZE_T" bytes, wrote "F_SIZE_T".\n", towrite, written);
      exit(1);
    }

    position += written;
  }
}


int
AS_UTL_safeRead(FILE *file, void *buffer, char *desc, size_t nbytes) {
  size_t  position = 0;
  size_t  length   = 32 * 1024 * 1024;
  size_t  toread   = 0;
  size_t  written  = 0;  //  readen?
  int     filedes  = fileno(file);

  while (position < nbytes) {
    toread = length;
    if (position + toread > nbytes)
      toread = nbytes - position;

    errno = 0;
    written = fread(((char *)buffer) + position, sizeof(char), toread, file);

    if (feof(file) || (written == 0))
      return(TRUE);

    if ((errno) && (errno != EINTR)) {
      fprintf(stderr, "safeRead()-- Read failure on %s: %s.\n", desc, strerror(errno));
      fprintf(stderr, "safeRead()-- Wanted to read "F_SIZE_T" bytes, read "F_SIZE_T".\n", toread, written);
      exit(1);
    }

    position += written;
  }

  return(FALSE);
}

