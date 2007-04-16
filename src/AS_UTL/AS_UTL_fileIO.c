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

//static char CM_ID[] = "$Id: AS_UTL_fileIO.c,v 1.9 2007-04-16 15:29:09 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
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
AS_UTL_safeWrite(FILE *file, const void *buffer, char *desc, size_t size, size_t nobj) {
  size_t  position = 0;
  size_t  length   = 32 * 1024 * 1024 / size;
  size_t  towrite  = 0;
  size_t  written  = 0;
  size_t  nbytes   = size * nobj;

  while (position < nobj) {
    towrite = length;
    if (position + towrite > nobj)
      towrite = nobj - position;

    errno = 0;
    written = fwrite(((char *)buffer) + position * size, size, towrite, file);

    if (errno) {
      fprintf(stderr, "safeWrite()-- Write failure on %s: %s\n", desc, strerror(errno));
      fprintf(stderr, "safeWrite()-- Wanted to write "F_SIZE_T" objects (size="F_SIZE_T"), wrote "F_SIZE_T".\n",
              towrite, size, written);
      exit(1);
    }

    position += written;
  }
}


size_t
AS_UTL_safeRead(FILE *file, void *buffer, char *desc, size_t size, size_t nobj) {
  size_t  position = 0;
  size_t  length   = 32 * 1024 * 1024 / size;
  size_t  toread   = 0;
  size_t  written  = 0;  //  readen?

  while (position < nobj) {
    toread = length;
    if (position + toread > nobj)
      toread = nobj - position;

    errno = 0;
    written = fread(((char *)buffer) + position * size, size, toread, file);
    position += written;

    if (feof(file) || (written == 0))
      return(position);

    if ((errno) && (errno != EINTR)) {
      fprintf(stderr, "safeRead()-- Read failure on %s: %s.\n", desc, strerror(errno));
      fprintf(stderr, "safeRead()-- Wanted to read "F_SIZE_T" objects (size="F_SIZE_T"), read "F_SIZE_T".\n",
              toread, size, written);
      exit(1);
    }
  }

  return(position);
}


//  Ensure that directory 'dirname' exists.  Returns true if the
//  directory needed to be created, false if it already exists.
int
AS_UTL_mkdir(const char *dirname) {
  struct stat  st;

  errno = 0;
  stat(dirname, &st);
  if (errno == 0) {
    if (S_ISDIR(st.st_mode))
      return(0);

    fprintf(stderr, "AS_UTL_mkdir()--  ERROR!  '%s' is a file, and not a directory.\n", dirname);
    exit(1);
  }

  if (errno != ENOENT) {
    fprintf(stderr, "AS_UTL_mkdir()--  Couldn't stat '%s': %s\n", dirname, strerror(errno));
    exit(1);
  }

  errno = 0;
  mkdir(dirname, S_IRWXU | S_IRWXG | S_IRWXO);
  if (errno) {
    fprintf(stderr, "AS_UTL_mkdir()--  Couldn't create directory '%s': %s\n", dirname, strerror(errno));
    exit(1);
  }

  return(1);
}

