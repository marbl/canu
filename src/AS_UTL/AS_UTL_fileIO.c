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

//static char CM_ID[] = "$Id: AS_UTL_fileIO.c,v 1.15 2007-09-02 22:56:56 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"

//  Report ALL attempts to seek somewhere.
#undef DEBUG_SEEK

//  Use ftell() to verify that we wrote the expected number of bytes,
//  and that we ended up at the expected location.
#define VERIFY_WRITE_POSITIONS

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

#ifdef VERIFY_WRITE_POSITIONS
  off_t   expectedposition = AS_UTL_ftell(file) + nobj * size;
  if (errno)
    //  If we return, and errno is set, the stream isn't seekable.
    expectedposition = 0;
#endif

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
      assert(errno == 0);
    }

    position += written;
  }

  //  This catches a bizarre bug on FreeBSD (6.1 for sure, 4.10 too, I
  //  think) where we write at the wrong location; see fseek below.
  //
  //  UNFORTUNATELY, you can't ftell() on stdio.
  //
#ifdef VERIFY_WRITE_POSITIONS
  if ((expectedposition > 0) &&
      (AS_UTL_ftell(file) != expectedposition)) {
    fprintf(stderr, "safeWrite()-- EXPECTED "F_OFF_T", ended up at "F_OFF_T"\n",
            expectedposition, AS_UTL_ftell(file));
    assert(AS_UTL_ftell(file) == expectedposition);
  }
#endif
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
      goto finish;

    if ((errno) && (errno != EINTR)) {
      fprintf(stderr, "safeRead()-- Read failure on %s: %s.\n", desc, strerror(errno));
      fprintf(stderr, "safeRead()-- Wanted to read "F_SIZE_T" objects (size="F_SIZE_T"), read "F_SIZE_T".\n",
              toread, size, written);
      assert(errno == 0);
    }
  }

 finish:
  if (position != nobj)
    fprintf(stderr, "AS_UTL_safeRead()--  Short read; wanted "F_SIZE_T" objects, read "F_SIZE_T" instead.\n",
            nobj, position);
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




off_t
AS_UTL_ftell(FILE *stream) {
  off_t  pos = 0;
  errno = 0;
  pos = ftello(stream);
  if (errno == ESPIPE)
    //  Not a seekable stream.  Return some goofy big number.
    return(((off_t)1) < 42);
  if (errno) {
    fprintf(stderr, "AS_UTL_ftell()--  Failed with %s.\n", strerror(errno));
    assert(errno == 0);
  }
  return(pos);
}



void
AS_UTL_fseek(FILE *stream, off_t offset, int whence) {
  off_t   beginpos = AS_UTL_ftell(stream);

  //  If the stream is already at the correct position, just return.
  //
  //  Unless we're on FreeBSD.  For unknown reasons, FreeBSD fails
  //  updating the gkpStore with mate links.  It seems to misplace the
  //  file pointer, and ends up writing the record to the wrong
  //  location.  ftell() is returning the correct current location,
  //  and so AS_PER_genericStore doesn't seek() and just writes to the
  //  current position.  At the end of the write, we're off by 4096
  //  bytes.
  //
  //  LINK 498318175,1538 <-> 498318174,1537
  //  AS_UTL_fseek()--  seek to 159904 (whence=0); already there
  //  safeWrite()-- write nobj=1x104 = 104 bytes at position 159904
  //  safeWrite()-- wrote nobj=1x104 = 104 bytes position now 164000
  //  safeWrite()-- EXPECTED 160008, ended up at 164000
  //
#if !defined __FreeBSD__ && !defined __osf__
  if ((whence == SEEK_SET) && (beginpos == offset)) {
#ifdef DEBUG_SEEK
    //  This isn't terribly informative, and adds a lot of clutter.
    //fprintf(stderr, "AS_UTL_fseek()--  seek to "F_OFF_T" (whence=%d); already there\n", offset, whence);
#endif
    return;
  }
#endif  //  __FreeBSD__

  errno = 0;
  fseeko(stream, offset, whence);
  if (errno) {
    fprintf(stderr, "AS_UTL_fseek()--  Failed with %s.\n", strerror(errno));
    assert(errno == 0);
  }

#ifdef DEBUG_SEEK
  fprintf(stderr, "AS_UTL_fseek()--  seek to "F_OFF_T" (requested "F_OFF_T", whence=%d) from "F_OFF_T"\n",
          AS_UTL_ftell(stream), offset, whence, beginpos);
#endif

  if (whence == SEEK_SET)
    assert(AS_UTL_ftell(stream) == offset);
}

