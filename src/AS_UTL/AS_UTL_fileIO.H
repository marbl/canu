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

#ifndef AS_UTL_FILEIO_H
#define AS_UTL_FILEIO_H

static const char *rcsid_AS_UTL_FILEIO_H = "$Id$";

#include "AS_global.H"

//  Provides a safe and reliable mechanism for reading / writing
//  binary data.
//
//  Split writes/reads into smaller pieces, check the result of each
//  piece.  Really needed by OSF1 (V5.1), useful on other platforms to
//  be a little more friendly (big writes are usually not
//  interruptable).

void    AS_UTL_safeWrite(FILE *file, const void *buffer, const char *desc, size_t size, size_t nobj);
size_t  AS_UTL_safeRead (FILE *file, void *buffer,       const char *desc, size_t size, size_t nobj);

int     AS_UTL_mkdir(const char *dirname);
int     AS_UTL_unlink(const char *filename);

int     AS_UTL_fileExists(const char *path, int directory, int readwrite);
off_t   AS_UTL_sizeOfFile(const char *path);

off_t   AS_UTL_ftell(FILE *stream);
void    AS_UTL_fseek(FILE *stream, off_t offset, int whence);

class compressedFileReader {
public:
  compressedFileReader(char const *filename);
  ~compressedFileReader();

  FILE *operator*(void)     {  return(_file);  };
  FILE *file(void)          {  return(_file);  };

  bool  isCompressed(void)  {  return(_pipe);  };

private:
  char   _name[FILENAME_MAX];
  FILE  *_file;
  bool   _pipe;
  bool   _stdi;
};

class compressedFileWriter {
public:
  compressedFileWriter(char const *filename, int32 level=1);
  ~compressedFileWriter();

  FILE *operator*(void)     {  return(_file);  };
  FILE *file(void)          {  return(_file);  };

  bool  isCompressed(void)  {  return(_pipe);  };

private:
  char   _name[FILENAME_MAX];
  FILE  *_file;
  bool   _pipe;
  bool   _stdi;
};

#endif  //  AS_UTL_FILEIO_H
