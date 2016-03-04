
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2005-SEP-26 to 2013-SEP-27
 *      are Copyright 2005,2007-2009,2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_UTL_FILEIO_H
#define AS_UTL_FILEIO_H

#include "AS_global.H"

#include <vector>

using namespace std;


//  Provides a safe and reliable mechanism for reading / writing
//  binary data.
//
//  Split writes/reads into smaller pieces, check the result of each
//  piece.  Really needed by OSF1 (V5.1), useful on other platforms to
//  be a little more friendly (big writes are usually not
//  interruptable).

#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif

void    AS_UTL_safeWrite(FILE *file, const void *buffer, const char *desc, size_t size, size_t nobj);
size_t  AS_UTL_safeRead (FILE *file, void *buffer,       const char *desc, size_t size, size_t nobj);

int     AS_UTL_mkdir(const char *dirname);
int     AS_UTL_unlink(const char *filename);

int     AS_UTL_fileExists(const char *path, int directory=false, int readwrite=false);
off_t   AS_UTL_sizeOfFile(const char *path);

off_t   AS_UTL_ftell(FILE *stream);
void    AS_UTL_fseek(FILE *stream, off_t offset, int whence);

//  Read a file-of-files into a vector
void    AS_UTL_loadFileList(char *fileName, vector<char *> &fileList);

class compressedFileReader {
public:
  compressedFileReader(char const *filename);
  ~compressedFileReader();

  FILE *operator*(void)     {  return(_file);  };
  FILE *file(void)          {  return(_file);  };

  bool  isCompressed(void)  {  return(_pipe);  };

private:
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
  FILE  *_file;
  bool   _pipe;
  bool   _stdi;
};

#endif  //  AS_UTL_FILEIO_H
