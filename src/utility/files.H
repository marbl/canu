
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
 *  This file is derived from:
 *
 *    src/utility/AS_UTL_fasta.H
 *    src/utility/AS_UTL_fileIO.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef FILES_H
#define FILES_H

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

void    AS_UTL_findBaseFileName(char *basename, const char *filename);

bool    AS_UTL_readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F);

void    AS_UTL_mkdir(const char *dirname);
void    AS_UTL_rmdir(const char *dirname);

void    AS_UTL_symlink(const char *pathToFile, const char *pathToLink);

void    AS_UTL_unlink(const char *prefix, char separator='.', char const *suffix=NULL);

void    AS_UTL_rename(const char *oldname, const char *newname);

//  Tests on existence of paths, files and directories.
bool    pathExists     (const char *path);                        //  true if file or dir 'path' exists
bool    fileExists     (const char *path, bool writable=false);   //  true if file        'path' exists (and is writable)
bool    directoryExists(const char *path);                        //  true if         dir 'path' exists

off_t   AS_UTL_sizeOfFile(const char *path);
off_t   AS_UTL_sizeOfFile(FILE *file);

off_t   AS_UTL_ftell(FILE *stream);
void    AS_UTL_fseek(FILE *stream, off_t offset, int whence);

//  Search CANU_INSTALL_PATH, MERYL_INSTALL_PATH and PATH for a data file.
char   *findSharedFile(char *relpath, char *filename);

//  Read a file-of-files into a vector
void    AS_UTL_loadFileList(char *fileName, vector<char *> &FILE);

FILE   *AS_UTL_openInputFile (char const *prefix, char separator='.', char const *suffix=NULL, bool doOpen=true);
FILE   *AS_UTL_openOutputFile(char const *prefix, char separator='.', char const *suffix=NULL, bool doOpen=true);

void    AS_UTL_closeFile(FILE *&F, const char *filename, char separator, char const *suffix, bool critical=true);
void    AS_UTL_closeFile(FILE *&F, const char *filename=NULL, bool critical=true);

void    AS_UTL_createEmptyFile(char const *prefix, char separator='.', char const *suffix=NULL);

template<typename OBJ>
void    AS_UTL_loadFile(char const *prefix, char separator, char const *suffix, OBJ *objects, uint64 numberToLoad) {
  FILE    *file   = AS_UTL_openInputFile(prefix, separator, suffix);
  uint64   length = AS_UTL_sizeOfFile(file);

  if (numberToLoad == 0)
    numberToLoad = length / sizeof(OBJ);

  if (length / sizeof(OBJ) < numberToLoad) {
    if (suffix)
      fprintf(stderr, "AS_UTL_loadFile()-- File '%s%c%s' contains " F_U64 " objects, but asked to load " F_U64 ".\n",
              prefix, separator, suffix, length / sizeof(OBJ), numberToLoad);
    else
      fprintf(stderr, "AS_UTL_loadFile()-- File '%s' contains " F_U64 " objects, but asked to load " F_U64 ".\n",
              prefix, length / sizeof(OBJ), numberToLoad);
    exit(1);
  }

  loadFromFile(objects, "loadFile", numberToLoad, file);

  AS_UTL_closeFile(file, prefix, separator, suffix);
}

template<typename OBJ>
void    AS_UTL_loadFile(char const *prefix, OBJ *objects, uint64 numberToLoad) {
  AS_UTL_loadFile(prefix, '.', NULL, objects, numberToLoad);
}

template<typename OBJ>
void    AS_UTL_saveFile(char const *prefix, char separator, char const *suffix, OBJ *objects, uint64 numberToWrite) {
  FILE    *file = AS_UTL_openOutputFile(prefix, separator, suffix);

  writeToFile(objects, "saveFile", numberToWrite, file);

  AS_UTL_closeFile(file, prefix, separator, suffix);
}

template<typename OBJ>
void    AS_UTL_saveFile(char const *prefix, OBJ *objects, uint64 numberToWrite) {
  AS_UTL_saveFile(prefix, '.', NULL, objects, numberToWrite);
}







//  Read or write an object from or to a disk file.  Two forms are provided:
//
//    loadFromFile(POINTER, description, NUMBER_OF_OBJECTS, file);
//    loadFileFile(OBJECT,  description,                    file);
//
//  The first will load some number of objects into an array, and expects a
//  pointer.  The second will load one object into an instance, and expects
//  an actual object.
//
//  The second form hsa two implementations, one for OBJECTS and one for
//  POINTERS.  The OBJECT variant loads into the object directly, while the
//  POINTER varient loads into the pointed to object.  If you try to write a
//  pointer itself, you'll get unexpected results (but you shouldn't be
//  writing pointers to files anyway!).
//
uint64
loadFromFile(void        *objects,
             const char  *description,
             uint64       objectSize,
             uint64       nObjects,
             FILE        *file,
             bool         exact=true);

template<typename OBJ>
uint64
loadFromFile(OBJ *array, const char *description, uint64 nObjects, FILE *file, bool exact=true) {
  return(loadFromFile(array, description, sizeof(OBJ), nObjects, file, exact));
}

template<typename OBJ>
uint64
loadFromFile(OBJ &object, const char *description, FILE *file, bool exact=true) {
  return(loadFromFile(&object, description, sizeof(OBJ), 1, file, exact));
}

template<typename OBJ>
uint64
loadFromFile(OBJ *object, const char *description, FILE *file, bool exact=true) {
  return(loadFromFile(object, description, sizeof(OBJ), 1, file, exact));
}




void
writeToFile(void        *objects,
            const char  *description,
            uint64       objectSize,
            uint64       nObjects,
            FILE        *file);

template<typename OBJ>
void
writeToFile(OBJ *array, const char *description, uint64 nObjects, FILE *file) {
  writeToFile(array, description, sizeof(OBJ), nObjects, file);
}

template<typename OBJ>
void
writeToFile(OBJ &object, const char *description, FILE *file) {
  writeToFile(&object, description, sizeof(OBJ), 1, file);
}

template<typename OBJ>
void
writeToFile(OBJ *object, const char *description, FILE *file) {
  writeToFile(object, description, sizeof(OBJ), 1, file);
}







//  Writes sequence as fasta, with at most 'bl' letters per line (unlimited if 0).
void
AS_UTL_writeFastA(FILE *f,
                  char *s, int sl, int bl,
                  char *h, ...);

//  Writes FastQ, with Sanger QVs.
void
AS_UTL_writeFastQ(FILE  *f,
                  char  *s, int sl,
                  char  *q, int ql,   //  As Sanger QV, already encoded
                  char  *h, ...);

void
AS_UTL_writeFastQ(FILE  *f,
                  char  *s, int sl,
                  uint8 *q, int ql,   //  As Sanger QV, from integer values
                  char  *h, ...);


#include "files-compressed.H"
#include "files-buffered.H"
#include "files-buffered-implementation.H"
#include "files-memoryMapped.H"


#endif  //  FILES_H
