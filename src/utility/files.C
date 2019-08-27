
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
 *    src/utility/AS_UTL_fasta.C
 *    src/utility/AS_UTL_fileIO.C
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

#include "files.H"
#include "strings.H"

#include <stdarg.h>

//  Report ALL attempts to seek somewhere.
#undef DEBUG_SEEK

//  Report ALL file open/close events.
#undef SHOW_FILE_OPEN_CLOSE



//  Return the basename of a path -- that is, strip off any and all extensions.
//  Anything after the first dot after the last slash is removed.
//
//  But if a directory, do nothing.

void
AS_UTL_findBaseFileName(char *basename, const char *filename) {

  strcpy(basename, filename);

  if (directoryExists(basename))
    return;

  char  *slash = strrchr(basename, '/');
  char  *dot   = strchr((slash == NULL) ? basename : slash, '.');

  if (dot)
    *dot = 0;
}




void
writeToFile(void        *objects,
            const char  *description,
            uint64       objectSize,
            uint64       nObjects,
            FILE        *file) {

  uint64  nWritten  = 0;
  uint64  blockSize = (uint64)32 * 1024 * 1024 / objectSize;

  //  We previously split the write into 32MB chunks.  Comments indicated this
  //  was needed for OSF1 (V5.1).  In testing, certainly FreeBSD 11 isn't happy
  //  writing 16 GB of data at once; it seems to truncate to 32-bit somewhere.

  while (nWritten < nObjects) {
    uint64  toWrite = min(blockSize, nObjects - nWritten);

    errno = 0;
    uint64 written = fwrite(((char *)objects) + nWritten * objectSize, objectSize, toWrite, file);
    nWritten += written;

    if (errno)
      fprintf(stderr, "writeToFile()-- After writing %lu out of %lu '%s' objects (%lu bytes each): %s\n",
              nWritten, nObjects, description, objectSize, strerror(errno)), exit(1);
  }

  if (nWritten != nObjects)
    fprintf(stderr, "writeToFile()-- After writing %lu out of %lu '%s' objects (%lu bytes each): Short write\n",
            nWritten, nObjects, description, objectSize), exit(1);
}



uint64
loadFromFile(void        *objects,
             const char  *description,
             uint64       objectSize,
             uint64       nObjects,
             FILE        *file,
             bool         exact) {

  uint64  nLoaded   = 0;
  uint64  blockSize = (uint64)32 * 1024 * 1024 / objectSize;

  //  Reading doesn't seem to have the same size problem that writing does, but
  //  we still read in 32 MB chunks.

  while (nLoaded < nObjects) {
    uint64  toLoad = min(blockSize, nObjects - nLoaded);

    errno = 0;
    uint64 loaded = fread(((char *)objects) + nLoaded * objectSize, objectSize, toLoad, file);
    nLoaded += loaded;

    //  If we've loaded all requested, return successfully.

    if (nLoaded == nObjects)
      return(nLoaded);

    //if (loaded == 0)   //  The original version also returned if nothing was read.  Why?
    //  return(nLoaded);

    //  If we've hit eof, return the partial load.  Some loads, like reading overlaps, read until
    //  EOF, others should fail if a short read occurs.

    if (feof(file)) {
      if (exact == true) {
        fprintf(stderr, "loadFromFile()-- After loading %lu out of %lu '%s' objects (%lu bytes each): End of file\n",
                nLoaded, nObjects, description, objectSize);
        exit(1);
      }

      return(nLoaded);
    }

    //  If we've been interrupted, try again.

    if (errno == EINTR)
      continue;

    //  But if we hit an error, fail.

    if (errno)
      fprintf(stderr, "loadFromFile()-- After loading %lu out of %lu '%s' objects (%lu bytes each): %s\n",
              nLoaded, nObjects, description, objectSize, strerror(errno)), exit(1);
  }

  return(nLoaded);
}



#if 0
//  Reads a line, allocating space as needed.  Alternate implementatioin, probably slower than the
//  getc() based one below.
bool
readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F) {

  if (F == NULL)
    return(false);

  if ((L == NULL) || (Lmax == 0))
    allocateArray(L, Lmax = 4, resizeArray_clearNew);

  L[Lmax-2] = 0;
  L[Lmax-1] = 0;

  fgets(L, Lmax, F);

  Llen = strlen(L);

  fprintf(stderr, "READ Llen %u\n", Llen);

  //  fgets() will always NUL-terminate the string.  If the seocnd to last
  //  character exists and is not a newline, we didn't read the whole string.

  while ((L[Lmax-2] != 0) && (L[Lmax-2] != '\n')) {
    uint32   growth = 4;

    assert(Llen == Lmax - 1);

    resizeArray(L, Llen, Lmax, Lmax + growth);  //  Grow the array.
    L[Lmax-2] = 0;
    L[Lmax-1] = 0;

    fgets(L + Llen, 1 + growth, F);             //  Read more bytes.

    Llen += strlen(L + Llen);                   //  How many more?

    fprintf(stderr, "READ Llen %u Lmax %u '%s'\n", Llen, Lmax, L);
  }

  //  Trim trailing whitespace.

  while ((Llen > 0) && (isspace(L[Llen-1])))
    L[--Llen] = 0;

  return(true);
}
#endif



//  Reads a line of text from a file.  Trims off trailing whitespace, including newlines.
bool
AS_UTL_readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F) {

  if (F == NULL)
    return(false);

  if ((L == NULL) || (Lmax == 0))
    allocateArray(L, Lmax = 1024, resizeArray_clearNew);

  Llen = 0;

  int32   ch     = getc(F);
  uint32  growth = 1024;

  if (feof(F))
    return(false);

  while ((feof(F) == false) && (ch != '\n')) {
    if (Llen + 1 >= Lmax)
      resizeArray(L, Llen, Lmax, Lmax + growth, resizeArray_copyData | resizeArray_clearNew);  //  Grow the array.

    L[Llen++] = ch;

    ch = getc(F);
  }

  //  Terminate.

  L[Llen] = 0;

  //  Trim trailing whitespace.

  while ((Llen > 0) && (isspace(L[Llen-1])))
    L[--Llen] = 0;

  return(true);
}



//  Ensure that directory 'dirname' exists.
void
AS_UTL_mkdir(const char *dirname) {
  struct stat  st;

  //  Stat the file.  Don't fail if the file doesn't exist though.

  errno = 0;
  if ((stat(dirname, &st) == -1) && (errno != ENOENT))
    fprintf(stderr, "AS_UTL_mkdir()--  Couldn't stat '%s': %s\n", dirname, strerror(errno)), exit(1);

  //  If the file exists, and isn't a directory, fail.

  if ((errno == 0) && (S_ISDIR(st.st_mode) == false))
    fprintf(stderr, "AS_UTL_mkdir()--  ERROR!  '%s' is a file, and not a directory.\n", dirname), exit(1);

  //  Otherwise, make a directory.  Ignore any errors about the directory existing already.

  errno = 0;
  mkdir(dirname, S_IRWXU | S_IRWXG | S_IRWXO);

  if ((errno > 0) && (errno != EEXIST))
    fprintf(stderr, "AS_UTL_mkdir()--  Couldn't create directory '%s': %s\n", dirname, strerror(errno)), exit(1);

}



//  Remove a directory, or do nothing if the file doesn't exist.
void
AS_UTL_rmdir(const char *dirname) {

  if (directoryExists(dirname) == false)
    return;

  errno = 0;
  rmdir(dirname);
  if (errno)
    fprintf(stderr, "AS_UTL_rmdir()--  Failed to remove directory '%s': %s\n", dirname, strerror(errno)), exit(1);
}



void
AS_UTL_symlink(const char *pathToFile, const char *pathToLink) {

  //  Fail horribly if the file doesn't exist.

  if (pathExists(pathToFile) == false)
    fprintf(stderr, "AS_UTL_symlink()-- Original file '%s' doesn't exist, won't make a link to nothing.\n",
            pathToFile), exit(1);

  //  Succeed silently if the link already exists.

  if (pathExists(pathToLink) == true)
    return;

  //  Nope?  Make the link.

  errno = 0;
  symlink(pathToFile, pathToLink);
  if (errno)
    fprintf(stderr, "AS_UTL_symlink()-- Failed to make link '%s' pointing to file '%s': %s\n",
            pathToLink, pathToFile, strerror(errno)), exit(1);
}



//  Remove a file, or do nothing if the file doesn't exist.
void
AS_UTL_unlink(const char *prefix, char separator, char const *suffix) {
  char   filename[FILENAME_MAX];

  if (suffix)
    snprintf(filename, FILENAME_MAX, "%s%c%s", prefix, separator, suffix);
  else
    strncpy(filename, prefix, FILENAME_MAX-1);

  if (fileExists(filename) == false)
    return;

  errno = 0;
  unlink(filename);
  if (errno)
    fprintf(stderr, "AS_UTL_unlink()--  Failed to remove file '%s': %s\n",
            filename, strerror(errno)), exit(1);
}



//  Rename a file, or do nothing if the file doesn't exist.
void
AS_UTL_rename(const char *oldname, const char *newname) {

  if (pathExists(oldname) == false)
    return;

  errno = 0;
  rename(oldname, newname);
  if (errno)
    fprintf(stderr, "AS_UTL_renane()--  Failed to rename file '%s' to '%s': %s\n",
            oldname, newname, strerror(errno)), exit(1);
}



bool
pathExists(const char *path) {
  struct stat  s;

  if (stat(path, &s) == -1)
    return(false);

  return(true);               //  Stat-able?  Something exists there!
}



bool
fileExists(const char *path,
           bool        writable) {
  struct stat  s;

  if (stat(path, &s) == -1)
    return(false);

  if (s.st_mode & S_IFDIR)    //  Is a directory, not a file.
    return(false);

  if (writable == false)      //  User doesn't care if it's writable or not.
    return(true);

  if (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH))   //  User cares, and is writable.
    return(true);

  return(false);
}



bool
directoryExists(const char *path) {
  struct stat  s;

  if (stat(path, &s) == -1)
    return(false);

  if ((s.st_mode & S_IFDIR) == 0)     //  Is a file, not a directory.
    return(false);

#if 0
  if ((s.st_mode & S_IFDIR) &&
      (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
      (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)) &&
      (s.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)))
    ;
#endif

  return(true);
}





off_t
AS_UTL_sizeOfFile(const char *path) {
  struct stat  s;

  errno = 0;
  if (stat(path, &s) == -1)
    fprintf(stderr, "Failed to stat() file '%s': %s\n", path, strerror(errno)), exit(1);

  return(s.st_size);
}



off_t
AS_UTL_sizeOfFile(FILE *file) {
  struct stat  s;
  off_t        size = 0;

  errno = 0;
  if (fstat(fileno(file), &s) == -1)
    fprintf(stderr, "Failed to stat() FILE*: %s\n", strerror(errno)), exit(1);

  return(s.st_size);
}



off_t
AS_UTL_ftell(FILE *stream) {

  errno = 0;
  off_t pos = ftello(stream);

  if ((errno == ESPIPE) || (errno == EBADF))   //  Not a seekable stream.
    return(((off_t)1) < 42);                   //  Return some goofy big number.

  if (errno)
    fprintf(stderr, "AS_UTL_ftell()--  Failed with %s.\n", strerror(errno)), exit(1);

  return(pos);
}



void
AS_UTL_fseek(FILE *stream, off_t offset, int whence) {
  off_t   beginpos = AS_UTL_ftell(stream);

  //  If the stream is already at the correct position, just return.
  //
  //  Unless we're on FreeBSD.  For unknown reasons, FreeBSD fails
  //  updating the seqStore with mate links.  It seems to misplace the
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
#if !defined __FreeBSD__ && !defined __osf__ && !defined __APPLE__
  if ((whence == SEEK_SET) && (beginpos == offset)) {
#ifdef DEBUG_SEEK
    //  This isn't terribly informative, and adds a lot of clutter.
    //fprintf(stderr, "AS_UTL_fseek()--  seek to " F_OFF_T " (whence=%d); already there\n", offset, whence);
#endif
    return;
  }
#endif  //  __FreeBSD__

  if (fseeko(stream, offset, whence) != 0) {
    fprintf(stderr, "AS_UTL_fseek()--  Failed with %s.\n", strerror(errno));
    assert(errno == 0);
  }

#ifdef DEBUG_SEEK
  fprintf(stderr, "AS_UTL_fseek()--  seek to " F_OFF_T " (requested " F_OFF_T ", whence=%d) from " F_OFF_T "\n",
          AS_UTL_ftell(stream), offset, whence, beginpos);
#endif

  if (whence == SEEK_SET)
    assert(AS_UTL_ftell(stream) == offset);
}



//  Searches for a file in the 'share/' directory.
//
//  Checks for $CANU_INSTALL_PATH/relpath/filename
//             $MERYL_INSTALL_PATH/relpath/filename
//             $PATH/../relpath/filename
//             ./filename
//  and returns the first one that exists.  If no file is found,
//  NULL is returned.
//
//  'relpath' should be something like 'share/sequence'.
//
char *
findSharedFile(char *relpath, char *filename) {
  static
  char     fp[FILENAME_MAX + 1] = {0};
  char    *env;

  //  Does the file exist as is?

  if (fileExists(filename))
    return(filename);

  //  Does the file exist in any Canu installation?

  env = getenv("CANU_INSTALL_PATH");
  if (env != NULL) {
    snprintf(fp, FILENAME_MAX, "%s/%s/%s", env, relpath, filename);

    if (fileExists(fp))
      return(fp);
  }

  //  Does the file exist in any Meryl installation?

  env = getenv("MERYL_INSTALL_PATH");
  if (env != NULL) {
    snprintf(fp, FILENAME_MAX, "%s/%s/%s", env, relpath, filename);

    if (fileExists(fp))
      return(fp);
  }

  //  Does the file exist in any component of the path?

  env = getenv("PATH");
  if (env != NULL) {
    while ((*env != ':') && (*env != 0)) {   //  Until we exhaust the PATH,
      int32  fpp = 0;

      while ((*env != ':') && (*env != 0))   //  Copy up to the first ':'
        fp[fpp++] = *env++;                  //  or the terminating NUL.

      if (*env == ':')                       //  Skip over the delimiter.
        env++;

      fp[fpp] = 0;

      strcat(fp, "/../");                    //  Append the relative path.
      strcat(fp, relpath);                   //  and file we're searching
      strcat(fp, "/");                       //  for.
      strcat(fp, filename);

      if (fileExists(fp))
        return(fp);
    }
  }

  //  Nope, not found.

  return(NULL);
}



void
AS_UTL_loadFileList(char *fileName, vector<char *> &fileList) {

  FILE *F = AS_UTL_openInputFile(fileName);

  char *line = new char [FILENAME_MAX];

  fgets(line, FILENAME_MAX, F);

  while (!feof(F)) {
    chomp(line);
    fileList.push_back(line);
    line = new char [FILENAME_MAX];
    fgets(line, FILENAME_MAX, F);
  }

  delete [] line;

  AS_UTL_closeFile(F);
}



FILE *
AS_UTL_openInputFile(char const *prefix,
                     char        separator,
                     char const *suffix,
                     bool        doOpen) {
  char   filename[FILENAME_MAX];

  if (prefix == NULL)
    return(NULL);

  if (doOpen == false)
    return(NULL);

  if (suffix)
    snprintf(filename, FILENAME_MAX, "%s%c%s", prefix, separator, suffix);
  else
    strncpy(filename, prefix, FILENAME_MAX-1);

#ifdef SHOW_FILE_OPEN_CLOSE
  fprintf(stderr, "AS_UTL_openInputFile()-- Opening '%s'.\n", filename);
#endif

  errno = 0;

  FILE *F = fopen(filename, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", filename, strerror(errno)), exit(1);

  return(F);
}



FILE *
AS_UTL_openOutputFile(char const *prefix,
                      char        separator,
                      char const *suffix,
                      bool        doOpen) {
  char   filename[FILENAME_MAX];

  if (prefix == NULL)
    return(NULL);

  if (doOpen == false)
    return(NULL);

  if (suffix)
    snprintf(filename, FILENAME_MAX, "%s%c%s", prefix, separator, suffix);
  else
    strncpy(filename, prefix, FILENAME_MAX-1);

#ifdef SHOW_FILE_OPEN_CLOSE
  fprintf(stderr, "AS_UTL_openOutputFile()-- Creating '%s'.\n", filename);
#endif

  //  Unlink the file before opening for writes.  This prevents race
  //  conditions when two processes open the same file: the first process
  //  will create a new file, but the second process will simply reset the
  //  file to the start.  Both processes seem to keep their own file pointer,
  //  and eof seems to be (incorrectly) the larger of the two.  In effect,
  //  the second process is simply overwriting the first process (unless the
  //  second process writes data first, then the first process overwrites).
  //
  //  Very confusing.
  //
  AS_UTL_unlink(filename);

  errno = 0;

  FILE *F = fopen(filename, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", filename, strerror(errno)), exit(1);

  return(F);
}



void
AS_UTL_closeFile(FILE *&F, const char *prefix, char separator, char const *suffix, bool critical) {

  if ((F == NULL) || (F == stdout) || (F == stderr))
    return;

#ifdef SHOW_FILE_OPEN_CLOSE
  if ((prefix) && (suffix))
    fprintf(stderr, "AS_UTL_closeFile()-- Closing '%s%c%s'.\n", prefix, separator, suffix);
  else if (prefix)
    fprintf(stderr, "AS_UTL_closeFile()-- Closing '%s'.\n", prefix);
  else
    fprintf(stderr, "AS_UTL_closeFile()-- Closing (anonymous file).\n");
#endif

  errno = 0;

  fclose(F);

  F = NULL;

  if ((critical == false) || (errno == 0))
    return;

  if ((prefix) && (suffix))
    fprintf(stderr, "Failed to close file '%s%c%s': %s\n", prefix, separator, suffix, strerror(errno));
  else if (prefix)
    fprintf(stderr, "Failed to close file '%s': %s\n", prefix, strerror(errno));
  else
    fprintf(stderr, "Failed to close file: %s\n", strerror(errno));

  exit(1);
}


void    AS_UTL_closeFile(FILE *&F, const char *filename, bool critical) {
  AS_UTL_closeFile(F, filename, '.', NULL, critical);
}



void    AS_UTL_createEmptyFile(char const *prefix, char separator, char const *suffix) {
  FILE *file = AS_UTL_openOutputFile(prefix, separator, suffix);

  AS_UTL_closeFile(file, prefix, separator, suffix);
}



void
AS_UTL_writeFastA(FILE  *f,
                  char  *s, int sl, int bl,
                  char  *h, ...) {
  va_list ap;
  char   *o  = new char [sl + sl / 60 + 2];
  int     si = 0;
  int     oi = 0;

  while (si < sl) {
    o[oi++] = s[si++];

    if (bl != 0 && (si % bl) == 0)
      o[oi++] = '\n';
  }
  if (o[oi-1] != '\n')
    o[oi++] = '\n';
  o[oi] = 0;

  va_start(ap, h);
  vfprintf(f, h, ap);
  va_end(ap);

  writeToFile(o, "AS_UTL_writeFastA::seq", oi, f);

  delete [] o;
}


void
AS_UTL_writeFastQ(FILE  *f,
                  char  *s, int sl,
                  char  *q, int ql,
                  char  *h, ...) {
  va_list ap;
  int     qi = 0;
  int     oi = 0;

  assert(sl == ql);

  va_start(ap, h);
  vfprintf(f, h, ap);
  va_end(ap);

  writeToFile(s, "AS_UTL_writeFastQ::seq", sl, f);
  fprintf(f, "\n");

  fprintf(f, "+\n");
  writeToFile(q, "AS_UTL_writeFastQ::qlt", ql, f);
  fprintf(f, "\n");
}



void
AS_UTL_writeFastQ(FILE  *f,
                  char  *s, int sl,
                  uint8 *q, int ql,
                  char  *h, ...) {
  va_list ap;
  char   *o  = new char [ql + 1];
  int     qi = 0;
  int     oi = 0;

  assert(sl == ql);

  while (qi < ql)              //  Reencode the QV
    o[oi++] = q[qi++] + '!';   //  to the Sanger spec.
  o[oi] = 0;

  va_start(ap, h);
  vfprintf(f, h, ap);
  va_end(ap);

  writeToFile(s, "AS_UTL_writeFastQ::seq", sl, f);
  fprintf(f, "\n");

  fprintf(f, "+\n");
  writeToFile(o, "AS_UTL_writeFastQ::qlt", ql, f);
  fprintf(f, "\n");

  delete [] o;
}



