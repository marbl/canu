
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "files.H"
#include "strings.H"

#include <stdarg.h>

//  Report ALL attempts to seek somewhere.
#undef DEBUG_SEEK

//  Use ftell() to verify that we wrote the expected number of bytes,
//  and that we ended up at the expected location.
#undef VERIFY_WRITE_POSITIONS




//  Return the basename of a path -- that is, strip off any and all extensions.
//  Anything after the first dot after the last slash is removed.
//
//  But if a directory, do nothing.

void
AS_UTL_findBaseFileName(char *basename, const char *filename) {

  strcpy(basename, filename);

  if (AS_UTL_fileExists(basename, true, false))
    return;

  char  *slash = strrchr(basename, '/');
  char  *dot   = strchr((slash == NULL) ? basename : slash, '.');

  if (dot)
    *dot = 0;
}




//  Provides a safe and reliable mechanism for reading / writing
//  binary data.
//
//  Split writes/reads into smaller pieces, check the result of each
//  piece.  Really needed by OSF1 (V5.1), useful on other platforms to
//  be a little more friendly (big writes are usually not
//  interruptable).

void
AS_UTL_safeWrite(FILE *file, const void *buffer, const char *desc, size_t size, size_t nobj) {
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
      fprintf(stderr, "safeWrite()-- Wanted to write " F_SIZE_T " objects (size=" F_SIZE_T "), wrote " F_SIZE_T ".\n",
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
    fprintf(stderr, "safeWrite()-- EXPECTED " F_OFF_T ", ended up at " F_OFF_T "\n",
            expectedposition, AS_UTL_ftell(file));
    assert(AS_UTL_ftell(file) == expectedposition);
  }
#endif
}


size_t
AS_UTL_safeRead(FILE *file, void *buffer, const char *desc, size_t size, size_t nobj) {
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
      fprintf(stderr, "safeRead()-- Wanted to read " F_SIZE_T " objects (size=" F_SIZE_T "), read " F_SIZE_T ".\n",
              toread, size, written);
      assert(errno == 0);
    }
  }

 finish:
  //  Just annoys developers.  Stop it.
  //if (position != nobj)
  //  fprintf(stderr, "AS_UTL_safeRead()--  Short read; wanted " F_SIZE_T " objects, read " F_SIZE_T " instead.\n",
  //          nobj, position);
  return(position);
}



#if 0
//  Reads a line, allocating space as needed.  Alternate implementatioin, probably slower than the
//  getc() based one below.
bool
readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F) {

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
  stat(dirname, &st);

  if ((errno > 0) && (errno != ENOENT))
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

  if (AS_UTL_fileExists(dirname, false, false) == false)
    return;

  errno = 0;
  rmdir(dirname);
  if (errno)
    fprintf(stderr, "AS_UTL_rmdir()--  Failed to remove directory '%s': %s\n", dirname, strerror(errno)), exit(1);
}



void
AS_UTL_symlink(const char *pathToFile, const char *pathToLink) {

  //  Fail horribly if the file doesn't exist.

  if (AS_UTL_fileExists(pathToFile, false, false) == false)
    fprintf(stderr, "AS_UTL_symlink()-- Original file '%s' doesn't exist, won't make a link to nothing.\n",
            pathToFile), exit(1);

  //  Succeed silently if the link already exists.

  if (AS_UTL_fileExists(pathToLink, false, false) == true)
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
AS_UTL_unlink(const char *filename) {

  if (AS_UTL_fileExists(filename, false, false) == false)
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

  if (AS_UTL_fileExists(oldname, false, false) == false)
    return;

  errno = 0;
  rename(oldname, newname);
  if (errno)
    fprintf(stderr, "AS_UTL_renane()--  Failed to rename file '%s' to '%s': %s\n",
            oldname, newname, strerror(errno)), exit(1);
}




//  Returns true if the named file/directory exists, and permissions
//  allow us to read and/or write.
//
int
AS_UTL_fileExists(const char *path,
                  int directory,
                  int readwrite) {
  struct stat  s;

  errno = 0;
  stat(path, &s);
  if (errno)
    return(0);

  if ((directory == 1) &&
      (readwrite == 0) &&
      (s.st_mode & S_IFDIR) &&
      (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
      (s.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)))
    return(1);

  if ((directory == 1) &&
      (readwrite == 1) &&
      (s.st_mode & S_IFDIR) &&
      (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
      (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)) &&
      (s.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)))
    return(1);

  if ((directory == 0) &&
      (readwrite == 0) &&
      (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)))
    return(1);

  if ((directory == 0) &&
      (readwrite == 1) &&
      (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
      (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)))
    return(1);

  return(0);
}





off_t
AS_UTL_sizeOfFile(const char *path) {
  struct stat  s;
  off_t        size = 0;

  errno = 0;
  stat(path, &s);
  if (errno) {
    fprintf(stderr, "Failed to stat() file '%s': %s\n", path, strerror(errno));
    exit(1);
  }

  //  gzipped files contain a file contents list, which we can
  //  use to get the uncompressed size.
  //
  //  gzip -l <file>
  //  compressed        uncompressed  ratio uncompressed_name
  //       14444               71680  79.9% up.tar
  //
  //  bzipped files have no contents and we just guess.

  if        (strcasecmp(path+strlen(path)-3, ".gz") == 0) {
    char   cmd[FILENAME_MAX], *p = cmd;

    snprintf(cmd, FILENAME_MAX, "gzip -l '%s'", path);

    FILE *F = popen(cmd, "r");
    fgets(cmd, FILENAME_MAX, F);    //   compressed uncompressed  ratio uncompressed_name
    fgets(cmd, FILENAME_MAX, F);    //     30264891     43640320  30.6% file
    pclose(F);

    while (isspace(*p) == true)  p++;  //  Skip spaces at the start of the line
    while (isspace(*p) == false) p++;  //  Skip the compressed size
    while (isspace(*p) == true)  p++;  //  Skip spaces

    size = strtoull(p, NULL, 10);      //  Retain the uncompresssed size
  }

  else if (strcasecmp(path+strlen(path)-4, ".bz2") == 0) {
    size = s.st_size * 14 / 10;
  }

  else {
    size = s.st_size;
  }

  return(size);
}



off_t
AS_UTL_sizeOfFile(FILE *file) {
  struct stat  s;
  off_t        size = 0;

  errno = 0;

  fstat(fileno(file), &s);

  if (errno)
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

  AS_UTL_safeWrite(f, o, "AS_UTL_writeFastA", sizeof(char), oi);

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

  AS_UTL_safeWrite(f, s, "AS_UTL_writeFastQ", sizeof(char), sl);
  fprintf(f, "\n");

  fprintf(f, "+\n");
  AS_UTL_safeWrite(f, q, "AS_UTL_writeFastQ", sizeof(char), ql);
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

  AS_UTL_safeWrite(f, s, "AS_UTL_writeFastQ", sizeof(char), sl);
  fprintf(f, "\n");

  fprintf(f, "+\n");
  AS_UTL_safeWrite(f, o, "AS_UTL_writeFastQ", sizeof(char), ql);
  fprintf(f, "\n");

  delete [] o;
}



cftType
compressedFileType(char const *filename) {

  if ((filename == NULL) || (filename[0] == 0) || (strcmp(filename, "-") == 0))
    return(cftSTDIN);

  int32  len = strlen(filename);

  if      ((len > 3) && (strcasecmp(filename + len - 3, ".gz") == 0))
    return(cftGZ);

  else if ((len > 4) && (strcasecmp(filename + len - 4, ".bz2") == 0))
    return(cftBZ2);

  else if ((len > 3) && (strcasecmp(filename + len - 3, ".xz") == 0))
    return(cftXZ);

  else
    return(cftNONE);
}



compressedFileReader::compressedFileReader(const char *filename) {
  char    cmd[FILENAME_MAX];
  int32   len = 0;

  _file     = NULL;
  _filename = duplicateString(filename);
  _pipe     = false;
  _stdi     = false;

  cftType   ft = compressedFileType(_filename);

  if ((ft != cftSTDIN) && (AS_UTL_fileExists(_filename, false, false) == false))
    fprintf(stderr, "ERROR:  Failed to open input file '%s': %s\n", _filename, strerror(errno)), exit(1);

  errno = 0;

  switch (ft) {
    case cftGZ:
      snprintf(cmd, FILENAME_MAX, "gzip -dc '%s'", _filename);
      _file = popen(cmd, "r");
      _pipe = true;
      break;

    case cftBZ2:
      snprintf(cmd, FILENAME_MAX, "bzip2 -dc '%s'", _filename);
      _file = popen(cmd, "r");
      _pipe = true;
      break;

    case cftXZ:
      snprintf(cmd, FILENAME_MAX, "xz -dc '%s'", _filename);
      _file = popen(cmd, "r");
      _pipe = true;

      if (_file == NULL)    //  popen() returns NULL on error.  It does not reliably set errno.
        fprintf(stderr, "ERROR:  Failed to open input file '%s': popen() returned NULL\n", _filename), exit(1);

      errno = 0;
      break;

    case cftSTDIN:
      _file = stdin;
      _stdi = true;
      break;

    default:
      _file = fopen(_filename, "r");
      _pipe = false;
      break;
  }

  if (errno)
    fprintf(stderr, "ERROR:  Failed to open input file '%s': %s\n", _filename, strerror(errno)), exit(1);
}



compressedFileReader::~compressedFileReader() {

  if (_file == NULL)
    return;

  if (_stdi)
    return;

  if (_pipe)
    pclose(_file);
  else
    AS_UTL_closeFile(_file);

  delete [] _filename;
}



compressedFileWriter::compressedFileWriter(const char *filename, int32 level) {
  char   cmd[FILENAME_MAX];
  int32  len = 0;

  _file     = NULL;
  _filename = duplicateString(filename);
  _pipe     = false;
  _stdi     = false;

  cftType   ft = compressedFileType(_filename);

  errno = 0;

  switch (ft) {
    case cftGZ:
      snprintf(cmd, FILENAME_MAX, "gzip -%dc > '%s'", level, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftBZ2:
      snprintf(cmd, FILENAME_MAX, "bzip2 -%dc > '%s'", level, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftXZ:
      snprintf(cmd, FILENAME_MAX, "xz -%dc > '%s'", level, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftSTDIN:
      _file = stdout;
      _stdi = 1;
      break;

    default:
      _file = fopen(_filename, "w");
      _pipe = false;
      break;
  }

  if (errno)
    fprintf(stderr, "ERROR:  Failed to open output file '%s': %s\n", _filename, strerror(errno)), exit(1);
}


compressedFileWriter::~compressedFileWriter() {

  if (_file == NULL)
    return;

  if (_stdi)
    return;

  errno = 0;

  if (_pipe)
    pclose(_file);
  else
    AS_UTL_closeFile(_file);

  if (errno)
    fprintf(stderr, "ERROR:  Failed to cleanly close output file '%s': %s\n", _filename, strerror(errno)), exit(1);

  delete [] _filename;
}
