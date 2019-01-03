
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
 *    src/utility/files.C
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

  if ((ft != cftSTDIN) && (fileExists(_filename) == false))
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
