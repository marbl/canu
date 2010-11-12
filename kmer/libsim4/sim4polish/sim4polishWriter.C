#include "sim4polishWriter.H"

#include "util++.H"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

static
const
char base64[65] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-";

sim4polishWriter::sim4polishWriter(const char *name, sim4polishStyle style, bool hidden) {

  if (hidden) {
    //  We are supposed to be a hidden file.
    strcpy(_otName, "(hidden)");
    _otFile = makeTempFile(NULL);

  } else if ((name == 0L) || ((name[0] == '-') && (name[1] == 0))) {
    //  We are stdout.
    strcpy(_otName, "(stdout)");
    _otFile = stdout;

  } else {
    //  Nope, just a regular ol' file.
    if (strlen(name) > FILENAME_MAX)
      fprintf(stderr, "sim4polishWriter()-- Failed to open '%s' for writing: file name too long.\n",
              name), exit(1);

    strncpy(_otName, name, FILENAME_MAX);

    errno = 0;
    _otFile = fopen(name, "w");
    if (errno)
      fprintf(stderr, "sim4polishWriter()-- Failed to open '%s' for writing: %s\n",
              _otName, strerror(errno)), exit(1);
  }

  _style = style;

  switch (_style) {
    case sim4polishS4DB: s4p_putHeaderS4DB(); break; 
    case sim4polishGFF3: s4p_putHeaderGFF3(); break;
    case sim4polishATAC: s4p_putHeaderATAC(); break; 
  }

  memset(_sourceName,    0, sizeof(char) * 32);
  memset(_matchIDprefix, 0, sizeof(char) * 32);
  memset(_matchIDsalt,   0, sizeof(char) * 8);

  _matchID = 0;

  //  Construct a match ID salt based on the current time and process ID.  We make a 48-bit
  //  number from the combination of process ID and curent time, then convert that to base-64.

  u64bit  saltTime    = (u64bit)getTime();  //  returns a double, fraction of seconds
  u64bit  saltPID     = (u64bit)getpid();

  u64bit  saltInteger = (saltPID << 32) | (saltTime);
  u64bit  saltMask    = u64bitMASK(6);

  _matchIDsalt[0] = base64[saltInteger & saltMask];  saltInteger >>= 6;  //   6 bits
  _matchIDsalt[1] = base64[saltInteger & saltMask];  saltInteger >>= 6;  //  12 bits
  _matchIDsalt[2] = base64[saltInteger & saltMask];  saltInteger >>= 6;  //  18 bits
  _matchIDsalt[3] = base64[saltInteger & saltMask];  saltInteger >>= 6;  //  24 bits
  _matchIDsalt[4] = base64[saltInteger & saltMask];  saltInteger >>= 6;  //  30 bits
  _matchIDsalt[5] = base64[saltInteger & saltMask];  saltInteger >>= 6;  //  36 bits
  _matchIDsalt[6] = base64[saltInteger & saltMask];  saltInteger >>= 6;  //  42 bits
  _matchIDsalt[7] = 0;

#if DEBUG_WRITER
  fprintf(stderr, "SALT: "u64bitFMT" + "u64bitFMT" = %s\n",
          saltPID, saltTime, _matchIDsalt);
#endif
}

void
sim4polishWriter::s4p_putHeaderS4DB() {
  return;
}

void
sim4polishWriter::s4p_putHeaderATAC() {
  return;
}

void
sim4polishWriter::s4p_putHeaderGFF3() {
  fputs( "##gff-version 3\n", _otFile);

  return;
}

sim4polishWriter::~sim4polishWriter() {

  if (strcmp(_otName, "(hidden)") == 0) {
    if (_otFile)
      fprintf(stderr, "sim4polishWriter()-- WARNING:  Hidden output file was lost; surrenderToReader() never called.\n");
  } else {
    errno = 0;
    if (_otFile)
      fclose(_otFile);
    if (errno)
      fprintf(stderr, "sim4polishWriter()-- WARNING:  Failed to close '%s': %s\n",
              _otName, strerror(errno));
  }

  _otFile = NULL;
}


FILE *
sim4polishWriter::surrenderToReader(void) {
  FILE *retval = _otFile;

  _otFile = 0L;

  fflush(retval);
  rewind(retval);
  return(retval);
}


void
sim4polishWriter::setSourceName(const char *sourceName) {

  //  Find the last slash, if any.
  const char *lastSlash = strrchr(sourceName, '/');

  //  If found, advance one letter to the first letter in the name, otherwise
  //  reset lastSlash to the first letter in the sourceName.
  if (lastSlash)
    lastSlash++;
  else
    lastSlash = sourceName;

  if (lastSlash[0] == 0)
    fprintf(stderr, "sim4polishWriter()-- source name is empty, or ends in a '/'; no source name used.\n");

  if (strlen(lastSlash) > 32)
    fprintf(stderr, "sim4polishWriter()-- source name too long, truncating to 31 letters.\n");

  strncpy(_sourceName, lastSlash, 32);
  _sourceName[31] = 0;
}


void
sim4polishWriter::setMatchIDPrefix(const char *prefix) {

  if (strlen(prefix) > 32)
    fprintf(stderr, "sim4polishWriter()-- ID prefix too long, truncating to 31 letters.\n");

  strncpy(_matchIDprefix, prefix, 32);
  _matchIDprefix[31] = 0;
}


void
sim4polishWriter::writeAlignment(sim4polish *out) {
  char *str = 0L;

  switch (_style) {
    case sim4polishS4DB:
      str = out->s4p_polishToStringS4DB();
      break;
    case sim4polishGFF3:
      str = out->s4p_polishToStringGFF3();
      break;
    case sim4polishATAC:
      str = out->s4p_polishToStringATAC();
      break;
  }

  fputs(str, _otFile);

  delete [] str;
}
