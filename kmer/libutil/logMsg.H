#ifndef LOGMSG_H
#define LOGMSG_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>

#include "util.h"

class logMsg {
public:
  logMsg(bool toScreen=false) {
    _logLen       = 0;
    _logMax       = 8192;
    _log          = new char [_logMax];
    _toScreenToo  = toScreen;
  };


  ~logMsg() {
    delete [] _log;
  };


  //  Ensure that the string has at least 'moreSpace' available.
  //
  void     resize(u32bit moreSpace) {
    if (_logLen + moreSpace >= _logMax) {
      _logMax += _logMax + moreSpace + 1;
      char *ll = new char [_logMax];
      memcpy(ll, _log, sizeof(char) * _logLen);
      delete [] _log;
      _log = ll;
    }
  };


  //  Add a message to the log, assume the message is less than 8192
  //  bytes.  Would be nice to parse the fmt string (and any args) but
  //  that's a lot of work (and already done if you have vsnprintf.
  //
  //  It warns if you overwrote memory.
  //
  void     add(char const *fmt, ...) {
    va_list   ap;

    resize(8192);

    va_start(ap, fmt);

    if (_toScreenToo)
      vfprintf(stderr, fmt, ap);

    _logLen += vsprintf(_log + _logLen, fmt, ap);

    va_end(ap);

    if (_logLen > _logMax)
      fprintf(stderr,
              "logMsg::add()-- HEY!  I wrote "u32bitFMT" bytes beyond the end of the buffer!\n"
              "logMsg::add()-- This program will probably crash soon....\n",
              _logLen - _logMax);
  };


  //  Dump the message to a file, taking care of errors.
  //
  void     write(int file, char const *name=0L) {
    errno = 0;
    ::write(file, _log, sizeof(char) * _logLen);
    if (errno) {
      fprintf(stderr, "logMsg::write()-- Couldn't write to the log message file '%s': %s\n",
              name ? name : "(unknown)",
              strerror(errno));
      exit(1);
    }
  };

  void     fwrite(FILE *file, char const *name=0L) {
    errno = 0;
    ::fwrite(_log, sizeof(char), _logLen, file);
    if (errno) {
      fprintf(stderr, "logMsg::fwrite()-- Couldn't write to the log message file '%s': %s\n",
              name ? name : "(unknown)",
              strerror(errno));
      exit(1);
    }
  };


private:
  u32bit   _logLen;
  u32bit   _logMax;
  char    *_log;
  bool     _toScreenToo;
};


#endif //  LOGMSG_H
