#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "bio.h"

//  Define this to also get log messages to the screen
//
#define LOG_TO_SCREEN

class logMsg {
public:
  u32bit   theLogLen;
  u32bit   theLogMax;
  char    *theLog;

  logMsg() {
    theLogLen = 0;
    theLogMax = 8192;
    theLog    = new char [theLogMax];
  };

  ~logMsg() {
    delete [] theLog;
  };



  //  Ensure that the string has at least 'moreSpace' available.
  //
  void     resize(u32bit moreSpace) {
    if (theLogLen + moreSpace >= theLogMax) {
      theLogMax += theLogMax + moreSpace;
      char *ll = new char [theLogMax];
      memcpy(ll, theLog, sizeof(char) * theLogLen);
      delete [] theLog;
      theLog = ll;
    }
  };

  //  Add a message to the log, assume the message is less than 1024
  //  bytes.  Would be nice to parse the fmt string (and any args) but
  //  that's a lot of work (and already done if you have vsnprintf.
  //
  //  It warns if you overwrote memory.
  //
  void     add(char *fmt, ...) {
    va_list   ap;
    resize(1024);
    va_start(ap, fmt);
    int ret = vsprintf(theLog + theLogLen, fmt, ap);
    theLogLen += ret;
    if (theLogLen > theLogMax) {
      fprintf(stderr, "WARNING:  logMsg::add() wrote %d bytes beyond the end of the buffer!\n",
              theLogLen - theLogMax);
      fprintf(stderr, "WARNING:  This program will probably crash soon.\n");
    }
    va_end(ap);

#ifdef LOG_TO_SCREEN
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
#endif

  };
};
