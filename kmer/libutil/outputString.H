#ifndef OUTPUTSTRING_H
#define OUTPUTSTRING_H

//  Used by sim4-cleaner

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class outputString {
private:
  int          _outputPos;
  int          _outputMax;
  char        *_output;

  int          _bufferMax;
  char        *_buffer;

  //  Reallocate the output string (if needed) then append
  //  the string passed in.
  //
  void         _appendString(char *str, int strl) {
    if (strl + 1 + _outputPos > _outputMax) {
      _outputMax = strl + _outputPos + 10240;
      char *n = new char [_outputMax];
      strcpy(n, _output);
      delete [] _output;
      _output = n;
    }
    strcpy(_output + _outputPos, str);
    _outputPos += strl;
  };
public:
  outputString(int size=10240) {
    _outputPos  = 0;
    _outputMax  = size;
    _output     = new char [_outputMax];
    _bufferMax  = size;
    _buffer     = new char [_bufferMax];

    _output[0] = 0;
  };
  ~outputString() {
    delete [] _output;
    delete [] _buffer;
  };

  char *getBuffer(int size=0) {
    if (size && (_bufferMax < size)) {
      delete [] _buffer;
      _buffer = new char [_bufferMax];
    }
    return(_buffer);
  };

  void  appendBuffer(void) {
    _appendString(_buffer, strlen(_buffer));
  };

  void  appendString(char *str) {
    if ((str == 0L) || (str[0] == 0))
      return;
    _appendString(str, strlen(str));
  };

  char *getString(void) {
    return(_output);
  };
};

#endif  //  OUTPUTSTRING_H
