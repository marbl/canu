#ifndef FASTA_H
#define FASTA_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <errno.h>
#include <ctype.h>
#include "libbritypes.H"


//
//  We use this in reverseComplementSequence().  It's really defined
//  in libbri.H, but that also includes this file.
//
extern const unsigned char  complementSymbol[256];

typedef off_t  _idxfa_pos;
typedef u32bit _idxfa_len;

struct _idxfa_desc {
  _idxfa_pos   _headerStart;       //  Position of the header in the file
  _idxfa_len   _headerLen;         //  Length of the header, not counting the '\n'
  _idxfa_pos   _seqStart;          //  Position of the sequence in the file
  _idxfa_len   _seqLen;            //  Length of the sequence
};


class FastABuffer {
private:
  u32bit         _idx;

  u32bit         _headerLen;
  u32bit         _headerMax;
  unsigned char *_header;

  u32bit         _seqLen;
  u32bit         _seqMax;
  unsigned char *_seq;

  friend class FastA;
public:
  FastABuffer() {
    _idx       = u32bitMAX;

    _headerLen = 0;
    _headerMax = 1024;
    _header    = new unsigned char [_headerMax];
    _header[0] = 0;

    _seqLen = 0;
    _seqMax = 4 * 1024;
    _seq    = new unsigned char [_seqMax];
    _seq[0] = 0;
  };

  ~FastABuffer() {
    delete [] _header;
    delete [] _seq;
  };

  void               clear(void) {
    _idx = u32bitMAX;

    _headerLen = 0;
    _header[0] = 0;

    _seqLen = 0;
    _seq[0] = 0;
  };

  void               pushHeader(unsigned char c) {
    _header[_headerLen++] = c;

    if (_headerLen == (_headerMax-1)) {
      unsigned char *nh = new unsigned char [_headerMax*2];
      memcpy((void *)nh, (void *)_header, sizeof(unsigned char) * _headerMax);
      delete [] _header;
      _header     = nh;
      _headerMax *= 2;
    }
  };

  void               pushSequence(unsigned char c) {
    _seq[_seqLen++] = c;

    if (_seqLen == (_seqMax-1)) {
      unsigned char *ns = new unsigned char [_seqMax*2];
      memcpy((void *)ns, (void *)_seq, sizeof(unsigned char) * _seqMax);
      delete [] _seq;
      _seq     = ns;
      _seqMax *= 2;
    }
  };

  void               terminateSequence(void) {
    _seq[_seqLen] = 0;
  };

  void               reverseComplementSequence(void) {
    unsigned char  *s = _seq;
    unsigned char  *e = _seq + _seqLen - 1;
    unsigned char   t;
    unsigned int    c = _seqLen / 2;

    while (c--) {
      t = complementSymbol[*s];
      *(s++) = complementSymbol[*e];
      *(e--) = t;
    }

    if (s == e)
      *s = complementSymbol[*s];
  };


  unsigned char     *header(void)         { return(_header); };
  u32bit             headerLength(void)   { return(_headerLen); };

  unsigned char     *sequence(void)       { return(_seq); };
  u32bit             sequenceLength(void) { return(_seqLen); };

  u32bit             index(void)          { return(_idx); };
};





class FastA {
public:
  FastA(char *name, bool randomAccess=false, bool beVerbose=false);
  ~FastA();

  //  Sequential access methods
  //
  void          first(FastABuffer &b);
  bool          next(FastABuffer &b);
  bool          eof(void)   { return(_eof); };

  //  Random access methods
  //
  bool          seek(FastABuffer &b, u32bit id);
  u32bit        numberOfSequences(void) { return(_numSeqs); };

  u32bit        sequenceLength(u32bit id) {
    if (_numSeqs == 0) {
      fprintf(stderr, "%s: FastA::sequenceLength() called, but file is not random-access!\n", _fastaname);
      return(0);
    }
    return(_descr[id]._seqLen);
  };

  u32bit        headerLength(u32bit id) {
    if (_numSeqs == 0) {
      fprintf(stderr, "%s: FastA::headerLength() called, but file is not random-access!\n", _fastaname);
      return(0);
    }
    return(_descr[id]._headerLen);
  };
private:
  int            _file;
  bool           _eof;

  u32bit         _numSeqs;    //  If random-access, _numSeqs will be non-zero.
  u32bit         _maxSeqs;
  _idxfa_desc   *_descr;

  size_t         _bufferLen;
  size_t         _bufferPos;
  size_t         _bufferMax;
  char unsigned *_buffer;

  char          *_fastaname;

  void            fillBuffer(void) {
    if (_bufferPos >= _bufferLen) {
      errno = 0;
      _bufferPos = 0;
      errno = 0;
      _bufferLen = (u32bit)read(_file, _buffer, _bufferMax);
      if (errno) {
        fprintf(stderr, "%s: FastA::fillBuffer() got '%s'\n", _fastaname, strerror(errno));
        exit(1);
      }

      if (_bufferLen == 0)
        _eof = true;
    }
  };

  unsigned char   theCharacter(void) {
    return(_buffer[_bufferPos]);
  };

  void            nextCharacter(void) {
    _bufferPos++;
    fillBuffer();
  };
};

#endif  //  FASTA_H
