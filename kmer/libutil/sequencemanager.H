#ifndef SEQUENCE_MANAGER_H
#define SEQUENCE_MANAGER_H

#include "libbri.H"

//  The sequence manager will handle loading of sequence from a
//  multi-FastA.  It caches and is lazy.  (Always lazy)
//
class SequenceManager {
private:
  //  The FastA, and how many sequences are in it.
  //
  FastA            *_seqFile;
  u32bit            _max;

  //  In core buffers.
  //
  //  _inCore is a map from a FastA index to a buffer.  If the map is
  //  ~0, then it's not loaded.
  //
  u32bit            _numBuffers;
  FastABuffer      *_buf;
  u32bit           *_bufSeqID;
  u32bit            _nextAvail;
  u32bit           *_inCore;

  void           loadSequence(unsigned int i) {

    //  If it's already loaded, just return.
    //
    if (_inCore[i] < _max)
      return;

    //  If a sequence is in the next buffer, unload it.
    //
    if (_bufSeqID[_nextAvail] < _max)
      _inCore[ _bufSeqID[_nextAvail] ] = u32bitMAX;

    //  Assign the next avail buffer to this sequence.
    //
    _inCore[i]            = _nextAvail;
    _bufSeqID[_nextAvail] = i;

    //  Load the sequence into the buffer.
    //
    _seqFile->seek(_buf[_nextAvail], i);

    //  Move to the next buffer, maybe wrapping around
    //
    _nextAvail++;
    if (_nextAvail >= _numBuffers)
      _nextAvail = 0;
  };
public:
  SequenceManager(char *filename, unsigned int keep=128) {
    _seqFile    = new FastA(filename, true);
    _max        = _seqFile->numberOfSequences();

    _numBuffers = keep;
    _buf        = new FastABuffer [keep];
    _bufSeqID   = new u32bit [keep];
    _nextAvail  = 0;
    _inCore     = new u32bit [_max];

    u32bit i;
    for (i=0; i<_numBuffers; i++)
      _bufSeqID[i] = u32bitMAX;
    for (i=0; i<_max; i++)
      _inCore[i]   = u32bitMAX;
  };

  SequenceManager(unsigned char *seq, int low, int high) {
    _seqFile    = 0L;
    _max        = 1;

    _numBuffers = 1;
    _buf        = new FastABuffer [1];
    _bufSeqID   = 0L;
    _nextAvail  = 0;
    _inCore     = new u32bit [1];

    _inCore[0] = 0;  //  sequence 0 is in buffer 0

    if (high < low) {
      int t = low;
      low   = high;
      high  = t;
    }

    if (low < 0)
      low = 0;

    for (int i=low; seq[i] && i<high; i++)
      _buf->pushSequence(seq[i]);

    _buf->terminateSequence();
  }

  ~SequenceManager() {
    delete    _seqFile;
    delete [] _buf;
    delete [] _inCore;
    delete [] _bufSeqID;
  };

  unsigned char  *getSequence(unsigned int i) {
    if (i >= _max) {
      return(0L);
    } else {
      loadSequence(i);
      return(_buf[_inCore[i]].sequence());
    }
  };
  u32bit          getSequenceLength(unsigned int i) {
    if (i >= _max) {
      return(0);
    } else {
      loadSequence(i);
      return(_buf[_inCore[i]].sequenceLength());
    }
  };
  unsigned char  *getHeader(unsigned int i) {
    if (i >= _max) {
      return(0L);
    } else {
      loadSequence(i);
      return(_buf[_inCore[i]].header());
    }
  };
  u32bit          getHeaderLength(unsigned int i) {
    if (i >= _max) {
      return(0);
    } else {
      loadSequence(i);
      return(_buf[_inCore[i]].headerLength());
    }
  };

  unsigned int    getNumberOfSequences(void) {
    return(_max);
  };
};

#endif  //  SEQUENCE_MANAGER_H
