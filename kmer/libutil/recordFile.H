#ifndef RECORDFILE_H
#define RECORDFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

//  A file of fixed size records, with an optional header at the
//  start.  Derived from the bitPackedFile at SVN-1533, but heavily
//  modified.  Records can only be added, not updated (probably
//  trivial to fix).  Records must be dense (also probably trivial to
//  fix).

class recordFile {
public:
  recordFile(char const *name,
             u32bit      headerSize,
             u32bit      recordSize,
             bool        append=true);
  ~recordFile();

  void      *header(void) { return(_header); };

  u64bit     numRecords(void) { return(_numRecords); };

  bool       getRecord(void *record, u32bit numRecords=1);
  void       putRecord(void *record, u32bit numRecords=1);

  void       seek(u64bit rec, bool forced=false);

private:
  void       flushDirty(void);

  int       _file;
  char     *_name;

  u64bit    _numRecords;
  u32bit    _recordSize;

  u32bit    _headerSize;
  char     *_header;

  u64bit    _bfrmax;  //  Number of records in the buffer
  char     *_bfr;     //  A chunk of the bitPackedFile in core
  u64bit    _pos;     //  The location this chunk is from (in records)
  u64bit    _rec;     //  The record we're modifying relative to _pos

  bool      _bfrDirty;
  bool      _isReadOnly;
};


inline
bool
recordFile::getRecord(void *record, u32bit numRecords) {
  if ((_numRecords <= _pos + _rec) || (_numRecords < _pos + _rec + numRecords))
    return(false);
  if (_bfrmax < _rec + numRecords + 1)
    seek(_pos + _rec, true);
  memcpy(record, _bfr + _rec * _recordSize, _recordSize * numRecords);
  _rec += numRecords;
  return(true);
}


inline
void
recordFile::putRecord(void *record, u32bit numRecords) {
  if (_bfrmax < _rec + numRecords + 1)
    seek(_pos + _rec, true);
  memcpy(_bfr + _rec * _recordSize, record, _recordSize * numRecords);
  _rec        += numRecords;
  _numRecords += numRecords;
  _bfrDirty = true;
}


#endif  //  RECORDFILE_H



