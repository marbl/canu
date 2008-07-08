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

  bool       getRecord(void *record);
  void       putRecord(void *record);

  void       seek(u64bit rec);

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
recordFile::getRecord(void *record) {
  if (_numRecords <= _pos + _rec)
    return(false);
  if (_bfrmax < _rec + 2)
    seek(_pos + _rec);
  memcpy(record, _bfr + _rec * _recordSize, _recordSize);
  _rec++;
  return(true);
}


inline
void
recordFile::putRecord(void *record) {
  if (_bfrmax < _rec + 2)
    seek(_pos + _rec);
  memcpy(_bfr + _rec * _recordSize, record, _recordSize);
  _rec++;
  _numRecords++;
  _bfrDirty = true;
}


#endif  //  RECORDFILE_H



