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
             char        mode);
  ~recordFile();

  void      *header(void) { return(_header); };

  u64bit     numRecords(void) { return(_numRecords); };

  //  Read/write records.
  u32bit     getRecord(void *record, u32bit num=1);
  void       putRecord(void *record, u32bit num=1);

  //  Seek to record rec, optionally repositioning the buffer to that
  //  record.
  void       seek(u64bit rec, bool forced=false);

  //  Set an artificial EOF at record rec.
  void       limit(u64bit rec) { _limit = rec; };

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

  u64bit    _limit;   //  An artificial EOF

  u64bit    _pos;     //  The location this chunk is from (in records)
  u64bit    _rec;     //  The record we're modifying relative to _pos

  bool      _bfrDirty;
  bool      _isReadOnly;
};

#endif  //  RECORDFILE_H



