
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2008-JUL-08 to 2014-APR-11
 *      are Copyright 2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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
             uint32      headerSize,
             uint32      recordSize,
             char        mode);
  ~recordFile();

  void      *header(void) { return(_header); };

  uint64     numRecords(void) { return(_numRecords); };

  //  Read/write records.
  uint32     getRecord(void *record, uint32 num=1);
  void       putRecord(void *record, uint32 num=1);

  //  Seek to record rec, optionally repositioning the buffer to that
  //  record.
  void       seek(uint64 rec, bool forced=false);

  //  Set an artificial EOF at record rec.
  void       limit(uint64 rec) { _limit = rec; };

private:
  void       flushDirty(void);

  int       _file;
  char     *_name;

  uint64    _numRecords;
  uint32    _recordSize;

  uint32    _headerSize;
  char     *_header;

  uint64    _bfrmax;  //  Number of records in the buffer
  char     *_bfr;     //  A chunk of the bitPackedFile in core

  uint64    _limit;   //  An artificial EOF

  uint64    _pos;     //  The location this chunk is from (in records)
  uint64    _rec;     //  The record we're modifying relative to _pos

  bool      _bfrDirty;
  bool      _isReadOnly;
};

#endif  //  RECORDFILE_H



