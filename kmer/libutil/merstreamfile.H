#ifndef MERSTREAMFILE_H
#define MERSTREAMFILE_H

#include <sys/types.h>

#include "libbritypes.h"
#include "bitPackedFile.H"

//
//  Uses the merStream to build a compressed merStreamFile.  The
//  merStreamFile contains all the forward mers that the stream
//  returns.  It does not contain position information or deflines.
//  Just the mers, ma'am.
//

//  To get position, sequence number and defline:
//
//  Each mer block can store the position start, and the sequence
//  number.  If we then store the deflines, we get everything!


class merStreamFileBuilder {
public:
  merStreamFileBuilder(u32bit merSize, const char *inputFile, const char *outputFile);
  ~merStreamFileBuilder();

  void                   build(bool beVerbose=false);
  bool                   test(bool beVerbose=false);
private:
  u32bit                _merSize;
  const char           *_inputFile;
  const char           *_outputFile;
};


class merStreamFileReader {
public:
  merStreamFileReader(const char *filename);
  ~merStreamFileReader();

  u32bit                 merSize(void)      { return(_merSize); };
  u32bit                 numberOfMers(void) { return(_numMers); };

  u64bit                 theFMer(void) { return(_theFMer); };
  u64bit                 theRMer(void) { return(_theRMer); };

  u64bit                 thePosition(void)       { return(_thePosition); };
  u64bit                 theSequenceNumber(void) { return(_theSequence); };
  char const            *theDefLine(void)        { return(_theDefline);  };

  bool                   seekToMer(u64bit merNumber);
  bool                   seekToSequence(u64bit seqNumber);

  bool                   nextMer(u32bit skip=0);
private:
  const char           *_inputFile;

  u32bit                _merSize;
  u64bit                _numMers;
  u64bit                _numBlocks;
  u64bit                _numDefs;
  u64bit                _defLength;
  off_t                 _blockFileSize;
  off_t                 _deflineFileSize;
  off_t                 _streamFileSize;
  off_t                 _blkStart;
  off_t                 _defStart;
  off_t                 _strStart;

  u32bit               *_blockSize;
  u32bit               *_blockSequence;
  u64bit               *_blockPosition;
  char                 *_deflineStorage;
  char                **_deflines;

  bitPackedFileReader  *_streamFile;

  //  Things used to decode the stream
  //
  u64bit                _thisBlock;
  u64bit                _thisBlockSize;
  u64bit                _merMask;
  u64bit                _firstMer;

  //  Things returned by the accessors.  These are set/updated by
  //  nextMer().
  //
  u64bit                _theFMer;
  u64bit                _theRMer;
  u64bit                _thePosition;
  u64bit                _theSequence;
  char                 *_theDefline;
};

#endif  //  MERSTREAMFILE_H
