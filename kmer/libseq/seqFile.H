#ifndef SEQFILE_H
#define SEQFILE_H

#include "util.h"

//  General flow of the constructors is:
//    Clear all data
//    Open the file
//    Set _filename, _typename
//    Read/build the index structure
//    Position the file to the first read
//    Set _numberOfSequences (IMPORTANT, and subtle)

class seqFile {
protected:
  seqFile(const char *filename) {};
  seqFile() {};

public:
  virtual ~seqFile() {};

protected:
  virtual seqFile      *openFile(const char *filename) = 0;

public:
  virtual const char   *getSourceName(void)    { return(_filename); };
  virtual const char   *getFileTypeName(void)  { return(_typename); };

  virtual bool          randomAccessSupported(void) { return(_randomAccessSupported); };

  virtual uint32        getNumberOfSequences(void) { return(_numberOfSequences); };

public:
  virtual uint32        find(const char *sequencename) = 0;

  virtual uint32        getSequenceLength(uint32 id) = 0;
  virtual bool          getSequence(uint32 id,
                                    char *&h, uint32 &hLen, uint32 &hMax,
                                    char *&s, uint32 &sLen, uint32 &sMax) = 0;
  virtual bool          getSequence(uint32 iid,
                                    uint32 bgn, uint32 end, char *s) = 0;

protected:
  char                 _filename[FILENAME_MAX];
  char                 _typename[FILENAME_MAX];

  bool                 _randomAccessSupported;

  uint32               _numberOfSequences;

  friend class seqFactory;
};

#endif //  SEQFILE_H
