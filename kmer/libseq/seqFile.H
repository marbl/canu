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
  virtual u32bit        find(const char *sequencename) = 0;

  virtual const char   *getSourceName(void)    { return(_filename); };
  virtual const char   *getFileTypeName(void)  { return(_typename); };

  //  Needs to go.  It is unreliable for stdin fasta files.
  virtual u32bit        getNumberOfSequences(void) { return(_numberOfSequences); };

  virtual u32bit        getSequenceLength(u32bit id) = 0;
  virtual bool          getSequence(u32bit id,
                                    char *&h, u32bit &hLen, u32bit &hMax,
                                    char *&s, u32bit &sLen, u32bit &sMax) = 0;
  virtual bool          getSequence(u32bit iid,
                                    u32bit bgn, u32bit end, char *s) = 0;

protected:
  char                 _filename[FILENAME_MAX];
  char                 _typename[FILENAME_MAX];

  u32bit               _numberOfSequences;

  friend class seqFactory;
};

#endif //  SEQFILE_H
