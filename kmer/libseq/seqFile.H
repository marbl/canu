#ifndef SEQFILE_H
#define SEQFILE_H

#include "util.h"

class seqFile {
protected:
  seqFile(const char *filename) {};

public:
  seqFile() {};
  virtual ~seqFile() {};

  virtual seqFile      *openFile(const char *filename) = 0;

  virtual u32bit        find(const char *sequencename) = 0;

  virtual const char   *getSourceName(void)    { return(_filename); };
  virtual const char   *getFileTypeName(void)  { return(_typename); };

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
};

#endif //  SEQFILE_H
