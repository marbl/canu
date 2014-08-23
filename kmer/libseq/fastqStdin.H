#ifndef FASTQSTDIN_H
#define FASTQSTDIN_H

#include "util++.H"
#include "bio++.H"

#include "seqFile.H"


class fastqStdin : public seqFile {
protected:
  fastqStdin(const char *filename);
  fastqStdin();

public:
  ~fastqStdin();

protected:
  seqFile            *openFile(const char *filename);

public:
  uint32              getNumberOfSequences(void);

public:
  uint32              find(const char *sequencename);

  uint32              getSequenceLength(uint32 iid);
  bool                getSequence(uint32 iid,
                                  char *&h, uint32 &hLen, uint32 &hMax,
                                  char *&s, uint32 &sLen, uint32 &sMax);
  bool                getSequence(uint32 iid,
                                  uint32 bgn, uint32 end, char *s);

private:
  void                clear(void);
  bool                loadNextSequence(char *&h, uint32 &hLen, uint32 &hMax,
                                       char *&s, uint32 &sLen, uint32 &sMax);

  readBuffer        *_rb;
  uint32             _nextIID;

  FILE              *_pipe;

  char              *_header;
  uint32             _headerLen;
  uint32             _headerMax;

  char              *_sequence;
  uint32             _sequenceLen;
  uint32             _sequenceMax;

  char              *_quality;
  uint32             _qualityLen;
  uint32             _qualityMax;

  friend class seqFactory;
};


#endif  //  FASTQSTDIN_H
