#ifndef FASTASTDIN_H
#define FASTASTDIN_H

#include "util++.H"
#include "bio++.H"

#include "seqFile.H"


class fastaStdin : public seqFile {
protected:
  fastaStdin(const char *filename);
  fastaStdin();

public:
  ~fastaStdin();

protected:
  seqFile            *openFile(const char *filename);

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
  uint32             _thisIID;

  friend class seqFactory;
};


#endif  //  FASTASTDIN_H
