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
  u32bit              find(const char *sequencename);

  u32bit              getSequenceLength(u32bit iid);
  bool                getSequence(u32bit iid,
                                  char *&h, u32bit &hLen, u32bit &hMax,
                                  char *&s, u32bit &sLen, u32bit &sMax);
  bool                getSequence(u32bit iid,
                                  u32bit bgn, u32bit end, char *s);

private:
  void                clear(void);
  bool                loadNextSequence(char *&h, u32bit &hLen, u32bit &hMax,
                                       char *&s, u32bit &sLen, u32bit &sMax);

  readBuffer        *_rb;
  u32bit             _thisIID;

  friend class seqFactory;
};


#endif  //  FASTASTDIN_H
