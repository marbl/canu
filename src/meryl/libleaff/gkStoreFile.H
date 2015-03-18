#ifndef GKSTOREFILE_H
#define GKSTOREFILE_H

#include "seqFile.H"
#include "gkStore.H"

class gkStoreFile : public seqFile {
protected:
  gkStoreFile(const char *filename);
  gkStoreFile();

public:
  ~gkStoreFile();

protected:
  seqFile      *openFile(const char *name);

public:
  uint32        find(const char *sequencename) {
    assert(0);
    return(0);
  };

  uint32        getSequenceLength(uint32 iid) {
    return(gkp->gkStore_getRead(iid + 1)->gkRead_sequenceLength());
  };

  bool          getSequence(uint32 iid,
                            char *&h, uint32 &hLen, uint32 &hMax,
                            char *&s, uint32 &sLen, uint32 &sMax);
  bool          getSequence(uint32 iid,
                            uint32 bgn, uint32 end, char *s);

private:
  void          clear(void);

  gkStore      *gkp;
  gkReadData    readData;

  friend class seqFactory;
};


#endif  //  GKSTOREFILE_H
