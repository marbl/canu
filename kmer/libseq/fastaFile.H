#ifndef FASTAFILE_H
#define FASTAFILE_H

#include "util++.H"
#include "bio++.H"

#include "seqFile.H"

struct fastaFileHeader {
  uint64       _magic[2];
  uint32       _numberOfSequences;      //  Number of sequences in the file
  uint32       _namesLength;            //  Bytes in the names
  uint64       _fastaFileSize;          //  st_size  - size of file in bytes
  uint64       _fastaModificationTime;  //  st_mtime - time of last data modification
  uint64       _fastaCreationTime;      //  st_ctime - time of last file status change
};


struct fastaFileIndex {
  uint64       _seqPosition;       //  Position of the sequence in the file
  uint32       _seqLength;         //  Length of the sequence (no whitespace counted)
};


class fastaFile : public seqFile {
protected:
  fastaFile(const char *filename);
  fastaFile();

public:
  ~fastaFile();

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
  void                loadIndex(char *indexname);
  void                constructIndex(void);

  readBuffer        *_rb;

  fastaFileHeader    _header;
  fastaFileIndex    *_index;
  char              *_names;

  uint32             _nextID;         //  Next sequence in the read buffer

  uint32             _gs_iid;
  uint32             _gs_pos;

  friend class seqFactory;
};


#endif  //  FASTAFILE_H
