#ifndef FASTAFILE_H
#define FASTAFILE_H

#include "util++.H"
#include "bio++.H"

#include "seqFile.H"

struct fastaFileHeader {
  u64bit       _magic[2];
  u32bit       _numberOfSequences;      //  Number of sequences in the file
  u32bit       _namesLength;            //  Bytes in the names
  u64bit       _fastaFileSize;          //  st_size  - size of file in bytes
  u64bit       _fastaModificationTime;  //  st_mtime - time of last data modification
  u64bit       _fastaCreationTime;      //  st_ctime - time of last file status change
};


struct fastaFileIndex {
  u64bit       _seqPosition;       //  Position of the sequence in the file
  u32bit       _seqLength;         //  Length of the sequence (no whitespace counted)
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
  u32bit              find(const char *sequencename);

  u32bit              getSequenceLength(u32bit iid);
  bool                getSequence(u32bit iid,
                                  char *&h, u32bit &hLen, u32bit &hMax,
                                  char *&s, u32bit &sLen, u32bit &sMax);
  bool                getSequence(u32bit iid,
                                  u32bit bgn, u32bit end, char *s);

private:
  void                clear(void);
  void                loadIndex(char *indexname);
  void                constructIndex(void);

  readBuffer        *_rb;

  fastaFileHeader    _header;
  fastaFileIndex    *_index;
  char              *_names;

  u32bit             _nextID;         //  Next sequence in the read buffer

  u32bit             _gs_iid;
  u32bit             _gs_pos;

  friend class seqFactory;
};


#endif  //  FASTAFILE_H
