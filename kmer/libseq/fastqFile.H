#ifndef FASTQFILE_H
#define FASTQFILE_H

#include "util++.H"
#include "bio++.H"

#include "seqFile.H"

struct fastqFileHeader {
  u64bit       _magic[2];
  u32bit       _numberOfSequences;      //  Number of sequences in the file
  u32bit       _namesLength;            //  Bytes in the names
  u64bit       _fastqFileSize;          //  st_size  - size of file in bytes
  u64bit       _fastqModificationTime;  //  st_mtime - time of last data modification
  u64bit       _fastqCreationTime;      //  st_ctime - time of last file status change
};


struct fastqFileIndex {
  u64bit       _seqPosition;       //  Position of the sequence in the file
  u32bit       _seqLength;         //  Length of the sequence (no whitespace counted)
};


class fastqFile : public seqFile {
protected:
  fastqFile(const char *filename);
  fastqFile();

public:
  ~fastqFile();

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

  fastqFileHeader    _header;
  fastqFileIndex    *_index;
  char              *_names;

  u32bit             _nextID;         //  Next sequence in the read buffer

  u32bit             _gs_iid;
  u32bit             _gs_pos;

  friend class seqFactory;
};


#endif  //  FASTQFILE_H
