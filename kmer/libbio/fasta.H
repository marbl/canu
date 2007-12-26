#ifndef FASTA_H
#define FASTA_H

#include "seq.H"

class seqInCore;
class seqOnDisk;

#define FASTA_MAGICNUMBER   0x4abeadeddadababeULL
#define FASTA_VERSIONNUMBER 0x00000010

#define FASTA_INDEX_ANY                0x00000000
#define FASTA_INDEX_ONLY               0x00000001
#define FASTA_INDEX_PLUS_IDS           0x00000002
#define FASTA_INDEX_PLUS_DEFLINES      0x00000003

struct _idxfa_global {
  u64bit       _magic;
  u32bit       _version;
  u32bit       _indexType;
  u32bit       _numberOfSequences;      //  Number of sequences in the file
  u64bit       _fastaFileSize;          //  st_size  - size of file in bytes
  u64bit       _fastaModificationTime;  //  st_mtime - time of last data modification
  u64bit       _fastaCreationTime;      //  st_ctime - time of last file status change
  u32bit       _seqlineLength;          //  Length of the sequence lines, if fixed width
  u32bit       _seqendlLength;          //  Length of the sequence line terminators, if fixed width (really, "\n" or "\r\n"?)
  u32bit       _fixedWidth;             //  All (but last) lines are the same width, and there is no embedded space
  u32bit       _squeezedSequences;      //  No whitespace in the sequences
  char         _alphabet[256];          //  true/false for each character seen
};

struct _idxfa_desc {
  u64bit       _headerStart;       //  Position of the header in the file
  u64bit       _seqStart;          //  Position of the sequence in the file
  u32bit       _headerLen;         //  Length of the header, not counting the '\n'
  u32bit       _seqLen;            //  Length of the sequence
};


class fastaFile : public seqFile {
protected:
  fastaFile(char const *filename);

public:
  fastaFile();
  ~fastaFile();

  const char      *fileTypeName(void) { return("FastA"); };

  seqFile         *openFile(const char *name);

  void             openIndex(u32bit indextypetoload=FASTA_INDEX_ANY);

  bool             rewind(void)   { _curIID = 0; return(_filebuffer->seek(0)); };
  bool             eof(void)      { return(_filebuffer->eof()); };

  seqInCore       *getSequenceInCore(void);
  seqOnDisk       *getSequenceOnDisk(void);

  u64bit           timeStamp(void) { return(_fileTimeStamp); };

  u32bit           currentIID(void) { return(_curIID); };

  bool             find(seqIID  iid);
  bool             find(char   *uid);

  u32bit           sequenceLength(seqIID iid) {
    return((_isIndexed) ? _theSeqs[iid]._seqLen : 0);
  };

  u32bit           nameLength(seqIID iid) {
    return((_isIndexed) ? _theSeqs[iid]._headerLen : 0);
  };

  u32bit           getNumberOfSequences(void) {
    return((_isIndexed) ? _theDesc._numberOfSequences : 0);
  }

  char const      *getSourceName(void) {
    return(_filename);
  };

  void             printDescription(FILE *out, char *name);

  bool             isIndexed(void)  { return(_isIndexed); };
  bool             isSqueezed(void) { return(_theDesc._squeezedSequences); };

private:
  void             createIndex(u32bit indextype);

  //  probably should go away
  void             optimizeRandomAccess(void);

  _idxfa_global   _theDesc;

  u64bit          _fileTimeStamp;

  u32bit          _numSeqs;  //  for create
  u32bit          _maxSeqs;  //  for create
  _idxfa_desc    *_theSeqs;

  u32bit          _theNamesLen;
  u32bit          _theNamesMax;  //  for create
  char           *_theNames;

  readBuffer     *_filebuffer;
  char            _filename[FILENAME_MAX];
  char            _indexname[FILENAME_MAX];

  //u32bit          _curIID;

  bool            _isStreamInput;      //  Forced to be not random access
  bool            _isIndexed;          //  Set if an index is open
  bool            _isRandomAccessOpt;  //  Set if we resized the readBuffer
};


#endif  //  FASTA_H
