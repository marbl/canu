#ifndef FASTA_TYPES_H
#define FASTA_TYPES_H

//  Types used in the index
//
typedef off_t  _idxfa_pos;
typedef u32bit _idxfa_len;

//  INDEX TYPES
//
//  Types are currently subsets.  An INDEX_ONLY is a PLUS_IDS is a
//  PLUS_DEFLINES This is used during the openIndex().
//
#define FASTA_INDEX_ANY                0x00000000

#define FASTA_INDEX_ONLY               0x00000001
#define FASTA_INDEX_PLUS_IDS           0x00000002
#define FASTA_INDEX_PLUS_DEFLINES      0x00000003

#define FASTA_INDEX_MASK               0x000000ff

#define FASTA_INDEX_MD5                0x80000000


//  FASTA TYPES
//
#define FASTA_UNDECLARED    0x00000001


struct _idxfa_global {
  u64bit       _magic;
  u32bit       _version;
  u32bit       _indexType;
  u32bit       _fastaType;
  u32bit       _numberOfSequences;      //  Number of sequences in the file
  off_t        _fastaFileSize;          //  st_size  - size of file in bytes
  time_t       _fastaModificationTime;  //  st_mtime - time of last data modification
  time_t       _fastaCreationTime;      //  st_ctime - time of last file status change
  _idxfa_len   _seqlineLength;          //  Length of the sequence lines, if fixed width
  _idxfa_len   _seqendlLength;          //  Length of the sequence line terminators, if fixed width (really, "\n" or "\r\n"?)
  bool         _fixedWidth;             //  All (but last) lines are the same width, and there is no embedded space
  bool         _squeezedSequences;      //  No whitespace in the sequences
  char         _alphabet[256];          //  true/false for each character seen
};

//  20040508 - On the next major change, this structure needs to be
//  rearranged so that 8-byte values are 8-byte aligned (moving
//  _seqStart before _headerLen will do that).  This will fix an
//  optimization on OS-X that tries to align data to their natural
//  boundaries.
//
#ifdef __APPLE__

//  This is the "correct" version that we should use everywhere.
//
struct _idxfa_desc {
  _idxfa_pos   _headerStart;       //  Position of the header in the file
  _idxfa_pos   _seqStart;          //  Position of the sequence in the file
  _idxfa_len   _headerLen;         //  Length of the header, not counting the '\n'
  _idxfa_len   _seqLen;            //  Length of the sequence
};

#else

//  Ick!
//
struct _idxfa_desc {
  _idxfa_pos   _headerStart;       //  Position of the header in the file
  _idxfa_len   _headerLen;         //  Length of the header, not counting the '\n'
  _idxfa_pos   _seqStart;          //  Position of the sequence in the file
  _idxfa_len   _seqLen;            //  Length of the sequence
};

#endif


#endif  //  FASTA_TYPES_H
