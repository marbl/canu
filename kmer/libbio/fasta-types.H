#ifndef FASTA_TYPES_H
#define FASTA_TYPES_H

//  VERSION NUMBERS
//
//  versionNumber 1 contains the simple index, name or defline.  You
//  can't revert to this version.
//
//  versionNumber 2 extends that to contain md5 checksums.  A version
//  1 index is binary compatible with version 2.  You can revert to
//  this version.
//
//  versionNumber 3 changes the order of _idxfa_desc, so that data
//  is natuarlly aligned.  You can revert to this version.
//
//  versionNumber 4 changes flags in _idxfa_global from bool to u32bit
//  so that the size of the structure is more reliable across
//  platforms (in particular, _idxfa_global on Linux/gcc323 was
//  smaller than on OSF1/cxx65, even though the endianess was the
//  same).
//
//  versionNumber 5 adjusts for differences in off_t / time_t between
//  Linux and Linux64.


#define FASTA_MAGICNUMBER   0x4abeadeddadababeULL
#define FASTA_VERSIONNUMBER 5


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


typedef u32bit IID_t;
typedef off_t  _idxfa_pos;
typedef u32bit _idxfa_len;

#if FASTA_VERSIONNUMBER < 5
typedef time_t _idxfa_time
#else
typedef s64bit _idxfa_time;
#endif

#if FASTA_VERSIONNUMBER > 4
#define offtSwap       u64bitSwap
#define idxfaTimeSwap  u64bitSwap
#define idxfaPosSwap   u64bitSwap
#define idxfaLenSwap   u32bitSwap
#endif

struct _idxfa_global {
  u64bit       _magic;
  u32bit       _version;
  u32bit       _indexType;
  u32bit       _fastaType;
  u32bit       _numberOfSequences;      //  Number of sequences in the file
  _idxfa_pos   _fastaFileSize;          //  st_size  - size of file in bytes
  _idxfa_time  _fastaModificationTime;  //  st_mtime - time of last data modification
  _idxfa_time  _fastaCreationTime;      //  st_ctime - time of last file status change
  _idxfa_len   _seqlineLength;          //  Length of the sequence lines, if fixed width
  _idxfa_len   _seqendlLength;          //  Length of the sequence line terminators, if fixed width (really, "\n" or "\r\n"?)
#if FASTA_VERSIONNUMBER < 4
  bool         _fixedWidth;             //  All (but last) lines are the same width, and there is no embedded space
  bool         _squeezedSequences;      //  No whitespace in the sequences
#else
  u32bit       _fixedWidth;             //  All (but last) lines are the same width, and there is no embedded space
  u32bit       _squeezedSequences;      //  No whitespace in the sequences
#endif
  char         _alphabet[256];          //  true/false for each character seen
};



//  Version three rearranges so that 8-byte values are 8-byte aligned
//  (moving _seqStart before _headerLen will do that).  This will fix
//  an optimization on OS-X that tries to align data to their natural
//  boundaries.
//
struct _idxfa_desc {
  _idxfa_pos   _headerStart;       //  Position of the header in the file
#if FASTA_VERSIONNUMBER < 3
  _idxfa_len   _headerLen;         //  Length of the header, not counting the '\n'
  _idxfa_pos   _seqStart;          //  Position of the sequence in the file
#else
  _idxfa_pos   _seqStart;          //  Position of the sequence in the file
  _idxfa_len   _headerLen;         //  Length of the header, not counting the '\n'
#endif
  _idxfa_len   _seqLen;            //  Length of the sequence
};


#endif  //  FASTA_TYPES_H
