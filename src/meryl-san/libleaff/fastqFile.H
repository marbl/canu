
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef FASTQFILE_H
#define FASTQFILE_H

#include "seqFile.H"

struct fastqFileHeader {
  uint64       _magic[2];
  uint32       _numberOfSequences;      //  Number of sequences in the file
  uint32       _namesLength;            //  Bytes in the names
  uint64       _fastqFileSize;          //  st_size  - size of file in bytes
  uint64       _fastqModificationTime;  //  st_mtime - time of last data modification
  uint64       _fastqCreationTime;      //  st_ctime - time of last file status change
};


struct fastqFileIndex {
  uint64       _seqPosition;       //  Position of the sequence in the file
  uint32       _seqLength;         //  Length of the sequence (no whitespace counted)
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

  fastqFileHeader    _header;
  fastqFileIndex    *_index;
  char              *_names;

  uint32             _nextID;         //  Next sequence in the read buffer

  uint32             _gs_iid;
  uint32             _gs_pos;

  friend class seqFactory;
};


#endif  //  FASTQFILE_H
