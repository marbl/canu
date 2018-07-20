
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

#ifndef FASTQSTDIN_H
#define FASTQSTDIN_H

#include "seqFile.H"


class fastqStdin : public seqFile {
protected:
  fastqStdin(const char *filename);
  fastqStdin();

public:
  ~fastqStdin();

protected:
  seqFile            *openFile(const char *filename);

public:
  uint32              getNumberOfSequences(void);

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
  uint32             _nextIID;

  FILE              *_pipe;

  char              *_header;
  uint32             _headerLen;
  uint32             _headerMax;

  char              *_sequence;
  uint32             _sequenceLen;
  uint32             _sequenceMax;

  char              *_quality;
  uint32             _qualityLen;
  uint32             _qualityMax;

  friend class seqFactory;
};


#endif  //  FASTQSTDIN_H
