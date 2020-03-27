
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

class maThreadData {
public:
  maThreadData(trGlobalData *g, uint32 tid) {
    threadID = tid;

    bSeqsLen = 0;
    bSeqsMax = 0;
    bSeqs    = NULL;
  };

  ~maThreadData() {
    delete [] bSeqs;
  };


  uint32      threadID;

  dnaSeq      aSeq;

  uint32      bSeqsLen;
  uint32      bSeqsMax;
  dnaSeq    **bSeqs;
};
