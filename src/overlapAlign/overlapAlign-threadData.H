
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
 *    Brian P. Walenz beginning on 2019-APR-22
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
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
