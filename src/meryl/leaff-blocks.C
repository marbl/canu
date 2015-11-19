
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
 *  This file is derived from:
 *
 *    kmer/leaff/blocks.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2009-FEB-07 to 2014-APR-11
 *      are Copyright 2009,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "seqCache.H"

void
dumpBlocks(char *filename) {
seqCache      *F     = 0L;
  seqInCore   *S     = 0L;

  bool                  V[256] = {0};

  for (uint32 i=0; i<256; i++)
    V[i] = false;

  V['n'] = true;
  V['N'] = true;

  F = new seqCache(filename);

  for (uint32 s=0; s<F->getNumberOfSequences(); s++) {
    seqInCore *S = F->getSequenceInCore(s);

    uint32  len    = S->sequenceLength();
    char    begseq = S->sequence()[0];
    bool    nnn    = V[begseq];
    uint32  begpos = 0;
    uint32  pos    = 0;

    for (pos=0; pos<len; pos++) {
      char seq = S->sequence()[pos];

      if (nnn != V[seq]) {
        fprintf(stdout, "%c "F_U32" "F_U32" "F_U32" "F_U32"\n",
                begseq, s, begpos, pos, pos - begpos);
        nnn = V[seq];
        begpos = pos;
        begseq = seq;
      }
    }

    fprintf(stdout, "%c "F_U32" "F_U32" "F_U32" "F_U32"\n",
            begseq, s, begpos, pos, pos - begpos);
    fprintf(stdout, ". "F_U32" "F_U32" "F_U32"\n", s, pos, 0);

    delete S;
  }

  delete F;
}


