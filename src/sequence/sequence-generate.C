
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
 *    src/sequence/sequence.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sequence/sequence.H"

#include "files.H"
#include "mt19937ar.H"



void
doGenerate(generateParameters &genPar) {
  mtRandom   MT;

  uint64  nSeqs  = 0;
  uint64  nBases = 0;

  uint64  seqLen = 0;
  uint64  seqMax = 65536;
  char   *seq    = new char  [seqMax + 1];
  uint8  *qlt    = new uint8 [seqMax + 1];

  double Athresh = genPar.aFreq;
  double Cthresh = genPar.aFreq + genPar.cFreq;
  double Gthresh = genPar.aFreq + genPar.cFreq + genPar.gFreq;
  double Tthresh = genPar.aFreq + genPar.cFreq + genPar.gFreq + genPar.tFreq;

  while ((nSeqs  < genPar.nSeqs) &&
         (nBases < genPar.nBases)) {
    double   len = MT.mtRandomGaussian(genPar.gMean, genPar.gStdDev);

    while (len < -0.5)
      len = MT.mtRandomGaussian(genPar.gMean, genPar.gStdDev);

    seqLen = (uint64)round(len);

    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, 0, seqMax, seqLen+1, resizeArray_doNothing);

    for (uint64 ii=0; ii<seqLen; ii++) {
      double  bp = MT.mtRandomRealOpen();

      if        (bp < Athresh) {
        seq[ii] = 'A';
        qlt[ii] = 0;

      } else if (bp < Cthresh) {
        seq[ii] = 'C';
        qlt[ii] = 0;

      } else if (bp < Gthresh) {
        seq[ii] = 'G';
        qlt[ii] = 0;

      } else {
        seq[ii] = 'T';
        qlt[ii] = 0;
      }
    }

    seq[seqLen] = 0;
    qlt[seqLen] = 0;

    AS_UTL_writeFastA(stdout,
                      seq, seqLen, 0,
                      ">random%08 " F_U64P "\n", nSeqs);

    nSeqs  += 1;
    nBases += seqLen;
  }

  delete [] seq;
  delete [] qlt;
}

