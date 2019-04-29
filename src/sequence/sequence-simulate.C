
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

#include "utility/sequence.H"
#include "sequence/sequence.H"
#include "mt19937ar.H"

#include <vector>

using namespace std;



void
doSimulate_loadSequences(simulateParameters  &simPar,
                         vector<dnaSeq *>    &seqs,
                         uint64              &seqLen) {
  fprintf(stderr, "Loading sequences from '%s'\n", simPar.genomeName);

  dnaSeq           *seq    = new dnaSeq;
  dnaSeqFile       *sf     = new dnaSeqFile(simPar.genomeName);

  while (sf->loadSequence(*seq)) {
    //fprintf(stderr, "  %9lu '%s'\n", seq->length(), seq->name());

    seqLen += seq->length();

    seqs.push_back(seq);

    seq = new dnaSeq;
  }

  delete seq;
  delete sf;

  fprintf(stderr, "Loaded %lu sequences.\n", seqs.size());
}



void
doSimulate_extract(simulateParameters &simPar,
                   vector<dnaSeq *>   &seqs,
                   uint64              seqLen,        //  Total length of all sequences
                   mtRandom           &mt,
                   uint64              nReadsMax,
                   uint64              nBasesMax) {

  uint64  nReads = 0;
  uint64  nBases = 0;

  uint32  rLen = 0;
  uint32  rMax = 1048576;
  char   *r    = new char [rMax];

  while ((nReads < nReadsMax) &&
         (nBases < nBasesMax)) {

    //  Based on the input length distribution, generate a random read length.
    //
    uint32  readLength    = simPar.dist.getValue(mt.mtRandomRealOpen());

    resizeArray(r, 0, rMax, readLength+1, resizeArray_doNothing);

    //  For normal non-circular references, we cannot start a read in the
    //  last few bases of the sequence.  But we can for circular references.
    //  normalSeqDiff adjusts the reference length (at certain times) to
    //  account for this.
    //
    uint32  normalSeqDiff = (simPar.circular == true) ? (0) : (readLength);

    //  Compute a position in the sequences.
    //
    //  If we're circular, any position is valid.
    //
    //  If not, we cannot start a new read in the last N bases of the
    //  sequence.  And if the read length is longer than the reference
    //  length, we cannot get any read out of this sequence.  Thus, we need
    //  to explicitly sum the lengths of sequences for each read length.

    uint64  sl = 0;

    if (simPar.circular == false) {
      for (uint32 ss=0; ss < seqs.size(); ss++)
        if (readLength <= seqs[ss]->length())
          sl += seqs[ss]->length() - readLength + 1;
    } else {
      sl = seqLen;
    }

    uint64  position = (uint64)floor(mt.mtRandomRealOpen() * sl);

    //  Search for the sequence that has bases that start at 'position'.

    for (uint32 ss=0; ss < seqs.size(); ss++) {

      //  Skip to the next sequence if the desired read length is longer
      //  than the sequence length.
      if (seqs[ss]->length() < readLength) {
        continue;
      }

      //  Skip to the next sequence if the desired start position is after
      //  this sequence.
      if (seqs[ss]->length() - normalSeqDiff < position) {
        position -= seqs[ss]->length() - normalSeqDiff;
        continue;
      }

      //  Otherwise, we're in the correct sequence.  If the reference is
      //  linear, we're guaranteed to have all the bases for the read in a
      //  contiguous block.  But if circular, we might need to wrap around.

      if (position + readLength <= seqs[ss]->length()) {
        memcpy(r, seqs[ss]->bases() + position, sizeof(char) * readLength);
      }

      else {
        uint32  l1 = seqs[ss]->length() - position;
        uint32  l2 = readLength - l1;

        assert(simPar.circular == true);
        assert(position < seqs[ss]->length());
        assert(l2 <= seqs[ss]->length());

        memcpy(r,      seqs[ss]->bases() + position, sizeof(char) * l1);
        memcpy(r + l1, seqs[ss]->bases(),            sizeof(char) * l2);
      }

      r[readLength] = 0;

      fprintf(stdout, ">read=%u,position=%u,length=%u,%s\n", nReads+1, position, readLength, seqs[ss]->name());
      fprintf(stdout, "%s\n", r);

      break;
    }

    //  Account for the read we just emitted.

    nReads += 1;
    nBases += readLength;
  }
}



void
doSimulate(vector<char *>     &inputs,
           simulateParameters &simPar) {

  //  Load the genome sequences.

  vector<dnaSeq *>  seqs;
  uint64            seqLen = 0;

  doSimulate_loadSequences(simPar, seqs, seqLen);

  //  Decide how many reads or bases to make.

  mtRandom   mt;

  uint64  nReads = 0, nReadsMax = UINT64_MAX;
  uint64  nBases = 0, nBasesMax = UINT64_MAX;

  if (simPar.desiredCoverage > 0)
    nBasesMax = simPar.desiredCoverage * simPar.genomeSize;

  if (simPar.desiredNumReads > 0)
    nReadsMax = simPar.desiredNumReads;

  if (simPar.desiredNumBases > 0)
    nBasesMax = simPar.desiredNumBases;

  //  Make reads!

  //
  //  Can't support desiredMinLength and desiredMaxLength without
  //  changing sample.H
  //

  doSimulate_extract(simPar, seqs, seqLen, mt, nReadsMax, nBasesMax);
}
