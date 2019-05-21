
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

#include "utility/sequence.H"



uint64
doExtract_compress(char    *outputString,
                   uint64   outputStringLen) {

  if (outputStringLen == 0)
    return(0);

  uint64  cc = 0;
  uint64  rr = 1;

  while (rr < outputStringLen) {
    if (outputString[cc] == outputString[rr])
      rr++;
    else
      outputString[++cc] = outputString[rr++];
  }

  outputString[++cc] = 0;

  return(cc);

}


void
doExtract(vector<char *>    &inputs,
          extractParameters &extPar) {

  char            C[256] = {0};
  char            U[256] = {0};
  char            L[256] = {0};

  uint32          nameMax = 0;
  char           *name    = NULL;
  uint64          seqMax  = 0;
  char           *seq     = NULL;
  uint8          *qlt     = NULL;
  uint64          seqLen  = 0;

  uint64  outputStringLen = 0;
  uint64  outputStringMax = 0;
  char   *outputString    = NULL;

  //  Initialize complement. toUpper and toLower arrays.

  C['a'] = 't';  U['a'] = 'A';  L['a'] = 'a';
  C['c'] = 'g';  U['c'] = 'C';  L['c'] = 'c';
  C['g'] = 'c';  U['g'] = 'G';  L['g'] = 'g';
  C['t'] = 'a';  U['t'] = 'T';  L['t'] = 't';

  C['A'] = 'T';  U['A'] = 'A';  L['A'] = 'a';
  C['C'] = 'G';  U['C'] = 'C';  L['C'] = 'c';
  C['G'] = 'C';  U['G'] = 'G';  L['G'] = 'g';
  C['T'] = 'A';  U['T'] = 'T';  L['T'] = 't';



  for (uint32 fi=0; fi<inputs.size(); fi++) {
    dnaSeqFile  *sf   = new dnaSeqFile(inputs[fi], true);

    //  Allocate a string big enough to hold the largest output.
    //
    //  Later, maybe, we can analyze the bases to output and make this exactly the correct size.

    uint64  maxStringLength = 0;

    for (uint32 ss=0; ss<sf->numberOfSequences(); ss++)
      maxStringLength = max(maxStringLength, sf->sequenceLength(ss));

    resizeArray(outputString, 0, outputStringMax, maxStringLength + 1);

    //fprintf(stderr, "seqs - length %u first %u %u\n", extPar.seqsBgn.size(), extPar.seqsBgn[0], extPar.seqsEnd[0]);

    for (uint32 si=0; si<extPar.seqsBgn.size(); si++) {
      uint64  sbgn = extPar.seqsBgn[si];
      uint64  send = extPar.seqsEnd[si];

      sbgn = min(sbgn, sf->numberOfSequences());
      send = min(send, sf->numberOfSequences());

      //fprintf(stderr, "sbgn %u send %u\n", sbgn, send);

      for (uint32 ss=sbgn; ss<send; ss++) {
        uint64  seqLen = sf->sequenceLength(ss);

        //fprintf(stderr, "lens - length %u first %u %u\n", extPar.lensBgn.size(), extPar.lensBgn[0], extPar.lensEnd[0]);

        for (uint32 li=0; li<extPar.lensBgn.size(); li++) {
          uint64  lmin = extPar.lensBgn[li];
          uint64  lmax = extPar.lensEnd[li];

          if ((seqLen < lmin) ||
              (lmax < seqLen))
            seqLen = UINT64_MAX;
        }

        if (seqLen == UINT64_MAX)
          continue;

        if (sf->findSequence(ss) == false) {
          //fprintf(stderr, "Failed to find sequence #%u in file '%s'\n", ss, inputs[fi]);
          continue;
        }

        if (sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen) == false) {
          //fprintf(stderr, "Failed to load sequence #%u in file '%s'\n", ss, inputs[fi]);
          continue;
        }

        //fprintf(stderr, "base - length %u first %u %u\n", extPar.baseBgn.size(), extPar.baseBgn[0], extPar.baseEnd[0]);

        outputStringLen = 0;

        for (uint32 bi=0; bi<extPar.baseBgn.size(); bi++) {
          uint64  bbgn = extPar.baseBgn[bi];
          uint64  bend = extPar.baseEnd[bi];

          bbgn = min(bbgn, seqLen);
          bend = min(bend, seqLen);

          //fprintf(stderr, "base - seqLen %u base[%u] %u %u limited %u %u\n", seqLen, bi, extPar.baseBgn[bi], extPar.baseEnd[bi], bbgn, bend);

          if (bbgn == bend)
            continue;

          memcpy(outputString + outputStringLen, seq + bbgn, bend - bbgn);

          outputStringLen += bend - bbgn;
        }

        outputString[outputStringLen] = 0;

        if (extPar.asReverse)
          reverse(outputString, outputString + outputStringLen);

        if (extPar.asComplement)
          for (uint32 ii=0; ii<outputStringLen; ii++)
            outputString[ii] = C[outputString[ii]];

        if (extPar.asUpperCase)
          for (uint32 ii=0; ii<outputStringLen; ii++)
            outputString[ii] = U[outputString[ii]];

        if (extPar.asLowerCase)
          for (uint32 ii=0; ii<outputStringLen; ii++)
            outputString[ii] = L[outputString[ii]];

        if (extPar.asCompressed)
          outputStringLen = doExtract_compress(outputString, outputStringLen);

        fprintf(stdout, ">%s\n%s\n", name, outputString);
      }
    }

    //  Done with this file.  Get rid of it.

    delete sf;
  }

  //  Cleanup.

  delete [] name;
  delete [] seq;
  delete [] qlt;
  delete [] outputString;
}





