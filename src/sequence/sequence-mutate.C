
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sequence/sequence.H"

#include "utility//sequence.H"
#include "mt19937ar.H"


void
doMutate_substitute(double p,
                    char   base, uint8  qual,
                    char  *ns,   uint8 *nq,    uint32 &oo,
                    mutateParameters &mutPar) {

  for (uint32 xx=0; xx<256; xx++) {
    if (p < mutPar.pS[base][xx]) {
      //fprintf(stderr, "sub %c -> %c at pos %u\n", base, xx, oo);
      ns[oo] = xx;
      nq[oo] = 0;
      oo++;
      break;
    }

    p -= mutPar.pS[base][xx];
  }
}



void
doMutate_insert(double p,
                char   base, uint8  qual,
                char  *ns,   uint8 *nq,    uint32 &oo,
                mutateParameters &mutPar) {

  for (uint32 xx=0; xx<256; xx++) {
    if (p < mutPar.pI[xx]) {
      ns[oo] = xx;
      nq[oo] = 0;
      oo++;

      ns[oo] = base;
      nq[oo] = qual;
      oo++;

      break;
    }

    p -= mutPar.pI[xx];
  }
}



void
doMutate(vector<char *> &inputs, mutateParameters &mutPar) {
  mtRandom   MT;

  mutPar.finalize();

  dnaSeq     seq;

  uint32     nMax   = 0;
  char      *nBases = NULL;
  uint8     *nQuals = NULL;

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf     = new dnaSeqFile(inputs[ff]);
    uint64       num    = 0;

    while (sf->loadSequence(seq)) {
      uint32  iPos = 0;   //  Position in the input read
      uint32  oLen = 0;   //  Position in the output read

      //  Resize the output to fit the input string with the expected
      //  number of insert/delete made.

      uint32  expectedLen = seq.length() * (1 + 2 * mutPar.pInsert - mutPar.pDelete);

      resizeArrayPair(nBases, nQuals, 0, nMax, expectedLen);

      //  Over every base, randomly substitue, insert or delete bases.

      for (iPos=0, oLen=0; iPos<seq.length(); iPos++) {
        char   base = seq.bases()[iPos];
        uint8  qual = seq.quals()[iPos];
        double p    = 0.0;

        //  Whoops!  Resize again?
        if (oLen + 2 > nMax)
          resizeArrayPair(nBases, nQuals, oLen, nMax, nMax + 1000);

        //  If a chance of doing something, make a random number.
        if (mutPar.pSubstitute[base] + mutPar.pInsert + mutPar.pD[base] > 0.0)
          p = MT.mtRandomRealClosed();

        //  If a substitution, make it.
        if (p < mutPar.pSubstitute[base]) {
          doMutate_substitute(p, base, qual, nBases, nQuals, oLen, mutPar);
          continue;
        }
        p -= mutPar.pSubstitute[base];


        //  If an insertion, make it.
        if (p < mutPar.pInsert) {
          doMutate_insert(p, base, qual, nBases, nQuals, oLen, mutPar);
          continue;
        }
        p -= mutPar.pInsert;


        //  If a deletion, make it.  Hamm, nothing to do here,
        //  just don't add the base to the output.
        if (p < mutPar.pD[base]) {
          continue;
        }
        p -= mutPar.pD[base];


        //  Otherwise, no change.  Just copy the base.
        nBases[oLen] = base;
        nQuals[oLen] = qual;
        oLen++;
      }


      //  And one more for an inserting at the end of the read.
      double p = MT.mtRandomRealClosed();

      if (p < mutPar.pInsert)
        doMutate_insert(p, 0, 0, nBases, nQuals, oLen, mutPar);


      //  All done changing.  Output the modified read.

      AS_UTL_writeFastA(stdout, nBases, oLen, 0, ">%s\n", seq.name());
    }

    delete sf;
  }
}


