 
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

#include "utility/sequence.H"
#include "mt19937ar.H"



class seqEntry {
public:
  seqEntry(mtRandom &MT, uint64 pos_) {
    pos = pos_;
    rnd = MT.mtRandomRealOpen();
  };

  uint64   pos;    //  Position in the seqLengths vector.
  double   rnd;    //  A random number.
};


bool seqOrderRandom(const seqEntry &a, const seqEntry &b) {
  return(a.rnd < b.rnd);
}

bool seqOrderNormal(const seqEntry &a, const seqEntry &b) {
  return(a.pos < b.pos);
}



void
doSample_sample(sampleParameters &samPar,
                uint64            numSeqsTotal,
                uint64            numBasesTotal,
                vector<uint64>   &seqLengths,
                vector<seqEntry> &seqOrder) {

  //  Randomize the sequences.

  sort(seqOrder.begin(), seqOrder.end(), seqOrderRandom);

  //  Do some math to figure out what sequences to report.

  if (samPar.desiredCoverage > 0.0) {
    samPar.desiredNumBases = (uint64)ceil(samPar.desiredCoverage * samPar.genomeSize);
  }

  if (samPar.desiredNumReads > 0) {
    fprintf(stderr, "Emitting " F_U64 " reads.\n",
            samPar.desiredNumReads);

    for (uint64 ii=0; ii<numSeqsTotal; ii++)
      if (ii < samPar.desiredNumReads)
        ;
      else
        seqLengths[seqOrder[ii].pos] = 0;
  }

  if (samPar.desiredNumBases > 0) {
    if (samPar.desiredCoverage > 0)
      fprintf(stderr, "Emitting %.3fx coverage; " F_U64 " bases.\n",
              samPar.desiredCoverage,
              samPar.desiredNumBases);
    else
      fprintf(stderr, "Emitting " F_U64 " bases.\n",
              samPar.desiredNumBases);

    for (uint64 nbe=0, ii=0; ii<numSeqsTotal; ii++) {
      if (nbe < samPar.desiredNumBases)
        nbe += seqLengths[seqOrder[ii].pos];
      else
        seqLengths[seqOrder[ii].pos] = 0;
    }
  }

  if (samPar.desiredFraction > 0.0) {
    fprintf(stderr, "Emitting %.4f fraction of the reads.\n",
            samPar.desiredFraction);

    for (uint64 ii=0; ii<numSeqsTotal; ii++) {
      if (seqOrder[ii].rnd < samPar.desiredFraction)
        ;
      else
        seqLengths[seqOrder[ii].pos] = 0;
    }
  }

  //  Unrandomize the sequences.  Not needed for 2-pass, but needed for 1-pass.

  sort(seqOrder.begin(), seqOrder.end(), seqOrderNormal);
}



void
doSample_paired(vector<char *> &inputs, sampleParameters &samPar) {

  samPar.initialize();

  vector<uint64>    numSeqsPerFile;
  vector<uint64>    seqLengths;
  vector<seqEntry>  seqOrder;

  uint64            numSeqsTotal  = 0;
  uint64            numBasesTotal = 0;

  vector<char *>    names;
  vector<char *>    sequences;
  vector<char *>    qualities;

  mtRandom          MT;

  //  Open output files. If paired, replace #'s in the output names with 1 or 2.

  FILE *outFile1 = NULL;
  FILE *outFile2 = NULL;

  {
    char  *a = strrchr(samPar.output1, '#');
    char  *b = strrchr(samPar.output2, '#');

    if (a == NULL)
      fprintf(stderr, "ERROR: Failed to find '#' in output name '%s'\n", samPar.output1), exit(1);
    if (b == NULL)
      fprintf(stderr, "ERROR: Failed to find '#' in output name '%s'\n", samPar.output2), exit(1);

    *a = '1';
    *b = '2';

    outFile1 = AS_UTL_openOutputFile(samPar.output1);
    outFile2 = AS_UTL_openOutputFile(samPar.output2);
  }

  //  Scan the inputs, saving the number of sequences in each and the length of each sequence.

  dnaSeq   seq1;
  dnaSeq   seq2;

  for (uint32 ff=0; ff<inputs.size(); ff += 2) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff+0]);
    dnaSeqFile  *sf2 = new dnaSeqFile(inputs[ff+1]);
    uint64       num = 0;

    bool   sf1more = sf1->loadSequence(seq1);
    bool   sf2more = sf2->loadSequence(seq2);

    while ((sf1more == true) &&
           (sf2more == true)) {
      seqLengths.push_back(seq1.length() + seq2.length());
      seqOrder.push_back(seqEntry(MT, numSeqsTotal));

      numSeqsTotal  += 1;
      numBasesTotal += seq1.length() + seq2.length();

      num += 1;

      sf1more = sf1->loadSequence(seq1);
      sf2more = sf2->loadSequence(seq2);
    }

    numSeqsPerFile.push_back(num);

    delete sf1;
    delete sf2;
  }

  //  Figure out what to output.

  doSample_sample(samPar, numSeqsTotal, numBasesTotal, seqLengths, seqOrder);

  //  Scan the inputs again, this time emitting sequences if their saved length isn't zero.

  for (uint32 ff=0; ff<inputs.size(); ff += 2) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff+0]);
    dnaSeqFile  *sf2 = new dnaSeqFile(inputs[ff+1]);
    uint64       num = 0;

    bool   sf1more = sf1->loadSequence(seq1);
    bool   sf2more = sf2->loadSequence(seq2);

    while ((sf1more == true) &&
           (sf2more == true)) {
      if (seqLengths[num] > 0) {
        AS_UTL_writeFastA(outFile1, seq1.bases(), seq1.length(), 0, ">%s\n", seq1.name());
        AS_UTL_writeFastA(outFile2, seq2.bases(), seq2.length(), 0, ">%s\n", seq2.name());
      }

      num += 1;

      sf1more = sf1->loadSequence(seq1);
      sf2more = sf2->loadSequence(seq2);
    }

    delete sf1;
    delete sf2;
  }

  AS_UTL_closeFile(outFile1, samPar.output1);
  AS_UTL_closeFile(outFile2, samPar.output2);
}



void
doSample_single(vector<char *> &inputs, sampleParameters &samPar) {

  samPar.initialize();

  vector<uint64>    numSeqsPerFile;
  vector<uint64>    seqLengths;
  vector<seqEntry>  seqOrder;

  uint64            numSeqsTotal  = 0;
  uint64            numBasesTotal = 0;

  vector<char *>    names;
  vector<char *>    sequences;
  vector<char *>    qualities;

  mtRandom          MT;

  //  Open output files. If paired, replace #'s in the output names with 1 or 2.

  FILE *outFile1 = AS_UTL_openOutputFile(samPar.output1);

  //  Scan the inputs, saving the number of sequences in each and the length of each sequence.

  dnaSeq   seq1;

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff]);
    uint64       num = 0;

    while (sf1->loadSequence(seq1)) {
      seqLengths.push_back(seq1.length());
      seqOrder.push_back(seqEntry(MT, numSeqsTotal));

      numSeqsTotal  += 1;
      numBasesTotal += seq1.length();

      num += 1;
    }

    numSeqsPerFile.push_back(num);

    delete sf1;
  }

  //  Figure out what to output.

  doSample_sample(samPar, numSeqsTotal, numBasesTotal, seqLengths, seqOrder);

  //  Scan the inputs again, this time emitting sequences if their saved length isn't zero.

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff]);
    uint64       num = 0;

    while (sf1->loadSequence(seq1)) {
      if (seqLengths[num] > 0)
        AS_UTL_writeFastA(outFile1, seq1.bases(), seq1.length(), 0, ">%s\n", seq1.name());

      num += 1;
    }

    delete sf1;
  }

  AS_UTL_closeFile(outFile1, samPar.output1);
}



void
doSample(vector<char *> &inputs, sampleParameters &samPar) {

  if (samPar.isPaired == false)
    doSample_single(inputs, samPar);
  else
    doSample_paired(inputs, samPar);
}

