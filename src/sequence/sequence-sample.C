
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
#include "mt19937ar.H"



class seqEntry {
public:
  seqEntry(mtRandom &MT, uint64 pos_, uint32 len_) {
    rnd = MT.mtRandomRealOpen();
    pos = pos_;
    len = len_;
    out = UINT32_MAX;
  };

  double   rnd;    //  A random number.
  uint64   pos;    //  Position in the seqLengths vector.
  uint32   len;    //  Length of this sequence.
  uint32   out;    //  Output file index for this sequence.
};


bool seqOrderRandom(const seqEntry &a, const seqEntry &b) {
  return(a.rnd < b.rnd);
}

bool seqOrderNormal(const seqEntry &a, const seqEntry &b) {
  return(a.pos < b.pos);
}


static
void
doSample_sample_logNumReads(sampleParameters &samPar) {
  if (samPar.numCopies == 1)
    fprintf(stderr, "Emitting " F_U64 " reads.\n",
            samPar.desiredNumReads);
  else
    fprintf(stderr, "Emitting " F_U64 " reads into each of %u copies.\n",
            samPar.desiredNumReads, samPar.numCopies);
}


static
void
doSample_sample_logCoverage(sampleParameters &samPar) {
  if (samPar.numCopies == 1)
    fprintf(stderr, "Emitting %.3fx coverage; " F_U64 " bases.\n",
            samPar.desiredCoverage, samPar.desiredNumBases);
  else
    fprintf(stderr, "Emitting %.3fx coverage; " F_U64 " bases into each of %u copies.\n",
            samPar.desiredCoverage, samPar.desiredNumBases, samPar.numCopies);
}


static
void
doSample_sample_logNumBases(sampleParameters &samPar) {
  if (samPar.numCopies == 1)
    fprintf(stderr, "Emitting " F_U64 " bases.\n",
            samPar.desiredNumBases);
  else
    fprintf(stderr, "Emitting " F_U64 " bases into each of %u copies.\n",
            samPar.desiredNumBases, samPar.numCopies);
}


static
void
doSample_sample_logFraction(sampleParameters &samPar) {
  if (samPar.numCopies == 1)
    fprintf(stderr, "Emitting %.4f fraction of the reads.\n",
            samPar.desiredFraction);
  else
    fprintf(stderr, "Emitting %.4f fraction of the reads into each of %u copies.\n",
            samPar.desiredFraction, samPar.numCopies);
}


void
doSample_sample(sampleParameters &samPar,
                uint64            numSeqsTotal,
                uint64            numBasesTotal,
                vector<seqEntry> &seqOrder) {

  //  Randomize the sequences.

  sort(seqOrder.begin(), seqOrder.end(), seqOrderRandom);

  //  Do some math to figure out what sequences to report.

  if (samPar.desiredCoverage > 0.0) {
    samPar.desiredNumBases = (uint64)ceil(samPar.desiredCoverage * samPar.genomeSize);
  }

  if (samPar.desiredFraction > 0.0) {
    samPar.desiredNumBases = (uint64)ceil(samPar.desiredFraction * numBasesTotal);
  }

  //  Scan the randomized reads, assigning each to an output file,
  //  and moving to the next file when the current one is too big.

  uint64  nr  = 0;
  uint64  nbe = 0;

  uint32  of = 0;

  if (samPar.desiredNumReads > 0) {
    doSample_sample_logNumReads(samPar);

    for (uint64 ii=0; ii<numSeqsTotal; ii++) {
      if (of < samPar.desiredNumReads) {
        nr++;
        seqOrder[ii].out = of;
      } else {
        nr = 0;
        seqOrder[ii].out = ++of;
      }
    }
  }

  if (samPar.desiredNumBases > 0) {
    if (samPar.desiredCoverage > 0)
      doSample_sample_logCoverage(samPar);
    else if (samPar.desiredFraction > 0)
      doSample_sample_logFraction(samPar);
    else
      doSample_sample_logNumBases(samPar);

    for (uint64 nbe=0, ii=0; ii<numSeqsTotal; ii++) {
      if (nbe < samPar.desiredNumBases) {
        nbe += seqOrder[ii].len;
        seqOrder[ii].out = of;
      } else {
        nbe  = seqOrder[ii].len;
        seqOrder[ii].out = ++of;
      }
    }
  }

  //  Unrandomize the sequences.

  sort(seqOrder.begin(), seqOrder.end(), seqOrderNormal);
}



void
doSample_paired(vector<char *> &inputs, sampleParameters &samPar) {

  samPar.initialize();

  vector<uint64>    numSeqsPerFile;
  vector<seqEntry>  seqOrder;

  uint64            numSeqsTotal  = 0;
  uint64            numBasesTotal = 0;

  vector<char *>    names;
  vector<char *>    sequences;
  vector<char *>    qualities;

  mtRandom          MT;

  //  No support for multiple copies.

  if (samPar.numCopies != 1)
    fprintf(stderr, "ERROR: No support for -copies in paried-end mode.\n"), exit(1);

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
      seqOrder.push_back(seqEntry(MT, numSeqsTotal, seq1.length() + seq2.length()));

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

  doSample_sample(samPar, numSeqsTotal, numBasesTotal, seqOrder);

  //  Scan the inputs again, this time emitting sequences if their saved length isn't zero.

  for (uint32 ff=0; ff<inputs.size(); ff += 2) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff+0]);
    dnaSeqFile  *sf2 = new dnaSeqFile(inputs[ff+1]);
    uint64       num = 0;

    bool   sf1more = sf1->loadSequence(seq1);
    bool   sf2more = sf2->loadSequence(seq2);

    while ((sf1more == true) &&
           (sf2more == true)) {
      uint32 of = seqOrder[num].out;

      if (of < samPar.numCopies) {
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



FILE *
doSample_single_openOutput(sampleParameters &samPar, uint32 ii) {
  uint32  ap = 0;

  while ((samPar.output1[ap] != 0) &&
         (samPar.output1[ap] != '#'))
    ap++;

  //  If no #'s and multiple copies requested, fail.  But if only one copy, just open
  //  the file and return.

  if ((samPar.output1[ap] == 0) && (samPar.numCopies > 1))
    fprintf(stderr, "ERROR: Failed to find '#' in output name '%s', and asked to make multiple copies.\n", samPar.output1), exit(1);

  if  (samPar.output1[ap] == 0)
    return(AS_UTL_openOutputFile(samPar.output1));

  //  We've got #'s in the string.  We want to replace the last block of 'em
  //  with digits (ap found above is the start of the first block, sigh).
  //  Search backwards.

  while (samPar.output1[ap] != 0)                   //  Find the end of the string.
    ap++;

  while ((ap > 0) && (samPar.output1[ap] != '#'))   //  Find the last #.
    ap--;

  while ((ap > 0) && (samPar.output1[ap] == '#'))   //  Find the start of that run.
    ap--;

  if (samPar.output1[ap] != '#')                    //  Handle the stupid edge case.
    ap++;

  //  We're guaranteed to have some #'s in the name.  Count 'em.

  uint32 dp = 0;

  while (samPar.output1[ap + dp] == '#')
    dp++;

  //  Make a copy of the name, insert the digits, and return the file.
  //
  //  Suppose we have three #'s in the string; dp will be 4.
  //  We'll copy in digs[7-3 = 4]; dp will be 2 after this
  //                digs[7-2 = 5]; dp will be 1 after this
  //                digs[7-1 = 6]; dp will be 0 after this
  //                digs[7-0 = 7]; the NUL byte is not copied.

  char  name[FILENAME_MAX+1] = {0};
  char  digs[8];

  snprintf(digs, 8, "%07u", ii);

  for (uint32 ii=0; samPar.output1[ii]; ii++) {
    if ((ii < ap) || (dp == 0))
      name[ii] = samPar.output1[ii];
    else
      name[ii] = digs[7 - dp--];
  }

  return(AS_UTL_openOutputFile(name));
}



void
doSample_single(vector<char *> &inputs, sampleParameters &samPar) {

  samPar.initialize();

  vector<uint64>    numSeqsPerFile;
  vector<seqEntry>  seqOrder;

  uint64            numSeqsTotal  = 0;
  uint64            numBasesTotal = 0;

  vector<char *>    names;
  vector<char *>    sequences;
  vector<char *>    qualities;

  mtRandom          MT;

  //  Open output files. If paired, replace #'s in the output names with 1 or 2.

  FILE **outFiles = new FILE * [samPar.numCopies];

  for (uint32 ii=0; ii<samPar.numCopies; ii++)
    outFiles[ii] = doSample_single_openOutput(samPar, ii);

  //  Scan the inputs, saving the number of sequences in each and the length of each sequence.

  dnaSeq   seq1;

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff]);
    uint64       num = 0;

    while (sf1->loadSequence(seq1)) {
      seqOrder.push_back(seqEntry(MT, numSeqsTotal, seq1.length()));

      numSeqsTotal  += 1;
      numBasesTotal += seq1.length();

      num += 1;
    }

    numSeqsPerFile.push_back(num);

    delete sf1;
  }

  //  Figure out what to output.

  doSample_sample(samPar, numSeqsTotal, numBasesTotal, seqOrder);

  //  Scan the inputs again, this time emitting sequences if their saved length isn't zero.

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf1 = new dnaSeqFile(inputs[ff]);
    uint64       num = 0;

    while (sf1->loadSequence(seq1)) {
      uint32 of = seqOrder[num].out;

      if (of < samPar.numCopies)
        AS_UTL_writeFastA(outFiles[of], seq1.bases(), seq1.length(), 0, ">%s\n", seq1.name());

      num += 1;
    }

    delete sf1;
  }

  for (uint32 ii=0; ii<samPar.numCopies; ii++)
    AS_UTL_closeFile(outFiles[ii], samPar.output1);
}



void
doSample(vector<char *> &inputs, sampleParameters &samPar) {

  if (samPar.isPaired == false)
    doSample_single(inputs, samPar);
  else
    doSample_paired(inputs, samPar);
}

