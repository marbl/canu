
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



bool
doSummarize_loadSequence(dnaSeqFile  *sf,
                         bool         asSequences,
                         char       *&name,   uint32    &nameMax,
                         char       *&seq,
                         uint8      *&qlt,    uint64    &seqMax,
                         uint64      &seqLen) {

  if (asSequences)
    return(sf->loadSequence(name, nameMax, seq, qlt, seqMax, seqLen));

  //  Otherwise, piece it together from multiple calls to get bases.
  //  loadBases() returns true if bases were loaded, and sets endOfSeq
  //  if the block returned is to the end of a sequence.

  uint64   bufferMax = 23;
  uint64   bufferLen = 0;
  char    *buffer    = new char [bufferMax];
  bool     endOfSeq  = false;

  resizeArray(name, 0, nameMax, (uint32)1024);
  resizeArrayPair(seq, qlt, 0, seqMax, seqLen+1, resizeArray_doNothing);

  name[0] = 0;
  seq[0]  = 0;
  qlt[0]  = 0;

  seqLen = 0;

  while (sf->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
    if (seqLen + bufferLen >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 2 * (seqLen + bufferLen + 1));

    assert(seqLen + bufferLen + 1 < seqMax);

    memcpy(seq + seqLen, buffer, sizeof(char) * bufferLen);
    seqLen += bufferLen;

    seq[seqLen] = 0;

    if (endOfSeq)
      break;
  }

  //  We get here in two ways:
  //    loadBases() immediately hit EOF, it will return false and set endOfSeq = false.
  //    Otherwise, it found a sequence and endOfSeq (here) _must_ be true.
  //  So, we can just return endOfSeq to mean 'there might be another sequence to load'.

  delete [] buffer;

  return(endOfSeq);
}



void
doSummarize_lengthHistogramSimple(vector<uint64> lengths) {
  uint64  bgn = 0;
  uint64  end = 0;

  sort(lengths.begin(), lengths.end(), less<uint64>());

  while (bgn < lengths.size()) {
    end = bgn + 1;

    while ((end < lengths.size()) &&
           (lengths[bgn] == lengths[end]))
      end++;

    fprintf(stdout, "%lu\t%lu\n", lengths[bgn], end - bgn);

    bgn = end;
  }
}



void
doSummarize_dumpLengths(vector<uint64> lengths) {

  sort(lengths.begin(), lengths.end(), less<uint64>());

  for (uint64 ii=0; ii<lengths.size(); ii++)
    fprintf(stdout, "%lu\n", lengths[ii]);
}



void
doSummarize_lengthHistogram(vector<uint64> lengths,
                            uint64         genomeSize,
                            bool           limitTo1x) {

  sort(lengths.begin(), lengths.end(), greater<uint64>());

  uint32   nLines   = 0;                      //  Number of lines in the NG table.

  uint32   nCols    = 63;                     //  Magic number to make the histogram the same width as the trinucleotide list
  uint32   nRows    = 0;                      //  Height of the histogram; dynamically set.
  uint32   nRowsMin = 50;                     //  Nothing really magic, just fits on the screen.

  uint64   lSum  = 0;                         //  Sum of the lengths we've encountered so far

  uint32   nStep = 10;                        //  Step of each N report.
  uint32   nVal  = nStep;                     //  Index of the threshold we're next printing.
  uint64   nThr  = genomeSize * nVal / 100;   //  Threshold lenth; if sum is bigger, emit and move to the next threshold

  //  Count the number of lines we expect to get in the NG table.

  for (uint32 ii=0; ii<lengths.size(); ii++) {
    lSum += lengths[ii];

    while (lSum >= nThr) {
      nLines++;

      if      (nVal <    200)  nVal += nStep;
      else if (nVal <   2000)  nVal += nStep * 10;
      else if (nVal <  20000)  nVal += nStep * 100;
      else if (nVal < 200000)  nVal += nStep * 1000;
      else                     nVal += nStep * 10000;

      nThr  = genomeSize * nVal / 100;
    }
  }

  if (lengths.size() == 0)
    return;

  uint64   minLength = lengths[lengths.size()-1];
  uint64   maxLength = lengths[0];

  if (nLines < nRowsMin)                                //  If there are too few lines in the NG table, make the
    nRows = nRowsMin;                                   //  histogram plot some minimal size, otherwise, make it
  else                                                  //  the same size as the NG table.
    nRows = nLines;

  double   bucketSized = (double)(maxLength - minLength) / nRows;
  uint32   bucketSize  = (uint32)ceil(bucketSized);

  if (bucketSize == 0)                                  //  Happens when all data is the same length.
    bucketSize = 1;

  nRows = (maxLength - minLength) / bucketSize;         //  With new bucketSize set, compute number of rows.

  if (nRows > nRowsMin) {                               //  But if we get WAY too many rows, reset.
    nRows       = nRowsMin;
    bucketSized = (double)(maxLength - minLength) / nRows;
    bucketSize  = (uint32)ceil(bucketSized);
  }

  lSum  = 0;                                            //  Reset for actually generating the length histogram.
  nStep = 10;
  nVal  = nStep;
  nThr  = genomeSize * nVal / 100;

  //  Generate the length histogram.

  uint32  *nSeqPerLen = new uint32 [nRows + 1];

  for (uint32 rr=0; rr<nRows+1; rr++)                   //  Clear the histogram.
    nSeqPerLen[rr] = 0;

  for (uint32 ii=0; ii<lengths.size(); ii++) {          //  Count number of sequences per size range.
    uint32 r = (lengths[ii] - minLength) / bucketSize;

    assert(r < nRows+1);
    nSeqPerLen[r]++;
  }

  uint64  maxCount = 1;                                 //  Avoids divide-by-zero, even if zero can never actually occur.

  for (uint32 rr=0; rr<nRows+1; rr++)                   //  Find the maximum number of sequences in a size range.
    if (maxCount < nSeqPerLen[rr])
      maxCount = nSeqPerLen[rr];

  char **histPlot = new char * [nRows + 1];

  for (uint32 rr=0; rr<nRows+1; rr++)                   //  28 is a magic number based on the histPlot[] format string below.
    histPlot[rr] = new char [28 + nCols + 1];           //  28 = 9 + 1 + 9 + 1 + 7 + 1

  for (uint32 rr=0; rr<nRows+1; rr++) {                 //  Generate histogram in text.
    uint32  nn = (uint32)ceil(nSeqPerLen[rr] * nCols / (double)maxCount);
    uint64  lo = (rr+0) * bucketSize + minLength;
    uint64  hi = (rr+1) * bucketSize + minLength - 1;

    if (lo == hi)
      sprintf(histPlot[rr], "%9" F_U64P "           %7" F_U32P "|",
              lo, nSeqPerLen[rr]);
    else
      sprintf(histPlot[rr], "%9" F_U64P "-%-9" F_U64P " %7" F_U32P "|",
              lo, hi, nSeqPerLen[rr]);

    for (uint32 cc=0; cc<nn; cc++)
      histPlot[rr][28 + cc] = '-';

    histPlot[rr][28 + nn]   = 0;
  }

  //  Output N table, with length histogram appended at the end of each line.

  uint32  hp = 0;

  fprintf(stdout, "\n");
  fprintf(stdout, "G=%-12" F_U64P "                     sum of  ||               length     num\n", genomeSize);
  fprintf(stdout,   "NG         length     index       lengths  ||                range    seqs\n");
  fprintf(stdout,   "----- ------------ --------- ------------  ||  ------------------- -------\n");

  //  Write lines if we're showing all data, or if we're below 1x coverage.

  for (uint32 ii=0; ii<lengths.size(); ii++) {
    lSum += lengths[ii];

    while (lSum >= nThr) {
      if ((limitTo1x == false) ||
          (nVal <= 100)) {
        if (hp <= nRows)
          fprintf(stdout, "%05"    F_U32P " %12" F_U64P " %9" F_U32P " %12" F_U64P "  ||  %s\n",
                  nVal, lengths[ii], ii, lSum,
                  histPlot[hp++]);
        else
          fprintf(stdout, "%05"    F_U32P " %12" F_U64P " %9" F_U32P " %12" F_U64P "  ||\n",
                  nVal, lengths[ii], ii, lSum);
      }

      if      (nVal <    200)   nVal += nStep;
      else if (nVal <   2000)   nVal += nStep * 10;
      else if (nVal <  20000)   nVal += nStep * 100;
      else if (nVal < 200000)   nVal += nStep * 1000;
      else                      nVal += nStep * 10000;

      nThr  = genomeSize * nVal / 100;
    }
  }

  //  If we're displaying exactly 1x, write empty lines to get to there.

  if (limitTo1x == true) {
    while (nVal <= 100) {
      if (hp <= nRows)
        fprintf(stdout, "%05"    F_U32P " %12s %9s %12s  ||  %s\n",
                nVal, "-", "-", "-",
                histPlot[hp++]);
      else
        fprintf(stdout, "%05"    F_U32P " %12s %9s %12s  ||\n",
                nVal, "-", "-", "-");

      nVal += nStep;
    }
  }

  //  Now the final summary line.

  if (genomeSize == 0)
    fprintf(stdout, "%07.3fx           %9" F_U64P " %12" F_U64P "  ||  %s\n", 0.0, lengths.size(), lSum, histPlot[hp++]);   //  Occurs if only empty sequences in the input!
  else if (hp <= nRows)
    fprintf(stdout, "%07.3fx           %9" F_U64P " %12" F_U64P "  ||  %s\n", (double)lSum / genomeSize, lengths.size(), lSum, histPlot[hp++]);
  else
    fprintf(stdout, "%07.3fx           %9" F_U64P " %12" F_U64P "  ||\n",     (double)lSum / genomeSize, lengths.size(), lSum);

  while (hp <= nRows)
    fprintf(stdout, "                                           ||  %s\n", histPlot[hp++]);

  fprintf(stdout, "\n");


  //  Cleanup.

  for (uint32 rr=0; rr<nRows+1; rr++)
    delete [] histPlot[rr];

  delete [] histPlot;
  delete [] nSeqPerLen;
}



void
doSummarize(vector<char *>       &inputs,
            summarizeParameters  &sumPar) {

  vector<uint64>  lengths;

  uint64          nSeqs  = 0;
  uint64          nBases = 0;

  uint32          mer = 0;

  uint64          mn[4]     = {0};
  uint64          dn[4*4]   = {0};
  uint64          tn[4*4*4] = {0};

  double          nmn = 0;
  double          ndn = 0;
  double          ntn = 0;

  uint32          nameMax = 0;
  char           *name    = NULL;
  uint64          seqMax  = 0;
  char           *seq     = NULL;
  uint8          *qlt     = NULL;
  uint64          seqLen  = 0;

  for (uint32 ff=0; ff<inputs.size(); ff++) {
    dnaSeqFile  *sf = new dnaSeqFile(inputs[ff]);

    //  If sequence,
    //    Count mono-, di- and tri-nucleotides.
    //    Count number of mono-, di- and tri-nucleotides.
    //    Count number of sequences and total bases.
    //    Save the lengths of sequences.

    while (doSummarize_loadSequence(sf, sumPar.asSequences, name, nameMax, seq, qlt, seqMax, seqLen) == true) {
      uint64  pos = 0;
      uint64  bgn = 0;

      if (pos < seqLen) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
      }

      if (pos < seqLen) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
        dn[mer & 0x0f]++;
      }

      while (pos < seqLen) {
        mer = ((mer << 2) | ((seq[pos++] >> 1) & 0x03)) & 0x3f;
        mn[mer & 0x03]++;
        dn[mer & 0x0f]++;
        tn[mer & 0x3f]++;
      }

      nmn +=                    (seqLen-0);
      ndn += (seqLen < 2) ? 0 : (seqLen-1);
      ntn += (seqLen < 3) ? 0 : (seqLen-2);

      //  If we're NOT splitting on N, add one sequence of the given length.

      if (sumPar.breakAtN == false) {
        nSeqs  += 1;
        nBases += seqLen;

        lengths.push_back(seqLen);
        continue;
      }

      //  But if we ARE splitting on N, add multiple sequences.

      pos = 0;
      bgn = 0;

      while (pos < seqLen) {

        //  Skip any N's.
        while ((pos < seqLen) && ((seq[pos] == 'n') ||
                                  (seq[pos] == 'N')))
          pos++;

        //  Remember our start position.
        bgn = pos;

        //  Move ahead until the end of sequence or an N.
        while ((pos < seqLen) && ((seq[pos] != 'n') &&
                                  (seq[pos] != 'N')))
          pos++;

        //  If a sequence, increment stuff.
        if (pos - bgn > 0) {
          nSeqs  += 1;
          nBases += pos - bgn;

          lengths.push_back(pos - bgn);
        }
      }
    }

    //  All done!

    delete sf;
  }

  delete [] name;
  delete [] seq;
  delete [] qlt;

  if (sumPar.genomeSize == 0)      //  If no genome size supplied, set it to the sum of lengths.
    sumPar.genomeSize = nBases;

  //  If only a simple histogram of lengths is requested, dump and done.

  if (sumPar.asSimple == true) {
    doSummarize_lengthHistogramSimple(lengths);
  }

  //  If only the read lengths are requested, dump and done.

  else if (sumPar.asLength == true) {
    doSummarize_dumpLengths(lengths);
  }

  //  Otherwise, generate a fancy histogram plot.
  //  And finish with the mono-, di- and tri-nucleotide frequencies.

#define FMT "%12" F_U64P " %6.4f"
#define GC "%05.02f%%"

  else {
    doSummarize_lengthHistogram(lengths, sumPar.genomeSize, sumPar.limitTo1x);

    if (nmn == 0)  nmn = 1;   //  Avoid divide by zero.
    if (ndn == 0)  ndn = 1;
    if (ntn == 0)  ntn = 1;

    double gc = 100.0 * (mn[0x01] + mn[0x03]) / (mn[0x00] + mn[0x01] + mn[0x03] + mn[0x02]);

    fprintf(stdout, "--------------------- --------------------- ----------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "       mononucleotide          dinucleotide                                                                                  trinucleotide\n");
    fprintf(stdout, "--------------------- --------------------- ----------------------------------------------------------------------------------------------\n");
    fprintf(stdout, ""             FMT " A" FMT " AA" FMT " AAA " FMT " AAC " FMT " AAG " FMT " AAT\n", mn[0x00], mn[0x00] / nmn, dn[0x00], dn[0x00] / ndn, tn[0x00], tn[0x00] / ntn, tn[0x01], tn[0x01] / ntn, tn[0x03], tn[0x03] / ntn, tn[0x02], tn[0x02] / ntn);
    fprintf(stdout, ""             FMT " C" FMT " AC" FMT " ACA " FMT " ACC " FMT " ACG " FMT " ACT\n", mn[0x01], mn[0x01] / nmn, dn[0x01], dn[0x01] / ndn, tn[0x04], tn[0x04] / ntn, tn[0x05], tn[0x05] / ntn, tn[0x07], tn[0x07] / ntn, tn[0x06], tn[0x06] / ntn);
    fprintf(stdout, ""             FMT " G" FMT " AG" FMT " AGA " FMT " AGC " FMT " AGG " FMT " AGT\n", mn[0x03], mn[0x03] / nmn, dn[0x03], dn[0x03] / ndn, tn[0x0c], tn[0x0c] / ntn, tn[0x0d], tn[0x0d] / ntn, tn[0x0f], tn[0x0f] / ntn, tn[0x0e], tn[0x0e] / ntn);
    fprintf(stdout, ""             FMT " T" FMT " AT" FMT " ATA " FMT " ATC " FMT " ATG " FMT " ATT\n", mn[0x02], mn[0x02] / nmn, dn[0x02], dn[0x02] / ndn, tn[0x08], tn[0x08] / ntn, tn[0x09], tn[0x09] / ntn, tn[0x0b], tn[0x0b] / ntn, tn[0x0a], tn[0x0a] / ntn);
    fprintf(stdout, "                     " FMT " CA" FMT " CAA " FMT " CAC " FMT " CAG " FMT " CAT\n",                           dn[0x04], dn[0x04] / ndn, tn[0x10], tn[0x10] / ntn, tn[0x11], tn[0x11] / ntn, tn[0x13], tn[0x13] / ntn, tn[0x12], tn[0x12] / ntn);
    fprintf(stdout, "      --GC--  --AT-- " FMT " CC" FMT " CCA " FMT " CCC " FMT " CCG " FMT " CCT\n",                           dn[0x05], dn[0x05] / ndn, tn[0x14], tn[0x14] / ntn, tn[0x15], tn[0x15] / ntn, tn[0x17], tn[0x17] / ntn, tn[0x16], tn[0x16] / ntn);
    fprintf(stdout, "      " GC "  " GC " " FMT " CG" FMT " CGA " FMT " CGC " FMT " CGG " FMT " CGT\n", gc, 100 - gc,             dn[0x07], dn[0x07] / ndn, tn[0x1c], tn[0x1c] / ntn, tn[0x1d], tn[0x1d] / ntn, tn[0x1f], tn[0x1f] / ntn, tn[0x1e], tn[0x1e] / ntn);
    fprintf(stdout, "                     " FMT " CT" FMT " CTA " FMT " CTC " FMT " CTG " FMT " CTT\n",                           dn[0x06], dn[0x06] / ndn, tn[0x18], tn[0x18] / ntn, tn[0x19], tn[0x19] / ntn, tn[0x1b], tn[0x1b] / ntn, tn[0x1a], tn[0x1a] / ntn);
    fprintf(stdout, "                     " FMT " GA" FMT " GAA " FMT " GAC " FMT " GAG " FMT " GAT\n",                           dn[0x0c], dn[0x0c] / ndn, tn[0x30], tn[0x30] / ntn, tn[0x31], tn[0x31] / ntn, tn[0x33], tn[0x33] / ntn, tn[0x32], tn[0x32] / ntn);
    fprintf(stdout, "                     " FMT " GC" FMT " GCA " FMT " GCC " FMT " GCG " FMT " GCT\n",                           dn[0x0d], dn[0x0d] / ndn, tn[0x34], tn[0x34] / ntn, tn[0x35], tn[0x35] / ntn, tn[0x37], tn[0x37] / ntn, tn[0x36], tn[0x36] / ntn);
    fprintf(stdout, "                     " FMT " GG" FMT " GGA " FMT " GGC " FMT " GGG " FMT " GGT\n",                           dn[0x0f], dn[0x0f] / ndn, tn[0x3c], tn[0x3c] / ntn, tn[0x3d], tn[0x3d] / ntn, tn[0x3f], tn[0x3f] / ntn, tn[0x3e], tn[0x3e] / ntn);
    fprintf(stdout, "                     " FMT " GT" FMT " GTA " FMT " GTC " FMT " GTG " FMT " GTT\n",                           dn[0x0e], dn[0x0e] / ndn, tn[0x38], tn[0x38] / ntn, tn[0x39], tn[0x39] / ntn, tn[0x3b], tn[0x3b] / ntn, tn[0x3a], tn[0x3a] / ntn);
    fprintf(stdout, "                     " FMT " TA" FMT " TAA " FMT " TAC " FMT " TAG " FMT " TAT\n",                           dn[0x08], dn[0x08] / ndn, tn[0x20], tn[0x20] / ntn, tn[0x21], tn[0x21] / ntn, tn[0x23], tn[0x23] / ntn, tn[0x22], tn[0x22] / ntn);
    fprintf(stdout, "                     " FMT " TC" FMT " TCA " FMT " TCC " FMT " TCG " FMT " TCT\n",                           dn[0x09], dn[0x09] / ndn, tn[0x24], tn[0x24] / ntn, tn[0x25], tn[0x25] / ntn, tn[0x27], tn[0x27] / ntn, tn[0x26], tn[0x26] / ntn);
    fprintf(stdout, "                     " FMT " TG" FMT " TGA " FMT " TGC " FMT " TGG " FMT " TGT\n",                           dn[0x0b], dn[0x0b] / ndn, tn[0x2c], tn[0x2c] / ntn, tn[0x2d], tn[0x2d] / ntn, tn[0x2f], tn[0x2f] / ntn, tn[0x2e], tn[0x2e] / ntn);
    fprintf(stdout, "                     " FMT " TT" FMT " TTA " FMT " TTC " FMT " TTG " FMT " TTT\n",                           dn[0x0a], dn[0x0a] / ndn, tn[0x28], tn[0x28] / ntn, tn[0x29], tn[0x29] / ntn, tn[0x2b], tn[0x2b] / ntn, tn[0x2a], tn[0x2a] / ntn);
    fprintf(stdout, "\n");
  }
}

