
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
 *    Brian P. Walenz beginning on 2018-JUL-23
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sqStore.H"

#include "files.H"
#include "strings.H"
#include "kmers.H"

int
main(int argc, char **argv) {
  char    *seqStorePath = NULL;
  char    *merylPath    = NULL;

  uint32   bgnID        = 1;
  uint32   endID        = UINT32_MAX;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStorePath = argv[++arg];
    }

    else if   (strcmp(argv[arg], "-M") == 0) {
      merylPath = argv[++arg];
    }

    else if   (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], bgnID, endID);
    }

    else {
    }

    arg++;
  }

  if (seqStorePath == NULL)    err.push_back("No sequence store (-S option) supplied.\n");
  if (merylPath    == NULL)    err.push_back("No kmer data (-M option) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqPath -M merylData ...\n", argv[0]);
    fprintf(stderr, "\n");
    exit(1);
  }


  sqStore               *seqStore   = sqStore::sqStore_open(seqStorePath);

  bgnID = max(bgnID, (uint32)1);
  endID = min(endID, seqStore->sqStore_getNumReads() + 1);

  sqReadData           *readData    = new sqReadData;

  kmerCountFileReader  *merylReader = new kmerCountFileReader(merylPath, false, true);
  kmerCountExactLookup *merylLookup = new kmerCountExactLookup(merylReader);

  fprintf(stdout, "      ----kmers---- -------value----------\n");
  fprintf(stdout, "   ID  found tested    min    max  average\n");
  fprintf(stdout, "----- ------ ------ ------ ------ --------\n");


  for (uint32 ii=bgnID; ii<endID; ii++) {
    sqRead  *read  = seqStore->sqStore_getRead(ii);

    seqStore->sqStore_loadReadData(read, readData);

    char    *seq    = readData->sqReadData_getSequence();
    uint32   seqLen = read->sqRead_sequenceLength();

    uint32   kmersTested = 0;
    uint32   kmersFound  = 0;

    kmer     fmer;
    kmer     rmer;

    uint32   kmerLoad  = 0;
    uint32   kmerValid = fmer.merSize() - 1;

    uint64   average = 0;
    uint32   minV    = UINT32_MAX;
    uint32   maxV    = 0;

    for (uint32 ss=0; ss<seqLen; ss++) {
      if ((seq[ss] != 'A') && (seq[ss] != 'a') &&   //  If not valid DNA, don't
          (seq[ss] != 'C') && (seq[ss] != 'c') &&   //  make a kmer, and reset
          (seq[ss] != 'G') && (seq[ss] != 'g') &&   //  the count until the next
          (seq[ss] != 'T') && (seq[ss] != 't')) {   //  valid kmer is available.
        kmerLoad = 0;
        continue;
      }

      fmer.addR(seq[ss]);
      rmer.addL(seq[ss]);

      if (kmerLoad < kmerValid) {   //  If not a full kmer, increase the length we've
        kmerLoad++;                 //  got loaded, and keep going.
        continue;
      }

      uint32  value;

      if (fmer < rmer)
        value = merylLookup->value(fmer);
      else
        value = merylLookup->value(rmer);

      kmersTested++;

      if (value) {
        minV = min(minV, value);
        maxV = max(maxV, value);

        kmersFound++;
      }

      average += value;
    }

    if (kmersTested > 0)
      fprintf(stdout, "%5u %6u %6u %6u %6u %8.2f\n",
              ii, kmersFound, kmersTested, minV, maxV, (double)average / kmersTested);
  }


  delete merylLookup;
  delete merylReader;

  seqStore->sqStore_close();

  exit(0);
}
