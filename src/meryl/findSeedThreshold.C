
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



//  For bases [bgn..end), return a list of the frequencies of the kmer at each position.
//
//  profile[0] == kmer at positions[bgn..bgn+K]
//  profile[x] == 0 if the kmer at that position is invalid
//
void
findKmerProfile(uint32 *profile, char *seq, uint32 bgn, uint32 end, kmerCountExactLookup *merylLookup) {
  kmer     fmer;
  kmer     rmer;

  uint32   kmerLoad  = 0;
  uint32   kmerValid = fmer.merSize() - 1;

  for (uint32 ss=bgn; ss<end; ss++) {

    profile[ss] = 0;

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

    profile[ss - kmerValid] = (fmer < rmer) ? merylLookup->value(fmer) : merylLookup->value(rmer);
  }
}



void
dumpProfile(kmerCountExactLookup  *merylLookup,
            sqStore               *seqStore,
            uint32                 bgnID,
            uint32                 endID,
            char                  *outputPrefix) {
  char          profileName[FILENAME_MAX+1];
  uint32       *profileData = new uint32 [AS_MAX_READLEN];
  sqReadData   *readData    = new sqReadData;

  for (uint32 rr=bgnID; rr <= endID; rr++) {
    sqRead  *read  = seqStore->sqStore_getRead(rr);
    uint32   seqLen = read->sqRead_sequenceLength();

    if (seqLen == 0)
      continue;

    seqStore->sqStore_loadReadData(read, readData);


    findKmerProfile(profileData,
                    readData->sqReadData_getSequence(),
                    0, seqLen,
                    merylLookup);

    snprintf(profileName, FILENAME_MAX, "%s-%06u.dat", outputPrefix, rr);

    FILE *F = AS_UTL_openOutputFile(profileName);

    for (uint32 ii=0; ii<seqLen; ii++)
      fprintf(F, "%u\t%u\n", ii, profileData[ii]);
      
    AS_UTL_closeFile(F, profileName);
  }

  delete    readData;
  delete [] profileData;
}



uint32
pickOverlapLength(sqStore   *seqStore,
                  double     lengthFraction) {

  uint32    maxReadLen = 0;

  //  Discover the maximum read length.

  for (uint32 ii=1; ii <= seqStore->sqStore_getNumReads(); ii++)
    maxReadLen = max(maxReadLen, seqStore->sqStore_getRead(ii)->sqRead_sequenceLength());

  //  Allocate space to count the lengths of each read.

  uint32   *lengthHistogram = new uint32 [maxReadLen + 1];
  uint32    nValidReads     = 0;

  for (uint32 ii=0; ii<maxReadLen+1; ii++)
    lengthHistogram[ii] = 0;

  //  Count the number of reads at each length.

  for (uint32 ii=1; ii <= seqStore->sqStore_getNumReads(); ii++) {
    sqRead   *read    = seqStore->sqStore_getRead(ii);
    uint32    readLen = read->sqRead_sequenceLength();

    if (readLen == 0)
      continue;

    lengthHistogram[readLen]++;
    nValidReads++;
  }

  //  Find the length such that 'lengthFraction' of the reads are shorter than it.

  uint32    nReadsTarget = lengthFraction * nValidReads;
  uint32    nReadsSum    = 0;
  uint32    olapLength   = 0;

  for (uint32 ii=1; ii <= maxReadLen; ii++) {
    if (nReadsSum + lengthHistogram[ii] >= nReadsTarget) {
      olapLength = ii - 1;
      break;
    }

    nReadsSum += lengthHistogram[ii];
  }

  //  Done.

  fprintf(stderr, "\n");
  fprintf(stderr, "Num valid reads:  %u\n",   nValidReads);
  fprintf(stderr, "Max read length:  %u\n",   maxReadLen);
  fprintf(stderr, "Length fraction:  %.2f\n", lengthFraction);
  fprintf(stderr, "  target nReads:  %u\n",   nReadsTarget);
  fprintf(stderr, "\n");
  fprintf(stderr, "Length at target: %u\n",   olapLength);
  fprintf(stderr, "\n");

  return(olapLength);
}



//  Build a profile.
//
//  If we were to ignore errors, then the smallest kmer in the list would
//  be the smallest threshold we could use to generate overlaps for this
//  read.  E.g., if there is a kmer with count 3, then setting the
//  threshold to 4 would mean that kmer survives to seed an overlap.
//
//  What we're really looking for is how much of the read is covered
//  by kmers of some value.  If, say, 10% of the read is covered by kmers
//  smaller than X then we declare that the overlap can be found with
//  a threshold of X+1.
//
//  With the sorted kmer values, we then need to find the 100th
//  smallest value, and declare the overlap detectable for values below that.

//  Find the 100th non-zero profile value.  This is the threshold at which we would
//  mask out 90% of the kmers, so declare that we WOULD NOT FIND the overlap for
//  thresholds below this, and WOULD FIND the overlap for thresholds above this.

void
pickThreshold(kmerCountExactLookup  *merylLookup,
              sqStore               *seqStore,
              uint32                 bgnID,
              uint32                 endID,
              uint32                 olapLength,
              char                  *outputPrefix) {

  uint32         nReadsShort   = 0;
  uint32         nReadsZero    = 0;
  uint32         nReadsTested  = 0;
  uint32         nReads        = 0;

  uint32  MIN_COUNT_TRUSTED = 2;

  double  MIN_FRACTION      = 0.10;
  uint32  MIN_NUMBER        = 100;

  uint32  maxValue          = 50000;


  sqReadData   *readData    = new sqReadData;

  uint32       *kmerProfile5 = new uint32 [AS_MAX_READLEN];
  uint32       *kmerProfile3 = new uint32 [AS_MAX_READLEN];

  uint32       *missedOverlaps = new uint32 [maxValue];
  memset(missedOverlaps, 0, sizeof(uint32) * maxValue);




  for (uint32 ii=bgnID; ii<endID; ii++) {
    sqRead  *read  = seqStore->sqStore_getRead(ii);

    seqStore->sqStore_loadReadData(read, readData);

    char    *seq    = readData->sqReadData_getSequence();
    uint32   seqLen = read->sqRead_sequenceLength();

    if (seqLen == 0) {
      nReadsZero++;
      continue;
    }

    if (seqLen < olapLength) {
      nReadsShort++;
      continue;
    }

    //
    //  5' end
    //

    findKmerProfile(kmerProfile5, seq, 0, olapLength, merylLookup);
    sort(kmerProfile5, kmerProfile5 + olapLength);

    uint32  cutoff = maxValue;

    for (uint32 pp=0, nn=0; pp<olapLength; pp++) {
      if (kmerProfile5[pp] < MIN_COUNT_TRUSTED)
        continue;
      nn++;
      if ((nn < MIN_NUMBER) && (nn < MIN_FRACTION * seqLen))
        continue;
      cutoff = min(kmerProfile5[pp], cutoff);
      break;
    }

    for (uint32 xx=0; xx < cutoff; xx++)
      missedOverlaps[xx]++;

    //
    //  3' end
    //

    findKmerProfile(kmerProfile3, seq, seqLen-olapLength, seqLen, merylLookup);
    sort(kmerProfile3, kmerProfile3 + olapLength);

    cutoff = maxValue;

    for (uint32 pp=0, nn=0; pp<olapLength; pp++) {
      if (kmerProfile3[pp] < MIN_COUNT_TRUSTED)
        continue;
      nn++;
      if ((nn < MIN_NUMBER) && (nn < MIN_FRACTION * seqLen))
        continue;
      cutoff = min(kmerProfile3[pp], cutoff);
      break;
    }

    for (uint32 xx=0; xx < cutoff; xx++)
      missedOverlaps[xx]++;

    //
    //  Spurs
    //



    nReads += 2;
  }

  //  On output, we can invert the number so it now means 'the fraction of reads number of reads that definitely cannot have an overlap'
  //  at some kmer threshold X.

  FILE *F = fopen("profile.out", "w");
  for (uint32 xx=0; xx<maxValue; xx++)
    fprintf(F, "%u\t%7.3f\t%u\n", xx, (double)missedOverlaps[xx] / nReads, missedOverlaps[xx]);
  fclose(F);


  fprintf(stdout, "Total reads:   %u\n", seqStore->sqStore_getNumReads());
  fprintf(stderr, "Reads < %ubp %u\n",   olapLength, nReadsShort);
  fprintf(stderr, "Reads          %u\n", nReads);
  fprintf(stderr, "Totally Lost   %f\n", (double)missedOverlaps[maxValue-1] / nReads);
}





int
main(int argc, char **argv) {
  char    *seqStorePath     = NULL;
  char    *merylPath        = NULL;

  bool     doDumpProfile    = false;
  bool     doPickThreshold  = true;

  char    *outputPrefix     = NULL;

  uint32   bgnID            = 1;
  uint32   endID            = UINT32_MAX;

  double   lengthFraction   = 0.5;


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

    else if   (strcmp(argv[arg], "-k") == 0) {
      kmer::setSize(strtouint32(argv[++arg]));
    }

    else if   (strcmp(argv[arg], "-profile") == 0) {
      doDumpProfile   = true;
      doPickThreshold = false;
    }

    else if   (strcmp(argv[arg], "-threshold") == 0) {
      doDumpProfile   = false;
      doPickThreshold = true;
    }

    else if   (strcmp(argv[arg], "-length") == 0) {
      lengthFraction = strtodouble(argv[++arg]);
    }

    else if   (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], bgnID, endID);
    }

    else if   (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];
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


  sqStore               *seqStore    = sqStore::sqStore_open(seqStorePath);

  kmerCountFileReader  *merylReader  = new kmerCountFileReader(merylPath, false, true);
  kmerCountExactLookup *merylLookup  = new kmerCountExactLookup(merylReader);

  bgnID = max(bgnID, (uint32)1);
  endID = min(endID, seqStore->sqStore_getNumReads() + 1);

  if (doDumpProfile)
    dumpProfile(merylLookup, seqStore, bgnID, endID, outputPrefix);

  if (doPickThreshold) {
    uint32 olapLength = pickOverlapLength(seqStore, lengthFraction);
    pickThreshold(merylLookup, seqStore, bgnID, endID, olapLength, outputPrefix);
  }

  delete merylLookup;
  delete merylReader;

  seqStore->sqStore_close();

  fprintf(stderr, "Bye.\n");
  exit(0);
}





#if 0
uint32   kmersTested = 0;
uint32   kmersFound  = 0;

kmer     fmer;
kmer     rmer;

uint32   kmerLoad  = 0;
uint32   kmerValid = fmer.merSize() - 1;

uint64   average = 0;
uint32   minV    = UINT32_MAX;
uint32   maxV    = 0;


fprintf(stdout, "      ----kmers---- -------value----------\n");
fprintf(stdout, "   ID  found tested    min    max  average\n");
fprintf(stdout, "----- ------ ------ ------ ------ --------\n");


#endif

#if 0
if (kmersTested > 0)
  fprintf(stdout, "%5u %6u %6u %6u %6u %8.2f\n",
          ii, kmersFound, kmersTested, minV, maxV, (double)average / kmersTested);
#endif
