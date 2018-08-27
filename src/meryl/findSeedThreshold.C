
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
//  profile[0] == frequency of kmer at positions[bgn..bgn+K]
//  profile[x] == 0 if the kmer at that position is invalid
//
//  profile is NOT resized.
//

class kmerProfile {
public:
  kmerProfile() {
  };

  uint32   position;
  uint32   value;
};

class
sortByPosition {
public:
  bool operator()(const kmerProfile &a, const kmerProfile &b) {
    if (a.position < b.position)  return(true);
    if (a.position > b.position)  return(false);

    return(a.value < b.value);
  };
};

class
sortByValue {
public:
  bool operator()(const kmerProfile &a, const kmerProfile &b) {
    if (a.value < b.value)  return(true);
    if (a.value > b.value)  return(false);

    return(a.position < b.position);
  };
};




void
findKmerProfile(kmerProfile *profile, char *seq, uint32 bgn, uint32 end, kmerCountExactLookup *merylLookup) {
  kmer     fmer;
  kmer     rmer;

  uint32   kmerLoad  = 0;
  uint32   kmerValid = fmer.merSize() - 1;

  for (uint32 ss=bgn+kmerValid, pp=0; ss<=end; ss++, pp++) {
    profile[pp].position = ss;
    profile[pp].value    = 0;
  }

  for (uint32 ss=bgn, pp=0; ss<end; ss++, pp++) {
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

    uint64  fval = merylLookup->value(fmer);
    uint64  rval = merylLookup->value(rmer);

#if 0
    if ((fval > 0) ||
        (rval > 0)) {
      char  str1[128], str2[128];

      fprintf(stderr, "pos %6u %s %9lu -  %s %9lu\n",
              ss,
              fmer.toString(str1), fval,
              rmer.toString(str2), rval);
    }
#endif

    if (fmer != rmer) {
      if (fmer < rmer)
        assert(rval == 0);
      else
        assert(fval == 0);
    }

    if (fmer < rmer)
      profile[pp].value = fval;
    else
      profile[pp].value = rval;
  }
}



//  For each read bgnID..endID, creates a file with the kmer frequency profile.
//  Files called 'outputPrefix-######.dat', based on readID.
//  Lines are 'position \t frequency', frequency is the count of the kmer that begins at position.
//
void
dumpProfile(kmerCountExactLookup  *merylLookup,
            sqStore               *seqStore,
            uint32                 bgnID,
            uint32                 endID,
            char                  *outputPrefix) {
  char          profileName[FILENAME_MAX+1];
  kmerProfile  *profileData = new kmerProfile [AS_MAX_READLEN];
  sqReadData   *readData    = new sqReadData;

  for (uint32 rr=bgnID; rr <= endID; rr++) {
    sqRead  *read  = seqStore->sqStore_getRead(rr);
    uint32   seqLen = read->sqRead_sequenceLength();
    FILE    *F;

    if (seqLen == 0)
      continue;

    seqStore->sqStore_loadReadData(read, readData);


    findKmerProfile(profileData,
                    readData->sqReadData_getSequence(),
                    0, seqLen,
                    merylLookup);

    {
      snprintf(profileName, FILENAME_MAX, "%s-%06u.dat", outputPrefix, rr);

      F = AS_UTL_openOutputFile(profileName);

      for (uint32 ii=0; ii<seqLen - kmerTiny::merSize(); ii++)
        fprintf(F, "%u\t%u\n", profileData[ii].position, profileData[ii].value);

      AS_UTL_closeFile(F, profileName);
    }

    {
      snprintf(profileName, FILENAME_MAX, "%s-%06u.gp", outputPrefix, rr);

      F = AS_UTL_openOutputFile(profileName);

      fprintf(F, "set logscale y\n");
      fprintf(F, "set terminal png size 1280,800\n");
      fprintf(F, "set output '%s-%06u.png'\n", outputPrefix, rr);
      fprintf(F, "plot '%s-%06u.dat' using 1:2 with lines title 'read %u'\n", outputPrefix, rr, rr);

      AS_UTL_closeFile(F, profileName);

      snprintf(profileName, FILENAME_MAX, "gnuplot %s-%06u.gp > /dev/null 2>&1", outputPrefix, rr);

      system(profileName);
    }
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

  //  Generate a histogram of read lengths.

  uint32   *lengthHistogram = new uint32 [maxReadLen + 1];
  uint32    nValidReads     = 0;

  for (uint32 ii=0; ii<maxReadLen+1; ii++)
    lengthHistogram[ii] = 0;

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

class thrReadData {
public:
  thrReadData() {
    threshold5 = 500000;
    repeatLen5 = 0;

    threshold3 = 500000;
    repeatLen3 = 0;
  };

  uint32   threshold5;
  uint32   repeatLen5;

  uint32   threshold3;
  uint32   repeatLen3;
};



class pickThreshold_ThreadData {
public:
  pickThreshold_ThreadData() {
    readData    = new sqReadData;

    kmerProfile5 = new kmerProfile [AS_MAX_READLEN];
    kmerProfile3 = new kmerProfile [AS_MAX_READLEN];
  };
  ~pickThreshold_ThreadData() {
    delete    readData;
    delete [] kmerProfile5;
    delete [] kmerProfile3;
  };

  sqReadData   *readData;
  kmerProfile  *kmerProfile5;
  kmerProfile  *kmerProfile3;
};



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

  uint32  MIN_COUNT_TRUSTED = 4;

  double  MIN_FRACTION      = 0.25;
  uint32  MIN_NUMBER        = 1000;

  uint32  maxValue          = 500000;
  uint32  cutoff;

  pickThreshold_ThreadData  *threads = new pickThreshold_ThreadData [omp_get_max_threads()];

  thrReadData  *readProfile  = new thrReadData [endID - bgnID + 1];

  //
  //  For each read, find the threshold where there are at least
  //  MIN_NUMBER of MIN_FRACTION*seqLen kmers at or below the threshold.
  //

#pragma omp parallel for schedule(dynamic, 10000)
  for (uint32 ii=bgnID; ii<endID; ii++) {
    pickThreshold_ThreadData  *thread = threads + omp_get_thread_num();

    sqRead  *read  = seqStore->sqStore_getRead(ii);

    seqStore->sqStore_loadReadData(read, thread->readData);

    char    *seq    = thread->readData->sqReadData_getSequence();
    uint32   seqLen = read->sqRead_sequenceLength();

    if (seqLen == 0)          {  nReadsZero++;   continue;  }
    if (seqLen  < olapLength) {  nReadsShort++;  continue;  }

    {
      findKmerProfile(thread->kmerProfile5, seq, 0, olapLength, merylLookup);               //  5' end
      sort(thread->kmerProfile5, thread->kmerProfile5 + olapLength, sortByValue());

      uint32 cutoff   = maxValue;
      uint32 position = olapLength;  //  pos of first kmer below threshold

      for (uint32 pp=0, nn=0; pp<olapLength; pp++) {
        if (thread->kmerProfile5[pp].value < MIN_COUNT_TRUSTED)
          continue;
        nn++;
        if ((nn < MIN_NUMBER) && (nn < MIN_FRACTION * seqLen))
          continue;
        if (thread->kmerProfile5[pp].value < cutoff) {
          cutoff   =     thread->kmerProfile5[pp].value;
          position = min(thread->kmerProfile5[pp].position, position);
        }
        break;
      }

      readProfile[ii-bgnID].threshold5 = cutoff;
      readProfile[ii-bgnID].repeatLen5 = position;
    }

    {
      findKmerProfile(thread->kmerProfile3, seq, seqLen-olapLength, seqLen, merylLookup);   //  3' end
      sort(thread->kmerProfile3, thread->kmerProfile3 + olapLength, sortByValue());

      uint32 cutoff   = maxValue;
      uint32 position = seqLen - olapLength;  //  pos of last kmer below threshold

      for (uint32 pp=0, nn=0; pp<olapLength; pp++) {
        if (thread->kmerProfile3[pp].value < MIN_COUNT_TRUSTED)
          continue;
        nn++;
        if ((nn < MIN_NUMBER) && (nn < MIN_FRACTION * seqLen))
          continue;
        if (thread->kmerProfile5[pp].value < cutoff) {
          cutoff   =     thread->kmerProfile5[pp].value;
          position = max(thread->kmerProfile5[pp].position, position);
        }
        break;
      }

      readProfile[ii-bgnID].threshold3 = cutoff;
      readProfile[ii-bgnID].repeatLen3 = position;
    }

    //fprintf(stderr, "readID %9u 5'Thresh %6u 5'Score %6u  3'Thresh %6u 3'Score %6u\n",
    //        ii,
    //        readProfile[ii-bgnID].threshold5, readProfile[ii-bgnID].repeatLen5,
    //        readProfile[ii-bgnID].threshold3, readProfile[ii-bgnID].repeatLen3);
  }

  //
  //  Now, with those parameters computed, score each threshold with:
  //    the number of overlaps missed (threshold below read threshold)
  //    the number of overlaps found  (threshold above read threshold)
  //    the number of overlaps extra  (threshold above read threshold)
  //

  for (uint32 threshold=2; threshold<1000000; threshold++) {
    uint64  olapsMissed = 0;    //  Potential overlaps missed; threshold too low.
    uint64  olapsFound  = 0;    //  Good!  Threshold found overlaps.
    uint64  olapsExtra  = 0;    //  Bad!   Threshold too high, found repeats.

    //  Over all reads

    for (uint32 ii=0; ii<endID-bgnID+1; ii++) {

      //  Analyze the 5' end.

      if      (threshold < readProfile[ii].threshold5) {
        olapsMissed += readProfile[ii].threshold5;
      }
      else {
        olapsFound  +=             readProfile[ii].threshold5;
        olapsExtra  += threshold - readProfile[ii].threshold5;
      }

      //  Analyze the 3' end.

      if      (threshold < readProfile[ii].threshold3) {
        olapsMissed += readProfile[ii].threshold3;
      }
      else {
        olapsFound  +=             readProfile[ii].threshold3;
        olapsExtra  += threshold - readProfile[ii].threshold3;
      }
    }
  
    //  Report for this threshold

    fprintf(stdout, "threshold %u missed %lu found %lu extra %lu\n",
            threshold, olapsMissed, olapsFound, olapsExtra);

    if (olapsMissed == 0)
      break;
  }



  //  On output, we can invert the number so it now means 'the fraction of reads number of reads that definitely cannot have an overlap'
  //  at some kmer threshold X.

#if 0
  FILE *F = fopen("profile.out", "w");
  for (uint32 xx=0; xx<maxValue; xx++)
    fprintf(F, "%u\t%7.3f\t%u\n", xx, (double)missedOverlaps[xx] / nReads, missedOverlaps[xx]);
  fclose(F);


  fprintf(stdout, "Total reads:   %u\n", seqStore->sqStore_getNumReads());
  fprintf(stderr, "Reads < %ubp %u\n",   olapLength, nReadsShort);
  fprintf(stderr, "Reads          %u\n", nReads);
  fprintf(stderr, "Totally Lost   %f\n", (double)missedOverlaps[maxValue-1] / nReads);
#endif
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

  uint32   numThreads       = 1;


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


    else if   (strcmp(argv[arg], "-threads") == 0) {
      numThreads = strtouint32(argv[++arg]);
      omp_set_num_threads(numThreads);
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

  kmerCountFileReader  *merylReader  = new kmerCountFileReader(merylPath, true);
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
