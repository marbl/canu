#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libbri.H"
#include "mcBucket.H"
#include "mcDescription.H"

static mcDescription mcd;

void
dump(char *inputFile) {

  if (inputFile == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  //  Open the counted sequence files
  //
  char *inpath = new char [strlen(inputFile) + 17];

  sprintf(inpath, "%s.mcidx", inputFile);
  bitPackedFileReader *IDX = new bitPackedFileReader(inpath);

  sprintf(inpath, "%s.mcdat", inputFile);
  bitPackedFileReader *DAT = new bitPackedFileReader(inpath);

  delete [] inpath;

  //
  //  Read the parameters
  //
  mcd.read(DAT);
  mcd.print(stderr);

  //
  //  Determine the count for each mer
  //
  mcBucket *B = new mcBucket(IDX, DAT, &mcd);

  for (u32bit i=0; i<mcd._tableSizeInEntries; i++)
    B->dump(stdout);

  delete DAT;
  delete IDX;
}

void
dumpThreshold(char *inputFile, u32bit threshold) {

  if (inputFile == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  //  Open the counted sequence files
  //
  char *inpath = new char [strlen(inputFile) + 17];

  sprintf(inpath, "%s.mcidx", inputFile);
  bitPackedFileReader *IDX = new bitPackedFileReader(inpath);

  sprintf(inpath, "%s.mcdat", inputFile);
  bitPackedFileReader *DAT = new bitPackedFileReader(inpath);

  delete [] inpath;

  //
  //  Read the parameters
  //
  mcd.read(DAT);
  mcd.print(stderr);

  //
  //  Determine the count for each mer
  //
  mcBucket *B = new mcBucket(IDX, DAT, &mcd);

  u64bit mer;
  char   ms[33];

  for (u32bit i=0; i<mcd._tableSizeInEntries; i++) {
    for (u32bit i=0; i<B->_items; i++) {
      if (B->_counts[i] >= threshold) {
        mer = B->_bucketID << B->_chckBits | B->_checks[i];

        for (u32bit z=0; z<mcd._merSizeInBases; z++)
          ms[mcd._merSizeInBases-z-1] = decompressSymbol[(mer >> (2*z)) & 0x03];
        ms[mcd._merSizeInBases] = 0;

#ifdef TRUE64BIT
        fprintf(stdout, ">%lu\n%s\n", B->_counts[i], ms);
#else
        fprintf(stdout, ">%llu\n%s\n", B->_counts[i], ms);
#endif
      }
    }

    B->readBucket();
  }

  delete DAT;
  delete IDX;
}

void
countUnique(char *inputFile) {

  if (inputFile == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  char *inpath = new char [strlen(inputFile) + 17];

  sprintf(inpath, "%s.mcidx", inputFile);
  bitPackedFileReader *IDX = new bitPackedFileReader(inpath);

  sprintf(inpath, "%s.mcdat", inputFile);
  bitPackedFileReader *DAT = new bitPackedFileReader(inpath);

  delete [] inpath;

  mcd.read(DAT);

  u64bit numMers     = 0;
  u64bit numUnique   = 0;
  u64bit numDistinct = 0;

  mcBucket *B = new mcBucket(IDX, DAT, &mcd);

  for (u32bit i=0; i<mcd._tableSizeInEntries; i++) {
    numDistinct += B->_items;
    for (u32bit j=0; j<B->_items; j++) {
      numMers += B->_counts[j];
      if (B->_counts[j] == 0)
        fprintf(stderr, "Got something with zero count?\n");
      if (B->_counts[j] == 1)
        numUnique++;
    }
    B->readBucket();
  }

  delete DAT;
  delete IDX;

#ifdef TRUE64BIT
  fprintf(stderr, "Found %lu mers.\n", numMers);
  fprintf(stderr, "Found %lu distinct mers.\n", numDistinct);
  fprintf(stderr, "Found %lu unique mers.\n", numUnique);
#else
  fprintf(stderr, "Found %llu mers.\n", numMers);
  fprintf(stderr, "Found %llu distinct mers.\n", numDistinct);
  fprintf(stderr, "Found %llu unique mers.\n", numUnique);
#endif
}

#define MAX_DISTANCE (16 * 1024 * 1024)

void
plotDistanceBetweenMers(char *inputFile) {
  int     distances[MAX_DISTANCE] = { 0 };

  if (inputFile == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  char *inpath = new char [strlen(inputFile) + 17];

  sprintf(inpath, "%s.mcidx", inputFile);
  bitPackedFileReader *IDX = new bitPackedFileReader(inpath);

  sprintf(inpath, "%s.mcdat", inputFile);
  bitPackedFileReader *DAT = new bitPackedFileReader(inpath);

  delete [] inpath;

  mcd.read(DAT);

  mcBucket *B = new mcBucket(IDX, DAT, &mcd);

  u64bit  lastMer = 0;
  u64bit  thisMer = 0;

  for (u32bit i=0; i<mcd._tableSizeInEntries; i++) {

    for (u32bit j=0; j<B->_items; j++) {

      thisMer = B->_bucketID << B->_chckBits | B->_checks[j];

      if ((thisMer - lastMer) < MAX_DISTANCE)
        distances[thisMer - lastMer]++;
      else
#ifdef TRUE64BIT
        fprintf(stderr, "Too large!  %lu\n", thisMer - lastMer);
#else
        fprintf(stderr, "Too large!  %llu\n", thisMer - lastMer);
#endif
      lastMer = thisMer;
    }

    B->readBucket();
  }

  delete DAT;
  delete IDX;

  for (u32bit i=0; i< MAX_DISTANCE; i++)
    printf("%u\n", distances[i]);

}





void
plotHistogram(char *inputFile, char *outputFile) {
  u32bit   *H = new u32bit [64 * 1024 * 1024];

  if (inputFile == 0L) {
    fprintf(stderr, "ERROR - no counted database file specified.\n");
    exit(1);
  }

  char *inpath = new char [strlen(inputFile) + 17];

  sprintf(inpath, "%s.mcidx", inputFile);
  bitPackedFileReader *IDX = new bitPackedFileReader(inpath);

  sprintf(inpath, "%s.mcdat", inputFile);
  bitPackedFileReader *DAT = new bitPackedFileReader(inpath);

  delete [] inpath;

  mcd.read(DAT);

  u64bit numMers     = 0;
  u64bit numUnique   = 0;
  u64bit numDistinct = 0;
  u64bit numHuge     = 0;

  mcBucket *B = new mcBucket(IDX, DAT, &mcd);

  for (u32bit i=0; i<16*1024*1024; i++)
    H[i] = 0;

  for (u32bit i=0; i<mcd._tableSizeInEntries; i++) {
    numDistinct += B->_items;
    for (u32bit j=0; j<B->_items; j++) {
      numMers += B->_counts[j];
      if (B->_counts[j] == 0)
        fprintf(stderr, "Got something with zero count?\n");
      if (B->_counts[j] == 1)
        numUnique++;
      if (B->_counts[j] < 16 * 1024 * 1024)
        H[B->_counts[j]]++;
      else
        numHuge++;
    }
    B->readBucket();
  }

  delete DAT;
  delete IDX;

#ifdef TRUE64BIT
  fprintf(stderr, "Found %lu mers.\n", numMers);
  fprintf(stderr, "Found %lu distinct mers.\n", numDistinct);
  fprintf(stderr, "Found %lu unique mers.\n", numUnique);
  fprintf(stderr, "Found %lu huge mers.\n", numHuge);
#else
  fprintf(stderr, "Found %llu mers.\n", numMers);
  fprintf(stderr, "Found %llu distinct mers.\n", numDistinct);
  fprintf(stderr, "Found %llu unique mers.\n", numUnique);
  fprintf(stderr, "Found %llu huge mers.\n", numHuge);
#endif

  FILE *F = fopen(outputFile, "w");
  for (u32bit i=0; i<16 * 1024 * 1024; i++)
#ifdef TRUE64BIT
    fprintf(F, "%u\n", H[i]);
#else
    fprintf(F, "%lu\n", H[i]);
#endif
  fclose(F);
}



#if 0
void
dumpStatsOnTwoFiles(char   **mergeFiles,
                    u32bit   mergeFilesLen) {
}
#endif




//  Dumps the shared mers between two files, such that the counts of
//  both are between lowCount and highCount
//
#if 0
void
dumpShared(char   **mergeFiles,
           u32bit   mergeFilesLen,
           u32bit   lowCount,
           u32bit   highCount) {
  if (mergeFilesLen != 2) {
    fprintf(stderr, "ERROR - I can only dumpShared on exactly two files!\n");
    exit(1);
  }

  //  Open all the input files.
  //
  bitPackedFileReader   **IDX = new bitPackedFileReader* [mergeFilesLen];
  bitPackedFileReader   **DAT = new bitPackedFileReader* [mergeFilesLen];

  for (u32bit i=0; i<mergeFilesLen; i++) {
    char *inpath = new char [strlen(mergeFiles[i]) + 17];

    sprintf(inpath, "%s.mcidx", mergeFiles[i]);
    IDX[i] = new bitPackedFileReader(inpath);

    sprintf(inpath, "%s.mcdat", mergeFiles[i]);
    DAT[i] = new bitPackedFileReader(inpath);

    delete [] inpath;
  }

  //  Read the parameters for each of the input files.  Check
  //  that the input files are compatable.
  //
  mcDescription *mcd  = new mcDescription [mergeFilesLen];
  for (u32bit i=0; i<mergeFilesLen; i++)
    mcd[i].read(DAT[i]);

  if (checkDescriptions(mcd, mergeFilesLen)) {
    fprintf(stderr, "ERROR:  Files are not compatable.\n");
    exit(1);
  }

  //  Create buckets
  //
  mcBucket **B = new mcBucket* [mergeFilesLen];

  for (u32bit i=0; i<mergeFilesLen; i++)
    B[i] = new mcBucket(IDX[i], DAT[i], &mcd[i]);




#if 0


  u64bit mer;
  char   ms[33];

  for (u32bit i=0; i<mcd[0]._tableSizeInEntries; i++) {

    //  The first bucket is read already

    bool notDone = true;
    while (notDone) {

      //  Find the first bucket with a value
      //
      u32bit firstBucket = 0;


      //  XXX:  Need to keep track of the position in each bucket!


      for (u32bit z=0; (z<mergeFilesLen) && (B; z++) {
        if (
      }
    }


    while (

    for (u32bit i=0; i<B->_items; i++) {
      if (B->_counts[i] >= threshold) {
        mer = B->_bucketID << B->_chckBits | B->_checks[i];

        for (u32bit z=0; z<mcd._merSizeInBases; z++)
          ms[mcd._merSizeInBases-z-1] = decompressSymbol[(mer >> (2*z)) & 0x03];
        ms[mcd._merSizeInBases] = 0;

#ifdef TRUE64BIT
        fprintf(stdout, ">%lu\n%s\n", B->_counts[i], ms);
#else
        fprintf(stdout, ">%llu\n%s\n", B->_counts[i], ms);
#endif
      }
    }

    B->readBucket();
  }



#endif

  for (u32bit i=0; i<mergeFilesLen; i++) {
    delete IDX[i];
    delete DAT[i];
  }
  delete [] IDX;
  delete [] DAT;
}

#endif
