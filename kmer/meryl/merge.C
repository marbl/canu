#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "meryl.H"
#include "mcBucket.H"
#include "mcDescription.H"
#include "outputMer.H"


//  Tests if the newCheck is in the mask file.  Returns true/false,
//  and modified the bucket position.
//
bool
decideIfMasked(mcBucket  *mmB,
               u32bit&    mmbucketPosition,
               u64bit     newCheck) {
  bool isMasked = false;

  if (mmB) {
    while ((mmbucketPosition < mmB->_items) &&
           (mmB->_checks[mmbucketPosition] < newCheck))
      mmbucketPosition++;

    if ((mmbucketPosition >= mmB->_items) ||
        (mmB->_checks[mmbucketPosition] != newCheck))
      isMasked = true;
  }

  return(isMasked);
}


u32bit
multipleWrite(char                   personality,
              bitPackedFileWriter   *oDAT,
              mcDescription         *mcd,
              u64bit                 b,
              u32bit                 mergeFilesLen,
              u64bit                 currentCheck,
              u32bit                 currentCount,
              u32bit                 currentTimes,
              u32bit                 itemsWritten) {

  //  Reset the count to zero, unless the mer is in all the databases.
  //
  if ((personality == PERSONALITY_MIN) &&
      (currentTimes != mergeFilesLen))
    currentCount = u32bitZERO;

  if (currentCount == u32bitZERO)
    return(itemsWritten);

  switch (personality) {
    case PERSONALITY_MIN:
    case PERSONALITY_MINEXIST:
    case PERSONALITY_MAX:
    case PERSONALITY_ADD:
      outputMer(oDAT, mcd[0], b, currentCheck, currentCount);
      itemsWritten++;
      break;
    case PERSONALITY_AND:
      if (currentTimes == mergeFilesLen) {
        outputMer(oDAT, mcd[0], b, currentCheck, u32bitONE);
        itemsWritten++;
      }
      break;
    case PERSONALITY_NAND:
      if (currentTimes != mergeFilesLen) {
        outputMer(oDAT, mcd[0], b, currentCheck, u32bitONE);
        itemsWritten++;
      }
      break;
    case PERSONALITY_OR:
      outputMer(oDAT, mcd[0], b, currentCheck, u32bitONE);
      itemsWritten++;
      break;
    case PERSONALITY_XOR:
      if ((currentTimes % 2) == 1) {
        outputMer(oDAT, mcd[0], b, currentCheck, u32bitONE);
        itemsWritten++;
      }
      break;
    default:
      fprintf(stderr, "ERROR - invalid personality in multipleOperations::write\n");
      fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
      exit(1);
      break;
  }

  return(itemsWritten);
}



void
multipleOperations(char     personality,
                   char    **mergeFiles,
                   u32bit   mergeFilesLen,
		   char    *maskFile,
                   char    *outputFile,
                   bool     beVerbose) {

  if (mergeFilesLen < 2) {
    fprintf(stderr, "ERROR - must have at least two databases!\n");
    exit(1);
  }
  if (outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }
  if ((personality != PERSONALITY_MIN) &&
      (personality != PERSONALITY_MINEXIST) &&
      (personality != PERSONALITY_MAX) &&
      (personality != PERSONALITY_ADD) &&
      (personality != PERSONALITY_AND) &&
      (personality != PERSONALITY_NAND) &&
      (personality != PERSONALITY_OR) &&
      (personality != PERSONALITY_XOR)) {
    fprintf(stderr, "ERROR - only personalities min, minexist, max, add, and, nand, or, xor\n");
    fprintf(stderr, "ERROR - are supported in multipleOperations().\n");
    fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
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

  //  Open the mask, if it exists
  //
  bitPackedFileReader *mIDX = 0L;
  bitPackedFileReader *mDAT = 0L;
  if (maskFile) {
    char *maskpath = new char [strlen(maskFile) + 17];

    sprintf(maskpath, "%s.mcidx", maskFile);
    mIDX = new bitPackedFileReader(maskpath);

    sprintf(maskpath, "%s.mcdat", maskFile);
    mDAT = new bitPackedFileReader(maskpath);

    delete [] maskpath;
  }


  //  Open the output file
  //
  char *outpath = new char [strlen(outputFile) + 17];

  sprintf(outpath, "%s.mcidx", outputFile);
  bitPackedFileWriter   *oIDX = new bitPackedFileWriter(outpath);

  sprintf(outpath, "%s.mcdat", outputFile);
  bitPackedFileWriter   *oDAT = new bitPackedFileWriter(outpath);

  delete [] outpath;


  //  Read the parameters for each of the input (and mask) files.
  //  Check that they are compatable.
  //
  mcDescription *mcd  = new mcDescription [mergeFilesLen];
  mcDescription  mmcd;

  for (u32bit i=0; i<mergeFilesLen; i++)
    mcd[i].read(DAT[i]);

  bool    fail = false;

  if (mDAT) {
    mmcd.read(mDAT);

    fprintf(stderr, "Checking files AND mask.\n");
    for (u32bit i=0; i<mergeFilesLen; i++)
      fail |= checkSingleDescription(mcd+i, mergeFiles[i], &mmcd, maskFile);
  } else {
    fprintf(stderr, "Checking files.\n");
    for (u32bit i=1; i<mergeFilesLen; i++)
      fail |= checkSingleDescription(mcd+i-1, mergeFiles[i-1], mcd+i, mergeFiles[i]);
  }

  if (fail) {
    fprintf(stderr, "ERROR:  Files are not compatable.\n");
    exit(1);
  }

  mcDescription  omcd(mcd[0]);


  //  Determine the number of mers in all the input files, and
  //  the number of bits needed in the hash table for these
  //  mers.
  //
#if 0

  //  XXX:  WHY?!

  omcd._actualNumberOfMers = 0;
  for (u32bit i=0; i<mergeFilesLen; i++)
    omcd._actualNumberOfMers += mcd[i]._actualNumberOfMers;

  omcd._hashWidth  = 1;
  while ((omcd._actualNumberOfMers+1) > (u64bitONE << omcd._hashWidth))
    omcd._hashWidth++;


  //
  //  XXX: if -H is given on the command line, override the setting of
  //  hashWidth.
  //
  //  Also needed in binaryOp
  //
  omcd._hashWidth = 20;
#endif

  //  Write the description to the output.
  //
  omcd.write(oDAT);


  //
  //  Read buckets from each file, merging them into the output
  //


  //  Create buckets
  //
  mcBucket **B = new mcBucket* [mergeFilesLen];
  for (u32bit i=0; i<mergeFilesLen; i++)
    B[i] = new mcBucket(IDX[i], DAT[i], &mcd[i]);

  mcBucket  *mmB = 0L;
  if (mDAT)
    mmB = new mcBucket(mIDX, mDAT, &mmcd);

  u64bit   currentCheck     = ~u64bitZERO;
  u32bit   currentCount     =  u32bitZERO;
  u32bit   currentTimes     =  u32bitZERO;

  u32bit  *bucketPosition   = new u32bit [mergeFilesLen];
  u32bit   mmbucketPosition = u32bitZERO;

  u32bit   newPosition      =  u32bitZERO;
  u64bit   newCheck         = ~u64bitZERO;
  u32bit   newCount         =  u32bitZERO;

  u64bit   totalItems       =  u64bitZERO;
  u32bit   itemsWritten     =  u32bitZERO;

  u64bit   maxBucket        = mcd[0]._tableSizeInEntries;

  //  For each bucket, build and output the merged bucket.
  //
  for (u64bit b=0; b<maxBucket; b++) {

    if ((beVerbose) && ((b & 0xfff) == 0)) {
#ifdef TRUE64BIT
      fprintf(stderr, "Bucket 0x%016lx\r", b);
#else
      fprintf(stderr, "Bucket 0x%016llx\r", b);
#endif
      fflush(stderr);
    }

    //  Clear stuff
    //
    currentCheck     = ~u64bitZERO;
    currentCount     =  u32bitZERO;
    currentTimes     =  u32bitZERO;

    for(u32bit i=0; i<mergeFilesLen; i++)
      bucketPosition[i] = u32bitZERO;

    mmbucketPosition = u32bitZERO;

    newCheck         = ~u64bitZERO;
    newCount         =  u32bitZERO;
    newPosition      =  u32bitZERO;

    totalItems       =  u64bitZERO;
    itemsWritten     =  u32bitZERO;

    //  How many items?  If none, skip the merge and just write the index
    //  (otherwise, we would always be writing a count)
    //
    for (u32bit i=0; i<mergeFilesLen; i++)
      totalItems += B[i]->_items;

    if (totalItems > 0) {
      bool isFirstEntryInBucket = true;

      currentCheck = ~u64bitZERO;
      currentCount =  u32bitZERO;
      currentTimes =  u32bitZERO;

      //  How many buckets still have stuff in them?
      //
      u32bit   notDone = u32bitZERO;
      for (u32bit i=0; i<mergeFilesLen; i++)
        if (bucketPosition[i] < B[i]->_items)
          notDone++;

      //  Now, until the buckets are exhausted, merge them.
      //
      while (notDone > 0) {

        //  Find the smallest check value
        //
        newCheck    = ~u64bitZERO;
        newCount    =  u32bitZERO;
        newPosition =  u32bitZERO;
        for (u32bit i=0; i<mergeFilesLen; i++) {
          if ((bucketPosition[i] < B[i]->_items) &&
              (newCheck >= B[i]->_checks[bucketPosition[i]])) {
            newCheck    = B[i]->_checks[bucketPosition[i]];
            newCount    = (u32bit)B[i]->_counts[bucketPosition[i]];
            newPosition = i;
          }
        }
        bucketPosition[newPosition]++;


        ///////////////////////////////////////
        //
        //  Output and reset
        //
        //  If this is a new check, write the old one to disk, and reset
        //  the currentCheck and currentCount.
        //
        if ((isFirstEntryInBucket == false) &&
            (decideIfMasked(mmB, mmbucketPosition, currentCheck) == false) &&
            (newCheck != currentCheck)) {
          itemsWritten = multipleWrite(personality,
                                       oDAT,
                                       mcd,
                                       b,
                                       mergeFilesLen,
                                       currentCheck, currentCount, currentTimes, itemsWritten);

          currentCount = u32bitZERO;
          currentTimes = u32bitZERO;
        }

        //  currentCheck is either the same as newCheck, or unset
        //
        currentCheck = newCheck;

        //fprintf(stderr, "current: 0x%016lx %u -- new 0x%016lx %u\n", currentCheck, currentCount, newCheck, newCount);


        //  Now we're not on the first entry in the bucket!
        //
        isFirstEntryInBucket = false;

        ///////////////////////////////////////
        //
        //  Perform the operation
        //
        switch (personality) {
          case PERSONALITY_MIN:
          case PERSONALITY_MINEXIST:
            if (currentTimes == 0) {
              currentCount = newCount;
            } else {
              if (currentCount > newCount)
                currentCount = newCount;
            }
            break;
          case PERSONALITY_MAX:
            if (currentCount < newCount)
              currentCount = newCount;
            break;
          case PERSONALITY_ADD:
            //fprintf(stderr, "ADD: %u += %u\n", currentCount, newCount);
            currentCount += newCount;
            break;
          case PERSONALITY_AND:
          case PERSONALITY_NAND:
          case PERSONALITY_OR:
          case PERSONALITY_XOR:
            currentCount = 1;
            break;
          default:
            fprintf(stderr, "ERROR - invalid personality in multipleOperations::operate\n");
            fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
            exit(1);
            break;
        }

        currentTimes++;

        //  Decrement the semaphore if the bucket is done.  This eliminates
        //  the bucket from the smallest check value code above.
        //
        if (bucketPosition[newPosition] >= B[newPosition]->_items)
          notDone--;
      }

      //  write the final entry
      //
      if (decideIfMasked(mmB, mmbucketPosition, currentCheck) == false)
        itemsWritten = multipleWrite(personality,
                                     oDAT,
                                     mcd,
                                     b,
                                     mergeFilesLen,
                                     currentCheck, currentCount, currentTimes, itemsWritten);
    }

    //  write the number of entries
    //
    oIDX->putBits(itemsWritten, 32);

    //  read the next set of buckets
    //
    for (u32bit i=0; i<mergeFilesLen; i++)
      B[i]->readBucket();

    if (mmB)
      mmB->readBucket();
  }

  delete [] bucketPosition;

  for (u32bit i=0; i<mergeFilesLen; i++)
    delete B[i];

  delete [] B;

  delete [] mcd;

  delete oIDX;
  delete oDAT;

  for (u32bit i=0; i<mergeFilesLen; i++) {
    delete IDX[i];
    delete DAT[i];
  }
  delete [] IDX;
  delete [] DAT;
}
