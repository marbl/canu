#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "meryl.H"
#include "mcBucket.H"
#include "mcDescription.H"
#include "outputMer.H"



///////////////////////////////////////
//
//  Implements binary operations
//    sub
//    abs
//
void
unaryOperations(char      personality,
                char    **mergeFiles,
                u32bit    mergeFilesLen,
                u32bit    desiredCount,
                char     *outputFile,
                bool      beVerbose) {

  if (mergeFilesLen != 1) {
    fprintf(stderr, "ERROR - must have exactly one file!\n");
    exit(1);
  }
  if (outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }
  if ((personality != PERSONALITY_LEQ) &&
      (personality != PERSONALITY_GEQ) &&
      (personality != PERSONALITY_EQ)) {
    fprintf(stderr, "ERROR - only personalities lessthan, lessthanorequal,\n");
    fprintf(stderr, "ERROR - greaterthan, greaterthanorequal, and equal\n");
    fprintf(stderr, "ERROR - are supported in unaryOperations().\n");
    fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
    exit(1);
  }

  //  Open all the input files.
  //
  char *inpath = new char [strlen(mergeFiles[0]) + 17];

  sprintf(inpath, "%s.mcidx", mergeFiles[0]);
  bitPackedFileReader   *IDX = new bitPackedFileReader(inpath);

  sprintf(inpath, "%s.mcdat", mergeFiles[0]);
  bitPackedFileReader   *DAT = new bitPackedFileReader(inpath);

  delete [] inpath;


  //  Open the output file
  //
  char *outpath = new char [strlen(outputFile) + 17];

  sprintf(outpath, "%s.mcidx", outputFile);
  bitPackedFileWriter   *oIDX = new bitPackedFileWriter(outpath);

  sprintf(outpath, "%s.mcdat", outputFile);
  bitPackedFileWriter   *oDAT = new bitPackedFileWriter(outpath);

  delete [] outpath;


  //  Read the parameters for the input file, and write them into the output file.
  //
  mcDescription  mcd;
  mcd.read(DAT);
  mcd.write(oDAT);


  //
  //  Read buckets from each file, merging them into the output
  //

  //  Create buckets
  //
  mcBucket *B = new mcBucket(IDX, DAT, &mcd);

  u32bit   itemsWritten    =  u32bitZERO;
  u64bit   maxBucket       = mcd._tableSizeInEntries;

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

    //  We'll count the number of things we have written.
    //
    itemsWritten    =  u32bitZERO;

    switch (personality) {
      case PERSONALITY_LEQ:
        for (u32bit pos=0; pos < B->_items; pos++) {
          if (B->_counts[pos] <= desiredCount) {
            outputMer(oDAT, mcd, b, B->_checks[pos], (u32bit)B->_counts[pos]);
            itemsWritten++;
          }
        }
        break;
      case PERSONALITY_GEQ:
        for (u32bit pos=0; pos < B->_items; pos++) {
          if (B->_counts[pos] >= desiredCount) {
            outputMer(oDAT, mcd, b, B->_checks[pos], (u32bit)B->_counts[pos]);
            itemsWritten++;
          }
        }
        break;
      case PERSONALITY_EQ:
        for (u32bit pos=0; pos < B->_items; pos++) {
          if (B->_counts[pos] == desiredCount) {
            outputMer(oDAT, mcd, b, B->_checks[pos], (u32bit)B->_counts[pos]);
            itemsWritten++;
          }
        }
        break;
    }

    //  write the number of entries
    //
    oIDX->putBits(itemsWritten, 32);

    //  read the next bucket
    //
    B->readBucket();
  }

  delete B;
  delete oIDX;
  delete oDAT;
  delete IDX;
  delete DAT;
}
