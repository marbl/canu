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
binaryOperations(merylArgs *args) {

  if (args->mergeFilesLen != 2) {
    fprintf(stderr, "ERROR - must have exactly two files to use reduce!\n");
    exit(1);
  }
  if (args->outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }
  if ((args->personality != PERSONALITY_SUB) &&
      (args->personality != PERSONALITY_ABS)) {
    fprintf(stderr, "ERROR - only personalities sub and abs\n");
    fprintf(stderr, "ERROR - are supported in binaryOperations().\n");
    fprintf(stderr, "ERROR - this is a coding error, not a user error.\n");
    exit(1);
  }

  //  Open all the input files.
  //
  bitPackedFileReader   **IDX = new bitPackedFileReader* [args->mergeFilesLen];
  bitPackedFileReader   **DAT = new bitPackedFileReader* [args->mergeFilesLen];

  for (u32bit i=0; i<args->mergeFilesLen; i++) {
    char *inpath = new char [strlen(args->mergeFiles[i]) + 17];

    sprintf(inpath, "%s.mcidx", args->mergeFiles[i]);
    IDX[i] = new bitPackedFileReader(inpath);

    sprintf(inpath, "%s.mcdat", args->mergeFiles[i]);
    DAT[i] = new bitPackedFileReader(inpath);

    delete [] inpath;
  }


  //  Open the output file
  //
  char *outpath = new char [strlen(args->outputFile) + 17];

  sprintf(outpath, "%s.mcidx", args->outputFile);
  bitPackedFileWriter   *oIDX = new bitPackedFileWriter(outpath);

  sprintf(outpath, "%s.mcdat", args->outputFile);
  bitPackedFileWriter   *oDAT = new bitPackedFileWriter(outpath);

  delete [] outpath;


  //  Read the parameters for each of the input files.  Check
  //  that the input files are compatable.
  //
  mcDescription *mcd  = new mcDescription [args->mergeFilesLen];
  for (u32bit i=0; i<args->mergeFilesLen; i++)
    mcd[i].read(DAT[i]);

  if (checkSingleDescription(mcd+0, args->mergeFiles[0], mcd+1, args->mergeFiles[1])) {
    fprintf(stderr, "ERROR:  Files are not compatable.\n");
    exit(1);
  }

  mcDescription  omcd(mcd[0]);


  //  Determine the number of mers in all the input files, and
  //  the number of bits needed in the hash table for these
  //  mers.
  //
#if 0
  omcd._actualNumberOfMers = 0;
  for (u32bit i=0; i<args->mergeFilesLen; i++)
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
  mcBucket **B = new mcBucket* [args->mergeFilesLen];

  for (u32bit i=0; i<args->mergeFilesLen; i++)
    B[i] = new mcBucket(IDX[i], DAT[i], &mcd[i]);

  u32bit   itemsWritten    =  u32bitZERO;
  u64bit   maxBucket       = mcd[0]._tableSizeInEntries;

  //  For each bucket, build and output the merged bucket.
  //
  for (u64bit b=0; b<maxBucket; b++) {

    if ((args->beVerbose) && ((b & 0xfff) == 0)) {
      fprintf(stderr, "Bucket "u64bitHEX"\r", b);
      fflush(stderr);
    }

    //  We'll count the number of things we have written.
    //
    itemsWritten    =  u32bitZERO;

    //  If the second bucket isn't the current bucket, there are no counts
    //  to write.  We'll print a warning.
    //
    //  If both buckets are the current bucket, we have to subtract.
    //
    //  If only the first bucket is the current, we just dump it.
    //
    //  If neither bucket is the current, we do nothing until we get
    //  to the current.

    //  If the first bucketID is the current bucket, but the second isn't,
    //  there is nothing to subtract, and we just need to write out the
    //  first bucket.
    //
    if ((B[0]->_bucketID == b) && (B[1]->_bucketID != b)) {
      for (u32bit i=0; i<B[0]->_items; i++) {
        outputMer(oDAT, &mcd[0], b, B[0]->_checks[i], (u32bit)B[0]->_counts[i]);
        itemsWritten++;
      }
    }

    //  Scream and should if the first bucket is not current, but the
    //  second one is.
    //
    if ((B[0]->_bucketID != b) && (B[1]->_bucketID == b)) {
      fprintf(stderr, "WARNING: Negative count for bck="u64bitHEX"; "u64bitFMT" mers with zero output!\n",
              B[1]->_bucketID, B[1]->_items);
    }

    //  If the buckets are the same, we need to subtract the counts.
    //
    if ((B[0]->_bucketID == b) && (B[1]->_bucketID == b)) {

      //  The size of each bucket can be different, and we're probably
      //  not going to get the same mers in each.
      //
      u32bit pos1 = 0;
      u32bit pos2 = 0;

      while ((pos1 < B[0]->_items) && (pos2 < B[1]->_items)) {


          switch (args->personality) {
            case PERSONALITY_SUB:
              if        (B[0]->_checks[pos1] == B[1]->_checks[pos2]) {
                if (B[0]->_counts[pos1] < B[1]->_counts[pos2]) {
                  fprintf(stderr, "WARNING: Negative count for bck="u64bitHEX" chk="u64bitHEX", counts = "u64bitFMT" - "u64bitFMT"; zero output instead\n",
                          B[0]->_bucketID, B[0]->_checks[pos1],
                          B[0]->_counts[pos1], B[1]->_counts[pos2]);
                  outputMer(oDAT, &mcd[0], b, B[0]->_checks[pos1], 0);
                  itemsWritten++;
                }
                if (B[0]->_counts[pos1] > B[1]->_counts[pos2]) {
                  outputMer(oDAT, &mcd[0], b, B[0]->_checks[pos1], (u32bit)B[0]->_counts[pos1] - (u32bit)B[1]->_counts[pos2]);
                  itemsWritten++;
                }
                pos1++;
                pos2++;
              } else if (B[0]->_checks[pos1] < B[1]->_checks[pos2]) {
                outputMer(oDAT, &mcd[0], b, B[0]->_checks[pos1], (u32bit)B[0]->_counts[pos1]);
                itemsWritten++;
                pos1++;
              } else {
                fprintf(stderr, "WARNING: Negative count for bck="u64bitHEX" chk="u64bitHEX", counts = "u64bitFMT" - "u64bitFMT"; zero output instead\n",
                        B[0]->_bucketID, B[0]->_checks[pos1],
                        B[0]->_counts[pos1], B[1]->_counts[pos2]);
                outputMer(oDAT, &mcd[0], b, B[0]->_checks[pos1], u32bitZERO);
                itemsWritten++;
                pos2++;
              }
              break;
            case PERSONALITY_ABS:
              if        (B[0]->_checks[pos1] == B[1]->_checks[pos2]) {
                if (B[0]->_counts[pos1] < B[1]->_counts[pos2]) {
                  outputMer(oDAT, &mcd[0], b, B[0]->_checks[pos1], (u32bit)B[1]->_counts[pos2] - (u32bit)B[0]->_counts[pos1]);
                  itemsWritten++;
                }
                if (B[0]->_counts[pos1] > B[1]->_counts[pos2]) {
                  outputMer(oDAT, &mcd[0], b, B[0]->_checks[pos1], (u32bit)B[0]->_counts[pos1] - (u32bit)B[1]->_counts[pos2]);
                  itemsWritten++;
                }
                pos1++;
                pos2++;
              } else if (B[0]->_checks[pos1] < B[1]->_checks[pos2]) {
                outputMer(oDAT, &mcd[0], b, B[0]->_checks[pos1], (u32bit)B[0]->_counts[pos1]);
                itemsWritten++;
                pos1++;
              } else {
                outputMer(oDAT, &mcd[0], b, B[1]->_checks[pos2], (u32bit)B[1]->_counts[pos2]);
                itemsWritten++;
                pos2++;
              }
              break;
#if 0
            case PERSONALITY_MASK:
              if        (B[0]->_checks[pos1] == B[1]->_checks[pos2]) {
                outputMer(oDAT, &mcd[0], b, B[1]->_checks[pos2], (u32bit)B[1]->_counts[pos2]);
                itemsWritten++;
                pos1++;
                pos2++;
              } else if (B[0]->_checks[pos1] < B[1]->_checks[pos2]) {
                pos1++;
              } else {
                pos2++;
              }
              break;
#endif
          }
      }
    }

    //  write the number of entries
    //
    oIDX->putBits(itemsWritten, 32);

    //  read the next set of buckets
    //
    if (B[0]->_bucketID == b)
      B[0]->readBucket();

    if (B[1]->_bucketID == b)
      B[1]->readBucket();
  }

  for (u32bit i=0; i<args->mergeFilesLen; i++)
    delete B[i];

  delete [] B;

  delete [] mcd;

  delete oIDX;
  delete oDAT;

  for (u32bit i=0; i<args->mergeFilesLen; i++) {
    delete IDX[i];
    delete DAT[i];
  }
  delete [] IDX;
  delete [] DAT;
}
