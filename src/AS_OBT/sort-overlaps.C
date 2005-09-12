#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>

#include "util++.H"
#include "overlap.H"
#include "constants.H"

//  Reads the output of overlap used to find trim points.  Bucket
//  sorts.

//  XXXX: OPTIMIZE: If one batch, skip writing the disk, and just load
//  things into core.  This will probably be a completely different
//  path than the bucket sort here.


int
overlap_t_sort(const void *a, const void *b) {
  overlap_t const *A = (overlap_t const *)a;
  overlap_t const *B = (overlap_t const *)b;
  if (A->Aiid < B->Aiid)  return(-1);
  if (A->Aiid > B->Aiid)  return(1);
  if (A->Biid < B->Biid)  return(-1);
  if (A->Biid > B->Biid)  return(1);
  return(0);
}



int
main(int argc, char **argv) {
  u64bit  memoryLimit = 512 * 1024 * 1024;
  u64bit  maxIID      = 0;
  u32bit  fileListLen = 0;
  u32bit  fileListMax = 10 * 1024;  //  If you run more than 10,000 overlapper jobs, you'll die.
  char  **fileList    = new char * [fileListMax];
  char   *tmpPath = ".";

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-memory") == 0) {
      memoryLimit = strtou64bit(argv[++arg], 0L) * 1024 * 1024;
    } else if (strcmp(argv[arg], "-maxiid") == 0) {
      //  Plus one, because it's all too easy to be given the
      //  last valid IID, rather than the number of IIDs.
      maxIID = strtou64bit(argv[++arg], 0L) + 1;
    } else if (strcmp(argv[arg], "-L") == 0) {
      //  The next arg is a file with the list of files to use
      arg++;
      errno = 0;
      FILE *F = fopen(argv[arg], "r");
      if (errno)
        fprintf(stderr, "Can't open '%s': %s\n", argv[arg], strerror(errno)), exit(1);

      //  Ugh, terrible.

      char  *line = new char [1024];
      fgets(line, 1024, F);
      while (!feof(F)) {
        chomp(line);
        fileList[fileListLen++] = line;
        //fprintf(stderr, "line=%s\n", line);
        line = new char [1024];
        fgets(line, 1024, F);
      }
      fclose(F);

      if (fileListLen >= fileListMax)
        fprintf(stderr, "Too many input files, increase fileListMax.\n"), exit(1);
    } else {
      //  Assume it's an input file
      //
      fileList[fileListLen++] = argv[arg];
    }
    arg++;
  }

  if (maxIID == 0) {
    fprintf(stderr, "usage: %s [-memory xMB] -maxiid # <ovl-file> <ovl-file> <...>\n");
    exit(1);
  }

  //  Based on filesize, guess the number of overlaps in the file.
  //  Unless someone changed the format string in overlap, there are
  //  exactly 57 bytes per line (per overlap).
  //
  u64bit  numOverlaps    = 0;

  for (u32bit i=0; i<fileListLen; i++) {
    numOverlaps += sizeOfFile(fileList[i]) / 57;
    //fprintf(stderr, "%s -- "u64bitFMT" -- %d\n", fileList[i], numOverlaps, (int)numOverlaps);
  }

  fprintf(stderr, "I think there are "u64bitFMT" overlaps in your input, so "u64bitFMT" overlaps to sort.\n", numOverlaps, 2*numOverlaps);
  fprintf(stderr, "You'll let me use "u64bitFMT" bytes of memory.\n", memoryLimit);

  u64bit  overlapsPerBatch   = memoryLimit / sizeof(overlap_t);
  fprintf(stderr, "overlapsPerBatch:    "u64bitFMT" (memory / sizeof(overlap))\n", overlapsPerBatch);

  //  We really need to work with twice the number of overlaps in the input.
  numOverlaps *= 2;

  //  +1 just to get around having fewer overlaps than IID's.
  u64bit  overlapsPerIID     = numOverlaps / maxIID + 1;
  fprintf(stderr, "maxIID:              "u64bitFMT"\n", maxIID);
  fprintf(stderr, "overlapsPerIID:      "u64bitFMT"\n", overlapsPerIID);

  u64bit  overlapIIDPerBatch = overlapsPerBatch / overlapsPerIID;
  fprintf(stderr, "overlapIIDPerBatch:  "u64bitFMT"\n", overlapIIDPerBatch);

  //  +1 again to get around having fewer overlaps than IID's
  u64bit  overlapBatches     = maxIID / overlapIIDPerBatch + 1;
  fprintf(stderr, "overlapBatches:      "u64bitFMT"\n", overlapBatches);

  //
  //  The first pass is to throw things into buckets.  Our buckets are
  //  completely on disk.  Better implementation would use a buffer.
  //

#ifdef RESTART
  overlapBatches = 10;
#endif

  FILE       **dumpFiles    = new FILE * [overlapBatches];
  u32bit      *dumpFilesLen = new u32bit [overlapBatches];


  //  XXXXXXXX  SKIP THIS FOR A RESTART OF MACAQUE
#ifndef RESTART

  for (u32bit i=0; i<overlapBatches; i++) {
    char name[1024];
    sprintf(name, "%s/sort-overlap.dump."u32bitFMTW(03), tmpPath, i);
    errno = 0;
    dumpFiles[i]    = fopen(name, "wb");
    if (errno)
      fprintf(stderr, "Failed to create '%s': %s\n", name, strerror(errno)), exit(1);
    dumpFilesLen[i] = 0;
  }

  overlap_t    overlap;

#ifdef SPEEDCOUNTER_H
  speedCounter  *C = new speedCounter("%7.2f Moverlaps -- %5.2f Moverlaps/second\r",
                                      1000000.0, 0x7fff, true);
  C->enableLiner();
#endif

  char line[1024];

  for (u32bit i=0; i<fileListLen; i++) {
    fprintf(stderr, "\nWorking on %s\n", fileList[i]);
    errno = 0;
    FILE *inFile = fopen(fileList[i], "r");
    if (errno)
      fprintf(stderr, "Failed to open input '%s': %s\n", fileList[i], strerror(errno)), exit(1);

    fgets(line, 1024, inFile);
    while (!feof(inFile)) {
      overlap.decode(line, false);

      if ((overlap.Aiid >= maxIID) || (overlap.Biid >= maxIID))
        fprintf(stderr, "ERROR:  Too many IID's!  Input is:\n        %s\n", line), exit(1);
      if (overlap.Aiid / overlapIIDPerBatch >= overlapBatches)
        fprintf(stderr, "ERROR:  Aiid overflowed the batch!  Input is:\n        %s\n", line), exit(1);
      if (overlap.Biid / overlapIIDPerBatch >= overlapBatches)
        fprintf(stderr, "ERROR:  Biid overflowed the batch!  Input is:\n        %s\n", line), exit(1);


      if (overlap.acceptable()) {
        overlap.dump(dumpFiles[overlap.Aiid / overlapIIDPerBatch]);
        dumpFilesLen[overlap.Aiid / overlapIIDPerBatch]++;

        overlap.decode(line, true);
        overlap.dump(dumpFiles[overlap.Aiid / overlapIIDPerBatch]);
        dumpFilesLen[overlap.Aiid / overlapIIDPerBatch]++;
      }

      fgets(line, 1024, inFile);
#ifdef SPEEDCOUNTER_H
      C->tick();
#endif
    }

    fclose(inFile);
  }
#ifdef SPEEDCOUNTER_H
  delete C;
#endif

  for (u32bit i=0; i<overlapBatches; i++) {
    fclose(dumpFiles[i]);
  }

#endif

  //
  //  Read each bucket, sort it, and dump it to the output
  //

  FILE *mergeFile = stdout;
  //errno = 0;
  //FILE *mergeFile = fopen("/scratch/sort-overlap-out", "wb");
  //if (errno)
  //  fprintf(stderr, "Failed to create '%s': %s\n", name, strerror(errno)), exit(1);

  for (u32bit i=0; i<overlapBatches; i++) {
    char name[1024];
    sprintf(name, "%s/sort-overlap.dump."u32bitFMTW(03), tmpPath, i);

#ifdef RESTART
    dumpFilesLen[i] = sizeOfFile(name) / sizeof(overlap_t);
#endif

    overlap_t *overlapsort = new overlap_t [dumpFilesLen[i]];

    errno = 0;
    dumpFiles[i] = fopen(name, "rb");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for reading: %s\n", name, strerror(errno)), exit(1);

    fprintf(stderr, "Read %s at %f\n", name, getTime());
    fread(overlapsort, sizeof(overlap_t), dumpFilesLen[i], dumpFiles[i]);
    if (errno)
      fprintf(stderr, "Failed to read "u32bitFMT" overlaps from '%s': %s\n", dumpFilesLen[i], name, strerror(errno)), exit(1);

    fclose(dumpFiles[i]);

    fprintf(stderr, "Sort %s at %f\n", name, getTime());
    qsort(overlapsort, dumpFilesLen[i], sizeof(overlap_t), overlap_t_sort);

    fprintf(stderr, "Dump %s at %f\n", name, getTime());
#ifdef ASCII_OVERLAPS
    for (u32bit x=0; x<dumpFilesLen[i]; x++)
      overlapsort[x].print(mergeFile);
#else
    fwrite(overlapsort, sizeof(overlap_t), dumpFilesLen[i], mergeFile);
#endif

    delete [] overlapsort;

  }

  if (mergeFile != stdout)
    fclose(mergeFile);

  //  Clean up our temporary files
  //
  for (u32bit i=0; i<overlapBatches; i++) {
    char name[1024];
    sprintf(name, "/scratch/sort-overlap.dump."u32bitFMTW(03), i);
    unlink(name);
  }

  exit(0);
}

