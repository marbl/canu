#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "meryl.H"
#include "mcDescription.H"
#include "outputMer.H"
#include "britime.H"



void
adjustHeap(u64bit *M, s64bit i, s64bit n) {
  u64bit   m = M[i];
  s64bit    j = (i << 1) + 1;  //  let j be the left child

  while (j < n) {
    if (j<n-1 && M[j] < M[j+1])
      j++;                   //  j is the larger child

    if (m >= M[j])           //  a position for M[i] has been found
      break;

    M[(j-1)/2] = M[j];       //  Move larger child up a level

    j = (j << 1) + 1;
  }

  M[(j-1)/2] = m;
}



void
sortAndOutput(merylArgs      *args,
              mcDescription  *mcd,
              u64bit         *chck,
              u64bit         *hash) {

  u64bit   m     = u64bitONE << mcd->_tableSizeInBits;
  u32bit   count = 0;
  u32bit   items = 0;
  u64bit  *sortedList    = 0L;
  u32bit   sortedListMax = 0;
  u32bit   sortedListLen = 0;


  if (args->beVerbose)
    fprintf(stderr, " 6) Writing output.\n");

  //  Open the output files
  //
  char *outpath = new char [strlen(args->outputFile) + 17];

  sprintf(outpath, "%s.mcidx", args->outputFile);
  bitPackedFileWriter *IDX = new bitPackedFileWriter(outpath);

  sprintf(outpath, "%s.mcdat", args->outputFile);
  bitPackedFileWriter *DAT = new bitPackedFileWriter(outpath);

  delete [] outpath;


  //  Write the parameters to the DAT file.  It probably doesn't
  //  matter which file we write these to, as we can't really do
  //  random access anyway.
  //
  mcd->write(DAT);

  speedCounter  C(" %7.2f Mbuckets -- %5.2f Mbuckets/second\r", 1000000.0, 0x1fffff, args->beVerbose);

  //  For each bucket, sort it.  The output is done
  //  in the sort.
  //
  for (u64bit B=0, b=0; b<m; b++) {
    C.tick();

    u64bit st = getDecodedValue(hash, B, mcd->_hashWidth);
    B        += mcd->_hashWidth;
    u64bit ed = getDecodedValue(hash, B, mcd->_hashWidth);

    if (ed < st)
      fprintf(stderr, "ERROR: Bucket "u64bitFMT" ends before it starts!  start="u64bitFMT" end="u64bitFMT"\n", b, st, ed);

    sortedListLen = (u32bit)(ed - st);

    count = 0;
    items = 0;

    if (sortedListLen > 0) {

      //  Allocate more space, if we need to.
      //
      if (sortedListLen > sortedListMax) {
        delete [] sortedList;
        sortedList    = new u64bit [sortedListLen + 1];
        sortedListMax = sortedListLen;
      }

      //  Unpack the check values
      //
      for (u64bit i=st, J=st*mcd->_chckBits; i<ed; i++, J += mcd->_chckBits)
        sortedList[i-st] = getDecodedValue(chck, J, mcd->_chckBits);

      //  Sort if there is more than one item
      //
      if (sortedListLen > 1) {

        //  Create the heap of lines.
        //
        for (s64bit t=(sortedListLen-2)/2; t>=0; t--)
          adjustHeap(sortedList, t, sortedListLen);

        //  Interchange the new maximum with the element at the end of the tree
        //
        for (s64bit t=sortedListLen-1; t>0; t--) {
          u64bit           tv = sortedList[t];
          sortedList[t]      = sortedList[0];
          sortedList[0]      = tv;

          adjustHeap(sortedList, 0, t);
        }
      }


      //  Scan the list of sorted mers, counting them.  Whenever we 
      //  know the count, output it.
      //
      count = 1;
      if (sortedListLen > 0) {
        for (u32bit t=1; t<sortedListLen; t++) {
          if (sortedList[t] != sortedList[t-1]) {
            if ((args->lowCount <= count) && (count <= args->highCount)) {
              outputMer(DAT, mcd, b, sortedList[t-1], count);
              items++;
            }
            count = 0;
          }

          count++;
        }

        if ((args->lowCount <= count) && (count <= args->highCount)) {
          outputMer(DAT, mcd, b, sortedList[sortedListLen-1], count);
          items++;
        }
      }
    }

    //  Output the index
    //
    IDX->putBits(items, 32);
  }

  delete DAT;
  delete IDX;

  if (args->beVerbose)
    fprintf(stderr, "\n");
}

