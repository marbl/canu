#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "util++.H"
#include "overlap.H"

//  This used to be sort-overlap-trim.C -- which would read raw overlaps
//  from overlap, sort them, and then do everything here.

//  Reads the output of overlap used to find trim points.  Bucket
//  sorts.
//
//  An overlap is
//    fragA fragB  ori  leftA rightA lenA  leftB rightB lenB  qual
//
//  Discards overlap if it is less than 100bp
//  Discards overlap if percent error is more than 1.5
//  (But see the code for what it really does)
//
//  The output contains:
//      fragment id
//      the list of 5' and 3' points
//        mode
//        max/min
//        max/min with more than one occurance
//

#define MAX_OVERLAPS_PER_FRAG   (16 * 1024 * 1024)


//  sort the position values on the 5' end -- this sorts increasingly
int
position_compare5(const void *a, const void *b) {
  u32bit  A = (u32bit)*((u32bit *)a);
  u32bit  B = (u32bit)*((u32bit *)b);

  if (A < B)  return(-1);
  if (A > B)  return(1);
  return(0);
}

//  sort the position values on the 3' end -- this sorts decreasingly
int
position_compare3(const void *a, const void *b) {
  u32bit  A = (u32bit)*((u32bit *)a);
  u32bit  B = (u32bit)*((u32bit *)b);

  if (A < B)  return(1);
  if (A > B)  return(-1);
  return(0);
}




void
sortAndOutput(u32bit fid, u32bit numOvl, u32bit *left, u32bit *right) {

  qsort(left,  numOvl, sizeof(u32bit), position_compare5);
  qsort(right, numOvl, sizeof(u32bit), position_compare3);

  //  XXX:  Print the 5' and 3' stuff
  //
  //  We might as well find the
  //    min/max
  //    min/max with more than one hit
  //    mode
  //  since we have everything here
  //
  //  minN    -- minimum value we've ever seen
  //  minmN   -- minimum value we've seen more than once
  //  minmNc  -- number of times we've seen minm
  //  modeN   -- mode
  //  modeNc  -- number of times we've seen the mode
  //  modeNt  -- temp copy of the mode
  //  modeNtc -- temp copy of the mode, number of times
  //
  u32bit  min5, minm5, minm5c,  mode5, mode5c,  mode5t, mode5tc;
  u32bit  max3, maxm3, maxm3c,  mode3, mode3c,  mode3t, mode3tc;

  min5 = left[0];
  max3 = right[0];

  minm5 = 9999;       minm5c  = 0;
  maxm3 = 9999;       maxm3c  = 0;

  mode5 = left[0];    mode5c  = 1;
  mode3 = right[0];   mode3c  = 1;

  mode5t = left[0];   mode5tc = 1;
  mode3t = right[0];  mode3tc = 1;

  for (u32bit i=1; i<numOvl; i++) {

    //  5' end.  We scan the list, remembering the best mode
    //  we've seen so far.  When a better one arrives, we copy
    //  it to the saved one -- and keep copying it as it gets
    //  better.
    //
    if (mode5t == left[i]) {  //  Same mode?  Count.
      mode5tc++;
    } else {
      mode5t  = left[i];  //  Different mode, restart.
      mode5tc = 1;
    }
    if (mode5tc > mode5c) {  //  Bigger mode?  Save it.
      mode5  = mode5t;
      mode5c = mode5tc;
    }

    //  If our mode is more than one and we've not seen a multiple hit before
    //  save this position.
    //
    if ((mode5c > 1) && (minm5 == 9999))
      minm5  = mode5;
    if (minm5 == mode5)
      minm5c = mode5c;


    //  Do it all again for the 3' -- remember that we've
    //  sorted this decreasingly.


    if (mode3t == right[i]) {  //  Same mode?  Count.
      mode3tc++;
    } else {
      mode3t  = right[i];  //  Different mode, restart.
      mode3tc = 1;
    }
    if (mode3tc > mode3c) {  //  Bigger mode?  Save it.
      mode3  = mode3t;
      mode3c = mode3tc;
    }
    if ((mode3c > 1) && (maxm3 == 9999))
      maxm3  = mode3;
    if (maxm3 == mode3)
      maxm3c = mode3c;
  }


  //  Output!
  //
  fprintf(stdout, u32bitFMT"  "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"  "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT" "u32bitFMT"",
          fid,
          min5, minm5, minm5c, mode5, mode5c,
          max3, maxm3, maxm3c, mode3, mode3c);

#if 0
  //  Save all the overlaps too
  //
  fprintf(stdout, "  "u32bitFMT"", numOvl);
  for (u32bit i=0; i<numOvl; i++)
    fprintf(stdout, "  "u32bitFMT" "u32bitFMT"",
            left[i], right[i]);
#endif

  fprintf(stdout, "\n");
}



bool
readOverlap(FILE *file, overlap_t &ovl) {
  static char line[1024];

  if (feof(file))
    return(false);

#ifdef ASCII_OVERLAPS
  fgets(line, 1024, stdin);
  ovl.decode(line, false);
#else
  ovl.load(file);
#endif

  if (feof(file))
    return(false);

  return(true);
}



int
main(int argc, char **argv) {

  if (argc != 1) {
    fprintf(stderr, "usage: %s < overlap-trim-results\n", argv[0]);
    exit(1);
  }

  u32bit   idAlast     = 0;
  u32bit   numOverlaps = 0;
  u32bit  *left        = new u32bit [MAX_OVERLAPS_PER_FRAG];
  u32bit  *right       = new u32bit [MAX_OVERLAPS_PER_FRAG];

  speedCounter  *C = new speedCounter("%7.2f Moverlaps -- %5.2f Moverlaps/second\r",
                                      1000000.0, 0x7fff, true);
  C->enableLiner();

  overlap_t ovl;
  while (readOverlap(stdin, ovl)) {
    u64bit   wa, wb;

    //  We should get rid of these aliases...
    u32bit idA    = ovl.Aiid;
    //u32bit idB    = strtou32bit(W[1], 0L);
    char   ori    = (ovl.ori) ? 'f' : 'r';
    u32bit leftA  = ovl.Abeg;
    u32bit rightA = ovl.Aend;
    //u32bit lenA   = strtou32bit(W[5], 0L);
    u32bit leftB  = ovl.Bbeg;
    u32bit rightB = ovl.Bend;
    //u32bit lenB   = strtou32bit(W[8], 0L);
    double  error  = ovl.erate;

    //  If we see a different idA than we had last time, process
    //  the previous read.
    //
    if ((idAlast != idA) && (numOverlaps > 0)) {
      sortAndOutput(idAlast, numOverlaps, left, right);
      numOverlaps = 0;
    }

    idAlast     = idA;

    //  Save the location of the overlap in the A fragment if it's
    //  an acceptable 
    //
    if (ovl.acceptable()) {
      left[numOverlaps]  = leftA;
      right[numOverlaps] = rightA;
      numOverlaps++;
    }

    C->tick();
  }

  //  Don't forget to do the last batch!
  if (numOverlaps > 0) {
    sortAndOutput(idAlast, numOverlaps, left, right);
    numOverlaps = 0;
  }

  delete [] left;
  delete [] right;
}
