#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"

//  Attempts to look for query that are chimeric.  It is assumed that
//  your query have been mapped to a target reference genome, such
//  that little pieces will be mapped.  The heuristic used is simple.
//  The mapping intervals are merged together, and if there are two
//  blocks that do not overlap, then it is chimeric.  Intervals are
//  decreased by 3bp before merging.

#define QUERY_LENGTH     256

int
main(int argc, char **argv) {
  bool   beVerbose      = false;
  u32bit chimeraOverlap = 5;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
    } else if (strncmp(argv[arg], "-o", 2) == 0) {
      chimeraOverlap = strtou32bit(argv[++arg], 0L);
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
    }
    arg++;
  }

  intervalList   IL;
  u32bit         ILid = 0;
  char           lastdefline[1024] = { 0 };

  u32bit         numPts = 0;
  u32bit         maxPts = 1024;
  u32bit        *begPt = new u32bit [maxPts];
  u32bit        *endPt = new u32bit [maxPts];

  u32bit         queryLength = 0;

  char   spaces[QUERY_LENGTH];
  char   lines[QUERY_LENGTH];
  char   equals[QUERY_LENGTH];

  for (u32bit i=0; i<QUERY_LENGTH; i++) {
    spaces[i] = ' ';
    lines[i]  = '-';
    equals[i]  = '=';
  }

  while (!feof(stdin)) {
    sim4polish *p = s4p_readPolish(stdin);

    if ((p == 0L) || (p->estID != ILid)) {
      if (lastdefline[0]) {
        IL.merge();
        if (IL.numberOfIntervals() > 1) {
          fprintf(stdout, "%s\n", lastdefline);

          equals[queryLength] = 0;
          fprintf(stdout, "        %s\n", equals);
          equals[queryLength] = '=';

          //  Bubble sort the positions.
          //
          for (u32bit a=0; a<numPts; a++) {
            for (u32bit b=a+1; b<numPts; b++) {
              if ((begPt[a] > begPt[b]) ||
                  ((begPt[a] == begPt[b]) && (endPt[a] > endPt[b]))) {
                u32bit x = begPt[a];
                u32bit y = endPt[a];
                begPt[a] = begPt[b];
                endPt[a] = endPt[b];
                begPt[b] = x;
                endPt[b] = y;
              }
            }
          }


          for (u32bit i=0; i<numPts && i<maxPts; i++) {
            if (begPt[i] > 255) {
              fprintf(stdout, "WARNING:  Next line (begin) truncated to %d positions!\n", QUERY_LENGTH);
              begPt[i] = QUERY_LENGTH-1;
            }
            if (endPt[i] > 255) {
              fprintf(stdout, "WARNING:  Next line (end) truncated to %d positions!\n", QUERY_LENGTH);
              endPt[i] = QUERY_LENGTH-1;
            }


            spaces[begPt[i]] = 0;
            lines[endPt[i] - begPt[i]] = 0;
            fprintf(stdout, u32bitFMTW(3)"-"u32bitFMTW(3)" %s%s\n", begPt[i], endPt[i], spaces, lines);
            spaces[begPt[i]] = ' ';
            lines[endPt[i] - begPt[i]] = '-';
          }

          fprintf(stdout, "\n\n");
        }  //  end of chimera detected
      }

      IL.clear();
      numPts = 0;
    }

    if (p != 0L) {
      strcpy(lastdefline, p->estDefLine);
      ILid = p->estID;

      queryLength = p->estLen;

      u32bit  beg = p->exons[0].estFrom - 1;
      u32bit  end = p->exons[p->numExons-1].estTo;

      if (numPts == maxPts) {
        fprintf(stdout, "Wow!  The next guy is a deep mapping!  I'm only showing the\n");
        fprintf(stdout, "first "u32bitFMT" alignments.\n", maxPts);
      } else if (numPts < maxPts) {
        begPt[numPts] = beg;
        endPt[numPts] = end;
      }
      numPts++;

      beg += chimeraOverlap;
      end -= chimeraOverlap;

      if (beg < end)
        IL.add(beg, end-beg);

      s4p_destroyPolish(p);
    }
  }

  return(0);
}

