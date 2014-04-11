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

#define QUERY_LENGTH     2048

int
main(int argc, char **argv) {
  bool   beVerbose      = false;
  uint32 chimeraOverlap = 5;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
    } else if (strncmp(argv[arg], "-o", 2) == 0) {
      chimeraOverlap = strtouint32(argv[++arg], 0L);
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
    }
    arg++;
  }

  intervalList   IL;
  intervalList   ILfull;
  uint32         ILid = 0;
  char           lastdefline[1024] = { 0 };

  uint32         numPts = 0;
  uint32         maxPts = 1024;
  uint32        *begPt = new uint32 [maxPts];
  uint32        *endPt = new uint32 [maxPts];

  uint32        *genBeg = new uint32 [maxPts];
  uint32        *genEnd = new uint32 [maxPts];

  uint32         queryLength = 0;

  char   spaces[QUERY_LENGTH+1];
  char   lines[QUERY_LENGTH+1];
  char   equals[QUERY_LENGTH+1];

  for (uint32 i=0; i<QUERY_LENGTH; i++) {
    spaces[i]  = ' ';
    lines[i]   = '-';
    equals[i]  = '=';
  }
  spaces[QUERY_LENGTH] = 0;
  lines[QUERY_LENGTH]  = 0;
  equals[QUERY_LENGTH] = 0;

  sim4polishReader *R = new sim4polishReader("-");
  sim4polish       *p = 0L;

  while (R->nextAlignment(p)) {
    if ((p->_estID != ILid) &&
        (lastdefline[0])) {

#if 0
      fprintf(stdout, "\n\n");

      fprintf(stdout, "IL "uint32FMT"\n", IL.numberOfIntervals());
      for (uint32 i=0; i<IL.numberOfIntervals(); i++)
        fprintf(stderr, "IL["uint32FMTW(3)"] "uint64FMT" "uint64FMT"\n", i, IL.lo(i), IL.hi(i));

      fprintf(stdout, "ILfull "uint32FMT"\n", ILfull.numberOfIntervals());
      for (uint32 i=0; i<ILfull.numberOfIntervals(); i++)
        fprintf(stderr, "ILfull["uint32FMTW(3)"] "uint64FMT" "uint64FMT"\n", i, ILfull.lo(i), ILfull.hi(i));
#endif

      IL.merge();
      ILfull.merge();

      if ((IL.numberOfIntervals() > 1) &&
          (ILfull.sumOfLengths() >= 0.9 * queryLength)) {
        fprintf(stdout, "%s\n", lastdefline);

        equals[queryLength] = 0;
        fprintf(stdout, "        %s\n", equals);
        equals[queryLength] = '=';

        //  Bubble sort the positions.
        //
        for (uint32 a=0; a<numPts; a++) {
          for (uint32 b=a+1; b<numPts; b++) {
            if ((begPt[a] > begPt[b]) ||
                ((begPt[a] == begPt[b]) && (endPt[a] > endPt[b]))) {
              uint32 x = begPt[a];
              uint32 y = endPt[a];
              begPt[a] = begPt[b];
              endPt[a] = endPt[b];
              begPt[b] = x;
              endPt[b] = y;

              x         = genBeg[a];
              y         = genEnd[a];
              genBeg[a] = genBeg[b];
              genEnd[a] = genEnd[b];
              genBeg[b] = x;
              genEnd[b] = y;
            }
          }
        }


        for (uint32 i=0; i<numPts && i<maxPts; i++) {
          if (begPt[i] >= QUERY_LENGTH) {
            fprintf(stdout, "WARNING:  Next line (begin) truncated to %d positions!\n", QUERY_LENGTH);
            begPt[i] = QUERY_LENGTH-1;
          }
          if (endPt[i] >= QUERY_LENGTH) {
            fprintf(stdout, "WARNING:  Next line (end) truncated to %d positions!\n", QUERY_LENGTH);
            endPt[i] = QUERY_LENGTH-1;
          }


          spaces[begPt[i]] = 0;
          lines[endPt[i] - begPt[i]] = 0;
          fprintf(stdout, uint32FMTW(3)"-"uint32FMTW(3)" %s%s ("uint32FMT","uint32FMT")\n",
                  begPt[i], endPt[i], spaces, lines, genBeg[i], genEnd[i]);
          spaces[begPt[i]] = ' ';
          lines[endPt[i] - begPt[i]] = '-';
        }

        fprintf(stdout, "\n\n");
      }  //  end of chimera detected

      IL.clear();
      ILfull.clear();
      numPts = 0;
    }

    strcpy(lastdefline, p->_estDefLine);
    ILid = p->_estID;

    queryLength = p->_estLen;

    uint32  beg = p->_exons[0]._estFrom - 1;
    uint32  end = p->_exons[p->_numExons-1]._estTo;

    if (numPts == maxPts) {
      fprintf(stdout, "Wow!  The next guy is a deep mapping!  I'm only showing the\n");
      fprintf(stdout, "first "uint32FMT" alignments.\n", maxPts);
    } else if (numPts < maxPts) {
      begPt[numPts] = beg;
      endPt[numPts] = end;
      genBeg[numPts] = p->_exons[0]._genFrom - 1;
      genEnd[numPts] = p->_exons[p->_numExons-1]._genTo;
    }
    numPts++;

    //fprintf(stdout, "beg,end = %d,%d\n", (int)beg, (int)end);

    if (end - beg > 2 * chimeraOverlap) {
      IL.add(beg + chimeraOverlap, end - beg - 2 * chimeraOverlap);
      ILfull.add(beg, end - beg);
    }
  }

  return(0);
}

