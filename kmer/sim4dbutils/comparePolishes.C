#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "bri++.H"
#include "sim4polish.h"
#include "sim4polishBuilder.H"
#include "sim4polishFile.H"

//
//  Matches two sets of polishes to each other using a simple overlap
//  heuristic.
//
//  Output is:
//    AmatchIID BmatchIID  ESTlen  A%id A%cov  B%id B%cov  overlap
//
//  Arguments:
//   
//    -i min-percent-id (default 95)
//    -c min-percent-coverage (default 50)
//    -a polishes-file-1
//    -b polishes-file-1
//
//    -O output-file
//    -A output-polishes-from-file-a
//    -B output-polishes-from-file-b
//    -L log-output-file
//
//  The default is to output on stdout.  -A, -B write the overlapped
//  matches in order.  -L reports a log of weird stuff.
//

u32bit  findOverlap(sim4polish *A, sim4polish *B);

int
main(int argc, char **argv) {
  u32bit           minI = 95;
  u32bit           minC = 50;
  sim4polishFile  *Afile = 0L;
  sim4polishFile  *Bfile = 0L;
  FILE            *outfile = stdout;
  FILE            *Aout = 0L;
  FILE            *Bout = 0L;

  //  Stats
  u32bit           goodOverlap = 0;
  u32bit           multipleInA = 0;
  u32bit           multipleInB = 0;

  int arg=1;
  while(arg < argc) {
    if        (strcmp(argv[arg], "-i") == 0) {
      minI = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      minC = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-a") == 0) {
      Afile = new sim4polishFile(argv[++arg]);
    } else if (strcmp(argv[arg], "-b") == 0) {
      Bfile = new sim4polishFile(argv[++arg]);
    } else if (strcmp(argv[arg], "-O") == 0) {
      errno = 0L;
      outfile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for output: %s\n", argv[arg], strerror(errno)), exit(1);
    } else if (strcmp(argv[arg], "-A") == 0) {
      errno = 0L;
      Aout = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for output: %s\n", argv[arg], strerror(errno)), exit(1);
    } else if (strcmp(argv[arg], "-B") == 0) {
      errno = 0L;
      Bout = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for output: %s\n", argv[arg], strerror(errno)), exit(1);
    }
    arg++;
  }

  if ((Afile == 0L) || (Bfile == 0L)) {
    fprintf(stderr, "usage: %s [read-the-code]\n");
    exit(1);
  }

  //  Force index builds
  //
  Afile->setPosition(0);
  Bfile->setPosition(0);


  //  Ask both for the largest EST iid seen, then iterate over those.
  //
  u32bit  largestIID = Afile->maxIID();
  if (largestIID < Bfile->maxIID())
    largestIID = Bfile->maxIID();

  for (u32bit iid=0; iid<largestIID; iid++) {
    sim4polishList *A = Afile->getEST(iid);
    sim4polishList *B = Bfile->getEST(iid);

    //  Filter by quality
    A->filterByQuality(minI, minC);
    B->filterByQuality(minI, minC);

    //  fill out the overlap matrix

    u32bit  **overlap = new u32bit* [A->length()];
    overlap[0] = new u32bit [A->length() * B->length()];
    for (u32bit i=1; i<A->length(); i++)
      overlap[i] = overlap[i-1] + B->length();

    for (u32bit a=0; a<A->length(); a++)
      for (u32bit b=0; b<B->length(); b++) {
        overlap[a][b] = findOverlap((*A)[a], (*B)[b]);

#if 0
        if (overlap[a][b] > 0) {
          fprintf(stderr, "overlap[%2d][%2d] = "u32bitFMT"\n", a, b, overlap[a][b]);
          //s4p_printPolish(stderr, (*A)[a], S4P_PRINTPOLISH_NORMALIZED | S4P_PRINTPOLISH_MINIMAL);
          //s4p_printPolish(stderr, (*B)[b], S4P_PRINTPOLISH_NORMALIZED | S4P_PRINTPOLISH_MINIMAL);
        }
#endif
      }

    //  For each match in A, find the one in B that overlaps it.  If
    //  there is more than one, record an error, and merge the B's
    //  into one match.
    //

    for (u32bit a=0; a<A->length(); a++) {

      //  Count the number of things we overlap.  We count the number
      //  of times that this A overlaps any B, and if exactly once,
      //  then we count the number of times that B overlaps any A.  If
      //  also once, we declare a match.

      u32bit Aoverlaps = 0;
      u32bit theAovl   = 0;
      u32bit Boverlaps = 0;
      u32bit theBovl   = 0;

      for (u32bit b=0; b<B->length(); b++) {
        if (overlap[a][b]) {
          Aoverlaps++;
          theAovl = b;
        }
      }

      //  Continue only if we found an overlap in B.  If we didn't there is no match, period.

      if (Aoverlaps > 0) {

        for (u32bit i=0; i<A->length(); i++) {
          if (overlap[i][theAovl]) {
            Boverlaps++;
            theBovl = i;
          }
        }

        //  By definition if there is an A overlap, there will be a(t
        //  least one) B overlap.

        if ((Aoverlaps == 1) && (Boverlaps == 1)) {
          //  Got a good overlap!  match 'theAovl' and match 'theBovl'
          //  are probably the same.
          //
          goodOverlap++;
        } else if ((Aoverlaps > 1) || (Boverlaps > 1)) {
          //  Dang!  Somebody matched more than once.
          //
          if (Aoverlaps > 1)
            multipleInB++;
          if (Boverlaps > 1)
            multipleInA++;

#if 0
          //  This will print out details for the multiple match.

          fprintf(stderr, "Multiple overlaps!\n");

          fprintf(stderr, "A="u32bitFMT" with overlaps of ", a);
          for (u32bit x=0; x<B->length(); x++)
            fprintf(stderr, " B="u32bitFMT":"u32bitFMT, x, overlap[a][x]);
          fprintf(stderr, "\n");

          s4p_printPolish(stderr, (*A)[a], S4P_PRINTPOLISH_NORMALIZED | S4P_PRINTPOLISH_MINIMAL);

          for (u32bit x=0; x<B->length(); x++) {
            if (overlap[a][x]) {
              fprintf(stderr, "B="u32bitFMT" with overlaps of ", x);
              for (u32bit y=0; y<A->length(); y++)
                fprintf(stderr, " A="u32bitFMT":"u32bitFMT, y, overlap[y][x]);
              fprintf(stderr, "\n");

              s4p_printPolish(stderr, (*B)[x], S4P_PRINTPOLISH_NORMALIZED | S4P_PRINTPOLISH_MINIMAL);
            }
          }
#endif


        } else {
          fprintf(stderr, "Huh?!  I shouldn't be printing this!\n");
          exit(1);
        }
      }
    }

    fprintf(stderr, "IID:"u32bitFMTW(8)" A:"u32bitFMTW(4)" B:"u32bitFMTW(4)"  multiInA:"u32bitFMT" multiInB:"u32bitFMT" good:"u32bitFMT"\n",
            iid,
            A->length(), B->length(),
            multipleInA, multipleInB, goodOverlap);

    //  Do it again for B -> A, iterate until it's stable.

    //  Have now a match in A and a match in B.  Decide if they are
    //  close enough to be counted as the same.  They need to overlap
    //  by at least XX% of the longer guy.



    delete [] overlap[0];
    delete [] overlap;

    delete A;
    delete B;
  }

  delete Afile;
  delete Bfile;
}


