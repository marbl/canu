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
//  Arguments (are horrible, whatcha gonna do about it?):
//
//    -i min percent id (default 95)
//    -c min percent coverage (default 50)
//    -a polishes file 1
//    -b polishes file 2
//
//    -O output file
//    -A polishes from file a
//    -B polishes from file b
//    -L log output file
//
//  The default is to output on stdout.  -A, -B write the overlapped
//  matches in the same order as the output-file.  -L reports a log of
//  weird stuff.
//
//  Output is tab-delimited:
//    ESTiid ESTlen  overlap  A%id A%cov #cdnagaps #exons  B%id B%cov #cdnagaps #exons
//

//  e,b good:113408 Anovel:1547 Amulti: 216 Bnovel:83220 Bmulti:  13 hairy: 270
//  b,e good:113408 Bnovel:1547 Bmulti:  76 Anovel:83220 Amulti:  69 hairy: 129

//  b,e good:113408 Anovel:1547 Amulti: 216 Bnovel:83220 Bmulti:  13 hairy:  54
//  e,b good:113408 Bnovel:1547 Bmulti:  76 Anovel:83220 Amulti:  69 hairy:  60


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

  //  goodOverlap   -- match in A maps uniquely to B and likewise.
  //
  //  novelInA      -- a match in A has no counterpart in B.
  //  novelInB      -- similar for B.
  //
  //  multipleInA   -- a match in B maps to multiple things in A.
  //  multipleInB   -- similar for A.
  //
  //  multipleInA requires that the matches in A map only to the single
  //  match in B.
  //
  //  hairyOverlap  -- multiple matches in both.
  //
  u32bit           goodOverlap  = 0;  //  the number of lines in the output
  u32bit           novelInA     = 0;
  u32bit           novelInB     = 0;
  u32bit           multipleInA  = 0;
  u32bit           multipleInB  = 0;
  u32bit           hairyOverlap = 0;

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


  //  Find the largest IID
  //
  u32bit  largestIID = Afile->maxIID();
  if (largestIID < Bfile->maxIID())
    largestIID = Bfile->maxIID();


  //  Iterate over all the ESTs.

  for (u32bit iid=0; iid<largestIID; iid++) {
    sim4polishList *A = Afile->getEST(iid);
    sim4polishList *B = Bfile->getEST(iid);

    //  Filter by quality.
    A->filterByQuality(minI, minC);
    B->filterByQuality(minI, minC);

    //  fill out the overlap matrix

    u32bit  **overlap = new u32bit* [A->length()];
    overlap[0] = new u32bit [A->length() * B->length()];
    for (u32bit i=1; i<A->length(); i++)
      overlap[i] = overlap[i-1] + B->length();

    for (u32bit a=0; a<A->length(); a++)
      for (u32bit b=0; b<B->length(); b++)
        overlap[a][b] = findOverlap((*A)[a], (*B)[b]);



    //  For each match in A, find the one in B that overlaps it.

    for (u32bit a=0; a<A->length(); a++) {

      u32bit Boverlaps = 0;
      u32bit theBovl   = 0;

      //  Count the number of things we overlap in B.
      for (u32bit b=0; b<B->length(); b++) {
        if (overlap[a][b]) {
          Boverlaps++;
          theBovl = b;
        }
      }

      //  No overlaps --> novelInA
      //
      if (Boverlaps == 0) {
        novelInA++;
          //  XXX OUTPUT
      }

      //  One overlap?  Possibly a goodOverlap, possibly multipleInA.
      //
      else if (Boverlaps == 1) {

        //  Count the number of overlaps the guy in B has with A.  If
        //  1, it's a goodOverlap, else it's a multipleInA.

        u32bit Aoverlaps = 0;
        for (u32bit x=0; x<A->length(); x++)
          if (overlap[x][theBovl])
            Aoverlaps++;

        if (Aoverlaps == 1) {
          goodOverlap++;
          //  XXX OUTPUT
        } else {
          multipleInA++;
          //  XXX OUTPUT
        }
      }

      //  More than one overlap?  Possibly multipleInB, possilby hairyOverlap.
      //
      else if (Boverlaps > 1) {

        //  Examine the overlaps the guys in B have with A.  If any
        //  are to something other than a, it's hairy, otherwise it's
        //  multipleInB.

        bool hairy = false;

        for (u32bit b=0; b<B->length(); b++)
          if (overlap[a][b]) {

            //  Found an overlap, check all of A and see if this b overlaps
            //  anything other than this a

            for (u32bit x=0; x<A->length(); x++)
              if ((overlap[x][b]) && (a != x))
                hairy = true;
          }

        if (hairy == false) {
          multipleInB++;
          //  XXX OUTPUT
        } else {
          hairyOverlap++;
          //  XXX OUTPUT
        }
      }
    }

    //  Look for novelInB separately.
    //
    for (u32bit b=0; b<B->length(); b++) {
      u32bit Boverlaps = 0;

      for (u32bit a=0; a<A->length(); a++)
        if (overlap[a][b])
          Boverlaps++;

      if (Boverlaps == 0)
        novelInB++;
    }


    if ((iid % 1234) == 0) {
      fprintf(stderr, "IID:"u32bitFMTW(8)" A:"u32bitFMTW(4)" B:"u32bitFMTW(4)"  good:"u32bitFMTW(4)" Anovel:"u32bitFMTW(4)" Amulti:"u32bitFMTW(4)" Bnovel:"u32bitFMTW(4)" Bmulti:"u32bitFMTW(4)" hairy:"u32bitFMTW(4)"\r",
              iid,
              A->length(), B->length(),
              goodOverlap, novelInA, multipleInA, novelInB, multipleInB, hairyOverlap);
      fflush(stderr);
    }

    //  Do it again for B -> A, iterate until it's stable.

    //  Have now a match in A and a match in B.  Decide if they are
    //  close enough to be counted as the same.  They need to overlap
    //  by at least XX% of the longer guy.



    delete [] overlap[0];
    delete [] overlap;

    delete A;
    delete B;
  }

  fprintf(stderr, "\ngood:"u32bitFMTW(4)" Anovel:"u32bitFMTW(4)" Amulti:"u32bitFMTW(4)" Bnovel:"u32bitFMTW(4)" Bmulti:"u32bitFMTW(4)" hairy:"u32bitFMTW(4)"\n",
          goodOverlap, novelInA, multipleInA, novelInB, multipleInB, hairyOverlap);


  delete Afile;
  delete Bfile;
}
