#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "bri++.H"
#include "sim4polish.h"
#include "sim4polishBuilder.H"
#include "sim4polishFile.H"


sim4polishList*
filterByQuality(u32bit minI, u32bit minC, sim4polishList *A) {
  sim4polishList *l = new sim4polishList;

  for (u32bit i=0; i<A->length(); i++)
    if (((*A)[i]->percentIdentity  >= minI) &&
        ((*A)[i]->querySeqIdentity >= minC))
      l->push( s4p_copyPolish((*A)[i]) );

  delete A;

  return(l);
}


//  Build an interval list with all exons (from both guys), merge
//  overlapping regions, compute the length, subtract from the
//  total.
//
u32bit
findOverlap(sim4polish *A, sim4polish *B) {

  if ((A->genID != B->genID) || (A->matchOrientation != B->matchOrientation))
    return(0);

  u32bit        length = 0;
  u32bit        total  = 0;
  intervalList  IL;

  for (u32bit i=0; i<A->numExons; i++) {
    length = A->exons[i].genTo - A->exons[i].genFrom + 1;
    total  += length;
    IL.add(A->genLo + A->exons[i].genTo, length);
  }

  for (u32bit i=0; i<B->numExons; i++) {
    length = B->exons[i].genTo - B->exons[i].genFrom + 1;
    total  += length;
    IL.add(B->genLo + B->exons[i].genTo, length);
  }

  IL.merge();

  return(total - IL.sumOfLengths());
}


void
comparePolishFiles(int argc, char **argv) {
  u32bit      minI  = atoi(argv[1]);
  u32bit      minC  = atoi(argv[2]);
  char       *Apath = argv[3];
  char       *Bpath = argv[4];

  sim4polishFile *Afile = new sim4polishFile(Apath);
  sim4polishFile *Bfile = new sim4polishFile(Bpath);

  //  Force index builds
  //
  Afile->setPosition(0);
  Bfile->setPosition(0);


  u32bit  goodOverlap = 0;
  u32bit  multipleInA = 0;
  u32bit  multipleInB = 0;


  //  Ask both for the largest EST iid seen, then iterate over those.
  //
  u32bit  largestIID = 100000;

  for (u32bit iid=0; iid<largestIID; iid++) {
    sim4polishList *A = Afile->getEST(iid);
    sim4polishList *B = Bfile->getEST(iid);

    //  Filter by quality
    A = filterByQuality(minI, minC, A);
    B = filterByQuality(minI, minC, B);

    //  fill out the overlap matrix

    u32bit  **overlap = new u32bit* [A->length()];
    overlap[0] = new u32bit [A->length() * B->length()];
    for (u32bit i=1; i<A->length(); i++)
      overlap[i] = overlap[i-1] + B->length();

    for (u32bit a=0; a<A->length(); a++)
      for (u32bit b=0; b<B->length(); b++) {
        overlap[a][b] = findOverlap((*A)[a], (*B)[b]);
        //fprintf(stderr, "overlap[%2d][%2d] = "u32bitFMT"\n", a, b, overlap[a][b]);
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


int
main(int argc, char **argv) {

  if (argc != 5) {
    fprintf(stderr, "       %s <minid> <mincov> <polishes-file-1> <polishes-file-2>\n", argv[0]);
    exit(1);
  }

  comparePolishFiles(argc, argv);
}


