#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "bri++.H"
#include "sim4polish.h"
#include "sim4polishBuilder.H"
#include "sim4polishFile.H"

//  Remove redundant polishes from an input set.
//
//  Redundancy is defined as two polishes that overlap on the genome.
//  Any amount of overlap is redundant.
//
//  The longest of the overlapping matches is saved.



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
}


int
main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "usage: %s <polishes-file>\n", argv[0]);
    fprintf(stderr, "(yes, you _must_ give it a file.  stdin is not possible.)\n");
    exit(1);
  }

  u32bit  ESTaffected    = 0;
  u32bit  matchesRemoved = 0;

  //  Open a polishFile and force the index to build
  //
  sim4polishFile *Afile = new sim4polishFile(argv[1]);
  Afile->setPosition(0);

  //  Ask both for the largest EST iid seen, then iterate over those.
  //
  u32bit  largestIID = Afile->maxIID();

  for (u32bit iid=0; iid<largestIID; iid++) {
    sim4polishList *A = Afile->getEST(iid);

    if (A->length() > 0) {
      bool affected = false;

      //  fill out the overlap matrix

      u32bit  **overlap = new u32bit* [A->length()];
      overlap[0] = new u32bit [A->length() * A->length()];
      for (u32bit i=1; i<A->length(); i++)
        overlap[i] = overlap[i-1] + A->length();

      for (u32bit a=0; a<A->length(); a++)
        for (u32bit b=0; b<A->length(); b++)
          overlap[a][b] = findOverlap((*A)[a], (*A)[b]);

      //  For each match, find the ones that overlap.

      for (u32bit a=0; a<A->length(); a++) {
        u32bit longest = 0;

        for (u32bit b=0; b<A->length(); b++) {
          if (overlap[a][longest] < overlap[a][b])
            longest = b;
        }

        //  Emit the match if we are the longest.

        if (a == longest)
          s4p_printPolish(stdout, (*A)[a], 0);
        else {
          affected = true;
          matchesRemoved++;
        }
      }

      if (affected)
        ESTaffected++;

      delete [] overlap[0];
      delete [] overlap;
    }

    delete A;
  }

  fprintf(stderr, "Removed "u32bitFMT" matches, affecting "u32bitFMT" ESTs\n", matchesRemoved, ESTaffected);

  delete Afile;
}


