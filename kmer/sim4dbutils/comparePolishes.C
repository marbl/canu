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



//  For cDNA gaps larger than GAP_MINIMUM, count it as a gap only if
//  the genomic gap is within GAP_DIFFERENCE of the cDNA gap.
//
//  XXX This needs some tweaking!
//
#define GAP_MINIMUM     10
#define GAP_DIFFERENCE   4



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
    sim4polishList *A  = Afile->getEST(iid);
    sim4polishList *B  = Bfile->getEST(iid);
    sim4polishList *Ta = 0L;
    sim4polishList *Tb = 0L;

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



    //  Find and remove those matches that are unique to either set.
    //  Removing is a big pain, because we either have to know
    //  something about the removal process, or we need to rebuild the
    //  overlap matrix after each removal.  Instead, we build a new set.

    bool *removeA = new bool [A->length()];
    bool *removeB = new bool [B->length()];

    for (u32bit a=0; a<A->length(); a++)
      removeA[a] = false;

    for (u32bit b=0; b<B->length(); b++)
      removeB[b] = false;


    for (u32bit a=0; a<A->length(); a++) {
      u32bit ovl = 0;

      for (u32bit b=0; b<B->length(); b++)
        if (overlap[a][b])
          ovl++;

      if (ovl == 0) {
        removeA[a] = true;
        novelInA++;
        //  XXX OUTPUT
      }
    }

    for (u32bit b=0; b<B->length(); b++) {
      u32bit ovl = 0;

      for (u32bit a=0; a<A->length(); a++)
        if (overlap[a][b])
          ovl++;

      if (ovl == 0) {
        removeB[b] = true;
        novelInB++;
        //  XXX OUTPUT
      }
    }    

    //
    //  Now find all those that are perfect matches.  Yeah, yeah, we
    //  could ignore those that we already marked for removal.
    //

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

      //  If exactly one overlap, we just need to check if the guy in B
      //  also has one overlap with anybody in A.

      if (Boverlaps == 1) {

        //  Count the number of overlaps the guy in B has with A.  If
        //  1, it's a goodOverlap, else it's a multipleInA.

        u32bit b = theBovl;

        u32bit Aoverlaps = 0;
        for (u32bit x=0; x<A->length(); x++)
          if (overlap[x][b])
            Aoverlaps++;

        if (Aoverlaps == 1) {
          removeA[a] = true;
          removeB[b] = true;
          goodOverlap++;

          //    ESTiid ESTlen  overlap  A%id A%cov AgenLen #exons #cdnagaps  B%id B%cov BgenLen #exons #cdnagaps

          u32bit AgenLen = 0, BgenLen = 0;
          u32bit Agaps   = 0, Bgaps   = 0;

          for (u32bit x=0; x < (*A)[a]->numExons; x++)
            AgenLen += (*A)[a]->exons[x].genTo - (*A)[a]->exons[x].genFrom + 1;

          for (u32bit x=0; x < (*B)[b]->numExons; x++)
            BgenLen += (*B)[b]->exons[x].genTo - (*B)[b]->exons[x].genFrom + 1;

#ifdef GAP_MINIMUM
          for (u32bit x=1; x < (*A)[a]->numExons; x++) {
            int egap = (*A)[a]->exons[x].estFrom - (*A)[a]->exons[x-1].estTo;
            int ggap = (*A)[a]->exons[x].genFrom - (*A)[a]->exons[x-1].genTo;
            int dgap = 0;

            if (egap > ggap)
              dgap = egap - ggap;
            else
              dgap = ggap - egap;

            if ((egap > GAP_MINIMUM) &&
                (dgap < GAP_DIFFERENCE))
              Agaps++;
          }

          for (u32bit x=1; x < (*B)[b]->numExons; x++) {
            int egap = (*B)[b]->exons[x].estFrom - (*B)[b]->exons[x-1].estTo;
            int ggap = (*B)[b]->exons[x].genFrom - (*B)[b]->exons[x-1].genTo;
            int dgap = 0;

            if (egap > ggap)
              dgap = egap - ggap;
            else
              dgap = ggap - egap;

            if ((egap > GAP_MINIMUM) &&
                (dgap < GAP_DIFFERENCE))
              Bgaps++;
          }
#else
          for (u32bit x=1; x < (*A)[a]->numExons; x++)
            if ( (*A)[a]->exons[x].estFrom - (*A)[a]->exons[x-1].estTo != 1 )
              Agaps++;

          for (u32bit x=1; x < (*B)[b]->numExons; x++)
            if ( (*B)[b]->exons[x].estFrom - (*B)[b]->exons[x-1].estTo != 1 )
              Bgaps++;
#endif

          double score = 0;
          if (AgenLen > BgenLen)
            score = (double)overlap[a][b] / (double)BgenLen;
          else
            score = (double)overlap[a][b] / (double)AgenLen;

          fprintf(outfile, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%f\t%8.3f\t%8.3f\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%8.3f\t%8.3f\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
                  iid,
                  (*A)[a]->estLen,
                  overlap[a][b],
                  score,
                  s4p_percentIdentity( (*A)[a] ),  //->percentIdentity,
                  100.0 * (double)((*A)[a]->numCovered) / (double)((*A)[a]->estLen - (*A)[a]->estPolyA - (*A)[a]->estPolyT),
                  //(*A)[a]->querySeqIdentity,
                  AgenLen, (*A)[a]->numExons, Agaps,
                  s4p_percentIdentity( (*B)[b] ),  //->percentIdentity,
                  100.0 * (double)((*B)[b]->numCovered) / (double)((*B)[b]->estLen - (*B)[b]->estPolyA - (*B)[b]->estPolyT),
                  //(*B)[b]->querySeqIdentity,
                  BgenLen, (*B)[b]->numExons, Bgaps);
                  
          //  XXX OUTPUT
        }
      }
    }

    //
    //  Rebuild
    //

    Ta = new sim4polishList;
    Tb = new sim4polishList;

    for (u32bit a=0; a<A->length(); a++)
      if (removeA[a] == false)
        Ta->push( s4p_copyPolish((*A)[a]) );

    for (u32bit b=0; b<B->length(); b++)
      if (removeB[b] == false)
        Tb->push( s4p_copyPolish((*B)[b]) );

    delete A;
    delete B;
    A = Ta;
    B = Tb;
    Ta = Tb = 0L;

    //  Rebuild overlaps
    //
    for (u32bit a=0; a<A->length(); a++)
      for (u32bit b=0; b<B->length(); b++)
        overlap[a][b] = findOverlap((*A)[a], (*B)[b]);


    //
    //  And now all we're left with is a bunch of intersecting crud.
    //


    //  Grab the first match in A.  Find all the overlaps with things
    //  in B.  For each of those, find the overlaps in A.  Repeat
    //  until nothing changes.  Generate a report.  Remove all those
    //  matches.  Do it all again until there are no more matches.

    while (A->length()) {
      for (u32bit a=0; a<A->length(); a++)
        removeA[a] = false;

      for (u32bit b=0; b<B->length(); b++)
        removeB[b] = false;

      removeA[0] = true;

      bool keepGoing = true;

      while (keepGoing) {
        keepGoing = false;

        //  For all of A, if we have something marked for removal, see if we
        //  overlap with anything in B.  If that b is not marked for removal,
        //  mark it, and keep going.
        //
        for (u32bit a=0; a<A->length(); a++) {
          if (removeA[a]) {
            for (u32bit b=0; b<B->length(); b++) {
              if ((overlap[a][b]) && (removeB[b] == false)) {
                removeB[b] = true;
                keepGoing = true;
              }
            }
          }
        }

        //  Same thing, but for B.
        // 
        for (u32bit b=0; b<B->length(); b++) {
          if (removeB[b]) {
            for (u32bit a=0; a<A->length(); a++) {
              if ((overlap[a][b]) && (removeA[a] == false)) {
                removeA[a] = true;
                keepGoing = true;
              }
            }
          }
        }
      }

      //  Found a component.  Output it.

      u32bit inA = 0;
      u32bit inB = 0;

      for (u32bit a=0; a<A->length(); a++)
        if (removeA[a])
          inA++;
      for (u32bit b=0; b<B->length(); b++)
        if (removeB[b])
          inB++;

      if        ((inA  > 1) && (inB  > 1)) {
        hairyOverlap++;
      } else if ((inA == 1) && (inB  > 1)) {
        multipleInB++;
      } else if ((inA  > 1) && (inB == 1)) {
        multipleInA++;
      } else {
        fprintf(stderr, "ERROR!  inA="u32bitFMT" inB="u32bitFMT"\n", inA, inB);
      }

      //
      //  Rebuild
      //

      Ta = new sim4polishList;
      Tb = new sim4polishList;

      for (u32bit a=0; a<A->length(); a++)
        if (removeA[a] == false)
          Ta->push( s4p_copyPolish((*A)[a]) );

      for (u32bit b=0; b<B->length(); b++)
        if (removeB[b] == false)
          Tb->push( s4p_copyPolish((*B)[b]) );
    
      delete A;
      delete B;
      A = Ta;
      B = Tb;
      Ta = Tb = 0L;

      //  Rebuild overlaps
      //
      for (u32bit a=0; a<A->length(); a++)
        for (u32bit b=0; b<B->length(); b++)
          overlap[a][b] = findOverlap((*A)[a], (*B)[b]);
    }

    if ((iid % 1234) == 0) {
      fprintf(stderr, "IID:"u32bitFMTW(8)"  good:"u32bitFMTW(4)" Anovel:"u32bitFMTW(4)" Amulti:"u32bitFMTW(4)" Bnovel:"u32bitFMTW(4)" Bmulti:"u32bitFMTW(4)" hairy:"u32bitFMTW(4)"\r",
              iid,
              goodOverlap, novelInA, multipleInA, novelInB, multipleInB, hairyOverlap);
      fflush(stderr);
    }

    delete [] overlap[0];
    delete [] overlap;

    delete [] removeA;
    delete [] removeB;

    delete A;
    delete B;
  }

  fprintf(stderr, "\ngood:"u32bitFMTW(4)" Anovel:"u32bitFMTW(4)" Amulti:"u32bitFMTW(4)" Bnovel:"u32bitFMTW(4)" Bmulti:"u32bitFMTW(4)" hairy:"u32bitFMTW(4)"\n",
          goodOverlap, novelInA, multipleInA, novelInB, multipleInB, hairyOverlap);


  delete Afile;
  delete Bfile;
}
