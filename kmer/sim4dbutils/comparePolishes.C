#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bio++.H"
#include "sim4.H"
#include "s4p_overlap.H"

//  Matches two sets of polishes to each other using a simple overlap
//  heuristic.
//
//  Arguments (are horrible, whatcha gonna do about it?):
//
//    -i min percent id (default 95)
//    -c min percent coverage (default 50)
//    -a polishes input file 1
//    -b polishes input file 2
//    -gff3 write output as GFF3
//
//    Output is on standard out, and is tab-delimited.  It reports
//    stuff about the 'same' matches:
//
//    ESTiid ESTlen  overlap  A%id A%cov #cdnagaps #exons  B%id B%cov #cdnagaps #exons

//  Try to analyze cDNA gaps.
//
//  For cDNA gaps larger than GAP_MINIMUM, count it as a gap only if
//  the genomic gap is within GAP_DIFFERENCE of the cDNA gap.
//
//  XXX This needs some tweaking!
//
#define GAP_MINIMUM     10
#define GAP_DIFFERENCE   4


sim4polishWriter *
openOutput(const char *prefix, const char *suffix, sim4polishStyle style) {
  char name[FILENAME_MAX];
  sprintf(name, "%s.%s", prefix, suffix);
  return(new sim4polishWriter(name, style));
}



int
main(int argc, char **argv) {
  u32bit           minI   = 95;
  u32bit           minC   = 50;
  const char      *prefix = "comparePolishes";
  sim4polishFile  *Afile  = 0L;
  sim4polishFile  *Bfile  = 0L;

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

  bool             doGFF3;

  sim4polishStyle Astyle = sim4polishStyleDefault;
  sim4polishStyle Bstyle = sim4polishStyleDefault;

  sim4polishStyle style  = sim4polishStyleDefault;


  int arg=1;
  while(arg < argc) {
    if        (strcmp(argv[arg], "-i") == 0) {
      minI = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-c") == 0) {
      minC = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-a") == 0) {
      //  Ugly hack to obtain the style of the input files, but can be fixed later

      sim4polishReader *AR = new sim4polishReader(argv[++arg]);
      Astyle = AR->getsim4polishStyle();
      delete   AR;     

      Afile = new sim4polishFile(argv[arg], Astyle);
    } else if (strcmp(argv[arg], "-b") == 0) {
      //  Ugly hack to obtain the style of the input files, but can be fixed later

      sim4polishReader *BR = new sim4polishReader(argv[++arg]);
      Bstyle = BR->getsim4polishStyle();
      delete   BR;

      Bfile = new sim4polishFile(argv[arg], Astyle);
    } else if (strcmp(argv[arg], "-p") == 0) {
      prefix = argv[++arg];
    } else if (strcmp(argv[arg], "-gff3") == 0) {
      doGFF3 = true;
      style = sim4polishGFF3;
    }
    arg++;
  }

  if ((Afile == 0L) || (Bfile == 0L)) {
    fprintf(stderr, "usage: %s [-i percent-identity] [-c percent-coverage] -a input-set-a -b input-set-b [-p output-prefix] [-gff3]\n", argv[0]);
    fprintf(stderr, "only -a and -b are mandatory, but you should give all anyway\n");
    exit(1);
  }

  //  Open the output files
  //
  sim4polishWriter *fasame  = openOutput(prefix, "a-same", style);
  sim4polishWriter *fbsame  = openOutput(prefix, "b-same", style);
  sim4polishWriter *fanovel = openOutput(prefix, "a-novel", style);
  sim4polishWriter *fbnovel = openOutput(prefix, "b-novel", style);
  sim4polishWriter *famulti = openOutput(prefix, "a-multi", style);
  sim4polishWriter *fbmulti = openOutput(prefix, "b-multi", style);
  sim4polishWriter *fhairy  = openOutput(prefix, "hairy", style);

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

    olap_t  **overlap = new olap_t* [A->length()];
    overlap[0] = new olap_t [A->length() * B->length()];
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

        if (fanovel)
          fanovel->writeAlignment((*A)[a]);
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

        if (fbnovel)
          fbnovel->writeAlignment((*B)[b]);
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

          for (u32bit x=0; x < (*A)[a]->_numExons; x++)
            AgenLen += (*A)[a]->_exons[x]._genTo - (*A)[a]->_exons[x]._genFrom + 1;

          for (u32bit x=0; x < (*B)[b]->_numExons; x++)
            BgenLen += (*B)[b]->_exons[x]._genTo - (*B)[b]->_exons[x]._genFrom + 1;

#ifdef GAP_MINIMUM
          for (u32bit x=1; x < (*A)[a]->_numExons; x++) {
            int egap = (*A)[a]->_exons[x]._estFrom - (*A)[a]->_exons[x-1]._estTo;
            int ggap = (*A)[a]->_exons[x]._genFrom - (*A)[a]->_exons[x-1]._genTo;
            int dgap = 0;

            if (egap > ggap)
              dgap = egap - ggap;
            else
              dgap = ggap - egap;

            if ((egap > GAP_MINIMUM) &&
                (dgap < GAP_DIFFERENCE))
              Agaps++;
          }

          for (u32bit x=1; x < (*B)[b]->_numExons; x++) {
            int egap = (*B)[b]->_exons[x]._estFrom - (*B)[b]->_exons[x-1]._estTo;
            int ggap = (*B)[b]->_exons[x]._genFrom - (*B)[b]->_exons[x-1]._genTo;
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
          for (u32bit x=1; x < (*A)[a]->_numExons; x++)
            if ( (*A)[a]->_exons[x]._estFrom - (*A)[a]->_exons[x-1]._estTo != 1 )
              Agaps++;

          for (u32bit x=1; x < (*B)[b]->_numExons; x++)
            if ( (*B)[b]->_exons[x]._estFrom - (*B)[b]->_exons[x-1]._estTo != 1 )
              Bgaps++;
#endif

          double score = 0;
          if (AgenLen > BgenLen)
            score = (double)overlap[a][b] / (double)BgenLen;
          else
            score = (double)overlap[a][b] / (double)AgenLen;

          fprintf(stdout, u32bitFMT"\t"u32bitFMT"\t"OLAPTFMT"\t%f\t%8.3f\t%8.3f\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t%8.3f\t%8.3f\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
                  iid,
                  (*A)[a]->_estLen,
                  overlap[a][b],
                  score,
                  (*A)[a]->s4p_percentIdentityExact(),
                  (*A)[a]->s4p_percentCoverageExact(),
                  AgenLen, (*A)[a]->_numExons, Agaps,
                  (*B)[b]->s4p_percentIdentityExact(),
                  (*B)[b]->s4p_percentCoverageExact(),
                  BgenLen, (*B)[b]->_numExons, Bgaps);
                  
          if (fasame)
            fasame->writeAlignment((*A)[a]);
          if (fbsame)
            fbsame->writeAlignment((*B)[b]);
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
        Ta->push(new sim4polish((*A)[a]));

    for (u32bit b=0; b<B->length(); b++)
      if (removeB[b] == false)
        Tb->push(new sim4polish((*B)[b]));

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

        //fprintf(fhairy, "EST="u32bitFMT" "u32bitFMT" "u32bitFMT"\n", (*A)[0]->_estID, inA, inB);
        for (u32bit a=0; a<A->length(); a++)
          if (removeA[a])
            fhairy->writeAlignment((*A)[a]);
        for (u32bit b=0; b<B->length(); b++)
          if (removeB[b])
            fhairy->writeAlignment((*B)[b]);
      } else if ((inA == 1) && (inB  > 1)) {
        multipleInB++;

        //fprintf(fbmulti, "EST="u32bitFMT" "u32bitFMT" "u32bitFMT"\n", (*A)[0]->_estID, inA, inB);
        for (u32bit a=0; a<A->length(); a++)
          if (removeA[a])
            fbmulti->writeAlignment((*A)[a]);
        for (u32bit b=0; b<B->length(); b++)
          if (removeB[b])
            fbmulti->writeAlignment((*B)[b]);
      } else if ((inA  > 1) && (inB == 1)) {
        multipleInA++;

        //fprintf(famulti, "EST="u32bitFMT" "u32bitFMT" "u32bitFMT"\n", (*A)[0]->_estID, inA, inB);
        for (u32bit a=0; a<A->length(); a++)
          if (removeA[a])
            famulti->writeAlignment((*A)[a]);
        for (u32bit b=0; b<B->length(); b++)
          if (removeB[b])
            famulti->writeAlignment((*B)[b]);
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
          Ta->push(new sim4polish((*A)[a]));

      for (u32bit b=0; b<B->length(); b++)
        if (removeB[b] == false)
          Tb->push(new sim4polish((*B)[b]));
    
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

    if ((iid % 100) == 0) {
      fprintf(stderr, "IID:"u32bitFMTW(8)"  good:"u32bitFMTW(4)" Anovel:"u32bitFMTW(4)" Amulti:"u32bitFMTW(4)" Bnovel:"u32bitFMTW(4)" Bmulti:"u32bitFMTW(4)" hairy:"u32bitFMTW(4)"\r",
              iid,
              goodOverlap, novelInA, multipleInA, novelInB, multipleInB, hairyOverlap);
      fflush(stderr);
    }

#if 0
    if ((iid % 1234) == 0) {
      fprintf(stderr, "IID:"u32bitFMTW(8)"  good:"u32bitFMTW(4)" Anovel:"u32bitFMTW(4)" Amulti:"u32bitFMTW(4)" Bnovel:"u32bitFMTW(4)" Bmulti:"u32bitFMTW(4)" hairy:"u32bitFMTW(4)"\r",
              iid,
              goodOverlap, novelInA, multipleInA, novelInB, multipleInB, hairyOverlap);
      fflush(stderr);
    }
#endif

    delete [] overlap[0];
    delete [] overlap;

    delete [] removeA;
    delete [] removeB;

    delete A;
    delete B;
  }

  delete fasame;
  delete fbsame;
  delete fanovel;
  delete fbnovel;
  delete famulti;
  delete fbmulti;
  delete fhairy;

  delete Afile;
  delete Bfile;

  fprintf(stderr, "\ngood:"u32bitFMTW(4)" Anovel:"u32bitFMTW(4)" Amulti:"u32bitFMTW(4)" Bnovel:"u32bitFMTW(4)" Bmulti:"u32bitFMTW(4)" hairy:"u32bitFMTW(4)"\n",
          goodOverlap, novelInA, multipleInA, novelInB, multipleInB, hairyOverlap);

  exit(0);
}
